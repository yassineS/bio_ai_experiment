package samtools

import (
	"bytes"
	"io"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// unsortedSAM has records mixed across (chr1 / chr2) with positions out of
// order; sorting by coordinate must produce a (refID, pos)-ordered output.
const unsortedSAM = `@HD	VN:1.6	SO:unsorted
@SQ	SN:chr1	LN:1000
@SQ	SN:chr2	LN:1000
r10	0	chr2	50	60	5M	*	0	0	ACGTA	IIIII	NM:i:3
r2	0	chr1	200	60	5M	*	0	0	ACGTA	IIIII	NM:i:1
r1	0	chr1	100	60	5M	*	0	0	ACGTA	IIIII	NM:i:2
r20	0	chr1	150	60	5M	*	0	0	ACGTA	IIIII	NM:i:9
u1	4	*	0	0	*	*	0	0	*	*
`

// readBAMRecords drains a BAM-encoded buffer into a []*sam.Record.
func readBAMRecords(t *testing.T, raw []byte) []*sam.Record {
	t.Helper()
	br, err := sam.NewBAMReader(bytes.NewReader(raw))
	if err != nil {
		t.Fatalf("NewBAMReader: %v", err)
	}
	var out []*sam.Record
	for {
		rec, err := br.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("BAM Read: %v", err)
		}
		out = append(out, rec)
	}
	return out
}

func TestSortByCoordinate(t *testing.T) {
	var out bytes.Buffer
	if err := Sort(strings.NewReader(unsortedSAM), &out, SortOptions{Order: SortCoordinate, OutputBAM: true}); err != nil {
		t.Fatalf("Sort: %v", err)
	}
	recs := readBAMRecords(t, out.Bytes())
	if len(recs) != 5 {
		t.Fatalf("record count: got %d, want 5", len(recs))
	}
	// Expected order: chr1 100 (r1), chr1 150 (r20), chr1 200 (r2),
	// chr2 50 (r10), unmapped u1 last.
	wantQNames := []string{"r1", "r20", "r2", "r10", "u1"}
	for i, w := range wantQNames {
		if recs[i].QName != w {
			t.Errorf("position %d: got %q, want %q", i, recs[i].QName, w)
		}
	}
}

func TestSortByCoordinateSAMOutput(t *testing.T) {
	var out bytes.Buffer
	if err := Sort(strings.NewReader(unsortedSAM), &out, SortOptions{Order: SortCoordinate, OutputSAM: true}); err != nil {
		t.Fatalf("Sort: %v", err)
	}
	got := out.String()
	if !strings.Contains(got, "@HD") {
		t.Errorf("expected @HD in SAM output: %q", got[:50])
	}
	if !strings.Contains(got, "SO:coordinate") {
		t.Errorf("expected SO:coordinate in @HD: %q", got[:80])
	}
	// First record line (after the header) should be r1.
	lines := strings.Split(strings.TrimRight(got, "\n"), "\n")
	firstBody := ""
	for _, l := range lines {
		if !strings.HasPrefix(l, "@") {
			firstBody = l
			break
		}
	}
	if !strings.HasPrefix(firstBody, "r1\t") {
		t.Errorf("first body line: got %q, want r1...", firstBody)
	}
}

func TestSortByName(t *testing.T) {
	var out bytes.Buffer
	if err := Sort(strings.NewReader(unsortedSAM), &out, SortOptions{Order: SortByName, OutputBAM: true}); err != nil {
		t.Fatalf("Sort: %v", err)
	}
	recs := readBAMRecords(t, out.Bytes())
	got := make([]string, len(recs))
	for i, r := range recs {
		got[i] = r.QName
	}
	// Plain lexicographic: r1 < r10 < r2 < r20 < u1.
	want := []string{"r1", "r10", "r2", "r20", "u1"}
	for i := range want {
		if got[i] != want[i] {
			t.Errorf("position %d: got %q, want %q", i, got[i], want[i])
		}
	}
}

func TestSortByNaturalName(t *testing.T) {
	var out bytes.Buffer
	if err := Sort(strings.NewReader(unsortedSAM), &out, SortOptions{Order: SortByNameNatural, OutputBAM: true}); err != nil {
		t.Fatalf("Sort: %v", err)
	}
	recs := readBAMRecords(t, out.Bytes())
	got := make([]string, len(recs))
	for i, r := range recs {
		got[i] = r.QName
	}
	// Natural ordering: r1 < r2 < r10 < r20 < u1.
	want := []string{"r1", "r2", "r10", "r20", "u1"}
	for i := range want {
		if got[i] != want[i] {
			t.Errorf("position %d: got %q, want %q", i, got[i], want[i])
		}
	}
}

func TestSortByTag(t *testing.T) {
	var out bytes.Buffer
	if err := Sort(strings.NewReader(unsortedSAM), &out, SortOptions{Order: SortByTag, Tag: "NM", OutputBAM: true}); err != nil {
		t.Fatalf("Sort: %v", err)
	}
	recs := readBAMRecords(t, out.Bytes())
	// NM values: r10=3, r2=1, r1=2, r20=9, u1=missing. Upstream sorts the
	// tagless record first, then by NM value: u1, r2(1), r1(2), r10(3),
	// r20(9). Verified against `samtools sort -t NM`.
	want := []string{"u1", "r2", "r1", "r10", "r20"}
	for i, w := range want {
		if recs[i].QName != w {
			t.Errorf("position %d: got %q, want %q", i, recs[i].QName, w)
		}
	}
}

func TestSortByTagMissing(t *testing.T) {
	// Sorting by a tag absent on every record falls back to the coordinate
	// secondary key (refID, pos), matching upstream `samtools sort -t ZZ`:
	// r1(chr1:100), r20(chr1:150), r2(chr1:200), r10(chr2:50), u1(unmapped).
	var out bytes.Buffer
	if err := Sort(strings.NewReader(unsortedSAM), &out, SortOptions{Order: SortByTag, Tag: "ZZ", OutputBAM: true}); err != nil {
		t.Fatalf("Sort: %v", err)
	}
	recs := readBAMRecords(t, out.Bytes())
	want := []string{"r1", "r20", "r2", "r10", "u1"}
	for i, w := range want {
		if recs[i].QName != w {
			t.Errorf("position %d on missing tag: got %q, want %q", i, recs[i].QName, w)
		}
	}
}

// TestSortMultiShard forces external-merge sort by setting an absurdly tiny
// memory budget so every record spills to its own shard, then verifies the
// final order is correct.
func TestSortMultiShard(t *testing.T) {
	var out bytes.Buffer
	opts := SortOptions{
		Order:       SortCoordinate,
		OutputBAM:   true,
		MaxMemBytes: 1, // forces a flush after every record
		TmpPrefix:   t.TempDir(),
	}
	if err := Sort(strings.NewReader(unsortedSAM), &out, opts); err != nil {
		t.Fatalf("Sort: %v", err)
	}
	recs := readBAMRecords(t, out.Bytes())
	want := []string{"r1", "r20", "r2", "r10", "u1"}
	if len(recs) != len(want) {
		t.Fatalf("record count: got %d, want %d", len(recs), len(want))
	}
	for i, w := range want {
		if recs[i].QName != w {
			t.Errorf("position %d: got %q, want %q", i, recs[i].QName, w)
		}
	}
}

// TestSortPackedMatchesDecode is the byte-identity gate for the packed-spill
// fast path: for every sort mode (coordinate, -N, -n, -t) and both the in-memory
// and the every-record-spills configurations, sorting a BAM (packed path) must
// produce the SAME final BAM bytes as sorting the equivalent SAM (decode path).
// A comparator or spill/merge mismatch surfaces here as a byte diff.
func TestSortPackedMatchesDecode(t *testing.T) {
	bamIn := samToBAM(t, unsortedSAM)

	modes := []struct {
		name string
		opts SortOptions
	}{
		{"coordinate", SortOptions{Order: SortCoordinate, OutputBAM: true}},
		{"name-lex", SortOptions{Order: SortByName, OutputBAM: true}},
		{"name-natural", SortOptions{Order: SortByNameNatural, OutputBAM: true}},
		{"tag-NM", SortOptions{Order: SortByTag, Tag: "NM", OutputBAM: true}},
		{"tag-missing", SortOptions{Order: SortByTag, Tag: "ZZ", OutputBAM: true}},
	}
	for _, m := range modes {
		for _, spill := range []bool{false, true} {
			t.Run(m.name+map[bool]string{false: "/in-mem", true: "/spill"}[spill], func(t *testing.T) {
				decOpts := m.opts
				pkOpts := m.opts
				if spill {
					decOpts.MaxMemBytes = 1
					decOpts.TmpPrefix = t.TempDir()
					pkOpts.MaxMemBytes = 1
					pkOpts.TmpPrefix = t.TempDir()
				}
				var decOut, pkOut bytes.Buffer
				if err := Sort(strings.NewReader(unsortedSAM), &decOut, decOpts); err != nil {
					t.Fatalf("decode-path Sort: %v", err)
				}
				if err := Sort(bytes.NewReader(bamIn), &pkOut, pkOpts); err != nil {
					t.Fatalf("packed-path Sort: %v", err)
				}
				if !bytes.Equal(decOut.Bytes(), pkOut.Bytes()) {
					// Fall back to a record-level diff for a readable failure.
					dr := readBAMRecords(t, decOut.Bytes())
					pr := readBAMRecords(t, pkOut.Bytes())
					if len(dr) != len(pr) {
						t.Fatalf("record count: decode=%d packed=%d", len(dr), len(pr))
					}
					for i := range dr {
						if dr[i].QName != pr[i].QName || dr[i].RName != pr[i].RName || dr[i].Pos != pr[i].Pos || dr[i].Flag != pr[i].Flag {
							t.Fatalf("record %d differs: decode=%s/%s/%d packed=%s/%s/%d",
								i, dr[i].QName, dr[i].RName, dr[i].Pos, pr[i].QName, pr[i].RName, pr[i].Pos)
						}
					}
					t.Fatalf("packed BAM bytes differ from decode path (records match field-wise; header/encoding diff)")
				}
			})
		}
	}
}

// TestSortPackedSAMOutput verifies the packed path's SAM-output branch (raw
// bodies decoded on the way out) matches the decode path's SAM output exactly.
func TestSortPackedSAMOutput(t *testing.T) {
	bamIn := samToBAM(t, unsortedSAM)
	for _, order := range []SortOrder{SortCoordinate, SortByName, SortByNameNatural} {
		var decOut, pkOut bytes.Buffer
		opts := SortOptions{Order: order, OutputSAM: true}
		if err := Sort(strings.NewReader(unsortedSAM), &decOut, opts); err != nil {
			t.Fatalf("decode SAM Sort: %v", err)
		}
		if err := Sort(bytes.NewReader(bamIn), &pkOut, opts); err != nil {
			t.Fatalf("packed SAM Sort: %v", err)
		}
		if decOut.String() != pkOut.String() {
			t.Errorf("order %d: packed SAM output differs from decode path:\n--decode--\n%s\n--packed--\n%s", order, decOut.String(), pkOut.String())
		}
	}
}

func TestSortEmptyTag(t *testing.T) {
	if err := Sort(strings.NewReader(unsortedSAM), &bytes.Buffer{}, SortOptions{Order: SortByTag, Tag: ""}); err != ErrEmptyTag {
		t.Errorf("expected ErrEmptyTag, got %v", err)
	}
	// Tag with wrong length is also rejected.
	if err := Sort(strings.NewReader(unsortedSAM), &bytes.Buffer{}, SortOptions{Order: SortByTag, Tag: "X"}); err != ErrEmptyTag {
		t.Errorf("expected ErrEmptyTag for short tag, got %v", err)
	}
}

func TestParseSortOrder(t *testing.T) {
	cases := []struct {
		in   string
		want SortOrder
		bad  bool
	}{
		{"", SortCoordinate, false},
		{"coordinate", SortCoordinate, false},
		{"queryname", SortByName, false},
		{"name", SortByName, false},
		{"natural", SortByNameNatural, false},
		{"tag", SortByTag, false},
		{"xyz", 0, true},
	}
	for _, c := range cases {
		got, err := ParseSortOrder(c.in)
		if c.bad {
			if err == nil {
				t.Errorf("ParseSortOrder(%q): expected error", c.in)
			}
			continue
		}
		if err != nil {
			t.Errorf("ParseSortOrder(%q): %v", c.in, err)
		}
		if got != c.want {
			t.Errorf("ParseSortOrder(%q): got %v, want %v", c.in, got, c.want)
		}
	}
}

func TestParseMemBudget(t *testing.T) {
	cases := []struct {
		in   string
		want int64
		bad  bool
	}{
		{"", 0, false},
		{"100", 100, false},
		{"4K", 4 << 10, false},
		{"4k", 4 << 10, false},
		{"8M", 8 << 20, false},
		{"2G", 2 << 30, false},
		{"-1", 0, true},
		{"abc", 0, true},
	}
	for _, c := range cases {
		got, err := ParseMemBudget(c.in)
		if c.bad {
			if err == nil {
				t.Errorf("ParseMemBudget(%q): expected error", c.in)
			}
			continue
		}
		if err != nil {
			t.Errorf("ParseMemBudget(%q): %v", c.in, err)
		}
		if got != c.want {
			t.Errorf("ParseMemBudget(%q): got %d, want %d", c.in, got, c.want)
		}
	}
}

// emptyRefIndex is the reference map used by the tagLess unit tests; the
// records they build are unmapped, so the coordinate secondary key resolves
// every one to the same unmapped sentinel.
var emptyRefIndex = map[string]int{"chr1": 0}

// TestTagLessFloat exercises the float branch of tagLess: floats compare
// numerically and equal floats fall through to the coordinate secondary key.
func TestTagLessFloat(t *testing.T) {
	a := &sam.Record{QName: "a", RName: "chr1", Pos: 10, Aux: []sam.Aux{{Tag: "XF", Type: 'f', Value: 1.0}}}
	b := &sam.Record{QName: "b", RName: "chr1", Pos: 20, Aux: []sam.Aux{{Tag: "XF", Type: 'f', Value: 2.0}}}
	if !tagLess(a, b, "XF", emptyRefIndex) {
		t.Error("XF=1.0 should be less than XF=2.0")
	}
	if tagLess(b, a, "XF", emptyRefIndex) {
		t.Error("XF=2.0 should NOT be less than XF=1.0")
	}
	// Equal floats fall back to the coordinate key (pos 10 < pos 30).
	c := &sam.Record{QName: "z", RName: "chr1", Pos: 30, Aux: []sam.Aux{{Tag: "XF", Type: 'f', Value: 1.0}}}
	if !tagLess(a, c, "XF", emptyRefIndex) {
		t.Error("tie on XF=1.0 should defer to coordinate (pos 10 < 30)")
	}
}

// TestTagLessString exercises the string branch of tagLess.
func TestTagLessString(t *testing.T) {
	a := &sam.Record{QName: "a", RName: "chr1", Pos: 10, Aux: []sam.Aux{{Tag: "RG", Type: 'Z', Value: "rg-aa"}}}
	b := &sam.Record{QName: "b", RName: "chr1", Pos: 20, Aux: []sam.Aux{{Tag: "RG", Type: 'Z', Value: "rg-bb"}}}
	if !tagLess(a, b, "RG", emptyRefIndex) {
		t.Error("rg-aa < rg-bb expected")
	}
	// Equal strings fall back to the coordinate key (pos 10 < pos 30).
	c := &sam.Record{QName: "z", RName: "chr1", Pos: 30, Aux: []sam.Aux{{Tag: "RG", Type: 'Z', Value: "rg-aa"}}}
	if !tagLess(a, c, "RG", emptyRefIndex) {
		t.Error("tie on rg-aa should defer to coordinate (pos 10 < 30)")
	}
}

// TestTagLessMissingOneSide verifies upstream's rule that records lacking the
// sort tag sort *before* records that carry it (bam1_cmp_by_tag returns -1
// when only the left record's tag is NULL).
func TestTagLessMissingOneSide(t *testing.T) {
	withTag := &sam.Record{QName: "a", Aux: []sam.Aux{{Tag: "NM", Type: 'i', Value: int64(1)}}}
	without := &sam.Record{QName: "b"}
	// The tagless record sorts first.
	if !tagLess(without, withTag, "NM", emptyRefIndex) {
		t.Error("record without tag should sort before record with tag")
	}
	if tagLess(withTag, without, "NM", emptyRefIndex) {
		t.Error("record with tag should NOT sort before record without tag")
	}
}

// TestSortWriteShardError forces writeShard to fail by pointing TmpPrefix
// at a path the OS refuses to create files in.
func TestSortWriteShardErrorPath(t *testing.T) {
	// Use a path under a non-existent directory; writeShard's os.Create
	// must fail.
	bad := "/this/does/not/exist/sort"
	opts := SortOptions{
		Order:       SortCoordinate,
		OutputBAM:   true,
		MaxMemBytes: 1,
		TmpPrefix:   bad,
	}
	err := Sort(strings.NewReader(unsortedSAM), &bytes.Buffer{}, opts)
	if err == nil {
		t.Error("expected error from unwritable TmpPrefix")
	}
}

// TestSortBadInput rejects input that cannot be parsed as SAM.
func TestSortBadInput(t *testing.T) {
	bad := "@SQ\tbroken\n"
	if err := Sort(strings.NewReader(bad), &bytes.Buffer{}, SortOptions{}); err == nil {
		t.Error("expected error on bad input")
	}
}

func TestNaturalLess(t *testing.T) {
	cases := []struct {
		a, b string
		want bool
	}{
		{"r1", "r2", true},
		{"r2", "r10", true},
		{"r10", "r2", false},
		{"abc", "abd", true},
		{"a1b", "a1c", true},
		{"a01", "a1", false}, // same numeric value: shorter literal run sorts first
		{"", "x", true},
	}
	for _, c := range cases {
		if got := naturalLess(c.a, c.b); got != c.want {
			t.Errorf("naturalLess(%q, %q): got %v, want %v", c.a, c.b, got, c.want)
		}
	}
}

// TestStrnumCmp checks the three-way strnum_cmp port in both natural and
// plain modes, including symmetry and the equal-value/equal-string cases.
func TestStrnumCmp(t *testing.T) {
	sign := func(n int) int {
		switch {
		case n < 0:
			return -1
		case n > 0:
			return 1
		default:
			return 0
		}
	}
	cases := []struct {
		a, b    string
		natural bool
		want    int
	}{
		// Natural mode: numeric runs compare by value.
		{"r1", "r2", true, -1},
		{"r2", "r10", true, -1},
		{"r10", "r2", true, 1},
		{"r10", "r10", true, 0},
		{"a01", "a1", true, 0}, // equal numeric value
		{"abc", "abd", true, -1},
		{"", "x", true, -1},
		{"x", "", true, 1},
		// Plain (lexicographic) mode: q10 sorts before q2.
		{"q10", "q2", false, -1},
		{"q2", "q10", false, 1},
		{"q2", "q2", false, 0},
	}
	for _, c := range cases {
		got := sign(strnumCmp(c.a, c.b, c.natural))
		if got != c.want {
			t.Errorf("strnumCmp(%q, %q, natural=%v): got %d, want %d", c.a, c.b, c.natural, got, c.want)
		}
		// The comparator must be antisymmetric.
		rev := sign(strnumCmp(c.b, c.a, c.natural))
		if rev != -c.want {
			t.Errorf("strnumCmp(%q, %q, natural=%v) not antisymmetric: got %d, want %d", c.b, c.a, c.natural, rev, -c.want)
		}
	}
}

// TestFlagSortKey pins the exact integer keys the FLAG-to-key transform
// produces (a direct port of upstream's bit shuffle). The resulting numeric
// ordering across the canonical FLAG categories is
// primary < supplementary < secondary < read1 < read2, which is what the
// name-sort secondary key relies on.
func TestFlagSortKey(t *testing.T) {
	cases := []struct {
		flag uint16
		want int
	}{
		{0, 0},                         // primary, neither read1 nor read2
		{sam.FlagSupplementary, 0x100}, // 0x800 >> 3
		{sam.FlagSecondary, 0x800},     // 0x100 << 3
		{sam.FlagRead1, 0x4000},        // 0x40 << 8
		{sam.FlagRead2, 0x8000},        // 0x80 << 8
		{sam.FlagRead1 | sam.FlagRead2, 0xc000},
	}
	for _, c := range cases {
		if got := flagSortKey(c.flag); got != c.want {
			t.Errorf("flagSortKey(0x%x): got 0x%x, want 0x%x", c.flag, got, c.want)
		}
	}
	// The numeric ordering must be primary < supp < sec < read1 < read2.
	order := []uint16{0, sam.FlagSupplementary, sam.FlagSecondary, sam.FlagRead1, sam.FlagRead2}
	for i := 1; i < len(order); i++ {
		if !(flagSortKey(order[i-1]) < flagSortKey(order[i])) {
			t.Errorf("FLAG key order broken at index %d", i)
		}
	}
	// Strand and other non-category bits must not perturb the key.
	if flagSortKey(sam.FlagRead1) != flagSortKey(sam.FlagRead1|sam.FlagReverse|sam.FlagPaired) {
		t.Error("non-category flag bits must not affect the name-sort key")
	}
}

// TestNameLessFlagTiebreak exercises nameLess's FLAG secondary key on a
// QNAME tie. The numeric key ordering means a primary read sorts before a
// supplementary, which sorts before a secondary, and read1 sorts before
// read2.
func TestNameLessFlagTiebreak(t *testing.T) {
	primary := &sam.Record{QName: "q", Flag: 0}
	supp := &sam.Record{QName: "q", Flag: sam.FlagSupplementary}
	sec := &sam.Record{QName: "q", Flag: sam.FlagSecondary}
	read1 := &sam.Record{QName: "q", Flag: sam.FlagPaired | sam.FlagRead1}
	read2 := &sam.Record{QName: "q", Flag: sam.FlagPaired | sam.FlagRead2}
	chain := []*sam.Record{primary, supp, sec, read1, read2}
	for i := 1; i < len(chain); i++ {
		if !nameLess(chain[i-1], chain[i], true) {
			t.Errorf("nameLess: record %d (flag 0x%x) should sort before record %d (flag 0x%x)",
				i-1, chain[i-1].Flag, i, chain[i].Flag)
		}
		if nameLess(chain[i], chain[i-1], true) {
			t.Errorf("nameLess not antisymmetric at index %d", i)
		}
	}
	// Distinct QNAMEs ignore the FLAG key entirely.
	other := &sam.Record{QName: "p", Flag: sam.FlagSecondary}
	if !nameLess(other, read1, true) {
		t.Error("QNAME p should sort before q regardless of FLAG")
	}
}

// TestUnitCoordSortKeyTieBreak pins bug #3's fix: the coordinate sort key is
// (refID, pos, reverse-strand-flag) — NOT a QNAME tie-break. It mirrors
// upstream bam_sort.c bam1_cmp_core's non-QueryName branch: forward-strand
// reads sort before reverse-strand reads at the same (rname, pos), and
// records identical in all three keys compare equal (the stable sort then
// preserves input order, matching htslib's observed output). This is a pure
// helper test with no external binary.
func TestUnitCoordSortKeyTieBreak(t *testing.T) {
	refIndex := map[string]int{"chr1": 0, "chr2": 1}

	fwd := &sam.Record{QName: "z", RName: "chr1", Pos: 100, Flag: 0}
	rev := &sam.Record{QName: "a", RName: "chr1", Pos: 100, Flag: sam.FlagReverse}

	// Same (rname, pos): forward sorts before reverse, regardless of QNAME
	// (here the reverse read's QNAME "a" sorts lexically *before* "z", which
	// would have won under the old buggy QNAME tie-break).
	if !coordLess(fwd, rev, refIndex) {
		t.Error("forward read should sort before reverse read at equal (rname,pos)")
	}
	if coordLess(rev, fwd, refIndex) {
		t.Error("coordLess not antisymmetric on the reverse-strand tie-break")
	}

	// QNAME must NOT influence ordering when (rname, pos, strand) are equal.
	x := &sam.Record{QName: "zzz", RName: "chr1", Pos: 100, Flag: 0}
	y := &sam.Record{QName: "aaa", RName: "chr1", Pos: 100, Flag: 0}
	if coordLess(x, y, refIndex) || coordLess(y, x, refIndex) {
		t.Error("records equal in (rname,pos,strand) must compare equal regardless of QNAME")
	}
	if coordCmp(x, y, refIndex) != 0 {
		t.Errorf("coordCmp of two equal-key records = %d, want 0", coordCmp(x, y, refIndex))
	}

	// Primary key ordering: refID then pos.
	earlierPos := &sam.Record{QName: "p", RName: "chr1", Pos: 50, Flag: sam.FlagReverse}
	if !coordLess(earlierPos, fwd, refIndex) {
		t.Error("smaller pos must sort first even if it is reverse-strand")
	}
	chr2 := &sam.Record{QName: "p", RName: "chr2", Pos: 1, Flag: 0}
	if !coordLess(fwd, chr2, refIndex) {
		t.Error("chr1 must sort before chr2")
	}

	// Unmapped (no RNAME) sorts after every mapped record.
	unmapped := &sam.Record{QName: "u", RName: "*", Pos: 0, Flag: sam.FlagUnmapped}
	if !coordLess(chr2, unmapped, refIndex) {
		t.Error("mapped record must sort before unmapped")
	}
}
