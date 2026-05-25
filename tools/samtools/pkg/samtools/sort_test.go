package samtools

import (
	"bytes"
	"io"
	"strconv"
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
	if err := Sort(strings.NewReader(unsortedSAM), &out, SortOptions{
		Order:           SortByTag,
		Tag:             "NM",
		SecondaryByName: true,
		OutputBAM:       true,
	}); err != nil {
		t.Fatalf("Sort: %v", err)
	}
	recs := readBAMRecords(t, out.Bytes())
	// Upstream `samtools sort -n -t NM`: records missing the tag sort
	// FIRST (per bam1_cmp_by_tag), then the rest by tag value with a
	// qname+FLAG secondary key.
	// NM: r10=3, r2=1, r1=2, r20=9; u1 missing.
	want := []string{"u1", "r2", "r1", "r10", "r20"}
	for i, w := range want {
		if recs[i].QName != w {
			t.Errorf("position %d: got %q, want %q", i, recs[i].QName, w)
		}
	}
}

func TestSortByTagMissing(t *testing.T) {
	// Sorting by a tag absent on every record falls back to the secondary
	// comparator. With SecondaryByName=true the order is qname-natural,
	// so r1 lands first.
	var out bytes.Buffer
	if err := Sort(strings.NewReader(unsortedSAM), &out, SortOptions{
		Order:           SortByTag,
		Tag:             "ZZ",
		SecondaryByName: true,
		OutputBAM:       true,
	}); err != nil {
		t.Fatalf("Sort: %v", err)
	}
	recs := readBAMRecords(t, out.Bytes())
	if recs[0].QName != "r1" {
		t.Errorf("first record on missing tag: got %q, want %q", recs[0].QName, "r1")
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

// TestTagLessFloat exercises the float branch of tagLess (qname-secondary).
func TestTagLessFloat(t *testing.T) {
	a := &sam.Record{QName: "a", Aux: []sam.Aux{{Tag: "XF", Type: 'f', Value: 1.0}}}
	b := &sam.Record{QName: "b", Aux: []sam.Aux{{Tag: "XF", Type: 'f', Value: 2.0}}}
	if !tagLess(a, b, "XF", nil, true) {
		t.Error("XF=1.0 should be less than XF=2.0")
	}
	if tagLess(b, a, "XF", nil, true) {
		t.Error("XF=2.0 should NOT be less than XF=1.0")
	}
	// Equal floats fall back to QName.
	c := &sam.Record{QName: "z", Aux: []sam.Aux{{Tag: "XF", Type: 'f', Value: 1.0}}}
	if !tagLess(a, c, "XF", nil, true) {
		t.Error("tie on XF=1.0 should defer to QName (a < z)")
	}
}

func TestTagLessString(t *testing.T) {
	a := &sam.Record{QName: "a", Aux: []sam.Aux{{Tag: "RG", Type: 'Z', Value: "rg-aa"}}}
	b := &sam.Record{QName: "b", Aux: []sam.Aux{{Tag: "RG", Type: 'Z', Value: "rg-bb"}}}
	if !tagLess(a, b, "RG", nil, true) {
		t.Error("rg-aa < rg-bb expected")
	}
	c := &sam.Record{QName: "z", Aux: []sam.Aux{{Tag: "RG", Type: 'Z', Value: "rg-aa"}}}
	if !tagLess(a, c, "RG", nil, true) {
		t.Error("tie on rg-aa should defer to QName")
	}
}

// TestTagLessMissingOneSide mirrors upstream `bam1_cmp_by_tag`: the record
// MISSING the sort tag sorts FIRST. This used to be inverted; the change
// brings the port in line with upstream samtools sort -t TAG.
func TestTagLessMissingOneSide(t *testing.T) {
	withTag := &sam.Record{QName: "a", Aux: []sam.Aux{{Tag: "NM", Type: 'i', Value: int64(1)}}}
	without := &sam.Record{QName: "b"}
	if !tagLess(without, withTag, "NM", nil, true) {
		t.Error("missing-tag record should sort before tag-bearing record (upstream parity)")
	}
	if tagLess(withTag, without, "NM", nil, true) {
		t.Error("tag-bearing record should NOT sort before missing-tag record")
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

// largeUnsortedSAM produces enough records (with byte sizes exceeding the
// tight memory budget used by the determinism tests) to force the
// external-merge path. Records cycle across two references with positions
// chosen to scramble both the (refID, pos) and natural-QName orderings,
// guaranteeing multiple non-trivial shards.
func largeUnsortedSAM(t *testing.T, n int) string {
	t.Helper()
	var b strings.Builder
	b.WriteString("@HD\tVN:1.6\tSO:unsorted\n")
	b.WriteString("@SQ\tSN:chr1\tLN:100000\n")
	b.WriteString("@SQ\tSN:chr2\tLN:100000\n")
	for i := 0; i < n; i++ {
		// Mix references and positions deterministically but out of order.
		ref := "chr1"
		if i%3 == 0 {
			ref = "chr2"
		}
		pos := ((i*9277)%99000 + 1) // pseudo-random but reproducible
		// Pad QName with the iteration number so every record is unique
		// and the bytes payload (Seq + Qual) is non-trivial; this gives
		// the memory budget something real to count against.
		qname := "read_" + strings.Repeat("x", i%17) + strconv.Itoa(i)
		seq := strings.Repeat("ACGT", 12)
		qual := strings.Repeat("I", len(seq))
		b.WriteString(qname)
		b.WriteByte('\t')
		b.WriteString("0\t")
		b.WriteString(ref)
		b.WriteByte('\t')
		b.WriteString(strconv.Itoa(pos))
		b.WriteString("\t60\t48M\t*\t0\t0\t")
		b.WriteString(seq)
		b.WriteByte('\t')
		b.WriteString(qual)
		b.WriteByte('\n')
	}
	return b.String()
}

// TestSortParallelDeterminism asserts byte-equal BAM output across the
// serial path (Threads=0) and several parallel widths (Threads=2, 4, 8).
// This is the load-bearing parity test for the -@/--threads feature: any
// non-determinism in the worker pool (out-of-order shard collection,
// shared mutable state, racy merge tie-breaks) would surface here.
func TestSortParallelDeterminism(t *testing.T) {
	// 4 KiB budget on records of ~150-200 bytes each forces ~30-40 shards
	// at n=2000, exercising the parallel submit/drain/merge pipeline.
	const memBudget int64 = 4 * 1024
	const n = 2000
	input := largeUnsortedSAM(t, n)

	run := func(workers int) []byte {
		var out bytes.Buffer
		opts := SortOptions{
			Order:       SortCoordinate,
			OutputBAM:   true,
			MaxMemBytes: memBudget,
			TmpPrefix:   t.TempDir(),
			Threads:     workers,
		}
		if err := Sort(strings.NewReader(input), &out, opts); err != nil {
			t.Fatalf("Sort(Threads=%d): %v", workers, err)
		}
		return out.Bytes()
	}

	want := run(0) // serial baseline
	// Sanity: the baseline must actually contain records.
	if recs := readBAMRecords(t, want); len(recs) != n {
		t.Fatalf("baseline record count: got %d, want %d", len(recs), n)
	}
	for _, w := range []int{1, 2, 4, 8} {
		got := run(w)
		if !bytes.Equal(got, want) {
			t.Errorf("Sort with Threads=%d produced bytes that differ from serial baseline (len got=%d, want=%d)", w, len(got), len(want))
		}
	}
}

// TestSortParallelDeterminismByName covers the same parity property for
// the queryname sort path — its `cmp` and tie-break logic differ from
// coordinate sort, so the parallel pool's ordering needs to be verified
// separately.
func TestSortParallelDeterminismByName(t *testing.T) {
	const memBudget int64 = 4 * 1024
	const n = 1500
	input := largeUnsortedSAM(t, n)

	run := func(workers int) []byte {
		var out bytes.Buffer
		opts := SortOptions{
			Order:       SortByName,
			OutputBAM:   true,
			MaxMemBytes: memBudget,
			TmpPrefix:   t.TempDir(),
			Threads:     workers,
		}
		if err := Sort(strings.NewReader(input), &out, opts); err != nil {
			t.Fatalf("Sort(Threads=%d): %v", workers, err)
		}
		return out.Bytes()
	}
	want := run(0)
	for _, w := range []int{1, 2, 4, 8} {
		got := run(w)
		if !bytes.Equal(got, want) {
			t.Errorf("Sort byName with Threads=%d differs from serial baseline (len got=%d, want=%d)", w, len(got), len(want))
		}
	}
}

// TestSortParallelFastPath exercises the parallel-path branch where no
// shards spill (everything fits in memory). The pool is built and torn
// down without doing any work; the output must still match the serial
// path byte-for-byte.
func TestSortParallelFastPath(t *testing.T) {
	var serial bytes.Buffer
	if err := Sort(strings.NewReader(unsortedSAM), &serial, SortOptions{Order: SortCoordinate, OutputBAM: true}); err != nil {
		t.Fatalf("serial sort: %v", err)
	}
	var par bytes.Buffer
	if err := Sort(strings.NewReader(unsortedSAM), &par, SortOptions{Order: SortCoordinate, OutputBAM: true, Threads: 4}); err != nil {
		t.Fatalf("parallel sort: %v", err)
	}
	if !bytes.Equal(serial.Bytes(), par.Bytes()) {
		t.Errorf("fast-path parallel output differs from serial (lens %d vs %d)", par.Len(), serial.Len())
	}
}
