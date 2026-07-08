package mosdepth

import (
	"bytes"
	"compress/gzip"
	"io"
	"os"
	"path/filepath"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// TestSortedGuardInvalidation proves the covAccum.sorted fast-path guard is
// correctly invalidated by every event-appending path. sortEvents must be a
// no-op only while no new event has been appended since the last sort; any
// add / overlap-correction append must force the next sortEvents to re-sort.
func TestSortedGuardInvalidation(t *testing.T) {
	a := newCovAccum(100)
	// Insert two runs out of order so the slice genuinely needs sorting.
	a.add(50, 60)
	a.add(10, 20)
	if a.sorted {
		t.Fatalf("sorted should be false after add() before any sortEvents")
	}
	a.sortEvents()
	if !a.sorted {
		t.Fatalf("sorted should be true after sortEvents")
	}
	if !isSortedByPos(a.events) {
		t.Fatalf("events not sorted after sortEvents: %+v", a.events)
	}
	// A second sortEvents while sorted==true must not disturb the ordering and
	// must be a genuine no-op (the fast path). Capture the slice header to make
	// sure nothing is reallocated.
	before := append([]covEvent(nil), a.events...)
	a.sortEvents()
	if !equalEvents(before, a.events) {
		t.Fatalf("no-op sortEvents changed events")
	}
	// Appending a new low-position event must clear sorted and the next
	// sortEvents must actually re-sort it into place.
	a.add(0, 5)
	if a.sorted {
		t.Fatalf("add() after a sort must clear sorted")
	}
	a.sortEvents()
	if !isSortedByPos(a.events) {
		t.Fatalf("events not re-sorted after post-sort add: %+v", a.events)
	}
}

// TestRegionCursorResetOnSort proves the resumable region cursor is invalidated
// whenever the event slice is re-sorted, so a stale cursor can never seed a
// regionStats call with an index into a differently-ordered slice.
func TestRegionCursorResetOnSort(t *testing.T) {
	a := newCovAccum(100)
	a.add(10, 90)
	// Seed the cursor via a regionStats call (which sorts and checkpoints).
	a.regionStats(20, 30, nil, nil, 0)
	if !a.curValid {
		t.Fatalf("cursor should be valid after regionStats")
	}
	// Appending an event and re-sorting must reset the cursor.
	a.add(0, 100)
	a.sortEvents()
	if a.curValid {
		t.Fatalf("cursor should be reset after a re-sort")
	}
}

// TestStreamingMatchesBuffered is the end-to-end determinism guard for the
// streaming (SO:coordinate) driver: for the same records it must produce output
// byte-for-byte identical to the buffered fallback path. We run the pipeline
// twice over a multi-chromosome record set — once with an @HD SO:coordinate
// header line (streaming) and once without it (buffered) — and require every
// output file to match exactly.
func TestStreamingMatchesBuffered(t *testing.T) {
	refs := []sam.Reference{
		{Name: "chr1", Length: 60},
		{Name: "chr2", Length: 40},
		{Name: "chr3", Length: 30}, // zero-coverage reference (no reads).
		{Name: "chr4", Length: 50},
	}
	mk := func(name, ref string, pos int32, cigar string, flag uint16) *sam.Record {
		c := mustParseCigar(t, cigar)
		return &sam.Record{
			QName: name, RName: ref, Pos: int64(pos), MapQ: 60, Flag: flag,
			Cigar: c,
			Seq:   repeatA(c.QueryLength()),
		}
	}
	// Coordinate order: chr1 reads, then chr2, then (skipping zero-coverage
	// chr3) chr4. Includes a proper-pair overlap on chr1 to exercise the
	// overlap-correction append path under streaming.
	recs := []*sam.Record{
		mkPair(t, "p1", "chr1", 5, 12, "10M"),  // overlapping mate pair
		mkPair2(t, "p1", "chr1", 12, 5, "10M"), // mate of p1
		mk("a", "chr1", 30, "5M", 0),
		mk("b", "chr2", 3, "8M", 0),
		mk("c", "chr2", 20, "6M", 0),
		mk("d", "chr4", 10, "15M", 0),
	}

	stream := runToDir(t, refs, recs, true)
	buffer := runToDir(t, refs, recs, false)

	names := []string{
		"x.mosdepth.global.dist.txt",
		"x.mosdepth.summary.txt",
		"x.per-base.bed.gz",
	}
	for _, n := range names {
		sb := mustRead(t, filepath.Join(stream, n))
		bb := mustRead(t, filepath.Join(buffer, n))
		if n[len(n)-3:] == ".gz" {
			sb = mustGunzip(t, sb)
			bb = mustGunzip(t, bb)
		}
		if !bytes.Equal(sb, bb) {
			t.Fatalf("streaming vs buffered differ for %s:\nstream=%q\nbuffer=%q", n, sb, bb)
		}
	}
}

// TestStreamingMatchesBufferedByRegions is the --by counterpart to
// TestStreamingMatchesBuffered: it proves the streaming (SO:coordinate) and
// buffered drivers produce byte-identical region outputs when --by selects a
// mix of BED intervals AND a fixed window. This exercises the resumable region
// cursor and the sorted-guard under the exact multi-region-per-chromosome
// workload the perf fix targets. The regions.bed.gz (decompressed),
// region.dist.txt, and summary outputs must match exactly between the two
// paths.
func TestStreamingMatchesBufferedByRegions(t *testing.T) {
	refs := []sam.Reference{
		{Name: "chr1", Length: 60},
		{Name: "chr2", Length: 40},
		{Name: "chr3", Length: 30}, // zero-coverage reference (no reads).
		{Name: "chr4", Length: 50},
	}
	mk := func(name, ref string, pos int32, cigar string) *sam.Record {
		c := mustParseCigar(t, cigar)
		return &sam.Record{
			QName: name, RName: ref, Pos: int64(pos), MapQ: 60, Flag: 0,
			Cigar: c,
			Seq:   repeatA(c.QueryLength()),
		}
	}
	recs := []*sam.Record{
		mkPair(t, "p1", "chr1", 5, 12, "10M"),  // overlapping mate pair
		mkPair2(t, "p1", "chr1", 12, 5, "10M"), // mate of p1
		mk("a", "chr1", 30, "5M"),
		mk("b", "chr2", 3, "8M"),
		mk("c", "chr2", 20, "6M"),
		mk("d", "chr4", 10, "15M"),
	}

	// A BED with multiple, deliberately overlapping regions per chromosome so
	// the region cursor and sorted-guard are exercised heavily on both paths.
	bed := "" +
		"chr1\t0\t20\tr1a\n" +
		"chr1\t10\t40\tr1b\n" +
		"chr1\t5\t15\tr1c\n" +
		"chr2\t0\t10\tr2a\n" +
		"chr2\t15\t30\tr2b\n" +
		"chr4\t5\t25\tr4a\n"

	runBy := func(coordSorted bool) string {
		t.Helper()
		dir := t.TempDir()
		bedPath := filepath.Join(dir, "regions.bed")
		if err := os.WriteFile(bedPath, []byte(bed), 0o644); err != nil {
			t.Fatalf("write bed: %v", err)
		}
		bam := makeBAMSorted(t, refs, recs, coordSorted)
		opts := Options{
			Prefix:      filepath.Join(dir, "x"),
			ExcludeFlag: DefaultExcludeFlag,
			ByBED:       bedPath,
			Thresholds:  []int{1, 5},
		}
		if err := Run(bytes.NewReader(bam), opts); err != nil {
			t.Fatalf("Run(--by, coordSorted=%v): %v", coordSorted, err)
		}
		return dir
	}

	stream := runBy(true)
	buffer := runBy(false)

	names := []string{
		"x.mosdepth.global.dist.txt",
		"x.mosdepth.region.dist.txt",
		"x.mosdepth.summary.txt",
		"x.regions.bed.gz",
		"x.thresholds.bed.gz",
	}
	for _, n := range names {
		sb := mustRead(t, filepath.Join(stream, n))
		bb := mustRead(t, filepath.Join(buffer, n))
		if n[len(n)-3:] == ".gz" {
			sb = mustGunzip(t, sb)
			bb = mustGunzip(t, bb)
		}
		if !bytes.Equal(sb, bb) {
			t.Fatalf("streaming vs buffered --by differ for %s:\nstream=%q\nbuffer=%q", n, sb, bb)
		}
	}
}

// TestRegionCursorOutOfOrderMatchesFreshAccum locks the resumable region
// cursor's rewind-on-earlier-region behaviour: for a batch of deliberately
// overlapping and descending-begin (out-of-order) regions sharing one
// accumulator, each region's coverage sum must equal what a fresh-from-zero
// accumulator computes independently (a brute-force oracle). If the cursor ever
// failed to rewind when a region begins before curPos, the shared-accumulator
// sums would diverge from the fresh ones.
func TestRegionCursorOutOfOrderMatchesFreshAccum(t *testing.T) {
	const refLen = 200
	// A spread of overlapping coverage runs so depth varies across the reference.
	runs := [][2]int{
		{0, 50}, {10, 40}, {20, 120}, {30, 35}, {60, 180}, {90, 100}, {150, 200},
	}
	build := func() *covAccum {
		a := newCovAccum(refLen)
		for _, r := range runs {
			a.add(r[0], r[1])
		}
		return a
	}
	// Deliberately overlapping and out-of-order (descending, then jumping back
	// forward) so the cursor must both resume and rewind.
	regions := [][2]int{
		{100, 160}, // forward
		{20, 80},   // jumps backwards -> cursor must rewind
		{50, 130},  // overlaps previous, begins after
		{0, 30},    // jumps back again
		{140, 200}, // forward
		{10, 190},  // wide, begins before curPos
		{10, 190},  // exact repeat
	}

	shared := build()
	for _, reg := range regions {
		gotSum, _, gotMin, gotMax, _ := shared.regionStats(reg[0], reg[1], nil, nil, 0)
		// Fresh, independent oracle: a brand-new accumulator scanned once.
		fresh := build()
		wantSum, _, wantMin, wantMax, _ := fresh.regionStats(reg[0], reg[1], nil, nil, 0)
		if gotSum != wantSum {
			t.Fatalf("region %v: cursor sum=%d, fresh oracle=%d", reg, gotSum, wantSum)
		}
		if gotMin != wantMin || gotMax != wantMax {
			t.Fatalf("region %v: cursor min/max=%d/%d, fresh oracle=%d/%d",
				reg, gotMin, gotMax, wantMin, wantMax)
		}
	}
}

// runToDir writes a BAM from refs+recs (optionally SO:coordinate to select the
// streaming driver) and runs the default mosdepth pipeline into a fresh temp
// dir, returning that dir. The per-base output is enabled so the coverage
// track is compared too.
func runToDir(t *testing.T, refs []sam.Reference, recs []*sam.Record, coordSorted bool) string {
	t.Helper()
	dir := t.TempDir()
	bam := makeBAMSorted(t, refs, recs, coordSorted)
	opts := Options{
		Prefix:      filepath.Join(dir, "x"),
		ExcludeFlag: DefaultExcludeFlag,
	}
	if err := Run(bytes.NewReader(bam), opts); err != nil {
		t.Fatalf("Run(coordSorted=%v): %v", coordSorted, err)
	}
	return dir
}

// makeBAMSorted is makeBAM with an optional @HD SO:coordinate line so the test
// can drive either the streaming or the buffered path from the same records.
func makeBAMSorted(t *testing.T, refs []sam.Reference, recs []*sam.Record, coordSorted bool) []byte {
	t.Helper()
	hdr := &sam.Header{Refs: refs}
	if coordSorted {
		hd := sam.HeaderLine{Tag: "HD", Fields: []sam.HeaderField{
			{Tag: "VN", Value: "1.6"},
			{Tag: "SO", Value: "coordinate"},
		}}
		hdr.Lines = append(hdr.Lines, hd)
		hdr.HDFields = hd.Fields
	}
	for _, r := range refs {
		hdr.Lines = append(hdr.Lines, sam.HeaderLine{
			Tag: "SQ",
			Fields: []sam.HeaderField{
				{Tag: "SN", Value: r.Name},
				{Tag: "LN", Value: intToStr(int(r.Length))},
			},
		})
	}
	var buf bytes.Buffer
	bw := sam.NewBAMWriter(&buf)
	if err := bw.WriteHeader(hdr); err != nil {
		t.Fatalf("WriteHeader: %v", err)
	}
	for _, rec := range recs {
		if err := bw.Write(rec); err != nil {
			t.Fatalf("Write rec: %v", err)
		}
	}
	if err := bw.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	return buf.Bytes()
}

// mkPair / mkPair2 build the two mates of a proper, overlapping read pair so the
// streaming path exercises the addOverlapCorrection event append (which must
// also clear the sorted guard).
func mkPair(t *testing.T, name, ref string, pos, matePos int32, cigar string) *sam.Record {
	t.Helper()
	c := mustParseCigar(t, cigar)
	return &sam.Record{
		QName: name, RName: ref, Pos: int64(pos), MapQ: 60,
		Flag:  sam.FlagPaired | sam.FlagProperPair | sam.FlagRead1,
		Cigar: c, RNext: "=", PNext: int64(matePos), TLen: 17,
		Seq: repeatA(c.QueryLength()),
	}
}

func mkPair2(t *testing.T, name, ref string, pos, matePos int32, cigar string) *sam.Record {
	t.Helper()
	c := mustParseCigar(t, cigar)
	return &sam.Record{
		QName: name, RName: ref, Pos: int64(pos), MapQ: 60,
		Flag:  sam.FlagPaired | sam.FlagProperPair | sam.FlagRead2,
		Cigar: c, RNext: "=", PNext: int64(matePos), TLen: -17,
		Seq: repeatA(c.QueryLength()),
	}
}

func repeatA(n int) string {
	b := make([]byte, n)
	for i := range b {
		b[i] = 'A'
	}
	return string(b)
}

func isSortedByPos(ev []covEvent) bool {
	for i := 1; i < len(ev); i++ {
		if ev[i-1].pos > ev[i].pos {
			return false
		}
	}
	return true
}

func equalEvents(a, b []covEvent) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

func mustRead(t *testing.T, path string) []byte {
	t.Helper()
	b, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("read %s: %v", path, err)
	}
	return b
}

func mustGunzip(t *testing.T, b []byte) []byte {
	t.Helper()
	gr, err := gzip.NewReader(bytes.NewReader(b))
	if err != nil {
		t.Fatalf("gzip.NewReader: %v", err)
	}
	gr.Multistream(true)
	out, err := io.ReadAll(gr)
	if err != nil {
		t.Fatalf("gunzip: %v", err)
	}
	return out
}
