package bedmulticov

import (
	"bytes"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// upstreamFixture loads a file under testdata/parity relative to this file.
func upstreamFixture(t *testing.T, name string) []byte {
	t.Helper()
	p := filepath.Join("..", "..", "testdata", "parity", name)
	data, err := os.ReadFile(p)
	if err != nil {
		t.Fatalf("read fixture %s: %v", name, err)
	}
	return data
}

// openFixture returns the fixture as a reader (kept open until test ends).
func openFixture(t *testing.T, name string) io.Reader {
	t.Helper()
	return bytes.NewReader(upstreamFixture(t, name))
}

// Parity.t1: BED-input mirror of upstream multicov.t1 (one_block.bam vs
// multicov.bed). The single read on chr1:0-30 overlaps all four A
// intervals on chr1:15-{20,27}, regardless of strand.
func TestParity_T1_DefaultOverlap(t *testing.T) {
	want := upstreamFixture(t, "multicov.t1.expected")
	var got bytes.Buffer
	if _, err := Run(openFixture(t, "multicov.bed"),
		[]io.Reader{openFixture(t, "one_block.bed")}, &got, Options{}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	if !bytes.Equal(want, got.Bytes()) {
		t.Fatalf("output mismatch:\n got:\n%s\nwant:\n%s", got.String(), string(want))
	}
}

// Parity.t2: -s same-strand. The read is '-'; only A.{a3,a4} should match.
func TestParity_T2_SameStrand(t *testing.T) {
	want := upstreamFixture(t, "multicov.t2.expected")
	var got bytes.Buffer
	if _, err := Run(openFixture(t, "multicov.bed"),
		[]io.Reader{openFixture(t, "one_block.bed")}, &got,
		Options{SameStrand: true}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	if !bytes.Equal(want, got.Bytes()) {
		t.Fatalf("output mismatch:\n got:\n%s\nwant:\n%s", got.String(), string(want))
	}
}

// Parity.t3: -S opposite-strand. The read is '-'; only A.{a1,a2} should match.
func TestParity_T3_OppositeStrand(t *testing.T) {
	want := upstreamFixture(t, "multicov.t3.expected")
	var got bytes.Buffer
	if _, err := Run(openFixture(t, "multicov.bed"),
		[]io.Reader{openFixture(t, "one_block.bed")}, &got,
		Options{OppositeStrand: true}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	if !bytes.Equal(want, got.Bytes()) {
		t.Fatalf("output mismatch:\n got:\n%s\nwant:\n%s", got.String(), string(want))
	}
}

// Parity.t10: multi-input case. Two BED inputs against test-multi.bed:
// each file contributes its 4 records to its disjoint window.
func TestParity_T10_MultiInput(t *testing.T) {
	want := upstreamFixture(t, "multicov.t10.expected")
	var got bytes.Buffer
	if _, err := Run(openFixture(t, "test-multi.bed"),
		[]io.Reader{
			openFixture(t, "test-multi.1.bed"),
			openFixture(t, "test-multi.2.bed"),
		}, &got, Options{}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	if !bytes.Equal(want, got.Bytes()) {
		t.Fatalf("output mismatch:\n got:\n%s\nwant:\n%s", got.String(), string(want))
	}
}

// Parity.t1-BAM mirrors upstream multicov.t1 exactly: one BAM alignment
// (`one_block.bam`) on chr1:1-30 ('-' strand, MAPQ 40, CIGAR 30M) should
// overlap all four A intervals.
func TestParity_T1_BAMInput(t *testing.T) {
	want := upstreamFixture(t, "multicov.t1.expected")
	bam := makeBAM(t, []bamAln{{rname: "chr1", pos: 1, mapq: 40, cigar: "30M", flag: 16}})
	var got bytes.Buffer
	if _, err := RunSources(openFixture(t, "multicov.bed"),
		[]Source{{Reader: bytes.NewReader(bam), Kind: SourceBAM}},
		&got, Options{}); err != nil {
		t.Fatalf("RunSources: %v", err)
	}
	if !bytes.Equal(want, got.Bytes()) {
		t.Fatalf("BAM-t1 mismatch:\n got:\n%s\nwant:\n%s", got.String(), string(want))
	}
}

// Parity.t2-BAM mirrors upstream multicov.t2 (BAM + -s).
func TestParity_T2_BAMInput_SameStrand(t *testing.T) {
	want := upstreamFixture(t, "multicov.t2.expected")
	bam := makeBAM(t, []bamAln{{rname: "chr1", pos: 1, mapq: 40, cigar: "30M", flag: 16}})
	var got bytes.Buffer
	if _, err := RunSources(openFixture(t, "multicov.bed"),
		[]Source{{Reader: bytes.NewReader(bam), Kind: SourceBAM}},
		&got, Options{SameStrand: true}); err != nil {
		t.Fatalf("RunSources: %v", err)
	}
	if !bytes.Equal(want, got.Bytes()) {
		t.Fatalf("BAM-t2 mismatch:\n got:\n%s\nwant:\n%s", got.String(), string(want))
	}
}

// Parity.t3-BAM mirrors upstream multicov.t3 (BAM + -S).
func TestParity_T3_BAMInput_OppositeStrand(t *testing.T) {
	want := upstreamFixture(t, "multicov.t3.expected")
	bam := makeBAM(t, []bamAln{{rname: "chr1", pos: 1, mapq: 40, cigar: "30M", flag: 16}})
	var got bytes.Buffer
	if _, err := RunSources(openFixture(t, "multicov.bed"),
		[]Source{{Reader: bytes.NewReader(bam), Kind: SourceBAM}},
		&got, Options{OppositeStrand: true}); err != nil {
		t.Fatalf("RunSources: %v", err)
	}
	if !bytes.Equal(want, got.Bytes()) {
		t.Fatalf("BAM-t3 mismatch:\n got:\n%s\nwant:\n%s", got.String(), string(want))
	}
}

// Parity.t4-BAM mirrors upstream multicov.t4: a single `two_blocks` BAM
// alignment (CIGAR `15M10N15M`, pos 1, '-' strand). Without `-split` the
// alignment's full reference footprint [0,40) covers all four A intervals
// — identical to t1.
func TestParity_T4_BAMInput_TwoBlocks_NoSplit(t *testing.T) {
	// Upstream's exp differs from t1's only in the input BAM; the expected
	// output is byte-for-byte identical to t1.expected. We assert against
	// t1.expected for clarity.
	want := upstreamFixture(t, "multicov.t1.expected")
	bam := makeBAM(t, []bamAln{{rname: "chr1", pos: 1, mapq: 40, cigar: "15M10N15M", flag: 16}})
	var got bytes.Buffer
	if _, err := RunSources(openFixture(t, "multicov.bed"),
		[]Source{{Reader: bytes.NewReader(bam), Kind: SourceBAM}},
		&got, Options{}); err != nil {
		t.Fatalf("RunSources: %v", err)
	}
	if !bytes.Equal(want, got.Bytes()) {
		t.Fatalf("BAM-t4 mismatch:\n got:\n%s\nwant:\n%s", got.String(), string(want))
	}
}

// Parity.t5: -split alone. The `two_blocks` alignment (CIGAR 15M10N15M,
// pos 1, '-' strand) decomposes into [0,15) and [25,40). A.{a1,a3}=[15,20)
// no longer overlap any block; A.{a2,a4}=[15,27) still pick up the second
// block via [25,27).
func TestParity_T5_BAMInput_Split(t *testing.T) {
	want := upstreamFixture(t, "multicov.t5.expected")
	bam := makeBAM(t, []bamAln{{rname: "chr1", pos: 1, mapq: 40, cigar: "15M10N15M", flag: 16}})
	var got bytes.Buffer
	if _, err := RunSources(openFixture(t, "multicov.bed"),
		[]Source{{Reader: bytes.NewReader(bam), Kind: SourceBAM}},
		&got, Options{Split: true}); err != nil {
		t.Fatalf("RunSources: %v", err)
	}
	if !bytes.Equal(want, got.Bytes()) {
		t.Fatalf("BAM-t5 mismatch:\n got:\n%s\nwant:\n%s", got.String(), string(want))
	}
}

// Parity.t6: -split + -s same-strand. Alignment is '-', so only the two
// '-' A intervals are candidates; only A.a4 has a positive block overlap.
func TestParity_T6_BAMInput_Split_SameStrand(t *testing.T) {
	want := upstreamFixture(t, "multicov.t6.expected")
	bam := makeBAM(t, []bamAln{{rname: "chr1", pos: 1, mapq: 40, cigar: "15M10N15M", flag: 16}})
	var got bytes.Buffer
	if _, err := RunSources(openFixture(t, "multicov.bed"),
		[]Source{{Reader: bytes.NewReader(bam), Kind: SourceBAM}},
		&got, Options{Split: true, SameStrand: true}); err != nil {
		t.Fatalf("RunSources: %v", err)
	}
	if !bytes.Equal(want, got.Bytes()) {
		t.Fatalf("BAM-t6 mismatch:\n got:\n%s\nwant:\n%s", got.String(), string(want))
	}
}

// Parity.t7: -split + -S opposite-strand. Alignment is '-', so only the
// two '+' A intervals are candidates; only A.a2 has a positive block
// overlap (block [25,40) vs [15,27)).
func TestParity_T7_BAMInput_Split_OppositeStrand(t *testing.T) {
	want := upstreamFixture(t, "multicov.t7.expected")
	bam := makeBAM(t, []bamAln{{rname: "chr1", pos: 1, mapq: 40, cigar: "15M10N15M", flag: 16}})
	var got bytes.Buffer
	if _, err := RunSources(openFixture(t, "multicov.bed"),
		[]Source{{Reader: bytes.NewReader(bam), Kind: SourceBAM}},
		&got, Options{Split: true, OppositeStrand: true}); err != nil {
		t.Fatalf("RunSources: %v", err)
	}
	if !bytes.Equal(want, got.Bytes()) {
		t.Fatalf("BAM-t7 mismatch:\n got:\n%s\nwant:\n%s", got.String(), string(want))
	}
}

// Parity.t8: -split + -f 0.01. With overlap=2 on a2/a4 and the BAM block
// footprint = 30, upstream's check is 2/30 = 0.0667 > 0.01 → pass for
// a2/a4 (a1/a3 still don't overlap any block).
func TestParity_T8_BAMInput_Split_FractionA01(t *testing.T) {
	want := upstreamFixture(t, "multicov.t8.expected")
	bam := makeBAM(t, []bamAln{{rname: "chr1", pos: 1, mapq: 40, cigar: "15M10N15M", flag: 16}})
	var got bytes.Buffer
	if _, err := RunSources(openFixture(t, "multicov.bed"),
		[]Source{{Reader: bytes.NewReader(bam), Kind: SourceBAM}},
		&got, Options{Split: true, FractionA: 0.01}); err != nil {
		t.Fatalf("RunSources: %v", err)
	}
	if !bytes.Equal(want, got.Bytes()) {
		t.Fatalf("BAM-t8 mismatch:\n got:\n%s\nwant:\n%s", got.String(), string(want))
	}
}

// Parity.t9: -split + -f 0.10. Same arithmetic as t8 but the threshold
// is now 0.10 and 2/30 = 0.0667 is NOT > 0.10 → all counts become 0.
// This exercises the upstream-specific strict-`>` divide-by-footprint
// semantics that this port faithfully mirrors.
func TestParity_T9_BAMInput_Split_FractionA10(t *testing.T) {
	want := upstreamFixture(t, "multicov.t9.expected")
	bam := makeBAM(t, []bamAln{{rname: "chr1", pos: 1, mapq: 40, cigar: "15M10N15M", flag: 16}})
	var got bytes.Buffer
	if _, err := RunSources(openFixture(t, "multicov.bed"),
		[]Source{{Reader: bytes.NewReader(bam), Kind: SourceBAM}},
		&got, Options{Split: true, FractionA: 0.10}); err != nil {
		t.Fatalf("RunSources: %v", err)
	}
	if !bytes.Equal(want, got.Bytes()) {
		t.Fatalf("BAM-t9 mismatch:\n got:\n%s\nwant:\n%s", got.String(), string(want))
	}
}

// Parity.t10-BAM mirrors upstream multicov.t10 exactly: two BAM inputs
// with disjoint chr1 windows ([10,250) and [500,1030)) against the two A
// intervals in test-multi.bed.
func TestParity_T10_BAMInput_MultiFile(t *testing.T) {
	want := upstreamFixture(t, "multicov.t10.expected")
	// test-multi.bam: four 30M alignments at chr1 pos 10, 100, 120, 200.
	// All overlap A.1 (chr1:0-250), none overlap A.2 (chr1:500-1000).
	bam1 := makeBAM(t, []bamAln{
		{rname: "chr1", pos: 10, mapq: 1, cigar: "30M", flag: 1},
		{rname: "chr1", pos: 100, mapq: 1, cigar: "30M", flag: 1},
		{rname: "chr1", pos: 120, mapq: 1, cigar: "30M", flag: 1},
		{rname: "chr1", pos: 200, mapq: 1, cigar: "30M", flag: 1},
	})
	// test-multi.2.bam: four 30M alignments at chr1 pos 510, 520, 600, 1000.
	// None overlap A.1, all overlap A.2.
	bam2 := makeBAM(t, []bamAln{
		{rname: "chr1", pos: 510, mapq: 1, cigar: "30M", flag: 1},
		{rname: "chr1", pos: 520, mapq: 1, cigar: "30M", flag: 1},
		{rname: "chr1", pos: 600, mapq: 1, cigar: "30M", flag: 1},
		{rname: "chr1", pos: 1000, mapq: 1, cigar: "30M", flag: 1},
	})
	var got bytes.Buffer
	if _, err := RunSources(openFixture(t, "test-multi.bed"),
		[]Source{
			{Reader: bytes.NewReader(bam1), Kind: SourceBAM},
			{Reader: bytes.NewReader(bam2), Kind: SourceBAM},
		},
		&got, Options{}); err != nil {
		t.Fatalf("RunSources: %v", err)
	}
	if !bytes.Equal(want, got.Bytes()) {
		t.Fatalf("BAM-t10 mismatch:\n got:\n%s\nwant:\n%s", got.String(), string(want))
	}
}

// TestRunSources_CRAMInput exercises the CRAM path through RunSources using
// the real htslib-produced a.cram fixture (no upstream binary needed). Both
// of its mapped reads overlap chr1:10000-10100; neither overlaps chr1:9000-9500.
func TestRunSources_CRAMInput(t *testing.T) {
	cR := openFixture(t, "a.cram")
	a := "chr1\t10000\t10100\tregion1\nchr1\t9000\t9500\tregion2\n"
	var got bytes.Buffer
	if _, err := RunSources(strings.NewReader(a),
		[]Source{{Reader: cR, Kind: SourceCRAM}}, &got, Options{}); err != nil {
		t.Fatalf("RunSources(CRAM): %v", err)
	}
	want := "chr1\t10000\t10100\tregion1\t2\nchr1\t9000\t9500\tregion2\t0\n"
	if got.String() != want {
		t.Fatalf("CRAM mismatch:\n got: %q\nwant: %q", got.String(), want)
	}
}
