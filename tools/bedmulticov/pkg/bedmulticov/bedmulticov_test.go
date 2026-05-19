package bedmulticov

import (
	"bytes"
	"io"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// readers turns N input strings into N io.Readers — convenient for tests.
func readers(in ...string) []io.Reader {
	out := make([]io.Reader, len(in))
	for i, s := range in {
		out[i] = strings.NewReader(s)
	}
	return out
}

// Hand-computed: A has 2 intervals on chr1; two B files. B1 has 2 overlaps
// with A.row1 and 1 with A.row2; B2 has 0 with A.row1 and 2 with A.row2.
func TestRun_BasicTwoFiles(t *testing.T) {
	a := "chr1\t0\t100\tA1\n" +
		"chr1\t200\t300\tA2\n"
	b1 := "chr1\t10\t20\n" +
		"chr1\t50\t60\n" +
		"chr1\t250\t260\n"
	b2 := "chr1\t210\t220\n" +
		"chr1\t290\t295\n"
	var out bytes.Buffer
	n, err := Run(strings.NewReader(a), readers(b1, b2), &out, Options{})
	if err != nil {
		t.Fatalf("Run: %v", err)
	}
	if n != 2 {
		t.Fatalf("expected 2 A records, got %d", n)
	}
	want := "chr1\t0\t100\tA1\t2\t0\n" +
		"chr1\t200\t300\tA2\t1\t2\n"
	if got := out.String(); got != want {
		t.Fatalf("output mismatch:\n got:\n%s\nwant:\n%s", got, want)
	}
}

// Hand-computed: strand filter -s drops opposite-strand B records. A is
// BED6 ('+'). B1 has one + and one - overlap; the - should not count.
// With -S the result inverts.
func TestRun_StrandFilters(t *testing.T) {
	a := "chr1\t100\t200\tA1\t0\t+\n"
	b1 := "chr1\t110\t120\t.\t0\t+\n" +
		"chr1\t150\t160\t.\t0\t-\n"
	var out bytes.Buffer
	if _, err := Run(strings.NewReader(a), readers(b1), &out, Options{SameStrand: true}); err != nil {
		t.Fatalf("Run -s: %v", err)
	}
	if got, want := out.String(), "chr1\t100\t200\tA1\t0\t+\t1\n"; got != want {
		t.Fatalf("-s mismatch:\n got %q\nwant %q", got, want)
	}

	out.Reset()
	if _, err := Run(strings.NewReader(a), readers(b1), &out, Options{OppositeStrand: true}); err != nil {
		t.Fatalf("Run -S: %v", err)
	}
	if got, want := out.String(), "chr1\t100\t200\tA1\t0\t+\t1\n"; got != want {
		t.Fatalf("-S mismatch:\n got %q\nwant %q", got, want)
	}
}

// Hand-computed: -f 0.5 requires >= 50% of A covered by a single B record.
// A is 100bp; B1 has a 40bp hit (fails) and a 60bp hit (passes); B2 has
// two 20bp hits (each fails). Reciprocal: -r adds the same threshold on B.
func TestRun_FractionAndReciprocal(t *testing.T) {
	a := "chr1\t0\t100\n"
	b1 := "chr1\t0\t40\n" + // 40% of A
		"chr1\t10\t70\n" // 60% of A — passes -f 0.5
	b2 := "chr1\t0\t20\n" + // 20% of A
		"chr1\t30\t50\n" // 20% of A
	var out bytes.Buffer
	if _, err := Run(strings.NewReader(a), readers(b1, b2), &out, Options{FractionA: 0.5}); err != nil {
		t.Fatalf("Run -f: %v", err)
	}
	if got, want := out.String(), "chr1\t0\t100\t1\t0\n"; got != want {
		t.Fatalf("-f mismatch:\n got %q\nwant %q", got, want)
	}

	// -r with -f 0.5 also requires 50% of B covered: the 60bp B hit
	// extends from 10..70 and A overlap is 60bp = 100% of B, so it still
	// passes; but with -f 0.7 the 60% A coverage now fails entirely.
	out.Reset()
	if _, err := Run(strings.NewReader(a), readers(b1, b2), &out,
		Options{FractionA: 0.7, Reciprocal: true}); err != nil {
		t.Fatalf("Run -r: %v", err)
	}
	if got, want := out.String(), "chr1\t0\t100\t0\t0\n"; got != want {
		t.Fatalf("-r mismatch:\n got %q\nwant %q", got, want)
	}
}

// Conflicting strand flags should error out.
func TestRun_ConflictingStrandFlags(t *testing.T) {
	var out bytes.Buffer
	_, err := Run(strings.NewReader(""), nil, &out,
		Options{SameStrand: true, OppositeStrand: true})
	if err == nil {
		t.Fatal("expected error for -s + -S, got nil")
	}
}

// -r without -f is a user error.
func TestRun_ReciprocalWithoutFraction(t *testing.T) {
	var out bytes.Buffer
	_, err := Run(strings.NewReader(""), nil, &out, Options{Reciprocal: true})
	if err == nil {
		t.Fatal("expected error for -r without -f, got nil")
	}
}

// Range checks on -f / -F should surface user errors.
func TestRun_FractionRangeValidation(t *testing.T) {
	var out bytes.Buffer
	if _, err := Run(strings.NewReader(""), nil, &out, Options{FractionA: 1.5}); err == nil {
		t.Fatal("expected -f range error")
	}
	out.Reset()
	if _, err := Run(strings.NewReader(""), nil, &out, Options{FractionB: -0.1}); err == nil {
		t.Fatal("expected -F range error")
	}
}

// BED comments / track / browser lines and short rows should be handled
// gracefully — the former skipped, the latter surfaced as an error.
func TestRun_CommentsAndMalformed(t *testing.T) {
	a := "# header line\n" +
		"track name=foo\n" +
		"browser hide all\n" +
		"\n" +
		"chr1\t0\t10\n"
	var out bytes.Buffer
	if _, err := Run(strings.NewReader(a), nil, &out, Options{}); err != nil {
		t.Fatalf("Run on comments: %v", err)
	}
	if got, want := out.String(), "chr1\t0\t10\n"; got != want {
		t.Fatalf("comment handling mismatch:\n got %q\nwant %q", got, want)
	}

	// Malformed: 2 columns.
	out.Reset()
	if _, err := Run(strings.NewReader("chr1\t0\n"), nil, &out, Options{}); err == nil {
		t.Fatal("expected error on short record")
	}
	// Malformed: non-int start.
	out.Reset()
	if _, err := Run(strings.NewReader("chr1\tBAD\t10\n"), nil, &out, Options{}); err == nil {
		t.Fatal("expected error on non-int start")
	}
	// Malformed: non-int end.
	out.Reset()
	if _, err := Run(strings.NewReader("chr1\t0\tBAD\n"), nil, &out, Options{}); err == nil {
		t.Fatal("expected error on non-int end")
	}
}

// strand "" on either side under -s / -S should be treated as "no match".
func TestRun_MissingStrandUnderFilter(t *testing.T) {
	a := "chr1\t0\t100\tA1\t0\t+\n"
	b := "chr1\t10\t20\n" // no strand column
	var out bytes.Buffer
	if _, err := Run(strings.NewReader(a), readers(b), &out, Options{SameStrand: true}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	if got, want := out.String(), "chr1\t0\t100\tA1\t0\t+\t0\n"; got != want {
		t.Fatalf("missing-strand mismatch:\n got %q\nwant %q", got, want)
	}
	out.Reset()
	if _, err := Run(strings.NewReader(a), readers(b), &out, Options{OppositeStrand: true}); err != nil {
		t.Fatalf("Run -S: %v", err)
	}
	if got, want := out.String(), "chr1\t0\t100\tA1\t0\t+\t0\n"; got != want {
		t.Fatalf("missing-strand -S mismatch:\n got %q\nwant %q", got, want)
	}
}

// bamAln is a compact alignment description used to build test BAM streams.
type bamAln struct {
	rname string
	pos   int32 // 1-based per SAM
	mapq  uint8
	cigar string
	flag  uint16
}

// makeBAM builds an in-memory BGZF-wrapped BAM containing the given
// alignments. Header @SQ lines are inferred from the alignments' RNAMEs.
func makeBAM(t *testing.T, alns []bamAln) []byte {
	t.Helper()
	hdr := &sam.Header{}
	seen := map[string]bool{}
	for _, a := range alns {
		if seen[a.rname] {
			continue
		}
		seen[a.rname] = true
		hdr.Refs = append(hdr.Refs, sam.Reference{Name: a.rname, Length: 1000000})
	}
	var buf bytes.Buffer
	bw := sam.NewBAMWriter(&buf)
	if err := bw.WriteHeader(hdr); err != nil {
		t.Fatalf("WriteHeader: %v", err)
	}
	for _, a := range alns {
		cig, err := sam.ParseCigar(a.cigar)
		if err != nil {
			t.Fatalf("ParseCigar(%q): %v", a.cigar, err)
		}
		// Seq length must match CIGAR query length, but the BAM writer
		// accepts "*" (empty Seq). We use a literal of the right length so
		// the round-trip is well-formed.
		ql := cig.QueryLength()
		seq := strings.Repeat("A", ql)
		rec := &sam.Record{
			QName: "r",
			Flag:  a.flag,
			RName: a.rname,
			Pos:   a.pos,
			MapQ:  a.mapq,
			Cigar: cig,
			RNext: "*",
			Seq:   seq,
		}
		if err := bw.Write(rec); err != nil {
			t.Fatalf("Write: %v", err)
		}
	}
	if err := bw.Close(); err != nil {
		t.Fatalf("Close BAM: %v", err)
	}
	return buf.Bytes()
}

// TestRun_BAMInput exercises the BAM path end-to-end through RunSources.
// One alignment on chr1:1-30 must overlap all four A intervals (regardless
// of strand). Mirrors upstream's multicov.t1.
func TestRun_BAMInput(t *testing.T) {
	bam := makeBAM(t, []bamAln{{rname: "chr1", pos: 1, mapq: 40, cigar: "30M", flag: 16}})
	a := "chr1\t15\t20\ta1\t1\t+\n" +
		"chr1\t15\t27\ta2\t2\t+\n" +
		"chr1\t15\t20\ta3\t3\t-\n" +
		"chr1\t15\t27\ta4\t4\t-\n"
	var got bytes.Buffer
	if _, err := RunSources(strings.NewReader(a),
		[]Source{{Reader: bytes.NewReader(bam), Kind: SourceBAM}},
		&got, Options{}); err != nil {
		t.Fatalf("RunSources: %v", err)
	}
	want := "chr1\t15\t20\ta1\t1\t+\t1\n" +
		"chr1\t15\t27\ta2\t2\t+\t1\n" +
		"chr1\t15\t20\ta3\t3\t-\t1\n" +
		"chr1\t15\t27\ta4\t4\t-\t1\n"
	if got.String() != want {
		t.Fatalf("BAM-input mismatch:\n got:\n%s\nwant:\n%s", got.String(), want)
	}
}

// TestRun_BAMInput_SameStrand: under -s, the '-' alignment only matches A
// intervals on '-' strand (a3, a4). Mirrors upstream's multicov.t2.
func TestRun_BAMInput_SameStrand(t *testing.T) {
	bam := makeBAM(t, []bamAln{{rname: "chr1", pos: 1, mapq: 40, cigar: "30M", flag: 16}})
	a := "chr1\t15\t20\ta1\t1\t+\n" +
		"chr1\t15\t27\ta2\t2\t+\n" +
		"chr1\t15\t20\ta3\t3\t-\n" +
		"chr1\t15\t27\ta4\t4\t-\n"
	var got bytes.Buffer
	if _, err := RunSources(strings.NewReader(a),
		[]Source{{Reader: bytes.NewReader(bam), Kind: SourceBAM}},
		&got, Options{SameStrand: true}); err != nil {
		t.Fatalf("RunSources: %v", err)
	}
	want := "chr1\t15\t20\ta1\t1\t+\t0\n" +
		"chr1\t15\t27\ta2\t2\t+\t0\n" +
		"chr1\t15\t20\ta3\t3\t-\t1\n" +
		"chr1\t15\t27\ta4\t4\t-\t1\n"
	if got.String() != want {
		t.Fatalf("BAM-input -s mismatch:\n got:\n%s\nwant:\n%s", got.String(), want)
	}
}

// TestRun_BAMInput_OppositeStrand: under -S, the '-' alignment only matches
// '+'-strand A intervals (a1, a2). Mirrors upstream's multicov.t3.
func TestRun_BAMInput_OppositeStrand(t *testing.T) {
	bam := makeBAM(t, []bamAln{{rname: "chr1", pos: 1, mapq: 40, cigar: "30M", flag: 16}})
	a := "chr1\t15\t20\ta1\t1\t+\n" +
		"chr1\t15\t27\ta2\t2\t+\n" +
		"chr1\t15\t20\ta3\t3\t-\n" +
		"chr1\t15\t27\ta4\t4\t-\n"
	var got bytes.Buffer
	if _, err := RunSources(strings.NewReader(a),
		[]Source{{Reader: bytes.NewReader(bam), Kind: SourceBAM}},
		&got, Options{OppositeStrand: true}); err != nil {
		t.Fatalf("RunSources: %v", err)
	}
	want := "chr1\t15\t20\ta1\t1\t+\t1\n" +
		"chr1\t15\t27\ta2\t2\t+\t1\n" +
		"chr1\t15\t20\ta3\t3\t-\t0\n" +
		"chr1\t15\t27\ta4\t4\t-\t0\n"
	if got.String() != want {
		t.Fatalf("BAM-input -S mismatch:\n got:\n%s\nwant:\n%s", got.String(), want)
	}
}

// TestRun_BAMInput_MAPQFilter: alignments below the -q threshold are
// dropped during indexing.
func TestRun_BAMInput_MAPQFilter(t *testing.T) {
	bam := makeBAM(t, []bamAln{
		{rname: "chr1", pos: 1, mapq: 5, cigar: "30M", flag: 0},
		{rname: "chr1", pos: 1, mapq: 40, cigar: "30M", flag: 0},
	})
	a := "chr1\t15\t20\n"
	var got bytes.Buffer
	if _, err := RunSources(strings.NewReader(a),
		[]Source{{Reader: bytes.NewReader(bam), Kind: SourceBAM}},
		&got, Options{MinMAPQ: 20}); err != nil {
		t.Fatalf("RunSources: %v", err)
	}
	if got.String() != "chr1\t15\t20\t1\n" {
		t.Fatalf("MAPQ filter mismatch: got %q", got.String())
	}
}

// TestRun_BAMInput_MaxDepthCap: -D caps the count reported per A interval
// per BAM input.
func TestRun_BAMInput_MaxDepthCap(t *testing.T) {
	// Stage 5 identical alignments overlapping the single A interval.
	alns := make([]bamAln, 5)
	for i := range alns {
		alns[i] = bamAln{rname: "chr1", pos: 1, mapq: 40, cigar: "30M", flag: 0}
	}
	bam := makeBAM(t, alns)
	a := "chr1\t15\t20\n"
	var got bytes.Buffer
	if _, err := RunSources(strings.NewReader(a),
		[]Source{{Reader: bytes.NewReader(bam), Kind: SourceBAM}},
		&got, Options{MaxDepth: 3}); err != nil {
		t.Fatalf("RunSources: %v", err)
	}
	if got.String() != "chr1\t15\t20\t3\n" {
		t.Fatalf("-D cap mismatch: got %q (want capped at 3)", got.String())
	}
	// MaxDepth=0 disables the cap.
	got.Reset()
	if _, err := RunSources(strings.NewReader(a),
		[]Source{{Reader: bytes.NewReader(bam), Kind: SourceBAM}},
		&got, Options{MaxDepth: 0}); err != nil {
		t.Fatalf("RunSources: %v", err)
	}
	if got.String() != "chr1\t15\t20\t5\n" {
		t.Fatalf("-D=0 disable mismatch: got %q (want 5)", got.String())
	}
}

// TestRun_BAMInput_SkipsUnmappedSecondaryDup: records flagged unmapped /
// secondary / supplementary / duplicate / QC-fail must not be indexed.
func TestRun_BAMInput_SkipsUnmappedSecondaryDup(t *testing.T) {
	alns := []bamAln{
		{rname: "chr1", pos: 1, mapq: 40, cigar: "30M", flag: 4},    // unmapped
		{rname: "chr1", pos: 1, mapq: 40, cigar: "30M", flag: 256},  // secondary
		{rname: "chr1", pos: 1, mapq: 40, cigar: "30M", flag: 2048}, // supplementary
		{rname: "chr1", pos: 1, mapq: 40, cigar: "30M", flag: 1024}, // duplicate
		{rname: "chr1", pos: 1, mapq: 40, cigar: "30M", flag: 512},  // QC fail
		{rname: "chr1", pos: 1, mapq: 40, cigar: "30M", flag: 0},    // primary, counted
	}
	bam := makeBAM(t, alns)
	a := "chr1\t15\t20\n"
	var got bytes.Buffer
	if _, err := RunSources(strings.NewReader(a),
		[]Source{{Reader: bytes.NewReader(bam), Kind: SourceBAM}},
		&got, Options{}); err != nil {
		t.Fatalf("RunSources: %v", err)
	}
	if got.String() != "chr1\t15\t20\t1\n" {
		t.Fatalf("flag filter mismatch: got %q (want exactly 1)", got.String())
	}
}

// TestRunSources_RejectsUnknownKind: defensive — an out-of-range SourceKind
// surfaces a clear error.
func TestRunSources_RejectsUnknownKind(t *testing.T) {
	var got bytes.Buffer
	_, err := RunSources(strings.NewReader(""),
		[]Source{{Reader: strings.NewReader(""), Kind: SourceKind(99)}},
		&got, Options{})
	if err == nil {
		t.Fatal("expected error for unknown SourceKind")
	}
}

// TestRunSources_RejectsBadMAPQ / MaxDepth: option validation.
func TestRunSources_OptionValidation(t *testing.T) {
	var got bytes.Buffer
	if _, err := RunSources(strings.NewReader(""), nil, &got, Options{MinMAPQ: -1}); err == nil {
		t.Fatal("expected error for negative -q")
	}
	if _, err := RunSources(strings.NewReader(""), nil, &got, Options{MinMAPQ: 300}); err == nil {
		t.Fatal("expected error for out-of-range -q")
	}
	if _, err := RunSources(strings.NewReader(""), nil, &got, Options{MaxDepth: -1}); err == nil {
		t.Fatal("expected error for negative -D")
	}
}

// TestSplit_Reciprocal exercises the -split + -r path: with Reciprocal
// set, the upstream check additionally requires
// total_overlap/lenA > FractionA. A short A interval that satisfies the
// per-footprint check but fails the per-A check should NOT count.
func TestSplit_Reciprocal(t *testing.T) {
	a := "chr1\t0\t200\n" // long A: 1/200 < 0.1 fails reciprocal at A side.
	// One read with CIGAR 1M99N1M starting at pos 1: blocks [0,1) and
	// [100,101). Footprint=2. Total overlap with A=[0,200): 2.
	// 2/2 = 1.0 > 0.1 (footprint side, pass).
	// 2/200 = 0.01 < 0.1 (A side, fail under -r).
	bam := makeBAM(t, []bamAln{{rname: "chr1", pos: 1, mapq: 0, cigar: "1M99N1M", flag: 0}})
	var got bytes.Buffer
	if _, err := RunSources(strings.NewReader(a),
		[]Source{{Reader: bytes.NewReader(bam), Kind: SourceBAM}},
		&got, Options{Split: true, FractionA: 0.1, Reciprocal: true}); err != nil {
		t.Fatalf("RunSources: %v", err)
	}
	if got.String() != "chr1\t0\t200\t0\n" {
		t.Fatalf("expected 0 under -r; got %q", got.String())
	}

	// Now drop Reciprocal: the footprint check (1.0 > 0.1) passes, so
	// the alignment counts once.
	var got2 bytes.Buffer
	if _, err := RunSources(strings.NewReader(a),
		[]Source{{Reader: bytes.NewReader(bam), Kind: SourceBAM}},
		&got2, Options{Split: true, FractionA: 0.1}); err != nil {
		t.Fatalf("RunSources (no -r): %v", err)
	}
	if got2.String() != "chr1\t0\t200\t1\n" {
		t.Fatalf("expected 1 without -r; got %q", got2.String())
	}
}

// TestSplit_OncePerAlignment: two blocks of one alignment both overlap
// the same A interval — the alignment must count ONCE (the whole point
// of `-split`'s de-duplication contract).
func TestSplit_OncePerAlignment(t *testing.T) {
	a := "chr1\t0\t100\n" // covers both blocks
	bam := makeBAM(t, []bamAln{{rname: "chr1", pos: 1, mapq: 0, cigar: "10M20N10M", flag: 0}})
	var got bytes.Buffer
	if _, err := RunSources(strings.NewReader(a),
		[]Source{{Reader: bytes.NewReader(bam), Kind: SourceBAM}},
		&got, Options{Split: true}); err != nil {
		t.Fatalf("RunSources: %v", err)
	}
	if got.String() != "chr1\t0\t100\t1\n" {
		t.Fatalf("expected 1 (once-per-aln dedupe); got %q", got.String())
	}
}

// TestSplit_DeletionExtendsBlock: D ops extend the current block under
// multicov's `-split` semantics (breakOnDeletionOps=false in upstream).
func TestSplit_DeletionExtendsBlock(t *testing.T) {
	a := "chr1\t0\t100\n"
	// 5M2D5M produces ONE block of length 12 starting at pos 0 (no N).
	bam := makeBAM(t, []bamAln{{rname: "chr1", pos: 1, mapq: 0, cigar: "5M2D5M", flag: 0}})
	var got bytes.Buffer
	if _, err := RunSources(strings.NewReader(a),
		[]Source{{Reader: bytes.NewReader(bam), Kind: SourceBAM}},
		&got, Options{Split: true}); err != nil {
		t.Fatalf("RunSources: %v", err)
	}
	if got.String() != "chr1\t0\t100\t1\n" {
		t.Fatalf("expected D to extend block; got %q", got.String())
	}
}

// TestSplit_AllSkippedAlignments: alignments that consume zero
// reference bases produce no blocks and contribute no counts.
func TestSplit_AllSkippedAlignments(t *testing.T) {
	a := "chr1\t0\t100\n"
	// Pure soft-clip — no reference advance.
	bam := makeBAM(t, []bamAln{{rname: "chr1", pos: 1, mapq: 0, cigar: "30S", flag: 0}})
	var got bytes.Buffer
	if _, err := RunSources(strings.NewReader(a),
		[]Source{{Reader: bytes.NewReader(bam), Kind: SourceBAM}},
		&got, Options{Split: true}); err != nil {
		t.Fatalf("RunSources: %v", err)
	}
	if got.String() != "chr1\t0\t100\t0\n" {
		t.Fatalf("expected 0 for zero-block alignment; got %q", got.String())
	}
}

// Multi-chrom inputs and intervals that span chrom gaps should still
// report a 0 for the file that has no records on a given chrom.
func TestRun_MultiChromAndMissingChrom(t *testing.T) {
	a := "chr1\t0\t100\nchr2\t0\t100\n"
	b1 := "chr1\t10\t20\n" // chr1 only
	b2 := "chr2\t10\t20\n" // chr2 only
	var out bytes.Buffer
	if _, err := Run(strings.NewReader(a), readers(b1, b2), &out, Options{}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	want := "chr1\t0\t100\t1\t0\n" +
		"chr2\t0\t100\t0\t1\n"
	if got := out.String(); got != want {
		t.Fatalf("output mismatch:\n got:\n%s\nwant:\n%s", got, want)
	}
}
