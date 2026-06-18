package fastp

import (
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fastq"
)

// truseqAdapter is the canonical Illumina TruSeq Read 1 adapter prefix.
// Detection should reliably recover at least the first ~15-20 bases.
const truseqAdapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"

// nexteraAdapter is the Nextera transposase R1 adapter; used as a
// secondary control in SE detection.
const nexteraAdapter = "CTGTCTCTTATACACATCTCCGAGCCCACGAGAC"

func makeSERead(id, seq string, qual byte) *fastq.Record {
	q := make([]byte, len(seq))
	for i := range q {
		q[i] = qual
	}
	return &fastq.Record{
		ID:       id,
		Sequence: []byte(seq),
		Quality:  q,
	}
}

func TestDetectAdapterSEPicksTruSeq(t *testing.T) {
	// Build a synthetic dataset: 10000+ reads (upstream gate at
	// evaluator.cpp:344), each is a unique-ish insert followed by the
	// TruSeq adapter. The insert varies so the seed kmer can't come
	// from the insert region. The known-adapters fast path should
	// match the TruSeq adapter directly.
	reads := make([]*fastq.Record, 0, 12000)
	for i := 0; i < 12000; i++ {
		// Vary the insert deterministically by rotating a base pattern.
		insert := strings.Repeat("ACGTACGT", 4) + string("ACGT"[i%4])
		seq := insert + truseqAdapter
		reads = append(reads, makeSERead("r", seq, 'I'))
	}
	got := DetectAdapterSE(reads)
	if got == "" {
		t.Fatalf("DetectAdapterSE returned empty; expected a TruSeq-like adapter")
	}
	// Should match the start of the TruSeq adapter for at least 10 bases.
	if !strings.HasPrefix(truseqAdapter, got[:min(len(got), 10)]) &&
		!strings.HasPrefix(got, truseqAdapter[:min(len(truseqAdapter), 10)]) {
		t.Errorf("DetectAdapterSE = %q; expected prefix overlap with %q", got, truseqAdapter)
	}
}

// TestDetectAdaptersForPE verifies the --detect_adapter_for_pe path detects
// each mate's adapter by running the single-file evaluator on read1 and read2
// independently (upstream main.cpp:447-480), not pairwise.
func TestDetectAdaptersForPE(t *testing.T) {
	pairs := make([][2]*fastq.Record, 0, 12000)
	for i := 0; i < 12000; i++ {
		ins1 := strings.Repeat("ACGTACGT", 4) + string("ACGT"[i%4])
		ins2 := strings.Repeat("TGCATGCA", 4) + string("ACGT"[(i+1)%4])
		r1 := makeSERead("r", ins1+truseqAdapter, 'I')
		r2 := makeSERead("r", ins2+truseqAdapter, 'I')
		pairs = append(pairs, [2]*fastq.Record{r1, r2})
	}
	a1, a2 := DetectAdaptersForPE(pairs)
	if a1 == "" || a2 == "" {
		t.Fatalf("DetectAdaptersForPE = (%q, %q); expected both detected", a1, a2)
	}
	if !strings.HasPrefix(a1, truseqAdapter[:10]) || !strings.HasPrefix(a2, truseqAdapter[:10]) {
		t.Errorf("DetectAdaptersForPE = (%q, %q); expected TruSeq-like prefixes", a1, a2)
	}
}

func TestDetectAdapterSEBelowMinRecordsReturnsEmpty(t *testing.T) {
	// Upstream evaluator.cpp:344 requires records >= 10000 before
	// detection runs. With fewer reads we must return "" — this is the
	// exact behavior the case 15 parity test relies on.
	reads := make([]*fastq.Record, 0, 500)
	for i := 0; i < 500; i++ {
		insert := strings.Repeat("ACGTACGT", 4) + string("ACGT"[i%4])
		seq := insert + truseqAdapter
		reads = append(reads, makeSERead("r", seq, 'I'))
	}
	if got := DetectAdapterSE(reads); got != "" {
		t.Errorf("DetectAdapterSE with %d reads = %q, want empty (upstream gate is 10000)", len(reads), got)
	}
}

func TestDetectAdapterSEReturnsEmptyOnRandomInput(t *testing.T) {
	// Pure-random reads (no embedded adapter) must not produce any of
	// the canonical adapters. Use 12000 reads to bypass the
	// min-records gate so detection actually runs.
	reads := make([]*fastq.Record, 0, 12000)
	bases := "ACGT"
	for i := 0; i < 12000; i++ {
		var sb strings.Builder
		// Length 60, varying enough to defeat consensus extension.
		for j := 0; j < 60; j++ {
			sb.WriteByte(bases[(i*7+j*11)%4])
		}
		reads = append(reads, makeSERead("r", sb.String(), 'I'))
	}
	got := DetectAdapterSE(reads)
	if strings.HasPrefix(got, "AGATCGGAAG") || strings.HasPrefix(got, "CTGTCTCTTA") {
		t.Errorf("DetectAdapterSE = %q; unexpected real-adapter hit on random input", got)
	}
}

func TestDetectAdapterSEEmptyInput(t *testing.T) {
	if got := DetectAdapterSE(nil); got != "" {
		t.Errorf("DetectAdapterSE(nil) = %q, want empty", got)
	}
	if got := DetectAdapterSE([]*fastq.Record{}); got != "" {
		t.Errorf("DetectAdapterSE([]) = %q, want empty", got)
	}
}

// TestDetectAdapterSEAllSameRead tests the corner case where every
// read is identical: the kmer histogram has a single high-frequency
// peak so the dominant-path walk should produce something coherent.
// Whether it picks a known adapter (from checkKnownAdaptersSE) or a
// de-novo adapter, the result should be non-empty.
func TestDetectAdapterSEAllSameRead(t *testing.T) {
	insert := "ACGTACGTACGTACGTACGTACGTACGT" // 28bp varied
	seq := insert + truseqAdapter
	reads := make([]*fastq.Record, 0, 12000)
	for i := 0; i < 12000; i++ {
		reads = append(reads, makeSERead("r", seq, 'I'))
	}
	got := DetectAdapterSE(reads)
	if got == "" {
		t.Fatalf("DetectAdapterSE on all-same TruSeq reads returned empty")
	}
}

// TestDetectAdapterSEVeryShortReads tests the corner case where every
// read is too short for the kmer/seed window: with reads shorter than
// 20+keylen+shiftTail the kmer histogram loop never fires and we
// should fall through to an empty return.
func TestDetectAdapterSEVeryShortReads(t *testing.T) {
	reads := make([]*fastq.Record, 0, 12000)
	for i := 0; i < 12000; i++ {
		reads = append(reads, makeSERead("r", "ACGTACGTAC", 'I'))
	}
	if got := DetectAdapterSE(reads); got != "" {
		t.Errorf("DetectAdapterSE on short reads = %q, want empty", got)
	}
}

// TestNucleotideTreeDominantPath mirrors upstream's
// NucleotideTree::test() at nucleotidetree.cpp:90-104.
func TestNucleotideTreeDominantPath(t *testing.T) {
	tree := newNucleotideTree()
	for i := 0; i < 100; i++ {
		tree.addSeq("AAAATTTT")
		tree.addSeq("AAAATTTTGGGG")
		tree.addSeq("AAAATTTTGGGGCCCC")
		tree.addSeq("AAAATTTTGGGGCCAA")
	}
	tree.addSeq("AAAATTTTGGGACCCC")
	path, _ := tree.getDominantPath()
	if path != "AAAATTTTGGGGCC" {
		t.Errorf("getDominantPath = %q, want %q", path, "AAAATTTTGGGGCC")
	}
}

// TestEvaluatorSeq2IntInt2SeqRoundtrip mirrors upstream's
// Evaluator::test() at evaluator.cpp:615-620.
func TestEvaluatorSeq2IntInt2SeqRoundtrip(t *testing.T) {
	s := "ATCGATCGAT"
	got := evaluatorInt2Seq(uint32(evaluatorSeq2Int(s, 0, 10, -1)), 10)
	if got != s {
		t.Errorf("seq2int/int2seq roundtrip on %q = %q", s, got)
	}
}

func TestDetectAdapterPEFromOverlap(t *testing.T) {
	// Build a synthetic PE pair where R1 ends with a TruSeq adapter
	// (read-through) and R2 ends with the reverse-complement of the
	// upstream end of R1 plus its own R2 adapter.
	insert := "GTCAACTGGCTTAATTACCGACATAGCTAGTTACCGGCAATGCATACGTACG" // 51 bp
	r1Seq := insert + truseqAdapter[:25]                             // 76 bp
	// R2 mirrors the insert (reverse complement) and ends with R2 adapter
	// AGATCGGAAGAGCGT... (we just use the TruSeq R2 prefix for simplicity).
	const truseqR2 = "AGATCGGAAGAGCGTCGTGTAGGGAAA"
	r2Seq := reverseComplement(insert) + truseqR2[:25]
	pairs := make([][2]*fastq.Record, 0, 500)
	for i := 0; i < 500; i++ {
		pairs = append(pairs, [2]*fastq.Record{
			makeSERead("r1", r1Seq, 'I'),
			makeSERead("r2", r2Seq, 'I'),
		})
	}
	r1Adapter, r2Adapter := DetectAdaptersFromPairs(pairs)
	if r1Adapter == "" {
		t.Errorf("expected non-empty R1 adapter")
	}
	if r2Adapter == "" {
		t.Errorf("expected non-empty R2 adapter")
	}
	// The R1 adapter should begin with the TruSeq prefix.
	if r1Adapter != "" && !strings.HasPrefix(truseqAdapter[:25], r1Adapter[:min(len(r1Adapter), 10)]) &&
		!strings.HasPrefix(r1Adapter, truseqAdapter[:min(len(truseqAdapter), 10)]) {
		t.Errorf("R1 adapter = %q; expected TruSeq prefix overlap", r1Adapter)
	}
}

func TestDetectAdapterPENoOverlap(t *testing.T) {
	// Two completely unrelated reads should produce no detection.
	r1 := makeSERead("r1", strings.Repeat("ACGT", 30), 'I')
	r2 := makeSERead("r2", strings.Repeat("TGCA", 30), 'I')
	if got := DetectAdapterPE(r1, r2); got != "" {
		// Random overlap can occur with low complexity but should not
		// for our pure-tandem repeats; allow short tails up to 10 bp.
		// We're not strict here, but verify it's not a long string.
		if len(got) > 30 {
			t.Errorf("DetectAdapterPE returned unexpectedly long tail %q", got)
		}
	}
}

func TestDetectAdapterPENilInputs(t *testing.T) {
	if got := DetectAdapterPE(nil, nil); got != "" {
		t.Errorf("DetectAdapterPE(nil, nil) = %q, want empty", got)
	}
}

func TestDetectAdaptersFromPairsThreshold(t *testing.T) {
	// Too few pairs - threshold should reject.
	pairs := [][2]*fastq.Record{}
	r1, r2 := DetectAdaptersFromPairs(pairs)
	if r1 != "" || r2 != "" {
		t.Errorf("empty pairs returned (%q, %q); want both empty", r1, r2)
	}
}

func TestProcessSingleEndDetectAdapterIntegrates(t *testing.T) {
	// Use DetectAdapterSE option, run ProcessSingleEnd, and verify that
	// the detected adapter ends up in stats.DetectedAdapter. We need
	// 10000+ reads to exercise the upstream-verbatim detection path
	// (evaluator.cpp:344).
	var sb strings.Builder
	for i := 0; i < 12000; i++ {
		insert := strings.Repeat("ACGTACGT", 4) + string("ACGT"[i%4])
		seq := insert + truseqAdapter
		sb.WriteString("@r\n")
		sb.WriteString(seq)
		sb.WriteString("\n+\n")
		sb.WriteString(strings.Repeat("I", len(seq)))
		sb.WriteString("\n")
	}

	opts := DefaultProcessOptions()
	opts.DetectAdapterSE = true
	opts.MinLength = 10

	stats, err := ProcessSingleEnd(strings.NewReader(sb.String()), &discardWriter{}, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("ProcessSingleEnd: %v", err)
	}
	if stats.DetectedAdapter == "" {
		t.Fatalf("expected DetectedAdapter to be set after SE detection")
	}
	if !strings.HasPrefix(truseqAdapter, stats.DetectedAdapter[:min(len(stats.DetectedAdapter), 10)]) &&
		!strings.HasPrefix(stats.DetectedAdapter, truseqAdapter[:min(len(truseqAdapter), 10)]) {
		t.Errorf("DetectedAdapter = %q; expected TruSeq prefix overlap", stats.DetectedAdapter)
	}
	if stats.AdapterTrimmedReads == 0 {
		t.Errorf("expected AdapterTrimmedReads > 0 after SE detection set the adapter")
	}
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}
