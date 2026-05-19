package skewer

// Unit tests for the two algorithms ported from
// reference_code/skewer/src/matrix.cpp:
//
//  - findAdapterWithQual: the quality-weighted "SW-with-tail-penalty"
//    matcher (cAdapter::align + cMatrix::penalty[]).
//  - calcRevCompScore + detectPairedTrim: the matrix-mode paired-end
//    overlap gate (cMatrix::findAdapterWithPE + CalcRevCompScore).
//
// Each test maps to a specific edge case described in the upstream code
// or in tools/PARITY_VALIDATION.md.

import (
	"bytes"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fastq"
)

// makeQual returns a Phred-33 quality buffer of length n, all set to the
// given quality byte. Helper for the find/RC tests.
func makeQual(n int, q byte) []byte {
	out := make([]byte, n)
	for i := range out {
		out[i] = q
	}
	return out
}

// TestFindAdapterWithQual_ExactMatch confirms the matcher locates a perfect
// adapter and returns its leftmost index.
func TestFindAdapterWithQual_ExactMatch(t *testing.T) {
	seq := "ACGTACGTACGTAGATCGGAAGAGCNN"
	adp := "AGATCGGAAGAGC"
	q := makeQual(len(seq), 'I')
	got := findAdapterWithQual(seq, adp, q, 3, 0.1)
	if got != 12 {
		t.Fatalf("exact match pos: got %d, want 12", got)
	}
}

// TestFindAdapterWithQual_MismatchHighQualRejected reproduces case05: a
// single mismatch at Phred 'I' (Q40) is enough to exceed dMaxPenalty when
// adapter is 13 bp (0.1 * MEAN_PENALTY * 13 ≈ 3.22 < MAX_PENALTY=4.477).
func TestFindAdapterWithQual_MismatchHighQualRejected(t *testing.T) {
	// Adapter has G at index 9; read has C there (1 mismatch at Q40).
	seq := "ACGTACGTACGTAGATCGGAACAGC"
	adp := "AGATCGGAAGAGC"
	q := makeQual(len(seq), 'I')
	got := findAdapterWithQual(seq, adp, q, 3, 0.1)
	if got != -1 {
		t.Fatalf("Q40 mismatch should be rejected: got pos=%d, want -1", got)
	}
}

// TestFindAdapterWithQual_MismatchLowQualAccepted: at Q5 (chr 38) the per-
// mismatch penalty is minPenalty+0.5 = 0.977, well under dMaxPenalty=3.22,
// so a 1-mismatch match should be accepted.
func TestFindAdapterWithQual_MismatchLowQualAccepted(t *testing.T) {
	seq := "ACGTACGTACGTAGATCGGAACAGC"
	adp := "AGATCGGAAGAGC"
	q := makeQual(len(seq), '&') // Q5
	got := findAdapterWithQual(seq, adp, q, 3, 0.1)
	if got != 12 {
		t.Fatalf("low-quality mismatch should be accepted: got pos=%d, want 12", got)
	}
}

// TestFindAdapterWithQual_MismatchAtStart: mismatch at the first base of
// the adapter, with Q40 quality — should still be rejected (single high-Q
// mismatch exceeds the threshold regardless of position).
func TestFindAdapterWithQual_MismatchAtStart(t *testing.T) {
	// Adapter is AGATCGGAAGAGC; read has TGATCGGAAGAGC at pos 12.
	seq := "ACGTACGTACGTTGATCGGAAGAGC"
	adp := "AGATCGGAAGAGC"
	q := makeQual(len(seq), 'I')
	got := findAdapterWithQual(seq, adp, q, 3, 0.1)
	if got != -1 {
		t.Fatalf("start-mismatch at Q40: got pos=%d, want -1", got)
	}
}

// TestFindAdapterWithQual_TwoMismatchesRejected: two mismatches at any
// quality exceed the threshold by an even wider margin.
func TestFindAdapterWithQual_TwoMismatchesRejected(t *testing.T) {
	// Adapter AGATCGGAAGAGC; read has two diffs at idx 9 and 11.
	seq := "ACGTACGTACGTAGATCGGAACATC"
	adp := "AGATCGGAAGAGC"
	q := makeQual(len(seq), 'I')
	got := findAdapterWithQual(seq, adp, q, 3, 0.1)
	if got != -1 {
		t.Fatalf("two Q40 mismatches: got pos=%d, want -1", got)
	}
}

// TestFindAdapterWithQual_EmptyAdapter and minOverlap guard behaviour.
func TestFindAdapterWithQual_EmptyAndShort(t *testing.T) {
	got := findAdapterWithQual("ACGT", "", nil, 3, 0.1)
	if got != -1 {
		t.Errorf("empty adapter: got %d, want -1", got)
	}
	got = findAdapterWithQual("AC", "AGATC", nil, 3, 0.1)
	if got != -1 {
		t.Errorf("seq shorter than minOverlap: got %d, want -1", got)
	}
}

// TestFindAdapterWithQual_PartialTail: a 5-bp partial adapter at the very
// 3' end is acceptable when minOverlap=3 and the bases all match.
func TestFindAdapterWithQual_PartialTail(t *testing.T) {
	seq := "ACGTACGTACGTAGATC"
	adp := "AGATCGGAAGAGC"
	q := makeQual(len(seq), 'I')
	got := findAdapterWithQual(seq, adp, q, 3, 0.1)
	if got != 12 {
		t.Fatalf("partial tail match: got pos=%d, want 12", got)
	}
}

// TestCalcRevCompScore_Overlap: R1 prefix is the literal reverse-complement
// of R2 prefix → CalcRevCompScore should return true.
func TestCalcRevCompScore_Overlap(t *testing.T) {
	// R1 prefix AAAAA, R2 prefix TTTTT → revcomp(R2)=AAAAA, matches R1.
	seq1 := "AAAAAGATCGGAAGAGC"
	seq2 := "TTTTTGGGGGGGGGGGG"
	q := makeQual(len(seq1), 'I')
	if !calcRevCompScore(seq1, seq2, q, q, 5, 0.1) {
		t.Fatal("identical RC prefixes should pass")
	}
}

// TestCalcRevCompScore_NoOverlap: R1 and R2 prefixes are unrelated (every
// position mismatches under RC) → must return false.
func TestCalcRevCompScore_NoOverlap(t *testing.T) {
	seq1 := "ACGTACGT"
	seq2 := "GCGCGCGC"
	q := makeQual(len(seq1), 'I')
	if calcRevCompScore(seq1, seq2, q, q, 8, 0.1) {
		t.Fatal("unrelated prefixes should fail RC check")
	}
}

// TestCalcRevCompScore_PartialOverlap: case04 scenario — R1 prefix is
// palindromic ACGTACGTACGT, R2 prefix is TGCATGCATGCA (also palindromic),
// but they are NOT reverse-complements of each other, so the RC check
// must fail (which is what makes case04 a pass-through in upstream).
func TestCalcRevCompScore_PartialOverlap(t *testing.T) {
	r1 := "ACGTACGTACGTAGATCGGAAGAGCNN"
	r2 := "TGCATGCATGCAAGATCGGAAGAGCNN"
	q := makeQual(len(r1), 'I')
	if calcRevCompScore(r1, r2, q, q, 12, 0.1) {
		t.Fatal("case04 prefixes should NOT pass RC check (this is what makes it untrimmed upstream)")
	}
}

// TestDetectPairedTrim_BothNoAdapter: when neither mate sees the adapter,
// the gate returns (-1,-1) and the mates pass through.
func TestDetectPairedTrim_BothNoAdapter(t *testing.T) {
	tp1, tp2 := detectPairedTrim("ACGT", "TGCA", nil, nil, -1, -1, 0.1)
	if tp1 != -1 || tp2 != -1 {
		t.Errorf("no adapter: got (%d,%d), want (-1,-1)", tp1, tp2)
	}
}

// TestDetectPairedTrim_OneAdapterRCSupport: R1 has adapter at pos 5 and the
// 5-base prefix is the RC of R2's 5-base prefix. Both mates should be
// trimmed at 5.
func TestDetectPairedTrim_OneAdapterRCSupport(t *testing.T) {
	r1 := "AAAAAGATCGGAAGAGC"
	r2 := "TTTTTGGGGGGGGGGGG"
	q := makeQual(len(r1), 'I')
	tp1, tp2 := detectPairedTrim(r1, r2, q, q, 5, -1, 0.1)
	if tp1 != 5 || tp2 != 5 {
		t.Errorf("RC-supported trim: got (%d,%d), want (5,5)", tp1, tp2)
	}
}

// TestDetectPairedTrim_OneAdapterNoRCSupport: R1 has adapter at 12 but the
// prefixes aren't reverse-complements (case04 second pair). Gate returns
// (-1,-1) so both mates stay untrimmed. R1 prefix is `ACGTACGTACGT` and R2
// prefix is `TGCATGCATGCA` — both happen to be palindromes individually,
// but reverse-complement of R2 prefix is `TGCATGCATGCA` (not equal to R1's
// `ACGTACGTACGT`), so the gate rejects.
func TestDetectPairedTrim_OneAdapterNoRCSupport(t *testing.T) {
	r1 := "ACGTACGTACGTAGATCGGAAGAGCNN"
	r2 := "TGCATGCATGCAAGATCGGAAGAGCNN"
	q := makeQual(len(r1), 'I')
	tp1, tp2 := detectPairedTrim(r1, r2, q, q, 12, 12, 0.1)
	if tp1 != -1 || tp2 != -1 {
		t.Errorf("no RC support: got (%d,%d), want (-1,-1)", tp1, tp2)
	}
}

// TestTrimPairedEnd_MatrixModeGate is an end-to-end check that matrix mode
// gates trimming on the RC overlap. With PEMatrixMode=true and the case04
// inputs, both mates should pass through unchanged.
func TestTrimPairedEnd_MatrixModeGate(t *testing.T) {
	in1 := "@r1/1\nACGTACGTACGTAGATCGGAAGAGCNN\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
	in2 := "@r1/2\nTGCATGCATGCAAGATCGGAAGAGCNN\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
	var out1, out2 bytes.Buffer
	opts := TrimOptions{
		Adapter3:     "AGATCGGAAGAGC",
		MinLength:    8,
		MinOverlap:   3,
		ErrorRate:    0.1,
		PEMatrixMode: true,
	}
	_, err := TrimPairedEnd(strings.NewReader(in1), strings.NewReader(in2), &out1, &out2, nil, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("TrimPairedEnd: %v", err)
	}
	if !bytes.Contains(out1.Bytes(), []byte("ACGTACGTACGTAGATCGGAAGAGCNN")) {
		t.Errorf("matrix mode should leave R1 untrimmed; got:\n%s", out1.String())
	}
	if !bytes.Contains(out2.Bytes(), []byte("TGCATGCATGCAAGATCGGAAGAGCNN")) {
		t.Errorf("matrix mode should leave R2 untrimmed; got:\n%s", out2.String())
	}
}

// TestTrimPairedEnd_MatrixModeAllowsTrim: when the RC overlap holds (R1
// prefix AAAA, R2 prefix TTTT), matrix mode should permit the trim.
// Note: the adapter `AGATCGGAAGAGC` starts at index 4 (not 5) because the
// read's 5th char 'A' is the first 'A' of the adapter — so the prefix is
// 4 chars long, not 5.
func TestTrimPairedEnd_MatrixModeAllowsTrim(t *testing.T) {
	in1 := "@p/1\nAAAAAGATCGGAAGAGC\n+\nIIIIIIIIIIIIIIIII\n"
	in2 := "@p/2\nTTTTTGGGGGGGGGGGG\n+\nIIIIIIIIIIIIIIIII\n"
	var out1, out2 bytes.Buffer
	opts := TrimOptions{
		Adapter3:     "AGATCGGAAGAGC",
		MinLength:    3,
		MinOverlap:   3,
		ErrorRate:    0.1,
		PEMatrixMode: true,
	}
	_, err := TrimPairedEnd(strings.NewReader(in1), strings.NewReader(in2), &out1, &out2, nil, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("TrimPairedEnd: %v", err)
	}
	if !bytes.Contains(out1.Bytes(), []byte("AAAA\n")) {
		t.Errorf("matrix mode with RC support should trim R1 to 4 bases; got:\n%s", out1.String())
	}
	if !bytes.Contains(out2.Bytes(), []byte("TTTT\n")) {
		t.Errorf("matrix mode with RC support should trim R2 to 4 bases; got:\n%s", out2.String())
	}
}

// TestMismatchPenalty_QualityRamp confirms the penalty[] table mirrors
// upstream's matrix.cpp:547-556: Q<=0 → minPenalty, Q40+ → maxPenalty,
// linear 0.1 increments between.
func TestMismatchPenalty_QualityRamp(t *testing.T) {
	cases := []struct {
		qual byte
		want float64
	}{
		{'!', minPenalty}, // chr 33, Q0
		{'\'' /* chr 39 */, minPenalty + 0.6},
		{'I', maxPenalty}, // chr 73, Q40+ → max
		{'J', maxPenalty}, // chr 74 also max
	}
	for _, c := range cases {
		got := mismatchPenalty([]byte{c.qual}, 0)
		if got != c.want {
			t.Errorf("mismatchPenalty(%q)=%g, want %g", c.qual, got, c.want)
		}
	}
	if got := mismatchPenalty(nil, 0); got != meanPenalty {
		t.Errorf("nil qual fallback: got %g, want %g", got, meanPenalty)
	}
}

// TestComplementBase covers ACGT (case-insensitive) and N for unknowns.
func TestComplementBase(t *testing.T) {
	cases := map[byte]byte{
		'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
		'a': 'T', 't': 'A', 'c': 'G', 'g': 'C',
		'N': 'N', 'X': 'N',
	}
	for in, want := range cases {
		if got := complementBase(in); got != want {
			t.Errorf("complementBase(%q)=%q, want %q", in, got, want)
		}
	}
}
