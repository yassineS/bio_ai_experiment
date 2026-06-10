package fastp

import (
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fastq"
)

// TestCorrectByOverlapAnalysis_UpstreamFixture mirrors upstream's
// BaseCorrector::test (basecorrector.cpp:85-107): a hand-built pair where
// exactly one base in each mate is corrected from the other, and the result
// sequences/qualities are known.
func TestCorrectByOverlapAnalysis_UpstreamFixture(t *testing.T) {
	r1 := &fastq.Record{
		Sequence: []byte("TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAATTTTCCACGGGG"),
		Quality:  []byte("EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE/EEEEE"),
	}
	r2 := &fastq.Record{
		Sequence: []byte("AAAAAAAAAACCCCGGGGAAAATTTTAAAATTGGGGGGGGGGTGGGGGGGGGGGGG"),
		Quality:  []byte("EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE/EEEEEEEEEEEEE"),
	}
	rcSeq2 := reverseComplement(string(r2.Sequence))
	ov := analyzeOverlapPair(string(r1.Sequence), rcSeq2, 5, 30, 0.2)
	if !ov.Overlapped {
		t.Fatalf("expected overlap to be detected")
	}
	corrected, reads := correctByOverlapAnalysis(r1, r2, ov, fastq.Phred33)
	if corrected != 2 {
		t.Fatalf("corrected bases = %d, want 2", corrected)
	}
	if reads != 2 {
		t.Fatalf("corrected reads = %d, want 2", reads)
	}
	wantR1 := "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAATTTTCCCCGGGG"
	wantR2 := "AAAAAAAAAACCCCGGGGAAAATTTTAAAATTGGGGGGGGGGGGGGGGGGGGGGGG"
	if string(r1.Sequence) != wantR1 {
		t.Fatalf("r1 = %s\nwant %s", r1.Sequence, wantR1)
	}
	if string(r2.Sequence) != wantR2 {
		t.Fatalf("r2 = %s\nwant %s", r2.Sequence, wantR2)
	}
	wantQual := "EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE"
	if string(r1.Quality) != wantQual || string(r2.Quality) != wantQual {
		t.Fatalf("qualities not normalised:\n r1q=%s\n r2q=%s", r1.Quality, r2.Quality)
	}
}

// TestCorrectByOverlapAnalysis_NoDiff verifies that a perfectly matching
// overlap (Diff == 0) is left untouched, matching upstream's early return.
func TestCorrectByOverlapAnalysis_NoDiff(t *testing.T) {
	seq := "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
	r1 := &fastq.Record{Sequence: []byte(seq), Quality: []byte(rep('I', len(seq)))}
	r2 := &fastq.Record{Sequence: []byte(reverseComplement(seq)), Quality: []byte(rep('I', len(seq)))}
	ov := analyzeOverlapPair(seq, reverseComplement(string(r2.Sequence)), 5, 30, 0.2)
	corrected, reads := correctByOverlapAnalysis(r1, r2, ov, fastq.Phred33)
	if corrected != 0 || reads != 0 {
		t.Fatalf("expected no correction on identical overlap, got %d/%d", corrected, reads)
	}
}

func TestComplementBase(t *testing.T) {
	cases := map[byte]byte{'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'a': 'T', 'N': 'N', 'X': 'N'}
	for in, want := range cases {
		if got := complementBase(in); got != want {
			t.Errorf("complementBase(%c) = %c, want %c", in, got, want)
		}
	}
}

func rep(b byte, n int) string {
	out := make([]byte, n)
	for i := range out {
		out[i] = b
	}
	return string(out)
}
