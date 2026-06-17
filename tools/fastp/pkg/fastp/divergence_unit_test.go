package fastp

// Binary-free UNIT tests for the three formerly-divergent fastp algorithms.
//
// These pin the pure helpers directly, with NO dependency on the upstream
// fastp binary or the reference_code submodule, so they pass on a bare CI
// checkout. The expected values are taken from the upstream C++ in-source
// self-tests (polyx.cpp::test, filter.cpp::test, adaptertrimmer.cpp::test)
// or are hand-derived from the documented upstream algorithm — never from a
// frozen Go golden.
//
// The matching live byte-parity / similarity assertions against the upstream
// binary live in parity_test.go (Cases 12, 13, 14, 16).

import (
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fastq"
)

// TestUnitTrimPolyG checks the mismatch-tolerant 3' poly-G scan
// (polyx.cpp:16-42): 1 mismatch allowed per 8 bases scanned, capped at 5,
// anchoring the trim on the last-G position. trimPolyG returns the index at
// which to truncate (the new length).
func TestUnitTrimPolyG(t *testing.T) {
	cases := []struct {
		name       string
		seq        string
		compareReq int
		wantLen    int
	}{
		{
			// Pure poly-G tail of 20 after a G-free insert. The trim anchors
			// on the leftmost G of the run; since the insert "ACTACTACTA" has
			// no G, the run starts exactly at index 10.
			name:       "clean_polyg_tail",
			seq:        "ACTACTACTA" + strings.Repeat("G", 20),
			compareReq: 10,
			wantLen:    10,
		},
		{
			// One embedded mismatch (a T) inside the poly-G tail is tolerated
			// (1 mismatch over >=8 bases scanned). The whole tail back to the
			// first G is still trimmed.
			name:       "one_mismatch_in_tail",
			seq:        "ACTACTACTA" + "GGGGGGGGTGGGGGGGGGGG",
			compareReq: 10,
			wantLen:    10,
		},
		{
			// No poly-G tail at all: nothing trimmed.
			name:       "no_polyg",
			seq:        "ACGTACGTACGTACGTACGTACGTACGT",
			compareReq: 10,
			wantLen:    28,
		},
		{
			// Tail shorter than compareReq is not trimmed (i < compareReq).
			name:       "short_tail",
			seq:        "ACGTACGTACGTACGT" + "GGGG",
			compareReq: 10,
			wantLen:    20,
		},
	}
	for _, tc := range cases {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			got := trimPolyG(tc.seq, tc.compareReq)
			if got != tc.wantLen {
				t.Fatalf("trimPolyG(%q, %d) = %d, want %d", tc.seq, tc.compareReq, got, tc.wantLen)
			}
		})
	}
}

// TestUnitTrimPolyX reproduces upstream's PolyX::test() golden exactly
// (polyx.cpp:118-130): the documented input must trim to "ATTTT" with 51
// bases removed. We also check a clean poly-A tail and a no-poly read.
func TestUnitTrimPolyX(t *testing.T) {
	// Upstream golden (polyx.cpp:120-129).
	const upstreamSeq = "ATTTTAAAAAAAAAATAAAAAAAAAAAAACAAAAAAAAAAAAAAAAAAAAAAAAAT"
	newLen, trimmed := trimPolyXUpstream(upstreamSeq, 10)
	if got := upstreamSeq[:newLen]; got != "ATTTT" {
		t.Fatalf("trimPolyXUpstream upstream golden = %q (len %d), want %q", got, newLen, "ATTTT")
	}
	if trimmed != 51 {
		t.Fatalf("trimPolyXUpstream upstream golden trimmed = %d, want 51", trimmed)
	}

	// Clean poly-A tail.
	seq := "ACGTACGTACGT" + strings.Repeat("A", 20)
	newLen, trimmed = trimPolyXUpstream(seq, 10)
	if newLen != 12 || trimmed != 20 {
		t.Fatalf("clean poly-A: newLen=%d trimmed=%d, want 12/20", newLen, trimmed)
	}

	// No poly tail: no trim.
	seq = "ACGTACGTACGTACGTACGTACGT"
	newLen, trimmed = trimPolyXUpstream(seq, 10)
	if newLen != len(seq) || trimmed != 0 {
		t.Fatalf("no-poly: newLen=%d trimmed=%d, want %d/0", newLen, trimmed, len(seq))
	}
}

// TestUnitSlidingWindowCut reproduces upstream's Filter::test() golden
// (filter.cpp:260-279) for cut_front + cut_tail with window=4, quality=20,
// then checks cut_right and the N-skip behaviour. slidingWindowCut returns
// the [lo, hi) range to keep; here we always pass the full read (the caller
// pre-clips front/tail), so we emulate the upstream "tail=1" golden by
// dropping the last base before calling.
func TestUnitSlidingWindowCut(t *testing.T) {
	enc := fastq.Phred33

	t.Run("cut_front_tail_upstream_golden", func(t *testing.T) {
		// Upstream golden runs trimAndCut(&r, front=0, tail=1). Emulate the
		// tail=1 hard-trim by removing the final base up front.
		seq := "TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTT"
		qual := "/////CCCCCCCCCCCC////CCCCCCCCCCCCCC////E"
		seq = seq[:len(seq)-1]
		qual = qual[:len(qual)-1]
		opts := DefaultProcessOptions()
		opts.CutFront = true
		opts.CutTail = true
		opts.CutWindowSize = 4
		opts.CutMeanQuality = 20
		lo, hi := slidingWindowCut([]byte(seq), []byte(qual), enc, opts)
		gotSeq := seq[lo:hi]
		gotQual := qual[lo:hi]
		if gotSeq != "CCCCCCCCCCCCCCCCCCCCCCCCCCCC" {
			t.Fatalf("cut_front+tail seq = %q, want %q", gotSeq, "CCCCCCCCCCCCCCCCCCCCCCCCCCCC")
		}
		if gotQual != "CCCCCCCCCCC////CCCCCCCCCCCCC" {
			t.Fatalf("cut_front+tail qual = %q, want %q", gotQual, "CCCCCCCCCCC////CCCCCCCCCCCCC")
		}
	})

	t.Run("cut_right_keeps_high_q_prefix", func(t *testing.T) {
		// 30 high-Q bases then a low-Q tail. cut_right should find the first
		// low window and keep the leading high-Q bases of that window.
		seq := strings.Repeat("A", 40)
		qual := strings.Repeat("I", 30) + strings.Repeat("#", 10) // I=40, #=2
		opts := DefaultProcessOptions()
		opts.CutRight = true
		opts.CutWindowSize = 4
		opts.CutMeanQuality = 20
		lo, hi := slidingWindowCut([]byte(seq), []byte(qual), enc, opts)
		if lo != 0 {
			t.Fatalf("cut_right lo = %d, want 0", lo)
		}
		// The high-Q prefix (30 'I') is kept; the low-Q tail is cut.
		if hi != 30 {
			t.Fatalf("cut_right hi = %d, want 30", hi)
		}
	})

	t.Run("cut_front_skips_leading_N", func(t *testing.T) {
		// Low-Q front, then an N at the window boundary, then high-Q. The
		// front advance must skip the N (filter.cpp:138-139).
		seq := "NNNN" + strings.Repeat("A", 36)
		qual := strings.Repeat("#", 4) + strings.Repeat("I", 36)
		opts := DefaultProcessOptions()
		opts.CutFront = true
		opts.CutWindowSize = 4
		opts.CutMeanQuality = 20
		lo, hi := slidingWindowCut([]byte(seq), []byte(qual), enc, opts)
		if lo < 4 {
			t.Fatalf("cut_front lo = %d, expected to skip the leading Ns (>=4)", lo)
		}
		if seq[lo] == 'N' {
			t.Fatalf("cut_front kept a leading N at lo=%d", lo)
		}
		_ = hi
	})
}

// TestUnitDetectAdapterSE checks the SE kmer/nucleotide-tree auto-detector
// (evaluator.cpp) WITHOUT the upstream binary: a synthetic set of >=10000
// reads carrying the TruSeq adapter must recover it; below the 10000-record
// gate it must return ""; random reads must not produce a false adapter.
func TestUnitDetectAdapterSE(t *testing.T) {
	const adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"

	// Deterministic LCG so the test needs no external seed/source.
	var state uint64 = 0x9E3779B97F4A7C15
	next := func() uint64 {
		state = state*6364136223846793005 + 1442695040888963407
		return state >> 33
	}
	bases := []byte("ACGT")
	randSeq := func(n int) string {
		b := make([]byte, n)
		for i := range b {
			b[i] = bases[next()%4]
		}
		return string(b)
	}
	mkReads := func(n int, withAdapter bool) []*fastq.Record {
		out := make([]*fastq.Record, 0, n)
		for i := 0; i < n; i++ {
			var seq string
			if withAdapter && i%10 != 0 {
				insertLen := 20 + int(next()%50)
				seq = randSeq(insertLen) + adapter
				if len(seq) < 100 {
					seq += randSeq(100 - len(seq))
				}
				seq = seq[:100]
			} else {
				seq = randSeq(100)
			}
			q := make([]byte, len(seq))
			for j := range q {
				q[j] = 'I'
			}
			out = append(out, &fastq.Record{ID: "r", Sequence: []byte(seq), Quality: q})
		}
		return out
	}

	t.Run("recovers_truseq_above_gate", func(t *testing.T) {
		got := DetectAdapterSE(mkReads(12000, true))
		if got == "" {
			t.Fatalf("DetectAdapterSE returned empty on 12000 adapter reads")
		}
		// The detected string must be the TruSeq adapter or a prefix of it
		// (the de-novo walk can stop a base or two early).
		if !strings.HasPrefix(adapter, got) && !strings.HasPrefix(got, "AGATCGGAAGAGC") {
			t.Fatalf("detected %q is not a TruSeq prefix", got)
		}
	})

	t.Run("below_gate_returns_empty", func(t *testing.T) {
		// 5000 < 10000 record gate (evaluator.cpp:344).
		if got := DetectAdapterSE(mkReads(5000, true)); got != "" {
			t.Fatalf("DetectAdapterSE = %q below the 10000-record gate, want \"\"", got)
		}
	})

	t.Run("random_reads_no_false_adapter", func(t *testing.T) {
		if got := DetectAdapterSE(mkReads(12000, false)); got != "" {
			t.Fatalf("DetectAdapterSE = %q on random reads, want \"\"", got)
		}
	})
}
