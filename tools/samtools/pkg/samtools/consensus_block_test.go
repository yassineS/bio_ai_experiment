package samtools

import (
	"bytes"
	"fmt"
	"os"
	"strings"
	"testing"
)

// makeMultiBlockSAM builds a coordinate-spanning SAM whose single contig is
// long enough to cross several consensusBlockWidth boundaries, with reads
// placed in the first block, straddling the block-0/block-1 boundary, and in
// later blocks. The straddling read (its [start,end) crosses an exact multiple
// of consensusBlockWidth) is the case the block gather must carry into the next
// block; if the carry or readIdx rebase were wrong, the streaming output would
// diverge from the whole-buffer walk.
func makeMultiBlockSAM() string {
	const seq30 = "ACGTACGTACGTACGTACGTACGTACGTAC"
	const qual30 = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
	// Contig spanning > 2 blocks so at least two interior block edges exist.
	refLen := 2*consensusBlockWidth + 4096
	var b strings.Builder
	fmt.Fprintf(&b, "@HD\tVN:1.6\tSO:unsorted\n")
	fmt.Fprintf(&b, "@SQ\tSN:chr1\tLN:%d\n", refLen)
	// rec emits one 30M read at 1-based pos p (with an MD tag so the bayesian
	// NM-halo path is exercised, matching the existing indexed fixture).
	rec := func(name string, pos int) {
		fmt.Fprintf(&b, "%s\t0\tchr1\t%d\t60\t30M\t*\t0\t0\t%s\t%s\tNM:i:0\tMD:Z:30\n",
			name, pos, seq30, qual30)
	}
	// Block 0 interior coverage (overlapping reads so calls are non-trivial).
	rec("b0_a", 1000)
	rec("b0_b", 1010)
	rec("b0_c", 1025)
	// Reads straddling the first block boundary (1-based pos chosen so the read
	// covers positions on both sides of the block-0/block-1 edge).
	edge1 := consensusBlockWidth // 0-based edge; 1-based start a touch before it
	rec("edge1_a", edge1-10)     // covers [edge1-11 .. edge1+18], crosses edge1
	rec("edge1_b", edge1-5)
	rec("edge1_c", edge1+3)
	// Block 1 interior coverage.
	rec("b1_a", consensusBlockWidth+5000)
	rec("b1_b", consensusBlockWidth+5012)
	// Reads straddling the second block boundary.
	edge2 := 2 * consensusBlockWidth
	rec("edge2_a", edge2-12)
	rec("edge2_b", edge2-2)
	// Block 2 interior coverage.
	rec("b2_a", 2*consensusBlockWidth+1000)
	rec("b2_b", 2*consensusBlockWidth+1015)
	return b.String()
}

// TestConsensusBlockChunkedMatchesWholeBuffer pins the task-#50 invariant: the
// coordinate-blocked, streaming consensus engine (the indexed `consensus -r`
// fast path, which processes the contig in consensusBlockWidth blocks and only
// keeps ~one block's reads resident) emits output byte-identical to the
// whole-buffer linear walk (Consensus -> bucketByChrom, which holds every read
// of the contig at once). The synthetic BAM spans multiple block boundaries and
// includes reads that straddle a block edge, so the block carry and the
// block-local readIdx rebase are both exercised. If either were wrong the two
// outputs would differ.
func TestConsensusBlockChunkedMatchesWholeBuffer(t *testing.T) {
	bam := makeIndexedBAM(t, writeSAMFile(t, makeMultiBlockSAM()))

	base := func() ConsensusOptions { return ConsensusOptions{Input: bam} }
	cases := []struct {
		name   string
		mutate func(*ConsensusOptions)
	}{
		{"fasta-simple-wholechrom", func(o *ConsensusOptions) {
			o.Format, o.Mode, o.Regions = ConsensusFASTA, ConsensusModeSimple, []string{"chr1"}
		}},
		{"pileup-simple-wholechrom", func(o *ConsensusOptions) {
			o.Format, o.Mode, o.Regions = ConsensusPileup, ConsensusModeSimple, []string{"chr1"}
		}},
		{"fasta-bayesian-wholechrom", func(o *ConsensusOptions) {
			o.Format, o.Mode, o.Regions = ConsensusFASTA, ConsensusModeBayesian, []string{"chr1"}
		}},
		{"pileup-bayesian-wholechrom", func(o *ConsensusOptions) {
			o.Format, o.Mode, o.Regions = ConsensusPileup, ConsensusModeBayesian, []string{"chr1"}
		}},
		{"fasta-bayesian-spanning-region", func(o *ConsensusOptions) {
			// A sub-region whose window spans a block boundary; exercises the
			// streamRecSource drainRest path (window end < contig length) too.
			o.Format, o.Mode = ConsensusFASTA, ConsensusModeBayesian
			o.Regions = []string{fmt.Sprintf("chr1:%d-%d",
				consensusBlockWidth-100, consensusBlockWidth+100)}
		}},
	}

	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			// Indexed/streaming fast path via ConsensusFile.
			idxOpts := base()
			tc.mutate(&idxOpts)
			var got bytes.Buffer
			if err := ConsensusFile(idxOpts, &got, nil); err != nil {
				t.Fatalf("ConsensusFile (block-chunked indexed): %v", err)
			}

			// Whole-buffer linear oracle: Consensus over the whole file, which
			// buckets every read of the contig at once (bucketByChrom).
			f, err := os.Open(bam)
			if err != nil {
				t.Fatal(err)
			}
			linOpts := base()
			tc.mutate(&linOpts)
			var want bytes.Buffer
			if err := Consensus(f, &want, linOpts); err != nil {
				_ = f.Close()
				t.Fatalf("Consensus (whole-buffer linear): %v", err)
			}
			_ = f.Close()

			if got.String() != want.String() {
				t.Errorf("block-chunked output differs from whole-buffer\n--- chunked (len %d) ---\n%s\n--- whole (len %d) ---\n%s",
					got.Len(), trunc(got.String()), want.Len(), trunc(want.String()))
			}
		})
	}
}

// trunc shortens a possibly-large output for a readable test failure message.
func trunc(s string) string {
	const max = 2000
	if len(s) <= max {
		return s
	}
	return s[:max] + "\n...[truncated]..."
}
