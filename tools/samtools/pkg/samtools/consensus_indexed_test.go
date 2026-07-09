package samtools

import (
	"bytes"
	"os"
	"runtime/debug"
	"testing"
)

// consensusIndexedSAM is a small coordinate-sortable SAM exercising the indexed
// `consensus -r` fast path. It deliberately includes reads that START before a
// sub-region but span INTO it (read_a at pos 5, read_b at pos 12), so the test
// also pins that the index fast path keeps overlap reads exactly as htslib's
// iterator does — the same reads the linear post-filter retains. Several reads
// carry an MD aux tag so the bayesian NM-halo path (which reads MD before the
// per-record aux payload is released) is exercised too.
const consensusIndexedSAM = `@HD	VN:1.6	SO:unsorted
@SQ	SN:chr1	LN:120
read_a	0	chr1	5	60	30M	*	0	0	ACGTACGTACGTACGTACGTACGTACGTAC	IIIIIIIIIIIIIIIIIIIIIIIIIIIIII	NM:i:0	MD:Z:30
read_b	0	chr1	12	60	30M	*	0	0	GTACGTACGTACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIIIIIIIIIIIII	NM:i:1	MD:Z:5A24
read_c	0	chr1	20	60	25M	*	0	0	ACGTACGTACGTACGTACGTACGTA	IIIIIIIIIIIIIIIIIIIIIIIII	NM:i:0	MD:Z:25
read_d	0	chr1	22	60	28M	*	0	0	CGTACGTACGTACGTACGTACGTACGTA	IIIIIIIIIIIIIIIIIIIIIIIIIIII	NM:i:0	MD:Z:28
read_e	16	chr1	35	40	20M	*	0	0	ACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIII	NM:i:0	MD:Z:20
read_f	0	chr1	60	60	30M	*	0	0	TACGTACGTACGTACGTACGTACGTACGTA	IIIIIIIIIIIIIIIIIIIIIIIIIIIIII	NM:i:0	MD:Z:30
read_g	0	chr1	90	60	20M	*	0	0	ACGTACGTACGTACGTACGT	IIIIIIIIIIIIIIIIIIII	NM:i:0	MD:Z:20
`

// TestConsensusIndexedMatchesLinear pins the core invariant of the indexed
// region fast path added for task #45: for a single coordinate-sorted, indexed
// BAM, `consensus -r <region>` taking the BAI seek path (ConsensusFile) must
// emit output byte-identical to the linear walk over the same reads (Consensus,
// which applies the region as a post-filter). The index only changes which BGZF
// blocks are inflated, never which reads reach the consensus engine, so every
// format/mode/region pair must agree exactly. This locks ours-indexed ==
// ours-linear without needing the multi-GB GIAB BAM.
func TestConsensusIndexedMatchesLinear(t *testing.T) {
	bam := makeIndexedBAM(t, writeSAMFile(t, consensusIndexedSAM))

	base := func() ConsensusOptions { return ConsensusOptions{Input: bam} }
	cases := []struct {
		name   string
		mutate func(*ConsensusOptions)
	}{
		{"fasta-simple-sub", func(o *ConsensusOptions) {
			o.Format, o.Mode, o.Regions = ConsensusFASTA, ConsensusModeSimple, []string{"chr1:15-45"}
		}},
		{"fastq-simple-sub", func(o *ConsensusOptions) {
			o.Format, o.Mode, o.Regions = ConsensusFASTQ, ConsensusModeSimple, []string{"chr1:15-45"}
		}},
		{"pileup-simple-sub", func(o *ConsensusOptions) {
			o.Format, o.Mode, o.Regions = ConsensusPileup, ConsensusModeSimple, []string{"chr1:15-45"}
		}},
		{"fasta-bayesian-sub", func(o *ConsensusOptions) {
			o.Format, o.Mode, o.Regions = ConsensusFASTA, ConsensusModeBayesian, []string{"chr1:15-45"}
		}},
		{"pileup-bayesian-sub", func(o *ConsensusOptions) {
			o.Format, o.Mode, o.Regions = ConsensusPileup, ConsensusModeBayesian, []string{"chr1:15-45"}
		}},
		{"fasta-simple-wholechrom", func(o *ConsensusOptions) {
			o.Format, o.Mode, o.Regions = ConsensusFASTA, ConsensusModeSimple, []string{"chr1"}
		}},
		{"pileup-bayesian-wholechrom", func(o *ConsensusOptions) {
			o.Format, o.Mode, o.Regions = ConsensusPileup, ConsensusModeBayesian, []string{"chr1"}
		}},
		{"fasta-simple-deepstart", func(o *ConsensusOptions) {
			// A region that begins past the leading reads exercises the
			// startIdx lower bound skipping reads that ended before the tile.
			o.Format, o.Mode, o.Regions = ConsensusFASTA, ConsensusModeSimple, []string{"chr1:55-100"}
		}},
	}

	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			// Indexed fast path via ConsensusFile.
			idxOpts := base()
			tc.mutate(&idxOpts)
			var got bytes.Buffer
			if err := ConsensusFile(idxOpts, &got, nil); err != nil {
				t.Fatalf("ConsensusFile (indexed): %v", err)
			}

			// Linear oracle: the same options fed through Consensus over the
			// whole file, which applies the region as a post-filter. This is the
			// "ours-linear" reference the fast path must reproduce byte-for-byte.
			f, err := os.Open(bam)
			if err != nil {
				t.Fatal(err)
			}
			linOpts := base()
			tc.mutate(&linOpts)
			var want bytes.Buffer
			if err := Consensus(f, &want, linOpts); err != nil {
				_ = f.Close()
				t.Fatalf("Consensus (linear): %v", err)
			}
			_ = f.Close()

			if got.String() != want.String() {
				t.Errorf("indexed output differs from linear\n--- indexed ---\n%s\n--- linear ---\n%s",
					got.String(), want.String())
			}
		})
	}
}

// TestConsensusFallsBackWithoutIndex confirms that when the sibling .bai/.csi is
// absent, ConsensusFile drops to the unchanged linear path and still produces
// output identical to the linear Consensus walk — the fast-path guard is purely
// an optimisation, never a behaviour change.
func TestConsensusFallsBackWithoutIndex(t *testing.T) {
	bam := makeIndexedBAM(t, writeSAMFile(t, consensusIndexedSAM))
	if err := os.Remove(bam + ".bai"); err != nil {
		t.Fatalf("remove bai: %v", err)
	}

	opts := ConsensusOptions{
		Input:   bam,
		Format:  ConsensusFASTA,
		Mode:    ConsensusModeSimple,
		Regions: []string{"chr1:15-45"},
	}

	var got bytes.Buffer
	if err := ConsensusFile(opts, &got, nil); err != nil {
		t.Fatalf("ConsensusFile (no index): %v", err)
	}

	f, err := os.Open(bam)
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	var want bytes.Buffer
	if err := Consensus(f, &want, opts); err != nil {
		t.Fatalf("Consensus: %v", err)
	}

	if got.String() != want.String() {
		t.Errorf("fallback output differs from linear\n got: %q\nwant: %q", got.String(), want.String())
	}
}

// TestConsensusFileGCNeutral verifies the debug.SetGCPercent tuning inside
// ConsensusFile is output-neutral and correctly scoped: the process GC target
// is restored to its previous value on return, and running the same input
// under different starting GC targets yields byte-identical output. This locks
// the memory-only nature of the GC-headroom change added for the consensus
// perf work (RSS lever only; must never perturb the emitted bytes).
func TestConsensusFileGCNeutral(t *testing.T) {
	bam := makeIndexedBAM(t, writeSAMFile(t, consensusIndexedSAM))
	run := func() string {
		var buf bytes.Buffer
		opts := ConsensusOptions{
			Input:   bam,
			Format:  ConsensusPileup,
			Mode:    ConsensusModeBayesian,
			Regions: []string{"chr1"},
		}
		if err := ConsensusFile(opts, &buf, nil); err != nil {
			t.Fatalf("ConsensusFile: %v", err)
		}
		return buf.String()
	}

	// GC target must be restored after ConsensusFile returns.
	sentinel := debug.SetGCPercent(123)
	_ = run()
	restored := debug.SetGCPercent(sentinel)
	if restored != 123 {
		t.Errorf("ConsensusFile did not restore GC percent: got %d, want 123", restored)
	}

	// Output must be identical regardless of the caller's starting GC target.
	debug.SetGCPercent(400)
	a := run()
	debug.SetGCPercent(10)
	b := run()
	debug.SetGCPercent(sentinel)
	if a != b {
		t.Errorf("ConsensusFile output depends on starting GC target:\n gc=400: %q\n gc=10:  %q", a, b)
	}
}
