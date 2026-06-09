package samtools

import (
	"os"
	"os/exec"
	"path/filepath"
	"testing"
)

// runUpstreamConsensus invokes the upstream samtools `consensus`
// subcommand on samPath and returns its stdout. samtools writes
// progress/warnings to stderr, which is discarded. The upstream binary is
// built once on demand via the shared upstreamSamtools helper
// (upstream_test.go); a genuine build failure is fatal, never skipped.
func runUpstreamConsensus(t *testing.T, bin, samPath string, args ...string) string {
	t.Helper()
	full := append([]string{"consensus"}, args...)
	full = append(full, samPath)
	cmd := exec.Command(bin, full...)
	out, err := cmd.Output()
	if err != nil {
		t.Fatalf("upstream samtools %v failed: %v", full, err)
	}
	return string(out)
}

// upstreamBugSAM mixes heterozygous and homozygous positions: four reads
// over chr1:1-3 give a het flank (A/C at pos1 and pos3) and a homozygous
// middle (all-G at pos2). Because it contains BOTH kinds of position, any
// genuine --het-only filter (like ours) visibly changes the output by
// dropping the homozygous column — which is exactly what lets the test
// below distinguish "filter applied" from "flag ignored".
const upstreamBugSAM = `@HD	VN:1.6
@SQ	SN:chr1	LN:3
r1	0	chr1	1	60	3M	*	0	0	AGA	III
r2	0	chr1	1	60	3M	*	0	0	AGA	III
r3	0	chr1	1	60	3M	*	0	0	CGC	III
r4	0	chr1	1	60	3M	*	0	0	CGC	III
`

// TestConsensus_HetOnlyUpstreamBug is the LIVE check that documents the
// upstream dead-option bug and our intentional, correct divergence from
// it. It builds (or reuses) the vendored upstream samtools binary and
// proves two facts by execution:
//
//  1. UPSTREAM IGNORES --het-only (the bug). samtools parses --het-only
//     into opts.het_only (bam_consensus.c) but never reads that variable
//     anywhere in the calling/output path, so its consensus output is
//     byte-for-byte IDENTICAL with and without the flag — even though the
//     fixture contains a homozygous position the flag's name says should
//     be suppressed. Documented in docs/UPSTREAM_BUGS.md.
//  2. OUR Go --het-only DIFFERS from our own output without the flag: it
//     suppresses the homozygous column (rendered 'N' in FASTA, omitted in
//     pileup), which is the behaviour the flag implies. This is the
//     intentional divergence — we fix upstream's dead option.
//
// No committed golden/snapshot file is used: upstream and Go outputs are
// produced live and compared directly. If the upstream binary cannot be
// built, the test fails hard (it never skips), per the project rules.
func TestConsensus_HetOnlyUpstreamBug(t *testing.T) {
	if testing.Short() {
		t.Skip("skipping live upstream build/parity test in -short mode")
	}
	bin := upstreamSamtools(t)

	// Write the fixture to a tmp SAM file for the upstream binary.
	samPath := filepath.Join(t.TempDir(), "bug.sam")
	if err := os.WriteFile(samPath, []byte(upstreamBugSAM), 0o600); err != nil {
		t.Fatalf("write fixture: %v", err)
	}

	cases := []struct {
		name      string
		cliFormat string
		goFormat  ConsensusFormat
	}{
		{"fasta", "fasta", ConsensusFASTA},
		{"pileup", "pileup", ConsensusPileup},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			// Fact 1: upstream --het-only is a no-op (the bug). Output
			// is identical with and without the flag, despite the
			// fixture's homozygous pos2 that a real filter would drop.
			// We run with -A/--ambig so heterozygous positions are
			// actually emitted as IUPAC codes; that makes the homozygous
			// vs het distinction visible and guarantees a genuine filter
			// would change the output.
			upNoFlag := runUpstreamConsensus(t, bin, samPath,
				"-m", "simple", "-A", "-f", c.cliFormat)
			upWithFlag := runUpstreamConsensus(t, bin, samPath,
				"-m", "simple", "-A", "-f", c.cliFormat, "--het-only")
			if upNoFlag != upWithFlag {
				t.Fatalf("upstream --het-only changed %s output, expected the "+
					"documented no-op (bug):\n--- without ---\n%s\n--- with ---\n%s",
					c.name, upNoFlag, upWithFlag)
			}

			// Fact 2: OUR --het-only differs from our own no-flag output
			// — we suppress the homozygous column. This is the
			// intentional divergence that fixes the upstream dead option.
			goNoFlag := runConsensusOnSAM(t, upstreamBugSAM,
				ConsensusOptions{Format: c.goFormat, AmbigCodes: true})
			goWithFlag := runConsensusOnSAM(t, upstreamBugSAM,
				ConsensusOptions{Format: c.goFormat, AmbigCodes: true, HetOnly: true})
			if goNoFlag == goWithFlag {
				t.Fatalf("Go --het-only %s output did not change (expected the "+
					"homozygous column to be suppressed):\n--- without ---\n%s",
					c.name, goNoFlag)
			}
			// And our no-flag output must still match upstream's no-flag
			// output (we only diverge WHEN the flag is set), confirming we
			// are not accidentally changing the baseline consensus.
			if goNoFlag != upNoFlag {
				t.Fatalf("Go baseline (no --het-only) %s output differs from "+
					"upstream:\n--- go ---\n%s\n--- upstream ---\n%s",
					c.name, goNoFlag, upNoFlag)
			}
		})
	}
}
