package samtools

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"testing"
)

// consensusRefFixtureSAM covers chr1:5-9 only (leading gap 1-4, trailing gap
// 10-20) and leaves chr2 wholly uncovered, so the -T substitution exercises
// leading, trailing and whole-contig gap fills.
const consensusRefFixtureSAM = `@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:20
@SQ	SN:chr2	LN:20
r1	0	chr1	5	60	5M	*	0	0	ACGTA	IIIII
r2	0	chr1	5	60	5M	*	0	0	ACGTA	IIIII
`

// consensusRefFixtureGapSAM covers chr1:1-3 and chr1:8-10, leaving an INTERNAL
// gap (4-7) that the -T substitution must fill from the reference.
const consensusRefFixtureGapSAM = `@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:20
@SQ	SN:chr2	LN:20
r1	0	chr1	1	60	3M	*	0	0	ACG	III
r3	0	chr1	8	60	3M	*	0	0	TAC	III
`

// consensusRefFASTA is the reference used for the -T substitution. Its bases
// are distinct from the read bases at the gap positions so a missing
// substitution is detectable.
const consensusRefFASTA = ">chr1\nACGTACGTACGTACGTACGT\n>chr2\nTTTTGGGGCCCCAAAATTTT\n"

// writeConsensusRefFASTA writes the reference and builds its .fai via the
// upstream samtools binary (so the .fai matches what upstream consensus reads).
func writeConsensusRefFASTA(t *testing.T, bin string) string {
	t.Helper()
	dir := t.TempDir()
	refPath := filepath.Join(dir, "ref.fa")
	if err := os.WriteFile(refPath, []byte(consensusRefFASTA), 0o600); err != nil {
		t.Fatalf("write ref: %v", err)
	}
	if out, err := exec.Command(bin, "faidx", refPath).CombinedOutput(); err != nil {
		t.Fatalf("samtools faidx: %v\n%s", err, out)
	}
	return refPath
}

// TestConsensusReferenceUpstreamParity is the LIVE parity check for
// `samtools consensus -T/--reference`. It builds a sorted+indexed BAM and a
// reference FASTA, then compares upstream `consensus ... -T ref.fa` against the
// Go port byte-for-byte across the pileup/FASTA/FASTQ formats, the simple and
// bayesian modes, the -a / -aa fills, an internal gap, and a non-default
// --ref-qual. The upstream binary is built on demand; a build failure is
// fatal, never skipped.
func TestConsensusReferenceUpstreamParity(t *testing.T) {
	bin := upstreamSamtools(t)
	refPath := writeConsensusRefFASTA(t, bin)
	bamPath := buildSortedIndexedBAM(t, bin, consensusRefFixtureSAM)
	gapBamPath := buildSortedIndexedBAM(t, bin, consensusRefFixtureGapSAM)

	// Each case fixes the consensus mode explicitly (upstream's default is
	// bayesian; the Go ConsensusOptions zero value is simple), so the CLI mode
	// string and the Go Mode always agree.
	cases := []struct {
		name    string
		bam     string
		modeCLI string
		mode    ConsensusMode
		cliArgs []string
		apply   func(*ConsensusOptions)
	}{
		{"pileup_a", bamPath, "simple", ConsensusModeSimple, []string{"-a", "--format", "pileup"},
			func(o *ConsensusOptions) { o.AllPositions = true; o.Format = ConsensusPileup }},
		{"pileup_aa", bamPath, "simple", ConsensusModeSimple, []string{"-aa", "--format", "pileup"},
			func(o *ConsensusOptions) { o.AllPositions = true; o.AllContigs = true; o.Format = ConsensusPileup }},
		{"fasta_default", bamPath, "simple", ConsensusModeSimple, []string{"-f", "fasta"},
			func(o *ConsensusOptions) { o.Format = ConsensusFASTA }},
		{"fasta_a", bamPath, "simple", ConsensusModeSimple, []string{"-a", "-f", "fasta"},
			func(o *ConsensusOptions) { o.AllPositions = true; o.Format = ConsensusFASTA }},
		{"fastq_a", bamPath, "simple", ConsensusModeSimple, []string{"-a", "-f", "fastq"},
			func(o *ConsensusOptions) { o.AllPositions = true; o.Format = ConsensusFASTQ }},
		{"fastq_a_refqual", bamPath, "simple", ConsensusModeSimple, []string{"-a", "-f", "fastq", "--ref-qual", "20"},
			func(o *ConsensusOptions) { o.AllPositions = true; o.Format = ConsensusFASTQ; o.RefQual = 20 }},
		{"fasta_internal_gap", gapBamPath, "simple", ConsensusModeSimple, []string{"-f", "fasta"},
			func(o *ConsensusOptions) { o.Format = ConsensusFASTA }},
		{"fasta_aa_empty_contig", gapBamPath, "simple", ConsensusModeSimple, []string{"-aa", "-f", "fasta"},
			func(o *ConsensusOptions) { o.AllPositions = true; o.AllContigs = true; o.Format = ConsensusFASTA }},
		{"pileup_a_bayesian", bamPath, "bayesian", ConsensusModeBayesian, []string{"-a", "--format", "pileup"},
			func(o *ConsensusOptions) { o.AllPositions = true; o.Format = ConsensusPileup }},
		{"fastq_a_bayesian", bamPath, "bayesian", ConsensusModeBayesian, []string{"-a", "-f", "fastq"},
			func(o *ConsensusOptions) { o.AllPositions = true; o.Format = ConsensusFASTQ }},
	}

	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			args := append([]string{"consensus", "-m", tc.modeCLI}, tc.cliArgs...)
			args = append(args, "-T", refPath, tc.bam)
			upCmd := exec.Command(bin, args...)
			var upOut, upErr bytes.Buffer
			upCmd.Stdout = &upOut
			upCmd.Stderr = &upErr
			if err := upCmd.Run(); err != nil {
				t.Fatalf("upstream consensus %v: %v\n%s", tc.cliArgs, err, upErr.String())
			}

			opts := ConsensusOptions{Reference: refPath, Mode: tc.mode}
			tc.apply(&opts)
			in, err := os.Open(tc.bam)
			if err != nil {
				t.Fatalf("open bam: %v", err)
			}
			defer in.Close()
			var goOut bytes.Buffer
			if err := Consensus(in, &goOut, opts); err != nil {
				t.Fatalf("Consensus: %v", err)
			}

			if goOut.String() != upOut.String() {
				t.Fatalf("consensus -T parity differs (%s)\n--- upstream ---\n%s\n--- go ---\n%s",
					tc.name, upOut.String(), goOut.String())
			}
		})
	}
}
