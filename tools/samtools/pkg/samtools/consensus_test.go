package samtools

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

// runConsensusOnSAM feeds one SAM-text input through Consensus() with
// the supplied options and returns the emitted text.
func runConsensusOnSAM(t *testing.T, samText string, opts ConsensusOptions) string {
	t.Helper()
	var buf bytes.Buffer
	if err := Consensus(strings.NewReader(samText), &buf, opts); err != nil {
		t.Fatalf("Consensus: %v", err)
	}
	return buf.String()
}

// allMatchSAM has three reads spelling ACGTA at chr1:1-5 with high quals.
const allMatchSAM = `@HD	VN:1.6
@SQ	SN:chr1	LN:8
r1	0	chr1	1	60	5M	*	0	0	ACGTA	IIIII
r2	0	chr1	1	60	5M	*	0	0	ACGTA	IIIII
r3	0	chr1	1	60	5M	*	0	0	ACGTA	IIIII
`

func TestConsensus_FASTA_AllMatch(t *testing.T) {
	out := runConsensusOnSAM(t, allMatchSAM, ConsensusOptions{
		Format: ConsensusFASTA,
	})
	want := ">chr1\nACGTA\n"
	if out != want {
		t.Errorf("FASTA: got %q want %q", out, want)
	}
}

func TestConsensus_FASTQ_AllMatch(t *testing.T) {
	out := runConsensusOnSAM(t, allMatchSAM, ConsensusOptions{
		Format: ConsensusFASTQ,
	})
	// Frequency-only (UseQual=false default): three reads vote for A,
	// freq[A]=3 with weight 8 -> score[A] = 24, tscore = 24, so qual =
	// 100 -> clamped to 93 -> phred byte = 33+93 = 126 = '~'.
	want := "@chr1\nACGTA\n+\n~~~~~\n"
	if out != want {
		t.Errorf("FASTQ: got %q want %q", out, want)
	}
}

func TestConsensus_Pileup_AllMatch(t *testing.T) {
	out := runConsensusOnSAM(t, allMatchSAM, ConsensusOptions{
		Format: ConsensusPileup,
	})
	lines := strings.Split(strings.TrimRight(out, "\n"), "\n")
	if len(lines) != 5 {
		t.Fatalf("want 5 lines, got %d:\n%s", len(lines), out)
	}
	if lines[0] != "chr1\t1\t0\t3\tA\t100\tAAA\tIII" {
		t.Errorf("pileup row 0 = %q", lines[0])
	}
	if lines[4] != "chr1\t5\t0\t3\tA\t100\tAAA\tIII" {
		t.Errorf("pileup row 4 = %q", lines[4])
	}
}

// mixedSAM exercises positions where one read disagrees out of four.
// pos 2: 3xC,1xT  -> usedScore=24, tscore=32, 24/32=0.75 == call_fract default
// pos 4: 3xT,1xA  -> usedScore=24, tscore=32, 24/32=0.75 == call_fract default
const mixedSAM = `@HD	VN:1.6
@SQ	SN:chr1	LN:8
r1	0	chr1	1	60	5M	*	0	0	ACGTA	IIIII
r2	0	chr1	1	60	5M	*	0	0	ACGTA	IIIII
r3	0	chr1	1	60	5M	*	0	0	ACGTA	IIIII
r4	0	chr1	1	60	5M	*	0	0	ATGAA	IIIII
`

func TestConsensus_FASTA_Mixed(t *testing.T) {
	// Default MinCallFraction=0.75. 3/4 = 0.75 -> threshold met
	// (`used_score >= call_fract * tscore` is `>=` upstream).
	out := runConsensusOnSAM(t, mixedSAM, ConsensusOptions{
		Format: ConsensusFASTA,
	})
	want := ">chr1\nACGTA\n"
	if out != want {
		t.Errorf("mixed FASTA: got %q want %q", out, want)
	}
}

func TestConsensus_FASTA_TightCallFraction_BecomesN(t *testing.T) {
	// Bump MinCallFraction past 0.75 so positions 2 and 4 (with 0.75
	// best/total) become N.
	out := runConsensusOnSAM(t, mixedSAM, ConsensusOptions{
		Format:          ConsensusFASTA,
		MinCallFraction: 0.9,
	})
	want := ">chr1\nANGNA\n"
	if out != want {
		t.Errorf("tight call-fract FASTA: got %q want %q", out, want)
	}
}

func TestConsensus_FASTA_AllPositions_FillsN(t *testing.T) {
	// Reads cover only positions 1..5 on an LN=8 contig; -a should
	// pad positions 6..8 with N.
	out := runConsensusOnSAM(t, allMatchSAM, ConsensusOptions{
		Format:       ConsensusFASTA,
		AllPositions: true,
	})
	want := ">chr1\nACGTANNN\n"
	if out != want {
		t.Errorf("-a FASTA: got %q want %q", out, want)
	}
}

func TestConsensus_NoCoverage_NoEmit(t *testing.T) {
	// No reads -> no output (without -a).
	sam := "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:5\n"
	out := runConsensusOnSAM(t, sam, ConsensusOptions{Format: ConsensusFASTA})
	if out != "" {
		t.Errorf("expected empty output, got %q", out)
	}
}

func TestConsensus_NoCoverage_AllPos_EmitsAllN(t *testing.T) {
	sam := "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:4\n"
	// -a alone does not emit fully-uncovered contigs (matches upstream:
	// the *T.out 31/32 fixtures show empty contigs only under -aa).
	out := runConsensusOnSAM(t, sam, ConsensusOptions{
		Format:       ConsensusFASTA,
		AllPositions: true,
	})
	if out != "" {
		t.Errorf("-a empty contig: got %q want empty", out)
	}
	// -aa emits the empty contig as all-N.
	out = runConsensusOnSAM(t, sam, ConsensusOptions{
		Format:       ConsensusFASTA,
		AllPositions: true,
		AllContigs:   true,
	})
	want := ">chr1\nNNNN\n"
	if out != want {
		t.Errorf("-aa empty contig: got %q want %q", out, want)
	}
}

// multiContigSAM covers two contigs.
const multiContigSAM = `@HD	VN:1.6
@SQ	SN:chr1	LN:4
@SQ	SN:chr2	LN:4
r1	0	chr1	1	60	4M	*	0	0	AAAA	IIII
r2	0	chr2	1	60	4M	*	0	0	GGGG	IIII
`

func TestConsensus_FASTA_MultiContig(t *testing.T) {
	out := runConsensusOnSAM(t, multiContigSAM, ConsensusOptions{Format: ConsensusFASTA})
	want := ">chr1\nAAAA\n>chr2\nGGGG\n"
	if out != want {
		t.Errorf("multi-contig: got %q want %q", out, want)
	}
}

// minDepthSAM exercises -d: two reads covering pos 1..3 but only one at pos 4.
const minDepthSAM = `@HD	VN:1.6
@SQ	SN:chr1	LN:4
r1	0	chr1	1	60	4M	*	0	0	ACGT	IIII
r2	0	chr1	1	60	3M	*	0	0	ACG	III
`

func TestConsensus_FASTA_MinDepth_BecomesN(t *testing.T) {
	// With min-depth=2, position 4 (only r1) -> N.
	out := runConsensusOnSAM(t, minDepthSAM, ConsensusOptions{
		Format:   ConsensusFASTA,
		MinDepth: 2,
	})
	want := ">chr1\nACGN\n"
	if out != want {
		t.Errorf("min-depth: got %q want %q", out, want)
	}
}

// delSAM has three reads agreeing on a 1bp deletion at position 3.
const delSAM = `@HD	VN:1.6
@SQ	SN:chr1	LN:5
r1	0	chr1	1	60	2M1D2M	*	0	0	ACTA	IIII
r2	0	chr1	1	60	2M1D2M	*	0	0	ACTA	IIII
r3	0	chr1	1	60	2M1D2M	*	0	0	ACTA	IIII
`

func TestConsensus_FASTA_Deletion_HiddenByDefault(t *testing.T) {
	out := runConsensusOnSAM(t, delSAM, ConsensusOptions{Format: ConsensusFASTA})
	// ShowDel defaults to false -> the deleted column is omitted.
	want := ">chr1\nACTA\n"
	if out != want {
		t.Errorf("del-default: got %q want %q", out, want)
	}
}

func TestConsensus_FASTA_Deletion_ShownWithShowDel(t *testing.T) {
	out := runConsensusOnSAM(t, delSAM, ConsensusOptions{
		Format:  ConsensusFASTA,
		ShowDel: true,
	})
	want := ">chr1\nAC*TA\n"
	if out != want {
		t.Errorf("show-del: got %q want %q", out, want)
	}
}

// TestConsensus_Pileup_ShowDelNo_OmitsStarRow covers reviewer
// correctness finding #5: --show-del no should also drop pileup rows
// whose call is '*' (matches upstream bam_consensus.c:2244).
func TestConsensus_Pileup_ShowDelNo_OmitsStarRow(t *testing.T) {
	out := runConsensusOnSAM(t, delSAM, ConsensusOptions{
		Format: ConsensusPileup,
		// ShowDel deliberately left false (upstream default).
	})
	lines := strings.Split(strings.TrimRight(out, "\n"), "\n")
	// We should see four rows (pos 1, 2, 4, 5) — the deleted column
	// at pos 3 must be suppressed.
	if len(lines) != 4 {
		t.Fatalf("expected 4 lines (no '*' row), got %d:\n%s", len(lines), out)
	}
	for _, ln := range lines {
		if strings.Contains(ln, "\t*\t") {
			t.Errorf("found '*' call row with ShowDel=false: %q", ln)
		}
	}
}

func TestConsensus_Pileup_ShowDelYes_KeepsStarRow(t *testing.T) {
	out := runConsensusOnSAM(t, delSAM, ConsensusOptions{
		Format:  ConsensusPileup,
		ShowDel: true,
	})
	if !strings.Contains(out, "chr1\t3\t0\t3\t*\t") {
		t.Errorf("expected '*' pileup row at pos 3 with ShowDel=true, got:\n%s", out)
	}
}

// insSAM has three reads with a 1bp insertion T between ref positions 2 and 3.
const insSAM = `@HD	VN:1.6
@SQ	SN:chr1	LN:4
r1	0	chr1	1	60	2M1I2M	*	0	0	ACTGA	IIIII
r2	0	chr1	1	60	2M1I2M	*	0	0	ACTGA	IIIII
r3	0	chr1	1	60	2M1I2M	*	0	0	ACTGA	IIIII
`

func TestConsensus_FASTA_Insertion_Included(t *testing.T) {
	out := runConsensusOnSAM(t, insSAM, ConsensusOptions{Format: ConsensusFASTA})
	// Reference has 4 bases (ACGA); insertion T after pos 2 -> ACTGA
	// when ShowIns=true (the default).
	want := ">chr1\nACTGA\n"
	if out != want {
		t.Errorf("insertion-include: got %q want %q", out, want)
	}
}

func TestConsensus_FASTA_Insertion_Suppressed(t *testing.T) {
	out := runConsensusOnSAM(t, insSAM, ConsensusOptions{
		Format:    ConsensusFASTA,
		NoShowIns: true,
	})
	want := ">chr1\nACGA\n"
	if out != want {
		t.Errorf("insertion-suppress: got %q want %q", out, want)
	}
}

func TestConsensus_FASTA_Insertion_MarkIns(t *testing.T) {
	out := runConsensusOnSAM(t, insSAM, ConsensusOptions{
		Format:  ConsensusFASTA,
		MarkIns: true,
	})
	// MarkIns prepends '_' before each inserted column, matching upstream
	// bam_consensus.c:2409-2412 (the marker byte is '_', not '+').
	want := ">chr1\nAC_TGA\n"
	if out != want {
		t.Errorf("mark-ins: got %q want %q", out, want)
	}
}

func TestConsensus_LineLen_Wrapping(t *testing.T) {
	out := runConsensusOnSAM(t, allMatchSAM, ConsensusOptions{
		Format:  ConsensusFASTA,
		LineLen: 2,
	})
	want := ">chr1\nAC\nGT\nA\n"
	if out != want {
		t.Errorf("line-len wrap: got %q want %q", out, want)
	}
}

func TestConsensus_Pileup_AllPositions_ZeroDepthRows(t *testing.T) {
	// allMatchSAM covers 1..5 on LN=8; -a should emit rows 6..8 as
	// zero depth, calling N with cq=0 and "*" for seq/qual.
	out := runConsensusOnSAM(t, allMatchSAM, ConsensusOptions{
		Format:       ConsensusPileup,
		AllPositions: true,
	})
	lines := strings.Split(strings.TrimRight(out, "\n"), "\n")
	if len(lines) != 8 {
		t.Fatalf("want 8 lines, got %d:\n%s", len(lines), out)
	}
	if lines[7] != "chr1\t8\t0\t0\tN\t0\t*\t*" {
		t.Errorf("zero-depth row: got %q", lines[7])
	}
}

// TestConsensus_AmbigCodes_50_30_20_LandsOnM is the canonical reviewer
// correctness fixture (finding #1): a 50/30/20 mix of A/C/G with -A
// (ambig) MUST land on 'M' (the A+C IUPAC code), NOT 'N'.
//
//	freq A=5, C=3, G=2 -> score[A]=40, score[C]=24, score[G]=16
//	tscore = 80
//	s1=40 (A), s2=24 (C)
//	het: s2 >= 0.5*s1 -> 24 >= 20 -> TRUE -> used = 1|2 = 3 -> 'M'
//	usedScore = 40+24 = 64; call_fract gate: 64 >= 0.75*80 = 60 -> PASS.
//
// Without the bogus "min fraction on dominant base alone" gate (the
// reviewer finding #1 said NOT to implement) we land on M; with it
// we'd see N because A alone is 40/80 = 0.5 < 0.6.
func TestConsensus_AmbigCodes_50_30_20_LandsOnM(t *testing.T) {
	var b strings.Builder
	b.WriteString("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1\n")
	for i := 0; i < 5; i++ {
		b.WriteString("r")
		b.WriteByte(byte('A' + i))
		b.WriteString("\t0\tchr1\t1\t60\t1M\t*\t0\t0\tA\tI\n")
	}
	for i := 0; i < 3; i++ {
		b.WriteString("c")
		b.WriteByte(byte('A' + i))
		b.WriteString("\t0\tchr1\t1\t60\t1M\t*\t0\t0\tC\tI\n")
	}
	for i := 0; i < 2; i++ {
		b.WriteString("g")
		b.WriteByte(byte('A' + i))
		b.WriteString("\t0\tchr1\t1\t60\t1M\t*\t0\t0\tG\tI\n")
	}
	out := runConsensusOnSAM(t, b.String(), ConsensusOptions{
		Format:     ConsensusFASTA,
		AmbigCodes: true,
	})
	// Must contain 'M', must NOT contain 'N'.
	if !strings.Contains(out, "M\n") {
		t.Errorf("50/30/20 + ambig must land on M, got %q", out)
	}
	if strings.Contains(out, "N\n") {
		t.Errorf("50/30/20 + ambig must NOT be N, got %q", out)
	}
}

// TestConsensus_AmbigCodes_Het_Simple keeps the simpler 50/50 fixture
// for completeness.
func TestConsensus_AmbigCodes_Het_Simple(t *testing.T) {
	sam := `@HD	VN:1.6
@SQ	SN:chr1	LN:1
r1	0	chr1	1	60	1M	*	0	0	A	I
r2	0	chr1	1	60	1M	*	0	0	A	I
r3	0	chr1	1	60	1M	*	0	0	C	I
r4	0	chr1	1	60	1M	*	0	0	C	I
`
	out := runConsensusOnSAM(t, sam, ConsensusOptions{
		Format:     ConsensusFASTA,
		AmbigCodes: true,
	})
	if !strings.Contains(out, "M\n") {
		t.Errorf("expected M ambig code, got %q", out)
	}
}

// TestConsensus_FrequencyOnlyByDefault asserts that with UseQual=false
// (upstream default), the score doesn't depend on per-base quality —
// a low-Q read is worth the same as a high-Q read (reviewer finding
// #3).
func TestConsensus_FrequencyOnlyByDefault(t *testing.T) {
	// Same fixture twice, one with all I (qual 40), one with all !
	// (qual 0). With UseQual=false the result must be identical.
	high := `@HD	VN:1.6
@SQ	SN:chr1	LN:3
r1	0	chr1	1	60	3M	*	0	0	ACG	III
r2	0	chr1	1	60	3M	*	0	0	ACG	III
`
	low := `@HD	VN:1.6
@SQ	SN:chr1	LN:3
r1	0	chr1	1	60	3M	*	0	0	ACG	!!!
r2	0	chr1	1	60	3M	*	0	0	ACG	!!!
`
	hiOut := runConsensusOnSAM(t, high, ConsensusOptions{Format: ConsensusFASTA})
	loOut := runConsensusOnSAM(t, low, ConsensusOptions{Format: ConsensusFASTA})
	if hiOut != loOut {
		t.Errorf("frequency-only mode (UseQual=false) should ignore quality:\nhigh-Q: %q\nlow-Q:  %q", hiOut, loOut)
	}
	// With FASTQ output the call-confidence must also be the same
	// in both cases since tscore/usedScore are independent of qual.
	hiOut = runConsensusOnSAM(t, high, ConsensusOptions{Format: ConsensusFASTQ})
	loOut = runConsensusOnSAM(t, low, ConsensusOptions{Format: ConsensusFASTQ})
	if hiOut != loOut {
		t.Errorf("FASTQ frequency-only mode mismatch:\nhigh-Q: %q\nlow-Q:  %q", hiOut, loOut)
	}
}

// TestConsensus_UseQualOn_PreferHighQ asserts that with UseQual=true,
// per-base quality DOES bias the call — a single high-Q read beats
// three low-Q reads on quality-weighted score. This tightens reviewer
// finding #3: the multiplier is gated on UseQual, not always-on.
func TestConsensus_UseQualOn_PreferHighQ(t *testing.T) {
	// One A at qual 40 (I) + three C at qual 1 (").
	// Frequency-only: scoreA=8, scoreC=24, tscore=32 -> 24/32=0.75
	// passes call-fract (>=0.75) -> call C.
	// UseQual on:     scoreA=40*8=320, scoreC=3*1*8=24,
	// tscore=344, s1=320 (A) -> 320/344=0.93 -> call A.
	sam := `@HD	VN:1.6
@SQ	SN:chr1	LN:1
r1	0	chr1	1	60	1M	*	0	0	A	I
r2	0	chr1	1	60	1M	*	0	0	C	"
r3	0	chr1	1	60	1M	*	0	0	C	"
r4	0	chr1	1	60	1M	*	0	0	C	"
`
	freqOut := runConsensusOnSAM(t, sam, ConsensusOptions{Format: ConsensusFASTA})
	qualOut := runConsensusOnSAM(t, sam, ConsensusOptions{
		Format:  ConsensusFASTA,
		UseQual: true,
	})
	if !strings.Contains(freqOut, "C\n") {
		t.Errorf("frequency-only should pick C (3 votes), got %q", freqOut)
	}
	if !strings.Contains(qualOut, "A\n") {
		t.Errorf("UseQual=true should pick A (40*8 > 3*1*8), got %q", qualOut)
	}
}

// TestConsensus_BayesianDefault_FromFile confirms the default
// invocation (Mode=Bayesian) runs the real Gap5 bayesian caller with no
// fallback warning on stderr.
func TestConsensus_BayesianDefault_FromFile(t *testing.T) {
	tmpDir := t.TempDir()
	bamPath := tmpDir + "/in.sam"
	if err := writeStringFile(bamPath, allMatchSAM); err != nil {
		t.Fatalf("write tmp SAM: %v", err)
	}
	var sout, serr bytes.Buffer
	opts := ConsensusOptions{
		Input:  bamPath,
		Format: ConsensusFASTA,
		Mode:   ConsensusModeBayesian,
	}
	if err := ConsensusFile(opts, &sout, &serr); err != nil {
		t.Fatalf("ConsensusFile(bayesian): %v", err)
	}
	if serr.Len() != 0 {
		t.Errorf("expected no stderr warning, got %q", serr.String())
	}
	if !strings.Contains(sout.String(), "ACGTA") {
		t.Errorf("bayesian default should emit a consensus, got %q", sout.String())
	}
}

// TestConsensus_BayesianDefault_LibraryCall confirms Consensus() with
// Mode=Bayesian runs the bayesian caller directly.
func TestConsensus_BayesianDefault_LibraryCall(t *testing.T) {
	var sout bytes.Buffer
	opts := ConsensusOptions{
		Format: ConsensusFASTA,
		Mode:   ConsensusModeBayesian,
	}
	if err := Consensus(strings.NewReader(allMatchSAM), &sout, opts); err != nil {
		t.Fatalf("Consensus(bayesian): %v", err)
	}
	if !strings.Contains(sout.String(), "ACGTA") {
		t.Errorf("bayesian default should emit a consensus, got %q", sout.String())
	}
}

func TestParseConsensusFormat(t *testing.T) {
	cases := []struct {
		in   string
		want ConsensusFormat
		err  bool
	}{
		{"", ConsensusFASTA, false},
		{"fasta", ConsensusFASTA, false},
		{"FASTA", ConsensusFASTA, false},
		{"fa", ConsensusFASTA, false},
		{"fastq", ConsensusFASTQ, false},
		{"FQ", ConsensusFASTQ, false},
		{"pileup", ConsensusPileup, false},
		{"junk", 0, true},
	}
	for _, c := range cases {
		got, err := ParseConsensusFormat(c.in)
		if (err != nil) != c.err {
			t.Errorf("ParseConsensusFormat(%q) err=%v want=%v", c.in, err, c.err)
			continue
		}
		if err == nil && got != c.want {
			t.Errorf("ParseConsensusFormat(%q) = %d, want %d", c.in, got, c.want)
		}
	}
}

func TestParseConsensusMode(t *testing.T) {
	cases := []struct {
		in   string
		want ConsensusMode
		err  bool
	}{
		{"", ConsensusModeSimple, false},
		{"simple", ConsensusModeSimple, false},
		{"SIMPLE", ConsensusModeSimple, false},
		{"bayesian", ConsensusModeBayesian, false},
		{"bayesian_r", ConsensusModeBayesian, false},
		{"bayesian_m", ConsensusModeBayesian, false},
		{"bayesian_p", ConsensusModeBayesian, false},
		{"bayesian_116", ConsensusModeBayesian, false},
		{"junk", 0, true},
	}
	for _, c := range cases {
		got, err := ParseConsensusMode(c.in)
		if (err != nil) != c.err {
			t.Errorf("ParseConsensusMode(%q) err=%v want=%v", c.in, err, c.err)
			continue
		}
		if err == nil && got != c.want {
			t.Errorf("ParseConsensusMode(%q) = %d, want %d", c.in, got, c.want)
		}
	}
}

// TestConsensus_BayesianUpstreamParity runs BOTH the live upstream
// `samtools consensus` binary and the Go port on the same vendored fixture
// and asserts the bayesian output is identical, replacing the former
// golden-file comparison. The upstream binary is built on demand and a
// build failure is fatal, never skipped. Each case's `args` is the exact
// invocation from reference_code/samtools/test/consensus/consensus.reg.
// Fixtures requiring a reference FASTA (-T, the *T.out goldens) are out of
// scope for v1 and are not exercised here; see docs/PARITY_ROADMAP.md.
func TestConsensus_BayesianUpstreamParity(t *testing.T) {
	bin := upstreamSamtools(t)
	const fixDir = "../../../../reference_code/samtools/test/consensus"
	cases := []struct {
		name  string
		input string
		args  []string // upstream `samtools consensus` args (sans input file)
		opts  ConsensusOptions
	}{
		{"18q", "consen1.sam", []string{"-f", "fastq", "--no-use-MQ", "-C", "0", "-m", "bayesian"}, consBayes(ConsensusFASTQ, 0, false, false, "", "")},
		{"19q", "consen1.sam", []string{"-f", "fastq", "--no-use-MQ", "-C", "19", "-m", "bayesian"}, consBayes(ConsensusFASTQ, 19, false, false, "", "")},
		{"18p", "consen1.sam", []string{"-f", "pileup", "--no-use-MQ", "-C", "0", "-m", "bayesian"}, consBayes(ConsensusPileup, 0, false, false, "", "")},
		{"19p", "consen1.sam", []string{"-f", "pileup", "--no-use-MQ", "-C", "19", "-m", "bayesian"}, consBayes(ConsensusPileup, 19, false, false, "", "")},
		{"20p", "consen1.sam", []string{"-f", "pileup", "--no-use-MQ", "-C", "30", "-A", "-m", "bayesian"}, consBayes(ConsensusPileup, 30, false, true, "", "")},
		{"21p", "consen1.sam", []string{"-f", "pileup", "--no-use-MQ", "-C", "31", "-A", "-m", "bayesian"}, consBayes(ConsensusPileup, 31, false, true, "", "")},
		{"30", "consen1c.sam", []string{"-f", "fastq", "--show-del", "yes", "--show-ins", "no", "-m", "bayesian", "-C0"}, consBayes(ConsensusFASTQ, 0, true, false, "yes", "no")},
		{"31", "consen1c.sam", []string{"-a", "-f", "fastq", "--show-del", "yes", "--show-ins", "no", "-m", "bayesian", "-C0"}, consBayesA(ConsensusFASTQ, 0, true, false, "yes", "no", 1)},
		{"32", "consen1c.sam", []string{"-aa", "-f", "fastq", "--show-del", "yes", "--show-ins", "no", "-m", "bayesian", "-C0"}, consBayesA(ConsensusFASTQ, 0, true, false, "yes", "no", 2)},
		{"40", "consen1c.sam", []string{"-f", "pileup", "--show-del", "yes", "--show-ins", "no", "-m", "bayesian", "-C0"}, consBayes(ConsensusPileup, 0, true, false, "yes", "no")},
		{"41", "consen1c.sam", []string{"-a", "-f", "pileup", "--show-del", "yes", "--show-ins", "no", "-m", "bayesian", "-C0"}, consBayesA(ConsensusPileup, 0, true, false, "yes", "no", 1)},
		{"42", "consen1c.sam", []string{"-aa", "-f", "pileup", "--show-del", "yes", "--show-ins", "no", "-m", "bayesian", "-C0"}, consBayesA(ConsensusPileup, 0, true, false, "yes", "no", 2)},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			input := filepath.Join(fixDir, c.input)
			samBytes, err := os.ReadFile(input)
			if err != nil {
				t.Fatalf("read input %s: %v", c.input, err)
			}

			// Live upstream invocation.
			cmd := exec.Command(bin, append(append([]string{"consensus"}, c.args...), input)...)
			var upOut, upErr bytes.Buffer
			cmd.Stdout = &upOut
			cmd.Stderr = &upErr
			if err := cmd.Run(); err != nil {
				t.Fatalf("upstream samtools consensus %v: %v\n%s", c.args, err, upErr.String())
			}

			// Go port.
			got := runConsensusOnSAM(t, string(samBytes), c.opts)
			if got != upOut.String() {
				t.Errorf("%s mismatch:\n--- got (go) ---\n%s\n--- want (upstream) ---\n%s", c.name, got, upOut.String())
			}
		})
	}
}

// consBayes builds a bayesian-mode ConsensusOptions for the parity test.
func consBayes(f ConsensusFormat, cutoff int, useMQ, ambig bool, showDel, showIns string) ConsensusOptions {
	return consBayesA(f, cutoff, useMQ, ambig, showDel, showIns, 0)
}

// consBayesA is consBayes with an explicit -a level (0/1/2).
func consBayesA(f ConsensusFormat, cutoff int, useMQ, ambig bool, showDel, showIns string, allLevel int) ConsensusOptions {
	o := ConsensusOptions{
		Format:        f,
		Mode:          ConsensusModeBayesian,
		ConsCutoff:    cutoff,
		ConsCutoffSet: true,
		AmbigCodes:    ambig,
		UseMQual:      useMQ,
		UseMQualSet:   true,
		AllPositions:  allLevel >= 1,
		AllContigs:    allLevel >= 2,
	}
	if showDel == "yes" {
		o.ShowDel = true
	}
	if showIns == "no" {
		o.NoShowIns = true
	}
	return o
}

// Table-driven smoke test across all three formats sharing a fixture.
func TestConsensus_FormatTable(t *testing.T) {
	cases := []struct {
		name string
		fmt  ConsensusFormat
		want string
	}{
		{
			name: "fasta",
			fmt:  ConsensusFASTA,
			want: ">chr1\nACGTA\n",
		},
		{
			name: "fastq",
			fmt:  ConsensusFASTQ,
			want: "@chr1\nACGTA\n+\n~~~~~\n",
		},
		{
			name: "pileup",
			fmt:  ConsensusPileup,
			want: "chr1\t1\t0\t3\tA\t100\tAAA\tIII\nchr1\t2\t0\t3\tC\t100\tCCC\tIII\nchr1\t3\t0\t3\tG\t100\tGGG\tIII\nchr1\t4\t0\t3\tT\t100\tTTT\tIII\nchr1\t5\t0\t3\tA\t100\tAAA\tIII\n",
		},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			out := runConsensusOnSAM(t, allMatchSAM, ConsensusOptions{Format: c.fmt})
			if out != c.want {
				t.Errorf("got %q\nwant %q", out, c.want)
			}
		})
	}
}

// writeStringFile is a tiny helper for the tmpdir SAM fixture used by
// the bayesian-fallback test (which needs ConsensusFile, which opens
// a path on disk).
func writeStringFile(path, body string) error {
	return os.WriteFile(path, []byte(body), 0o600)
}

// hetOnlySAM has four reads spanning chr1:1-3 with a heterozygous
// flanking pair (A/C at pos1 and pos3) and a homozygous middle (G at
// pos2). It deliberately mixes het and homozygous positions so that
// --het-only's suppression of the homozygous column is visible.
const hetOnlySAM = `@HD	VN:1.6
@SQ	SN:chr1	LN:3
r1	0	chr1	1	60	3M	*	0	0	AGA	III
r2	0	chr1	1	60	3M	*	0	0	AGA	III
r3	0	chr1	1	60	3M	*	0	0	CGC	III
r4	0	chr1	1	60	3M	*	0	0	CGC	III
`

// allHomozSAM has three identical reads spelling ACGTA — every position
// is homozygous, so --het-only must suppress all of them.
const allHomozSAM = `@HD	VN:1.6
@SQ	SN:chr1	LN:5
r1	0	chr1	1	60	5M	*	0	0	ACGTA	IIIII
r2	0	chr1	1	60	5M	*	0	0	ACGTA	IIIII
r3	0	chr1	1	60	5M	*	0	0	ACGTA	IIIII
`

// TestConsensus_HetOnly_SuppressesHomozygous covers the implemented
// --het-only behaviour across both calling modes (simple + bayesian) and
// both the --ambig and non-ambig paths. The homozygous middle position
// (pos2, all-G) must be suppressed: rendered 'N' in FASTA (coordinates
// preserved) and omitted entirely in pileup. The flanking heterozygous
// positions (pos1, pos3) survive.
func TestConsensus_HetOnly_SuppressesHomozygous(t *testing.T) {
	for _, mode := range []ConsensusMode{ConsensusModeSimple, ConsensusModeBayesian} {
		for _, ambig := range []bool{false, true} {
			modeName := "simple"
			if mode == ConsensusModeBayesian {
				modeName = "bayesian"
			}
			t.Run(modeName+"_ambig", func(t *testing.T) {
				// FASTA: homozygous middle becomes N, flanks kept.
				faOpts := ConsensusOptions{Format: ConsensusFASTA, Mode: mode, AmbigCodes: ambig, HetOnly: true}
				fa := runConsensusOnSAM(t, hetOnlySAM, faOpts)
				// The middle character must always be 'N' (homozygous
				// suppressed). With --ambig the flanks are the het code
				// 'M'; without --ambig the flanks cannot be represented
				// as a single base and fall to 'N', but they are still
				// emitted (not deleted), so the record stays length 3.
				wantFA := ">chr1\nMNM\n"
				if !ambig {
					wantFA = ">chr1\nNNN\n"
				}
				if fa != wantFA {
					t.Errorf("FASTA het-only (ambig=%v): got %q want %q", ambig, fa, wantFA)
				}

				// Pileup: homozygous middle row omitted; only pos1 and
				// pos3 rows remain.
				puOpts := faOpts
				puOpts.Format = ConsensusPileup
				pu := runConsensusOnSAM(t, hetOnlySAM, puOpts)
				if strings.Contains(pu, "\t2\t") {
					t.Errorf("pileup het-only kept homozygous pos2 row: %q", pu)
				}
				lines := strings.Split(strings.TrimRight(pu, "\n"), "\n")
				if len(lines) != 2 {
					t.Fatalf("pileup het-only: want 2 rows (pos1,pos3), got %d: %q", len(lines), pu)
				}
				for _, l := range lines {
					if !strings.HasPrefix(l, "chr1\t1\t") && !strings.HasPrefix(l, "chr1\t3\t") {
						t.Errorf("pileup het-only: unexpected row %q", l)
					}
				}
			})
		}
	}
}

// TestConsensus_HetOnly_AllHomozygous verifies that an input with no
// heterozygous positions yields empty FASTA/pileup output (the covered
// span trims away) and an all-'N' record under -a, for both calling
// modes.
func TestConsensus_HetOnly_AllHomozygous(t *testing.T) {
	for _, mode := range []ConsensusMode{ConsensusModeSimple, ConsensusModeBayesian} {
		fa := runConsensusOnSAM(t, allHomozSAM, ConsensusOptions{
			Format: ConsensusFASTA, Mode: mode, AmbigCodes: true, HetOnly: true})
		if fa != "" {
			t.Errorf("mode %v: all-homozygous FASTA het-only: want empty, got %q", mode, fa)
		}
		pu := runConsensusOnSAM(t, allHomozSAM, ConsensusOptions{
			Format: ConsensusPileup, Mode: mode, AmbigCodes: true, HetOnly: true})
		if pu != "" {
			t.Errorf("mode %v: all-homozygous pileup het-only: want empty, got %q", mode, pu)
		}
		// -a forces full-length emission: every position masked to 'N'.
		faAll := runConsensusOnSAM(t, allHomozSAM, ConsensusOptions{
			Format: ConsensusFASTA, Mode: mode, AmbigCodes: true, HetOnly: true, AllPositions: true})
		if faAll != ">chr1\nNNNNN\n" {
			t.Errorf("mode %v: all-homozygous FASTA het-only -a: got %q want %q",
				mode, faAll, ">chr1\nNNNNN\n")
		}
	}
}

// TestConsensus_HetOnly_OffIsUnaffected sanity-checks that without
// HetOnly the homozygous middle position is present (so the suppression
// tests above are exercising a real difference, not a degenerate input).
func TestConsensus_HetOnly_OffIsUnaffected(t *testing.T) {
	pu := runConsensusOnSAM(t, hetOnlySAM, ConsensusOptions{
		Format: ConsensusPileup, AmbigCodes: true})
	if !strings.Contains(pu, "chr1\t2\t") {
		t.Errorf("without het-only, expected homozygous pos2 row present: %q", pu)
	}
}
