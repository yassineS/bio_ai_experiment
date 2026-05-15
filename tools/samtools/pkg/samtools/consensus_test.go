package samtools

import (
	"bytes"
	"os"
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
	out := runConsensusOnSAM(t, sam, ConsensusOptions{
		Format:       ConsensusFASTA,
		AllPositions: true,
	})
	want := ">chr1\nNNNN\n"
	if out != want {
		t.Errorf("-a empty contig: got %q want %q", out, want)
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
	// MarkIns prepends '+' before the inserted base.
	want := ">chr1\nAC+TGA\n"
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

// TestConsensus_BayesianFallback_FromFile validates reviewer
// correctness finding #2: default invocation (Mode=Bayesian) must
// emit a stderr warning and fall back to simple.
func TestConsensus_BayesianFallback_FromFile(t *testing.T) {
	// We use ConsensusFile so that the errOut writer gets the warning
	// — Consensus() alone doesn't take a stderr handle.
	// Need a real file path on disk because ConsensusFile opens it.
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
		t.Fatalf("ConsensusFile(bayesian fallback): %v", err)
	}
	if !strings.Contains(serr.String(), "bayesian") {
		t.Errorf("expected stderr warning mentioning bayesian, got %q", serr.String())
	}
	if !strings.Contains(serr.String(), "falling back to simple") {
		t.Errorf("expected stderr warning mentioning 'falling back to simple', got %q", serr.String())
	}
	if !strings.Contains(sout.String(), "ACGTA") {
		t.Errorf("bayesian fallback should still emit consensus, got %q", sout.String())
	}
}

// TestConsensus_BayesianFallback_LibraryCall confirms that calling
// Consensus() directly with Mode=Bayesian also falls back (no warning,
// because there's no stderr handle, but no crash either — the
// behaviour is consistent with simple mode).
func TestConsensus_BayesianFallback_LibraryCall(t *testing.T) {
	var sout bytes.Buffer
	opts := ConsensusOptions{
		Format: ConsensusFASTA,
		Mode:   ConsensusModeBayesian,
	}
	if err := Consensus(strings.NewReader(allMatchSAM), &sout, opts); err != nil {
		t.Fatalf("Consensus(bayesian fallback): %v", err)
	}
	if !strings.Contains(sout.String(), "ACGTA") {
		t.Errorf("bayesian fallback should still emit consensus, got %q", sout.String())
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

// TestConsensus_UpstreamParityStub is a deliberate skip for the
// upstream `samtools test/consensus/` corpus. Re-enable once bayesian
// mode lands. Tracked in docs/PARITY_ROADMAP.md#samtools.
func TestConsensus_UpstreamParityStub(t *testing.T) {
	t.Skip("upstream samtools/test/consensus/ corpus exercises bayesian-mode by default; simple-mode parity covered above. Re-enable once bayesian mode lands.")
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
