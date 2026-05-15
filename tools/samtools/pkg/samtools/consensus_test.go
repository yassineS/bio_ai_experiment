package samtools

import (
	"bytes"
	"strings"
	"testing"
)

// runConsensusOnSAM feeds one SAM-text input through Consensus() with the
// supplied options and returns the emitted text.
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
	// With qual 'I' = 40 and three reads agreeing, the score for A is
	// 3*40=120 and tscore is the same (only one base). So qual byte
	// is 100 -> capped at 93 -> phred byte = 33+93 = 126 = '~'.
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
	// First row: chr1\t1\t0\t3\tA\t100\tAAA\tIII
	if lines[0] != "chr1\t1\t0\t3\tA\t100\tAAA\tIII" {
		t.Errorf("pileup row 0 = %q", lines[0])
	}
	if lines[4] != "chr1\t5\t0\t3\tA\t100\tAAA\tIII" {
		t.Errorf("pileup row 4 = %q", lines[4])
	}
}

// mixedSAM exercises positions where one read disagrees out of four:
//
//	pos 1: A,A,A,A  -> A (4/4)
//	pos 2: C,C,C,T  -> C (3/4 = 0.75 = MinCallFraction default)
//	pos 3: G,G,G,G  -> G (4/4)
//	pos 4: T,T,T,A  -> T (3/4)
//	pos 5: A,A,A,A  -> A (4/4)
const mixedSAM = `@HD	VN:1.6
@SQ	SN:chr1	LN:8
r1	0	chr1	1	60	5M	*	0	0	ACGTA	IIIII
r2	0	chr1	1	60	5M	*	0	0	ACGTA	IIIII
r3	0	chr1	1	60	5M	*	0	0	ACGTA	IIIII
r4	0	chr1	1	60	5M	*	0	0	ATGAA	IIIII
`

func TestConsensus_FASTA_Mixed(t *testing.T) {
	// Defaults: MinCallFraction=0.75, MinConsensusFraction=0.6.
	// 3/4 = 0.75 -> both fractions pass at positions 2 and 4.
	out := runConsensusOnSAM(t, mixedSAM, ConsensusOptions{
		Format: ConsensusFASTA,
	})
	want := ">chr1\nACGTA\n"
	if out != want {
		t.Errorf("mixed FASTA: got %q want %q", out, want)
	}
}

func TestConsensus_FASTA_BelowMinFraction_BecomesN(t *testing.T) {
	// Bump MinConsensusFraction past 0.75 so positions 2 and 4 (with
	// 3/4 = 0.75 majority) become N.
	out := runConsensusOnSAM(t, mixedSAM, ConsensusOptions{
		Format:               ConsensusFASTA,
		MinConsensusFraction: 0.9,
	})
	want := ">chr1\nANGNA\n"
	if out != want {
		t.Errorf("low-fraction FASTA: got %q want %q", out, want)
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

// insSAM has three reads with a 1bp insertion T between ref positions 2 and 3.
const insSAM = `@HD	VN:1.6
@SQ	SN:chr1	LN:4
r1	0	chr1	1	60	2M1I2M	*	0	0	ACTGA	IIIII
r2	0	chr1	1	60	2M1I2M	*	0	0	ACTGA	IIIII
r3	0	chr1	1	60	2M1I2M	*	0	0	ACTGA	IIIII
`

func TestConsensus_FASTA_Insertion_Included(t *testing.T) {
	out := runConsensusOnSAM(t, insSAM, ConsensusOptions{Format: ConsensusFASTA})
	// Reference has 4 bases (ACGA); insertion T after pos 2 brings us
	// to ACTGA when ShowIns=true (the default).
	want := ">chr1\nACTGA\n"
	if out != want {
		t.Errorf("insertion-include: got %q want %q", out, want)
	}
}

func TestConsensus_FASTA_Insertion_Suppressed(t *testing.T) {
	out := runConsensusOnSAM(t, insSAM, ConsensusOptions{
		Format:    ConsensusFASTA,
		NoShowIns: true, // override default include-insertions
	})
	want := ">chr1\nACGA\n"
	if out != want {
		t.Errorf("insertion-suppress: got %q want %q", out, want)
	}
}

func TestConsensus_LineLen_Wrapping(t *testing.T) {
	// Force a 2-column wrap.
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
	// allMatchSAM covers 1..5 on LN=8; -a should emit rows 6..8 as zero
	// depth (no reads), which call N with cq=0 and "*" for seq/qual.
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

func TestConsensus_AmbigCodes_Het(t *testing.T) {
	// Two reads A, two reads C at pos 1; with -A and high MinHetFraction
	// the call should be M (A+C ambiguity code).
	sam := `@HD	VN:1.6
@SQ	SN:chr1	LN:1
r1	0	chr1	1	60	1M	*	0	0	A	I
r2	0	chr1	1	60	1M	*	0	0	A	I
r3	0	chr1	1	60	1M	*	0	0	C	I
r4	0	chr1	1	60	1M	*	0	0	C	I
`
	out := runConsensusOnSAM(t, sam, ConsensusOptions{
		Format:               ConsensusFASTA,
		AmbigCodes:           true,
		MinHetFraction:       0.5,
		MinConsensusFraction: 0.4, // dominant alone is 0.5
		MinCallFraction:      0.9,
	})
	if !strings.Contains(out, "M\n") {
		t.Errorf("expected M ambig code, got %q", out)
	}
}

func TestConsensus_BayesianFallback(t *testing.T) {
	// Mode=bayesian should warn on stderr and fall back to simple.
	var sout, serr bytes.Buffer
	opts := ConsensusOptions{
		Format: ConsensusFASTA,
		Mode:   ConsensusModeBayesian,
	}
	// We use the streaming entrypoint directly to capture stderr via
	// ConsensusFile-like routing — but Consensus doesn't take a stderr.
	// Instead we verify the fallback path doesn't panic and produces
	// the same result as simple.
	if err := Consensus(strings.NewReader(allMatchSAM), &sout, opts); err != nil {
		t.Fatalf("Consensus(bayesian fallback): %v", err)
	}
	if !strings.Contains(sout.String(), "ACGTA") {
		t.Errorf("bayesian fallback should still emit consensus, got %q", sout.String())
	}
	_ = serr
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

// TestConsensus_UpstreamParityStub is a deliberate skip placeholder for
// the upstream `samtools test/consensus/` corpus. The v1 implementation
// only ships simple-mode FASTA/FASTQ/pileup; running the bayesian-mode
// corpus is meaningless until the Gap5 posterior caller lands.
//
// Tracked in docs/PARITY_ROADMAP.md#samtools (consensus bayesian-mode).
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
