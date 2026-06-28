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

// delRunMinQualSAM has a single read with a 2bp deletion (3M2D3M) where the
// PRE-gap base (query index 2, ref pos 3) has a LOW quality (Phred 10 = '+')
// and the POST-gap base (query index 3, ref pos 6) has a HIGH quality
// (Phred 40 = 'I'). The deletion '*' placeholders fall at ref positions 4 and 5.
//
//	ref pos:  1  2  3 | 4  5 | 6  7  8
//	op:       M  M  M | D  D | M  M  M
//	qry idx:  0  1  2 |      | 3  4  5
//	qual:     I  I  + |      | I  I  I   (Phred 40,40,10,40,40,40)
const delRunMinQualSAM = `@HD	VN:1.6
@SQ	SN:chr1	LN:20
r1	0	chr1	1	60	3M2D3M	*	0	0	ACGTAC	II+III
`

// TestConsensus_Pileup_DeletionRunningMinQual locks upstream's RUNNING-minimum
// deletion-quality rule for `samtools consensus -f pileup --mode simple`
// (consensus_pileup.c:195-202): each '*' placeholder gets MIN(pre-gap base
// qual, post-gap base qual). Here MIN(10, 40) = 10, so the '*' rows render
// quality byte '+' (10+33=43), NOT the post-gap 'I' (40+33=73) that the raw
// per-read placeholder quality would give. Crucially, the SAME read fed
// through mpileup must keep the post-gap quality 'I' for its '*' placeholders
// (mpileup's bam_plp engine renders the single post-gap base, no running min) —
// the consensus running-min must NOT leak into the mpileup path.
func TestConsensus_Pileup_DeletionRunningMinQual(t *testing.T) {
	// (1) consensus -f pileup --mode simple: '*' qual == MIN(pre,post) = '+'.
	out := runConsensusOnSAM(t, delRunMinQualSAM, ConsensusOptions{
		Format:  ConsensusPileup,
		Mode:    ConsensusModeSimple,
		ShowDel: true,
	})
	lines := strings.Split(strings.TrimRight(out, "\n"), "\n")
	starRows := 0
	for _, ln := range lines {
		f := strings.Split(ln, "\t")
		if len(f) < 8 || f[6] != "*" {
			continue
		}
		starRows++
		// Column 8 (index 7) is the quality string; one '*' read here.
		if f[7] != "+" {
			t.Errorf("pos %s: deletion '*' qual = %q, want %q (running MIN(pre=10,post=40)+33); full row %q",
				f[1], f[7], "+", ln)
		}
	}
	if starRows != 2 {
		t.Fatalf("expected 2 '*' deletion rows (ref pos 4,5), got %d:\n%s", starRows, out)
	}

	// (2) mpileup on the identical read: '*' placeholders keep the POST-gap
	// quality 'I' (40+33). The consensus running-min override is render-scoped
	// to the consensus pileup writer and must not perturb mpileup.
	mp := runMpileupOnSAM(t, []string{delRunMinQualSAM}, MpileupOptions{}, nil, nil)
	mpLines := strings.Split(strings.TrimRight(mp, "\n"), "\n")
	mpStar := 0
	for _, ln := range mpLines {
		f := strings.Split(ln, "\t")
		// mpileup columns: chrom pos ref depth bases quals ...
		if len(f) < 6 || !strings.Contains(f[4], "*") {
			continue
		}
		mpStar++
		if f[5] != "I" {
			t.Errorf("mpileup pos %s: '*' qual = %q, want %q (post-gap 40, unchanged); full row %q",
				f[1], f[5], "I", ln)
		}
	}
	if mpStar != 2 {
		t.Fatalf("expected 2 mpileup '*' rows (ref pos 4,5), got %d:\n%s", mpStar, mp)
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

// TestConsensus_Pileup_Insertion_PadRunningMin pins task #49: in
// `consensus -f pileup --mode simple` the '*' pad emitted for a read that
// LACKS an inserted base in an nth>0 insertion column must carry upstream's
// STATEFUL running-minimum quality, MIN(carried, b_qual[seq_offset+1]) with
// seq_offset pinned at the read's last consumed base (consensus_pileup.c:
// 183-191), not the read's M-base quality verbatim.
//
// Here rIns inserts "TT" after ref pos 2 (CIGAR 2M2I3M); rPad has a plain 5M
// alignment, so in the nth=1/nth=2 insertion columns rPad pads with '*'. The
// base at ref pos 2 in rPad has quality 'I' (40), but the NEXT base
// (seq_offset+1, ref pos 3) has quality '5' (20). Upstream therefore renders
// the pad '*' qual as the running minimum 20 ('5'), while the real inserted
// base of rIns keeps its own quality 'I' (40). Without the running min the pad
// would wrongly print 'I' (the pre-task-#49 behaviour).
func TestConsensus_Pileup_Insertion_PadRunningMin(t *testing.T) {
	const sam = `@HD	VN:1.6
@SQ	SN:chr1	LN:10
rIns	0	chr1	1	60	2M2I3M	*	0	0	ACTTGAA	IIIIIII
rPad	0	chr1	1	60	5M	*	0	0	ACGAA	II5II
`
	out := runConsensusOnSAM(t, sam, ConsensusOptions{
		Format: ConsensusPileup,
		Mode:   ConsensusModeSimple,
	})
	lines := strings.Split(strings.TrimRight(out, "\n"), "\n")
	// Collect the two nth>0 insertion rows (chrom pos 2, nth 1 and 2).
	var ins []string
	for _, ln := range lines {
		f := strings.Split(ln, "\t")
		if len(f) >= 3 && f[1] == "2" && (f[2] == "1" || f[2] == "2") {
			ins = append(ins, ln)
		}
	}
	if len(ins) != 2 {
		t.Fatalf("expected 2 insertion rows (pos 2, nth 1/2), got %d:\n%s", len(ins), out)
	}
	// Stable event order is rIns then rPad, so the seq is "T*" and the qual
	// column must be the inserted base's 'I' followed by the pad's running
	// min '5' (NOT 'I'). Assert both nth rows.
	for _, ln := range ins {
		f := strings.Split(ln, "\t")
		seq, qual := f[6], f[7]
		if seq != "T*" {
			t.Errorf("insertion row %q: seq = %q, want %q", ln, seq, "T*")
		}
		if qual != "I5" {
			t.Errorf("insertion row %q: qual = %q, want %q (inserted base 'I' unchanged, pad '*' = running MIN '5')", ln, qual, "I5")
		}
		// Guard against a silent regression to the M-base quality: the pad's
		// qual byte must be the running minimum, strictly below the M-base.
		if qual[1] != '5' {
			t.Errorf("insertion row %q: pad '*' qual = %q, want running min '5'; got the M-base qual instead", ln, string(qual[1]))
		}
	}
}

// TestConsensus_Pileup_DeletionPadRunningMin pins task #57: a read whose '*'
// deletion placeholder is padded into ANOTHER read's nth>0 insertion column must
// carry upstream's running-minimum deletion quality (the MIN(pre-gap, post-gap)
// computed at the nth==0 D column), with seq_offset PINNED at the pre-gap base —
// NOT seeded from the post-gap base index.
//
// Geometry (mirrors the GIAB 20:18390227 / 20:22642911 residual in miniature):
//
//	rIns 3M2I3M inserts "TT" after ref pos 3, opening nth=1/nth=2 columns there.
//	rDel 2M1D5M has its deletion '*' AT ref pos 3 (the same base column the
//	     insertion follows), so in the nth>0 columns rDel pads with '*'.
//
// rDel's quals are I I I 5 I I I (Phred 40,40,40,20,40,40,40). The deletion is at
// query offset 2: the PRE-gap base (query 1) and POST-gap base (query 2) are both
// 40, so the running minimum is 40 ('I'). The base one past the post-gap (query 3)
// is 20 ('5'). Upstream pins seq_offset at the pre-gap base and MINs only against
// the post-gap base, so every '*' (the nth==0 D row AND the nth=1/2 pad rows)
// renders 'I'. The pre-fix code seeded the carried seq_offset from the post-gap
// base index, so the FIRST pad column wrongly MIN'd against query 3 (20), printing
// '5' in the insertion rows. The nth==0 D row was already correct (task #48), so
// this asserts specifically the insertion-pad propagation.
func TestConsensus_Pileup_DeletionPadRunningMin(t *testing.T) {
	const sam = `@HD	VN:1.6
@SQ	SN:chr1	LN:20
rIns	0	chr1	1	60	3M2I3M	*	0	0	ACGTTACG	IIIIIIII
rDel	0	chr1	1	60	2M1D5M	*	0	0	ACTACGT	III5III
`
	out := runConsensusOnSAM(t, sam, ConsensusOptions{
		Format:  ConsensusPileup,
		Mode:    ConsensusModeSimple,
		ShowDel: true,
	})
	lines := strings.Split(strings.TrimRight(out, "\n"), "\n")
	// Stable event order is rIns then rDel, so at ref pos 3 the seq is the
	// inserted/aligned base of rIns followed by rDel's '*'. The qual column's
	// second byte is rDel's deletion placeholder: it must be 'I' (the running
	// min 40), never '5' (the pre-fix base-one-past-post-gap value 20).
	var nth0, padRows int
	for _, ln := range lines {
		f := strings.Split(ln, "\t")
		if len(f) < 8 || f[1] != "3" {
			continue
		}
		seq, qual := f[6], f[7]
		if len(seq) < 2 || seq[1] != '*' {
			t.Fatalf("pos 3 nth %s: expected rDel '*' as the 2nd seq byte, got seq %q (row %q)", f[2], seq, ln)
		}
		if len(qual) < 2 {
			t.Fatalf("pos 3 nth %s: short qual %q (row %q)", f[2], qual, ln)
		}
		if qual[1] != 'I' {
			t.Errorf("pos 3 nth %s: deletion-pad '*' qual = %q, want 'I' (running MIN(pre=40,post=40)); the '5' value means the pad read one base past the post-gap (row %q)",
				f[2], string(qual[1]), ln)
		}
		if f[2] == "0" {
			nth0++
		} else {
			padRows++
		}
	}
	if nth0 != 1 {
		t.Fatalf("expected 1 nth==0 D row at pos 3, got %d:\n%s", nth0, out)
	}
	if padRows != 2 {
		t.Fatalf("expected 2 insertion-pad rows (pos 3, nth 1/2), got %d:\n%s", padRows, out)
	}

	// Render-scope guard: the SAME read through mpileup must keep its post-gap
	// placeholder quality (the bam_plp engine renders the single post-gap base,
	// no consensus running min), proving the task-#57 fix did not leak into the
	// shared pileupEvent.qual the mpileup path consumes.
	mp := runMpileupOnSAM(t, []string{sam}, MpileupOptions{}, nil, nil)
	for _, ln := range strings.Split(strings.TrimRight(mp, "\n"), "\n") {
		f := strings.Split(ln, "\t")
		if len(f) < 6 || f[1] != "3" || !strings.Contains(f[4], "*") {
			continue
		}
		// rDel's '*' at pos 3 borrows the post-gap base quality 'I' (40);
		// the nth>0 insertion columns never appear in mpileup output.
		if !strings.Contains(f[5], "I") {
			t.Errorf("mpileup pos 3: '*' qual column = %q, want to contain post-gap 'I'; row %q", f[5], ln)
		}
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

// delHomopolymerRunMinSAM exercises the task #55 regression: the bayesian
// consensus caller must feed a deletion '*' placeholder the RUNNING-MINIMUM
// quality MIN(pre-gap base qual, post-gap base qual) (consensus_pileup.c:
// 195-202), NOT the post-gap base quality alone. The reference run is a 6-C
// homopolymer at chr1:2-7 (ref ACCCCCCGT... at 1-12). Some reads carry a 1bp
// deletion at the run's trailing edge (CIGAR 6M1D2M dropping the 7th ref base,
// the last C at pos 7); for those reads the PRE-gap base (the read's 6th base,
// quality '#'=2) is far lower than the POST-gap base (the G at pos 8, quality
// 'I'=40). Upstream MIN(2,40)=2, so each '*' contributes only quality 2 to the
// deletion-vs-base posterior at pos 7. With the (buggy) post-gap quality 40 the
// deletions look highly confident and outvote the C-supporting reads, dropping
// the run-edge base and frameshifting the consensus to a 5-C run. With the
// running minimum the deletions are weak and the 6th C is correctly called.
const delHomopolymerRunMinSAM = `@HD	VN:1.6
@SQ	SN:chr1	LN:12
c1	0	chr1	1	60	9M	*	0	0	ACCCCCCGT	IIIIIIIII
c2	0	chr1	1	60	9M	*	0	0	ACCCCCCGT	IIIIIIIII
d1	0	chr1	1	60	6M1D2M	*	0	0	ACCCCCGT	IIIII#II
d2	0	chr1	1	60	6M1D2M	*	0	0	ACCCCCGT	IIIII#II
d3	0	chr1	1	60	6M1D2M	*	0	0	ACCCCCGT	IIIII#II
d4	0	chr1	1	60	6M1D2M	*	0	0	ACCCCCGT	IIIII#II
d5	0	chr1	1	60	6M1D2M	*	0	0	ACCCCCGT	IIIII#II
`

// TestConsensus_FASTA_Bayesian_DeletionRunMinQual is the task #55 regression.
// It asserts the bayesian -f fasta caller CALLS the homopolymer run-edge base
// (a full 6-C run, "ACCCCCCGT...") rather than dropping it to a 5-C run. The
// deletion reads carry a HIGH post-gap quality but a LOW pre-gap quality; the
// call is correct only when the '*' placeholder uses the running minimum. This
// test FAILS when the del path does `bp.qual = e.qual` (post-gap only) and
// PASSES once it reads the running-min e.delPileupQual.
func TestConsensus_FASTA_Bayesian_DeletionRunMinQual(t *testing.T) {
	out := runConsensusOnSAM(t, delHomopolymerRunMinSAM, ConsensusOptions{
		Format: ConsensusFASTA,
		Mode:   ConsensusModeBayesian,
	})
	// Extract the unwrapped consensus sequence (skip the '>' header line).
	var seq string
	for _, ln := range strings.Split(strings.TrimRight(out, "\n"), "\n") {
		if strings.HasPrefix(ln, ">") {
			continue
		}
		seq += ln
	}
	cCount := strings.Count(seq, "C")
	if cCount != 6 {
		t.Errorf("bayesian -f fasta dropped a homopolymer run-edge base: got %q "+
			"(C-count %d), want a full 6-C run (the running-min deletion quality "+
			"MIN(pre=2,post=40) must keep the trailing C called); full output %q",
			seq, cCount, out)
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
// The reference-FASTA (-T) cases (the *T.out goldens) ARE now exercised: the
// no-coverage / gap ref-base substitution is implemented and byte-validated.
func TestConsensus_BayesianUpstreamParity(t *testing.T) {
	bin := upstreamSamtools(t)
	const fixDir = "../../../../reference_code/samtools/test/consensus"
	// refFixture is the per-case -T reference FASTA (empty when no -T).
	cases := []struct {
		name  string
		input string
		ref   string   // reference FASTA basename for -T, or "" for none
		args  []string // upstream `samtools consensus` args (sans input file)
		opts  ConsensusOptions
	}{
		{"18q", "consen1.sam", "", []string{"-f", "fastq", "--no-use-MQ", "-C", "0", "-m", "bayesian"}, consBayes(ConsensusFASTQ, 0, false, false, "", "")},
		{"19q", "consen1.sam", "", []string{"-f", "fastq", "--no-use-MQ", "-C", "19", "-m", "bayesian"}, consBayes(ConsensusFASTQ, 19, false, false, "", "")},
		{"18p", "consen1.sam", "", []string{"-f", "pileup", "--no-use-MQ", "-C", "0", "-m", "bayesian"}, consBayes(ConsensusPileup, 0, false, false, "", "")},
		{"19p", "consen1.sam", "", []string{"-f", "pileup", "--no-use-MQ", "-C", "19", "-m", "bayesian"}, consBayes(ConsensusPileup, 19, false, false, "", "")},
		{"20p", "consen1.sam", "", []string{"-f", "pileup", "--no-use-MQ", "-C", "30", "-A", "-m", "bayesian"}, consBayes(ConsensusPileup, 30, false, true, "", "")},
		{"21p", "consen1.sam", "", []string{"-f", "pileup", "--no-use-MQ", "-C", "31", "-A", "-m", "bayesian"}, consBayes(ConsensusPileup, 31, false, true, "", "")},
		{"30", "consen1c.sam", "", []string{"-f", "fastq", "--show-del", "yes", "--show-ins", "no", "-m", "bayesian", "-C0"}, consBayes(ConsensusFASTQ, 0, true, false, "yes", "no")},
		{"31", "consen1c.sam", "", []string{"-a", "-f", "fastq", "--show-del", "yes", "--show-ins", "no", "-m", "bayesian", "-C0"}, consBayesA(ConsensusFASTQ, 0, true, false, "yes", "no", 1)},
		{"32", "consen1c.sam", "", []string{"-aa", "-f", "fastq", "--show-del", "yes", "--show-ins", "no", "-m", "bayesian", "-C0"}, consBayesA(ConsensusFASTQ, 0, true, false, "yes", "no", 2)},
		{"40", "consen1c.sam", "", []string{"-f", "pileup", "--show-del", "yes", "--show-ins", "no", "-m", "bayesian", "-C0"}, consBayes(ConsensusPileup, 0, true, false, "yes", "no")},
		{"41", "consen1c.sam", "", []string{"-a", "-f", "pileup", "--show-del", "yes", "--show-ins", "no", "-m", "bayesian", "-C0"}, consBayesA(ConsensusPileup, 0, true, false, "yes", "no", 1)},
		{"42", "consen1c.sam", "", []string{"-aa", "-f", "pileup", "--show-del", "yes", "--show-ins", "no", "-m", "bayesian", "-C0"}, consBayesA(ConsensusPileup, 0, true, false, "yes", "no", 2)},
		// -T/--reference cases (consensus.reg 30T..42T): the no-coverage / gap
		// positions substitute consen1c.fa reference bases at --ref-qual 20.
		{"30T", "consen1c.sam", "consen1c.fa", []string{"-f", "fastq", "--show-del", "yes", "--show-ins", "no", "-m", "bayesian", "-C0", "--ref-qual", "20"}, consBayesRef(ConsensusFASTQ, 0, true, false, "yes", "no", 0, 20)},
		{"31T", "consen1c.sam", "consen1c.fa", []string{"-a", "-f", "fastq", "--show-del", "yes", "--show-ins", "no", "-m", "bayesian", "-C0", "--ref-qual", "20"}, consBayesRef(ConsensusFASTQ, 0, true, false, "yes", "no", 1, 20)},
		{"32T", "consen1c.sam", "consen1c.fa", []string{"-aa", "-f", "fastq", "--show-del", "yes", "--show-ins", "no", "-m", "bayesian", "-C0", "--ref-qual", "20"}, consBayesRef(ConsensusFASTQ, 0, true, false, "yes", "no", 2, 20)},
		{"40T", "consen1c.sam", "consen1c.fa", []string{"-f", "pileup", "--show-del", "yes", "--show-ins", "no", "-m", "bayesian", "-C0", "--ref-qual", "20"}, consBayesRef(ConsensusPileup, 0, true, false, "yes", "no", 0, 20)},
		{"41T", "consen1c.sam", "consen1c.fa", []string{"-a", "-f", "pileup", "--show-del", "yes", "--show-ins", "no", "-m", "bayesian", "-C0", "--ref-qual", "20"}, consBayesRef(ConsensusPileup, 0, true, false, "yes", "no", 1, 20)},
		{"42T", "consen1c.sam", "consen1c.fa", []string{"-aa", "-f", "pileup", "--show-del", "yes", "--show-ins", "no", "-m", "bayesian", "-C0", "--ref-qual", "20"}, consBayesRef(ConsensusPileup, 0, true, false, "yes", "no", 2, 20)},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			input := filepath.Join(fixDir, c.input)
			samBytes, err := os.ReadFile(input)
			if err != nil {
				t.Fatalf("read input %s: %v", c.input, err)
			}

			// Live upstream invocation. -T cases append "-T <ref>".
			upArgs := append([]string{"consensus"}, c.args...)
			if c.ref != "" {
				upArgs = append(upArgs, "-T", filepath.Join(fixDir, c.ref))
			}
			upArgs = append(upArgs, input)
			cmd := exec.Command(bin, upArgs...)
			var upOut, upErr bytes.Buffer
			cmd.Stdout = &upOut
			cmd.Stderr = &upErr
			if err := cmd.Run(); err != nil {
				t.Fatalf("upstream samtools consensus %v: %v\n%s", c.args, err, upErr.String())
			}

			// Go port.
			opts := c.opts
			if c.ref != "" {
				opts.Reference = filepath.Join(fixDir, c.ref)
			}
			got := runConsensusOnSAM(t, string(samBytes), opts)
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

// consBayesRef is consBayesA plus a --ref-qual setting for the -T cases. The
// Reference path itself is filled in by the test runner from the fixture dir.
func consBayesRef(f ConsensusFormat, cutoff int, useMQ, ambig bool, showDel, showIns string, allLevel, refQual int) ConsensusOptions {
	o := consBayesA(f, cutoff, useMQ, ambig, showDel, showIns, allLevel)
	o.RefQual = refQual
	return o
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

// delGapSAM places a 3bp deletion (5M3D5M at pos 3 → deletions at 8,9,10) on
// three reads of a 16bp contig, covering positions 3-7 and 11-15, with a
// leading (1-2) and a trailing (16) zero-coverage gap. It is the deterministic
// counterpart of the live
// all-positions parity sweep, pinning the placeholder-row format and the
// upstream duplicate-row quirk at deletion sites without needing the upstream
// binary.
const delGapSAM = `@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:16
d1	0	chr1	3	60	5M3D5M	*	0	0	ACGTAACGTA	IIIIIIIIII
d2	0	chr1	3	60	5M3D5M	*	0	0	ACGTAACGTA	IIIIIIIIII
d3	0	chr1	3	60	5M3D5M	*	0	0	ACGTAACGTA	IIIIIIIIII
`

// TestConsensus_AllPositionsPileup_DeletionAndGaps pins the exact -a pileup
// output: leading/trailing zero-coverage placeholder rows, and the deletion
// run 8,9,10 emitted as placeholders with upstream's duplicate-row quirk
// (last_pos is not advanced for a suppressed '*' column, so each deletion
// position is re-emitted by every following column: 8×3, 9×2, 10×1).
func TestConsensus_AllPositionsPileup_DeletionAndGaps(t *testing.T) {
	got := runConsensusOnSAM(t, delGapSAM, ConsensusOptions{
		Mode: ConsensusModeSimple, Format: ConsensusPileup, AllPositions: true})
	want := "chr1\t1\t0\t0\tN\t0\t*\t*\n" +
		"chr1\t2\t0\t0\tN\t0\t*\t*\n" +
		"chr1\t3\t0\t3\tA\t100\tAAA\tIII\n" +
		"chr1\t4\t0\t3\tC\t100\tCCC\tIII\n" +
		"chr1\t5\t0\t3\tG\t100\tGGG\tIII\n" +
		"chr1\t6\t0\t3\tT\t100\tTTT\tIII\n" +
		"chr1\t7\t0\t3\tA\t100\tAAA\tIII\n" +
		// Deletion run 8,9,10 with the duplicate quirk.
		"chr1\t8\t0\t0\tN\t0\t*\t*\n" +
		"chr1\t8\t0\t0\tN\t0\t*\t*\n" +
		"chr1\t9\t0\t0\tN\t0\t*\t*\n" +
		"chr1\t8\t0\t0\tN\t0\t*\t*\n" +
		"chr1\t9\t0\t0\tN\t0\t*\t*\n" +
		"chr1\t10\t0\t0\tN\t0\t*\t*\n" +
		"chr1\t11\t0\t3\tA\t100\tAAA\tIII\n" +
		"chr1\t12\t0\t3\tC\t100\tCCC\tIII\n" +
		"chr1\t13\t0\t3\tG\t100\tGGG\tIII\n" +
		"chr1\t14\t0\t3\tT\t100\tTTT\tIII\n" +
		"chr1\t15\t0\t3\tA\t100\tAAA\tIII\n" +
		// Trailing zero-coverage tail fill.
		"chr1\t16\t0\t0\tN\t0\t*\t*\n"
	if got != want {
		t.Fatalf("all-positions pileup mismatch:\n--- got ---\n%s\n--- want ---\n%s", got, want)
	}
}

// TestConsensus_AllPositionsPileup_ShowDelKeepsDeletionRows checks that with
// --show-del the deletion columns are emitted as genuine '*' rows (advancing
// last_pos), so the duplicate placeholder quirk disappears.
func TestConsensus_AllPositionsPileup_ShowDelKeepsDeletionRows(t *testing.T) {
	got := runConsensusOnSAM(t, delGapSAM, ConsensusOptions{
		Format: ConsensusPileup, AllPositions: true, ShowDel: true})
	// With show-del the deletion columns are real rows: depth 3, call '*',
	// each emitted exactly once and no duplicate placeholders.
	for _, pos := range []string{"chr1\t8\t0\t3\t*\t", "chr1\t9\t0\t3\t*\t", "chr1\t10\t0\t3\t*\t"} {
		if strings.Count(got, pos) != 1 {
			t.Fatalf("expected exactly one show-del deletion row %q, got %d:\n%s",
				pos, strings.Count(got, pos), got)
		}
	}
	if strings.Contains(got, "chr1\t8\t0\t0\tN") {
		t.Fatalf("show-del should not emit placeholder rows at deletion sites:\n%s", got)
	}
}

// TestConsensus_AllPositionsPileup_BEDFilterScopesPlaceholders verifies the
// -l/BED position filter (opts.BEDPath) confines the -a placeholder fill to
// the selected positions: rows for excluded coordinates are never emitted,
// not even as depth-0 placeholders.
func TestConsensus_AllPositionsPileup_BEDFilterScopesPlaceholders(t *testing.T) {
	dir := t.TempDir()
	bed := filepath.Join(dir, "sel.bed")
	// BED is 0-based half-open: chr1 1..4 selects 1-based positions 2,3,4.
	if err := os.WriteFile(bed, []byte("chr1\t1\t4\n"), 0o600); err != nil {
		t.Fatalf("write bed: %v", err)
	}
	got := runConsensusOnSAM(t, delGapSAM, ConsensusOptions{
		Mode: ConsensusModeSimple, Format: ConsensusPileup, AllPositions: true, BEDPath: bed})
	want := "chr1\t2\t0\t0\tN\t0\t*\t*\n" +
		"chr1\t3\t0\t3\tA\t100\tAAA\tIII\n" +
		"chr1\t4\t0\t3\tC\t100\tCCC\tIII\n"
	if got != want {
		t.Fatalf("BED-filtered all-positions pileup mismatch:\n--- got ---\n%s\n--- want ---\n%s", got, want)
	}
}

// nIslandGapSAM reproduces the whole-contig N-island gap that task #56 closes.
// Two reads cover chr1:1-5 ("ACGTA"); then positions 6-14 have ZERO coverage
// (an internal gap that the FASTA/FASTQ accumulator N-fills inline); then two
// reads at chr1:15-19 open with an 'N' base ("NCGTA"), so the default
// (bayesian) caller renders column 15 a deletion '*'. Upstream's basic_fasta
// returns early for that suppressed '*' BEFORE its gap-fill loop runs
// (bam_consensus.c:2403-2407), so the immediately-preceding zero-coverage gap
// is NEVER written — the byte-exact upstream output (verified against the
// samtools binary) is:
//
//	show-del off (default): >chr1\nACGTACGTA      (gap + '*' swallowed)
//	show-del yes:           >chr1\nACGTANNNNNNNNN*CGTA  ('*' emitted, gap kept)
const nIslandGapSAM = `@HD	VN:1.6
@SQ	SN:chr1	LN:30
r1	0	chr1	1	60	5M	*	0	0	ACGTA	IIIII
r2	0	chr1	1	60	5M	*	0	0	ACGTA	IIIII
r3	0	chr1	15	60	5M	*	0	0	NCGTA	IIIII
r4	0	chr1	15	60	5M	*	0	0	NCGTA	IIIII
`

// TestConsensus_FASTA_NIslandGap_SwallowedWhenShowDelOff is the task #56
// regression: with --show-del off (the default) the zero-coverage gap that
// immediately precedes a suppressed '*' deletion column must be SWALLOWED, not
// left as an N-run. It FAILS before the gapRunStart rollback (the pre-fix
// output is "ACGTANNNNNNNNNCGTA", the nine gap N's leaking through) and PASSES
// after.
func TestConsensus_FASTA_NIslandGap_SwallowedWhenShowDelOff(t *testing.T) {
	// Bayesian mode is the upstream `samtools consensus` default and is what
	// renders the post-gap 'N' base as a deletion '*' (simple mode would call
	// it 'N', not exercising the suppressed-'*' rollback path under test).
	out := runConsensusOnSAM(t, nIslandGapSAM, ConsensusOptions{
		Mode: ConsensusModeBayesian, Format: ConsensusFASTA})
	want := ">chr1\nACGTACGTA\n"
	if out != want {
		t.Errorf("N-island gap not swallowed: got %q want %q", out, want)
	}
	if strings.Contains(out, "N") {
		t.Errorf("gap N-run leaked into show-del-off output: %q", out)
	}
}

// TestConsensus_FASTQ_NIslandGap_QualLockstep verifies the rollback truncates
// qualBuf in lockstep with seqBuf for -f fastq: the swallowed gap must leave no
// orphaned quality bytes (seq and qual must be equal length and match upstream).
func TestConsensus_FASTQ_NIslandGap_QualLockstep(t *testing.T) {
	out := runConsensusOnSAM(t, nIslandGapSAM, ConsensusOptions{
		Mode: ConsensusModeBayesian, Format: ConsensusFASTQ})
	want := "@chr1\nACGTACGTA\n+\nGGGGGGGGG\n"
	if out != want {
		t.Errorf("fastq gap lockstep mismatch: got %q want %q", out, want)
	}
}

// TestConsensus_FASTA_NIslandGap_EmittedWhenShowDelYes is the negative control:
// with --show-del yes the '*' IS emitted and the rollback must NOT fire, so the
// gap N-run is kept exactly as upstream emits it.
func TestConsensus_FASTA_NIslandGap_EmittedWhenShowDelYes(t *testing.T) {
	out := runConsensusOnSAM(t, nIslandGapSAM, ConsensusOptions{
		Mode: ConsensusModeBayesian, Format: ConsensusFASTA, ShowDel: true})
	want := ">chr1\nACGTANNNNNNNNN*CGTA\n"
	if out != want {
		t.Errorf("show-del-yes must keep gap + '*': got %q want %q", out, want)
	}
}
