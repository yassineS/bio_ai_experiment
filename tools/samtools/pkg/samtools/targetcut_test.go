package samtools

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/errmod"
)

// ----- Simple-mode tests (legacy behaviour, behind SimpleMode opt-in) -----

// TestTargetcutSimple_Basic exercises soft-clipped flanks, insertions
// (kept), deletions (no query bases consumed), and unmapped /
// secondary record skipping in the legacy aligned-slice FASTA mode.
func TestTargetcutSimple_Basic(t *testing.T) {
	samText := strings.Join([]string{
		"@HD\tVN:1.6\tSO:coordinate",
		"@SQ\tSN:chr1\tLN:100",
		"r_basic\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII",
		"r_softclip\t0\tchr1\t5\t60\t2S6M2S\t*\t0\t0\tNNACGTACNN\tIIIIIIIIII",
		"r_insertion\t0\tchr1\t10\t60\t3M2I3M\t*\t0\t0\tACGTTACG\tIIIIIIII",
		"r_deletion\t0\tchr1\t20\t60\t3M2D3M\t*\t0\t0\tACGTAC\tIIIIII",
		"r_unmapped\t4\t*\t0\t0\t*\t*\t0\t0\tACGTAC\tIIIIII",
		"r_secondary\t256\tchr1\t30\t60\t6M\t*\t0\t0\tACGTAC\tIIIIII",
		"r_lowQ\t0\tchr1\t40\t60\t6M\t*\t0\t0\tACGTAC\t!!!!!!",
		"r_starseq\t0\tchr1\t50\t60\t6M\t*\t0\t0\t*\t*",
	}, "\n") + "\n"

	var buf bytes.Buffer
	n, err := Targetcut(strings.NewReader(samText), &buf, TargetcutOptions{
		MinBaseQ:   DefaultTargetcutMinBaseQ,
		SimpleMode: true,
	})
	if err != nil {
		t.Fatalf("Targetcut: %v", err)
	}
	if n != 4 {
		t.Fatalf("emitted = %d, want 4; full output:\n%s", n, buf.String())
	}

	want := strings.Join([]string{
		">r_basic",
		"ACGTACGTAC",
		">r_softclip",
		"ACGTAC",
		">r_insertion",
		"ACGTTACG",
		">r_deletion",
		"ACGTAC",
		"",
	}, "\n")
	if got := buf.String(); got != want {
		t.Errorf("Targetcut simple output mismatch\ngot:\n%s\nwant:\n%s", got, want)
	}
}

// TestTargetcutSimple_MinBaseQDropsBases checks per-base quality
// filtering in simple mode.
func TestTargetcutSimple_MinBaseQDropsBases(t *testing.T) {
	samText := strings.Join([]string{
		"@HD\tVN:1.6\tSO:coordinate",
		"@SQ\tSN:chr1\tLN:100",
		"r_mixed\t0\tchr1\t1\t60\t6M\t*\t0\t0\tACGTAC\tI!I!I!",
	}, "\n") + "\n"

	var buf bytes.Buffer
	if _, err := Targetcut(strings.NewReader(samText), &buf, TargetcutOptions{
		MinBaseQ:   30,
		SimpleMode: true,
	}); err != nil {
		t.Fatalf("Targetcut: %v", err)
	}
	want := ">r_mixed\nAGA\n"
	if got := buf.String(); got != want {
		t.Errorf("Targetcut min-baseq output:\ngot:  %q\nwant: %q", got, want)
	}
}

// TestTargetcut_DefaultExportedConst guards against drift in the
// upstream-default base-quality cutoff.
func TestTargetcut_DefaultExportedConst(t *testing.T) {
	if DefaultTargetcutMinBaseQ != 13 {
		t.Errorf("DefaultTargetcutMinBaseQ = %d, want 13 (upstream default)", DefaultTargetcutMinBaseQ)
	}
	if DefaultTargetcutEntryPenalty != 14000 {
		t.Errorf("DefaultTargetcutEntryPenalty = %d, want 14000 (upstream default)", DefaultTargetcutEntryPenalty)
	}
}

// ----- HMM-mode tests (the default; faithful port of cut_target.c) -----

// TestTargetcutHMM_SimpleCoverage exercises the HMM consensus mode on
// a uniform high-coverage region. The 2-state Viterbi enters state 1
// for free at position 0, self-loops through the no-coverage prefix
// (emission -4/pos) and the high-coverage core (emission +6/pos),
// then exits to state 0 only when the cumulative state-1 score drops
// below 0. Hand-derived from cut_target.c. Reference length is 60;
// the 20-bp core boost (+120) is exhausted by 30 trailing no-info
// positions (-120) so state 1 exits around pos 30. The exact region
// boundary depends on the upstream tie-break and the back-pointer
// loop that skips position 0 — confirmed against the C source.
func TestTargetcutHMM_SimpleCoverage(t *testing.T) {
	// 8 identical 20M reads at pos 11 (1-based) on a 60-bp reference.
	// Flanks 1..10 and 31..60 are uncovered.
	lines := []string{
		"@HD\tVN:1.6\tSO:coordinate",
		"@SQ\tSN:chr1\tLN:60",
	}
	const seq = "ACGTACGTACGTACGTACGT"
	const qual = "IIIIIIIIIIIIIIIIIIII"
	for i := 0; i < 8; i++ {
		lines = append(lines,
			"r"+targetcutItoa(i)+"\t0\tchr1\t11\t60\t20M\t*\t0\t0\t"+seq+"\t"+qual)
	}
	samText := strings.Join(lines, "\n") + "\n"

	var buf bytes.Buffer
	n, err := Targetcut(strings.NewReader(samText), &buf, TargetcutOptions{})
	if err != nil {
		t.Fatalf("Targetcut HMM: %v", err)
	}
	if n != 1 {
		t.Fatalf("emitted = %d, want 1; output:\n%s", n, buf.String())
	}
	cols := strings.Split(strings.TrimRight(buf.String(), "\n"), "\t")
	if len(cols) != 11 {
		t.Fatalf("expected 11 SAM columns, got %d: %q", len(cols), cols)
	}
	// Upstream's backtrack loop skips position 0, so the region's
	// 0-based start is index 1 (1-based "2"). The region ends at the
	// pos-30 (0-based) transition where state-0 starts winning, which
	// in upstream printf form yields "chr1:2-30" (1-based inclusive end).
	if cols[0] != "chr1:2-30" {
		t.Errorf("QNAME = %q, want chr1:2-30", cols[0])
	}
	if cols[1] != "0" {
		t.Errorf("FLAG = %q, want 0", cols[1])
	}
	if cols[2] != "chr1" {
		t.Errorf("RNAME = %q, want chr1", cols[2])
	}
	if cols[3] != "2" {
		t.Errorf("POS = %q, want 2", cols[3])
	}
	if cols[4] != "60" {
		t.Errorf("MAPQ = %q, want 60", cols[4])
	}
	if cols[5] != "29M" {
		t.Errorf("CIGAR = %q, want 29M", cols[5])
	}
	if cols[6] != "*" {
		t.Errorf("MRNM = %q, want *", cols[6])
	}
	if cols[7] != "0" {
		t.Errorf("MPOS = %q, want 0", cols[7])
	}
	if cols[8] != "0" {
		t.Errorf("ISIZE = %q, want 0", cols[8])
	}
	// SEQ shape: 9 leading 'N' (uncovered prefix from 0-based pos 1..9,
	// which is the trailing 9 of the 10-bp leading flank — position 0
	// was excluded by the upstream backtrack quirk) + 20 callable
	// bases. Trailing flank is not included because the HMM exits to
	// state 0.
	want := strings.Repeat("N", 9) + seq
	if cols[9] != want {
		t.Errorf("SEQ = %q\nwant %q", cols[9], want)
	}
	if len(cols[10]) != 29 {
		t.Errorf("QUAL length = %d, want 29", len(cols[10]))
	}
}

// TestTargetcutHMM_NoCoverage verifies that an entirely empty chrom
// produces no output (the Viterbi never enters state 1).
func TestTargetcutHMM_NoCoverage(t *testing.T) {
	samText := strings.Join([]string{
		"@HD\tVN:1.6\tSO:coordinate",
		"@SQ\tSN:chr1\tLN:100",
	}, "\n") + "\n"
	var buf bytes.Buffer
	n, err := Targetcut(strings.NewReader(samText), &buf, TargetcutOptions{})
	if err != nil {
		t.Fatalf("Targetcut HMM: %v", err)
	}
	if n != 0 {
		t.Fatalf("emitted = %d, want 0; output:\n%s", n, buf.String())
	}
	if buf.Len() != 0 {
		t.Errorf("expected empty output, got %q", buf.String())
	}
}

// TestTargetcutHMM_FiltersOutFlags verifies the upstream read-filter:
// unmapped / secondary / qcfail / duplicate records contribute nothing
// to the consensus.
func TestTargetcutHMM_FiltersOutFlags(t *testing.T) {
	// All "bad" reads at chr1:1 have qual 'I' (passes minBQ). If the
	// filter is broken they'd build a callable region; with the filter
	// applied the output must be empty.
	samText := strings.Join([]string{
		"@HD\tVN:1.6\tSO:coordinate",
		"@SQ\tSN:chr1\tLN:50",
		"r_unmap\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII",
		"r_secondary\t256\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII",
		"r_qcfail\t512\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII",
		"r_dup\t1024\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII",
	}, "\n") + "\n"

	var buf bytes.Buffer
	n, err := Targetcut(strings.NewReader(samText), &buf, TargetcutOptions{})
	if err != nil {
		t.Fatalf("Targetcut HMM: %v", err)
	}
	if n != 0 || buf.Len() != 0 {
		t.Errorf("expected zero records; got %d records: %q", n, buf.String())
	}
}

// TestTargetcutHMM_MinBaseQDropsCallableBase verifies that bumping
// MinBaseQ above every read's per-base quality renders all positions
// "no info" (cell == 0) and the HMM emits nothing.
func TestTargetcutHMM_MinBaseQDropsCallableBase(t *testing.T) {
	// 4 stacked 10M reads with quality 'I' = Phred 40. MinBaseQ = 50
	// excludes every base, so contribs[pos] is empty and gencns
	// returns 0 — cell class 0 (no info) at every position.
	lines := []string{
		"@HD\tVN:1.6\tSO:coordinate",
		"@SQ\tSN:chr1\tLN:40",
	}
	for i := 0; i < 4; i++ {
		lines = append(lines,
			"r"+targetcutItoa(i)+"\t0\tchr1\t5\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII")
	}
	samText := strings.Join(lines, "\n") + "\n"
	var buf bytes.Buffer
	n, err := Targetcut(strings.NewReader(samText), &buf, TargetcutOptions{MinBaseQ: 50})
	if err != nil {
		t.Fatalf("Targetcut HMM: %v", err)
	}
	if n != 0 || buf.Len() != 0 {
		t.Errorf("MinBaseQ=50 should have suppressed all output, got %d records: %q", n, buf.String())
	}
}

// TestTargetcutHMM_SmallEntryPenaltySplitsRegions probes the -i flag.
// With the default entry penalty (14000) the HMM never separates the
// two coverage blocks across a small gap (the inter-block emission
// cost of -4/cell at 80 cells = -320 stays well above -14000). Drop
// the entry penalty to 50 and the HMM will prefer two narrow
// "inside" runs over one wide one — emitting two regions.
//
// Hand-derived from cut_target.c: with -i 50, the cost of entering
// the second region (50) is cheaper than the 80 positions of inside-
// state self-loop with emission -4 each (= 320). So the optimal path
// is state-0 from pos 11..90 instead of state-1.
func TestTargetcutHMM_SmallEntryPenaltySplitsRegions(t *testing.T) {
	lines := []string{
		"@HD\tVN:1.6\tSO:coordinate",
		"@SQ\tSN:chr1\tLN:100",
	}
	for i := 0; i < 8; i++ {
		lines = append(lines,
			"r_left"+targetcutItoa(i)+"\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII")
	}
	for i := 0; i < 8; i++ {
		lines = append(lines,
			"r_right"+targetcutItoa(i)+"\t0\tchr1\t91\t60\t10M\t*\t0\t0\tTGCATGCATG\tIIIIIIIIII")
	}
	samText := strings.Join(lines, "\n") + "\n"

	var buf bytes.Buffer
	n, err := Targetcut(strings.NewReader(samText), &buf, TargetcutOptions{
		EntryPenalty: 50,
	})
	if err != nil {
		t.Fatalf("Targetcut HMM: %v", err)
	}
	if n != 2 {
		t.Fatalf("expected 2 separate regions, got %d; output:\n%s", n, buf.String())
	}
	got := buf.String()
	// First region starts at 0-based pos 1 (upstream backtrack quirk
	// skips pos 0); second region starts at 0-based pos 90.
	if !strings.HasPrefix(got, "chr1:2-10\t") {
		t.Errorf("first region not chr1:2-10; full output:\n%s", got)
	}
	if !strings.Contains(got, "\nchr1:91-100\t") {
		t.Errorf("second region not chr1:91-100; full output:\n%s", got)
	}
}

// TestTargetcutHMM_ConsensusBaseFromMajorityVote checks that gencns
// picks the majority allele when reads disagree at a single position
// inside an otherwise-uniform block.
func TestTargetcutHMM_ConsensusBaseFromMajorityVote(t *testing.T) {
	// 6 reads agree on 'A' at chr1:5, 2 reads vote 'C'. Consensus must
	// be 'A'. Bases at the other 9 positions are uniform.
	lines := []string{
		"@HD\tVN:1.6\tSO:coordinate",
		"@SQ\tSN:chr1\tLN:30",
	}
	// All reads are 10M starting at chr1:1.
	for i := 0; i < 6; i++ {
		lines = append(lines,
			"r_maj"+targetcutItoa(i)+"\t0\tchr1\t1\t60\t10M\t*\t0\t0\tCCCCACCCCC\tIIIIIIIIII")
	}
	for i := 0; i < 2; i++ {
		lines = append(lines,
			"r_min"+targetcutItoa(i)+"\t0\tchr1\t1\t60\t10M\t*\t0\t0\tCCCCCCCCCC\tIIIIIIIIII")
	}
	samText := strings.Join(lines, "\n") + "\n"

	var buf bytes.Buffer
	n, err := Targetcut(strings.NewReader(samText), &buf, TargetcutOptions{})
	if err != nil {
		t.Fatalf("Targetcut HMM: %v", err)
	}
	if n != 1 {
		t.Fatalf("expected 1 region, got %d; output:\n%s", n, buf.String())
	}
	cols := strings.Split(strings.TrimRight(buf.String(), "\n"), "\t")
	if len(cols) < 10 {
		t.Fatalf("malformed SAM line: %q", buf.String())
	}
	// The region starts at 0-based pos 1 (upstream backtrack quirk),
	// so SEQ index k corresponds to reference 0-based pos k+1. The
	// majority-vote position is reference 0-based pos 4, i.e. SEQ
	// index 3.
	seq := cols[9]
	if !strings.HasPrefix(cols[0], "chr1:2-") {
		t.Errorf("region QNAME = %q; want it to start at chr1:2-...", cols[0])
	}
	if len(seq) < 4 || seq[3] != 'A' {
		t.Errorf("consensus at majority-vote position = %q (full SEQ %q); want 'A' from 6-of-8 vote",
			string(seq[3:4]), seq)
	}
}

// TestErrModelSelfConsistent guards the shared errmod port through the
// targetcut call path: when every observed base agrees on one allele
// at qual 40, the diagonal entry for that allele must be the minimum
// (best) score and the other homozygous diagonals must score worse.
func TestErrModelSelfConsistent(t *testing.T) {
	em := errmod.Init(errModDepCorr)
	// 10 bases of 'A' (base = 0), forward strand, quality 40.
	bases := make([]uint16, 10)
	for i := range bases {
		bases[i] = uint16(40)<<5 | uint16(0)<<4 | uint16(0)
	}
	q := make([]float32, 16)
	em.Cal(bases, 4, q)
	for j := 1; j < 4; j++ {
		if q[j*4+j] <= q[0] {
			t.Errorf("homozygous q[%d,%d] = %g should be > q[A,A] = %g (10 A's)",
				j, j, q[j*4+j], q[0])
		}
	}
}

// TestTargetcutHMM_BAQReferenceChangesConsensus exercises the `-f`
// (BAQ) path: with a reference supplied, pkg/htsgo/baq.SamProbRealn
// lowers the qualities of bases flanked by mismatches, which shifts
// gencns' per-position consensus call and/or its quality. The test
// runs the same SAM input twice — once without `-f` and once with —
// and asserts that the per-position QUAL bytes differ. The
// `--simple` mode is unaffected, validated by a separate run.
func TestTargetcutHMM_BAQReferenceChangesConsensus(t *testing.T) {
	// 60-bp reference of all A's. Reads stack 10-deep at chr1:1 and
	// each read carries one of two mismatches near the centre. Without
	// BAQ every base is q40 and the consensus QUAL is uniformly high.
	// With BAQ the mismatch flanks get downweighted, so the consensus
	// QUAL bytes drop at those positions and the output differs.
	const refLen = 60
	refSeq := strings.Repeat("A", refLen)

	dir := t.TempDir()
	faPath := filepath.Join(dir, "ref.fa")
	if err := os.WriteFile(faPath, []byte(">chr1\n"+refSeq+"\n"), 0o644); err != nil {
		t.Fatalf("write ref: %v", err)
	}

	// Build a SAM with 10 stacked 30M reads at chr1:1. Each read has
	// two mismatches at positions 10 and 20 (1-based bases 'C') —
	// enough internal disagreement to trigger BAQ quality reduction
	// in the flanking match positions when a reference is supplied.
	lines := []string{
		"@HD\tVN:1.6\tSO:coordinate",
		"@SQ\tSN:chr1\tLN:60",
	}
	readSeq := strings.Repeat("A", 9) + "C" + strings.Repeat("A", 9) + "C" + strings.Repeat("A", 10)
	readQual := strings.Repeat("I", 30) // Phred 40 everywhere
	for i := 0; i < 10; i++ {
		lines = append(lines,
			"r"+targetcutItoa(i)+"\t0\tchr1\t1\t60\t30M\t*\t0\t0\t"+readSeq+"\t"+readQual)
	}
	samText := strings.Join(lines, "\n") + "\n"

	run := func(fasta string) string {
		var buf bytes.Buffer
		if _, err := Targetcut(strings.NewReader(samText), &buf, TargetcutOptions{
			FastaRef: fasta,
		}); err != nil {
			t.Fatalf("Targetcut(fasta=%q): %v", fasta, err)
		}
		return buf.String()
	}

	without := run("")
	with := run(faPath)
	if without == "" {
		t.Fatal("no-ref run produced no output (test fixture broken)")
	}
	if without == with {
		t.Fatalf("BAQ should have changed at least one output byte but outputs are identical:\n%s",
			without)
	}

	// And simple mode is unaffected by `-f`: same byte output.
	runSimple := func(fasta string) string {
		var buf bytes.Buffer
		if _, err := Targetcut(strings.NewReader(samText), &buf, TargetcutOptions{
			FastaRef:   fasta,
			SimpleMode: true,
		}); err != nil {
			t.Fatalf("Targetcut(simple, fasta=%q): %v", fasta, err)
		}
		return buf.String()
	}
	if runSimple("") != runSimple(faPath) {
		t.Errorf("simple-mode output changed with `-f`; it should be unaffected")
	}
}

// targetcutItoa is a local digit-only helper so the test doesn't
// pull strconv into the test-name format above. Three-digit indices
// suffice for every fixture in this file.
func targetcutItoa(i int) string {
	if i < 10 {
		return string('0' + byte(i))
	}
	if i < 100 {
		return string('0'+byte(i/10)) + string('0'+byte(i%10))
	}
	return string('0'+byte(i/100)) + string('0'+byte((i/10)%10)) + string('0'+byte(i%10))
}
