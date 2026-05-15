package samtools

import (
	"bytes"
	"strings"
	"testing"
)

// TestPhase_TwoHetsConsistentChain builds a tiny SAM with two het
// SNP positions and exercises the greedy chaining. Two cohorts of
// three reads each carry different (G,G) vs (T,C) allele
// combinations across the two sites. The chainer should phase both
// hets into one block; whether the second label is 1 or 2 depends on
// the deterministic allele-ordering policy in callHets (ACGT
// iteration + stable-sort by count). Here pos-3 alleles tie at 3:3
// → stable-sort keeps the ACGT order so allele0=G at pos 3; pos-7
// alleles also tie at 3:3 so allele0=C at pos 7. The reads that
// supported allele0 at pos 3 (G-bearing r_a/b/c) support allele1 at
// pos 7 (G), so the chainer flips the label and emits "2".
func TestPhase_TwoHetsConsistentChain(t *testing.T) {
	// Pileup design (chr1 positions 1..10):
	//   pos: 1 2 3 4 5 6 7 8 9 10
	//   ref: A C G T A C G T A C
	// Het sites are at pos 3 (G/T) and pos 7 (G/C).
	// r_a, r_b, r_c: ACGTACGTAC (G at 3, G at 7)
	// r_d, r_e, r_f: ACTTACCTAC (T at 3, C at 7)
	samText := strings.Join([]string{
		"@HD\tVN:1.6\tSO:coordinate",
		"@SQ\tSN:chr1\tLN:100",
		"r_a\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII",
		"r_b\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII",
		"r_c\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII",
		"r_d\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACTTACCTAC\tIIIIIIIIII",
		"r_e\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACTTACCTAC\tIIIIIIIIII",
		"r_f\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACTTACCTAC\tIIIIIIIIII",
	}, "\n") + "\n"

	var buf bytes.Buffer
	n, err := Phase(strings.NewReader(samText), &buf, PhaseOptions{})
	if err != nil {
		t.Fatalf("Phase: %v", err)
	}
	if n != 2 {
		t.Fatalf("emitted = %d, want 2; output:\n%s", n, buf.String())
	}
	want := strings.Join([]string{
		"PS\tchr1\t3\t1",
		"PS\tchr1\t7\t2",
		"",
	}, "\n")
	if got := buf.String(); got != want {
		t.Errorf("Phase output mismatch\ngot:\n%s\nwant:\n%s", got, want)
	}
}

// TestPhase_SecondHetLabelOne mirrors TestPhase_TwoHetsConsistentChain
// but with the second het's allele counts arranged so the chainer
// emits label 1 (no flip). At pos 7 the reads that supported allele0
// (G) at pos 3 carry C at pos 7. callHets's stable sort picks allele0
// = C (ACGT iteration places C before G when both have 3 reads), so
// the same readAssign mapping survives → same=6, opposite=0 → label 1.
func TestPhase_SecondHetLabelOne(t *testing.T) {
	// r_a/b/c: G@3 + C@7
	// r_d/e/f: T@3 + G@7
	samText := strings.Join([]string{
		"@HD\tVN:1.6\tSO:coordinate",
		"@SQ\tSN:chr1\tLN:100",
		"r_a\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACCTAC\tIIIIIIIIII",
		"r_b\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACCTAC\tIIIIIIIIII",
		"r_c\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACCTAC\tIIIIIIIIII",
		"r_d\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACTTACGTAC\tIIIIIIIIII",
		"r_e\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACTTACGTAC\tIIIIIIIIII",
		"r_f\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACTTACGTAC\tIIIIIIIIII",
	}, "\n") + "\n"

	var buf bytes.Buffer
	n, err := Phase(strings.NewReader(samText), &buf, PhaseOptions{})
	if err != nil {
		t.Fatalf("Phase: %v", err)
	}
	if n != 2 {
		t.Fatalf("emitted = %d, want 2; output:\n%s", n, buf.String())
	}
	got := buf.String()
	want := strings.Join([]string{
		"PS\tchr1\t3\t1",
		"PS\tchr1\t7\t1",
		"",
	}, "\n")
	if got != want {
		t.Errorf("Phase output mismatch\ngot:\n%s\nwant:\n%s", got, want)
	}
}

// TestPhase_AmbiguousLabel covers a het site where no reads connect
// it to the previous successfully-phased site. Expected label is 0.
func TestPhase_AmbiguousLabel(t *testing.T) {
	// First het block: reads r_a..r_f all spanning pos 3 (G/T).
	// Second het at pos 50 has its own reads (s_a..s_f) that do NOT
	// span the first het — so the chainer should restart the block
	// and label the first het of the new block 1.
	//
	// Actually: in the v1 implementation, "no overlap" yields label
	// 0 (ambiguous) and increments the consecUnphased counter; once
	// the counter exceeds BlockWindow a new block starts. With the
	// default window of 13, two consecutive hets in different read
	// groups → second one labels 0.
	samText := strings.Join([]string{
		"@HD\tVN:1.6\tSO:coordinate",
		"@SQ\tSN:chr1\tLN:200",
		// Het 1 reads at pos 1..10
		"r_a\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII",
		"r_b\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII",
		"r_c\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACTTACGTAC\tIIIIIIIIII",
		"r_d\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACTTACGTAC\tIIIIIIIIII",
		// Het 2 reads at pos 50..59 — disjoint from het 1 reads.
		"s_a\t0\tchr1\t50\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII",
		"s_b\t0\tchr1\t50\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII",
		"s_c\t0\tchr1\t50\t60\t10M\t*\t0\t0\tACTTACGTAC\tIIIIIIIIII",
		"s_d\t0\tchr1\t50\t60\t10M\t*\t0\t0\tACTTACGTAC\tIIIIIIIIII",
	}, "\n") + "\n"

	var buf bytes.Buffer
	if _, err := Phase(strings.NewReader(samText), &buf, PhaseOptions{}); err != nil {
		t.Fatalf("Phase: %v", err)
	}
	got := buf.String()
	// Het at pos 3: first het in block → label 1.
	// Het at pos 52: no read overlaps with het 1's reads → label 0.
	want := strings.Join([]string{
		"PS\tchr1\t3\t1",
		"PS\tchr1\t52\t0",
		"",
	}, "\n")
	if got != want {
		t.Errorf("Phase ambiguous output mismatch\ngot:\n%s\nwant:\n%s", got, want)
	}
}

// TestPhase_LowMAPQSkipped verifies the -q filter excludes records
// before they contribute to het calling.
func TestPhase_LowMAPQSkipped(t *testing.T) {
	// Two reads have allele G, two have allele T — but the T ones
	// are MAPQ=5 (below default 13) so they are dropped. With only
	// one allele surviving, no het is called → no output.
	samText := strings.Join([]string{
		"@HD\tVN:1.6\tSO:coordinate",
		"@SQ\tSN:chr1\tLN:100",
		"r_a\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII",
		"r_b\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII",
		"r_c\t0\tchr1\t1\t5\t10M\t*\t0\t0\tACTTACGTAC\tIIIIIIIIII",
		"r_d\t0\tchr1\t1\t5\t10M\t*\t0\t0\tACTTACGTAC\tIIIIIIIIII",
	}, "\n") + "\n"

	var buf bytes.Buffer
	n, err := Phase(strings.NewReader(samText), &buf, PhaseOptions{MinMAPQ: DefaultPhaseMinMAPQ})
	if err != nil {
		t.Fatalf("Phase: %v", err)
	}
	if n != 0 {
		t.Errorf("emitted = %d, want 0 (all T-bearing reads MAPQ-filtered); output:\n%s", n, buf.String())
	}
}

// TestPhase_DefaultExportedConsts guards against drift in the
// upstream-default parameter values.
func TestPhase_DefaultExportedConsts(t *testing.T) {
	if DefaultPhaseBlockWindow != 13 {
		t.Errorf("DefaultPhaseBlockWindow = %d, want 13", DefaultPhaseBlockWindow)
	}
	if DefaultPhaseMinMAPQ != 13 {
		t.Errorf("DefaultPhaseMinMAPQ = %d, want 13", DefaultPhaseMinMAPQ)
	}
	if DefaultPhaseMinBaseQ != 13 {
		t.Errorf("DefaultPhaseMinBaseQ = %d, want 13", DefaultPhaseMinBaseQ)
	}
	if DefaultPhaseMaxDepth != 256 {
		t.Errorf("DefaultPhaseMaxDepth = %d, want 256", DefaultPhaseMaxDepth)
	}
}
