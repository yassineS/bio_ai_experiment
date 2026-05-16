package seqtk

import (
	"bytes"
	"strings"
	"testing"
)

// makeTeloFASTA stitches together a synthetic FASTA record with a 5'
// telomere, middle filler, and 3' reverse-complement telomere. The
// defaults (60 × 6 = 360 bp at each end with default penalty=1) keep
// the max running score above the default min_score=300 so the
// produced fixture exercises both 5' and 3' BED rows.
func makeTeloFASTA(name string, fivePrime, mid, threePrime string, fiveRepeats, midLen, threeRepeats int) string {
	var b strings.Builder
	b.WriteString(">")
	b.WriteString(name)
	b.WriteByte('\n')
	for i := 0; i < fiveRepeats; i++ {
		b.WriteString(fivePrime)
	}
	for i := 0; i < midLen; i++ {
		b.WriteByte(mid[i%len(mid)])
	}
	for i := 0; i < threeRepeats; i++ {
		b.WriteString(threePrime)
	}
	b.WriteByte('\n')
	return b.String()
}

// TestTelo_DefaultBothEnds asserts the two-row BED output for a
// synthetic chromosome with CCCTAA at the 5' end and TTAGGG (the
// reverse complement) at the 3' end. With default options the 5' scan
// emits "<name>\t0\t<end>\t<len>" and the 3' scan emits
// "<name>\t<start>\t<len>\t<len>".
func TestTelo_DefaultBothEnds(t *testing.T) {
	in := makeTeloFASTA("chr1", "CCCTAA", "ACGTACGT", "TTAGGG", 60, 800, 60)
	var stdout, stderr bytes.Buffer
	opts := TeloOptions{
		Motif:    DefaultTeloMotif,
		Penalty:  DefaultTeloPenalty,
		MaxDrop:  DefaultTeloMaxDrop,
		MinScore: DefaultTeloMinScore,
	}
	if err := Telo(strings.NewReader(in), &stdout, &stderr, opts); err != nil {
		t.Fatalf("Telo: %v", err)
	}
	wantStdout := "chr1\t0\t360\t1520\nchr1\t1160\t1520\t1520\n"
	if stdout.String() != wantStdout {
		t.Fatalf("stdout mismatch\ngot:\n%s\nwant:\n%s", stdout.String(), wantStdout)
	}
	wantStderr := "720\t1520\n"
	if stderr.String() != wantStderr {
		t.Fatalf("stderr mismatch\ngot:\n%s\nwant:\n%s", stderr.String(), wantStderr)
	}
}

// TestTelo_TelomereAtStartOnly covers the edge case where only the
// 5' end has a telomere and the 3' scan should not emit anything.
// This pins the `st = max_i + 1` cursor so the 3' scan walks back
// only as far as the start of the 5' telomere.
func TestTelo_TelomereAtStartOnly(t *testing.T) {
	in := makeTeloFASTA("five_only", "CCCTAA", "ACGTACGT", "ACGTACGT", 60, 800, 0)
	var stdout, stderr bytes.Buffer
	opts := TeloOptions{
		Motif:    DefaultTeloMotif,
		Penalty:  DefaultTeloPenalty,
		MaxDrop:  DefaultTeloMaxDrop,
		MinScore: DefaultTeloMinScore,
	}
	if err := Telo(strings.NewReader(in), &stdout, &stderr, opts); err != nil {
		t.Fatalf("Telo: %v", err)
	}
	wantStdout := "five_only\t0\t360\t1160\n"
	if stdout.String() != wantStdout {
		t.Fatalf("stdout mismatch\ngot:\n%s\nwant:\n%s", stdout.String(), wantStdout)
	}
	wantStderr := "360\t1160\n"
	if stderr.String() != wantStderr {
		t.Fatalf("stderr mismatch\ngot:\n%s\nwant:\n%s", stderr.String(), wantStderr)
	}
}

// TestTelo_TelomereAtEndOnly covers the symmetric case: filler at the
// 5' end, telomere at the 3' end. The 5' scan should produce nothing.
func TestTelo_TelomereAtEndOnly(t *testing.T) {
	in := makeTeloFASTA("three_only", "ACGTACGT", "ACGTACGT", "TTAGGG", 0, 800, 60)
	var stdout, stderr bytes.Buffer
	opts := TeloOptions{
		Motif:    DefaultTeloMotif,
		Penalty:  DefaultTeloPenalty,
		MaxDrop:  DefaultTeloMaxDrop,
		MinScore: DefaultTeloMinScore,
	}
	if err := Telo(strings.NewReader(in), &stdout, &stderr, opts); err != nil {
		t.Fatalf("Telo: %v", err)
	}
	wantStdout := "three_only\t800\t1160\t1160\n"
	if stdout.String() != wantStdout {
		t.Fatalf("stdout mismatch\ngot:\n%s\nwant:\n%s", stdout.String(), wantStdout)
	}
	wantStderr := "360\t1160\n"
	if stderr.String() != wantStderr {
		t.Fatalf("stderr mismatch\ngot:\n%s\nwant:\n%s", stderr.String(), wantStderr)
	}
}

// TestTelo_CustomMotifSwapsStrands feeds a sequence whose 5' end has
// TTAGGG and 3' end has CCCTAA, then passes -m TTAGGG. Upstream's
// 3' scan checks the reverse-complement-encoded motif rotations using
// the same hash set, so swapping the motif effectively swaps which
// end matches forward and which matches reverse.
func TestTelo_CustomMotifSwapsStrands(t *testing.T) {
	in := makeTeloFASTA("chr1", "TTAGGG", "ACGTACGT", "CCCTAA", 60, 800, 60)
	var stdout, stderr bytes.Buffer
	opts := TeloOptions{
		Motif:    "TTAGGG",
		Penalty:  DefaultTeloPenalty,
		MaxDrop:  DefaultTeloMaxDrop,
		MinScore: DefaultTeloMinScore,
	}
	if err := Telo(strings.NewReader(in), &stdout, &stderr, opts); err != nil {
		t.Fatalf("Telo: %v", err)
	}
	wantStdout := "chr1\t0\t360\t1520\nchr1\t1160\t1520\t1520\n"
	if stdout.String() != wantStdout {
		t.Fatalf("stdout mismatch\ngot:\n%s\nwant:\n%s", stdout.String(), wantStdout)
	}
}

// TestTelo_ProfileMode runs the same FASTA in -P mode and asserts the
// first / last profile rows match the expected shape. Full byte-parity
// against upstream is in TestParity_Seqtk_Telo_Profile.
func TestTelo_ProfileMode(t *testing.T) {
	in := makeTeloFASTA("chr1", "CCCTAA", "ACGTACGT", "TTAGGG", 60, 800, 60)
	var stdout, stderr bytes.Buffer
	opts := TeloOptions{
		Motif:       DefaultTeloMotif,
		Penalty:     DefaultTeloPenalty,
		MaxDrop:     DefaultTeloMaxDrop,
		MinScore:    DefaultTeloMinScore,
		ShowProfile: true,
	}
	if err := Telo(strings.NewReader(in), &stdout, &stderr, opts); err != nil {
		t.Fatalf("Telo: %v", err)
	}
	lines := strings.Split(strings.TrimRight(stdout.String(), "\n"), "\n")
	if len(lines) < 2 {
		t.Fatalf("expected many profile lines, got %d", len(lines))
	}
	if !strings.HasPrefix(lines[0], "P\tchr1\t0\t0\t0") {
		t.Fatalf("first profile line shape unexpected: %q", lines[0])
	}
	// The very last line must be a Q row.
	if !strings.HasPrefix(lines[len(lines)-1], "Q\tchr1\t") {
		t.Fatalf("last profile line shape unexpected: %q", lines[len(lines)-1])
	}
}

// TestTelo_EmptyAndShortRecords pins the early-exit behaviour for
// records that cannot reach the default min_score. The 5'/3' loops
// must not produce BED rows, but the summary still accumulates input
// length (sum_telo stays at zero).
func TestTelo_EmptyAndShortRecords(t *testing.T) {
	const in = ">empty\n\n>short\nACGT\n"
	var stdout, stderr bytes.Buffer
	opts := TeloOptions{
		Motif:    DefaultTeloMotif,
		Penalty:  DefaultTeloPenalty,
		MaxDrop:  DefaultTeloMaxDrop,
		MinScore: DefaultTeloMinScore,
	}
	if err := Telo(strings.NewReader(in), &stdout, &stderr, opts); err != nil {
		t.Fatalf("Telo: %v", err)
	}
	if stdout.Len() != 0 {
		t.Fatalf("expected empty stdout for sub-min-score records, got %q", stdout.String())
	}
	if stderr.String() != "0\t4\n" {
		t.Fatalf("stderr summary mismatch: got %q", stderr.String())
	}
}

// TestTelo_NegativePenaltyFlipped pins upstream's
// `if (penalty < 0) penalty = -penalty;` quirk: a negative -p value
// is silently flipped to its absolute value. Without the flip, every
// non-hit base would INCREASE the score and any input would emit a
// telomere row.
func TestTelo_NegativePenaltyFlipped(t *testing.T) {
	in := makeTeloFASTA("chr1", "CCCTAA", "ACGTACGT", "TTAGGG", 60, 800, 60)
	var stdoutPos, stderrPos bytes.Buffer
	var stdoutNeg, stderrNeg bytes.Buffer
	opts := TeloOptions{
		Motif:    DefaultTeloMotif,
		Penalty:  1,
		MaxDrop:  DefaultTeloMaxDrop,
		MinScore: DefaultTeloMinScore,
	}
	if err := Telo(strings.NewReader(in), &stdoutPos, &stderrPos, opts); err != nil {
		t.Fatalf("Telo +1: %v", err)
	}
	opts.Penalty = -1
	if err := Telo(strings.NewReader(in), &stdoutNeg, &stderrNeg, opts); err != nil {
		t.Fatalf("Telo -1: %v", err)
	}
	if stdoutPos.String() != stdoutNeg.String() || stderrPos.String() != stderrNeg.String() {
		t.Fatalf("negative-penalty path diverges from positive:\npos stdout=%q\nneg stdout=%q",
			stdoutPos.String(), stdoutNeg.String())
	}
}

// TestTelo_InvalidMotif rejects empty and non-ACGT motifs with a
// clean error (upstream's assert() would abort).
func TestTelo_InvalidMotif(t *testing.T) {
	cases := []string{"", "NNNN", "ACNT", "RYSW"}
	for _, m := range cases {
		var stdout, stderr bytes.Buffer
		err := Telo(strings.NewReader(">x\nACGT\n"), &stdout, &stderr, TeloOptions{Motif: m})
		if err == nil {
			t.Fatalf("expected error for motif %q", m)
		}
	}
}
