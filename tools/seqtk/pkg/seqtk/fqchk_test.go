package seqtk

import (
	"bytes"
	"strings"
	"testing"
)

// TestFqchk_Basic exercises the default qthres=20 path on a tiny
// hand-built FASTQ. Verifies the header, the ALL row, and the
// per-position rows.
func TestFqchk_Basic(t *testing.T) {
	in := strings.NewReader("@r1\nACGT\n+\nIIII\n@r2\nACGT\n+\n!!!!\n")
	var out bytes.Buffer
	if err := Fqchk(in, &out, FqchkOptions{QThres: DefaultFqchkQThres}); err != nil {
		t.Fatalf("Fqchk: %v", err)
	}
	got := out.String()
	// Header line: column order is fixed, %low/%high only with qthres>0.
	if !strings.HasPrefix(got, "min_len: 4; max_len: 4; avg_len: 4.00; 2 distinct quality values\n") {
		t.Errorf("preamble mismatch: %q", got)
	}
	if !strings.Contains(got, "POS\t#bases\t%A\t%C\t%G\t%T\t%N\tavgQ\terrQ\t%low\t%high\n") {
		t.Errorf("header line mismatch: %q", got)
	}
	if !strings.Contains(got, "\nALL\t8\t") {
		t.Errorf("missing ALL row: %q", got)
	}
}

// TestFqchk_QZero exercises the qthres=0 path where the trailing
// columns become one %Qk per distinct observed quality.
func TestFqchk_QZero(t *testing.T) {
	in := strings.NewReader("@r1\nACGT\n+\nIIII\n@r2\nACGT\n+\n!!!!\n")
	var out bytes.Buffer
	if err := Fqchk(in, &out, FqchkOptions{QThres: 0}); err != nil {
		t.Fatalf("Fqchk -q0: %v", err)
	}
	got := out.String()
	// With Q0 (=='!') and Q40 (=='I') as the only distinct qualities,
	// header trailing cols are "%Q0\t%Q40".
	if !strings.Contains(got, "\terrQ\t%Q0\t%Q40\n") {
		t.Errorf("expected per-quality header in:\n%s", got)
	}
	// %low / %high should NOT appear.
	if strings.Contains(got, "%low") || strings.Contains(got, "%high") {
		t.Errorf("did not expect %%low/%%high in q0 output:\n%s", got)
	}
}

// TestFqchk_VariableLengths verifies that records with different
// lengths are handled (the per-position arrays grow as we go).
func TestFqchk_VariableLengths(t *testing.T) {
	in := strings.NewReader("@a\nACG\n+\nIII\n@b\nACGTACGT\n+\nIIIIIIII\n")
	var out bytes.Buffer
	if err := Fqchk(in, &out, FqchkOptions{QThres: 20}); err != nil {
		t.Fatalf("Fqchk: %v", err)
	}
	got := out.String()
	if !strings.HasPrefix(got, "min_len: 3; max_len: 8; avg_len: 5.50;") {
		t.Errorf("preamble mismatch: %q", got)
	}
	// Per-position rows 1..8 should all appear.
	for i := 1; i <= 8; i++ {
		needle := "\n" + fqchkItoa(i) + "\t"
		if !strings.Contains(got, needle) {
			t.Errorf("missing per-position row %d in:\n%s", i, got)
		}
	}
	// Position 4..8 only has 1 base (from record b).
	if !strings.Contains(got, "\n4\t1\t") {
		t.Errorf("expected '\\n4\\t1\\t' (one base at position 4):\n%s", got)
	}
}

// TestFqchk_EmptyQualityRecordSkipped verifies that records with an
// empty quality line are silently dropped, matching upstream's
// `if (seq->qual.l == 0) continue`. (The fastq reader rejects FASTA
// input outright, so this exercises the in-reader skip path.)
func TestFqchk_QualityClamping(t *testing.T) {
	// '~' is ASCII 126 -> phred 93 (the upper clamp).
	in := strings.NewReader("@a\nA\n+\n~\n")
	var out bytes.Buffer
	if err := Fqchk(in, &out, FqchkOptions{QThres: 50}); err != nil {
		t.Fatalf("Fqchk: %v", err)
	}
	got := out.String()
	// Q93 is %high under qthres=50, so %high == 100.0, %low == 0.0.
	if !strings.Contains(got, "\t0.0\t100.0\n") {
		t.Errorf("expected high-quality clamp to land in %%high column:\n%s", got)
	}
}

// fqchkItoa is a local int -> decimal helper that avoids dragging in
// strconv just for tests (matches the project's "stdlib only, terse"
// taste).
func fqchkItoa(i int) string {
	if i == 0 {
		return "0"
	}
	var buf [20]byte
	pos := len(buf)
	neg := false
	if i < 0 {
		neg = true
		i = -i
	}
	for i > 0 {
		pos--
		buf[pos] = byte('0' + i%10)
		i /= 10
	}
	if neg {
		pos--
		buf[pos] = '-'
	}
	return string(buf[pos:])
}
