package seqtk

import (
	"bytes"
	"strings"
	"testing"
)

// TestSamePairName covers the upstream-compatible name-equality
// predicate used by Dropse: identical bytes modulo a trailing
// "/<digit>" suffix.
func TestSamePairName(t *testing.T) {
	cases := []struct {
		name string
		p, q string
		want bool
	}{
		{"identical", "read1", "read1", true},
		{"slash 1 vs 2", "read1/1", "read1/2", true},
		{"slash same digit", "read1/1", "read1/1", true},
		{"different stem", "read1/1", "read2/2", false},
		{"different length", "read1", "read12", false},
		{"no slash mismatch", "read1A", "read1B", false},
		{"min length: 3 bytes with /digit", "x/1", "x/2", true},
		{"too short: 2 bytes /digit", "/1", "/2", false}, // l > 2 fails
		{"underscore not honoured", "read1_1", "read1_2", false},
		{"slash with non-digit", "read/A", "read/B", false},
		{"empty", "", "", true},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			got := samePairName(tc.p, tc.q)
			if got != tc.want {
				t.Errorf("samePairName(%q, %q) = %v, want %v", tc.p, tc.q, got, tc.want)
			}
		})
	}
}

// TestDropse_FastqDropsOrphans verifies that an interleaved FASTQ with
// orphan reads emits only the genuinely paired records.
func TestDropse_FastqDropsOrphans(t *testing.T) {
	input := strings.Join([]string{
		"@read1/1", "ACGT", "+", "IIII",
		"@read1/2", "TGCA", "+", "IIII",
		"@orphan/1", "NNNN", "+", "!!!!",
		"@read2/1", "AAAA", "+", "####",
		"@read2/2", "TTTT", "+", "$$$$",
		"@trail/1", "GGGG", "+", "%%%%",
		"",
	}, "\n")
	want := strings.Join([]string{
		"@read1/1", "ACGT", "+", "IIII",
		"@read1/2", "TGCA", "+", "IIII",
		"@read2/1", "AAAA", "+", "####",
		"@read2/2", "TTTT", "+", "$$$$",
		"",
	}, "\n")

	var out bytes.Buffer
	if err := Dropse(strings.NewReader(input), &out); err != nil {
		t.Fatalf("Dropse: %v", err)
	}
	if got := out.String(); got != want {
		t.Errorf("Dropse output:\n--want--\n%s\n--got--\n%s", want, got)
	}
}

// TestDropse_FastaDropsOrphans is the same test on a FASTA input.
func TestDropse_FastaDropsOrphans(t *testing.T) {
	input := strings.Join([]string{
		">read1/1", "ACGT",
		">read1/2", "TGCA",
		">orphan/1", "NNNN",
		">read2/1", "AAAA",
		">read2/2", "TTTT",
		">trail/1", "GGGG",
		"",
	}, "\n")
	want := strings.Join([]string{
		">read1/1", "ACGT",
		">read1/2", "TGCA",
		">read2/1", "AAAA",
		">read2/2", "TTTT",
		"",
	}, "\n")

	var out bytes.Buffer
	if err := Dropse(strings.NewReader(input), &out); err != nil {
		t.Fatalf("Dropse (fasta): %v", err)
	}
	if got := out.String(); got != want {
		t.Errorf("Dropse fasta output:\n--want--\n%s\n--got--\n%s", want, got)
	}
}

// TestDropse_AllPaired verifies a fully-paired interleaved file is
// emitted unchanged.
func TestDropse_AllPaired(t *testing.T) {
	input := strings.Join([]string{
		"@a/1", "ACGT", "+", "IIII",
		"@a/2", "TGCA", "+", "IIII",
		"@b/1", "AAAA", "+", "####",
		"@b/2", "TTTT", "+", "$$$$",
		"",
	}, "\n")
	var out bytes.Buffer
	if err := Dropse(strings.NewReader(input), &out); err != nil {
		t.Fatalf("Dropse: %v", err)
	}
	if got := out.String(); got != input {
		t.Errorf("Dropse should pass paired input through unchanged.\nwant:\n%s\ngot:\n%s", input, got)
	}
}

// TestDropse_EmptyInput verifies an empty stream produces no output.
func TestDropse_EmptyInput(t *testing.T) {
	var out bytes.Buffer
	if err := Dropse(strings.NewReader(""), &out); err != nil {
		t.Fatalf("Dropse on empty: %v", err)
	}
	if out.Len() != 0 {
		t.Errorf("Dropse on empty input emitted %q", out.String())
	}
}

// TestDropse_RunOfOrphans ensures consecutive non-matching records
// only keep the most recent as "last" (mirroring upstream's
// replacement behaviour).
func TestDropse_RunOfOrphans(t *testing.T) {
	// a, b, c, c/2: only c+c/2 should be emitted.
	input := strings.Join([]string{
		">a/1", "AAA",
		">b/1", "BBB",
		">c/1", "CCC",
		">c/2", "ccc",
		"",
	}, "\n")
	want := strings.Join([]string{
		">c/1", "CCC",
		">c/2", "ccc",
		"",
	}, "\n")
	var out bytes.Buffer
	if err := Dropse(strings.NewReader(input), &out); err != nil {
		t.Fatalf("Dropse: %v", err)
	}
	if got := out.String(); got != want {
		t.Errorf("run-of-orphans:\n--want--\n%s\n--got--\n%s", want, got)
	}
}
