package seqtk

import (
	"bytes"
	"strings"
	"testing"
)

// TestPair_FastqRoundTrip verifies that "pair" exactly inverts
// "mergepe" on a synthetic interleaved FASTQ stream.
func TestPair_FastqRoundTrip(t *testing.T) {
	r1Want := strings.Join([]string{
		"@a/1", "AAAA", "+", "IIII",
		"@b/1", "BBBB", "+", "JJJJ",
		"@c/1", "CCCC", "+", "KKKK",
		"",
	}, "\n")
	r2Want := strings.Join([]string{
		"@a/2", "tttt", "+", "iiii",
		"@b/2", "uuuu", "+", "jjjj",
		"@c/2", "vvvv", "+", "kkkk",
		"",
	}, "\n")

	// Build the interleaved input by interleaving the wanted outputs.
	var interleaved bytes.Buffer
	if err := MergePE(strings.NewReader(r1Want), strings.NewReader(r2Want), &interleaved); err != nil {
		t.Fatalf("MergePE setup: %v", err)
	}

	var out1, out2 bytes.Buffer
	if err := Pair(&interleaved, &out1, &out2); err != nil {
		t.Fatalf("Pair: %v", err)
	}
	if out1.String() != r1Want {
		t.Errorf("out1 mismatch:\n--want--\n%s\n--got--\n%s", r1Want, out1.String())
	}
	if out2.String() != r2Want {
		t.Errorf("out2 mismatch:\n--want--\n%s\n--got--\n%s", r2Want, out2.String())
	}
}

// TestPair_FastaRoundTrip verifies the same for FASTA inputs.
func TestPair_FastaRoundTrip(t *testing.T) {
	r1Want := strings.Join([]string{
		">a/1", "AAAA",
		">b/1", "BBBB",
		"",
	}, "\n")
	r2Want := strings.Join([]string{
		">a/2", "tttt",
		">b/2", "uuuu",
		"",
	}, "\n")

	var interleaved bytes.Buffer
	if err := MergePE(strings.NewReader(r1Want), strings.NewReader(r2Want), &interleaved); err != nil {
		t.Fatalf("MergePE setup: %v", err)
	}

	var out1, out2 bytes.Buffer
	if err := Pair(&interleaved, &out1, &out2); err != nil {
		t.Fatalf("Pair (fasta): %v", err)
	}
	if out1.String() != r1Want {
		t.Errorf("fasta out1 mismatch:\n--want--\n%s\n--got--\n%s", r1Want, out1.String())
	}
	if out2.String() != r2Want {
		t.Errorf("fasta out2 mismatch:\n--want--\n%s\n--got--\n%s", r2Want, out2.String())
	}
}

// TestPair_OddRecordCount asserts an unpaired trailing record is
// reported as an error rather than silently emitted.
func TestPair_OddRecordCount(t *testing.T) {
	input := strings.Join([]string{
		"@a/1", "AAAA", "+", "IIII",
		"@a/2", "tttt", "+", "iiii",
		"@trailing/1", "CCCC", "+", "KKKK",
		"",
	}, "\n")

	var out1, out2 bytes.Buffer
	err := Pair(strings.NewReader(input), &out1, &out2)
	if err == nil {
		t.Fatalf("expected error on odd record count, got nil")
	}
	if !strings.Contains(err.Error(), "odd record count") {
		t.Errorf("expected 'odd record count' in error, got: %v", err)
	}
}

// TestPair_Empty asserts an empty input is a no-op (no error, no
// output).
func TestPair_Empty(t *testing.T) {
	var out1, out2 bytes.Buffer
	if err := Pair(strings.NewReader(""), &out1, &out2); err != nil {
		t.Fatalf("Pair empty: %v", err)
	}
	if out1.Len() != 0 || out2.Len() != 0 {
		t.Errorf("Pair empty produced output: out1=%q out2=%q", out1.String(), out2.String())
	}
}
