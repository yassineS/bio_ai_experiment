package seqtk

import (
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// readSplitFile reads <dir>/<prefix>.<5-digit i>.fa and returns its
// contents as a string. Fails the test if the file is missing.
func readSplitFile(t *testing.T, dir, prefix string, i int) string {
	t.Helper()
	path := filepath.Join(dir, fmt.Sprintf("%s.%05d.fa", prefix, i))
	b, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("read split output %s: %v", path, err)
	}
	return string(b)
}

// TestSplit_FastaRoundRobin verifies the upstream round-robin
// distribution: record i (0-based) goes to file <prefix>.<i%n+1>.fa.
func TestSplit_FastaRoundRobin(t *testing.T) {
	in := strings.Join([]string{
		">a", "AAAA",
		">b", "BBBB",
		">c", "CCCC",
		">d", "DDDD",
		">e", "EEEE",
		"",
	}, "\n")
	dir := t.TempDir()
	prefix := filepath.Join(dir, "x")
	opts := SplitOptions{N: 2, Prefix: prefix}
	if err := Split(strings.NewReader(in), opts); err != nil {
		t.Fatalf("Split: %v", err)
	}
	// File 1: indices 0, 2, 4  -> a, c, e
	want1 := ">a\nAAAA\n>c\nCCCC\n>e\nEEEE\n"
	want2 := ">b\nBBBB\n>d\nDDDD\n"
	if got := readSplitFile(t, dir, "x", 1); got != want1 {
		t.Errorf("file 1:\nwant:\n%s\ngot:\n%s", want1, got)
	}
	if got := readSplitFile(t, dir, "x", 2); got != want2 {
		t.Errorf("file 2:\nwant:\n%s\ngot:\n%s", want2, got)
	}
}

// TestSplit_FastqPreservesFormat verifies that FASTQ input produces
// FASTQ output (despite the ".fa" file-name suffix, which upstream
// keeps even for FASTQ input).
func TestSplit_FastqPreservesFormat(t *testing.T) {
	in := strings.Join([]string{
		"@r1", "ACGT", "+", "IIII",
		"@r2", "GGGG", "+", "####",
		"",
	}, "\n")
	dir := t.TempDir()
	prefix := filepath.Join(dir, "fq")
	opts := SplitOptions{N: 2, Prefix: prefix}
	if err := Split(strings.NewReader(in), opts); err != nil {
		t.Fatalf("Split fastq: %v", err)
	}
	want1 := "@r1\nACGT\n+\nIIII\n"
	want2 := "@r2\nGGGG\n+\n####\n"
	if got := readSplitFile(t, dir, "fq", 1); got != want1 {
		t.Errorf("fastq file 1:\nwant:\n%s\ngot:\n%s", want1, got)
	}
	if got := readSplitFile(t, dir, "fq", 2); got != want2 {
		t.Errorf("fastq file 2:\nwant:\n%s\ngot:\n%s", want2, got)
	}
}

// TestSplit_LineLengthWrapsFasta verifies the -l flag wraps the
// sequence at the requested width.
func TestSplit_LineLengthWrapsFasta(t *testing.T) {
	in := ">a\nACGTACGTAC\n"
	dir := t.TempDir()
	prefix := filepath.Join(dir, "w")
	opts := SplitOptions{N: 1, LineLen: 4, Prefix: prefix}
	if err := Split(strings.NewReader(in), opts); err != nil {
		t.Fatalf("Split: %v", err)
	}
	want := ">a\nACGT\nACGT\nAC\n"
	if got := readSplitFile(t, dir, "w", 1); got != want {
		t.Errorf("wrapped fasta:\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestSplit_LineLengthWrapsFastq verifies that -l also wraps quality
// lines for FASTQ input (matching upstream's stk_printseq behaviour
// of calling stk_printstr on both seq and qual with the same len).
func TestSplit_LineLengthWrapsFastq(t *testing.T) {
	in := "@r\nACGTACGT\n+\nIIIIIIII\n"
	dir := t.TempDir()
	prefix := filepath.Join(dir, "wq")
	opts := SplitOptions{N: 1, LineLen: 4, Prefix: prefix}
	if err := Split(strings.NewReader(in), opts); err != nil {
		t.Fatalf("Split: %v", err)
	}
	want := "@r\nACGT\nACGT\n+\nIIII\nIIII\n"
	if got := readSplitFile(t, dir, "wq", 1); got != want {
		t.Errorf("wrapped fastq:\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestSplit_AlwaysCreatesNFiles verifies that all N output files are
// created up front, even ones that end up empty (matching upstream's
// "open all N before reading any records" behaviour).
func TestSplit_AlwaysCreatesNFiles(t *testing.T) {
	in := ">solo\nAAAA\n"
	dir := t.TempDir()
	prefix := filepath.Join(dir, "e")
	opts := SplitOptions{N: 3, Prefix: prefix}
	if err := Split(strings.NewReader(in), opts); err != nil {
		t.Fatalf("Split: %v", err)
	}
	for i := 1; i <= 3; i++ {
		path := filepath.Join(dir, fmt.Sprintf("e.%05d.fa", i))
		if _, err := os.Stat(path); err != nil {
			t.Errorf("expected file %s to exist: %v", path, err)
		}
	}
	if got := readSplitFile(t, dir, "e", 1); got != ">solo\nAAAA\n" {
		t.Errorf("file 1 = %q", got)
	}
	if got := readSplitFile(t, dir, "e", 2); got != "" {
		t.Errorf("file 2 should be empty, got %q", got)
	}
	if got := readSplitFile(t, dir, "e", 3); got != "" {
		t.Errorf("file 3 should be empty, got %q", got)
	}
}

// TestSplit_RejectsBadOptions covers the input-validation paths.
func TestSplit_RejectsBadOptions(t *testing.T) {
	if err := Split(strings.NewReader(""), SplitOptions{N: 0, Prefix: "p"}); err == nil {
		t.Error("Split with N=0 should return an error")
	}
	if err := Split(strings.NewReader(""), SplitOptions{N: 1, Prefix: ""}); err == nil {
		t.Error("Split with empty prefix should return an error")
	}
}
