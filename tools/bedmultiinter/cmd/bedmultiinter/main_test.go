package main

import (
	"os"
	"path/filepath"
	"testing"
)

// writeFile stages a file under dir and returns its path.
func writeFile(t *testing.T, dir, name, content string) string {
	t.Helper()
	p := filepath.Join(dir, name)
	if err := os.WriteFile(p, []byte(content), 0o644); err != nil {
		t.Fatalf("write %s: %v", name, err)
	}
	return p
}

// TestCLI_VariadicNames verifies that `-names` consumes a space-separated list
// (one per `-i` file), matching upstream's `-names A B C` syntax rather than a
// comma-separated single token. The header row must carry the supplied labels
// and the `list` column must use them.
func TestCLI_VariadicNames(t *testing.T) {
	dir := t.TempDir()
	a := writeFile(t, dir, "a.bed", "chr1\t10\t20\n")
	b := writeFile(t, dir, "b.bed", "chr1\t15\t25\n")
	out := filepath.Join(dir, "out.txt")

	if err := run([]string{"-i", a, b, "-header", "-names", "AAA", "BBB", "-o", out}, os.Stdout, os.Stderr); err != nil {
		t.Fatalf("run: %v", err)
	}
	got, err := os.ReadFile(out)
	if err != nil {
		t.Fatalf("read out: %v", err)
	}
	want := "chrom\tstart\tend\tnum\tlist\tAAA\tBBB\n" +
		"chr1\t10\t15\t1\tAAA\t1\t0\n" +
		"chr1\t15\t20\t2\tAAA,BBB\t1\t1\n" +
		"chr1\t20\t25\t1\tBBB\t0\t1\n"
	if string(got) != want {
		t.Fatalf("variadic -names output mismatch:\n got: %q\nwant: %q", string(got), want)
	}
}

// TestCLI_NamesCountMismatch verifies the upstream-style error when the number
// of -names labels does not match the number of -i files.
func TestCLI_NamesCountMismatch(t *testing.T) {
	dir := t.TempDir()
	a := writeFile(t, dir, "a.bed", "chr1\t10\t20\n")
	b := writeFile(t, dir, "b.bed", "chr1\t15\t25\n")

	err := run([]string{"-i", a, b, "-names", "ONLYONE"}, os.Stdout, os.Stderr)
	if err == nil {
		t.Fatalf("expected error for mismatched -names count")
	}
}

// TestCLI_NamesBeforeFiles verifies the two variadic flags are order
// independent: `-names` may appear before `-i`.
func TestCLI_NamesBeforeFiles(t *testing.T) {
	dir := t.TempDir()
	a := writeFile(t, dir, "a.bed", "chr1\t10\t20\n")
	b := writeFile(t, dir, "b.bed", "chr1\t15\t25\n")
	out := filepath.Join(dir, "out.txt")

	if err := run([]string{"-names", "X", "Y", "-i", a, b, "-o", out}, os.Stdout, os.Stderr); err != nil {
		t.Fatalf("run: %v", err)
	}
	got, err := os.ReadFile(out)
	if err != nil {
		t.Fatalf("read out: %v", err)
	}
	want := "chr1\t10\t15\t1\tX\t1\t0\n" +
		"chr1\t15\t20\t2\tX,Y\t1\t1\n" +
		"chr1\t20\t25\t1\tY\t0\t1\n"
	if string(got) != want {
		t.Fatalf("names-before-files output mismatch:\n got: %q\nwant: %q", string(got), want)
	}
}
