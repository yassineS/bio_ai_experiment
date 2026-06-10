package main

import (
	"os"
	"path/filepath"
	"testing"
)

// TestCLI_AcceptsCRAM verifies the CLI routes a `.cram` path to the CRAM
// decoder (no longer the old "not supported" rejection) and produces the
// expected per-interval counts. It uses the real htslib-produced a.cram
// fixture (shared with the package-level parity tests): two mapped chr1
// alignments around chr1:10003..10143.
func TestCLI_AcceptsCRAM(t *testing.T) {
	cramSrc, err := filepath.Abs(filepath.Join("..", "..", "testdata", "parity", "a.cram"))
	if err != nil {
		t.Fatalf("abs cram path: %v", err)
	}
	data, err := os.ReadFile(cramSrc)
	if err != nil {
		t.Fatalf("read a.cram fixture: %v", err)
	}
	dir := t.TempDir()
	cramPath := filepath.Join(dir, "a.cram")
	if err := os.WriteFile(cramPath, data, 0o644); err != nil {
		t.Fatalf("stage cram: %v", err)
	}

	bedPath := filepath.Join(dir, "a.bed")
	if err := os.WriteFile(bedPath,
		[]byte("chr1\t10000\t10100\tregion1\nchr1\t9000\t9500\tregion2\n"), 0o644); err != nil {
		t.Fatalf("write bed: %v", err)
	}
	outPath := filepath.Join(dir, "out.txt")

	if err := run([]string{"-bed", bedPath, "-bams", cramPath, "-o", outPath}, os.Stdout, os.Stderr); err != nil {
		t.Fatalf("run: %v", err)
	}
	got, err := os.ReadFile(outPath)
	if err != nil {
		t.Fatalf("read out: %v", err)
	}
	// Both reads overlap region1 (chr1:10000-10100); neither overlaps region2.
	want := "chr1\t10000\t10100\tregion1\t2\nchr1\t9000\t9500\tregion2\t0\n"
	if string(got) != want {
		t.Fatalf("CLI CRAM output mismatch:\n got: %q\nwant: %q", string(got), want)
	}
}
