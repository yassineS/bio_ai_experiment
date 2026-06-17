package main

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

// TestCLI_IdentifiesMultipleFiles compiles the htsfile binary and runs
// it against two synthetic files (a plain VCF and a plain FASTA),
// asserting the per-file one-line summary lands on stdout in the
// expected `path: <description>` form. End-to-end coverage of the
// argument loop in main.go that the unit tests don't reach.
func TestCLI_IdentifiesMultipleFiles(t *testing.T) {
	dir := t.TempDir()
	vcfPath := filepath.Join(dir, "in.vcf")
	if err := os.WriteFile(vcfPath, []byte("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	faPath := filepath.Join(dir, "in.fa")
	if err := os.WriteFile(faPath, []byte(">chr1\nACGTACGT\n"), 0o644); err != nil {
		t.Fatal(err)
	}

	bin := filepath.Join(dir, "htsfile")
	build := exec.Command("go", "build", "-o", bin, ".")
	build.Stderr = os.Stderr
	if err := build.Run(); err != nil {
		t.Fatalf("build: %v", err)
	}

	var out bytes.Buffer
	cmd := exec.Command(bin, vcfPath, faPath)
	cmd.Stdout = &out
	cmd.Stderr = os.Stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("run: %v", err)
	}
	got := out.String()
	for _, want := range []string{
		vcfPath + ":\tVCF version 4.2 variant calling text",
		faPath + ":\tFASTA sequence text",
	} {
		if !strings.Contains(got, want) {
			t.Errorf("missing %q in:\n%s", want, got)
		}
	}
}

// TestCLI_MissingFile pins the exit-1 + stderr message path when an
// argument refers to a file that doesn't exist.
func TestCLI_MissingFile(t *testing.T) {
	dir := t.TempDir()
	bin := filepath.Join(dir, "htsfile")
	if err := exec.Command("go", "build", "-o", bin, ".").Run(); err != nil {
		t.Fatalf("build: %v", err)
	}
	var stderr bytes.Buffer
	cmd := exec.Command(bin, filepath.Join(dir, "nope.bam"))
	cmd.Stderr = &stderr
	err := cmd.Run()
	if err == nil {
		t.Fatal("expected non-zero exit for missing file")
	}
	if !strings.Contains(stderr.String(), "htsfile:") {
		t.Errorf("stderr should mention 'htsfile:' diagnostic; got:\n%s", stderr.String())
	}
}
