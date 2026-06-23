package main

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

// flatTestInput is a small FASTQ fixture with one read that survives a
// modest min_len filter and one short read that should be rejected.
const flatTestInput = "@r1\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n" +
	"@r2_short\nACGT\n+\nIIII\n" +
	"@r3\nGGGGCCCCAAAATTTTGGGG\n+\nIIIIIIIIIIIIIIIIIIII\n"

// TestFlatCLIDropIn asserts that the upstream flat-flag CLI is accepted as a
// drop-in: "prinseq -fastq IN -out_good PREFIX -out_bad null -min_len N ...".
// It checks the upstream output-naming convention ("<prefix>.fastq") and that
// the flat form yields exactly the same bytes as the equivalent subcommand
// invocation (the package-level parity tests own byte-for-byte parity vs the
// real prinseq-lite.pl; this test owns CLI-shim fidelity).
func TestFlatCLIDropIn(t *testing.T) {
	bin := buildBinary(t)
	dir := t.TempDir()
	in := filepath.Join(dir, "in.fastq")
	if err := os.WriteFile(in, []byte(flatTestInput), 0o644); err != nil {
		t.Fatalf("write input: %v", err)
	}

	// Flat upstream CLI: no subcommand, single-dash long options, -out_good
	// PREFIX writes to "<PREFIX>.fastq".
	goodPrefix := filepath.Join(dir, "flat_good")
	cmd := exec.Command(bin, "-fastq", in,
		"-out_good", goodPrefix, "-out_bad", "null", "-min_len", "10")
	var stderr bytes.Buffer
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("flat CLI exited non-zero: %v\nstderr=%s", err, stderr.String())
	}

	flatOut := goodPrefix + ".fastq"
	flatBytes, err := os.ReadFile(flatOut)
	if err != nil {
		t.Fatalf("flat CLI did not write %q: %v", flatOut, err)
	}
	// r1 and r3 pass min_len 10; r2_short (4 bp) is rejected.
	if got := strings.Count(string(flatBytes), "@r"); got != 2 {
		t.Errorf("flat good output: got %d reads, want 2\n%s", got, flatBytes)
	}
	if strings.Contains(string(flatBytes), "r2_short") {
		t.Errorf("flat good output unexpectedly contains the short read:\n%s", flatBytes)
	}

	// Equivalent subcommand form should produce identical bytes.
	subOut := filepath.Join(dir, "sub_good.fastq")
	subCmd := exec.Command(bin, "filter", "--fastq", "-i", in, "-o", subOut, "--min-length", "10")
	subCmd.Stderr = &stderr
	if err := subCmd.Run(); err != nil {
		t.Fatalf("subcommand exited non-zero: %v\nstderr=%s", err, stderr.String())
	}
	subBytes, err := os.ReadFile(subOut)
	if err != nil {
		t.Fatalf("read subcommand output: %v", err)
	}
	if !bytes.Equal(flatBytes, subBytes) {
		t.Errorf("flat CLI output differs from subcommand output\nflat:\n%s\nsub:\n%s", flatBytes, subBytes)
	}
}

// TestFlatCLIStdout asserts "-out_good stdout" streams good reads to STDOUT and
// "-out_bad null" suppresses the bad stream, without creating stray files.
func TestFlatCLIStdout(t *testing.T) {
	bin := buildBinary(t)
	dir := t.TempDir()
	in := filepath.Join(dir, "in.fastq")
	if err := os.WriteFile(in, []byte(flatTestInput), 0o644); err != nil {
		t.Fatalf("write input: %v", err)
	}

	cmd := exec.Command(bin, "-fastq", in, "-out_good", "stdout", "-out_bad", "null", "-min_len", "10")
	cmd.Dir = dir
	var stdout, stderr bytes.Buffer
	cmd.Stdout = &stdout
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("flat -out_good stdout exited non-zero: %v\nstderr=%s", err, stderr.String())
	}
	if got := strings.Count(stdout.String(), "@r"); got != 2 {
		t.Errorf("stdout: got %d reads, want 2\n%s", got, stdout.String())
	}
	// No stray "null"-named or ".fastq" files should be created in the run dir.
	entries, _ := os.ReadDir(dir)
	for _, e := range entries {
		if e.Name() == "null" || e.Name() == "stdout" || strings.HasSuffix(e.Name(), ".fastq") && e.Name() != "in.fastq" {
			t.Errorf("unexpected stray output file created: %q", e.Name())
		}
	}
}

// TestFlatCLITrimQual asserts the upstream -trim_qual_right flag is wired and
// trims low-quality 3' bases (the package owns the exact trimming algorithm;
// here we only confirm the flat flag reaches it).
func TestFlatCLITrimQual(t *testing.T) {
	bin := buildBinary(t)
	dir := t.TempDir()
	in := filepath.Join(dir, "in.fastq")
	// A read whose 3' half is low quality (below Phred 20) so trim_qual_right
	// must shorten it.
	if err := os.WriteFile(in, []byte("@r1\nACGTACGTACGTACGT\n+\nIIIIIIII!!!!!!!!\n"), 0o644); err != nil {
		t.Fatalf("write input: %v", err)
	}

	goodPrefix := filepath.Join(dir, "good")
	cmd := exec.Command(bin, "-fastq", in, "-out_good", goodPrefix,
		"-out_bad", "null", "-min_len", "1", "-trim_qual_right", "20")
	var stderr bytes.Buffer
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("flat -trim_qual_right exited non-zero: %v\nstderr=%s", err, stderr.String())
	}
	out, err := os.ReadFile(goodPrefix + ".fastq")
	if err != nil {
		t.Fatalf("read output: %v", err)
	}
	lines := strings.Split(strings.TrimRight(string(out), "\n"), "\n")
	if len(lines) != 4 {
		t.Fatalf("expected one FASTQ record (4 lines), got %d:\n%s", len(lines), out)
	}
	// The 8 low-quality 3' bases must have been trimmed off.
	if len(lines[1]) != 8 {
		t.Errorf("trim_qual_right: trimmed sequence = %q (len %d), want length 8", lines[1], len(lines[1]))
	}
}
