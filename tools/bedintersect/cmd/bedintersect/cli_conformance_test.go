package main

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

func buildBedintersectCLI(t *testing.T) string {
	t.Helper()
	bin := filepath.Join(t.TempDir(), "bedintersect")
	out, err := exec.Command("go", "build", "-o", bin, ".").CombinedOutput()
	if err != nil {
		t.Fatalf("go build bedintersect: %v\n%s", err, out)
	}
	return bin
}

func runBedint(t *testing.T, bin string, args ...string) (string, string, int) {
	t.Helper()
	cmd := exec.Command(bin, args...)
	var stdout, stderr bytes.Buffer
	cmd.Stdout = &stdout
	cmd.Stderr = &stderr
	err := cmd.Run()
	code := 0
	if err != nil {
		ee, ok := err.(*exec.ExitError)
		if !ok {
			t.Fatalf("running bedintersect %v: %v", args, err)
		}
		code = ee.ExitCode()
	}
	return stdout.String(), stderr.String(), code
}

// TestBedintersectHelpVersion asserts -h/--help exit 0 with output, --version
// exits 0 with the version banner, and that -v stays bound to --invert (upstream
// bedtools parity) rather than printing the version.
func TestBedintersectHelpVersion(t *testing.T) {
	bin := buildBedintersectCLI(t)

	for _, fl := range []string{"-h", "--help"} {
		out, errb, code := runBedint(t, bin, fl)
		if code != 0 || (out == "" && errb == "") {
			t.Errorf("bedintersect %s: exit=%d out=%q err=%q", fl, code, out, errb)
		}
	}

	out, _, code := runBedint(t, bin, "--version")
	if code != 0 || !strings.Contains(out, "bedintersect version") {
		t.Errorf("bedintersect --version: exit=%d stdout=%q", code, out)
	}

	// -v must NOT be version: it is --invert and requires -a/-b, so a bare -v
	// errors on the missing inputs (exit 1) and never prints a version banner.
	out, errb, code := runBedint(t, bin, "-v")
	if strings.Contains(out, "bedintersect version") {
		t.Errorf("-v printed a version banner; it must remain --invert for bedtools parity")
	}
	if code != 1 {
		t.Errorf("bedintersect -v (no -a/-b): exit=%d, want 1 (missing inputs); stderr=%s", code, errb)
	}
}

// TestBedintersectStdinDash asserts file A can be read from stdin via '-'.
func TestBedintersectStdinDash(t *testing.T) {
	bin := buildBedintersectCLI(t)
	dir := t.TempDir()
	fileB := filepath.Join(dir, "b.bed")
	if err := os.WriteFile(fileB, []byte("chr1\t12\t18\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	cmd := exec.Command(bin, "-a", "-", "-b", fileB)
	cmd.Stdin = strings.NewReader("chr1\t10\t20\n")
	var stdout, stderr bytes.Buffer
	cmd.Stdout = &stdout
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("bedintersect -a - : %v\nstderr=%s", err, stderr.String())
	}
	if !strings.Contains(stdout.String(), "chr1\t12\t18") {
		t.Errorf("bedintersect -a - output=%q, want intersection chr1 12-18", stdout.String())
	}
}
