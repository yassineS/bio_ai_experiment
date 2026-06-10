package main

import (
	"bytes"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

func buildBedmergeCLI(t *testing.T) string {
	t.Helper()
	bin := filepath.Join(t.TempDir(), "bedmerge")
	out, err := exec.Command("go", "build", "-o", bin, ".").CombinedOutput()
	if err != nil {
		t.Fatalf("go build bedmerge: %v\n%s", err, out)
	}
	return bin
}

func runBedmerge(t *testing.T, bin, stdin string, args ...string) (string, string, int) {
	t.Helper()
	cmd := exec.Command(bin, args...)
	if stdin != "" {
		cmd.Stdin = strings.NewReader(stdin)
	}
	var stdout, stderr bytes.Buffer
	cmd.Stdout = &stdout
	cmd.Stderr = &stderr
	err := cmd.Run()
	code := 0
	if err != nil {
		ee, ok := err.(*exec.ExitError)
		if !ok {
			t.Fatalf("running bedmerge %v: %v", args, err)
		}
		code = ee.ExitCode()
	}
	return stdout.String(), stderr.String(), code
}

// TestBedmergeHelpVersionStdin asserts -h/--help (exit 0), -v/--version (exit 0,
// banner), unknown flag (exit 2) and bare '-'/no-arg stdin.
func TestBedmergeHelpVersionStdin(t *testing.T) {
	bin := buildBedmergeCLI(t)

	for _, fl := range []string{"-h", "--help"} {
		out, errb, code := runBedmerge(t, bin, "", fl)
		if code != 0 || (out == "" && errb == "") {
			t.Errorf("bedmerge %s: exit=%d out=%q err=%q", fl, code, out, errb)
		}
	}
	for _, fl := range []string{"-v", "--version"} {
		out, _, code := runBedmerge(t, bin, "", fl)
		if code != 0 || !strings.Contains(out, "bedmerge version") {
			t.Errorf("bedmerge %s: exit=%d stdout=%q", fl, code, out)
		}
	}
	if _, _, code := runBedmerge(t, bin, "", "--no-such-flag"); code != 2 {
		t.Errorf("bedmerge --no-such-flag: exit=%d, want 2", code)
	}

	// Bare '-' input reads stdin.
	out, errb, code := runBedmerge(t, bin, "chr1\t10\t20\nchr1\t15\t25\n", "-")
	if code != 0 {
		t.Fatalf("bedmerge - (stdin): exit=%d stderr=%s", code, errb)
	}
	if !strings.Contains(out, "chr1\t10\t25") {
		t.Errorf("bedmerge merged stdin output=%q, want merged interval chr1 10-25", out)
	}
}
