package main

import (
	"bytes"
	"os/exec"
	"strings"
	"testing"
)

func runFastp(t *testing.T, bin string, args ...string) (string, string, int) {
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
			t.Fatalf("running fastp %v: %v", args, err)
		}
		code = ee.ExitCode()
	}
	return stdout.String(), stderr.String(), code
}

// TestFastpHelpVersion asserts fastp answers -h/--help (exit 0, output) and
// -v/--version (exit 0, version banner) and rejects unknown flags with exit 2.
func TestFastpHelpVersion(t *testing.T) {
	bin := buildFastpCLI(t)
	for _, fl := range []string{"-h", "--help"} {
		out, errb, code := runFastp(t, bin, fl)
		if code != 0 {
			t.Errorf("fastp %s: exit=%d, want 0", fl, code)
		}
		if out == "" && errb == "" {
			t.Errorf("fastp %s: no output", fl)
		}
	}
	for _, fl := range []string{"-v", "--version"} {
		out, _, code := runFastp(t, bin, fl)
		if code != 0 {
			t.Errorf("fastp %s: exit=%d, want 0", fl, code)
		}
		if !strings.Contains(out, "fastp version") {
			t.Errorf("fastp %s: stdout=%q want %q", fl, out, "fastp version")
		}
	}
	_, _, code := runFastp(t, bin, "--no-such-flag")
	if code != 2 {
		t.Errorf("fastp --no-such-flag: exit=%d, want 2", code)
	}
}
