package main

import (
	"bytes"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

// TestBedslopHelpVersion is a representative conformance check for an
// already-compliant bed* tool: -h/--help exit 0 with output, -v/--version exit
// 0 with the version banner, and an unknown flag exits 2.
func TestBedslopHelpVersion(t *testing.T) {
	bin := filepath.Join(t.TempDir(), "bedslop")
	if out, err := exec.Command("go", "build", "-o", bin, ".").CombinedOutput(); err != nil {
		t.Fatalf("go build bedslop: %v\n%s", err, out)
	}

	run := func(args ...string) (string, string, int) {
		cmd := exec.Command(bin, args...)
		var so, se bytes.Buffer
		cmd.Stdout = &so
		cmd.Stderr = &se
		err := cmd.Run()
		code := 0
		if err != nil {
			ee, ok := err.(*exec.ExitError)
			if !ok {
				t.Fatalf("running bedslop %v: %v", args, err)
			}
			code = ee.ExitCode()
		}
		return so.String(), se.String(), code
	}

	for _, fl := range []string{"-h", "--help"} {
		out, errb, code := run(fl)
		if code != 0 || (out == "" && errb == "") {
			t.Errorf("bedslop %s: exit=%d out=%q err=%q", fl, code, out, errb)
		}
	}
	for _, fl := range []string{"-v", "--version"} {
		out, _, code := run(fl)
		if code != 0 || !strings.Contains(out, "bedslop version") {
			t.Errorf("bedslop %s: exit=%d stdout=%q", fl, code, out)
		}
	}
	if _, _, code := run("--no-such-flag"); code != 2 {
		t.Errorf("bedslop --no-such-flag: exit=%d, want 2", code)
	}
}
