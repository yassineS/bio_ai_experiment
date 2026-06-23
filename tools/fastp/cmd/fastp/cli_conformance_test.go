package main

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
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

// TestFastpHelpVersion asserts fastp answers -?/--help (exit 0, output) and
// -v/--version (exit 0, version banner) and rejects unknown flags with exit 2.
// Note: -h is upstream's HTML report flag (not help), so help is reached via
// -? / --help to keep the CLI drop-in compatible with stock fastp.
func TestFastpHelpVersion(t *testing.T) {
	bin := buildFastpCLI(t)
	for _, fl := range []string{"-?", "--help"} {
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

// TestFastpUpstreamCLIDropIn asserts that a stock upstream fastp command line
// — -i in1, -I in2, -o out1, -O out2, -j json, -h html — is accepted by our
// CLI (drop-in compatibility). It exercises the upstream-exact short-flag
// remap (-i=read1, -I=read2, -o=out1, -O=out2, -j=json, -h=html) and verifies
// both paired outputs and the JSON report are written.
func TestFastpUpstreamCLIDropIn(t *testing.T) {
	bin := buildFastpCLI(t)
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fq")
	r2 := filepath.Join(dir, "r2.fq")
	// A short overlapping pair (R2 is the reverse complement of R1) so the
	// pipeline runs the PE overlap path without error.
	if err := os.WriteFile(r1, []byte("@p/1\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	if err := os.WriteFile(r2, []byte("@p/2\nACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIII\n"), 0o644); err != nil {
		t.Fatal(err)
	}
	o1 := filepath.Join(dir, "o1.fq")
	o2 := filepath.Join(dir, "o2.fq")
	js := filepath.Join(dir, "report.json")
	ht := filepath.Join(dir, "report.html")
	_, errb, code := runFastp(t, bin, "-i", r1, "-I", r2, "-o", o1, "-O", o2, "-j", js, "-h", ht)
	if code != 0 {
		t.Fatalf("upstream-style fastp CLI: exit=%d, want 0\nstderr:\n%s", code, errb)
	}
	for _, f := range []string{o1, o2, js, ht} {
		if fi, err := os.Stat(f); err != nil || fi.Size() == 0 {
			t.Errorf("expected non-empty output %s (err=%v)", f, err)
		}
	}
}
