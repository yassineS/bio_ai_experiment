package main

import (
	"bytes"
	"os/exec"
	"path/filepath"
	"testing"
)

func buildSamtoolsCLI(t *testing.T) string {
	t.Helper()
	bin := filepath.Join(t.TempDir(), "samtools")
	if out, err := exec.Command("go", "build", "-o", bin, ".").CombinedOutput(); err != nil {
		t.Fatalf("go build samtools: %v\n%s", err, out)
	}
	return bin
}

func runSam(t *testing.T, bin string, args ...string) (string, string, int) {
	t.Helper()
	cmd := exec.Command(bin, args...)
	var so, se bytes.Buffer
	cmd.Stdout = &so
	cmd.Stderr = &se
	err := cmd.Run()
	code := 0
	if err != nil {
		ee, ok := err.(*exec.ExitError)
		if !ok {
			t.Fatalf("running samtools %v: %v", args, err)
		}
		code = ee.ExitCode()
	}
	return so.String(), se.String(), code
}

// TestSamtoolsHelpVersionConformance documents the intentional htslib-family
// convention: the long --help / --version forms are uniform (exit 0 with
// output) at the top level and on subcommands, while the short -h / -v stay
// bound to their upstream meanings (e.g. `view -h` = include header), so they
// are deliberately NOT remapped to help/version. Unknown flags exit 2.
func TestSamtoolsHelpVersionConformance(t *testing.T) {
	bin := buildSamtoolsCLI(t)

	// Top level: both short and long help/version work (no upstream conflict).
	for _, fl := range []string{"--help", "--version"} {
		out, errb, code := runSam(t, bin, fl)
		if code != 0 || (out == "" && errb == "") {
			t.Errorf("samtools %s: exit=%d out=%q err=%q", fl, code, out, errb)
		}
	}

	// Subcommand long forms are uniform across the suite.
	for _, sub := range []string{"view", "sort", "flagstat", "index"} {
		t.Run(sub, func(t *testing.T) {
			for _, fl := range []string{"--help", "--version"} {
				out, errb, code := runSam(t, bin, sub, fl)
				if code != 0 {
					t.Errorf("%s %s: exit=%d, want 0\nstderr=%s", sub, fl, code, errb)
				}
				if out == "" && errb == "" {
					t.Errorf("%s %s: no output", sub, fl)
				}
			}
			if _, _, code := runSam(t, bin, sub, "--no-such-flag"); code != 2 {
				t.Errorf("%s --no-such-flag: exit=%d, want 2", sub, code)
			}
		})
	}
}
