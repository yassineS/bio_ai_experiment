package main

import (
	"bytes"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

func buildSickleCLIConf(t *testing.T) string {
	t.Helper()
	bin := filepath.Join(t.TempDir(), "sickle")
	out, err := exec.Command("go", "build", "-o", bin, ".").CombinedOutput()
	if err != nil {
		t.Fatalf("go build sickle: %v\n%s", err, out)
	}
	return bin
}

func runSickle(t *testing.T, bin string, args ...string) (string, string, int) {
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
			t.Fatalf("running sickle %v: %v", args, err)
		}
		code = ee.ExitCode()
	}
	return stdout.String(), stderr.String(), code
}

// TestSickleSubcommandHelpVersion asserts the se/pe/batch subcommands answer
// -h/--help (exit 0, output), -v/--version (exit 0, version banner), and reject
// an unknown flag with exit 2.
func TestSickleSubcommandHelpVersion(t *testing.T) {
	bin := buildSickleCLIConf(t)
	for _, sub := range []string{"se", "pe", "batch"} {
		t.Run(sub, func(t *testing.T) {
			for _, fl := range []string{"-h", "--help"} {
				out, errb, code := runSickle(t, bin, sub, fl)
				if code != 0 {
					t.Errorf("%s %s: exit=%d, want 0", sub, fl, code)
				}
				if out == "" && errb == "" {
					t.Errorf("%s %s: no output", sub, fl)
				}
			}
			for _, fl := range []string{"-v", "--version"} {
				out, _, code := runSickle(t, bin, sub, fl)
				if code != 0 {
					t.Errorf("%s %s: exit=%d, want 0", sub, fl, code)
				}
				if !strings.Contains(out, "sickle version") {
					t.Errorf("%s %s: stdout=%q want %q", sub, fl, out, "sickle version")
				}
			}
			_, _, code := runSickle(t, bin, sub, "--nope-not-a-flag")
			if code != 2 {
				t.Errorf("%s --nope-not-a-flag: exit=%d, want 2", sub, code)
			}
		})
	}
}
