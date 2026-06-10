package main

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

func buildPrinseqCLIConf(t *testing.T) string {
	t.Helper()
	bin := filepath.Join(t.TempDir(), "prinseq")
	out, err := exec.Command("go", "build", "-o", bin, ".").CombinedOutput()
	if err != nil {
		t.Fatalf("go build prinseq: %v\n%s", err, out)
	}
	return bin
}

func runPrinseq(t *testing.T, bin, stdin string, args ...string) (string, string, int) {
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
			t.Fatalf("running prinseq %v: %v", args, err)
		}
		code = ee.ExitCode()
	}
	return stdout.String(), stderr.String(), code
}

// TestPrinseqSubcommandHelpVersion asserts the subcommands answer -h/--help
// (exit 0) and -v/--version (exit 0, version banner) and reject unknown flags
// with exit 2.
func TestPrinseqSubcommandHelpVersion(t *testing.T) {
	bin := buildPrinseqCLIConf(t)
	for _, sub := range []string{"stats", "filter", "graph", "report", "benchmark", "api", "batch"} {
		t.Run(sub, func(t *testing.T) {
			for _, fl := range []string{"-h", "--help"} {
				out, errb, code := runPrinseq(t, bin, "", sub, fl)
				if code != 0 {
					t.Errorf("%s %s: exit=%d, want 0", sub, fl, code)
				}
				if out == "" && errb == "" {
					t.Errorf("%s %s: no output", sub, fl)
				}
			}
			for _, fl := range []string{"-v", "--version"} {
				out, _, code := runPrinseq(t, bin, "", sub, fl)
				if code != 0 {
					t.Errorf("%s %s: exit=%d, want 0", sub, fl, code)
				}
				if !strings.Contains(out, "prinseq version") {
					t.Errorf("%s %s: stdout=%q want %q", sub, fl, out, "prinseq version")
				}
			}
			_, _, code := runPrinseq(t, bin, "", sub, "--no-such-flag")
			if code != 2 {
				t.Errorf("%s --no-such-flag: exit=%d, want 2", sub, code)
			}
		})
	}
}

// TestPrinseqFilterStdinAndDashOut asserts filter reads '-' as stdin and writes
// '-' as stdout rather than creating a file literally named "-".
func TestPrinseqFilterStdinAndDashOut(t *testing.T) {
	bin := buildPrinseqCLIConf(t)
	dir := t.TempDir()
	// Run from a temp dir so any stray "-" file would land there.
	cmd := exec.Command(bin, "filter", "--fasta", "-i", "-", "-o", "-", "-l", "1")
	cmd.Dir = dir
	cmd.Stdin = strings.NewReader(">s1\nACGTACGT\n")
	var stdout, stderr bytes.Buffer
	cmd.Stdout = &stdout
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("prinseq filter -i - -o -: %v\nstderr=%s", err, stderr.String())
	}
	if !strings.Contains(stdout.String(), "ACGTACGT") {
		t.Errorf("expected sequence on stdout, got %q", stdout.String())
	}
	if _, err := os.Stat(filepath.Join(dir, "-")); err == nil {
		t.Errorf("filter -o - created a file literally named '-' instead of writing to stdout")
	}
}
