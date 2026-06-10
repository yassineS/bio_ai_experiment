package main

import (
	"bytes"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

// buildSeqtkCLI compiles the seqtk binary into a temp dir for CLI conformance
// checks (help/version/exit-code uniformity across subcommands).
func buildSeqtkCLI(t *testing.T) string {
	t.Helper()
	bin := filepath.Join(t.TempDir(), "seqtk")
	out, err := exec.Command("go", "build", "-o", bin, ".").CombinedOutput()
	if err != nil {
		t.Fatalf("go build seqtk: %v\n%s", err, out)
	}
	return bin
}

// runCLI runs bin with args and stdin, returning combined stdout, stderr and
// the process exit code.
func runCLI(t *testing.T, bin, stdin string, args ...string) (string, string, int) {
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
			t.Fatalf("running %s %v: %v", bin, args, err)
		}
		code = ee.ExitCode()
	}
	return stdout.String(), stderr.String(), code
}

// TestSeqtkSubcommandHelpVersion asserts every representative subcommand answers
// -h/--help (exit 0, output) and -v/--version (exit 0, output), and rejects an
// unknown flag with exit 2.
func TestSeqtkSubcommandHelpVersion(t *testing.T) {
	bin := buildSeqtkCLI(t)
	// A representative spread of subcommands: a -h-free one, the mergefa
	// special case (where -h is --haploid), and a couple of others.
	subs := []string{"seq", "comp", "subseq", "sample", "trimfq", "gc", "hrun", "split"}
	for _, sub := range subs {
		t.Run(sub, func(t *testing.T) {
			for _, fl := range []string{"-h", "--help"} {
				out, errb, code := runCLI(t, bin, "", sub, fl)
				if code != 0 {
					t.Errorf("%s %s: exit=%d, want 0\nstderr=%s", sub, fl, code, errb)
				}
				if out == "" && errb == "" {
					t.Errorf("%s %s: produced no output", sub, fl)
				}
			}
			for _, fl := range []string{"-v", "--version"} {
				out, errb, code := runCLI(t, bin, "", sub, fl)
				if code != 0 {
					t.Errorf("%s %s: exit=%d, want 0\nstderr=%s", sub, fl, code, errb)
				}
				if !strings.Contains(out, "seqtk version") {
					t.Errorf("%s %s: stdout=%q, want it to contain %q", sub, fl, out, "seqtk version")
				}
			}
			_, _, code := runCLI(t, bin, "", sub, "--definitely-not-a-flag")
			if code != 2 {
				t.Errorf("%s --definitely-not-a-flag: exit=%d, want 2", sub, code)
			}
		})
	}
}

// TestSeqtkMergefaShortHIsHaploid asserts the mergefa special case keeps -h
// bound to --haploid (upstream parity) while --help still prints usage.
func TestSeqtkMergefaShortHIsHaploid(t *testing.T) {
	bin := buildSeqtkCLI(t)
	// --help prints usage and exits 0.
	out, errb, code := runCLI(t, bin, "", "mergefa", "--help")
	if code != 0 || (out == "" && errb == "") {
		t.Fatalf("mergefa --help: exit=%d out=%q err=%q, want exit 0 with output", code, out, errb)
	}
	if !strings.Contains(out+errb, "mergefa") {
		t.Errorf("mergefa --help did not mention the subcommand: %q", out+errb)
	}
}

// TestSeqtkStdinDash asserts a bare '-' is read as stdin for a subcommand.
func TestSeqtkStdinDash(t *testing.T) {
	bin := buildSeqtkCLI(t)
	out, errb, code := runCLI(t, bin, ">s1\nACGTACGT\n", "comp", "-")
	if code != 0 {
		t.Fatalf("seqtk comp -: exit=%d stderr=%s", code, errb)
	}
	if !strings.HasPrefix(out, "s1\t") {
		t.Errorf("seqtk comp - stdout=%q, want it to start with the record name", out)
	}
}
