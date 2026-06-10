package main

import (
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

// writeFixture writes content to a temp file and returns its path.
func writeFixture(t *testing.T, name, content string) string {
	t.Helper()
	p := filepath.Join(t.TempDir(), name)
	if err := os.WriteFile(p, []byte(content), 0o644); err != nil {
		t.Fatal(err)
	}
	return p
}

// TestTransformFlagsCLI exercises the new transform/misc knobs through the
// real command-line wiring.
func TestTransformFlagsCLI(t *testing.T) {
	bin := buildBinary(t)
	fasta := writeFixture(t, "in.fasta", ">read1 comment\nacgtACGT\n")

	cmd := exec.Command(bin, "filter", "--fasta", "-i", fasta,
		"--seq_case", "upper", "--dna_rna", "rna", "--seq_id", "S_", "--line_width", "4")
	out, err := cmd.CombinedOutput()
	if err != nil {
		t.Fatalf("filter: %v\n%s", err, out)
	}
	want := ">S_1 comment\nACGU\nACGU\n"
	if string(out) != want {
		t.Fatalf("CLI transform output mismatch\nwant: %q\ngot:  %q", want, string(out))
	}
}

// TestDefaultWrapCLI verifies the CLI wraps FASTA output at 60 columns by
// default (no --line_width), matching upstream's $linelen default, and that
// --line_width 0 disables wrapping.
func TestDefaultWrapCLI(t *testing.T) {
	bin := buildBinary(t)
	long := strings.Repeat("ACGT", 40) // 160 bases
	fasta := writeFixture(t, "in.fasta", ">a\n"+long+"\n")

	// Default: wrapped at 60.
	cmd := exec.Command(bin, "filter", "--fasta", "-i", fasta, "-l", "1")
	out, err := cmd.CombinedOutput()
	if err != nil {
		t.Fatalf("filter: %v\n%s", err, out)
	}
	lines := strings.Split(strings.TrimRight(string(out), "\n"), "\n")
	if len(lines) != 4 || len(lines[1]) != 60 || len(lines[3]) != 40 {
		t.Fatalf("default wrap not 60-col: %q", string(out))
	}

	// --line_width 0: single line.
	cmd = exec.Command(bin, "filter", "--fasta", "-i", fasta, "-l", "1", "--line_width", "0")
	out, err = cmd.CombinedOutput()
	if err != nil {
		t.Fatalf("filter: %v\n%s", err, out)
	}
	if string(out) != ">a\n"+long+"\n" {
		t.Fatalf("line_width 0 should not wrap: %q", string(out))
	}
}

// TestParamsFileCLI verifies that a -params file supplies flag values when
// the same flag is not given on the command line.
func TestParamsFileCLI(t *testing.T) {
	bin := buildBinary(t)
	fasta := writeFixture(t, "in.fasta", ">a comment\nacgt\n")
	params := writeFixture(t, "params.txt", "# transforms\n-seq_case upper\n-rm_header\n")

	cmd := exec.Command(bin, "filter", "--fasta", "-i", fasta, "--params", params)
	out, err := cmd.CombinedOutput()
	if err != nil {
		t.Fatalf("filter: %v\n%s", err, out)
	}
	// seq_case upper + rm_header (drops comment).
	want := ">a\nACGT\n"
	if string(out) != want {
		t.Fatalf("params-file output mismatch\nwant: %q\ngot:  %q", want, string(out))
	}
}

// TestParamsFileCLIPrecedence verifies that an explicit command-line flag wins
// over the params file.
func TestParamsFileCLIPrecedence(t *testing.T) {
	bin := buildBinary(t)
	fasta := writeFixture(t, "in.fasta", ">a\nacgt\n")
	params := writeFixture(t, "params.txt", "-seq_case upper\n")

	cmd := exec.Command(bin, "filter", "--fasta", "-i", fasta, "--params", params, "--seq_case", "lower")
	out, err := cmd.CombinedOutput()
	if err != nil {
		t.Fatalf("filter: %v\n%s", err, out)
	}
	if !strings.Contains(string(out), "\nacgt\n") {
		t.Fatalf("expected CLI seq_case lower to win: %q", string(out))
	}
}

// TestCustomParamsRepeatCLI verifies the repeat-rule rejection path through
// the CLI.
func TestCustomParamsRepeatCLI(t *testing.T) {
	bin := buildBinary(t)
	fasta := writeFixture(t, "in.fasta", ">keep\nACGTACGTAC\n>drop\nAAAAAAAAAA\n")

	cmd := exec.Command(bin, "filter", "--fasta", "-i", fasta, "--custom_params", "A 5")
	out, err := cmd.CombinedOutput()
	if err != nil {
		t.Fatalf("filter: %v\n%s", err, out)
	}
	got := string(out)
	if !strings.Contains(got, ">keep") || strings.Contains(got, ">drop") {
		t.Fatalf("custom_params repeat rule failed: %q", got)
	}
}

// TestTransformValidationErrors checks the domain/coupling validations exit
// non-zero with the expected messages.
func TestTransformValidationErrors(t *testing.T) {
	bin := buildBinary(t)
	fasta := writeFixture(t, "in.fasta", ">a\nACGT\n")

	cases := []struct {
		name   string
		args   []string
		errSub string
	}{
		{"bad_seq_case", []string{"--seq_case", "weird"}, "seq_case"},
		{"bad_dna_rna", []string{"--dna_rna", "xna"}, "dna_rna"},
		{"no_qual_header_fasta", []string{"--out_format", "1", "--no_qual_header"}, "no_qual_header"},
		// Upstream rejects --no_qual_header for any out_format != 3, including
		// the FASTQ+FASTA(+QUAL) modes 4 and 5 (prinseq-lite.pl:792-793).
		{"no_qual_header_fmt4", []string{"--out_format", "4", "--output", "x", "--no_qual_header"}, "no_qual_header"},
		{"no_qual_header_fmt5", []string{"--out_format", "5", "--output", "x", "--no_qual_header"}, "no_qual_header"},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			args := append([]string{"filter", "--fasta", "-i", fasta}, tc.args...)
			cmd := exec.Command(bin, args...)
			out, err := cmd.CombinedOutput()
			if err == nil {
				t.Fatalf("expected non-zero exit; output: %s", out)
			}
			if !strings.Contains(string(out), tc.errSub) {
				t.Fatalf("expected error mentioning %q, got: %s", tc.errSub, out)
			}
		})
	}
}
