package main

import (
	"os"
	"os/exec"
	"path/filepath"
	"testing"
)

// TestFilterTrimRangeFlagsWiring exercises the new quality-trim window/step/
// rule, trim_to_len, and range_len/range_gc flags through the real CLI so a
// regression in flag registration is caught. It checks behaviour, not exact
// upstream bytes (the package-level parity tests own byte-for-byte parity).
func TestFilterTrimRangeFlagsWiring(t *testing.T) {
	bin := buildBinary(t)
	dir := t.TempDir()

	in := filepath.Join(dir, "in.fastq")
	if err := os.WriteFile(in, []byte(
		"@r1\nACGTACGTACGT\n+\n!!!!IIIIIIII\n"+
			"@r2\nGGGGCCCCAAAA\n+\nIIIIIIIIIIII\n"), 0o644); err != nil {
		t.Fatalf("write input: %v", err)
	}

	cases := []struct {
		name string
		args []string
	}{
		{"trim_qual_window", []string{"--trim-qual-left", "20", "--trim_qual_window", "4", "--trim_qual_step", "1", "--trim_qual_type", "mean"}},
		{"trim_qual_rule", []string{"--trim-qual-right", "20", "--trim_qual_rule", "lt"}},
		{"trim_to_len", []string{"--trim_to_len", "6"}},
		{"range_len", []string{"--range_len", "12-12"}},
		{"range_gc", []string{"--range_gc", "0-100"}},
		// Hyphenated aliases must also resolve.
		{"hyphen_aliases", []string{"--trim-qual-left", "20", "--trim-qual-window", "2", "--trim-to-len", "8", "--range-gc", "0-100"}},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			out := filepath.Join(dir, tc.name+".out.fastq")
			args := append([]string{"filter", "--fastq", "-i", in, "-o", out}, tc.args...)
			cmd := exec.Command(bin, args...)
			cmd.Stderr = os.Stderr
			if err := cmd.Run(); err != nil {
				t.Fatalf("filter %s: %v", tc.name, err)
			}
			if _, err := os.Stat(out); err != nil {
				t.Fatalf("expected output file for %s: %v", tc.name, err)
			}
		})
	}
}

// TestFilterTrimQualRuleValidation confirms the CLI rejects an invalid
// --trim_qual_rule value, matching upstream's "invalid value" error.
func TestFilterTrimQualRuleValidation(t *testing.T) {
	bin := buildBinary(t)
	dir := t.TempDir()
	in := filepath.Join(dir, "in.fastq")
	if err := os.WriteFile(in, []byte("@r1\nACGT\n+\nIIII\n"), 0o644); err != nil {
		t.Fatalf("write input: %v", err)
	}
	cmd := exec.Command(bin, "filter", "--fastq", "-i", in, "--trim-qual-left", "20", "--trim_qual_rule", "bogus")
	cmd.Stdout = os.Stdout
	if err := cmd.Run(); err == nil {
		t.Fatal("expected non-zero exit for invalid --trim_qual_rule")
	}
}
