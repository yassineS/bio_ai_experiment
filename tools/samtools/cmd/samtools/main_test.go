package main

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// TestParseOutputFmtOptions checks the `--output-fmt-option` KEY=VALUE
// parser: a valid qbin value is extracted, last-wins on repeats, and a
// malformed entry or unknown key is an error.
func TestParseOutputFmtOptions(t *testing.T) {
	if got, err := parseOutputFmtOptions(nil); err != nil || got != "" {
		t.Errorf("parseOutputFmtOptions(nil) = (%q, %v), want (\"\", nil)", got, err)
	}
	if got, err := parseOutputFmtOptions(multiString{"qbin=8"}); err != nil || got != "8" {
		t.Errorf("parseOutputFmtOptions(qbin=8) = (%q, %v), want (\"8\", nil)", got, err)
	}
	// Last value wins for a repeated key.
	if got, err := parseOutputFmtOptions(multiString{"qbin=8", "qbin=4"}); err != nil || got != "4" {
		t.Errorf("repeated qbin = (%q, %v), want (\"4\", nil)", got, err)
	}
	if _, err := parseOutputFmtOptions(multiString{"qbin"}); err == nil {
		t.Error("parseOutputFmtOptions(\"qbin\") should fail (no '=')")
	}
	if _, err := parseOutputFmtOptions(multiString{"bogus=1"}); err == nil {
		t.Error("parseOutputFmtOptions with unknown key should fail")
	}
}

// TestRunMpileupBCFWiring drives the CLI `mpileup -g` path end-to-end:
// it parses real args, runs the genotype-likelihood caller, and asserts a
// BGZF-wrapped BCF file is produced on disk. `-u` is exercised for the
// uncompressed container. This pins the cmd-level wiring that maps the
// samtools mpileup flags onto the shared bcftools mpileup engine.
func TestRunMpileupBCFWiring(t *testing.T) {
	dir := t.TempDir()
	fa := filepath.Join(dir, "ref.fa")
	if err := os.WriteFile(fa, []byte(">chr1\n"+strings.Repeat("A", 20)+"\n"), 0o644); err != nil {
		t.Fatalf("write fasta: %v", err)
	}
	if err := os.WriteFile(fa+".fai", []byte("chr1\t20\t6\t20\t21\n"), 0o644); err != nil {
		t.Fatalf("write fai: %v", err)
	}
	samPath := filepath.Join(dir, "in.sam")
	sam := strings.Join([]string{
		"@HD\tVN:1.6\tSO:coordinate",
		"@SQ\tSN:chr1\tLN:20",
		"@RG\tID:rg1\tSM:s1",
		"r1\t0\tchr1\t1\t60\t8M\t*\t0\t0\tAAAACAAA\t????????\tRG:Z:rg1",
		"",
	}, "\n")
	if err := os.WriteFile(samPath, []byte(sam), 0o644); err != nil {
		t.Fatalf("write sam: %v", err)
	}

	cases := []struct {
		name      string
		flag      string
		wantMagic []byte
	}{
		{"bcf_g", "-g", []byte{0x1f, 0x8b}},     // BGZF/gzip magic
		{"ubcf_u", "-u", []byte{'B', 'C', 'F'}}, // uncompressed BCF magic
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			outPath := filepath.Join(dir, tc.name+".bcf")
			rc := runMpileup([]string{tc.flag, "-f", fa, "-o", outPath, samPath})
			if rc != 0 {
				t.Fatalf("runMpileup %s exit=%d, want 0", tc.flag, rc)
			}
			data, err := os.ReadFile(outPath)
			if err != nil {
				t.Fatalf("read output: %v", err)
			}
			if len(data) < len(tc.wantMagic) || !bytesHasPrefix(data, tc.wantMagic) {
				t.Fatalf("%s output missing magic %q, got % x", tc.flag, tc.wantMagic, data[:min(4, len(data))])
			}
		})
	}
}

func bytesHasPrefix(b, prefix []byte) bool {
	if len(b) < len(prefix) {
		return false
	}
	for i := range prefix {
		if b[i] != prefix[i] {
			return false
		}
	}
	return true
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}
