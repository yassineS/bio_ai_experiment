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

// TestUnitDepthFlagThresholdMapping is a binary-free unit test pinning bug
// #1's fix: `samtools depth -q` is the BASE-quality floor (per-base) and `-Q`
// is the MAPPING-quality floor (whole read), matching upstream bam2depth.c
// (case 'q' -> min_qual, case 'Q' -> min_mqual). It drives runDepth end to
// end (writing to a -o file, so no external binary is involved) on a fixture
// where the two knobs select genuinely different subsets.
func TestUnitDepthFlagThresholdMapping(t *testing.T) {
	dir := t.TempDir()
	// r1: MAPQ 60, all bases Phred 40 ('I').
	// r2: MAPQ 20, first two bases Phred 0 ('!'), rest 40.
	sam := "@HD\tVN:1.6\tSO:coordinate\n" +
		"@SQ\tSN:chr1\tLN:50\n" +
		"r1\t0\tchr1\t10\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n" +
		"r2\t0\tchr1\t10\t20\t5M\t*\t0\t0\tACGTA\t!!III\n"
	samPath := filepath.Join(dir, "in.sam")
	if err := os.WriteFile(samPath, []byte(sam), 0o644); err != nil {
		t.Fatal(err)
	}

	runDepthTo := func(t *testing.T, extra ...string) string {
		t.Helper()
		outPath := filepath.Join(dir, "out.txt")
		args := append([]string{"-o", outPath}, extra...)
		args = append(args, samPath)
		if rc := runDepth(args); rc != 0 {
			t.Fatalf("runDepth %v exit code %d", args, rc)
		}
		data, err := os.ReadFile(outPath)
		if err != nil {
			t.Fatalf("read depth output: %v", err)
		}
		return string(data)
	}

	// -q 30 (base quality): r2's first two bases ('!' = 0) are dropped, so
	// chr1:10 has depth 1 (only r1) and chr1:11 has depth 1.
	gotBase := runDepthTo(t, "-q", "30")
	wantBase := "chr1\t10\t1\nchr1\t11\t1\nchr1\t12\t2\nchr1\t13\t2\nchr1\t14\t2\n"
	if gotBase != wantBase {
		t.Errorf("depth -q 30 (base quality):\n got %q\nwant %q", gotBase, wantBase)
	}

	// -Q 30 (mapping quality): r2 (MAPQ 20) is dropped wholesale, leaving
	// only r1 across its 5 columns.
	gotMap := runDepthTo(t, "-Q", "30")
	wantMap := "chr1\t10\t1\nchr1\t11\t1\nchr1\t12\t1\nchr1\t13\t1\nchr1\t14\t1\n"
	if gotMap != wantMap {
		t.Errorf("depth -Q 30 (mapping quality):\n got %q\nwant %q", gotMap, wantMap)
	}

	// The long spellings --min-BQ / --min-MQ alias the same knobs.
	if got := runDepthTo(t, "--min-BQ", "30"); got != wantBase {
		t.Errorf("depth --min-BQ 30 should equal -q 30:\n got %q\nwant %q", got, wantBase)
	}
	if got := runDepthTo(t, "--min-MQ", "30"); got != wantMap {
		t.Errorf("depth --min-MQ 30 should equal -Q 30:\n got %q\nwant %q", got, wantMap)
	}
}
