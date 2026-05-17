// Wave-16 CLI-binary tests for the upstream-canonical flag aliases added
// to close gaps versus reference_code/vcftools/src/cpp/parameters.cpp:
//
//   - `-c` — upstream's short alias for `--stdout` (parameters.cpp:194).
//   - `--recode-INFO TAG` — upstream's canonical name for the repeatable
//     recode-INFO-column selector (parameters.cpp:319). The port's
//     `--keep-INFO` already implements this semantic; we add the
//     canonical spelling as a synonym, mirroring the existing
//     `--keep-INFO-all` ↔ `--recode-INFO-all` pattern (main.go).
//
// Both aliases are wired in main.go and have no `pkg/vcftools` Params
// equivalent of their own (they OR into the same fields as the long
// spellings), so the tests live here next to the CLI to exercise the
// actual argv → flag.Parse path rather than the package boundary.
package main

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

// buildVcftools compiles the CLI into a temp file. The pattern mirrors
// tools/prinseq/cmd/prinseq/graphdata_cli_test.go.
func buildVcftools(t *testing.T) string {
	t.Helper()
	dir := t.TempDir()
	bin := filepath.Join(dir, "vcftools")
	cmd := exec.Command("go", "build", "-o", bin, ".")
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("go build: %v", err)
	}
	return bin
}

const aliasVCF = "##fileformat=VCFv4.2\n" +
	"##INFO=<ID=AF,Number=1,Type=Float,Description=\"\">\n" +
	"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"\">\n" +
	"##INFO=<ID=AA,Number=1,Type=String,Description=\"\">\n" +
	"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">\n" +
	"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n" +
	"1\t100\t.\tA\tG\t.\tPASS\tAF=0.5;DP=10;AA=T\tGT\t0/0\n"

// runWithStdin runs the binary with the given argv and stdin payload.
func runWithStdin(t *testing.T, bin string, stdin string, args ...string) (string, string) {
	t.Helper()
	cmd := exec.Command(bin, args...)
	cmd.Stdin = strings.NewReader(stdin)
	var stdout, stderr bytes.Buffer
	cmd.Stdout = &stdout
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("run %v: %v; stderr=%s", args, err, stderr.String())
	}
	return stdout.String(), stderr.String()
}

// dataRowINFO returns the INFO field (column 8) of the first data row
// found in a VCF-format body, or "" if none is found.
func dataRowINFO(body string) string {
	for _, line := range strings.Split(body, "\n") {
		if strings.HasPrefix(line, "1\t") {
			fields := strings.Split(line, "\t")
			if len(fields) >= 8 {
				return fields[7]
			}
		}
	}
	return ""
}

// TestCLI_ShortStdoutFlag — `-c` is upstream's short alias for `--stdout`
// (parameters.cpp:194). It must route output through stdout exactly like
// `--stdout` does.
func TestCLI_ShortStdoutFlag(t *testing.T) {
	bin := buildVcftools(t)
	dir := t.TempDir()
	stdout, _ := runWithStdin(t, bin, aliasVCF,
		"--stdin", "--recode", "-c", "--out", filepath.Join(dir, "out"))
	if !strings.Contains(stdout, "##fileformat=VCFv4.2") {
		t.Fatalf("-c did not produce VCF on stdout; got:\n%s", stdout)
	}
	if dataRowINFO(stdout) == "" {
		t.Fatalf("-c stdout has no data row; got:\n%s", stdout)
	}
	// Confirm no .recode.vcf file was written to disk.
	if _, err := os.Stat(filepath.Join(dir, "out.recode.vcf")); err == nil {
		t.Errorf("with -c, expected no on-disk .recode.vcf, but the file exists")
	}
}

// TestCLI_RecodeINFOAlias — `--recode-INFO TAG` must restrict the INFO
// column in `.recode.vcf` to the listed tags, exactly like the port's
// `--keep-INFO` (parameters.cpp:319 + port's existing
// pkg/vcftools/info_filters.go::filterRecodeInfo).
func TestCLI_RecodeINFOAlias(t *testing.T) {
	bin := buildVcftools(t)
	dir := t.TempDir()
	prefix := filepath.Join(dir, "out")
	_, _ = runWithStdin(t, bin, aliasVCF,
		"--stdin", "--recode", "--recode-INFO", "AF", "--out", prefix)

	data, err := os.ReadFile(prefix + ".recode.vcf")
	if err != nil {
		t.Fatalf("read recode output: %v", err)
	}
	got := dataRowINFO(string(data))
	if got != "AF=0.5" {
		t.Errorf("--recode-INFO AF: INFO = %q, want %q", got, "AF=0.5")
	}
}

// TestCLI_RecodeINFOAlias_Repeatable — upstream's `--recode-INFO` is
// repeatable (`recode_INFO_to_keep.insert(...)` on each occurrence,
// parameters.cpp:319). Two consecutive `--recode-INFO` arguments must
// union into the kept set.
func TestCLI_RecodeINFOAlias_Repeatable(t *testing.T) {
	bin := buildVcftools(t)
	dir := t.TempDir()
	prefix := filepath.Join(dir, "out")
	_, _ = runWithStdin(t, bin, aliasVCF,
		"--stdin", "--recode",
		"--recode-INFO", "AF",
		"--recode-INFO", "DP",
		"--out", prefix)

	data, err := os.ReadFile(prefix + ".recode.vcf")
	if err != nil {
		t.Fatalf("read recode output: %v", err)
	}
	got := dataRowINFO(string(data))
	// Recoded INFO ordering: vcftools.go does not propagate InfoOrder
	// onto the recoded variant; surviving keys land in formatInfo's
	// alphabetical-leftover branch (pkg/bioformats/vcf/vcf.go:394-409).
	// AF precedes DP alphabetically.
	if got != "AF=0.5;DP=10" {
		t.Errorf("--recode-INFO AF + DP: INFO = %q, want %q", got, "AF=0.5;DP=10")
	}
}

// TestCLI_RecodeINFOAlias_MixedWithKeepINFO — `--recode-INFO` and the
// port's pre-existing `--keep-INFO` flow into the same `keepINFOParts`
// slice in main.go; mixing them must union into a single set.
func TestCLI_RecodeINFOAlias_MixedWithKeepINFO(t *testing.T) {
	bin := buildVcftools(t)
	dir := t.TempDir()
	prefix := filepath.Join(dir, "out")
	_, _ = runWithStdin(t, bin, aliasVCF,
		"--stdin", "--recode",
		"--keep-INFO", "AF",
		"--recode-INFO", "AA",
		"--out", prefix)

	data, err := os.ReadFile(prefix + ".recode.vcf")
	if err != nil {
		t.Fatalf("read recode output: %v", err)
	}
	got := dataRowINFO(string(data))
	// Recoded INFO ordering: tools/vcftools/pkg/vcftools/vcftools.go does
	// not propagate InfoOrder onto the recoded variant, so the surviving
	// keys fall through to formatInfo's alphabetical-leftover branch
	// (pkg/bioformats/vcf/vcf.go:394-409). AA precedes AF alphabetically.
	if got != "AA=T;AF=0.5" {
		t.Errorf("mixed --keep-INFO AF + --recode-INFO AA: INFO = %q, want %q", got, "AA=T;AF=0.5")
	}
}
