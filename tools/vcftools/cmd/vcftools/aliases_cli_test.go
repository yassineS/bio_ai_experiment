// Wave-16/17 CLI-binary tests for upstream-canonical flag aliases:
//
//   - `-c` — upstream's short alias for `--stdout` (parameters.cpp:194).
//   - `--recode-INFO TAG` — upstream's canonical name for the repeatable
//     recode-INFO-column selector (parameters.cpp:319,
//     `recode_INFO_to_keep`). Wave 17 separated this from `--keep-INFO`
//     (which was previously routed to the same internal slice as a
//     misaligned synonym); `--recode-INFO` is now the only flag wired
//     to the recode-column selector path.
//   - `--keep-INFO TAG` — upstream's site filter
//     (parameters.cpp:266 → `site_INFO_flags_to_keep` →
//     entry_filters.cpp:1033). Pre-wave-17 the port routed this to the
//     recode-column selector; wave-17 swapped it to its upstream SITE
//     FILTER semantic. See pkg/vcftools/info_filters_test.go for
//     coverage of the new semantic; this file covers the CLI surface.
//
// All three flags are wired in main.go, so the tests live here next to
// the CLI to exercise the actual argv → flag.Parse path rather than the
// package boundary.
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
// column in `.recode.vcf` to the listed tags
// (parameters.cpp:319 → recode_INFO_to_keep).
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

// keepInfoSiteVCF is a small fixture with one Flag-type INFO key,
// one site with the flag set and one without. Used by the wave-17
// `--keep-INFO` CLI tests below.
const keepInfoSiteVCF = "##fileformat=VCFv4.2\n" +
	"##INFO=<ID=FLAG_A,Number=0,Type=Flag,Description=\"\">\n" +
	"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"\">\n" +
	"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"\">\n" +
	"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n" +
	"1\t100\t.\tA\tG\t.\tPASS\tFLAG_A;DP=10\tGT\t0/0\n" +
	"1\t200\t.\tA\tC\t.\tPASS\tDP=20\tGT\t0/1\n"

// TestCLI_KeepINFOSiteFilter — wave-17 fix-on-port. `--keep-INFO TAG`
// is now a SITE FILTER (upstream parameters.cpp:266 +
// entry_filters.cpp:1033-1063). It must drop sites where the named
// INFO Flag is not present.
func TestCLI_KeepINFOSiteFilter(t *testing.T) {
	bin := buildVcftools(t)
	dir := t.TempDir()
	prefix := filepath.Join(dir, "out")
	_, _ = runWithStdin(t, bin, keepInfoSiteVCF,
		"--stdin", "--recode", "--recode-INFO-all",
		"--keep-INFO", "FLAG_A",
		"--out", prefix)

	data, err := os.ReadFile(prefix + ".recode.vcf")
	if err != nil {
		t.Fatalf("read recode output: %v", err)
	}
	body := string(data)
	if !strings.Contains(body, "\t100\t") {
		t.Errorf("expected site 100 (FLAG_A) kept; got:\n%s", body)
	}
	if strings.Contains(body, "\t200\t") {
		t.Errorf("expected site 200 (no FLAG_A) dropped; got:\n%s", body)
	}
}

// TestCLI_KeepINFOAndRecodeINFOAreSeparate — the two flags drive
// independent Params fields (KeepINFO vs RecodeINFO) and may be
// combined: --keep-INFO filters sites, --recode-INFO restricts the
// INFO column of survivors.
func TestCLI_KeepINFOAndRecodeINFOAreSeparate(t *testing.T) {
	bin := buildVcftools(t)
	dir := t.TempDir()
	prefix := filepath.Join(dir, "out")
	_, _ = runWithStdin(t, bin, keepInfoSiteVCF,
		"--stdin", "--recode",
		"--keep-INFO", "FLAG_A",
		"--recode-INFO", "DP",
		"--out", prefix)

	data, err := os.ReadFile(prefix + ".recode.vcf")
	if err != nil {
		t.Fatalf("read recode output: %v", err)
	}
	body := string(data)
	if !strings.Contains(body, "\t100\t") {
		t.Errorf("expected site 100 kept; got:\n%s", body)
	}
	if strings.Contains(body, "\t200\t") {
		t.Errorf("expected site 200 dropped; got:\n%s", body)
	}
	got := dataRowINFO(body)
	if got != "DP=10" {
		t.Errorf("INFO column = %q, want %q (only DP survives via --recode-INFO)", got, "DP=10")
	}
}

// TestCLI_RemoveINFOSiteFilter — wave-18 fix-on-port. `--remove-INFO TAG`
// is now a SITE FILTER (upstream parameters.cpp:328 +
// entry_filters.cpp:1068-1086). It must drop sites where the named
// INFO Flag IS present (polarity-inverted from --keep-INFO).
func TestCLI_RemoveINFOSiteFilter(t *testing.T) {
	bin := buildVcftools(t)
	dir := t.TempDir()
	prefix := filepath.Join(dir, "out")
	_, _ = runWithStdin(t, bin, keepInfoSiteVCF,
		"--stdin", "--recode", "--recode-INFO-all",
		"--remove-INFO", "FLAG_A",
		"--out", prefix)

	data, err := os.ReadFile(prefix + ".recode.vcf")
	if err != nil {
		t.Fatalf("read recode output: %v", err)
	}
	body := string(data)
	if strings.Contains(body, "\t100\t") {
		t.Errorf("expected site 100 (FLAG_A) dropped; got:\n%s", body)
	}
	if !strings.Contains(body, "\t200\t") {
		t.Errorf("expected site 200 (no FLAG_A) kept; got:\n%s", body)
	}
}
