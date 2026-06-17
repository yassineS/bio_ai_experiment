package main

import (
	"os"
	"path/filepath"
	"testing"
)

const genCLIVCF = `##fileformat=VCFv4.2
##contig=<ID=20,length=63025520>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
20	100	rs1	C	T	.	.	.	GT	0/0	0/1
`

// TestRunConvertGenSample drives the -g/-G dispatch end-to-end through the CLI
// runner, confirming the GEN/sample modes are wired up (no longer deferred).
func TestRunConvertGenSample(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "in.vcf")
	if err := os.WriteFile(in, []byte(genCLIVCF), 0o644); err != nil {
		t.Fatalf("write vcf: %v", err)
	}
	pfx := filepath.Join(dir, "out")
	if code := runConvert([]string{"-g", pfx, in}); code != 0 {
		t.Fatalf("convert -g exit=%d", code)
	}
	if _, err := os.Stat(pfx + ".gen.gz"); err != nil {
		t.Fatalf(".gen.gz not written: %v", err)
	}
	if _, err := os.Stat(pfx + ".samples"); err != nil {
		t.Fatalf(".samples not written: %v", err)
	}
	// Round-trip back to VCF via -G.
	out := filepath.Join(dir, "rt.vcf")
	if code := runConvert([]string{"-G", pfx, "-O", "v", "-o", out}); code != 0 {
		t.Fatalf("convert -G exit=%d", code)
	}
	if fi, err := os.Stat(out); err != nil || fi.Size() == 0 {
		t.Fatalf("round-trip VCF not written: err=%v", err)
	}
}

// TestRunConvertChromDeprecated locks in the upstream behaviour: --chrom is
// deprecated and exits non-zero pointing at --3N6.
func TestRunConvertChromDeprecated(t *testing.T) {
	if code := runConvert([]string{"--chrom", "-g", "x", "in.vcf"}); code == 0 {
		t.Fatalf("--chrom should be rejected (deprecated)")
	}
}

// TestUnitConvertGvcfAlias confirms that `--gvcf` is treated as the
// upstream prefix-abbreviation of `--gvcf2vcf` (no longer a deferred
// block-output mode). Upstream getopt_long resolves `--gvcf` to the
// no-argument `--gvcf2vcf` flag, which requires -f/--fasta-ref; our CLI
// reproduces that exact behaviour. Runs with no upstream binary.
func TestUnitConvertGvcfAlias(t *testing.T) {
	dir := t.TempDir()
	in := filepath.Join(dir, "in.vcf")
	if err := os.WriteFile(in, []byte(genCLIVCF), 0o644); err != nil {
		t.Fatalf("write vcf: %v", err)
	}
	// Without -f/--fasta-ref both --gvcf and --gvcf2vcf must fail with the
	// same "requires the --fasta-ref option" diagnostic (exit 1), proving
	// --gvcf dispatches into the gVCF->VCF expansion rather than a deferred
	// gate (which previously returned exit 2).
	if code := runConvert([]string{"--gvcf", in}); code != 1 {
		t.Fatalf("convert --gvcf without -f: exit=%d, want 1", code)
	}
	if code := runConvert([]string{"--gvcf2vcf", in}); code != 1 {
		t.Fatalf("convert --gvcf2vcf without -f: exit=%d, want 1", code)
	}
}
