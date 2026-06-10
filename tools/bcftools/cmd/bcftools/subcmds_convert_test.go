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

// TestCheckConvertDeferred locks in the upstream-flag-name surface that
// runConvert hard-rejects rather than silently accepting. Per the
// project's "every documented upstream flag must be recognised — either
// implemented or gracefully rejected with a pointer at PARITY_ROADMAP.md"
// rule (docs/PARITY_ROADMAP.md#definition-of-11), a future refactor that
// drops any of these from the rejection set without implementing the
// underlying behaviour is a regression.
func TestCheckConvertDeferred(t *testing.T) {
	if got := checkConvertDeferred(checkConvertDeferredInputs{}); got != "" {
		t.Fatalf("empty inputs: got deferred=%q, want \"\"", got)
	}
	cases := []struct {
		name string
		in   checkConvertDeferredInputs
		want string
	}{
		{"gvcf2vcf", checkConvertDeferredInputs{gvcf2vcf: true}, "--gvcf2vcf"},
		{"fasta-ref", checkConvertDeferredInputs{fastaRef: "ref.fa"}, "-f/--fasta-ref"},
		{"gvcf", checkConvertDeferredInputs{gvcfBlocks: "10,20"}, "--gvcf"},
		{"keep-duplicates", checkConvertDeferredInputs{keepDuplicates: true}, "--keep-duplicates"},
		{"hapsample", checkConvertDeferredInputs{hapsample: "x"}, "--hapsample"},
		{"hapsample2vcf", checkConvertDeferredInputs{hapsample2vcf: "x"}, "--hapsample2vcf"},
		{"haploid2diploid", checkConvertDeferredInputs{haploid2diploid: true}, "--haploid2diploid"},
		{"haplegendsample", checkConvertDeferredInputs{haplegendsample: "x"}, "--haplegendsample"},
		{"haplegendsample2vcf", checkConvertDeferredInputs{haplegendsample2vcf: "x"}, "--haplegendsample2vcf"},
		{"tsv2vcf", checkConvertDeferredInputs{tsv2vcf: "x"}, "--tsv2vcf"},
		{"columns", checkConvertDeferredInputs{columnsFlag: "CHROM,POS"}, "-c/--columns"},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			if got := checkConvertDeferred(tc.in); got != tc.want {
				t.Errorf("deferred(%s): got %q, want %q", tc.name, got, tc.want)
			}
		})
	}
}
