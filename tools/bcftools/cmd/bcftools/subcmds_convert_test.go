package main

import "testing"

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
		// --tsv2vcf, --gvcf2vcf and -c/--columns are now implemented and
		// dispatched in runConvert before this gate, so they are no longer
		// part of the deferred set. Bare -f/--fasta-ref (without those two
		// modes) still implies an unimplemented gen/hap export.
		{"fasta-ref", checkConvertDeferredInputs{fastaRef: "ref.fa"}, "-f/--fasta-ref"},
		{"gvcf", checkConvertDeferredInputs{gvcfBlocks: "10,20"}, "--gvcf"},
		{"gensample", checkConvertDeferredInputs{gensample: "x"}, "-g/--gensample"},
		{"gensample2vcf", checkConvertDeferredInputs{gensample2vcf: "x"}, "-G/--gensample2vcf"},
		{"3N6", checkConvertDeferredInputs{threeN6: true}, "--3N6"},
		{"tag", checkConvertDeferredInputs{tagFlag: "GT"}, "--tag"},
		{"chrom", checkConvertDeferredInputs{chromFlag: "chr1"}, "--chrom"},
		{"keep-duplicates", checkConvertDeferredInputs{keepDuplicates: true}, "--keep-duplicates"},
		{"sex", checkConvertDeferredInputs{sexFlag: "M"}, "--sex"},
		{"vcf-ids", checkConvertDeferredInputs{vcfIds: true}, "--vcf-ids"},
		{"hapsample", checkConvertDeferredInputs{hapsample: "x"}, "--hapsample"},
		{"hapsample2vcf", checkConvertDeferredInputs{hapsample2vcf: "x"}, "--hapsample2vcf"},
		{"haploid2diploid", checkConvertDeferredInputs{haploid2diploid: true}, "--haploid2diploid"},
		{"haplegendsample", checkConvertDeferredInputs{haplegendsample: "x"}, "--haplegendsample"},
		{"haplegendsample2vcf", checkConvertDeferredInputs{haplegendsample2vcf: "x"}, "--haplegendsample2vcf"},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			if got := checkConvertDeferred(tc.in); got != tc.want {
				t.Errorf("deferred(%s): got %q, want %q", tc.name, got, tc.want)
			}
		})
	}
}
