package vcftools

// Parity tests for --IMPUTE.
//
// Goldens were generated from upstream vcftools built from
// reference_code/vcftools (commit pinned by the submodule). The IMPUTE
// writer doesn't trigger the upstream stack-overflow bug (it doesn't use
// the LDhat-style temp files), but rebuilding with the FORTIFY_SOURCE
// workaround is harmless and matches the ldhat goldens' build:
//
//	cd reference_code/vcftools/src/cpp
//	make clean
//	make CXXFLAGS='-O0 -g -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=0'
//
// All fixtures and goldens live under tools/vcftools/testdata/parity/.

import (
	"path/filepath"
	"strings"
	"testing"
)

// TestParity_ImputePhased exercises --IMPUTE on a fully-phased 3-site /
// 3-sample fixture. All three sites pass the biallelic + no-missing + phased
// invariants, so the .impute.legend / .impute.hap / .impute.hap.indv outputs
// must match upstream byte-for-byte.
func TestParity_ImputePhased(t *testing.T) {
	prefix := runVcftoolsParity(t, "ldhat_phased.vcf", &Params{
		IMPUTE: true,
	})

	cases := []struct {
		gotPath  string
		wantPath string
	}{
		{prefix + ".impute.legend", "impute_phased.expected.impute.legend"},
		{prefix + ".impute.hap", "impute_phased.expected.impute.hap"},
		{prefix + ".impute.hap.indv", "impute_phased.expected.impute.hap.indv"},
	}
	for _, c := range cases {
		got := readFileBytes(t, c.gotPath)
		want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), c.wantPath))
		if string(got) != string(want) {
			t.Errorf("%s: byte mismatch\n--- got ---\n%s\n--- want ---\n%s",
				filepath.Base(c.gotPath), string(got), string(want))
		}
	}
}

// TestParity_ImputeMixedPhasing exercises the phased_only invariant: a
// fixture with two unphased sites among four total should yield a 2-site
// output bundle (positions 1000 and 2000 — the fully phased ones).
func TestParity_ImputeMixedPhasing(t *testing.T) {
	prefix := runVcftoolsParity(t, "phased_mixed.vcf", &Params{
		IMPUTE: true,
	})

	cases := []struct {
		gotPath  string
		wantPath string
	}{
		{prefix + ".impute.legend", "impute_mixed.expected.impute.legend"},
		{prefix + ".impute.hap", "impute_mixed.expected.impute.hap"},
		{prefix + ".impute.hap.indv", "impute_mixed.expected.impute.hap.indv"},
	}
	for _, c := range cases {
		got := readFileBytes(t, c.gotPath)
		want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), c.wantPath))
		if string(got) != string(want) {
			t.Errorf("%s: byte mismatch\n--- got ---\n%s\n--- want ---\n%s",
				filepath.Base(c.gotPath), string(got), string(want))
		}
	}
}

// TestParity_ImputeIDColumn exercises the ID-vs-default fallback at upstream
// lines 586-591: when the VCF ID column is ".", emit "CHROM-POS"; otherwise
// emit the ID verbatim. The fixture has one site with rs100 and one with ".".
func TestParity_ImputeIDColumn(t *testing.T) {
	prefix := runVcftoolsParity(t, "impute_id.vcf", &Params{
		IMPUTE: true,
	})

	cases := []struct {
		gotPath  string
		wantPath string
	}{
		{prefix + ".impute.legend", "impute_id.expected.impute.legend"},
		{prefix + ".impute.hap", "impute_id.expected.impute.hap"},
		{prefix + ".impute.hap.indv", "impute_id.expected.impute.hap.indv"},
	}
	for _, c := range cases {
		got := readFileBytes(t, c.gotPath)
		want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), c.wantPath))
		if string(got) != string(want) {
			t.Errorf("%s: byte mismatch\n--- got ---\n%s\n--- want ---\n%s",
				filepath.Base(c.gotPath), string(got), string(want))
		}
	}
}

// TestImpute_DropsMissingSite checks the inline per-site missing/unphased
// guard at upstream lines 556-584: any kept-sample missing or unphased GT
// disqualifies the site. The fixture has one all-phased biallelic site and
// one site with a missing GT — only the first should appear in the output.
func TestImpute_DropsMissingSite(t *testing.T) {
	vcfText := strings.Join([]string{
		"##fileformat=VCFv4.2",
		"##contig=<ID=20>",
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
		"20\t100\t.\tA\tG\t30\tPASS\t.\tGT\t0|0\t1|1",
		"20\t200\t.\tC\tT\t40\tPASS\t.\tGT\t0|1\t.|0",
		"",
	}, "\n")

	tmp := t.TempDir()
	in := strings.NewReader(vcfText)
	prefix := filepath.Join(tmp, "out")

	err := Run(in, &Params{
		OutPrefix:  prefix,
		IMPUTE:     true,
		MaxMissing: 0, // disable our generic missing filter to confirm
		// the IMPUTE-specific inline missing check fires.
	})
	if err != nil {
		t.Fatalf("Run: %v", err)
	}

	gotLegend := string(readFileBytes(t, prefix+".impute.legend"))
	wantLegend := "ID pos allele0 allele1\n20-100 100 A G\n"
	if gotLegend != wantLegend {
		t.Errorf("legend mismatch\n--- got ---\n%s\n--- want ---\n%s", gotLegend, wantLegend)
	}

	gotHap := string(readFileBytes(t, prefix+".impute.hap"))
	wantHap := "0 0 1 1\n"
	if gotHap != wantHap {
		t.Errorf("hap mismatch\n--- got ---\n%s\n--- want ---\n%s", gotHap, wantHap)
	}

	gotIndv := string(readFileBytes(t, prefix+".impute.hap.indv"))
	wantIndv := "S1\nS2\n"
	if gotIndv != wantIndv {
		t.Errorf("indv mismatch\n--- got ---\n%s\n--- want ---\n%s", gotIndv, wantIndv)
	}
}

// TestImpute_DropsUnphasedSite confirms that filter_sites_by_phase is
// applied (a site with one unphased GT is dropped before reaching the
// IMPUTE writer's inline check).
func TestImpute_DropsUnphasedSite(t *testing.T) {
	vcfText := strings.Join([]string{
		"##fileformat=VCFv4.2",
		"##contig=<ID=20>",
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
		"20\t100\t.\tA\tG\t30\tPASS\t.\tGT\t0|0\t1|1",
		"20\t200\t.\tC\tT\t40\tPASS\t.\tGT\t0|1\t1/0",
		"",
	}, "\n")

	tmp := t.TempDir()
	in := strings.NewReader(vcfText)
	prefix := filepath.Join(tmp, "out")

	err := Run(in, &Params{
		OutPrefix: prefix,
		IMPUTE:    true,
	})
	if err != nil {
		t.Fatalf("Run: %v", err)
	}

	gotLegend := string(readFileBytes(t, prefix+".impute.legend"))
	wantLegend := "ID pos allele0 allele1\n20-100 100 A G\n"
	if gotLegend != wantLegend {
		t.Errorf("legend mismatch\n--- got ---\n%s\n--- want ---\n%s", gotLegend, wantLegend)
	}
}

// TestImpute_MultiallelicSkipped checks that tri-allelic sites are silently
// dropped (with a single stderr warning, which we don't assert).
func TestImpute_MultiallelicSkipped(t *testing.T) {
	vcfText := strings.Join([]string{
		"##fileformat=VCFv4.2",
		"##contig=<ID=20>",
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
		"20\t100\t.\tA\tG\t30\tPASS\t.\tGT\t0|0\t1|1",
		"20\t200\t.\tC\tT,G\t40\tPASS\t.\tGT\t0|1\t1|2",
		"",
	}, "\n")

	tmp := t.TempDir()
	in := strings.NewReader(vcfText)
	prefix := filepath.Join(tmp, "out")

	// Override MinAlleles/MaxAlleles so the existing biallelic filter
	// doesn't drop the tri-allelic site before --IMPUTE sees it.
	err := Run(in, &Params{
		OutPrefix:  prefix,
		IMPUTE:     true,
		MinAlleles: 0,
		MaxAlleles: 0,
	})
	if err != nil {
		t.Fatalf("Run: %v", err)
	}

	gotLegend := string(readFileBytes(t, prefix+".impute.legend"))
	wantLegend := "ID pos allele0 allele1\n20-100 100 A G\n"
	if gotLegend != wantLegend {
		t.Errorf("legend mismatch\n--- got ---\n%s\n--- want ---\n%s", gotLegend, wantLegend)
	}
}
