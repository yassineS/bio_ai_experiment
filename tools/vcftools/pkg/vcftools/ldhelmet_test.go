package vcftools

// Parity tests for --ldhelmet.
//
// Goldens were generated from upstream vcftools built from
// reference_code/vcftools (commit pinned by the submodule) with the local
// CXXFLAGS workaround for the upstream stack-buffer overflow in
// variant_file_format_convert.cpp (`char tmpname[new_tmp.size()]` plus
// strcpy past the end). The workaround is:
//
//	cd reference_code/vcftools/src/cpp
//	make clean
//	make CXXFLAGS='-O0 -g -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=0'
//
// All fixtures and goldens live under tools/vcftools/testdata/parity/.

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// TestParity_LDhelmetPhased exercises --ldhelmet on a fully-phased 3-site /
// 3-sample fixture (ldhat_phased.vcf) restricted to chr 20. Both the
// .ldhelmet.snps and .ldhelmet.pos outputs must match upstream byte-for-byte.
func TestParity_LDhelmetPhased(t *testing.T) {
	prefix := runVcftoolsParity(t, "ldhat_phased.vcf", &Params{
		Chr:      "20",
		LDhelmet: true,
	})

	cases := []struct {
		gotPath  string
		wantPath string
	}{
		{prefix + ".ldhelmet.snps", "ldhelmet_phased.expected.ldhelmet.snps"},
		{prefix + ".ldhelmet.pos", "ldhelmet_phased.expected.ldhelmet.pos"},
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

// TestParity_LDhelmetMixedPhasing exercises the upstream phased_only
// invariant: a fixture with two unphased sites among four total should yield
// a 2-site output bundle (positions 1000 and 2000 — the fully phased ones).
func TestParity_LDhelmetMixedPhasing(t *testing.T) {
	prefix := runVcftoolsParity(t, "phased_mixed.vcf", &Params{
		Chr:      "20",
		LDhelmet: true,
	})

	cases := []struct {
		gotPath  string
		wantPath string
	}{
		{prefix + ".ldhelmet.snps", "ldhelmet_mixed.expected.ldhelmet.snps"},
		{prefix + ".ldhelmet.pos", "ldhelmet_mixed.expected.ldhelmet.pos"},
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

// TestParity_LDhelmetMissingGenotypes exercises the "N" fallback when an
// allele index is negative (missing). Sites with a phased ".|." or "0|."
// still pass filter_sites_by_phase (PHASE == '|' for both branches), so they
// reach the output writer; the writer then emits "N" for the missing
// haplotype. Matches upstream variant_file_format_convert.cpp:1110-1114.
func TestParity_LDhelmetMissingGenotypes(t *testing.T) {
	prefix := runVcftoolsParity(t, "ldhelmet_with_missing.vcf", &Params{
		Chr:      "20",
		LDhelmet: true,
		// MaxMissing default is 1 in our port, but the file has only
		// missing-allele entries in some genotypes (the GT separator
		// is '|' so the site as a whole is fully phased). Set
		// MaxMissing=0 to disable our generic per-site missing filter
		// so the missing-allele encoding path is reachable.
		MaxMissing: 0,
	})

	cases := []struct {
		gotPath  string
		wantPath string
	}{
		{prefix + ".ldhelmet.snps", "ldhelmet_missing.expected.ldhelmet.snps"},
		{prefix + ".ldhelmet.pos", "ldhelmet_missing.expected.ldhelmet.pos"},
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

// TestLDhelmet_RequiresChr verifies the upstream parameters.cpp:717
// invariant: --ldhelmet must be paired with --chr, otherwise the run fails
// with the upstream error string.
func TestLDhelmet_RequiresChr(t *testing.T) {
	tmp := t.TempDir()
	params := &Params{
		OutPrefix: filepath.Join(tmp, "out"),
		LDhelmet:  true,
	}
	in, err := os.Open(filepath.Join(vcftoolsFixtureDir(t), "ldhat_phased.vcf"))
	if err != nil {
		t.Fatalf("open vcf: %v", err)
	}
	defer in.Close()
	if err := Run(in, params); err == nil {
		t.Fatalf("expected error, got nil")
	} else if !strings.Contains(err.Error(), "Require a chromosome (--chr)") {
		t.Errorf("unexpected error: %v", err)
	}
}

// TestLDhelmet_ImpliesRemoveIndels confirms that --ldhelmet drops indel
// sites the way upstream's parameters.cpp:275 does (remove_indels = true).
func TestLDhelmet_ImpliesRemoveIndels(t *testing.T) {
	// Two sites: one SNP, one insertion (REF=A, ALT=AT).
	vcfText := strings.Join([]string{
		"##fileformat=VCFv4.2",
		"##contig=<ID=20>",
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
		"20\t100\t.\tA\tG\t30\tPASS\t.\tGT\t0|0\t1|1",
		"20\t200\t.\tA\tAT\t40\tPASS\t.\tGT\t0|1\t1|0",
		"20\t300\t.\tG\tC\t50\tPASS\t.\tGT\t0|1\t1|0",
		"",
	}, "\n")

	tmp := t.TempDir()
	in := strings.NewReader(vcfText)
	prefix := filepath.Join(tmp, "out")

	err := Run(in, &Params{
		OutPrefix: prefix,
		Chr:       "20",
		LDhelmet:  true,
	})
	if err != nil {
		t.Fatalf("Run: %v", err)
	}

	gotPos := string(readFileBytes(t, prefix+".ldhelmet.pos"))
	wantPos := "100\n300\n"
	if gotPos != wantPos {
		t.Errorf("pos mismatch\n--- got ---\n%s\n--- want ---\n%s", gotPos, wantPos)
	}

	gotSnps := string(readFileBytes(t, prefix+".ldhelmet.snps"))
	// Two samples, two SNP sites kept (insertion at pos 200 dropped).
	wantSnps := ">S1-0\nAG\n>S1-1\nAC\n>S2-0\nGC\n>S2-1\nGG\n"
	if gotSnps != wantSnps {
		t.Errorf("snps mismatch\n--- got ---\n%s\n--- want ---\n%s", gotSnps, wantSnps)
	}
}

// TestLDhelmetAllele unit-tests the per-allele lookup. Mirrors the upstream
// switch at variant_file_format_convert.cpp:1100-1114.
func TestLDhelmetAllele(t *testing.T) {
	alleles := []string{"A", "C", "G"}
	cases := []struct {
		idx  int
		want string
	}{
		{0, "A"},
		{1, "C"},
		{2, "G"},
		{-1, "N"},
		{-2, "N"},
		{3, "N"},
	}
	for _, c := range cases {
		if got := ldhelmetAllele(c.idx, alleles); got != c.want {
			t.Errorf("ldhelmetAllele(%d, %v) = %q; want %q", c.idx, alleles, got, c.want)
		}
	}
}
