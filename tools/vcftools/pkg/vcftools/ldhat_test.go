package vcftools

// Parity tests for --ldhat / --ldhat-geno / --phased.
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

// TestParity_LDhatPhased exercises --ldhat on a fully-phased 3-site / 3-sample
// fixture (ldhat_phased.vcf) restricted to chr 20. Both the .ldhat.sites and
// .ldhat.locs outputs must match upstream byte-for-byte.
func TestParity_LDhatPhased(t *testing.T) {
	prefix := runVcftoolsParity(t, "ldhat_phased.vcf", &Params{
		Chr:   "20",
		LDhat: true,
	})

	cases := []struct {
		gotPath  string
		wantPath string
	}{
		{prefix + ".ldhat.sites", "ldhat_phased.expected.ldhat.sites"},
		{prefix + ".ldhat.locs", "ldhat_phased.expected.ldhat.locs"},
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

// TestParity_LDhatGeno exercises --ldhat-geno on the same fully-phased
// fixture. The site row layout switches to one row per sample (N rows,
// state-count 2) and the per-genotype encoding follows upstream's
// variant_file_format_convert.cpp:903-936 switch.
func TestParity_LDhatGeno(t *testing.T) {
	prefix := runVcftoolsParity(t, "ldhat_phased.vcf", &Params{
		Chr:       "20",
		LDhatGeno: true,
	})

	cases := []struct {
		gotPath  string
		wantPath string
	}{
		{prefix + ".ldhat.sites", "ldhat_geno.expected.ldhat.sites"},
		{prefix + ".ldhat.locs", "ldhat_geno.expected.ldhat.locs"},
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

// TestParity_LDhatGenoHaploid runs --ldhat-geno on a chrX fixture that mixes
// diploid phased calls, haploid calls (e.g. "1"), and a haploid missing call
// ("."). Upstream maps:
//
//	"."  -> '?'    haploid missing
//	"0"  -> '0'    haploid ref
//	"1"  -> '1'    haploid alt
//
// We override MaxMissing to 0 so the existing "no missing allowed" default
// (which differs from upstream's --max-missing default 0) doesn't filter out
// the haploid-missing site.
func TestParity_LDhatGenoHaploid(t *testing.T) {
	prefix := runVcftoolsParity(t, "ldhat_haploid_x.vcf", &Params{
		Chr:        "X",
		LDhatGeno:  true,
		MaxMissing: 0,
	})

	cases := []struct {
		gotPath  string
		wantPath string
	}{
		{prefix + ".ldhat.sites", "ldhat_geno_haploid.expected.ldhat.sites"},
		{prefix + ".ldhat.locs", "ldhat_geno_haploid.expected.ldhat.locs"},
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

// TestParity_PhasedDropsUnphasedSites exercises --phased on a mixed-phasing
// fixture (sites with one unphased GT must be dropped). The expected output
// (matching upstream) is two of the four sites surviving. We compare the body
// lines of the recoded VCF rather than the full file because our header copy
// preserves the source meta-info verbatim and the kept rows must match.
func TestParity_PhasedDropsUnphasedSites(t *testing.T) {
	prefix := runVcftoolsParity(t, "phased_mixed.vcf", &Params{
		Phased: true,
		Recode: true,
	})

	gotPath := prefix + ".recode.vcf"
	got, err := os.ReadFile(gotPath)
	if err != nil {
		t.Fatalf("read %s: %v", gotPath, err)
	}

	// Extract non-meta lines.
	var bodyLines []string
	for _, line := range strings.Split(string(got), "\n") {
		if line == "" || strings.HasPrefix(line, "##") {
			continue
		}
		bodyLines = append(bodyLines, line)
	}

	want := []string{
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3",
		"20\t1000\t.\tA\tG\t30\tPASS\t.\tGT\t0|0\t0|1\t1|1",
		"20\t2000\t.\tG\tA\t50\tPASS\t.\tGT\t1|1\t0|0\t0|1",
	}
	assertLinesEqual(t, bodyLines, want)
}

// TestParity_LDhatImpliesPhased confirms that --ldhat composes with the
// site-level phase filter the way upstream does: a fixture with one unphased
// site has that site dropped from the --ldhat output even though the
// underlying records are biallelic.
func TestParity_LDhatImpliesPhased(t *testing.T) {
	// phased_mixed.vcf is on chr 20; --ldhat requires --chr. Two of its four
	// sites are fully phased, so both --ldhat output files should reflect a
	// 2-site bundle.
	prefix := runVcftoolsParity(t, "phased_mixed.vcf", &Params{
		Chr:   "20",
		LDhat: true,
	})

	// We expect 2 sites at positions 1000 and 2000, max pos = 2000.
	// locs header: "<n>\t<max/1000>\tL\n"
	wantLocs := "2\t2.0000\tL\n1.0000\n2.0000\n"
	gotLocs := string(readFileBytes(t, prefix+".ldhat.locs"))
	if gotLocs != wantLocs {
		t.Errorf("locs mismatch\n--- got ---\n%s\n--- want ---\n%s", gotLocs, wantLocs)
	}

	// sites header: "<2*N>\t<n_sites>\t1\n" with N=3 samples -> 6 haplotypes,
	// 2 sites. Per-haplotype rows below were captured byte-for-byte from
	// upstream vcftools (--ldhat on phased_mixed.vcf with --chr 20).
	wantSites := "6\t2\t1\n" +
		">S1-0\n01\n" +
		">S1-1\n01\n" +
		">S2-0\n00\n" +
		">S2-1\n10\n" +
		">S3-0\n10\n" +
		">S3-1\n11\n"
	gotSites := string(readFileBytes(t, prefix+".ldhat.sites"))
	if gotSites != wantSites {
		t.Errorf("sites mismatch\n--- got ---\n%s\n--- want ---\n%s", gotSites, wantSites)
	}
}

// TestLDhat_RequiresChr verifies the upstream parameters.cpp:717 invariant:
// --ldhat / --ldhat-geno must be paired with --chr, otherwise the run fails
// with the upstream error string.
func TestLDhat_RequiresChr(t *testing.T) {
	cases := []struct {
		name   string
		params *Params
	}{
		{"ldhat without --chr", &Params{LDhat: true}},
		{"ldhat-geno without --chr", &Params{LDhatGeno: true}},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			tmp := t.TempDir()
			c.params.OutPrefix = filepath.Join(tmp, "out")
			in, err := os.Open(filepath.Join(vcftoolsFixtureDir(t), "ldhat_phased.vcf"))
			if err != nil {
				t.Fatalf("open vcf: %v", err)
			}
			defer in.Close()
			err = Run(in, c.params)
			if err == nil {
				t.Fatalf("expected error, got nil")
			}
			if !strings.Contains(err.Error(), "Require a chromosome (--chr)") {
				t.Errorf("unexpected error: %v", err)
			}
		})
	}
}

// TestLDhat_MultiallelicSkipped checks the upstream one-off behaviour: only
// biallelic loci contribute to the LDhat output. A multi-allelic site at the
// same chrom should be silently dropped (with a single stderr warning, which
// we don't assert here).
func TestLDhat_MultiallelicSkipped(t *testing.T) {
	// Hand-built in-process: one biallelic site, one tri-allelic site.
	vcfText := strings.Join([]string{
		"##fileformat=VCFv4.2",
		"##contig=<ID=20>",
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2",
		"20\t100\t.\tA\tG\t30\tPASS\t.\tGT\t0|0\t1|1",
		"20\t200\t.\tC\tT,G\t40\tPASS\t.\tGT\t0|1\t1|2",
		"20\t300\t.\tG\tA\t50\tPASS\t.\tGT\t0|1\t1|0",
		"",
	}, "\n")

	tmp := t.TempDir()
	in := strings.NewReader(vcfText)
	prefix := filepath.Join(tmp, "out")

	// Override MinAlleles/MaxAlleles to 0 so the existing biallelic filter
	// doesn't drop the tri-allelic site before --ldhat sees it (otherwise we
	// can't tell whether the LDhat-level filter is doing its job).
	err := Run(in, &Params{
		OutPrefix:  prefix,
		Chr:        "20",
		LDhatGeno:  true,
		MinAlleles: 0,
		MaxAlleles: 0,
	})
	if err != nil {
		t.Fatalf("Run: %v", err)
	}

	gotSites := string(readFileBytes(t, prefix+".ldhat.sites"))
	// 2 samples, 2 biallelic sites kept (the tri-allelic at pos 200 dropped).
	// State count 2 for unphased layout.
	wantSites := "2\t2\t2\n>S1\n02\n>S2\n12\n"
	if gotSites != wantSites {
		t.Errorf("sites mismatch\n--- got ---\n%s\n--- want ---\n%s", gotSites, wantSites)
	}

	gotLocs := string(readFileBytes(t, prefix+".ldhat.locs"))
	// 2 sites; max pos 300 -> 0.3 kb.
	wantLocs := "2\t0.3000\tL\n0.1000\n0.3000\n"
	if gotLocs != wantLocs {
		t.Errorf("locs mismatch\n--- got ---\n%s\n--- want ---\n%s", gotLocs, wantLocs)
	}
}

// TestParseGTForLDhat is a table-driven unit test for the GT parser used by
// the LDhat writers. Mirrors upstream's
// vcf_entry::set_indv_GENOTYPE_and_PHASE (vcf_entry_setters.cpp:67-101).
func TestParseGTForLDhat(t *testing.T) {
	cases := []struct {
		gt             string
		wantA1, wantA2 int
		wantPhased     bool
	}{
		{"0|0", 0, 0, true},
		{"0|1", 0, 1, true},
		{"1|0", 1, 0, true},
		{"1|1", 1, 1, true},
		{"0/1", 0, 1, false},
		{"1/0", 1, 0, false},
		{"./.", -1, -1, false},
		{".|.", -1, -1, true},
		{".", -1, -2, true},   // haploid missing
		{"0", 0, -2, true},    // haploid ref
		{"1", 1, -2, true},    // haploid alt
		{"0|.", 0, -1, true},  // diploid with trailing missing
		{"./0", -1, 0, false}, // diploid with leading missing
		{"", -1, -2, true},    // empty -> haploid missing
		{"2|0", 2, 0, true},   // multi-allele index passes through
	}
	for _, c := range cases {
		t.Run(c.gt, func(t *testing.T) {
			a1, a2, p := parseGTForLDhat(c.gt)
			if a1 != c.wantA1 || a2 != c.wantA2 || p != c.wantPhased {
				t.Errorf("parseGTForLDhat(%q) = (%d,%d,%v); want (%d,%d,%v)",
					c.gt, a1, a2, p, c.wantA1, c.wantA2, c.wantPhased)
			}
		})
	}
}

// TestLDhatUnphasedChar enumerates the upstream switch in
// variant_file_format_convert.cpp:903-936 so we can detect a regression
// without round-tripping through Run.
func TestLDhatUnphasedChar(t *testing.T) {
	cases := []struct {
		a1, a2 int
		phased bool
		want   byte
	}{
		// Homozygous ref / het / homozygous alt
		{0, 0, false, '0'},
		{0, 0, true, '0'},
		{0, 1, false, '2'},
		{1, 0, false, '2'},
		{1, 1, false, '1'},
		// Phased haploid via trailing "." in diploid GT: alleles.second == -1
		// + PHASE='|'.
		{0, -1, true, '0'},
		{1, -1, true, '1'},
		// Unphased trailing-. — upstream falls through to '?'.
		{0, -1, false, '?'},
		{1, -1, false, '?'},
		// Truly haploid: alleles.second == -2 -> always mapped.
		{0, -2, true, '0'},
		{1, -2, true, '1'},
		{0, -2, false, '0'}, // upstream's mapping doesn't gate on phase here
		// Missing / out-of-range.
		{-1, -1, false, '?'},
		{-1, 0, true, '?'},
		{2, 0, true, '?'},
		{0, 2, true, '?'},
	}
	for _, c := range cases {
		got := ldhatUnphasedChar(c.a1, c.a2, c.phased)
		if got != c.want {
			t.Errorf("ldhatUnphasedChar(%d,%d,%v) = %q; want %q",
				c.a1, c.a2, c.phased, string(got), string(c.want))
		}
	}
}
