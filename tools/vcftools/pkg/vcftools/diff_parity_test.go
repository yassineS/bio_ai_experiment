package vcftools

// Parity tests for the --diff family long-tail extensions:
//
//   --diff-indv-map FILE        — rename file-2 sample IDs before matching
//   --diff-discordance-matrix   — 4x4 genotype-by-genotype counts
//
// Goldens were generated from upstream vcftools (built from the submodule at
// reference_code/vcftools) using the documented FORTIFY_SOURCE workaround
// (PARITY_ROADMAP.md#vcftools, near "Upstream build note for golden
// generation"):
//
//	cd reference_code/vcftools/src/cpp
//	make clean
//	make CXXFLAGS='-O0 -g -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=0'
//
// All fixtures and goldens live under tools/vcftools/testdata/parity/.

import (
	"path/filepath"
	"testing"
)

// TestParity_DiffDiscordanceMatrix_NoMap exercises --diff-discordance-matrix
// with no sample-ID mapping. The two fixtures share no sample names, so the
// matrix is entirely zero — but the header/labels/row count must still match
// upstream byte-for-byte.
func TestParity_DiffDiscordanceMatrix_NoMap(t *testing.T) {
	fix := vcftoolsFixtureDir(t)
	prefix := runVcftoolsParity(t, "diff_f1.vcf", &Params{
		Diff:                  filepath.Join(fix, "diff_f2.vcf"),
		DiffDiscordanceMatrix: true,
	})
	got := readFileBytes(t, prefix+".diff.discordance_matrix")
	want := readFileBytes(t, filepath.Join(fix, "diff_dm_nomap.expected.diff.discordance_matrix"))
	if string(got) != string(want) {
		t.Errorf("discordance_matrix byte mismatch:\n--- got ---\n%s\n--- want ---\n%s",
			string(got), string(want))
	}
}

// TestParity_DiffDiscordanceMatrix_WithMap exercises --diff-discordance-matrix
// together with --diff-indv-map. The map renames file-2's a1/a2 to file-1's
// s1/s2 so the matrix is populated by the two shared (post-rename) samples
// across the three shared sites.
func TestParity_DiffDiscordanceMatrix_WithMap(t *testing.T) {
	fix := vcftoolsFixtureDir(t)
	prefix := runVcftoolsParity(t, "diff_f1.vcf", &Params{
		Diff:                  filepath.Join(fix, "diff_f2.vcf"),
		DiffIndvMap:           filepath.Join(fix, "diff_indv_map.txt"),
		DiffDiscordanceMatrix: true,
	})
	got := readFileBytes(t, prefix+".diff.discordance_matrix")
	want := readFileBytes(t, filepath.Join(fix, "diff_dm_mapped.expected.diff.discordance_matrix"))
	if string(got) != string(want) {
		t.Errorf("discordance_matrix byte mismatch:\n--- got ---\n%s\n--- want ---\n%s",
			string(got), string(want))
	}
}

// TestParity_DiffIndvInFiles_WithMap exercises --diff-indv (the
// indv_in_files report) under a --diff-indv-map. Upstream renames file-2's
// a1/a2 to s1/s2 so the output should list s1, s2 as "B" rather than as 1-
// only and 2-only.
func TestParity_DiffIndvInFiles_WithMap(t *testing.T) {
	fix := vcftoolsFixtureDir(t)
	prefix := runVcftoolsParity(t, "diff_f1.vcf", &Params{
		Diff:        filepath.Join(fix, "diff_f2.vcf"),
		DiffIndvMap: filepath.Join(fix, "diff_indv_map.txt"),
		DiffIndv:    true,
	})
	got := readFileBytes(t, prefix+".diff.indv_in_files")
	want := readFileBytes(t, filepath.Join(fix, "diff_indv_mapped.expected.diff.indv_in_files"))
	if string(got) != string(want) {
		t.Errorf("indv_in_files byte mismatch:\n--- got ---\n%s\n--- want ---\n%s",
			string(got), string(want))
	}
}

// TestParity_DiffSiteDiscordance_NoMap exercises --diff-site-discordance
// without a sample-ID mapping. The two fixtures share no sample names, so
// every site has N_COMMON_CALLED = 0 and DISCORDANCE = -nan. The 7-column
// layout (CHROM, POS, FILES, MATCHING_ALLELES, N_COMMON_CALLED, N_DISCORD,
// DISCORDANCE) and the per-site MATCHING_ALLELES values still need to
// match upstream byte-for-byte.
func TestParity_DiffSiteDiscordance_NoMap(t *testing.T) {
	fix := vcftoolsFixtureDir(t)
	prefix := runVcftoolsParity(t, "diff_f1.vcf", &Params{
		Diff:                filepath.Join(fix, "diff_f2.vcf"),
		DiffSiteDiscordance: true,
	})
	got := readFileBytes(t, prefix+".diff.sites")
	want := readFileBytes(t, filepath.Join(fix, "diff_sites_nomap.expected.diff.sites"))
	if string(got) != string(want) {
		t.Errorf(".diff.sites byte mismatch:\n--- got ---\n%s\n--- want ---\n%s",
			string(got), string(want))
	}
}

// TestParity_DiffSiteDiscordance_WithMap exercises --diff-site-discordance
// together with --diff-indv-map. With a1→s1 and a2→s2 mapping, the two
// shared samples drive non-zero N_COMMON_CALLED and DISCORDANCE values
// (0.5, 0, 1 across the three shared sites).
func TestParity_DiffSiteDiscordance_WithMap(t *testing.T) {
	fix := vcftoolsFixtureDir(t)
	prefix := runVcftoolsParity(t, "diff_f1.vcf", &Params{
		Diff:                filepath.Join(fix, "diff_f2.vcf"),
		DiffIndvMap:         filepath.Join(fix, "diff_indv_map.txt"),
		DiffSiteDiscordance: true,
	})
	got := readFileBytes(t, prefix+".diff.sites")
	want := readFileBytes(t, filepath.Join(fix, "diff_sites_mapped.expected.diff.sites"))
	if string(got) != string(want) {
		t.Errorf(".diff.sites byte mismatch:\n--- got ---\n%s\n--- want ---\n%s",
			string(got), string(want))
	}
}

// TestParity_DiffIndvDiscordance_NoMap exercises --diff-indv-discordance
// without a sample-ID mapping. Without the map, file-1 (s1,s2,s3) and
// file-2 (a1,a2,a4) share no sample names, so every per-individual row
// is 0/0/-nan. The output lists the union of both files' samples in
// alphabetical order, matching upstream's std::map iteration order.
func TestParity_DiffIndvDiscordance_NoMap(t *testing.T) {
	fix := vcftoolsFixtureDir(t)
	prefix := runVcftoolsParity(t, "diff_f1.vcf", &Params{
		Diff:                filepath.Join(fix, "diff_f2.vcf"),
		DiffIndvDiscordance: true,
	})
	got := readFileBytes(t, prefix+".diff.indv")
	want := readFileBytes(t, filepath.Join(fix, "diff_indv_nomap.expected.diff.indv"))
	if string(got) != string(want) {
		t.Errorf(".diff.indv byte mismatch:\n--- got ---\n%s\n--- want ---\n%s",
			string(got), string(want))
	}
}

// TestParity_DiffSiteDiscordance_AltMismatch pins MATCHING_ALLELES=0 for
// a shared B-site whose ALT alleles differ between the two files (REF=A
// in both, ALT=G in file-1, ALT=C in file-2). Upstream's
// `alleles_match = (ALT1 == ALT2) && (REF1 == REF2)`
// (variant_file_diff.cpp:844) flags this row as a 0 — the discordance
// counts still accumulate because genotype IDs are looked up via the
// "alleles don't match" string-comparison branch
// (variant_file_diff.cpp:887-916). The same fixture also pins
// MATCHING_ALLELES=1 for a site where both REF and ALT line up.
func TestParity_DiffSiteDiscordance_AltMismatch(t *testing.T) {
	fix := vcftoolsFixtureDir(t)
	prefix := runVcftoolsParity(t, "diff_alt_mismatch_f1.vcf", &Params{
		Diff:                filepath.Join(fix, "diff_alt_mismatch_f2.vcf"),
		DiffSiteDiscordance: true,
	})
	got := readFileBytes(t, prefix+".diff.sites")
	want := readFileBytes(t, filepath.Join(fix, "diff_alt_mismatch.expected.diff.sites"))
	if string(got) != string(want) {
		t.Errorf(".diff.sites byte mismatch:\n--- got ---\n%s\n--- want ---\n%s",
			string(got), string(want))
	}
}

// TestParity_DiffIndvDiscordance_WithMap exercises --diff-indv-discordance
// under a --diff-indv-map. The mapping renames file-2's a1/a2 onto s1/s2
// so those samples carry real discordance counts (1/3, 2/3); a4 stays
// unmapped (file-2 only → -nan) and s3 is file-1 only (also -nan).
func TestParity_DiffIndvDiscordance_WithMap(t *testing.T) {
	fix := vcftoolsFixtureDir(t)
	prefix := runVcftoolsParity(t, "diff_f1.vcf", &Params{
		Diff:                filepath.Join(fix, "diff_f2.vcf"),
		DiffIndvMap:         filepath.Join(fix, "diff_indv_map.txt"),
		DiffIndvDiscordance: true,
	})
	got := readFileBytes(t, prefix+".diff.indv")
	want := readFileBytes(t, filepath.Join(fix, "diff_indv_mapped.expected.diff.indv"))
	if string(got) != string(want) {
		t.Errorf(".diff.indv byte mismatch:\n--- got ---\n%s\n--- want ---\n%s",
			string(got), string(want))
	}
}
