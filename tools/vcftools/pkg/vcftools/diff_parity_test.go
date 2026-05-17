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
