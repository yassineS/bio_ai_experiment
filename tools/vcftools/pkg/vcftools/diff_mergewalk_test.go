package vcftools

// Merge-walk parity tests for the --diff family. These pin the row ORDERING
// and the POS-equal REF-mismatch / REF-N-normalisation behaviour against the
// upstream two-pointer merge in reference_code/vcftools/src/cpp/
// variant_file_diff.cpp (output_sites_in_files / output_discordance_by_site).
//
// Upstream vcftools does not build on this darwin/arm64 host, so these
// expectations are hand-traced from the C++ source rather than captured from a
// live oracle; the byte-exact-vs-upstream gate for these paths runs in the
// Linux CI/container. Every expectation here is deterministic (ours-before/
// after stable) and matches the documented upstream algorithm line-for-line.

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// runDiffBodies runs the --diff family over the two in-memory VCFs and returns
// the .diff.sites_in_files, .diff.sites and .diff.indv bodies (header stripped).
func runDiffBodies(t *testing.T, f1, f2 string) (sitesInFiles, sites, indv string) {
	t.Helper()
	tmp := t.TempDir()
	f2path := filepath.Join(tmp, "f2.vcf")
	if err := os.WriteFile(f2path, []byte(f2), 0o644); err != nil {
		t.Fatalf("write f2: %v", err)
	}
	prefix := filepath.Join(tmp, "cmp")
	if err := Run(strings.NewReader(f1), &Params{
		OutPrefix:           prefix,
		Diff:                f2path,
		DiffSite:            true,
		DiffSiteDiscordance: true,
		DiffIndvDiscordance: true,
	}); err != nil {
		t.Fatalf("Run --diff: %v", err)
	}
	read := func(suffix string) string {
		body, err := os.ReadFile(prefix + suffix)
		if err != nil {
			t.Fatalf("read %s: %v", suffix, err)
		}
		// Strip the header line.
		_, rest, _ := strings.Cut(string(body), "\n")
		return rest
	}
	return read(".diff.sites_in_files"), read(".diff.sites"), read(".diff.indv")
}

const mergeHdr = `##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
`

// TestDiffMergeWalk_Interleave pins that file-2-only sites are emitted
// INTERLEAVED with file-1 sites in ascending position order, not appended
// after all file-1 rows. file-1 has 100,300; file-2 has 100,200,300 — the
// file-2-only 200 must appear BETWEEN 100 and 300.
func TestDiffMergeWalk_Interleave(t *testing.T) {
	f1 := mergeHdr + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n" +
		"chr1\t100\t.\tA\tG\t30\tPASS\t.\tGT\t0/0\n" +
		"chr1\t300\t.\tG\tA\t30\tPASS\t.\tGT\t0/1\n"
	f2 := mergeHdr + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n" +
		"chr1\t100\t.\tA\tG\t30\tPASS\t.\tGT\t0/0\n" +
		"chr1\t200\t.\tC\tT\t30\tPASS\t.\tGT\t0/1\n" +
		"chr1\t300\t.\tG\tA\t30\tPASS\t.\tGT\t0/1\n"

	sif, sites, _ := runDiffBodies(t, f1, f2)

	// 100: 0/0 vs 0/0 concordant; 300: 0/1 vs 0/1 concordant. 200 is
	// file-2-only and must land between 100 and 300.
	wantSites := "chr1\t100\tB\t1\t1\t0\t0\n" +
		"chr1\t200\t2\t0\t0\t0\t-nan\n" +
		"chr1\t300\tB\t1\t1\t0\t0\n"
	if sites != wantSites {
		t.Errorf(".diff.sites order mismatch:\n--- got ---\n%s\n--- want ---\n%s", sites, wantSites)
	}

	wantSIF := "chr1\t100\t100\tB\tA\tA\tG\tG\n" +
		"chr1\t.\t200\t2\t.\tC\t.\tT\n" +
		"chr1\t300\t300\tB\tG\tG\tA\tA\n"
	if sif != wantSIF {
		t.Errorf(".diff.sites_in_files order mismatch:\n--- got ---\n%s\n--- want ---\n%s", sif, wantSIF)
	}
}

// TestDiffMergeWalk_RefMismatchSkip pins upstream's POS-equal REF-mismatch
// handling: .diff.sites SKIPS the row entirely (no counts, one-off warning),
// while .diff.sites_in_files emits an "O" (overlap) row. The site's genotypes
// must NOT contribute to the per-individual discordance.
func TestDiffMergeWalk_RefMismatchSkip(t *testing.T) {
	// 100: shared, REF match (B). 200: shared but REF differs (A vs C) →
	// .diff.sites skip, sites_in_files O. 300: shared, REF match (B).
	f1 := mergeHdr + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n" +
		"chr1\t100\t.\tA\tG\t30\tPASS\t.\tGT\t0/1\n" +
		"chr1\t200\t.\tA\tG\t30\tPASS\t.\tGT\t1/1\n" +
		"chr1\t300\t.\tG\tA\t30\tPASS\t.\tGT\t0/1\n"
	f2 := mergeHdr + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n" +
		"chr1\t100\t.\tA\tG\t30\tPASS\t.\tGT\t0/1\n" +
		"chr1\t200\t.\tC\tG\t30\tPASS\t.\tGT\t0/0\n" +
		"chr1\t300\t.\tG\tA\t30\tPASS\t.\tGT\t0/0\n"

	sif, sites, indv := runDiffBodies(t, f1, f2)

	// .diff.sites: 200 is absent (skipped). 100 concordant (0/0 discord),
	// 300 discordant (0/1 vs 0/0).
	wantSites := "chr1\t100\tB\t1\t1\t0\t0\n" +
		"chr1\t300\tB\t1\t1\t1\t1\n"
	if sites != wantSites {
		t.Errorf(".diff.sites (REF-mismatch skip) mismatch:\n--- got ---\n%s\n--- want ---\n%s", sites, wantSites)
	}

	// .diff.sites_in_files: 200 becomes an O row.
	wantSIF := "chr1\t100\t100\tB\tA\tA\tG\tG\n" +
		"chr1\t200\t200\tO\tA\tC\tG\tG\n" +
		"chr1\t300\t300\tB\tG\tG\tA\tA\n"
	if sif != wantSIF {
		t.Errorf(".diff.sites_in_files (REF-mismatch O) mismatch:\n--- got ---\n%s\n--- want ---\n%s", sif, wantSIF)
	}

	// s1 discordance: only 100 (concordant) and 300 (discordant) count; the
	// skipped 200 must NOT contribute → 2 common, 1 discord.
	if !strings.Contains(indv, "s1\t2\t1\t0.5") {
		t.Errorf(".diff.indv should exclude the skipped REF-mismatch site, got:\n%s", indv)
	}
}

// TestDiffMergeWalk_RefNormalisation pins upstream's "N"/"."/"" REF
// normalisation at POS-equal sites: an "N" REF on one side is replaced by the
// other file's REF before the match test, so the site is a matching B row (not
// an O/skip) and MATCHING_ALLELES reflects the ALT comparison only.
func TestDiffMergeWalk_RefNormalisation(t *testing.T) {
	// file-1 REF=N at 100; file-2 REF=A. After normalisation REF matches, so
	// the site is a B row with MATCHING_ALLELES=1 (ALTs also match).
	f1 := mergeHdr + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n" +
		"chr1\t100\t.\tN\tG\t30\tPASS\t.\tGT\t0/1\n"
	f2 := mergeHdr + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n" +
		"chr1\t100\t.\tA\tG\t30\tPASS\t.\tGT\t0/1\n"

	sif, sites, _ := runDiffBodies(t, f1, f2)

	wantSites := "chr1\t100\tB\t1\t1\t0\t0\n"
	if sites != wantSites {
		t.Errorf(".diff.sites (REF-N norm) mismatch:\n--- got ---\n%s\n--- want ---\n%s", sites, wantSites)
	}
	// REF1 is normalised to REF2 (A) in the sites_in_files output.
	wantSIF := "chr1\t100\t100\tB\tA\tA\tG\tG\n"
	if sif != wantSIF {
		t.Errorf(".diff.sites_in_files (REF-N norm) mismatch:\n--- got ---\n%s\n--- want ---\n%s", sif, wantSIF)
	}
}

// TestDiffMergeWalk_IndvStricterRefSkip pins the divergence between the
// .diff.sites REF normalisation and the stricter .diff.indv one: at a POS-equal
// site with REF="." in file-1 and REF="A" in file-2, output_discordance_by_site
// normalises "."→"A" and emits a B row WITH genotype counts, but
// output_discordance_by_indv normalises only "N" so it SKIPS the site — the
// per-sample totals must therefore exclude it (variant_file_diff.cpp:508-517
// vs 780-790).
func TestDiffMergeWalk_IndvStricterRefSkip(t *testing.T) {
	// 100: REF "." vs "A", s1 discordant (0/1 vs 0/0). 200: normal B match,
	// s1 concordant. .diff.sites counts BOTH; .diff.indv counts only 200.
	f1 := mergeHdr + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n" +
		"chr1\t100\t.\t.\tG\t30\tPASS\t.\tGT\t0/1\n" +
		"chr1\t200\t.\tC\tT\t30\tPASS\t.\tGT\t0/0\n"
	f2 := mergeHdr + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n" +
		"chr1\t100\t.\tA\tG\t30\tPASS\t.\tGT\t0/0\n" +
		"chr1\t200\t.\tC\tT\t30\tPASS\t.\tGT\t0/0\n"

	_, sites, indv := runDiffBodies(t, f1, f2)

	// .diff.sites: 100 is a B row (REF "." normalised to A) with real counts
	// (s1 discordant → 1/1); 200 concordant (0/0 vs 0/0).
	wantSites := "chr1\t100\tB\t1\t1\t1\t1\n" +
		"chr1\t200\tB\t1\t1\t0\t0\n"
	if sites != wantSites {
		t.Errorf(".diff.sites (indv-skip divergence) mismatch:\n--- got ---\n%s\n--- want ---\n%s", sites, wantSites)
	}

	// .diff.indv: only site 200 contributes (site 100 skipped by the stricter
	// by_indv REF check) → s1 = 1 common, 0 discord.
	if !strings.Contains(indv, "s1\t1\t0\t0") {
		t.Errorf(".diff.indv should exclude the '.'-vs-REF site, got:\n%s", indv)
	}
}

// TestDiffMergeWalk_File2OnlyTail pins that file-2-only sites AFTER the last
// file-1 position are still emitted (the close() tail flush), in order.
func TestDiffMergeWalk_File2OnlyTail(t *testing.T) {
	f1 := mergeHdr + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n" +
		"chr1\t100\t.\tA\tG\t30\tPASS\t.\tGT\t0/1\n"
	f2 := mergeHdr + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n" +
		"chr1\t100\t.\tA\tG\t30\tPASS\t.\tGT\t0/1\n" +
		"chr1\t200\t.\tC\tT\t30\tPASS\t.\tGT\t0/1\n" +
		"chr1\t300\t.\tG\tA\t30\tPASS\t.\tGT\t0/1\n"

	_, sites, _ := runDiffBodies(t, f1, f2)
	wantSites := "chr1\t100\tB\t1\t1\t0\t0\n" +
		"chr1\t200\t2\t0\t0\t0\t-nan\n" +
		"chr1\t300\t2\t0\t0\t0\t-nan\n"
	if sites != wantSites {
		t.Errorf(".diff.sites tail mismatch:\n--- got ---\n%s\n--- want ---\n%s", sites, wantSites)
	}
}
