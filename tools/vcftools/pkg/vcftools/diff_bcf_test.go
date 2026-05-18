package vcftools

// Tests for `--diff-bcf` (BCF input for the --diff-* family). PR-G
// chunk 3, completing vcftools BCF I/O after wave 21 (--recode-bcf)
// and wave 22 (--bcf, --contigs).
//
// The flag accepts BGZF-compressed BCF v2.2 as the second comparison
// file; downstream --diff-* outputs are unchanged in shape (and
// pinned by the existing parity tests in diff_parity_test.go).

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// TestRun_DiffBCF_Roundtrip writes the second comparison file as a
// BCF via --recode-bcf, then runs the --diff-* family against the
// same source VCF on the file-1 side: every common site should be
// classified as "B" with matching alleles and zero discordance.
func TestRun_DiffBCF_Roundtrip(t *testing.T) {
	tmp := t.TempDir()

	// Step 1: write the file-2 BCF.
	bcfPrefix := filepath.Join(tmp, "f2")
	if err := Run(strings.NewReader(multiFormatBCFFixture), &Params{
		OutPrefix:     bcfPrefix,
		RecodeBCF:     true,
		RecodeInfoAll: true,
	}); err != nil {
		t.Fatalf("Run --recode-bcf: %v", err)
	}

	// Step 2: run --diff-bcf against the same source VCF on file-1.
	outPrefix := filepath.Join(tmp, "diff")
	if err := Run(strings.NewReader(multiFormatBCFFixture), &Params{
		OutPrefix:           outPrefix,
		DiffBCF:             bcfPrefix + ".recode.bcf",
		DiffSiteDiscordance: true,
	}); err != nil {
		t.Fatalf("Run --diff-bcf --diff-site-discordance: %v", err)
	}
	body, err := os.ReadFile(outPrefix + ".diff.sites")
	if err != nil {
		t.Fatal(err)
	}
	out := string(body)
	// Every site is in both files, alleles match, GTs are identical
	// (same source) → MATCHING_ALLELES=1, N_DISCORD=0, DISCORDANCE=0.
	for _, want := range []string{
		"chr1\t100\tB\t1\t2\t0\t0",
		"chr1\t200\tB\t1\t2\t0\t0",
	} {
		if !strings.Contains(out, want) {
			t.Errorf("output missing %q in:\n%s", want, out)
		}
	}
}

// TestRun_DiffBCF_ComposesWithMap exercises --diff-bcf together with
// --diff-indv-map, the most-used dual-file diff combination. The
// fixture has different sample names between the two files; the map
// renames file-2 samples onto file-1 sample names.
func TestRun_DiffBCF_ComposesWithMap(t *testing.T) {
	tmp := t.TempDir()

	// File-2 in BCF form with different sample names.
	const f2FixtureForBCF = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	a1	a2
chr1	100	.	A	G	30	PASS	DP=10	GT	0/0	0/1
chr1	200	.	C	T	30	PASS	DP=12	GT	0/1	1/1
`
	bcfPrefix := filepath.Join(tmp, "f2")
	if err := Run(strings.NewReader(f2FixtureForBCF), &Params{
		OutPrefix:     bcfPrefix,
		RecodeBCF:     true,
		RecodeInfoAll: true,
	}); err != nil {
		t.Fatalf("write f2 BCF: %v", err)
	}

	// File-1 VCF with names s1/s2.
	const f1Fixture = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1	s2
chr1	100	.	A	G	30	PASS	DP=15	GT	0/0	0/1
chr1	200	.	C	T	30	PASS	DP=18	GT	0/1	1/1
`
	mapPath := filepath.Join(tmp, "map.txt")
	if err := os.WriteFile(mapPath, []byte("a1\ts1\na2\ts2\n"), 0o644); err != nil {
		t.Fatal(err)
	}

	outPrefix := filepath.Join(tmp, "diff")
	if err := Run(strings.NewReader(f1Fixture), &Params{
		OutPrefix:           outPrefix,
		DiffBCF:             bcfPrefix + ".recode.bcf",
		DiffIndvMap:         mapPath,
		DiffSiteDiscordance: true,
		DiffIndvDiscordance: true,
	}); err != nil {
		t.Fatalf("Run --diff-bcf --diff-indv-map: %v", err)
	}
	sites, err := os.ReadFile(outPrefix + ".diff.sites")
	if err != nil {
		t.Fatal(err)
	}
	// Both sites have the same genotypes between file-1 and file-2
	// after the map renames a1→s1 and a2→s2, so DISCORDANCE = 0.
	if !strings.Contains(string(sites), "chr1\t100\tB\t1\t2\t0\t0") {
		t.Errorf(".diff.sites missing common+match line:\n%s", sites)
	}
	indv, err := os.ReadFile(outPrefix + ".diff.indv")
	if err != nil {
		t.Fatal(err)
	}
	if !strings.Contains(string(indv), "s1\t2\t0\t0") || !strings.Contains(string(indv), "s2\t2\t0\t0") {
		t.Errorf(".diff.indv missing zero-discord rows:\n%s", indv)
	}
}

// TestRun_DiffBCF_CompletesAllDiffOutputs verifies --diff-bcf with the
// full set of diff outputs to ensure none of them blow up on a BCF
// second file (regression coverage in case the wave-22 reader did
// something subtle that only showed up on certain outputs).
func TestRun_DiffBCF_CompletesAllDiffOutputs(t *testing.T) {
	tmp := t.TempDir()
	bcfPrefix := filepath.Join(tmp, "f2")
	if err := Run(strings.NewReader(multiFormatBCFFixture), &Params{
		OutPrefix:     bcfPrefix,
		RecodeBCF:     true,
		RecodeInfoAll: true,
	}); err != nil {
		t.Fatalf("write f2 BCF: %v", err)
	}
	outPrefix := filepath.Join(tmp, "diff")
	if err := Run(strings.NewReader(multiFormatBCFFixture), &Params{
		OutPrefix:             outPrefix,
		DiffBCF:               bcfPrefix + ".recode.bcf",
		DiffSite:              true,
		DiffIndv:              true,
		DiffSiteDiscordance:   true,
		DiffIndvDiscordance:   true,
		DiffDiscordanceMatrix: true,
	}); err != nil {
		t.Fatalf("Run --diff-bcf with all outputs: %v", err)
	}
	for _, suffix := range []string{
		".diff.sites_in_files",
		".diff.indv_in_files",
		".diff.sites",
		".diff.indv",
		".diff.discordance_matrix",
	} {
		path := outPrefix + suffix
		st, err := os.Stat(path)
		if err != nil {
			t.Errorf("missing %s: %v", suffix, err)
			continue
		}
		if st.Size() == 0 {
			t.Errorf("%s is empty", suffix)
		}
	}
}

// TestLoadDiffBCF_OpensAndReads is a focused unit test on the new
// `loadDiffBCF` loader: it should successfully open a BCF file, read
// every variant into the diffData map, and surface samples.
func TestLoadDiffBCF_OpensAndReads(t *testing.T) {
	tmp := t.TempDir()
	bcfPrefix := filepath.Join(tmp, "in")
	if err := Run(strings.NewReader(multiFormatBCFFixture), &Params{
		OutPrefix:     bcfPrefix,
		RecodeBCF:     true,
		RecodeInfoAll: true,
	}); err != nil {
		t.Fatalf("write BCF: %v", err)
	}
	d, err := loadDiffBCF(bcfPrefix + ".recode.bcf")
	if err != nil {
		t.Fatalf("loadDiffBCF: %v", err)
	}
	if got, want := d.samples, []string{"s1", "s2"}; !equalSS(got, want) {
		t.Errorf("samples: %v, want %v", got, want)
	}
	if _, ok := d.sites["chr1"][100]; !ok {
		t.Errorf("missing chr1:100 in sites map")
	}
	if _, ok := d.sites["chr1"][200]; !ok {
		t.Errorf("missing chr1:200 in sites map")
	}
}

// TestLoadDiffBCF_MissingFile pins the error path when the BCF file
// doesn't exist or isn't readable.
func TestLoadDiffBCF_MissingFile(t *testing.T) {
	if _, err := loadDiffBCF("/nope/missing.bcf"); err == nil {
		t.Fatal("expected error for missing BCF file")
	}
}
