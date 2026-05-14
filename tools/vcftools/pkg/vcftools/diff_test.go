package vcftools

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestCanonicalBiallelicGT(t *testing.T) {
	cases := []struct {
		gt   string
		a, b int
		miss bool
	}{
		{"0/0", 0, 0, false},
		{"0/1", 0, 1, false},
		{"1/0", 0, 1, false},
		{"1|1", 1, 1, false},
		{"./.", 0, 0, true},
		{".|.", 0, 0, true},
		{".", 0, 0, true},
		{"", 0, 0, true},
		{"0", 0, 0, true},   // haploid: not supported
		{"2/2", 0, 0, true}, // higher alt: treated as missing
		{"0/2", 0, 0, true},
		{"0/.", 0, 0, true},
	}
	for _, c := range cases {
		a, b, miss := canonicalBiallelicGT(c.gt)
		if miss != c.miss || a != c.a || b != c.b {
			t.Errorf("canonicalBiallelicGT(%q) = (%d,%d,%v), want (%d,%d,%v)",
				c.gt, a, b, miss, c.a, c.b, c.miss)
		}
	}
}

const diffVCF1 = `##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1	s2	s3
chr1	100	.	A	G	30	PASS	.	GT	0/0	0/1	1/1
chr1	200	.	C	T	30	PASS	.	GT	0/1	0/0	./.
chr1	400	.	G	A	30	PASS	.	GT	0/0	0/1	1/1
`

const diffVCF2 = `##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1	s2	s4
chr1	100	.	A	G	30	PASS	.	GT	0/0	1/1	1/1
chr1	200	.	C	T	30	PASS	.	GT	0/1	0/0	0/1
chr1	300	.	T	C	30	PASS	.	GT	0/1	1/1	0/0
`

func TestRunDiffFamily(t *testing.T) {
	tmp := t.TempDir()
	diff2 := filepath.Join(tmp, "f2.vcf")
	if err := os.WriteFile(diff2, []byte(diffVCF2), 0o644); err != nil {
		t.Fatalf("write f2: %v", err)
	}
	prefix := filepath.Join(tmp, "cmp")
	err := Run(strings.NewReader(diffVCF1), &Params{
		OutPrefix:           prefix,
		Diff:                diff2,
		DiffSite:            true,
		DiffIndv:            true,
		DiffSiteDiscordance: true,
		DiffIndvDiscordance: true,
	})
	if err != nil {
		t.Fatalf("Run --diff: %v", err)
	}

	mustRead := func(path string) string {
		body, err := os.ReadFile(path)
		if err != nil {
			t.Fatalf("read %s: %v", path, err)
		}
		return string(body)
	}

	// .diff.sites_in_files - chr1:100 and chr1:200 are in both; chr1:400 in 1
	// only; chr1:300 in 2 only.
	sif := mustRead(prefix + ".diff.sites_in_files")
	if !strings.Contains(sif, "CHROM\tPOS1\tPOS2\tIN_FILE\tREF1\tREF2\tALT1\tALT2") {
		t.Errorf("missing header: %s", sif)
	}
	expected := []string{
		"chr1\t100\t100\tB\tA\tA\tG\tG",
		"chr1\t200\t200\tB\tC\tC\tT\tT",
		"chr1\t400\t.\t1\tG\t.\tA\t.",
		"chr1\t.\t300\t2\t.\tT\t.\tC",
	}
	for _, line := range expected {
		if !strings.Contains(sif, line) {
			t.Errorf("sites_in_files missing line %q in:\n%s", line, sif)
		}
	}

	// .diff.indv_in_files - s1/s2 are in both, s3 in file1, s4 in file2.
	indv := mustRead(prefix + ".diff.indv_in_files")
	if !strings.Contains(indv, "INDV\tFILES") {
		t.Errorf("indv_in_files missing header: %s", indv)
	}
	for _, line := range []string{"s1\tB", "s2\tB", "s3\t1", "s4\t2"} {
		if !strings.Contains(indv, line) {
			t.Errorf("indv_in_files missing %q in:\n%s", line, indv)
		}
	}

	// .diff.sites: only chr1:100 and chr1:200 (the shared sites)
	//
	// chr1:100  s1=0/0 / 0/0 (concordant); s2=0/1 / 1/1 (discordant). s3
	// only in file1 → ignored. Common samples shared between files = {s1,s2}.
	// N_COMMON_CALLED = 2, N_DISCORD = 1.
	//
	// chr1:200 s1=0/1 / 0/1 (concordant); s2=0/0 / 0/0 (concordant).
	// N_COMMON_CALLED = 2, N_DISCORD = 0.
	sites := mustRead(prefix + ".diff.sites")
	for _, line := range []string{
		"chr1\t100\t2\t1",
		"chr1\t200\t2\t0",
	} {
		if !strings.Contains(sites, line) {
			t.Errorf("diff.sites missing %q in:\n%s", line, sites)
		}
	}
	if strings.Contains(sites, "chr1\t400") || strings.Contains(sites, "chr1\t300") {
		t.Errorf("diff.sites should only contain shared sites:\n%s", sites)
	}

	// .diff.indv: per-individual discordance over common sites.
	//   s1: chr1:100 concordant, chr1:200 concordant → 2 called, 0 discord
	//   s2: chr1:100 discordant, chr1:200 concordant → 2 called, 1 discord
	indvD := mustRead(prefix + ".diff.indv")
	for _, line := range []string{
		"s1\t2\t0",
		"s2\t2\t1",
	} {
		if !strings.Contains(indvD, line) {
			t.Errorf("diff.indv missing %q in:\n%s", line, indvD)
		}
	}
}

func TestLoadDiffVCFMissing(t *testing.T) {
	if _, err := loadDiffVCF(filepath.Join(t.TempDir(), "nope.vcf")); err == nil {
		t.Errorf("expected error for missing diff VCF")
	}
}

func TestNewDiffRunnerNoFlag(t *testing.T) {
	r, err := newDiffRunner(&Params{}, nil)
	if err != nil {
		t.Fatalf("err: %v", err)
	}
	if r != nil {
		t.Errorf("expected nil runner when --diff is empty")
	}
}

func TestNewDiffRunnerNoOutputFlags(t *testing.T) {
	// Diff set but no diff-* sub-flag → no-op (don't even open the file).
	r, err := newDiffRunner(&Params{Diff: "should-be-ignored"}, nil)
	if err != nil {
		t.Fatalf("err: %v", err)
	}
	if r != nil {
		t.Errorf("expected nil runner when no diff-* sub-flag")
	}
}

func TestDiffRunnerNilSafe(t *testing.T) {
	var r *diffRunner
	if err := r.addVariant(nil); err != nil {
		t.Errorf("nil addVariant: %v", err)
	}
	if err := r.close(); err != nil {
		t.Errorf("nil close: %v", err)
	}
}

func TestEnvDiffWarnNoop(t *testing.T) {
	// Just exercise the helper for coverage; it's a stderr writer.
	envDiffWarn(nil)
	envDiffWarn([]string{"a", "b"})
}
