package vcftools

import (
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func writeTempFile(t *testing.T, name, contents string) string {
	t.Helper()
	dir := t.TempDir()
	path := filepath.Join(dir, name)
	if err := os.WriteFile(path, []byte(contents), 0o644); err != nil {
		t.Fatalf("writing %s: %v", path, err)
	}
	return path
}

func TestMergeIntervals(t *testing.T) {
	cases := []struct {
		name string
		in   []bedInterval
		want []bedInterval
	}{
		{"empty", nil, nil},
		{
			"sorted disjoint",
			[]bedInterval{{0, 10}, {20, 30}},
			[]bedInterval{{0, 10}, {20, 30}},
		},
		{
			"out-of-order disjoint",
			[]bedInterval{{20, 30}, {0, 10}},
			[]bedInterval{{0, 10}, {20, 30}},
		},
		{
			"overlap",
			[]bedInterval{{0, 10}, {5, 15}},
			[]bedInterval{{0, 15}},
		},
		{
			"adjacent (touching)",
			[]bedInterval{{0, 10}, {10, 20}},
			[]bedInterval{{0, 20}},
		},
		{
			"nested",
			[]bedInterval{{0, 100}, {10, 20}, {30, 40}},
			[]bedInterval{{0, 100}},
		},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			got := mergeIntervals(append([]bedInterval(nil), c.in...))
			if len(got) != len(c.want) {
				t.Fatalf("len = %d, want %d (got %v)", len(got), len(c.want), got)
			}
			for i := range got {
				if got[i] != c.want[i] {
					t.Errorf("[%d] got %v, want %v", i, got[i], c.want[i])
				}
			}
		})
	}
}

func TestLoadBedRegionsAndContains(t *testing.T) {
	contents := strings.Join([]string{
		"# a comment line",
		"chr1\t100\t200",
		"chr1\t150\t250", // overlaps with previous → merge
		"chr1\t1000\t2000",
		"chr2\t0\t10",
		"track name=foo", // BED parser tolerates these
		"chr3\t5\t5",     // zero-length, should be skipped
		"",               // blank line
	}, "\n") + "\n"
	path := writeTempFile(t, "regions.bed", contents)

	r, err := loadBedRegions(path)
	if err != nil {
		t.Fatalf("loadBedRegions: %v", err)
	}

	// chr1 merged to [100,250) and [1000,2000); chr2 to [0,10).
	if got := len(r.byChrom["chr1"]); got != 2 {
		t.Errorf("chr1 intervals = %d, want 2; have %v", got, r.byChrom["chr1"])
	}
	if got := len(r.byChrom["chr3"]); got != 0 {
		t.Errorf("chr3 should be empty (zero-length skipped); got %v", r.byChrom["chr3"])
	}

	cases := []struct {
		chrom string
		pos   int
		want  bool
	}{
		// VCF POS is 1-based; BED is 0-based half-open. [100,200) maps to POS 101..200.
		{"chr1", 100, false},
		{"chr1", 101, true},
		{"chr1", 200, true},
		{"chr1", 250, true}, // upper edge after merge
		{"chr1", 251, false},
		{"chr1", 1000, false}, // gap
		{"chr1", 1001, true},
		{"chr1", 2000, true},
		{"chr1", 2001, false},
		{"chr2", 1, true},
		{"chr2", 10, true},
		{"chr2", 11, false},
		{"chrUnknown", 100, false},
	}
	for _, c := range cases {
		got := r.containsVCFPos(c.chrom, c.pos)
		if got != c.want {
			t.Errorf("containsVCFPos(%s, %d) = %v, want %v", c.chrom, c.pos, got, c.want)
		}
	}
}

func TestBedRegionsNilReceiver(t *testing.T) {
	var r *bedRegions
	if r.containsVCFPos("chr1", 100) {
		t.Errorf("nil receiver must always return false")
	}
}

func TestLoadBedRegionsMissingFile(t *testing.T) {
	if _, err := loadBedRegions(filepath.Join(t.TempDir(), "does-not-exist.bed")); err == nil {
		t.Errorf("expected error opening missing file")
	}
}

func TestRunBedFilters(t *testing.T) {
	const vcf = `##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1
chr1	100	.	A	G	30	PASS	.	GT	0/0
chr1	200	.	A	G	30	PASS	.	GT	0/0
chr1	300	.	A	G	30	PASS	.	GT	0/0
chr2	50	.	A	G	30	PASS	.	GT	0/0
`
	// [chr1 99 200) covers POS 100 and 200. [chr2 0 100) covers POS 50.
	bedContent := "chr1\t99\t200\nchr2\t0\t100\n"
	bedPath := writeTempFile(t, "include.bed", bedContent)
	tmp := t.TempDir()
	prefix := filepath.Join(tmp, "out")

	err := Run(strings.NewReader(vcf), &Params{
		OutPrefix: prefix,
		Bed:       bedPath,
		Recode:    true,
	})
	if err != nil {
		t.Fatalf("Run --bed: %v", err)
	}
	got, err := os.ReadFile(prefix + ".recode.vcf")
	if err != nil {
		t.Fatalf("read recode: %v", err)
	}
	body := string(got)
	if !strings.Contains(body, "chr1\t100\t") || !strings.Contains(body, "chr1\t200\t") || !strings.Contains(body, "chr2\t50\t") {
		t.Errorf("--bed missing expected sites: %s", body)
	}
	if strings.Contains(body, "chr1\t300\t") {
		t.Errorf("--bed should have filtered chr1 300: %s", body)
	}

	// --exclude-bed inverts.
	prefix2 := filepath.Join(tmp, "out2")
	err = Run(strings.NewReader(vcf), &Params{
		OutPrefix:  prefix2,
		ExcludeBed: bedPath,
		Recode:     true,
	})
	if err != nil {
		t.Fatalf("Run --exclude-bed: %v", err)
	}
	got2, err := os.ReadFile(prefix2 + ".recode.vcf")
	if err != nil {
		t.Fatalf("read recode 2: %v", err)
	}
	body2 := string(got2)
	if !strings.Contains(body2, "chr1\t300\t") {
		t.Errorf("--exclude-bed should keep chr1 300: %s", body2)
	}
	if strings.Contains(body2, "chr1\t100\t") || strings.Contains(body2, "chr2\t50\t") {
		t.Errorf("--exclude-bed should drop covered sites: %s", body2)
	}
}
