package vcftools

// Unit and parity tests for --mendel.
//
// Goldens generated from upstream vcftools (reference_code submodule) with the
// documented FORTIFY_SOURCE workaround. See diff_parity_test.go's header for
// the make incantation.

import (
	"os"
	"path/filepath"
	"testing"
)

// TestParity_Mendel exercises --mendel against upstream byte-for-byte. The
// fixture has a single trio (child / father / mother) and four sites: 100,
// 300 and 400 are Mendel errors; 200 is consistent.
func TestParity_Mendel(t *testing.T) {
	fix := vcftoolsFixtureDir(t)
	prefix := runVcftoolsParity(t, "mendel.vcf", &Params{
		MendelPedFile: filepath.Join(fix, "mendel.ped"),
	})
	got := readFileBytes(t, prefix+".mendel")
	want := readFileBytes(t, filepath.Join(fix, "mendel.expected.mendel"))
	if string(got) != string(want) {
		t.Errorf(".mendel byte mismatch:\n--- got ---\n%s\n--- want ---\n%s",
			string(got), string(want))
	}
}

// TestIsMendelConsistent walks the Mendelian truth table for biallelic
// genotypes. (m1,m2)×(f1,f2) yields four possible child genotypes; a child
// (c1,c2) is consistent iff (c1,c2) or (c2,c1) appears in that set.
func TestIsMendelConsistent(t *testing.T) {
	cases := []struct {
		name           string
		c1, c2, f1, f2 int
		m1, m2         int
		want           bool
	}{
		// (0/0) x (0/0) -> only 0/0 possible; 0/0 child consistent.
		{"hom00_hom00_hom00", 0, 0, 0, 0, 0, 0, true},
		// (0/0) x (0/0) -> 0/1 child impossible.
		{"hom00_hom00_het01", 0, 1, 0, 0, 0, 0, false},
		// (0/1) x (0/0) -> 0/0 or 0/1 possible.
		{"het01_hom00_het01", 0, 1, 0, 0, 0, 1, true},
		{"het01_hom00_hom11", 1, 1, 0, 0, 0, 1, false},
		// (1/1) x (0/0) -> only 0/1 possible.
		{"hom11_hom00_het01", 0, 1, 0, 0, 1, 1, true},
		{"hom11_hom00_hom11", 1, 1, 0, 0, 1, 1, false},
		// (0/1) x (0/1) -> 0/0, 0/1, 1/1 all possible.
		{"het01_het01_hom00", 0, 0, 0, 1, 0, 1, true},
		{"het01_het01_het01", 0, 1, 0, 1, 0, 1, true},
		{"het01_het01_hom11", 1, 1, 0, 1, 0, 1, true},
		// Swapped child ordering must still match.
		{"swapped_het10", 1, 0, 0, 0, 0, 1, true},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			got := isMendelConsistent(c.c1, c.c2, c.f1, c.f2, c.m1, c.m2)
			if got != c.want {
				t.Errorf("isMendelConsistent(c=%d/%d,f=%d/%d,m=%d/%d) = %v; want %v",
					c.c1, c.c2, c.f1, c.f2, c.m1, c.m2, got, c.want)
			}
		})
	}
}

// TestParseDiploidGT covers the diploid-only parser used by mendel.go.
// Haploid and any-slot-missing calls must return missing=true.
func TestParseDiploidGT(t *testing.T) {
	cases := []struct {
		gt      string
		a, b    int
		missing bool
	}{
		{"0/0", 0, 0, false},
		{"0/1", 0, 1, false},
		{"1|2", 1, 2, false},
		{"./.", 0, 0, true},
		{"0/.", 0, 0, true},
		{"./0", 0, 0, true},
		{"", 0, 0, true},
		{".", 0, 0, true},
		{"0", 0, 0, true},   // haploid → missing for Mendel
		{"abc", 0, 0, true}, // junk
	}
	for _, c := range cases {
		a, b, m := parseDiploidGT(c.gt)
		if a != c.a || b != c.b || m != c.missing {
			t.Errorf("parseDiploidGT(%q) = (%d,%d,%v); want (%d,%d,%v)",
				c.gt, a, b, m, c.a, c.b, c.missing)
		}
	}
}

// TestLoadMendelPED checks that:
//   - the first line is always skipped
//   - comments and blanks are skipped
//   - rows with "0" placeholders are dropped
//   - rows referencing samples not in the VCF are dropped
func TestLoadMendelPED(t *testing.T) {
	dir := t.TempDir()
	pedPath := filepath.Join(dir, "test.ped")
	const ped = `FAM IID FAT MOT
fam1 child father mother
# comment
fam2 c2 f2 m2
fam3 c3 0 m3
fam4 c4 f4 0
`
	if err := os.WriteFile(pedPath, []byte(ped), 0o600); err != nil {
		t.Fatalf("WriteFile: %v", err)
	}
	// VCF samples — only fam1's trio fully resolves.
	samples := []string{"child", "father", "mother", "c2"}
	trios, err := loadMendelPED(pedPath, samples)
	if err != nil {
		t.Fatalf("loadMendelPED: %v", err)
	}
	if len(trios) != 1 {
		t.Fatalf("got %d trios, want 1: %+v", len(trios), trios)
	}
	if trios[0].family != "child_father_mother" {
		t.Errorf("family = %q; want %q", trios[0].family, "child_father_mother")
	}
}
