package vcftools

import (
	"bytes"
	"path/filepath"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/vcf"
)

// TestParity_IndvBurden verifies --indv-burden byte-for-byte against
// upstream. The burden_fixture.vcf has 7 biallelic diploid sites and 4
// samples with a known mix of homozygous, heterozygous, and missing
// genotypes covering every branch of
// variant_file_output.cpp:471-483.
func TestParity_IndvBurden(t *testing.T) {
	prefix := runVcftoolsParity(t, "burden_fixture.vcf", &Params{
		IndvBurden: true,
		MinAlleles: 2,
	})
	got := readFileBytes(t, prefix+".iburden")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "burden.expected.iburden"))
	if !bytes.Equal(got, want) {
		t.Errorf(".iburden mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_IndvBurden_Derived verifies --indv-burden --derived. Two
// of the seven fixture sites carry INFO/AA (positions 600 and 700);
// the other five are dropped because AA is absent. Columns rename to
// N_HOM_ANC / N_HET / N_HOM_DER per upstream
// variant_file_output.cpp:402-403.
func TestParity_IndvBurden_Derived(t *testing.T) {
	prefix := runVcftoolsParity(t, "burden_fixture.vcf", &Params{
		IndvBurden: true,
		Derived:    true,
		MinAlleles: 2,
	})
	got := readFileBytes(t, prefix+".iburden")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "burden_derived.expected.iburden"))
	if !bytes.Equal(got, want) {
		t.Errorf(".iburden (derived) mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_IndvFreqBurden verifies --indv-freq-burden. The matrix is
// 4 (kept samples) × 9 (2*N+1) and exercises the
// double_count_hom_alt == 0 path where a hom-alt genotype contributes 2
// per-allele increments (upstream variant_file_output.cpp:603-609).
func TestParity_IndvFreqBurden(t *testing.T) {
	prefix := runVcftoolsParity(t, "burden_fixture.vcf", &Params{
		IndvFreqBurden: true,
		MinAlleles:     2,
	})
	got := readFileBytes(t, prefix+".ifreqburden")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "burden_freq.expected.ifreqburden"))
	if !bytes.Equal(got, want) {
		t.Errorf(".ifreqburden mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_IndvFreqBurden2 verifies --indv-freq-burden2 — same as
// --indv-freq-burden but with double_count_hom_alt=1, so a hom-alt
// genotype contributes 1 (not 2) to the corresponding allele-count bin.
// Mirrors upstream vcftools.cpp:64 + variant_file_output.cpp:605.
func TestParity_IndvFreqBurden2(t *testing.T) {
	prefix := runVcftoolsParity(t, "burden_fixture.vcf", &Params{
		IndvFreqBurden2: true,
		MinAlleles:      2,
	})
	got := readFileBytes(t, prefix+".ifreqburden")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "burden_freq2.expected.ifreqburden"))
	if !bytes.Equal(got, want) {
		t.Errorf(".ifreqburden (burden2) mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestParity_IndvFreqBurden_Derived verifies --indv-freq-burden
// --derived. Only the two AA-tagged sites contribute; the rest are
// dropped because AA cannot be resolved. The derived branch also picks
// the ancestral allele as the "skip" allele instead of REF (upstream
// variant_file_output.cpp:603, comparing against aa_idx).
func TestParity_IndvFreqBurden_Derived(t *testing.T) {
	prefix := runVcftoolsParity(t, "burden_fixture.vcf", &Params{
		IndvFreqBurden: true,
		Derived:        true,
		MinAlleles:     2,
	})
	got := readFileBytes(t, prefix+".ifreqburden")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "burden_freq_derived.expected.ifreqburden"))
	if !bytes.Equal(got, want) {
		t.Errorf(".ifreqburden (derived) mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestIndvFreqBurden_RemoveIndv_FixesLabelBug pins the FIX for the
// upstream label-index bug described in burden.go. Upstream
// variant_file_output.cpp:621 emits `meta_data.indv[indv_count]`
// instead of `meta_data.indv[ui]`, so removing S2 from
// [S1,S2,S3,S4] yields wrong labels `S1 S2 S3` next to S1/S3/S4
// burden values. Per CLAUDE.md ("don't replicate upstream bugs")
// we emit the CORRECT kept-sample labels. The golden fixture has
// been regenerated to match the fix.
func TestIndvFreqBurden_RemoveIndv_FixesLabelBug(t *testing.T) {
	prefix := runVcftoolsParity(t, "burden_fixture.vcf", &Params{
		IndvFreqBurden: true,
		RemoveIndvList: []string{"S2"},
		MinAlleles:     2,
	})
	got := readFileBytes(t, prefix+".ifreqburden")
	// Assert the first column is exactly the kept-sample IDs (S1, S3, S4)
	// — NOT the upstream-buggy (S1, S2, S3).
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "burden_freq_skip.expected.ifreqburden"))
	if !bytes.Equal(got, want) {
		t.Errorf(".ifreqburden mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestIndvBurden_SkipsNonDiploid pins that haploid (or partially-
// haploid) sites are skipped entirely by both burden runners, matching
// upstream's `if (e->is_diploid() == false) continue;` at
// variant_file_output.cpp:429-433 and :554-558. We use a small
// hand-built VCF mixing diploid and haploid rows; only the diploid row
// must contribute.
func TestIndvBurden_SkipsNonDiploid(t *testing.T) {
	vcfText := "##fileformat=VCFv4.0\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n" +
		"1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0/1\t1/1\n" + // diploid: contributes
		"1\t200\t.\tA\tC\t.\tPASS\t.\tGT\t1\t0\n" // haploid: skipped
	tmp := t.TempDir()
	in := bytes.NewBufferString(vcfText)
	prefix := filepath.Join(tmp, "out")
	if err := Run(in, &Params{OutPrefix: prefix, IndvBurden: true, MinAlleles: 2}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	got := string(readFileBytes(t, prefix+".iburden"))
	want := "INDV\tN_HOM_REF\tN_HET\tN_HOM_ALT\tN_MISS\n" +
		"A\t0\t1\t0\t0\n" +
		"B\t0\t0\t1\t0\n"
	if got != want {
		t.Errorf(".iburden mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestAncestralAlleleIndex covers the AA-resolution helper directly
// for the cases upstream calls out at variant_file_output.cpp:440-460:
// missing/`.`/`?` map to ok=false; case-insensitive REF/ALT matches
// return the right index; an AA that matches no allele yields
// ok=false.
func TestAncestralAlleleIndex(t *testing.T) {
	cases := []struct {
		name    string
		info    map[string]string
		ref     string
		alt     []string
		wantIdx int
		wantOk  bool
	}{
		{"missing", map[string]string{}, "A", []string{"C"}, 0, false},
		{"dot", map[string]string{"AA": "."}, "A", []string{"C"}, 0, false},
		{"question", map[string]string{"AA": "?"}, "A", []string{"C"}, 0, false},
		{"empty", map[string]string{"AA": ""}, "A", []string{"C"}, 0, false},
		{"ref-match-upper", map[string]string{"AA": "A"}, "A", []string{"C"}, 0, true},
		{"ref-match-lower", map[string]string{"AA": "a"}, "A", []string{"C"}, 0, true},
		{"alt1-match", map[string]string{"AA": "C"}, "A", []string{"C"}, 1, true},
		{"alt2-match", map[string]string{"AA": "G"}, "A", []string{"C", "G"}, 2, true},
		{"no-match", map[string]string{"AA": "N"}, "A", []string{"C"}, 0, false},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			v := &vcf.Variant{Ref: tc.ref, Alt: tc.alt, Info: tc.info}
			idx, ok := ancestralAlleleIndex(v)
			if ok != tc.wantOk || idx != tc.wantIdx {
				t.Errorf("ancestralAlleleIndex(...) = (%d,%v); want (%d,%v)", idx, ok, tc.wantIdx, tc.wantOk)
			}
		})
	}
}
