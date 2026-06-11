package main

import (
	"os"
	"path/filepath"
	"testing"

	"github.com/yassineS/bio_ai_experiment/tools/bcftools/pkg/bcftools"
)

// TestFilterMaskRequiresSoftFilter locks in upstream's rule that
// --mask / -M,--mask-file requires -s/--soft-filter (vcffilter.c:656):
// runFilter must reject a mask without a soft-filter name and accept it
// with one.
func TestFilterMaskRequiresSoftFilter(t *testing.T) {
	dir := t.TempDir()
	vcf := writeMainTempFile(t, dir, "f.vcf", "##fileformat=VCFv4.2\n"+
		"##contig=<ID=1>\n"+
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"+
		"1\t100\t.\tA\tT\t50\t.\t.\n")
	bed := writeMainTempFile(t, dir, "m.bed", "1\t50\t150\n")

	// Mask without -s must fail (rc 2).
	if rc := runFilter([]string{"-M", bed, "--no-version", vcf}); rc != 2 {
		t.Fatalf("mask without soft-filter: rc=%d, want 2", rc)
	}
	// Mask with -s must succeed (rc 0).
	out := filepath.Join(dir, "out.vcf")
	if rc := runFilter([]string{"-M", bed, "-s", "MASK", "-o", out, "--no-version", vcf}); rc != 0 {
		t.Fatalf("mask with soft-filter: rc=%d, want 0", rc)
	}
}

// writeMainTempFile writes content to dir/name and returns its path.
func writeMainTempFile(t *testing.T, dir, name, content string) string {
	t.Helper()
	p := filepath.Join(dir, name)
	if err := os.WriteFile(p, []byte(content), 0o644); err != nil {
		t.Fatalf("write %s: %v", p, err)
	}
	return p
}

// TestParseSnpGap accepts the upstream "INT[:TYPE]" form.
func TestParseSnpGap(t *testing.T) {
	cases := []struct {
		in   string
		want int
		err  bool
	}{
		{"", 0, false},
		{"5", 5, false},
		{"10:indel", 10, false},
		{"3:mnp,bnd,overlap", 3, false},
		{"abc", 0, true},
	}
	for _, tc := range cases {
		got, err := parseSnpGap(tc.in)
		if tc.err && err == nil {
			t.Errorf("parseSnpGap(%q) expected error, got %d", tc.in, got)
		}
		if !tc.err && err != nil {
			t.Errorf("parseSnpGap(%q) unexpected error: %v", tc.in, err)
		}
		if !tc.err && got != tc.want {
			t.Errorf("parseSnpGap(%q) = %d, want %d", tc.in, got, tc.want)
		}
	}
}

// TestParseHaplotypeNpIu checks the "NpIu" (phased-index / unphased-IUPAC)
// haplotype encoding parses into HapPhasedIUPAC with the right index, and
// that malformed variants are rejected.
func TestParseHaplotypeNpIu(t *testing.T) {
	cases := []struct {
		in      string
		wantSel bcftools.HaplotypeSelector
		wantIdx int
		wantErr bool
	}{
		{"2pIu", bcftools.HapPhasedIUPAC, 2, false},
		{"1pIu", bcftools.HapPhasedIUPAC, 1, false},
		{"10pIu", bcftools.HapPhasedIUPAC, 10, false},
		{"2PIU", bcftools.HapPhasedIUPAC, 2, false}, // case-insensitive suffix
		{"2", bcftools.HapIndex, 2, false},
		{"I", bcftools.HapIUPAC, 0, false},
		{"2pIuX", 0, 0, true},
		{"pIu", 0, 0, true},
		{"apIu", 0, 0, true},
	}
	for _, tc := range cases {
		sel, idx, err := bcftools.ParseHaplotypeSelector(tc.in)
		if tc.wantErr {
			if err == nil {
				t.Errorf("ParseHaplotypeSelector(%q): expected error, got sel=%v idx=%d", tc.in, sel, idx)
			}
			continue
		}
		if err != nil {
			t.Errorf("ParseHaplotypeSelector(%q): unexpected error %v", tc.in, err)
			continue
		}
		if sel != tc.wantSel || idx != tc.wantIdx {
			t.Errorf("ParseHaplotypeSelector(%q) = (%v,%d), want (%v,%d)", tc.in, sel, idx, tc.wantSel, tc.wantIdx)
		}
	}
}
