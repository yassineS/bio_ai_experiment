package main

// Live-upstream parity tests for the bedtools-gaps closure wave: BAM/VCF/GFF
// input support and the -c typed-writer column-drop fix. They build the real
// upstream `bedtools` binary from the vendored submodule and assert this
// port's output matches byte-for-byte. They t.Fatalf (never t.Skip), matching
// the project's parity-rig policy.

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"sync"
	"testing"
)

// upstreamBedtoolsGaps builds (once) and returns the path to the upstream
// `bedtools` binary. It is uniquely named so this suite's builder can be
// reasoned about independently of the sibling behavioral/compat builders.
var (
	upstreamBedtoolsGapsOnce sync.Once
	upstreamBedtoolsGapsPath string
	upstreamBedtoolsGapsErr  error
)

func upstreamBedtoolsGaps(t *testing.T) string {
	t.Helper()
	upstreamBedtoolsGapsOnce.Do(func() {
		root, err := repoRoot()
		if err != nil {
			upstreamBedtoolsGapsErr = err
			return
		}
		dir := filepath.Join(root, "reference_code", "bedtools")
		if _, statErr := os.Stat(filepath.Join(dir, "Makefile")); statErr != nil {
			upstreamBedtoolsGapsErr = statErr
			return
		}
		bin := filepath.Join(dir, "bin", "bedtools")
		if _, statErr := os.Stat(bin); statErr != nil {
			cmd := exec.Command("make", "-j", "4")
			cmd.Dir = dir
			if out, buildErr := cmd.CombinedOutput(); buildErr != nil {
				upstreamBedtoolsGapsErr = &buildError{buildErr, out}
				return
			}
		}
		upstreamBedtoolsGapsPath = bin
	})
	if upstreamBedtoolsGapsErr != nil {
		t.Fatalf("upstream bedtools unavailable: %v\n"+
			"run: git submodule update --init reference_code/bedtools && "+
			"(cd reference_code/bedtools && make -j\"$(nproc)\")", upstreamBedtoolsGapsErr)
	}
	return upstreamBedtoolsGapsPath
}

// TestGapsParity_VCFGFFInput asserts byte-for-byte parity for VCF and GFF
// inputs (on -a and -b) across the default, -wa, -wb, -c, -v, -wo, -wao and
// -loj modes. VCF records echo their full original line; GFF and BED records
// clip their coordinates in the default intersection output, exactly as
// upstream prints them.
func TestGapsParity_VCFGFFInput(t *testing.T) {
	bt := upstreamBedtoolsGaps(t)
	ours := buildOurs(t)

	const (
		vcfA = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
			"chr1\t10\t.\tACG\t.\t.\t.\t.\nchr1\t50\t.\tT\t.\t.\t.\t.\n" +
			"chr1\t120\t.\tATGCATGCAT\t.\t.\t.\t.\n"
		gffA = "chr1\ttest\tgene\t1\t100\t.\t+\t.\tID=big\n" +
			"chr1\ttest\tgene\t150\t160\t.\t-\t.\tID=small\n"
		bed = "chr1\t5\t12\tx\t0\t+\nchr1\t48\t60\ty\t0\t-\nchr1\t20\t40\tz\t0\t+\n"
	)

	modes := [][]string{
		{}, {"-wa"}, {"-wb"}, {"-c"}, {"-v"},
		{"-wo"}, {"-wao"}, {"-loj"},
	}

	dir := t.TempDir()
	vcfPath := writeFile(t, dir, "a.vcf", vcfA)
	gffPath := writeFile(t, dir, "a.gff", gffA)
	bedPath := writeFile(t, dir, "b.bed", bed)

	pairs := []struct {
		name, a, b string
	}{
		{"vcfA_bedB", vcfPath, bedPath},
		{"gffA_bedB", gffPath, bedPath},
		{"bedA_vcfB", bedPath, vcfPath},
		{"bedA_gffB", bedPath, gffPath},
	}

	for _, p := range pairs {
		for _, m := range modes {
			name := p.name + "_" + joinArgs(m)
			t.Run(name, func(t *testing.T) {
				up := append([]string{"intersect", "-a", p.a, "-b", p.b}, m...)
				want := runCapture(t, bt, up...)
				our := append([]string{"-a", p.a, "-b", p.b}, m...)
				got := runCapture(t, ours, our...)
				if !bytes.Equal(got, want) {
					t.Fatalf("mismatch %s\nupstream:\n%s\nours:\n%s", name, want, got)
				}
			})
		}
	}

	// GFF carries strand in column 7, so -s is meaningful and must match
	// upstream (VCF has no strand column, where upstream errors on -s).
	t.Run("gff_strand", func(t *testing.T) {
		sb := writeFile(t, dir, "sb.bed", "chr1\t5\t12\tx\t0\t+\nchr1\t150\t160\ty\t0\t-\n")
		for _, m := range [][]string{{"-s"}, {"-s", "-c"}} {
			up := append([]string{"intersect", "-a", gffPath, "-b", sb}, m...)
			want := runCapture(t, bt, up...)
			our := append([]string{"-a", gffPath, "-b", sb}, m...)
			got := runCapture(t, ours, our...)
			if !bytes.Equal(got, want) {
				t.Fatalf("gff_strand mismatch %v\nupstream:\n%s\nours:\n%s", m, want, got)
			}
		}
	})
}

// TestGapsParity_BAMInput asserts byte-for-byte parity for BAM input via the
// upstream test fixtures, converting each alignment into a BED12 line (the
// `-bed` output). It covers the single-block and CIGAR-N multi-block cases,
// across default / -wa / -wb / -c / -v / -wo / -loj and the -split variants.
func TestGapsParity_BAMInput(t *testing.T) {
	bt := upstreamBedtoolsGaps(t)
	ours := buildOurs(t)
	root, err := repoRoot()
	if err != nil {
		t.Fatalf("locate repo root: %v", err)
	}
	bamDir := filepath.Join(root, "reference_code", "bedtools", "test", "intersect")
	aBam := filepath.Join(bamDir, "a.bam")
	threeBam := filepath.Join(bamDir, "three_blocks_match.bam")
	if _, err := os.Stat(aBam); err != nil {
		t.Fatalf("BAM fixture missing: %v (run git submodule update --init reference_code/bedtools)", err)
	}

	dir := t.TempDir()
	b := writeFile(t, dir, "b.bed", "chr1\t90\t120\tb1\nchr1\t150\t250\tb2\nchr1\t0\t50\tb3\n")

	cases := []struct {
		name string
		bam  string
		args []string
	}{
		{"a_bed", aBam, []string{"-bed"}},
		{"a_bed_wa", aBam, []string{"-bed", "-wa"}},
		{"a_bed_wb", aBam, []string{"-bed", "-wb"}},
		{"a_bed_c", aBam, []string{"-bed", "-c"}},
		{"a_bed_v", aBam, []string{"-bed", "-v"}},
		{"a_bed_wo", aBam, []string{"-bed", "-wo"}},
		{"a_bed_loj", aBam, []string{"-bed", "-loj"}},
		{"three_bed", threeBam, []string{"-bed"}},
		{"three_bed_split", threeBam, []string{"-bed", "-split"}},
		{"three_bed_wo", threeBam, []string{"-bed", "-wo"}},
		{"three_bed_split_wo", threeBam, []string{"-bed", "-split", "-wo"}},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			up := append([]string{"intersect", "-a", tc.bam, "-b", b}, tc.args...)
			want := runCapture(t, bt, up...)
			our := append([]string{"-a", tc.bam, "-b", b}, tc.args...)
			got := runCapture(t, ours, our...)
			if !bytes.Equal(got, want) {
				t.Fatalf("mismatch %s\nupstream:\n%s\nours:\n%s", tc.name, want, got)
			}
		})
	}
}

// TestGapsParity_CountColumnDrop pins the -c column-drop fix: a record whose
// score is 0 (or which has a strand column) must echo every original column
// before the trailing count, matching upstream byte-for-byte. The previous
// typed-writer path dropped the score and strand columns whenever the score
// was 0.
func TestGapsParity_CountColumnDrop(t *testing.T) {
	bt := upstreamBedtoolsGaps(t)
	ours := buildOurs(t)
	dir := t.TempDir()

	cases := []struct {
		name, a, b string
	}{
		{"score0_strand", "chr1\t10\t20\tfeat\t0\t+\nchr1\t100\t200\tg2\t0\t-\n", "chr1\t12\t18\thit\nchr1\t150\t160\th2\n"},
		{"bed3", "chr1\t10\t20\n", "chr1\t12\t18\n"},
		{"bedplus", "chr1\t10\t20\tn\t0\t+\textra1\textra2\n", "chr1\t12\t18\n"},
		{"no_overlap_zero", "chr1\t10\t20\tn\t0\t+\n", "chr2\t12\t18\n"},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			a := writeFile(t, dir, tc.name+"_a.bed", tc.a)
			b := writeFile(t, dir, tc.name+"_b.bed", tc.b)
			want := runCapture(t, bt, "intersect", "-a", a, "-b", b, "-c")
			got := runCapture(t, ours, "-a", a, "-b", b, "-c")
			if !bytes.Equal(got, want) {
				t.Fatalf("mismatch %s\nupstream:\n%s\nours:\n%s", tc.name, want, got)
			}
		})
	}
}

// joinArgs renders a flag slice into a filesystem-safe subtest suffix.
func joinArgs(args []string) string {
	if len(args) == 0 {
		return "default"
	}
	out := ""
	for i, a := range args {
		if i > 0 {
			out += "_"
		}
		for _, r := range a {
			if r == '-' {
				continue
			}
			out += string(r)
		}
	}
	return out
}
