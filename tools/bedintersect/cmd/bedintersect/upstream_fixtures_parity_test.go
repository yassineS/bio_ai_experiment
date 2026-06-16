package main

// Live-upstream parity tests over the vendored bedtools test fixtures.
//
// These drive the gap-closing flag surface (-wa/-wb/-wo/-wao/-u/-c/-C/-loj,
// -f/-F/-r/-e, -s/-S, -split, -sorted/-g, -v, multiple -b with
// -names/-filenames/-sortout, BAM/VCF/GFF inputs, -header, -nonamecheck)
// against the live upstream `bedtools intersect` binary, comparing byte for
// byte over the real fixtures under reference_code/bedtools/test/intersect.
//
// Tests use t.Fatalf (never t.Skip): a missing/un-buildable upstream binary or
// any divergence is a hard failure. The upstream binary and our binary are each
// built once via the shared upstreamBedtools / buildOurs helpers.

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"testing"
)

// fixturesDir returns the absolute path to the upstream intersect fixtures,
// failing hard (with the submodule-init hint) when they are absent.
func fixturesDir(t *testing.T) string {
	t.Helper()
	root, err := repoRoot()
	if err != nil {
		t.Fatalf("locate repo root: %v", err)
	}
	dir := filepath.Join(root, "reference_code", "bedtools", "test", "intersect")
	if _, err := os.Stat(filepath.Join(dir, "a.bed")); err != nil {
		t.Fatalf("upstream intersect fixtures missing (%v)\n"+
			"run: git submodule update --init reference_code/bedtools", err)
	}
	return dir
}

// runOut runs a command in dir and returns stdout and stderr separately, never
// failing on a non-zero exit (several parity cases intentionally exercise
// upstream's error/exit paths).
func runOut(t *testing.T, dir, name string, args ...string) (stdout, stderr []byte) {
	t.Helper()
	cmd := exec.Command(name, args...)
	cmd.Dir = dir
	var so, se bytes.Buffer
	cmd.Stdout = &so
	cmd.Stderr = &se
	_ = cmd.Run()
	return so.Bytes(), se.Bytes()
}

// fixtureCases are the (name, flag-args) pairs exercised against both binaries
// over a fixed -a/-b fixture pair. Each set is run for its own fixtures below.
type fixtureCase struct {
	name string
	a    string
	b    []string // one or more B files
	args []string
}

// TestUpstreamFixtures_Intersect runs a broad matrix of flag combinations over
// the real upstream fixtures and asserts byte-for-byte stdout parity with the
// live `bedtools intersect`. stderr is compared as a line set for the cases
// that emit diagnostics, avoiding stdout/stderr interleaving artifacts.
func TestUpstreamFixtures_Intersect(t *testing.T) {
	bt := upstreamBedtools(t)
	ours := buildOurs(t)
	dir := fixturesDir(t)

	cases := []fixtureCase{
		// Core output modes over the canonical strand fixtures.
		{"default", "a.bed", []string{"b.bed"}, nil},
		{"wa", "a.bed", []string{"b.bed"}, []string{"-wa"}},
		{"wb", "a.bed", []string{"b.bed"}, []string{"-wb"}},
		{"wa_wb", "a.bed", []string{"b.bed"}, []string{"-wa", "-wb"}},
		{"wo", "a.bed", []string{"b.bed"}, []string{"-wo"}},
		{"wao", "a.bed", []string{"b.bed"}, []string{"-wao"}},
		{"loj", "a.bed", []string{"b.bed"}, []string{"-loj"}},
		{"count", "a.bed", []string{"b.bed"}, []string{"-c"}},
		{"countEach", "a.bed", []string{"b.bed"}, []string{"-C"}},
		{"unique", "a.bed", []string{"b.bed"}, []string{"-u"}},
		{"invert", "a.bed", []string{"b.bed"}, []string{"-v"}},

		// Strand.
		{"strand", "a.bed", []string{"b.bed"}, []string{"-s"}},
		{"opp_strand", "a.bed", []string{"b.bed"}, []string{"-S"}},
		{"strand_wb", "a.bed", []string{"b.bed"}, []string{"-s", "-wb"}},
		{"opp_strand_wb", "a.bed", []string{"b.bed"}, []string{"-S", "-wb"}},
		{"unique_strand", "a.bed", []string{"b.bed"}, []string{"-u", "-s"}},
		{"invert_opp_strand", "a.bed", []string{"b.bed"}, []string{"-v", "-S"}},
		{"countEach_strand", "a.bed", []string{"b.bed"}, []string{"-C", "-s"}},

		// Fractions: -f / -F / -r / -e.
		{"f_half_wo", "a.bed", []string{"b.bed"}, []string{"-f", "0.5", "-wo"}},
		{"e_f_F_a", "a.bed", []string{"b.bed"}, []string{"-f", "0.1", "-F", "0.5", "-e", "-wo"}},
		{"e_f_F_b", "a.bed", []string{"b.bed"}, []string{"-f", "0.9", "-F", "0.1", "-e", "-wo"}},
		{"reciprocal", "a.bed", []string{"b.bed"}, []string{"-f", "0.5", "-r", "-wo"}},
		{"countEach_f", "a.bed", []string{"b.bed"}, []string{"-C", "-f", "0.5"}},

		// x/y fraction fixtures (designed to stress -f/-F/-r/-e boundaries).
		{"xy_f", "x.bed", []string{"y.bed"}, []string{"-f", "0.2"}},
		{"xy_F_wawb", "x.bed", []string{"y.bed"}, []string{"-F", "0.21", "-wa", "-wb"}},
		{"xy_f_r_wawb", "x.bed", []string{"y.bed"}, []string{"-f", "0.19", "-r", "-wa", "-wb"}},
		{"xy_f_F_wawb", "x.bed", []string{"y.bed"}, []string{"-f", "0.19", "-F", "0.21", "-wa", "-wb"}},
		{"xy_f_F_e_wawb", "x.bed", []string{"y.bed"}, []string{"-f", "0.21", "-F", "0.21", "-e", "-wa", "-wb"}},

		// Bin-order / nesting fixtures (default hit order must match upstream).
		{"bug150", "bug150_a.bed", []string{"bug150_b.bed"}, nil},
		{"bug187_wao", "bug187_a.bed", []string{"bug187_b.bed"}, []string{"-wao"}},
		{"bug167_s", "bug167_strandSweep.bed", []string{"bug167_strandSweep.bed"}, []string{"-s"}},

		// Split (BED12) blocks.
		{"split_wo", "blocks.bed12", []string{"blocks.bed12"}, []string{"-split", "-wo"}},
		{"split_wa_wb", "blocks.bed12", []string{"blocks.bed12"}, []string{"-split", "-wa", "-wb"}},
		{"split_f", "blocks.bed12", []string{"blocks.bed12"}, []string{"-split", "-f", "0.1", "-wo"}},
		{"split_default", "blocks.bed12", []string{"blocks.bed12"}, []string{"-split"}},
		{"split_loj", "blocks.bed12", []string{"blocks.bed12"}, []string{"-split", "-loj"}},

		// Zero-length records.
		{"zerolen_wo", "a_testZeroLen.bed", []string{"b_testZeroLen.bed"}, []string{"-wo"}},
		{"zerolen_sorted_wo", "a_testZeroLen.bed", []string{"b_testZeroLen.bed"}, []string{"-wo", "-sorted"}},

		// GFF / VCF inputs (B record type drives the null-B placeholder shape).
		{"gff_self", "b.issue311.gff", []string{"b.issue311.gff"}, nil},
		{"gff_loj", "a.bed", []string{"b.issue311.gff"}, []string{"-loj"}},
		{"gff_wao", "a.bed", []string{"b.issue311.gff"}, []string{"-wao"}},
		{"vcf_self_loj", "b.issue311.vcf", []string{"b.issue311.vcf"}, []string{"-loj"}},
		{"vcf_sv_wo", "a_vcfSVtest.vcf", []string{"b_vcfSVtest.vcf"}, []string{"-wo"}},
		{"vcf_sv_wa_wb", "a_vcfSVtest.vcf", []string{"b_vcfSVtest.vcf"}, []string{"-wa", "-wb"}},
		{"vcf_sv_wa_wb_v", "a_vcfSVtest.vcf", []string{"b_vcfSVtest.vcf"}, []string{"-wa", "-wb", "-v"}},

		// BAM input (always rendered as BED here via -bed).
		{"bam_bed", "a.bam", []string{"b.bed"}, []string{"-bed"}},
		{"bam_bed_wo", "a.bam", []string{"b.bed"}, []string{"-bed", "-wo"}},
		{"bam_bed_wb", "a.bam", []string{"b.bed"}, []string{"-bed", "-wb"}},
		{"bam_bed_s", "a.bam", []string{"b.bed"}, []string{"-bed", "-s"}},
		{"bam_unmapped_v", "a_with_bothUnmapped.bam", []string{"b.bed"}, []string{"-bed", "-v"}},
		{"bam_oneUnmapped", "oneUnmapped.bam", []string{"j1.bed"}, []string{"-bed"}},
		{"bam_as_b", "a.bed", []string{"a.bam"}, nil},

		// Header echo (-header) over plain and compressed A files.
		{"header", "a_withLargeHeader.bed", []string{"b.bed"}, []string{"-header"}},
		{"header_gzip", "a_withLargeHeader_gzipped.bed.gz", []string{"b.bed"}, []string{"-header"}},
		{"no_header_default", "a_withLargeHeader.bed", []string{"b.bed"}, nil},

		// Multiple -b databases: DB-id column, -names, -filenames, -sortout.
		{"multi_wb", "null_a.bed", []string{"null_b.bed", "null_c.bed"}, []string{"-wb"}},
		{"multi_loj", "null_a.bed", []string{"null_b.bed", "null_c.bed"}, []string{"-loj"}},
		{"multi_wao", "null_a.bed", []string{"null_b.bed", "null_c.bed"}, []string{"-wao"}},
		{"multi_wo", "null_a.bed", []string{"null_b.bed", "null_c.bed"}, []string{"-wo"}},
		{"multi_c", "null_a.bed", []string{"null_b.bed", "null_c.bed"}, []string{"-c"}},
		{"multi_C", "null_a.bed", []string{"null_b.bed", "null_c.bed"}, []string{"-C"}},
		{"multi_u", "null_a.bed", []string{"null_b.bed", "null_c.bed"}, []string{"-u"}},
		{"multi_v", "null_a.bed", []string{"null_b.bed", "null_c.bed"}, []string{"-v"}},
		{"multi_loj_names", "null_a.bed", []string{"null_b.bed", "null_c.bed"}, []string{"-loj", "-names", "b", "c"}},
		{"multi_wao_names", "null_a.bed", []string{"null_b.bed", "null_c.bed"}, []string{"-wao", "-names", "b", "c"}},
		{"multi_C_names", "null_a.bed", []string{"null_b.bed", "null_c.bed"}, []string{"-C", "-names", "B", "C"}},
		{"multi_loj_filenames", "null_a.bed", []string{"null_b.bed", "null_c.bed"}, []string{"-loj", "-filenames"}},

		// -gzip'd / line-ending fixtures.
		{"dosLineChar", "dosLineChar_a.bed", []string{"dosLineCharWithExtraTab_b.bed"}, []string{"-v"}},

		// gzip'd VCF A with BED B.
		{"vcf_gz_wa_wb", "bug44_a.vcf.gz", []string{"bug44_b.bed"}, []string{"-wa", "-wb"}},
	}

	for _, tc := range cases {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			upArgs := append([]string{"intersect", "-a", tc.a, "-b"}, tc.b...)
			upArgs = append(upArgs, tc.args...)
			ourArgs := append([]string{"-a", tc.a, "-b"}, tc.b...)
			ourArgs = append(ourArgs, tc.args...)

			upOut, _ := runOut(t, dir, bt, upArgs...)
			ourOut, _ := runOut(t, dir, ours, ourArgs...)
			if !bytes.Equal(upOut, ourOut) {
				t.Fatalf("stdout mismatch for %s (args %v)\n--- upstream ---\n%s\n--- ours ---\n%s",
					tc.name, tc.args, upOut, ourOut)
			}
		})
	}
}

// TestUpstreamFixtures_SortedErrors verifies the -sorted / -g input-order
// validation messages match upstream verbatim (stderr), for both the
// out-of-order and genome-order-violation cases.
func TestUpstreamFixtures_SortedErrors(t *testing.T) {
	bt := upstreamBedtools(t)
	ours := buildOurs(t)
	dir := fixturesDir(t)

	cases := []struct {
		name string
		args []string
	}{
		{"chroms_out_of_order", []string{"-a", "chromsOutOfOrder.bed", "-b", "b.bed", "-sorted"}},
		{"genome_order", []string{"-a", "chromOrderA.bed", "-b", "chromOrderB.bed", "-sorted", "-g", "human.hg19.genome"}},
	}
	for _, tc := range cases {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			_, upErr := runOut(t, dir, bt, append([]string{"intersect"}, tc.args...)...)
			_, ourErr := runOut(t, dir, ours, tc.args...)
			if !bytes.Equal(upErr, ourErr) {
				t.Fatalf("stderr mismatch for %s\n--- upstream ---\n%s\n--- ours ---\n%s",
					tc.name, upErr, ourErr)
			}
		})
	}
}

// TestUpstreamFixtures_NameConventionWarning verifies the chromosome
// naming-convention WARNING (and its -nonamecheck suppression) matches upstream
// on stderr, for the unsorted, sorted, and suppressed cases.
func TestUpstreamFixtures_NameConventionWarning(t *testing.T) {
	bt := upstreamBedtools(t)
	ours := buildOurs(t)
	dir := fixturesDir(t)

	cases := []struct {
		name string
		args []string
	}{
		{"unsorted", []string{"-a", "nonamecheck_a.bed", "-b", "nonamecheck_b.bed"}},
		{"sorted", []string{"-a", "nonamecheck_a.bed", "-b", "nonamecheck_b.bed", "-sorted"}},
		{"suppressed", []string{"-a", "nonamecheck_a.bed", "-b", "nonamecheck_b.bed", "-sorted", "-nonamecheck"}},
	}
	for _, tc := range cases {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			_, upErr := runOut(t, dir, bt, append([]string{"intersect"}, tc.args...)...)
			_, ourErr := runOut(t, dir, ours, tc.args...)
			if !bytes.Equal(upErr, ourErr) {
				t.Fatalf("stderr mismatch for %s\n--- upstream ---\n%q\n--- ours ---\n%q",
					tc.name, upErr, ourErr)
			}
		})
	}
}

// TestUpstreamFixtures_GzipStdin verifies a gzip/BGZF-compressed BED piped to
// `-a -` is transparently decompressed and matches upstream.
func TestUpstreamFixtures_GzipStdin(t *testing.T) {
	bt := upstreamBedtools(t)
	ours := buildOurs(t)
	dir := fixturesDir(t)

	for _, gz := range []string{"a_gzipped.bed.gz", "a_bgzipped.bed.gz"} {
		gz := gz
		t.Run(gz, func(t *testing.T) {
			data, err := os.ReadFile(filepath.Join(dir, gz))
			if err != nil {
				t.Fatalf("read %s: %v", gz, err)
			}
			run := func(name string, args ...string) []byte {
				cmd := exec.Command(name, args...)
				cmd.Dir = dir
				cmd.Stdin = bytes.NewReader(data)
				var out bytes.Buffer
				cmd.Stdout = &out
				_ = cmd.Run()
				return out.Bytes()
			}
			up := run(bt, "intersect", "-a", "-", "-b", "b.bed")
			our := run(ours, "-a", "-", "-b", "b.bed")
			if !bytes.Equal(up, our) {
				t.Fatalf("gzip-stdin mismatch for %s\nupstream:\n%s\nours:\n%s", gz, up, our)
			}
		})
	}
}
