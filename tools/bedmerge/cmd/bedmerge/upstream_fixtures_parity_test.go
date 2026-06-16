package main

// Live-upstream byte-for-byte parity tests over the upstream bedtools merge
// fixtures (reference_code/bedtools/test/merge/*). Each case mirrors a scenario
// from upstream's test-merge.sh and asserts that bedmerge's stdout (and, for the
// error cases, stderr) matches the live `bedtools merge` binary exactly.
//
// The upstream binary is built once via sync.Once (see upstream_compat_test.go's
// upstreamBedtools); a missing binary is a t.Fatalf, never a t.Skip.

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"testing"
)

// mergeFixtureDir returns the absolute path to the upstream merge fixture
// directory, failing the test if the submodule is not checked out.
func mergeFixtureDir(t *testing.T) string {
	t.Helper()
	root, err := repoRoot()
	if err != nil {
		t.Fatalf("locate repo root: %v", err)
	}
	dir := filepath.Join(root, "reference_code", "bedtools", "test", "merge")
	if _, err := os.Stat(dir); err != nil {
		t.Fatalf("upstream merge fixtures unavailable: %v\n"+
			"run: git submodule update --init reference_code/bedtools", err)
	}
	return dir
}

// runStdoutStderr runs name with args and the fixture dir as the working
// directory, returning stdout and stderr separately. A non-zero exit is not a
// fatal error here because several upstream cases intentionally exit non-zero;
// the caller compares whichever stream it cares about.
func runStdoutStderr(t *testing.T, dir, name string, args ...string) (stdout, stderr []byte) {
	t.Helper()
	cmd := exec.Command(name, args...)
	cmd.Dir = dir
	var out, errb bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = &errb
	_ = cmd.Run()
	return out.Bytes(), errb.Bytes()
}

// TestUpstreamParity_MergeFixtures_Stdout asserts byte-for-byte stdout parity
// for every data-producing case in test-merge.sh.
func TestUpstreamParity_MergeFixtures_Stdout(t *testing.T) {
	bt := upstreamBedtools(t)
	ours := buildOurs(t)
	dir := mergeFixtureDir(t)

	cases := []struct {
		name string
		args []string // upstream args minus the leading "merge"
	}{
		{"t1_basic", []string{"-i", "a.bed"}},
		{"t3_count", []string{"-i", "a.bed", "-c", "1", "-o", "count"}},
		{"t5_collapse", []string{"-i", "a.names.bed", "-c", "4", "-o", "collapse"}},
		{"t6_collapse_sum", []string{"-i", "a.full.bed", "-c", "4,5", "-o", "collapse,sum"}},
		{"t7_count_sum", []string{"-i", "a.full.bed", "-c", "5", "-o", "count,sum"}},
		{"t8_three_ops", []string{"-i", "a.full.bed", "-c", "4,5,4", "-o", "collapse,sum,count"}},
		{"t9a_stranded", []string{"-i", "a.full.bed", "-s", "-c", "4,5,6", "-o", "collapse,sum,count"}},
		{"t9b_stranded_strand", []string{"-i", "a.full.bed", "-s", "-c", "4,5,6", "-o", "collapse,sum,collapse"}},
		{"t10_delim", []string{"-i", "a.names.bed", "-delim", "|", "-c", "4", "-o", "collapse"}},
		{"t13_vcf", []string{"-i", "testA.vcf"}},
		{"t14_gff", []string{"-i", "a.gff"}},
		{"t15_mixed_s", []string{"-i", "mixedStrands.bed", "-s"}},
		{"t16_S_plus", []string{"-i", "mixedStrands.bed", "-S", "+"}},
		{"t17_S_minus", []string{"-i", "mixedStrands.bed", "-S", "-"}},
		{"t20_chrom_change", []string{"-i", "b.bed"}},
		{"t21_bed3_from_full", []string{"-i", "a.full.bed"}},
		{"t22_mean_precision", []string{"-i", "precisionTest.bed", "-c", "5", "-o", "mean"}},
		{"t23a_sum_nonnumeric", []string{"-i", "a.names.bed", "-c", "4", "-o", "sum"}},
		{"t43_scientific_coords", []string{"-i", "expFormat.bed"}},
		{"t44a_vcf_sv_len", []string{"-i", "vcfSVtest.vcf"}},
		{"t46_sum_default_prec", []string{"-i", "precisionTest2.bed", "-c", "8", "-o", "sum"}},
		{"t47_sum_prec5", []string{"-i", "precisionTest2.bed", "-c", "8", "-o", "sum", "-prec", "5"}},
		{"t48_bug254_d", []string{"-i", "bug254_d.bed", "-s", "-d", "200"}},
		{"t49_bug254_e", []string{"-i", "bug254_e.bed", "-s", "-d", "200"}},
		{"t50_chained_gzip", []string{"-i", "chained.bed.gz"}},
		// BAM column operations (t24-t36).
		{"t24_bam_c1", []string{"-i", "fullFields.bam", "-c", "1", "-o", "collapse"}},
		{"t26_bam_c3", []string{"-i", "fullFields.bam", "-c", "3", "-o", "collapse"}},
		{"t27_bam_c4_mean", []string{"-i", "fullFields.bam", "-c", "4", "-o", "mean"}},
		{"t28_bam_c5_mean", []string{"-i", "fullFields.bam", "-c", "5", "-o", "mean"}},
		{"t29_bam_c6_cigar", []string{"-i", "fullFields.bam", "-c", "6", "-o", "collapse"}},
		{"t30_bam_c7_rnext", []string{"-i", "fullFields.bam", "-c", "7", "-o", "collapse"}},
		{"t31_bam_c8_mean", []string{"-i", "fullFields.bam", "-c", "8", "-o", "mean"}},
		{"t32_bam_c9_mean", []string{"-i", "fullFields.bam", "-c", "9", "-o", "mean"}},
		{"t33_bam_c10_seq", []string{"-i", "fullFields.bam", "-c", "10", "-o", "collapse"}},
		{"t34_bam_c11_qual", []string{"-i", "fullFields.bam", "-c", "11", "-o", "collapse"}},
		{"t36_bam_missing_mate", []string{"-i", "a.full.bam", "-c", "7", "-o", "collapse"}},
	}

	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			want, _ := runStdoutStderr(t, dir, bt, append([]string{"merge"}, tc.args...)...)
			got, _ := runStdoutStderr(t, dir, ours, tc.args...)
			if !bytes.Equal(got, want) {
				t.Fatalf("stdout mismatch for %v\nupstream:\n%q\nours:\n%q", tc.args, want, got)
			}
		})
	}
}

// TestUpstreamParity_MergeStdin asserts parity when reading from stdin (t45).
func TestUpstreamParity_MergeStdin(t *testing.T) {
	bt := upstreamBedtools(t)
	ours := buildOurs(t)
	dir := mergeFixtureDir(t)
	input, err := os.ReadFile(filepath.Join(dir, "a.bed"))
	if err != nil {
		t.Fatalf("read a.bed: %v", err)
	}
	run := func(bin string, args ...string) []byte {
		cmd := exec.Command(bin, args...)
		cmd.Dir = dir
		cmd.Stdin = bytes.NewReader(input)
		var out bytes.Buffer
		cmd.Stdout = &out
		_ = cmd.Run()
		return out.Bytes()
	}
	want := run(bt, "merge")
	got := run(ours)
	if !bytes.Equal(got, want) {
		t.Fatalf("stdin parity mismatch\nupstream:\n%q\nours:\n%q", want, got)
	}
}

// TestUpstreamParity_MergeErrors asserts byte-for-byte stderr parity for the
// documented error/diagnostic cases (deprecated flags, stranded VCF, bad -S,
// -iobuf validation, and BAM column constraints). Each compares the exact line
// the upstream test pulls from stderr.
func TestUpstreamParity_MergeErrors(t *testing.T) {
	bt := upstreamBedtools(t)
	ours := buildOurs(t)
	dir := mergeFixtureDir(t)

	// lastLine returns the final non-empty stderr line.
	lastLine := func(b []byte) []byte {
		lines := bytes.Split(bytes.TrimRight(b, "\n"), []byte("\n"))
		return lines[len(lines)-1]
	}
	// nthFromTopOfThree returns the middle of the first three stderr lines, which
	// is how the BAM-column upstream tests extract their message (head -3|tail -1).
	middleOfThree := func(b []byte) []byte {
		lines := bytes.Split(b, []byte("\n"))
		if len(lines) >= 3 {
			return lines[2]
		}
		return lastLine(b)
	}

	cases := []struct {
		name    string
		args    []string
		extract func([]byte) []byte
	}{
		{"t2_deprecated_n", []string{"-i", "a.bed", "-n"}, lastLine},
		{"t4_deprecated_nms", []string{"-i", "a.bed", "-nms"}, lastLine},
		{"t11_stranded_vcf", []string{"-i", "testA.vcf", "-s"}, lastLine},
		{"t18_bad_S_arg", []string{"-i", "mixedStrands.bed", "-S", ".", "-c", "6", "-o", "distinct"}, lastLine},
		{"t37_iobuf_missing", []string{"-i", "a.bed", "-iobuf"}, lastLine},
		{"t38_iobuf_bad_suffix", []string{"-i", "a.bed", "-iobuf", "20L"}, lastLine},
		{"t39_iobuf_too_small", []string{"-i", "a.bed", "-iobuf", "7"}, lastLine},
		{"t40_iobuf_nonnumeric", []string{"-i", "a.bed", "-iobuf", "beerM"}, lastLine},
		{"t25_bam_flags_col", []string{"-i", "fullFields.bam", "-c", "2", "-o", "collapse"}, middleOfThree},
		{"t35_bam_col_oob", []string{"-i", "fullFields.bam", "-c", "12", "-o", "sum"}, middleOfThree},
	}

	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			_, wantErr := runStdoutStderr(t, dir, bt, append([]string{"merge"}, tc.args...)...)
			_, gotErr := runStdoutStderr(t, dir, ours, tc.args...)
			want := tc.extract(wantErr)
			got := tc.extract(gotErr)
			if !bytes.Equal(got, want) {
				t.Fatalf("stderr mismatch for %v\nupstream: %q\nours:     %q", tc.args, want, got)
			}
		})
	}
}

// TestUpstreamParity_MergeWarning asserts the non-numeric WARNING lines on
// stderr match upstream byte-for-byte (t23b).
func TestUpstreamParity_MergeWarning(t *testing.T) {
	bt := upstreamBedtools(t)
	ours := buildOurs(t)
	dir := mergeFixtureDir(t)
	args := []string{"-i", "a.names.bed", "-c", "4", "-o", "sum"}
	_, wantErr := runStdoutStderr(t, dir, bt, append([]string{"merge"}, args...)...)
	_, gotErr := runStdoutStderr(t, dir, ours, args...)
	if !bytes.Equal(gotErr, wantErr) {
		t.Fatalf("warning stderr mismatch\nupstream:\n%q\nours:\n%q", wantErr, gotErr)
	}
}

// TestUpstreamParity_MergeHeader asserts `-header` echoes the input header
// before the merged output, byte-for-byte (covers the -header flag surface).
func TestUpstreamParity_MergeHeader(t *testing.T) {
	bt := upstreamBedtools(t)
	ours := buildOurs(t)
	dir := mergeFixtureDir(t)
	args := []string{"-i", "bug254_d.bed", "-header"}
	want, _ := runStdoutStderr(t, dir, bt, append([]string{"merge"}, args...)...)
	got, _ := runStdoutStderr(t, dir, ours, args...)
	if !bytes.Equal(got, want) {
		t.Fatalf("-header stdout mismatch\nupstream:\n%q\nours:\n%q", want, got)
	}
}

// TestUpstreamParity_MergeNewColumnOps asserts the extended KeyListOps
// vocabulary (absmin, absmax, stdev, sstdev, distinct_only, etc.) matches
// upstream over a crafted fixture with repeats and ties.
func TestUpstreamParity_MergeNewColumnOps(t *testing.T) {
	bt := upstreamBedtools(t)
	ours := buildOurs(t)
	dir := t.TempDir()
	in := writeFile(t, dir, "ops.bed",
		"chr1\t10\t20\t-3\n"+
			"chr1\t15\t25\t1\n"+
			"chr1\t24\t30\t-10\n"+
			"chr1\t29\t40\t3\n"+
			"chr1\t38\t50\t1\n")
	ops := []string{
		"absmin", "absmax", "stdev", "sstdev", "median",
		"distinct_sort_num", "distinct_sort_num_desc",
		"freqasc", "freqdesc", "count_distinct", "concat",
	}
	for _, op := range ops {
		t.Run(op, func(t *testing.T) {
			args := []string{"-i", in, "-c", "4", "-o", op}
			want := runCapture(t, bt, append([]string{"merge"}, args...)...)
			got := runCapture(t, ours, args...)
			if !bytes.Equal(got, want) {
				t.Fatalf("op %q mismatch\nupstream: %q\nours:     %q", op, want, got)
			}
		})
	}
}
