package vcftools

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// TestHapcount_CorrectBinTransitions verifies --hapcount produces the
// hand-traced, BUG-FIXED output for a multi-chromosome / multi-bin VCF.
//
// This is NOT a parity test against the upstream vcftools binary — that
// binary's `output_haplotype_count` has two well-known defects (see
// docs/UPSTREAM_BUGS.md and hapcount.go for the writeup). The expected
// fixture is hand-traced from the corrected semantics and pins:
//
//   - The within-chromosome bin-transition flush correctly attributes
//     SNP counts to the OLD bin (not the NEW bin's slot — upstream's
//     prev_bin_idx-shift bug).
//
//   - The last seen chromosome's bins are emitted at end-of-stream
//     (upstream's read-after-free at the post-loop flush block
//     silently drops them on a glibc-built binary).
//
//   - The BED file's first line is auto-detected as a header (starts
//     with '#') and dropped; subsequent data lines are read normally.
func TestHapcount_CorrectBinTransitions(t *testing.T) {
	prefix := runVcftoolsParity(t, "hapcount_fixture.vcf", &Params{
		HapcountBED: filepath.Join(vcftoolsFixtureDir(t), "hapcount_regions.bed"),
		MinAlleles:  2,
	})
	got := readFileBytes(t, prefix+".hapcount")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "hapcount.expected.hapcount"))
	if !bytes.Equal(got, want) {
		t.Errorf(".hapcount mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestHapcount_BEDFirstLineWithData verifies the FIX for the upstream
// BED-first-line silent-skip bug. A header-less BED file (first line is
// real data) must be fully consumed; otherwise the user silently loses
// chr-1 bin 0 from their output.
//
// We run against the same VCF fixture as the main test; since the BED
// content is identical (no '#' header in this variant), the expected
// .hapcount is byte-for-byte the same as the main test's fixture.
func TestHapcount_BEDFirstLineWithData(t *testing.T) {
	prefix := runVcftoolsParity(t, "hapcount_fixture.vcf", &Params{
		HapcountBED: filepath.Join(vcftoolsFixtureDir(t), "hapcount_regions_noheader.bed"),
		MinAlleles:  2,
	})
	got := readFileBytes(t, prefix+".hapcount")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "hapcount.expected.hapcount"))
	if !bytes.Equal(got, want) {
		t.Errorf(".hapcount (header-less BED) mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestHapcount_EndOfStreamFlush verifies the FIX for the upstream
// end-of-stream read-after-free bug. A multi-bin VCF whose LAST variant
// naturally ends mid-chromosome (no following chromosome to trigger the
// chrom-transition flush) must emit all bins of the last chromosome —
// including the final bin's accumulated counts.
//
// Without the fix, upstream silently drops the last chromosome's rows
// (have_data=true at EOF), or emits all-zero rows (have_data=false at
// EOF) — both losing real data.
func TestHapcount_EndOfStreamFlush(t *testing.T) {
	prefix := runVcftoolsParity(t, "hapcount_eof.vcf", &Params{
		HapcountBED: filepath.Join(vcftoolsFixtureDir(t), "hapcount_regions.bed"),
		MinAlleles:  2,
	})
	got := readFileBytes(t, prefix+".hapcount")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "hapcount_eof.expected.hapcount"))
	if !bytes.Equal(got, want) {
		t.Errorf(".hapcount (EOF flush) mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestHapcount_DropsUnphased pins that the implicit `phased_only`
// filter (upstream parameters.cpp:248) drops `0/1`-style sites before
// they reach the hapcount accumulator. The middle of three chr-1 sites
// is `0/1`; the output row for chr 1 bin 0 reflects only the two
// phased sites.
func TestHapcount_DropsUnphased(t *testing.T) {
	prefix := runVcftoolsParity(t, "hapcount_unphased.vcf", &Params{
		HapcountBED: filepath.Join(vcftoolsFixtureDir(t), "hapcount_regions.bed"),
		MinAlleles:  2,
	})
	got := readFileBytes(t, prefix+".hapcount")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "hapcount_unphased.expected.hapcount"))
	if !bytes.Equal(got, want) {
		t.Errorf(".hapcount (unphased) mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestHapcount_OverlappingBED rejects a BED file where two bins on the
// same chromosome overlap. Upstream variant_file_output.cpp:1212-1214
// errors with "BED file must be non-overlapping." (exit code 33). We
// surface a fatal Run-time error rather than producing garbage.
func TestHapcount_OverlappingBED(t *testing.T) {
	tmp := t.TempDir()
	bedPath := filepath.Join(tmp, "overlap.bed")
	if err := os.WriteFile(bedPath, []byte("#hdr\n1\t100\t300\n1\t200\t400\n"), 0o644); err != nil {
		t.Fatalf("write bed: %v", err)
	}
	in, err := os.Open(filepath.Join(vcftoolsFixtureDir(t), "hapcount_fixture.vcf"))
	if err != nil {
		t.Fatalf("open vcf: %v", err)
	}
	defer in.Close()
	err = Run(in, &Params{
		OutPrefix:   filepath.Join(tmp, "out"),
		HapcountBED: bedPath,
		MinAlleles:  2,
	})
	if err == nil {
		t.Fatal("expected Run() to error on overlapping BED, got nil")
	}
	if !strings.Contains(err.Error(), "non-overlapping") {
		t.Errorf("want error mentioning 'non-overlapping', got: %v", err)
	}
}

// TestHapcount_MissingBED reports a clear error when the BED file cannot
// be opened. Mirrors upstream variant_file_output.cpp:1178.
func TestHapcount_MissingBED(t *testing.T) {
	tmp := t.TempDir()
	in, err := os.Open(filepath.Join(vcftoolsFixtureDir(t), "hapcount_fixture.vcf"))
	if err != nil {
		t.Fatalf("open vcf: %v", err)
	}
	defer in.Close()
	err = Run(in, &Params{
		OutPrefix:   filepath.Join(tmp, "out"),
		HapcountBED: filepath.Join(tmp, "does-not-exist.bed"),
		MinAlleles:  2,
	})
	if err == nil {
		t.Fatal("expected Run() to error on missing BED file, got nil")
	}
	if !strings.Contains(err.Error(), "Could not open BED file") {
		t.Errorf("want error 'Could not open BED file', got: %v", err)
	}
}

// TestHapcount_EmptyVCFStillWritesHeader pins that the output file is
// created (header-only) even when no input variant ever touches a BED
// bin — useful as a defensive guarantee for downstream pipelines.
func TestHapcount_EmptyVCFStillWritesHeader(t *testing.T) {
	tmp := t.TempDir()
	vcfPath := filepath.Join(tmp, "empty.vcf")
	if err := os.WriteFile(vcfPath, []byte("##fileformat=VCFv4.0\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"), 0o644); err != nil {
		t.Fatalf("write vcf: %v", err)
	}
	in, err := os.Open(vcfPath)
	if err != nil {
		t.Fatalf("open vcf: %v", err)
	}
	defer in.Close()
	prefix := filepath.Join(tmp, "out")
	if err := Run(in, &Params{
		OutPrefix:   prefix,
		HapcountBED: filepath.Join(vcftoolsFixtureDir(t), "hapcount_regions.bed"),
		MinAlleles:  2,
	}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	got := readFileBytes(t, prefix+".hapcount")
	want := []byte("#CHROM\tBIN_START\tBIN_END\tN_SNP\tN_UNIQ_HAPS\tN_GROUPS\t{MULTIPLICITY:FREQ}\n")
	if !bytes.Equal(got, want) {
		t.Errorf("header-only output mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestShouldSkipBEDHeader directly exercises the header auto-detection
// helper since its behavior is the actual fix for upstream bug #3.
func TestShouldSkipBEDHeader(t *testing.T) {
	cases := []struct {
		line string
		skip bool
	}{
		{"", true},
		{"   ", true},
		{"#CHROM\tSTART\tEND", true},
		{"# any comment", true},
		{"track name=foo", true},
		{"browser position chr1:1-100", true},
		{"1\t100\t500", false},
		{"chr1\t100\t500\textra", false},
		{"X\t0\t1", false},
	}
	for _, c := range cases {
		if got := shouldSkipBEDHeader(c.line); got != c.skip {
			t.Errorf("shouldSkipBEDHeader(%q) = %v, want %v", c.line, got, c.skip)
		}
	}
}

// TestParams_TempDirAccepted is a smoke check that --temp is accepted
// without altering output. It mirrors upstream parameters.cpp:341
// (DIR is stored as the spill-file base path); we don't spill so the
// flag is observably a no-op.
func TestParams_TempDirAccepted(t *testing.T) {
	prefix := runVcftoolsParity(t, "burden_fixture.vcf", &Params{
		IndvBurden: true,
		TempDir:    "/some/path/we/never/touch",
		MinAlleles: 2,
	})
	got := readFileBytes(t, prefix+".iburden")
	want := readFileBytes(t, filepath.Join(vcftoolsFixtureDir(t), "burden.expected.iburden"))
	if !bytes.Equal(got, want) {
		t.Errorf(".iburden (with --temp) mismatch\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// TestGzdiffAliasesDiff smoke-tests that the --gzdiff CLI flag is wired
// as a plain alias for --diff (last-set wins). We don't exercise the
// CLI binary here — that's main.go's concern; instead we assert that
// Run() with Params.Diff set produces a .diff.sites_in_files row,
// which is the behaviour --gzdiff piggybacks on at the CLI layer.
func TestGzdiffAliasesDiff(t *testing.T) {
	tmp := t.TempDir()
	diff2 := filepath.Join(tmp, "f2.vcf")
	if err := os.WriteFile(diff2, []byte(diffVCF2), 0o644); err != nil {
		t.Fatalf("write f2: %v", err)
	}
	prefix := filepath.Join(tmp, "cmp")
	// Simulate `--gzdiff f2.vcf` (Params.Diff is the only slot; main.go
	// overwrites it when --gzdiff is set).
	err := Run(strings.NewReader(diffVCF1), &Params{
		OutPrefix: prefix,
		Diff:      diff2,
		DiffSite:  true,
	})
	if err != nil {
		t.Fatalf("Run --gzdiff (via Diff): %v", err)
	}
	body, err := os.ReadFile(prefix + ".diff.sites_in_files")
	if err != nil {
		t.Fatalf("read sites_in_files: %v", err)
	}
	// Must have at least the header line; full content is exercised
	// by TestRunDiffFamily.
	if !strings.Contains(string(body), "CHROM\tPOS1\tPOS2\tIN_FILE") {
		t.Errorf("missing diff sites_in_files header: %s", body)
	}
}
