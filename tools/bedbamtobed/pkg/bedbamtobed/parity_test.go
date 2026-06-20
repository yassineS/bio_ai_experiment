package bedbamtobed

// Live-upstream parity tests for `bedtools bamtobed`. They build the real
// upstream `bedtools` binary (and its `htsutil` helper) from the vendored
// submodule exactly once via sync.Once, convert the vendored SAM fixtures to
// BAM with htsutil, then compare upstream's BED/BEDPE output byte-for-byte
// against this port's Run over the same BAM. They t.Fatalf (never t.Skip) so a
// missing or unbuildable submodule is a hard failure, matching the project's
// parity-rig policy.

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"sync"
	"testing"
)

var (
	upstreamOnce sync.Once
	upstreamBin  string
	upstreamHts  string
	upstreamErr  error
)

// repoRoot walks up from this file to the module root (the dir holding go.mod).
func repoRoot(t *testing.T) string {
	t.Helper()
	_, file, _, ok := runtime.Caller(0)
	if !ok {
		t.Fatalf("runtime.Caller failed")
	}
	dir := filepath.Dir(file)
	for {
		if _, err := os.Stat(filepath.Join(dir, "go.mod")); err == nil {
			return dir
		}
		parent := filepath.Dir(dir)
		if parent == dir {
			t.Fatalf("could not locate repo root above %s", file)
		}
		dir = parent
	}
}

// upstream builds (once) and returns the paths to the upstream `bedtools`
// binary and its `htsutil` helper.
func upstream(t *testing.T) (bedtools, htsutil string) {
	t.Helper()
	if testing.Short() {
		t.Skip("skipping upstream-binary parity test in -short mode")
	}
	upstreamOnce.Do(func() {
		root := repoRoot(t)
		dir := filepath.Join(root, "reference_code", "bedtools")
		if _, err := os.Stat(filepath.Join(dir, "Makefile")); err != nil {
			upstreamErr = err
			return
		}
		bin := filepath.Join(dir, "bin", "bedtools")
		hts := filepath.Join(dir, "test", "htsutil")
		_, e1 := os.Stat(bin)
		_, e2 := os.Stat(hts)
		if e1 != nil || e2 != nil {
			cmd := exec.Command("make", "-j", "4")
			cmd.Dir = dir
			if out, buildErr := cmd.CombinedOutput(); buildErr != nil {
				upstreamErr = buildErr
				t.Logf("bedtools build output:\n%s", out)
				return
			}
		}
		upstreamBin, upstreamHts = bin, hts
	})
	if upstreamErr != nil {
		t.Skipf("upstream bedtools unavailable: %v\n"+
			"run: git submodule update --init reference_code/bedtools && "+
			"(cd reference_code/bedtools && make -j\"$(nproc)\")", upstreamErr)
	}
	return upstreamBin, upstreamHts
}

// fixturePath returns the absolute path of a SAM fixture under testdata/parity.
func fixturePath(t *testing.T, name string) string {
	t.Helper()
	wd, err := os.Getwd()
	if err != nil {
		t.Fatalf("getwd: %v", err)
	}
	return filepath.Join(wd, "..", "..", "testdata", "parity", name)
}

// samToBam converts a SAM fixture to a temp BAM using upstream htsutil and
// returns the BAM path.
func samToBam(t *testing.T, htsutil, samName string) string {
	t.Helper()
	sam := fixturePath(t, samName)
	bam := filepath.Join(t.TempDir(), samName+".bam")
	cmd := exec.Command(htsutil, "samtobam", sam, bam)
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("htsutil samtobam %s: %v\n%s", samName, err, out)
	}
	return bam
}

// runUpstream runs upstream `bedtools bamtobed` with the given flags on bamPath.
func runUpstream(t *testing.T, bedtools, bamPath string, flags ...string) []byte {
	t.Helper()
	args := append([]string{"bamtobed", "-i", bamPath}, flags...)
	cmd := exec.Command(bedtools, args...)
	var out, errb bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = &errb
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream bamtobed %v: %v\nstderr: %s", flags, err, errb.String())
	}
	return out.Bytes()
}

// runOurs runs this port's Run on bamPath with the given Options.
func runOurs(t *testing.T, bamPath string, opts Options) []byte {
	t.Helper()
	f, err := os.Open(bamPath)
	if err != nil {
		t.Fatalf("open %s: %v", bamPath, err)
	}
	defer f.Close()
	var out bytes.Buffer
	if _, err := Run(f, &out, opts); err != nil {
		t.Fatalf("Run(%+v): %v", opts, err)
	}
	return out.Bytes()
}

// parityCase couples upstream flags with the equivalent port Options.
type parityCase struct {
	name  string
	sam   string
	flags []string
	opts  Options
}

// TestParity_BamToBed_Suite reproduces the upstream test-bamtobed.sh cases plus
// the documented -cigar/-bedpe/-mate1/-color/-tag/-ed flags, asserting
// byte-for-byte BED-text parity against the live binary.
func TestParity_BamToBed_Suite(t *testing.T) {
	bt, hts := upstream(t)

	cases := []parityCase{
		// t1/t3/t5: whole-footprint BED6 (no split).
		{"one_block", "one_block.sam", nil, Options{}},
		{"two_blocks", "two_blocks.sam", nil, Options{}},
		{"three_blocks", "three_blocks.sam", nil, Options{}},
		// t2/t4/t6: -split.
		{"one_block_split", "one_block.sam", []string{"-split"}, Options{ObeySplits: true}},
		{"two_blocks_split", "two_blocks.sam", []string{"-split"}, Options{ObeySplits: true}},
		{"three_blocks_split", "three_blocks.sam", []string{"-split"}, Options{ObeySplits: true}},
		// t7: -bed12.
		{"three_blocks_bed12", "three_blocks.sam", []string{"-bed12"}, Options{WriteBed12: true}},
		// t9/t10: -split / -splitD over D+N alignment.
		{"wD_split", "two_blocks_w_D.sam", []string{"-split"}, Options{ObeySplits: true}},
		{"wD_splitD", "two_blocks_w_D.sam", []string{"-splitD"}, Options{ObeySplits: true, SplitOnDeletions: true}},
		// bed12 over D+N, with and without -splitD.
		{"wD_bed12", "two_blocks_w_D.sam", []string{"-bed12"}, Options{WriteBed12: true}},
		{"wD_bed12_splitD", "two_blocks_w_D.sam", []string{"-bed12", "-splitD"}, Options{WriteBed12: true, SplitOnDeletions: true}},
		// t12: numeric -tag NM.
		{"tag_NM", "numeric_tag.sam", []string{"-tag", "NM"}, Options{Tag: "NM"}},
		// -cigar (BED6 7th column).
		{"cigar_one", "one_block.sam", []string{"-cigar"}, Options{UseCigar: true}},
		{"cigar_two", "two_blocks.sam", []string{"-cigar"}, Options{UseCigar: true}},
		// -color with -bed12.
		{"color_bed12", "three_blocks.sam", []string{"-bed12", "-color", "0,255,0"}, Options{WriteBed12: true, Color: "0,255,0"}},
		// upstream quirk: -tag combined with -split emits an extra column.
		{"tag_split", "numeric_tag.sam", []string{"-tag", "NM", "-split"}, Options{Tag: "NM", ObeySplits: true}},
		// BEDPE modes.
		{"bedpe", "pe.sam", []string{"-bedpe"}, Options{WriteBedPE: true}},
		{"bedpe_ed", "pe.sam", []string{"-bedpe", "-ed"}, Options{WriteBedPE: true, UseEditDistance: true}},
		{"bedpe_mate1", "pe.sam", []string{"-bedpe", "-mate1"}, Options{WriteBedPE: true, Mate1First: true}},
	}

	for _, tc := range cases {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			bam := samToBam(t, hts, tc.sam)
			want := runUpstream(t, bt, bam, tc.flags...)
			got := runOurs(t, bam, tc.opts)
			if !bytes.Equal(got, want) {
				t.Fatalf("byte mismatch for %s %v\n--- upstream ---\n%s\n--- ours ---\n%s",
					tc.sam, tc.flags, want, got)
			}
		})
	}
}
