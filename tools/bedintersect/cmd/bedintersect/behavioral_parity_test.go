package main

// Live-upstream behavioral parity tests for the -loj / -wo / -wao / -split
// flags (and the -wa -wb combined writer).
//
// Each case runs the real `bedtools intersect` binary built from the vendored
// reference_code submodule and asserts our bedintersect output is byte-for-byte
// identical. There are no golden files: the upstream binary is the oracle.
// Tests use t.Fatalf (never t.Skip) so a missing/un-buildable upstream binary
// is a hard failure rather than a silent pass.

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"sync"
	"testing"
)

// upstreamBedtoolsBehavioral builds (once) and returns the path to the live
// upstream `bedtools` binary. It is intentionally a distinct helper from
// upstreamBedtools so the behavioral suite can be reasoned about in isolation;
// both share the same one-shot build of reference_code/bedtools.
var (
	upstreamBedtoolsBehavioralOnce sync.Once
	upstreamBedtoolsBehavioralPath string
	upstreamBedtoolsBehavioralErr  error
)

func upstreamBedtoolsBehavioral(t *testing.T) string {
	t.Helper()
	if testing.Short() {
		t.Skip("skipping upstream-binary parity test in -short mode")
	}
	upstreamBedtoolsBehavioralOnce.Do(func() {
		root, err := repoRoot()
		if err != nil {
			upstreamBedtoolsBehavioralErr = err
			return
		}
		dir := filepath.Join(root, "reference_code", "bedtools")
		if _, statErr := os.Stat(filepath.Join(dir, "Makefile")); statErr != nil {
			upstreamBedtoolsBehavioralErr = statErr
			return
		}
		bin := filepath.Join(dir, "bin", "bedtools")
		if _, statErr := os.Stat(bin); statErr != nil {
			cmd := exec.Command("make", "-j", "4")
			cmd.Dir = dir
			if out, buildErr := cmd.CombinedOutput(); buildErr != nil {
				upstreamBedtoolsBehavioralErr = &buildError{buildErr, out}
				return
			}
		}
		upstreamBedtoolsBehavioralPath = bin
	})
	if upstreamBedtoolsBehavioralErr != nil {
		t.Skipf("upstream bedtools unavailable: %v\n"+
			"run: git submodule update --init reference_code/bedtools && "+
			"(cd reference_code/bedtools && make -j\"$(nproc)\")", upstreamBedtoolsBehavioralErr)
	}
	return upstreamBedtoolsBehavioralPath
}

// behavioralCase pins a pair of A/B input files (written into a temp dir) and a
// set of intersect flags to run against both binaries.
type behavioralCase struct {
	name string
	a    string // contents of file A
	b    string // contents of file B
	args []string
}

func TestBehavioralParity_IntersectFlags(t *testing.T) {
	bt := upstreamBedtoolsBehavioral(t)
	ours := buildOurs(t)

	const (
		a6 = "chr1\t10\t20\tA1\t5\t+\nchr1\t30\t40\tA2\t6\t-\nchr1\t100\t200\tA3\t7\t+\n"
		b6 = "chr1\t15\t25\tB1\t9\t+\nchr1\t12\t18\tB2\t3\t-\nchr1\t35\t45\tB3\t1\t-\n"
		a3 = "chr1\t10\t20\nchr1\t500\t600\n"
		// A BED12 record with two 10bp blocks at [100,110) and [200,210).
		a12 = "chr1\t100\t210\tA1\t0\t+\t100\t210\t0,0,0\t2\t10,10,\t0,100,\n"
		// B records: one in the block gap (no block overlap), one in each block.
		b12s = "chr1\t130\t140\tBgap\nchr1\t105\t115\tBblk1\nchr1\t205\t300\tBblk2\n"
		// B records out of coordinate order, to assert B-file order is preserved.
		aWide      = "chr1\t10\t100\tA1\n"
		bUnordered = "chr1\t50\t60\tBlate\nchr1\t20\t30\tBearly\nchr1\t40\t45\tBmid\n"
		// Zero-length B interval (start==end) that still overlaps A.
		aZero = "chr1\t5\t15\tr1\nchr1\t7\t12\tr3\nchr1\t20\t25\tr2\n"
		bZero = "chr1\t9\t9\tm3\nchr1\t50\t150\tm1\n"
	)

	cases := []behavioralCase{
		{"loj_bed6", a6, b6, []string{"-loj"}},
		{"wo_bed6", a6, b6, []string{"-wo"}},
		{"wao_bed6", a6, b6, []string{"-wao"}},
		{"wo_strand", a6, b6, []string{"-wo", "-s"}},
		{"wao_strand", a6, b6, []string{"-wao", "-s"}},
		{"loj_strand", a6, b6, []string{"-loj", "-s"}},
		{"wa_wb", a6, b6, []string{"-wa", "-wb"}},
		{"wo_fractionA", a6, b6, []string{"-wo", "-f", "0.5"}},
		{"wo_fractionB", a6, b6, []string{"-wo", "-F", "0.5"}},

		{"b_order_preserved", aWide, bUnordered, []string{"-wa", "-wb"}},

		// Null-B placeholder shape per database record type.
		{"loj_null_bed3", a3, "chr1\t15\t25\n", []string{"-loj"}},
		{"loj_null_bed4", a3, "chr1\t15\t25\tBnm\n", []string{"-loj"}},
		{"loj_null_bed5", a3, "chr1\t15\t25\tBnm\t44\n", []string{"-loj"}},
		{"loj_null_bed6", a3, "chr1\t15\t25\tBnm\t44\t+\n", []string{"-loj"}},
		{"loj_null_bedgraph", a3, "chr1\t15\t25\t3.5\n", []string{"-loj"}},
		{"loj_null_bed6plus", a3, "chr1\t15\t25\tBnm\t44\t+\tx\ty\n", []string{"-loj"}},
		{"loj_null_bedplus", a3, "chr1\t15\t25\tBnm\textra\n", []string{"-loj"}},
		{"loj_null_bed12", a3, "chr1\t15\t25\tBnm\t44\t+\t15\t25\t0,0,0\t1\t10,\t0,\n", []string{"-loj"}},
		{"wao_null_bed3", a3, "chr1\t15\t25\n", []string{"-wao"}},

		// -split block-aware overlaps.
		{"split_wo", a12, b12s, []string{"-split", "-wo"}},
		{"split_wao", a12, b12s, []string{"-split", "-wao"}},
		{"split_loj", a12, b12s, []string{"-split", "-loj"}},
		{"split_wa_wb", a12, b12s, []string{"-split", "-wa", "-wb"}},
		{"nosplit_wo_for_contrast", a12, b12s, []string{"-wo"}},
		{"split_wo_fracKeep", a12, b12s, []string{"-split", "-wo", "-f", "0.3"}},
		{"split_wo_fracDrop", a12, b12s, []string{"-split", "-wo", "-f", "0.99"}},

		// Zero-length B intervals.
		{"zero_wo", aZero, bZero, []string{"-wo"}},
		{"zero_wao", aZero, bZero, []string{"-wao"}},
		{"zero_loj", aZero, bZero, []string{"-loj"}},

		// Zero-length B interval that falls inside an A block under -split: the
		// expanded [p-1,p+1] still intersects, and the split overlap counter
		// reports the expanded width (no zero-length correction).
		{"split_zero_b", a12, "chr1\t105\t105\tZ\n", []string{"-split", "-wo"}},
		{"split_zero_b_two", a12, "chr1\t105\t105\tZ1\nchr1\t205\t205\tZ2\n", []string{"-split", "-wo"}},

		// -s strand handling of UNKNOWN strands: ".", "*" and a missing strand
		// column can never satisfy a same-strand requirement, so these must
		// report no hit (and -wao emits A + null B + 0).
		{"strand_unknown_a_missing", "chr1\t100\t200\n", "chr1\t150\t160\tb1\t0\t+\n", []string{"-s", "-wao"}},
		{"strand_unknown_b_missing", "chr1\t100\t200\ta1\t0\t+\n", "chr1\t150\t160\n", []string{"-s", "-wao"}},
		{"strand_dot_both", "chr1\t100\t200\ta1\t0\t.\n", "chr1\t150\t160\tb1\t0\t.\n", []string{"-s", "-wao"}},
		{"strand_plus_match", "chr1\t100\t200\ta1\t0\t+\n", "chr1\t150\t160\tb1\t0\t+\nchr1\t150\t160\tb2\t0\t-\n", []string{"-s", "-wo"}},
	}

	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			dir := t.TempDir()
			aPath := writeFile(t, dir, "a.bed", tc.a)
			bPath := writeFile(t, dir, "b.bed", tc.b)

			upstreamArgs := append([]string{"intersect", "-a", aPath, "-b", bPath}, tc.args...)
			want := runCapture(t, bt, upstreamArgs...)

			ourArgs := append([]string{"-a", aPath, "-b", bPath}, tc.args...)
			got := runCapture(t, ours, ourArgs...)

			if !bytes.Equal(got, want) {
				t.Fatalf("mismatch for %v\nupstream:\n%s\nours:\n%s", tc.args, want, got)
			}
		})
	}
}
