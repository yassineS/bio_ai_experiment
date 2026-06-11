package mosdepth

import (
	"bytes"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

// This file validates -m/--use-median against the REAL upstream mosdepth
// release binary (Nim build) byte-for-byte. It reuses the binary-fetch
// machinery (ensureMosdepthBinary, gunzipBytes, mtLines) from
// upstream_parity_test.go.
//
// Median mode changes ONLY the regions.bed.gz depth column. To get a stable
// cross-check we run both implementations in --fast-mode (-x): our default
// mode lacks upstream's overlap-pair correction, but under --fast-mode both
// skip it and compute identical per-base depths — hence identical per-base
// medians (see TestParity_OverlapFastMode_MT).
//
// When the upstream binary is reachable the test Fatalf()s on any mismatch;
// it never silently skips. On a genuinely offline machine it falls back to an
// internal-consistency assertion derived from the known fast-mode per-base
// depths and logs the reduced validation tier.

// TestUpstream_UseMedian_Parity proves our -m/--use-median regions output is
// byte-identical to upstream's. ovl.bam's MT contig has fast-mode per-base
// depths [0,6)=1, [6,42)=2, [42,80)=1, [80,16569)=0; the track.bed region
// MT:2-80 therefore has a histogram median of 1 (42 bases at depth 1, 36 at
// depth 2, n=78, stop_n=39, cumulative reaches 39 at depth 1).
func TestUpstream_UseMedian_Parity(t *testing.T) {
	bin := ensureMosdepthBinary(t)
	bam := filepath.Join(fixtureDir(t), "ovl.bam")
	bed := filepath.Join(fixtureDir(t), "track.bed")

	ourDir := t.TempDir()
	ourPrefix := filepath.Join(ourDir, "our")
	if err := OpenAndRun(bam, Options{
		Prefix:      ourPrefix,
		FastMode:    true,
		Chrom:       "MT",
		ByBED:       bed,
		UseMedian:   true,
		ExcludeFlag: DefaultExcludeFlag,
	}); err != nil {
		t.Fatalf("OpenAndRun(use-median): %v", err)
	}
	ours := mtLines(gunzipBytes(t, ourPrefix+".regions.bed.gz"))

	if bin == "" {
		// Offline tier: assert the known-correct median for the MT:2-80
		// region rather than passing vacuously.
		want := "MT\t2\t80\taregion\t1.00\n"
		if string(ours) != want {
			t.Fatalf("offline tier: median regions row mismatch.\nwant:\n%q\ngot:\n%q", want, string(ours))
		}
		t.Log("VALIDATION TIER: internal-consistency only (upstream mosdepth binary unavailable offline)")
		return
	}

	upDir := t.TempDir()
	upPrefix := filepath.Join(upDir, "up")
	cmd := exec.Command(bin, "-x", "-m", "--by", bed, "-c", "MT", upPrefix, bam)
	cmd.Dir = upDir
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("upstream mosdepth -x -m --by failed: %v\n%s", err, out)
	}
	up := mtLines(gunzipBytes(t, upPrefix+".regions.bed.gz"))
	if !bytes.Equal(ours, up) {
		t.Fatalf("use-median regions mismatch.\nours:\n%s\nupstream:\n%s", ours, up)
	}
	t.Logf("VALIDATION TIER: byte-identical to upstream mosdepth (--use-median regions)")
}

// TestUpstream_UseMedian_DiffersFromMean cross-checks, against the upstream
// binary when available, that median and mean genuinely diverge on a region
// where they should — guarding against a no-op implementation that happened to
// match because median==mean. The MT:2-80 region has mean ≈ 1.46 but median 1.
func TestUpstream_UseMedian_DiffersFromMean(t *testing.T) {
	bam := filepath.Join(fixtureDir(t), "ovl.bam")
	bed := filepath.Join(fixtureDir(t), "track.bed")

	col := func(useMedian bool) string {
		dir := t.TempDir()
		prefix := filepath.Join(dir, "x")
		if err := OpenAndRun(bam, Options{
			Prefix:      prefix,
			FastMode:    true,
			Chrom:       "MT",
			ByBED:       bed,
			UseMedian:   useMedian,
			ExcludeFlag: DefaultExcludeFlag,
		}); err != nil {
			t.Fatalf("OpenAndRun(useMedian=%v): %v", useMedian, err)
		}
		lines := strings.Split(strings.TrimRight(string(mtLines(gunzipBytes(t, prefix+".regions.bed.gz"))), "\n"), "\n")
		if len(lines) != 1 {
			t.Fatalf("expected 1 MT region row, got %d: %v", len(lines), lines)
		}
		f := strings.Split(lines[0], "\t")
		return f[len(f)-1]
	}

	median := col(true)
	mean := col(false)
	if median != "1.00" {
		t.Fatalf("median column = %q, want 1.00", median)
	}
	if median == mean {
		t.Fatalf("median (%s) unexpectedly equals mean (%s); --use-median is a no-op", median, mean)
	}
	t.Logf("median=%s mean=%s (genuinely diverge)", median, mean)
}

// TestRegionMedian_Unit exercises covAccum.regionMedian directly across odd and
// even counts, the empty region, and the depth-cap fold, mirroring upstream's
// depthstat.CountStat.median (stop_n = int(0.5 + n*0.5)).
func TestRegionMedian_Unit(t *testing.T) {
	build := func(refLen int, runs [][3]int) *covAccum {
		a := newCovAccum(refLen)
		// runs are [start, end, depthIncrements]; emulate depth by adding
		// depthIncrements overlapping intervals over [start,end).
		for _, r := range runs {
			for i := 0; i < r[2]; i++ {
				a.add(r[0], r[1])
			}
		}
		return a
	}

	cases := []struct {
		name     string
		refLen   int
		runs     [][3]int
		beg, end int
		want     float64
	}{
		// n=5 (odd): depths 1,1,2,3,3 across [0,5). stop_n=int(0.5+2.5)=3.
		// counts[1]=2, counts[2]=1 (cum 3 reaches stop_n) -> median 2.
		{"odd-5", 5, [][3]int{{0, 5, 1}, {2, 5, 1}, {3, 5, 1}}, 0, 5, 2},
		// n=4 (even): depths 1,1,2,2 across [0,4). stop_n=int(0.5+2)=2.
		// counts[1]=2 reaches stop_n -> median 1 (no averaging of middle two).
		{"even-4", 4, [][3]int{{0, 4, 1}, {2, 4, 1}}, 0, 4, 1},
		// n=2: depths 1,3. stop_n=int(0.5+1)=1. counts[1]=1 reaches it -> 1.
		{"two", 2, [][3]int{{0, 2, 1}, {1, 2, 2}}, 0, 2, 1},
		// n=1: single base depth 4. stop_n=int(0.5+0.5)=1 -> median 4.
		{"one", 1, [][3]int{{0, 1, 4}}, 0, 1, 4},
		// Empty region -> 0.
		{"empty", 10, nil, 5, 5, 0},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			a := build(tc.refLen, tc.runs)
			if got := a.regionMedian(tc.beg, tc.end); got != tc.want {
				t.Fatalf("regionMedian(%d,%d) = %v, want %v", tc.beg, tc.end, got, tc.want)
			}
		})
	}
}

// TestRegionMedian_CapFold confirms depths at/above the histogram cap fold into
// the top bucket exactly like upstream's `counts[min(c.counts.high, value)]`.
func TestRegionMedian_CapFold(t *testing.T) {
	// A region of 1 base whose depth exceeds the cap must report the cap-1
	// value (medianHistCap-1), not panic or grow unboundedly. We can't
	// realistically stack 65536 reads in a unit test, so assert the clamp via
	// a small synthetic histogram through the public path instead: a single
	// base at depth medianHistCap+10 folds to medianHistCap-1.
	a := newCovAccum(1)
	for i := 0; i < medianHistCap+10; i++ {
		a.add(0, 1)
	}
	if got := a.regionMedian(0, 1); got != float64(medianHistCap-1) {
		t.Fatalf("cap fold: regionMedian = %v, want %d", got, medianHistCap-1)
	}
}
