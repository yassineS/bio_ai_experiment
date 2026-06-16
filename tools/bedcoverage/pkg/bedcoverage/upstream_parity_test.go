package bedcoverage

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"sync"
	"testing"
)

// Live-upstream parity tests for `bedtools coverage`.
//
// Unlike parity_test.go (which compares the library against vendored .expected
// snapshots), this file builds BOTH the real upstream `bedtools` binary from
// the vendored submodule AND this port's `bedcoverage` CLI, runs each on the
// same vendored fixtures, and asserts the two produce byte-for-byte identical
// output. Going through the CLI binary is what lets these tests cover the
// gaps that live in argument parsing / stderr text rather than the library:
//
//   - `-abam` BAM input for the A (query) side (coverage.t1 / t1b).
//   - `-sorted` (the chromsweep flag): we accept it and must match upstream's
//     identical output (coverage.t2b..t19b).
//   - the exact mutually-exclusive-modes stderr line (coverage.t14..t17).
//   - the covered-fraction column printed as a float32 with 7 decimals, so the
//     last digit matches upstream's `%0.7f` on a 32-bit float (coverage.t19,
//     e.g. 7/19 -> 0.3684210 not the float64-rounded 0.3684211).
//
// The tests t.Fatalf (never t.Skip): a missing/unbuildable submodule is a hard
// failure, matching the parity-rig policy used across the bed* tools.

var (
	upstreamBedtoolsOnce sync.Once
	upstreamBedtoolsPath string
	upstreamBedtoolsErr  error

	portBinOnce sync.Once
	portBinPath string
	portBinErr  error
)

// repoRoot walks up from this test file to the repo root (dir with go.mod).
func repoRoot(t *testing.T) string {
	t.Helper()
	_, here, _, ok := runtime.Caller(0)
	if !ok {
		t.Fatalf("runtime.Caller failed")
	}
	dir := filepath.Dir(here)
	for {
		if _, err := os.Stat(filepath.Join(dir, "go.mod")); err == nil {
			return dir
		}
		parent := filepath.Dir(dir)
		if parent == dir {
			t.Fatalf("could not locate repo root (go.mod) above %s", here)
		}
		dir = parent
	}
}

// upstreamBedtools builds (once) and returns the path to the upstream
// `bedtools` binary.
func upstreamBedtools(t *testing.T) string {
	t.Helper()
	upstreamBedtoolsOnce.Do(func() {
		root := repoRoot(t)
		dir := filepath.Join(root, "reference_code", "bedtools")
		bin := filepath.Join(dir, "bin", "bedtools")
		if _, err := os.Stat(bin); err == nil {
			upstreamBedtoolsPath = bin
			return
		}
		cmd := exec.Command("make", "-j", "4")
		cmd.Dir = dir
		if out, err := cmd.CombinedOutput(); err != nil {
			upstreamBedtoolsErr = err
			t.Logf("bedtools build output:\n%s", out)
			return
		}
		if _, err := os.Stat(bin); err != nil {
			upstreamBedtoolsErr = err
			return
		}
		upstreamBedtoolsPath = bin
	})
	if upstreamBedtoolsErr != nil {
		t.Fatalf("building upstream bedtools: %v (run `git submodule update --init reference_code/bedtools`)", upstreamBedtoolsErr)
	}
	if upstreamBedtoolsPath == "" {
		t.Fatalf("upstream bedtools binary not found after build")
	}
	return upstreamBedtoolsPath
}

// portBin builds (once) and returns the path to this port's bedcoverage CLI.
func portBin(t *testing.T) string {
	t.Helper()
	portBinOnce.Do(func() {
		root := repoRoot(t)
		// Build once into the OS temp dir so the binary outlives the per-test
		// t.TempDir()s and can be reused across every parity subtest.
		out := filepath.Join(os.TempDir(), "bedcoverage_parity_bin")
		cmd := exec.Command("go", "build", "-o", out, "./tools/bedcoverage/cmd/bedcoverage")
		cmd.Dir = root
		if combined, err := cmd.CombinedOutput(); err != nil {
			portBinErr = err
			t.Logf("port build output:\n%s", combined)
			return
		}
		portBinPath = out
	})
	if portBinErr != nil {
		t.Fatalf("building port bedcoverage: %v", portBinErr)
	}
	if portBinPath == "" {
		t.Fatalf("port bedcoverage binary not found after build")
	}
	return portBinPath
}

// fixtureAbs returns the absolute path to a testdata/parity fixture.
func fixtureAbs(t *testing.T, name string) string {
	t.Helper()
	p, err := filepath.Abs(filepath.Join("..", "..", "testdata", "parity", name))
	if err != nil {
		t.Fatalf("abs path for %s: %v", name, err)
	}
	return p
}

// run executes bin with args and returns (stdout, stderr) without failing on a
// non-zero exit (the error-text cases expect a non-zero exit).
func run(t *testing.T, bin string, args ...string) (stdout, stderr []byte) {
	t.Helper()
	cmd := exec.Command(bin, args...)
	var so, se bytes.Buffer
	cmd.Stdout = &so
	cmd.Stderr = &se
	_ = cmd.Run()
	return so.Bytes(), se.Bytes()
}

// assertStdoutParity runs upstream `bedtools coverage <args>` and the port's
// `bedcoverage <portArgs>` and asserts identical stdout. portArgs defaults to
// args when nil (the flag surface is shared for these cases).
func assertStdoutParity(t *testing.T, args []string, portArgs []string) {
	t.Helper()
	up := upstreamBedtools(t)
	pb := portBin(t)
	if portArgs == nil {
		portArgs = args
	}
	wantOut, _ := run(t, up, append([]string{"coverage"}, args...)...)
	gotOut, _ := run(t, pb, portArgs...)
	if !bytes.Equal(wantOut, gotOut) {
		t.Fatalf("stdout mismatch for args %v\nupstream:\n%s\nport:\n%s", args, wantOut, gotOut)
	}
}

// TestUpstreamParity_ABAM_QueryBAM covers coverage.t1: BAM input for A via the
// legacy -abam flag, plain BED on -b, default output mode. The expected line
// includes the full BED12 echo (trailing-comma block lists) of the alignment.
func TestUpstreamParity_ABAM_QueryBAM(t *testing.T) {
	a := fixtureAbs(t, "cov_a.bam")
	b := fixtureAbs(t, "cov_b.bed")
	assertStdoutParity(t,
		[]string{"-abam", a, "-b", b},
		[]string{"-abam", a, "-b", b},
	)
}

// TestUpstreamParity_ABAM_Sorted covers coverage.t1b: -abam plus -sorted. The
// port treats -sorted as a no-op and must match upstream's identical output.
func TestUpstreamParity_ABAM_Sorted(t *testing.T) {
	a := fixtureAbs(t, "cov_a.bam")
	b := fixtureAbs(t, "cov_b.bed")
	assertStdoutParity(t,
		[]string{"-abam", a, "-b", b, "-sorted"},
		[]string{"-abam", a, "-b", b, "-sorted"},
	)
}

// TestUpstreamParity_Sorted_AllModes covers coverage.t2b..t9b/t18b/t19b: the
// -sorted flag with every output mode and the strand/hist/depth/mean variants.
// Our interval-tree path must produce the exact bytes upstream's chromsweep does.
func TestUpstreamParity_Sorted_AllModes(t *testing.T) {
	a := fixtureAbs(t, "a.bed")
	b := fixtureAbs(t, "b.bed")
	x := fixtureAbs(t, "x.bed")
	y := fixtureAbs(t, "y.bed")
	cases := []struct {
		name string
		args []string
	}{
		{"default", []string{"-a", a, "-b", b, "-sorted"}},
		{"counts", []string{"-a", a, "-b", b, "-counts", "-sorted"}},
		{"hist", []string{"-a", a, "-b", b, "-hist", "-sorted"}},
		{"depth", []string{"-a", a, "-b", b, "-d", "-sorted"}},
		{"mean", []string{"-a", a, "-b", b, "-mean", "-sorted"}},
		{"sameStrand", []string{"-a", a, "-b", b, "-s", "-sorted"}},
		{"oppStrand", []string{"-a", a, "-b", b, "-S", "-sorted"}},
		{"depthLastNoOverlap", []string{"-a", x, "-b", y, "-d", "-sorted"}},
		{"histLastNoOverlap", []string{"-a", x, "-b", y, "-hist", "-sorted"}},
	}
	for _, tc := range cases {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			assertStdoutParity(t, tc.args, tc.args)
		})
	}
}

// TestUpstreamParity_FractionColumn_Float32 covers coverage.t19: the
// covered-fraction column must be computed as a float32 before printing with
// 7 decimals. 7/19 prints as 0.3684210 (float32) not 0.3684211 (float64).
func TestUpstreamParity_FractionColumn_Float32(t *testing.T) {
	x := fixtureAbs(t, "x.bed")
	y := fixtureAbs(t, "y.bed")
	assertStdoutParity(t,
		[]string{"-a", x, "-b", y, "-hist"},
		[]string{"-a", x, "-b", y, "-hist"},
	)
}

// TestUpstreamParity_MutuallyExclusiveError covers coverage.t14..t17: any pair
// of -counts/-d/-mean/-hist must produce upstream's exact stderr line. The test
// mirrors the upstream script's `2>&1 >/dev/null | tail -1` by comparing the
// last stderr line of each binary.
func TestUpstreamParity_MutuallyExclusiveError(t *testing.T) {
	a := fixtureAbs(t, "a.bed")
	b := fixtureAbs(t, "b.bed")
	up := upstreamBedtools(t)
	pb := portBin(t)
	combos := [][]string{
		{"-counts", "-hist"},
		{"-counts", "-d"},
		{"-hist", "-d"},
		{"-mean", "-d"},
	}
	for _, c := range combos {
		c := c
		t.Run(c[0]+c[1], func(t *testing.T) {
			args := append([]string{"-a", a, "-b", b}, c...)
			_, upErr := run(t, up, append([]string{"coverage"}, args...)...)
			_, pErr := run(t, pb, args...)
			wantLast := lastNonEmptyLine(upErr)
			gotLast := lastNonEmptyLine(pErr)
			if !bytes.Equal(wantLast, gotLast) {
				t.Fatalf("stderr last-line mismatch for %v\nupstream: %q\nport:     %q", c, wantLast, gotLast)
			}
		})
	}
}

// lastNonEmptyLine returns the final non-empty line of b without its newline,
// matching the upstream test harness's `tail -1` on the captured stderr.
func lastNonEmptyLine(b []byte) []byte {
	lines := bytes.Split(b, []byte("\n"))
	for i := len(lines) - 1; i >= 0; i-- {
		if len(bytes.TrimSpace(lines[i])) > 0 {
			return lines[i]
		}
	}
	return nil
}
