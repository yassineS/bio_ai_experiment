package bcftools

// Live upstream-binary parity for `bcftools view -x/--private` and
// `-X/--exclude-private`.
//
// Unlike the table-driven unit tests (view_test.go), this test runs the
// *actual upstream C binary* on the same fixture and compares its record
// selection against our Go port in-process. No committed golden/snapshot
// file is involved: the expected output is produced live by the binary the
// test locates (or builds) under reference_code/bcftools.
//
// By project rule the test must always run remotely; it therefore never
// t.Skip. If the upstream binary genuinely cannot be located or built it
// fails loudly via t.Fatalf so the gap is visible in CI.

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"sync"
	"testing"
	"time"
)

// upstreamBcftoolsOnce guards the build-or-locate work so a single test run
// builds the binary at most once even when invoked from several subtests.
var (
	upstreamBcftoolsOnce sync.Once
	upstreamBcftoolsPath string
	upstreamBcftoolsErr  error
)

// referenceCodeDir resolves the absolute path of the repo's reference_code/
// directory relative to this test file (tools/bcftools/pkg/bcftools).
func referenceCodeDir(t *testing.T) string {
	t.Helper()
	abs, err := filepath.Abs(filepath.Join("..", "..", "..", "..", "reference_code"))
	if err != nil {
		t.Fatalf("resolve reference_code dir: %v", err)
	}
	return abs
}

// runWithRetry runs name/args, retrying up to len(backoff)+1 times with the
// given backoff delays between attempts. It is used for the
// network-dependent `git submodule update` step so a transient fetch failure
// does not fail the whole parity build.
func runWithRetry(t *testing.T, backoff []time.Duration, name string, args ...string) error {
	t.Helper()
	var lastErr error
	for attempt := 0; attempt <= len(backoff); attempt++ {
		cmd := exec.Command(name, args...)
		out, err := cmd.CombinedOutput()
		if err == nil {
			return nil
		}
		lastErr = err
		t.Logf("attempt %d: %s %s failed: %v\n%s", attempt+1, name, strings.Join(args, " "), err, out)
		if attempt < len(backoff) {
			time.Sleep(backoff[attempt])
		}
	}
	return lastErr
}

// run executes name/args in dir and fails the test on error, surfacing the
// combined output for diagnosis.
func run(t *testing.T, dir, name string, args ...string) {
	t.Helper()
	cmd := exec.Command(name, args...)
	cmd.Dir = dir
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("%s %s (in %s) failed: %v\n%s", name, strings.Join(args, " "), dir, err, out)
	}
}

// ensureUpstreamBcftools returns the path to a usable upstream bcftools
// binary, building it from the vendored submodules if necessary. The work is
// done at most once per test process. On unrecoverable failure it returns a
// non-nil error and the caller t.Fatalf's — the project rule forbids skipping.
func ensureUpstreamBcftools(t *testing.T) (string, error) {
	t.Helper()
	upstreamBcftoolsOnce.Do(func() {
		refDir := referenceCodeDir(t)
		bcftoolsDir := filepath.Join(refDir, "bcftools")
		htslibDir := filepath.Join(refDir, "htslib")
		binPath := filepath.Join(bcftoolsDir, "bcftools")

		// 1. Fast path: a prebuilt binary is already present (cache reuse).
		if fi, err := os.Stat(binPath); err == nil && !fi.IsDir() {
			upstreamBcftoolsPath = binPath
			return
		}

		// 2. Ensure the submodule sources are checked out. The bcftools
		//    Makefile expects ../htslib, so init both. htslib in turn needs
		//    its own nested htscodecs submodule, hence --recursive.
		backoff := []time.Duration{2 * time.Second, 4 * time.Second, 8 * time.Second, 16 * time.Second}
		if _, err := os.Stat(filepath.Join(htslibDir, "Makefile")); err != nil {
			if err := runWithRetry(t, backoff, "git", "submodule", "update", "--init", "--recursive",
				htslibDir, bcftoolsDir); err != nil {
				upstreamBcftoolsErr = err
				return
			}
		}
		// htscodecs is a submodule of htslib; ensure it too (idempotent).
		if err := runWithRetry(t, backoff, "git", "-C", htslibDir, "submodule", "update", "--init", "--recursive"); err != nil {
			upstreamBcftoolsErr = err
			return
		}

		// 3. Build htslib (autoreconf + configure + make), then bcftools.
		run(t, htslibDir, "autoreconf", "-i")
		run(t, htslibDir, "./configure")
		run(t, htslibDir, "make", "-j")
		run(t, bcftoolsDir, "make", "-j")

		if fi, err := os.Stat(binPath); err != nil || fi.IsDir() {
			upstreamBcftoolsErr = err
			return
		}
		upstreamBcftoolsPath = binPath
	})
	return upstreamBcftoolsPath, upstreamBcftoolsErr
}

// runUpstreamView runs the upstream bcftools view on path with the given
// extra flags and returns stdout.
func runUpstreamView(t *testing.T, bin, path string, extraFlags ...string) []byte {
	t.Helper()
	args := append([]string{"view", "--no-version"}, extraFlags...)
	args = append(args, path)
	cmd := exec.Command(bin, args...)
	var stdout, stderr bytes.Buffer
	cmd.Stdout = &stdout
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream bcftools %v failed: %v\n%s", args, err, stderr.String())
	}
	return stdout.Bytes()
}

// TestView_PrivateUpstreamParity runs the upstream C binary and our Go port
// on the SAME fixture and asserts the record SELECTION (and every column the
// private filter governs) matches in-process — no committed snapshot.
//
// Upstream additionally recomputes INFO/AC/AN after sample subsetting (a
// separate, pre-existing gap tracked by TestParityView_SampleSubset). The
// private filter does not govern the INFO column, so we blank INFO (field 7)
// on both sides via dataRecordsStripINFO before comparing; every other
// column — CHROM/POS/ID/REF/ALT/QUAL/FILTER/FORMAT and the retained samples'
// genotypes — must match byte-for-byte.
func TestView_PrivateUpstreamParity(t *testing.T) {
	bin, err := ensureUpstreamBcftools(t)
	if err != nil || bin == "" {
		// Never t.Skip: the rule is this test must run remotely and in CI.
		t.Fatalf("could not locate or build the upstream bcftools binary: %v\n"+
			"the parity build needs the htslib + bcftools submodules and the htslib "+
			"build toolchain (autoconf/automake/gcc/zlib/bz2/lzma headers)", err)
	}

	fixture := parityPath(t, filepath.Join("view", "private.vcf"))
	in, readErr := os.ReadFile(fixture)
	if readErr != nil {
		t.Fatalf("read fixture %s: %v", fixture, readErr)
	}

	cases := []struct {
		name string
		flag string
		opts ViewOptions
	}{
		{"private", "-x", ViewOptions{Samples: []string{"S1", "S2"}, Private: true}},
		{"exclude-private", "-X", ViewOptions{Samples: []string{"S1", "S2"}, ExcludePrivate: true}},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			upstream := runUpstreamView(t, bin, fixture, "-s", "S1,S2", tc.flag)
			got := runParityView(t, in, tc.opts)

			wantRecs := dataRecordsStripINFO(string(upstream))
			gotRecs := dataRecordsStripINFO(string(got))
			if !equalStrings(gotRecs, wantRecs) {
				t.Fatalf("record selection mismatch vs live upstream bcftools.\nflag: %s\nwant: %v\ngot:  %v",
					tc.flag, wantRecs, gotRecs)
			}
		})
	}
}
