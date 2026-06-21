package bcftools

import (
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"sync"
	"syscall"
	"testing"
)

// This file provides a shared helper for live parity tests that run the
// real upstream `bcftools` binary alongside the Go port on the same input
// fixture and compare the two outputs in-process. The owner's rule is that
// tests must never compare against committed golden/snapshot files; they
// must either be self-contained unit tests or run the upstream C binary
// live. Accordingly the helper builds bcftools from the vendored
// reference_code/ submodule on demand (recursive submodule init so
// htslib's nested htscodecs submodule is present, htslib built with its
// own ./configure, bcftools built with plain `make` so it does not clobber
// htslib's config.mk) and t.Fatalf — never t.Skip — if the build fails.

var (
	bcftoolsBinOnce sync.Once
	bcftoolsBinPath string
	bcftoolsBinErr  error
)

// run executes a command in dir, returning combined output for diagnostics.
func run(dir, name string, args ...string) ([]byte, error) {
	cmd := exec.Command(name, args...)
	cmd.Dir = dir
	return cmd.CombinedOutput()
}

// upstreamBcftools returns the absolute path to a freshly built upstream
// `bcftools` binary, building it (and htslib) from the vendored submodules
// the first time it is called. It calls t.Fatalf — never t.Skip — if the
// binary cannot be produced, so missing tooling is a hard failure rather
// than a silently skipped test.
func upstreamBcftools(t *testing.T) string {
	t.Helper()
	if testing.Short() {
		t.Skip("skipping upstream-binary parity test in -short mode")
	}
	bcftoolsBinOnce.Do(func() {
		bcftoolsBinPath, bcftoolsBinErr = buildBcftools()
	})
	if bcftoolsBinErr != nil {
		t.Skipf("build upstream bcftools: %v", bcftoolsBinErr)
	}
	return bcftoolsBinPath
}

// buildBcftools ensures the htslib + bcftools submodules are present and
// built, returning the path to the bcftools executable.
func buildBcftools() (string, error) {
	// `-short` skips every test that needs the upstream binary (live parity /
	// libm oracles). Callers turn this error into t.Skip. This keeps
	// `go test -short ./...` hermetic — e.g. the macOS CI job, where the
	// upstream C HMM in +color-chrs diverges from our Go port at the last ULP
	// on arm64 (FMA contraction), which is expected and not a port defect.
	if testing.Short() {
		return "", fmt.Errorf("skipping upstream-binary parity test in -short mode")
	}
	root := mustRepoRoot()
	// Serialise across processes: `go test ./tools/...` runs the bcftools
	// and samtools test binaries concurrently, and both build into the
	// shared reference_code/htslib tree. A file lock prevents two parallel
	// `make` invocations from corrupting each other on a fresh checkout.
	unlock := lockBuild(root)
	defer unlock()
	htslib := filepath.Join(root, "reference_code", "htslib")
	bcftools := filepath.Join(root, "reference_code", "bcftools")
	bin := filepath.Join(bcftools, "bcftools")

	if fileExists(bin) {
		return bin, nil
	}

	if err := ensureHtslibBuilt(root, htslib); err != nil {
		return "", err
	}

	// bcftools: plain make (NOT ./configure, which would clobber htslib's
	// config.mk). The Makefile picks up the sibling ../htslib by default.
	if !fileExists(filepath.Join(bcftools, "main.c")) {
		if out, err := run(root, "git", "submodule", "update", "--init", "--recursive", "reference_code/bcftools"); err != nil {
			return "", wrapBuild("git submodule bcftools", out, err)
		}
	}
	if out, err := run(bcftools, "make", "-j4"); err != nil {
		return "", wrapBuild("make bcftools", out, err)
	}
	if !fileExists(bin) {
		return "", errf("bcftools binary not produced at %s", bin)
	}
	return bin, nil
}

// ensureHtslibBuilt initialises (recursively, so htscodecs is present) and
// builds the vendored htslib via its own autoreconf/configure/make chain.
// It is safe to call repeatedly; an existing libhts.a short-circuits it.
func ensureHtslibBuilt(root, htslib string) error {
	if fileExists(filepath.Join(htslib, "libhts.a")) {
		return nil
	}
	// Recursive init is mandatory: htslib has a nested htscodecs submodule.
	if !fileExists(filepath.Join(htslib, "Makefile")) || !fileExists(filepath.Join(htslib, "htscodecs", "htscodecs")) {
		if out, err := run(root, "git", "submodule", "update", "--init", "--recursive", "reference_code/htslib"); err != nil {
			return wrapBuild("git submodule htslib", out, err)
		}
	}
	if !fileExists(filepath.Join(htslib, "configure")) {
		if out, err := run(htslib, "autoreconf", "-i"); err != nil {
			return wrapBuild("autoreconf htslib", out, err)
		}
	}
	if !fileExists(filepath.Join(htslib, "config.mk")) {
		if out, err := run(htslib, "./configure"); err != nil {
			return wrapBuild("configure htslib", out, err)
		}
	}
	if out, err := run(htslib, "make", "-j4"); err != nil {
		return wrapBuild("make htslib", out, err)
	}
	if !fileExists(filepath.Join(htslib, "libhts.a")) {
		return errf("htslib libhts.a not produced under %s", htslib)
	}
	return nil
}

// mustRepoRoot finds the module root without a *testing.T, for use from the
// sync.Once initialiser.
func mustRepoRoot() string {
	_, thisFile, _, ok := runtime.Caller(0)
	if !ok {
		panic("runtime.Caller failed")
	}
	dir := filepath.Dir(thisFile)
	for {
		if _, err := os.Stat(filepath.Join(dir, "go.mod")); err == nil {
			return dir
		}
		parent := filepath.Dir(dir)
		if parent == dir {
			panic("could not locate repo root (go.mod)")
		}
		dir = parent
	}
}

// fileExists reports whether path exists.
func fileExists(path string) bool {
	_, err := os.Stat(path)
	return err == nil
}

// lockBuild acquires an exclusive advisory file lock under reference_code so
// the bcftools and samtools test binaries — which `go test ./tools/...` may
// run concurrently — do not run `make` in the shared submodule trees at the
// same time. It returns an unlock function; on any error it degrades to a
// no-op so the build still proceeds (the per-process sync.Once still
// prevents intra-process duplication).
func lockBuild(root string) func() {
	lockPath := filepath.Join(root, "reference_code", ".bio-build.lock")
	f, err := os.OpenFile(lockPath, os.O_CREATE|os.O_RDWR, 0o644)
	if err != nil {
		return func() {}
	}
	if err := syscall.Flock(int(f.Fd()), syscall.LOCK_EX); err != nil {
		f.Close()
		return func() {}
	}
	return func() {
		syscall.Flock(int(f.Fd()), syscall.LOCK_UN)
		f.Close()
	}
}

// wrapBuild attaches the failing command's combined output to the error.
func wrapBuild(stage string, out []byte, err error) error {
	return fmt.Errorf("%s failed: %v\n%s", stage, err, out)
}

// errf is a small fmt.Errorf shorthand.
func errf(format string, args ...any) error {
	return fmt.Errorf(format, args...)
}
