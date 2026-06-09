package samtools

import (
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"strings"
	"sync"
	"testing"
	"time"
)

// repoRootForTest returns the absolute path to the repository root,
// derived from this test file's own location (tools/samtools/pkg/samtools
// is four directories below the root). Deriving it from the source path
// rather than the working directory keeps the helper robust regardless of
// where `go test` is invoked from.
func repoRootForTest(t *testing.T) string {
	t.Helper()
	_, thisFile, _, ok := runtime.Caller(0)
	if !ok {
		t.Fatalf("runtime.Caller failed; cannot locate repo root")
	}
	// thisFile == <root>/tools/samtools/pkg/samtools/consensus_upstream_test.go
	root := filepath.Clean(filepath.Join(filepath.Dir(thisFile), "..", "..", "..", ".."))
	if _, err := os.Stat(filepath.Join(root, "go.mod")); err != nil {
		t.Fatalf("repo root %q has no go.mod: %v", root, err)
	}
	return root
}

// upstreamSamtoolsOnce memoises the (one-time, expensive) submodule
// init + htslib/samtools build across all tests in the package so that
// multiple upstream-parity tests share a single built binary.
var upstreamSamtoolsOnce struct {
	sync.Once
	path string
	err  error
}

// upstreamSamtoolsBinary locates, or builds from the vendored
// submodules, the upstream samtools binary and returns its absolute
// path. The build is performed at most once per test binary. It never
// calls t.Skip: per the project's testing rules the upstream parity
// check must actually execute, so a genuine inability to produce the
// binary is a hard failure (t.Fatalf), not a skip.
//
// Build steps mirror the docs and the CI upstream-parity job:
//  1. `git submodule update --init --recursive reference_code/htslib
//     reference_code/samtools` (only if the sources are missing), with
//     exponential backoff retries for transient network failures. The
//     --recursive flag is required because htslib carries its own nested
//     htscodecs submodule, without which its ./configure aborts.
//  2. `autoreconf -i && ./configure && make -j` in reference_code/htslib.
//  3. `make -j` in reference_code/samtools (its own ./configure is
//     skipped: it is optional and can clobber htslib's config.mk).
//
// An already-built reference_code/samtools/samtools is reused as-is.
func upstreamSamtoolsBinary(t *testing.T) string {
	t.Helper()
	root := repoRootForTest(t)
	upstreamSamtoolsOnce.Do(func() {
		upstreamSamtoolsOnce.path, upstreamSamtoolsOnce.err = buildUpstreamSamtools(root)
	})
	if upstreamSamtoolsOnce.err != nil {
		t.Fatalf("could not obtain upstream samtools binary: %v", upstreamSamtoolsOnce.err)
	}
	return upstreamSamtoolsOnce.path
}

// buildUpstreamSamtools performs the locate-or-build work for
// upstreamSamtoolsBinary, returning the binary path or an error.
func buildUpstreamSamtools(root string) (string, error) {
	samtoolsDir := filepath.Join(root, "reference_code", "samtools")
	htslibDir := filepath.Join(root, "reference_code", "htslib")
	bin := filepath.Join(samtoolsDir, "samtools")

	// Fast path: reuse an already-built binary.
	if fi, err := os.Stat(bin); err == nil && !fi.IsDir() {
		return bin, nil
	}

	// Ensure the submodule sources are present. --recursive pulls
	// htslib's nested htscodecs submodule, which its ./configure needs.
	htscodecsC := filepath.Join(htslibDir, "htscodecs", "htscodecs", "htscodecs.c")
	if _, err := os.Stat(filepath.Join(samtoolsDir, "bam_consensus.c")); err != nil {
		if err := runWithRetry(root, "git", "submodule", "update", "--init", "--recursive",
			"reference_code/htslib", "reference_code/samtools"); err != nil {
			return "", err
		}
	} else if _, err := os.Stat(htscodecsC); err != nil {
		// Sources present but the nested htscodecs submodule is missing
		// (e.g. a non-recursive prior init): pull it now.
		if err := runWithRetry(htslibDir, "git", "submodule", "update", "--init", "--recursive"); err != nil {
			return "", err
		}
	}

	// Build htslib: autoreconf -i (a git checkout ships no ./configure)
	// then ./configure && make. We always re-run ./configure so a stale
	// or error-stamped config.mk (samtools' own configure can overwrite
	// htslib's during bundled-htslib detection) is regenerated cleanly.
	if _, err := os.Stat(filepath.Join(htslibDir, "configure")); err != nil {
		if err := runCmd(htslibDir, "autoreconf", "-i"); err != nil {
			return "", err
		}
	}
	if err := runCmd(htslibDir, "./configure"); err != nil {
		return "", err
	}
	if err := runCmd(htslibDir, "make", "-j"); err != nil {
		return "", err
	}

	// Build samtools with a plain `make`. Per samtools' INSTALL, running
	// its own ./configure is optional and only diagnoses build problems;
	// the Makefile builds against the bundled ../htslib we just built.
	// We deliberately skip samtools' ./configure because it regenerates
	// (and can corrupt) htslib's config.mk.
	if err := runCmd(samtoolsDir, "make", "-j"); err != nil {
		return "", err
	}

	if fi, err := os.Stat(bin); err != nil || fi.IsDir() {
		return "", err
	}
	return bin, nil
}

// runCmd runs name+args in dir, returning a descriptive error (with
// captured combined output) on failure.
func runCmd(dir, name string, args ...string) error {
	cmd := exec.Command(name, args...)
	cmd.Dir = dir
	out, err := cmd.CombinedOutput()
	if err != nil {
		return &cmdError{name: name, args: args, dir: dir, out: out, err: err}
	}
	return nil
}

// runWithRetry runs name+args in dir, retrying up to four times with
// exponential backoff (2/4/8/16s) to absorb transient network failures
// (e.g. a flaky `git submodule` clone).
func runWithRetry(dir, name string, args ...string) error {
	backoff := []time.Duration{2 * time.Second, 4 * time.Second, 8 * time.Second, 16 * time.Second}
	var lastErr error
	for attempt := 0; attempt <= len(backoff); attempt++ {
		if attempt > 0 {
			time.Sleep(backoff[attempt-1])
		}
		if err := runCmd(dir, name, args...); err != nil {
			lastErr = err
			continue
		}
		return nil
	}
	return lastErr
}

// runUpstreamConsensus invokes the upstream samtools `consensus`
// subcommand on samPath and returns its stdout. samtools writes
// progress/warnings to stderr, which is discarded.
func runUpstreamConsensus(t *testing.T, bin, samPath string, args ...string) string {
	t.Helper()
	full := append([]string{"consensus"}, args...)
	full = append(full, samPath)
	cmd := exec.Command(bin, full...)
	out, err := cmd.Output()
	if err != nil {
		t.Fatalf("upstream samtools %v failed: %v", full, err)
	}
	return string(out)
}

// cmdError is a build/command failure carrying the captured output so
// the test log shows exactly what went wrong.
type cmdError struct {
	name string
	args []string
	dir  string
	out  []byte
	err  error
}

func (e *cmdError) Error() string {
	return strings.TrimSpace(
		e.name + " " + strings.Join(e.args, " ") + " (in " + e.dir + ") failed: " +
			e.err.Error() + "\n--- output ---\n" + string(e.out))
}

// upstreamParitySAM is a small, indel-free fixture spelling ACGTA at
// chr1:1-5 with three identical high-quality reads. Simple frequency
// mode makes the call deterministic and identical between the Go port
// and the upstream binary, so the comparison below is exact.
const upstreamParitySAM = `@HD	VN:1.6
@SQ	SN:chr1	LN:5
r1	0	chr1	1	60	5M	*	0	0	ACGTA	IIIII
r2	0	chr1	1	60	5M	*	0	0	ACGTA	IIIII
r3	0	chr1	1	60	5M	*	0	0	ACGTA	IIIII
`

// TestConsensus_HetOnlyUpstreamParity is the LIVE upstream-parity check
// for the --het-only no-op. It builds (or reuses) the vendored upstream
// samtools binary and proves, by execution, the two facts that justify
// our implementation:
//
//  1. Upstream ignores --het-only: its consensus output is byte-for-byte
//     identical with and without the flag (the C source parses it into
//     opts.het_only but never reads it). This is the 1:1 parity rationale
//     for accepting the flag as a no-op rather than erroring.
//  2. Our Go --het-only output matches upstream's output exactly, in both
//     FASTA and pileup modes, in simple calling mode.
//
// No committed golden/snapshot file is used: both upstream and Go outputs
// are produced live, in-process, and compared directly. If the upstream
// binary cannot be built, the test fails hard (it does not skip).
func TestConsensus_HetOnlyUpstreamParity(t *testing.T) {
	if testing.Short() {
		t.Skip("skipping live upstream build/parity test in -short mode")
	}
	bin := upstreamSamtoolsBinary(t)

	// Write the fixture to a tmp SAM file for the upstream binary.
	samPath := filepath.Join(t.TempDir(), "parity.sam")
	if err := os.WriteFile(samPath, []byte(upstreamParitySAM), 0o600); err != nil {
		t.Fatalf("write fixture: %v", err)
	}

	cases := []struct {
		name      string
		cliFormat string
		goFormat  ConsensusFormat
	}{
		{"fasta", "fasta", ConsensusFASTA},
		{"pileup", "pileup", ConsensusPileup},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			// Fact 1: upstream is a no-op w.r.t. --het-only.
			upNoFlag := runUpstreamConsensus(t, bin, samPath,
				"-m", "simple", "-f", c.cliFormat)
			upWithFlag := runUpstreamConsensus(t, bin, samPath,
				"-m", "simple", "-f", c.cliFormat, "--het-only")
			if upNoFlag != upWithFlag {
				t.Fatalf("upstream --het-only changed %s output (expected no-op):\n"+
					"--- without ---\n%s\n--- with ---\n%s", c.name, upNoFlag, upWithFlag)
			}

			// Fact 2: Go --het-only matches upstream exactly.
			goOpts := ConsensusOptions{Format: c.goFormat, HetOnly: true}
			goOut := runConsensusOnSAM(t, upstreamParitySAM, goOpts)
			if goOut != upWithFlag {
				t.Fatalf("Go --het-only %s output differs from upstream:\n"+
					"--- go ---\n%s\n--- upstream ---\n%s", c.name, goOut, upWithFlag)
			}
		})
	}
}
