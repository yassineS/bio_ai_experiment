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

// upstreamBugSAM mixes heterozygous and homozygous positions: four reads
// over chr1:1-3 give a het flank (A/C at pos1 and pos3) and a homozygous
// middle (all-G at pos2). Because it contains BOTH kinds of position, any
// genuine --het-only filter (like ours) visibly changes the output by
// dropping the homozygous column — which is exactly what lets the test
// below distinguish "filter applied" from "flag ignored".
const upstreamBugSAM = `@HD	VN:1.6
@SQ	SN:chr1	LN:3
r1	0	chr1	1	60	3M	*	0	0	AGA	III
r2	0	chr1	1	60	3M	*	0	0	AGA	III
r3	0	chr1	1	60	3M	*	0	0	CGC	III
r4	0	chr1	1	60	3M	*	0	0	CGC	III
`

// TestConsensus_HetOnlyUpstreamBug is the LIVE check that documents the
// upstream dead-option bug and our intentional, correct divergence from
// it. It builds (or reuses) the vendored upstream samtools binary and
// proves two facts by execution:
//
//  1. UPSTREAM IGNORES --het-only (the bug). samtools parses --het-only
//     into opts.het_only (bam_consensus.c) but never reads that variable
//     anywhere in the calling/output path, so its consensus output is
//     byte-for-byte IDENTICAL with and without the flag — even though the
//     fixture contains a homozygous position the flag's name says should
//     be suppressed. Documented in docs/UPSTREAM_BUGS.md.
//  2. OUR Go --het-only DIFFERS from our own output without the flag: it
//     suppresses the homozygous column (rendered 'N' in FASTA, omitted in
//     pileup), which is the behaviour the flag implies. This is the
//     intentional divergence — we fix upstream's dead option.
//
// No committed golden/snapshot file is used: upstream and Go outputs are
// produced live and compared directly. If the upstream binary cannot be
// built, the test fails hard (it never skips), per the project rules.
func TestConsensus_HetOnlyUpstreamBug(t *testing.T) {
	if testing.Short() {
		t.Skip("skipping live upstream build/parity test in -short mode")
	}
	bin := upstreamSamtoolsBinary(t)

	// Write the fixture to a tmp SAM file for the upstream binary.
	samPath := filepath.Join(t.TempDir(), "bug.sam")
	if err := os.WriteFile(samPath, []byte(upstreamBugSAM), 0o600); err != nil {
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
			// Fact 1: upstream --het-only is a no-op (the bug). Output
			// is identical with and without the flag, despite the
			// fixture's homozygous pos2 that a real filter would drop.
			// We run with -A/--ambig so heterozygous positions are
			// actually emitted as IUPAC codes; that makes the homozygous
			// vs het distinction visible and guarantees a genuine filter
			// would change the output.
			upNoFlag := runUpstreamConsensus(t, bin, samPath,
				"-m", "simple", "-A", "-f", c.cliFormat)
			upWithFlag := runUpstreamConsensus(t, bin, samPath,
				"-m", "simple", "-A", "-f", c.cliFormat, "--het-only")
			if upNoFlag != upWithFlag {
				t.Fatalf("upstream --het-only changed %s output, expected the "+
					"documented no-op (bug):\n--- without ---\n%s\n--- with ---\n%s",
					c.name, upNoFlag, upWithFlag)
			}

			// Fact 2: OUR --het-only differs from our own no-flag output
			// — we suppress the homozygous column. This is the
			// intentional divergence that fixes the upstream dead option.
			goNoFlag := runConsensusOnSAM(t, upstreamBugSAM,
				ConsensusOptions{Format: c.goFormat, AmbigCodes: true})
			goWithFlag := runConsensusOnSAM(t, upstreamBugSAM,
				ConsensusOptions{Format: c.goFormat, AmbigCodes: true, HetOnly: true})
			if goNoFlag == goWithFlag {
				t.Fatalf("Go --het-only %s output did not change (expected the "+
					"homozygous column to be suppressed):\n--- without ---\n%s",
					c.name, goNoFlag)
			}
			// And our no-flag output must still match upstream's no-flag
			// output (we only diverge WHEN the flag is set), confirming we
			// are not accidentally changing the baseline consensus.
			if goNoFlag != upNoFlag {
				t.Fatalf("Go baseline (no --het-only) %s output differs from "+
					"upstream:\n--- go ---\n%s\n--- upstream ---\n%s",
					c.name, goNoFlag, upNoFlag)
			}
		})
	}
}
