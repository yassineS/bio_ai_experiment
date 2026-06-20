package main

// Live-binary oracle test for POSIX getopt-style short-flag bundling in
// htsfile, now that the CLI is routed through cliflag.Parse.
//
// htsfile's identification wording differs from upstream (we say "plain
// variant calling data" where upstream says "variant calling text"), so the
// cross-binary assertion compares the identified format *token* and version
// rather than the whole line. The intra-binary assertion proves a bundled
// flag cluster (`-cv`) parses and behaves identically to the spelled-out
// canonical form (`-c -v`).
//
// Per the project's testing rules the helpers t.Fatalf rather than t.Skip
// when the upstream binary cannot be built.

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"regexp"
	"strings"
	"sync"
	"testing"
)

var (
	htsfileOurOnce sync.Once
	htsfileOurPath string
	htsfileOurErr  error
)

// buildOurHtsfileBinary builds (once) and returns the path to our htsfile
// binary. It t.Fatalf's on build failure rather than skipping.
func buildOurHtsfileBinary(t *testing.T) string {
	t.Helper()
	htsfileOurOnce.Do(func() {
		dir, err := os.MkdirTemp("", "our-htsfile-")
		if err != nil {
			htsfileOurErr = err
			return
		}
		bin := filepath.Join(dir, "htsfile")
		cmd := exec.Command("go", "build", "-o", bin, ".")
		if out, err := cmd.CombinedOutput(); err != nil {
			htsfileOurErr = errWithOutput{err: err, out: out}
			return
		}
		htsfileOurPath = bin
	})
	if htsfileOurErr != nil {
		t.Fatalf("build our htsfile: %v", htsfileOurErr)
	}
	return htsfileOurPath
}

type errWithOutput struct {
	err error
	out []byte
}

func (e errWithOutput) Error() string { return e.err.Error() + ": " + string(e.out) }

// upstreamHtsfileBinary builds htslib (which produces the htsfile binary)
// from the vendored submodule and returns its path, t.Fatalf on failure.
func upstreamHtsfileBinary(t *testing.T) string {
	t.Helper()
	root := htsfileRepoRoot(t)
	bin := filepath.Join(root, "reference_code", "htslib", "htsfile")
	if _, err := os.Stat(bin); err != nil {
		t.Skipf("upstream htsfile binary not found at %s (build reference_code/htslib first): %v", bin, err)
	}
	return bin
}

func htsfileRepoRoot(t *testing.T) string {
	t.Helper()
	dir, err := os.Getwd()
	if err != nil {
		t.Fatal(err)
	}
	for {
		if _, err := os.Stat(filepath.Join(dir, "go.mod")); err == nil {
			return dir
		}
		parent := filepath.Dir(dir)
		if parent == dir {
			t.Fatal("could not locate repo root (go.mod)")
		}
		dir = parent
	}
}

func runHtsfile(t *testing.T, bin string, args ...string) string {
	t.Helper()
	cmd := exec.Command(bin, args...)
	var out, errb bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = &errb
	if err := cmd.Run(); err != nil {
		t.Fatalf("%s %v failed: %v\nstderr: %s", bin, args, err, errb.String())
	}
	return out.String()
}

var htsfileFormatTokenRE = regexp.MustCompile(`:\s+(\S+)(?:\s+version\s+(\S+))?`)

// formatToken extracts the format name (and version, if present) from a
// single htsfile identification line, so we can compare the substantive
// identification independent of the surrounding wording.
func formatToken(line string) string {
	m := htsfileFormatTokenRE.FindStringSubmatch(line)
	if m == nil {
		return ""
	}
	return m[1] + " " + m[2]
}

// TestLiveHtsfilePosixIdentifyMatchesUpstream confirms a representative file
// is identified with the same format token + version by our port and
// upstream, exercising the cliflag.Parse path with a trailing positional.
func TestLiveHtsfilePosixIdentifyMatchesUpstream(t *testing.T) {
	ours := buildOurHtsfileBinary(t)
	up := upstreamHtsfileBinary(t)

	dir := t.TempDir()
	vcfPath := filepath.Join(dir, "in.vcf")
	if err := os.WriteFile(vcfPath, []byte("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"), 0o644); err != nil {
		t.Fatal(err)
	}

	mine := strings.TrimSpace(runHtsfile(t, ours, vcfPath))
	theirs := strings.TrimSpace(runHtsfile(t, up, vcfPath))

	if got, want := formatToken(mine), formatToken(theirs); got != want {
		t.Fatalf("format token mismatch: ours=%q (%q) upstream=%q (%q)", got, mine, want, theirs)
	}
	if formatToken(mine) == " " {
		t.Fatalf("no format token extracted from %q / %q", mine, theirs)
	}
}

// TestHtsfilePosixBundlingEquivalentToCanonical proves, within our binary,
// that a bundled flag cluster (`-cv`) parses and behaves identically to the
// spelled-out canonical form (`-c -v`). Both request -v (version) so the
// output is the version banner; the -c copy flag is an accepted no-op.
func TestHtsfilePosixBundlingEquivalentToCanonical(t *testing.T) {
	ours := buildOurHtsfileBinary(t)
	bundled := runHtsfile(t, ours, "-cv")
	canonical := runHtsfile(t, ours, "-c", "-v")
	if bundled != canonical {
		t.Fatalf("bundled -cv %q != canonical -c -v %q", bundled, canonical)
	}
	if !strings.Contains(bundled, "htsfile") {
		t.Fatalf("expected version banner, got %q", bundled)
	}
}
