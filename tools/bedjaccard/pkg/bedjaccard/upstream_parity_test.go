package bedjaccard

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"strings"
	"sync"
	"testing"
)

// This file holds live-upstream parity tests: it builds (or reuses) the real
// upstream `bedtools` binary from the vendored submodule, runs it, and compares
// the output byte for byte against this port. It closes the documented BAM
// input + `-bed` gap for `bedtools jaccard`.
//
// Upstream `jaccard` only emits the interval-metric summary for BAM input when
// `-bed` is supplied (without it, jaccard inherits BAM *output* mode and writes
// a garbage BAM blob to stdout). This port always renders BAM alignments as
// BED12 intervals, so `-bed` is accepted as a documented no-op; the test proves
// our default output equals upstream's `-bed` output byte for byte.
//
// The tests t.Fatalf (never t.Skip): a missing/unbuildable submodule is a hard
// failure, matching the parity-rig policy.

var (
	upstreamBedtoolsOnce sync.Once
	upstreamBedtoolsPath string
	upstreamBedtoolsErr  error
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
		// When running inside a git worktree the bedtools submodule may not be
		// checked out locally; the build's main checkout shares the same
		// submodule. Reuse it if present rather than failing.
		if _, err := os.Stat(filepath.Join(dir, "Makefile")); err != nil {
			if shared := sharedReferenceBedtools(); shared != "" {
				upstreamBedtoolsPath = shared
				return
			}
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

// sharedReferenceBedtools locates a pre-built upstream `bedtools` in the main
// (non-worktree) checkout. Git worktrees store their gitdir under the main
// repo's .git/worktrees/<name>; the path before "/.git/" is that main checkout,
// whose reference_code/bedtools submodule is the populated one. Returns "" if no
// such binary can be found.
func sharedReferenceBedtools() string {
	gd, err := os.ReadFile(".git")
	// In a worktree, the package dir's ancestor holds a `.git` *file*; but tests
	// run from the package dir, so walk up looking for it.
	if err != nil {
		dir, _ := os.Getwd()
		for {
			b, e := os.ReadFile(filepath.Join(dir, ".git"))
			if e == nil {
				gd = b
				break
			}
			parent := filepath.Dir(dir)
			if parent == dir {
				return ""
			}
			dir = parent
		}
	}
	line := strings.TrimSpace(string(gd))
	const prefix = "gitdir: "
	if !strings.HasPrefix(line, prefix) {
		return ""
	}
	gitdir := strings.TrimPrefix(line, prefix)
	idx := strings.Index(gitdir, "/.git/")
	if idx < 0 {
		return ""
	}
	mainRoot := gitdir[:idx]
	bin := filepath.Join(mainRoot, "reference_code", "bedtools", "bin", "bedtools")
	if _, err := os.Stat(bin); err == nil {
		return bin
	}
	return ""
}

// fixturePath returns the absolute path to a testdata/parity fixture.
func fixturePath(t *testing.T, name string) string {
	t.Helper()
	p, err := filepath.Abs(filepath.Join("..", "..", "testdata", "parity", name))
	if err != nil {
		t.Fatalf("abs path for %s: %v", name, err)
	}
	return p
}

// TestUpstreamParity_Jaccard_BAM_Bed mirrors jaccard.t09:
//
//	bedtools jaccard -a a.bam -b three_blocks_match.bam -bed
//
// It runs the real upstream binary with `-bed` and asserts this port's BAM
// path (which always treats BAM alignments as BED12 intervals) produces the
// identical `intersection\tunion\tjaccard\tn_intersections` table.
func TestUpstreamParity_Jaccard_BAM_Bed(t *testing.T) {
	bin := upstreamBedtools(t)
	aBam := fixturePath(t, "a.bam")
	bBam := fixturePath(t, "three_blocks_match.bam")

	cmd := exec.Command(bin, "jaccard", "-a", aBam, "-b", bBam, "-bed")
	var upStdout, upStderr bytes.Buffer
	cmd.Stdout = &upStdout
	cmd.Stderr = &upStderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream jaccard -bed: %v\nstderr: %s", err, upStderr.String())
	}
	want := upStdout.Bytes()

	// This port: BAMToBED mirrors the CLI `-bed` flag (a documented no-op since
	// BAM is always rendered as BED12 intervals).
	aData, err := os.ReadFile(aBam)
	if err != nil {
		t.Fatalf("read a.bam: %v", err)
	}
	bData, err := os.ReadFile(bBam)
	if err != nil {
		t.Fatalf("read three_blocks_match.bam: %v", err)
	}
	var got bytes.Buffer
	if _, err := Run(bytes.NewReader(aData), bytes.NewReader(bData), &got, Options{BAMToBED: true}); err != nil {
		t.Fatalf("port Run on BAM input: %v", err)
	}

	if !bytes.Equal(got.Bytes(), want) {
		t.Fatalf("BAM -bed parity mismatch.\nupstream:\n%s\nport:\n%s", want, got.Bytes())
	}
}
