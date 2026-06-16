package bedgroupby

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
// the output byte for byte against this port. It closes the documented
// scalar-numeric formatting gap for `bedtools groupby` (groupby.t20): a `min`
// over a large integer-like token must be rendered with std::setprecision(10),
// i.e. 7777788888899999 -> 7.777788889e+15, not preserved as a raw integer.
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

// sharedReferenceBedtools locates a pre-built upstream `bedtools` in the main
// (non-worktree) checkout. Git worktrees store their gitdir under the main
// repo's .git/worktrees/<name>; the path before "/.git/" is that main checkout,
// whose reference_code/bedtools submodule is the populated one. Returns "" if no
// such binary can be found.
func sharedReferenceBedtools() string {
	var gd []byte
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
		// In a git worktree the bedtools submodule may not be checked out
		// locally; the main checkout shares the same submodule. Reuse it.
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

// runUpstreamGroupby pipes stdin into `bedtools groupby <args>` and returns its
// stdout, failing hard on any error.
func runUpstreamGroupby(t *testing.T, stdin string, args ...string) []byte {
	t.Helper()
	bin := upstreamBedtools(t)
	cmd := exec.Command(bin, append([]string{"groupby"}, args...)...)
	cmd.Stdin = strings.NewReader(stdin)
	var stdout, stderr bytes.Buffer
	cmd.Stdout = &stdout
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream groupby %v: %v\nstderr: %s", args, err, stderr.String())
	}
	return stdout.Bytes()
}

// runPortGroupby runs this port's Group over the same stdin/options.
func runPortGroupby(t *testing.T, stdin string, opts Options) []byte {
	t.Helper()
	var out bytes.Buffer
	if _, err := Group(strings.NewReader(stdin), &out, opts); err != nil {
		t.Fatalf("port Group: %v", err)
	}
	return out.Bytes()
}

// TestUpstreamParity_Groupby_T20_LargeIntMin mirrors groupby.t20:
//
//	echo "a\t1253555555355577777777\t7777788888899999" |
//	  bedtools groupby -i - -g 1 -c 2,3 -o distinct,min
//
// The `distinct` on the 22-digit token must be preserved verbatim, while the
// `min` (a scalar numeric op) is rendered through setprecision(10), collapsing
// 7777788888899999 to 7.777788889e+15. Asserts byte-for-byte vs upstream.
func TestUpstreamParity_Groupby_T20_LargeIntMin(t *testing.T) {
	const input = "a\t1253555555355577777777\t7777788888899999\n"
	want := runUpstreamGroupby(t, input, "-i", "-", "-g", "1", "-c", "2,3", "-o", "distinct,min")
	got := runPortGroupby(t, input, Options{
		GroupCols: []int{1},
		AggCols:   []int{2, 3},
		Ops:       []string{"distinct", "min"},
	})
	if !bytes.Equal(got, want) {
		t.Fatalf("groupby.t20 parity mismatch.\nupstream:\n%q\nport:\n%q", want, got)
	}
}

// TestUpstreamParity_Groupby_ScalarNumericPrecision exercises a spread of
// scalar numeric ops on values that do and do not exceed precision-10
// significant digits, asserting this port matches upstream's setprecision(10)
// rendering byte-for-byte (small/exact values stay plain; oversized integers go
// scientific).
func TestUpstreamParity_Groupby_ScalarNumericPrecision(t *testing.T) {
	const input = "g\t12345678901234\t12345678901234\n" +
		"g\t99999999999999\t1\n"
	for _, op := range []string{"min", "max", "sum", "mean", "median", "absmin", "absmax"} {
		op := op
		t.Run(op, func(t *testing.T) {
			want := runUpstreamGroupby(t, input, "-i", "-", "-g", "1", "-c", "2", "-o", op)
			got := runPortGroupby(t, input, Options{
				GroupCols: []int{1},
				AggCols:   []int{2},
				Ops:       []string{op},
			})
			if !bytes.Equal(got, want) {
				t.Fatalf("op %s parity mismatch.\nupstream:\n%q\nport:\n%q", op, want, got)
			}
		})
	}
}
