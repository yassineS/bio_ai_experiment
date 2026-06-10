package main

// CLI-level retro-compatibility tests for bedintersect.
//
// These verify two things at once:
//
//  1. The *upstream* bedtools command line — which uses multi-character,
//     single-dash flag NAMES (-wa, -wb, -loj, ...) rather than POSIX
//     single-character bundling — parses correctly through Go's `flag`
//     package without any "flag provided but not defined" error. Go's flag
//     parser natively treats `-wa` as the flag named "wa", so these forms
//     must keep working WITHOUT routing through cliflag.Parse (which would
//     wrongly expand `-wa` into `-w -a`).
//
//  2. For the option subset bedintersect implements with the same semantics
//     as upstream, our output byte-for-byte matches the live upstream
//     `bedtools` binary built from reference_code/bedtools.
//
// The upstream binary is built once via the uniquely-named upstreamBedtools
// helper. Tests use t.Fatalf (never t.Skip) so a missing/un-buildable
// upstream binary is a hard failure that surfaces the regression.

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"sync"
	"testing"
)

// upstreamBedtools builds (once) and returns the path to the live upstream
// `bedtools` binary from the vendored reference_code submodule. It t.Fatalf's
// if the submodule is not checked out or the build fails.
var (
	upstreamBedtoolsOnce sync.Once
	upstreamBedtoolsPath string
	upstreamBedtoolsErr  error
)

func upstreamBedtools(t *testing.T) string {
	t.Helper()
	upstreamBedtoolsOnce.Do(func() {
		root, err := repoRoot()
		if err != nil {
			upstreamBedtoolsErr = err
			return
		}
		dir := filepath.Join(root, "reference_code", "bedtools")
		if _, statErr := os.Stat(filepath.Join(dir, "Makefile")); statErr != nil {
			upstreamBedtoolsErr = statErr
			return
		}
		bin := filepath.Join(dir, "bin", "bedtools")
		if _, statErr := os.Stat(bin); statErr != nil {
			cmd := exec.Command("make", "-j", "4")
			cmd.Dir = dir
			if out, buildErr := cmd.CombinedOutput(); buildErr != nil {
				upstreamBedtoolsErr = &buildError{buildErr, out}
				return
			}
		}
		upstreamBedtoolsPath = bin
	})
	if upstreamBedtoolsErr != nil {
		t.Fatalf("upstream bedtools unavailable: %v\n"+
			"run: git submodule update --init reference_code/bedtools && "+
			"(cd reference_code/bedtools && make -j\"$(nproc)\")", upstreamBedtoolsErr)
	}
	return upstreamBedtoolsPath
}

type buildError struct {
	err error
	out []byte
}

func (e *buildError) Error() string { return e.err.Error() + "\n" + string(e.out) }

// repoRoot walks up from this test file to the module root (the dir holding
// go.mod).
func repoRoot() (string, error) {
	_, file, _, _ := runtime.Caller(0)
	dir := filepath.Dir(file)
	for {
		if _, err := os.Stat(filepath.Join(dir, "go.mod")); err == nil {
			return dir, nil
		}
		parent := filepath.Dir(dir)
		if parent == dir {
			return "", os.ErrNotExist
		}
		dir = parent
	}
}

// buildOurs compiles the bedintersect binary under test into a temp dir.
func buildOurs(t *testing.T) string {
	t.Helper()
	root, err := repoRoot()
	if err != nil {
		t.Fatalf("locate repo root: %v", err)
	}
	bin := filepath.Join(t.TempDir(), "bedintersect")
	cmd := exec.Command("go", "build", "-o", bin,
		"github.com/yassineS/bio_ai_experiment/tools/bedintersect/cmd/bedintersect")
	cmd.Dir = root
	if out, buildErr := cmd.CombinedOutput(); buildErr != nil {
		t.Fatalf("build bedintersect: %v\n%s", buildErr, out)
	}
	return bin
}

func writeFile(t *testing.T, dir, name, content string) string {
	t.Helper()
	p := filepath.Join(dir, name)
	if err := os.WriteFile(p, []byte(content), 0o644); err != nil {
		t.Fatalf("write %s: %v", name, err)
	}
	return p
}

func runCapture(t *testing.T, name string, args ...string) []byte {
	t.Helper()
	cmd := exec.Command(name, args...)
	var stdout, stderr bytes.Buffer
	cmd.Stdout = &stdout
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("run %s %v failed: %v\nstderr:\n%s", name, args, err, stderr.String())
	}
	return stdout.Bytes()
}

// TestUpstreamCompat_Intersect_MultiCharSingleDash exercises the upstream
// multi-char single-dash flag NAMES against both binaries on 3-column BED
// (the column layout where bedintersect and upstream converge), asserting our
// output equals the live upstream binary's.
func TestUpstreamCompat_Intersect_MultiCharSingleDash(t *testing.T) {
	bt := upstreamBedtools(t)
	ours := buildOurs(t)
	dir := t.TempDir()
	a := writeFile(t, dir, "a.bed", "chr1\t10\t20\nchr1\t30\t40\nchr1\t60\t70\n")
	b := writeFile(t, dir, "b.bed", "chr1\t15\t25\nchr1\t35\t45\n")

	cases := []struct {
		name string
		args []string // flags only; -a/-b prepended
	}{
		{"default", nil},
		{"wa", []string{"-wa"}},
		{"count", []string{"-c"}},
		{"invert", []string{"-v"}},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			args := append([]string{"intersect", "-a", a, "-b", b}, tc.args...)
			want := runCapture(t, bt, args...)
			ourArgs := append([]string{"-a", a, "-b", b}, tc.args...)
			got := runCapture(t, ours, ourArgs...)
			if !bytes.Equal(got, want) {
				t.Fatalf("output mismatch for %v\nupstream:\n%s\nours:\n%s", tc.args, want, got)
			}
		})
	}
}

// TestUpstreamCompat_Intersect_DoubleDashAlsoWorks confirms the GNU-style
// double-dash form of the same multi-char names parses too (Go flag treats
// -wa and --wa identically). This is what guards against a future
// well-meaning switch to POSIX bundling that would break -wa.
func TestUpstreamCompat_Intersect_DoubleDashAlsoWorks(t *testing.T) {
	ours := buildOurs(t)
	dir := t.TempDir()
	a := writeFile(t, dir, "a.bed", "chr1\t10\t20\n")
	b := writeFile(t, dir, "b.bed", "chr1\t15\t25\n")
	single := runCapture(t, ours, "-a", a, "-b", b, "-wa")
	double := runCapture(t, ours, "-a", a, "-b", b, "--wa")
	if !bytes.Equal(single, double) {
		t.Fatalf("-wa and --wa diverged:\n-wa:\n%s\n--wa:\n%s", single, double)
	}
}

// TestUpstreamCompat_Intersect_StdinDash verifies `-` routes A from stdin,
// matching upstream's stdin convention.
func TestUpstreamCompat_Intersect_StdinDash(t *testing.T) {
	ours := buildOurs(t)
	dir := t.TempDir()
	b := writeFile(t, dir, "b.bed", "chr1\t15\t25\n")
	cmd := exec.Command(ours, "-a", "-", "-b", b, "-wa")
	cmd.Stdin = bytes.NewReader([]byte("chr1\t10\t20\n"))
	var out, errb bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = &errb
	if err := cmd.Run(); err != nil {
		t.Fatalf("stdin run failed: %v\nstderr:\n%s", err, errb.String())
	}
	if got := out.String(); got != "chr1\t10\t20\n" {
		t.Fatalf("stdin -a mismatch: got %q", got)
	}
}
