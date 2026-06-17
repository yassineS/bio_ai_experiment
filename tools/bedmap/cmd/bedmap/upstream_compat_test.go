package main

// CLI-level retro-compatibility tests for bedmap against the live upstream
// bedtools binary. Upstream `bedtools map` uses multi-char single-dash flag
// NAMES (-c, -o, -delim, -s, -f, ...), parsed natively by Go's `flag` package
// — no POSIX bundling, no cliflag.Parse. bedmap implements the column-mapping
// semantics with matching output for -c/-o aggregation, so these assert
// byte-for-byte equality with the live upstream binary. t.Fatalf, never
// t.Skip.

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"sync"
	"testing"
)

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

func buildOurs(t *testing.T) string {
	t.Helper()
	root, err := repoRoot()
	if err != nil {
		t.Fatalf("locate repo root: %v", err)
	}
	bin := filepath.Join(t.TempDir(), "bedmap")
	cmd := exec.Command("go", "build", "-o", bin,
		"github.com/yassineS/bio_ai_experiment/tools/bedmap/cmd/bedmap")
	cmd.Dir = root
	if out, buildErr := cmd.CombinedOutput(); buildErr != nil {
		t.Fatalf("build bedmap: %v\n%s", buildErr, out)
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

// TestUpstreamCompat_Map_MultiCharSingleDash drives both binaries with the
// upstream multi-char single-dash command lines for `bedtools map` and asserts
// our output equals the live upstream binary's.
func TestUpstreamCompat_Map_MultiCharSingleDash(t *testing.T) {
	bt := upstreamBedtools(t)
	ours := buildOurs(t)
	dir := t.TempDir()
	a := writeFile(t, dir, "a.bed", "chr1\t10\t20\tg1\nchr1\t50\t60\tg2\n")
	b := writeFile(t, dir, "b.bed", "chr1\t12\t18\tx\t5\nchr1\t15\t25\ty\t3\n")

	cases := []struct {
		name string
		args []string
	}{
		{"sum", []string{"-c", "5", "-o", "sum"}},
		{"mean", []string{"-c", "5", "-o", "mean"}},
		{"count", []string{"-c", "5", "-o", "count"}},
		{"collapse", []string{"-c", "4", "-o", "collapse"}},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			want := runCapture(t, bt, append([]string{"map", "-a", a, "-b", b}, tc.args...)...)
			got := runCapture(t, ours, append([]string{"-a", a, "-b", b}, tc.args...)...)
			if !bytes.Equal(got, want) {
				t.Fatalf("map mismatch for %v\nupstream:\n%s\nours:\n%s", tc.args, want, got)
			}
		})
	}
}

// TestUpstreamCompat_Map_TieBreakCollapseOrder asserts byte-for-byte parity
// with the live upstream binary for the order-sensitive ops when B has many
// equal-(chrom, start) records whose ends are NOT in ascending order. The
// collapsed/distinct values must follow upstream's stream (input) order, not
// chromEnd order. Order-independent ops (sum/count) are included as a guard
// against regressing the already-matching paths.
func TestUpstreamCompat_Map_TieBreakCollapseOrder(t *testing.T) {
	bt := upstreamBedtools(t)
	ours := buildOurs(t)
	dir := t.TempDir()

	a := writeFile(t, dir, "a.bed",
		"chr1\t0\t200\tA1\t0\t+\n"+
			"chr2\t0\t100\tA2\t0\t-\n")
	// Equal (chr1,10) B records, ends out of order, mixed strands; plus a second
	// chromosome with equal-start records.
	b := writeFile(t, dir, "b.bed",
		"chr1\t10\t100\ta\t1\t+\n"+
			"chr1\t10\t50\tb\t2\t-\n"+
			"chr1\t10\t75\tc\t3\t+\n"+
			"chr1\t10\t60\td\t4\t-\n"+
			"chr1\t20\t30\te\t5\t+\n"+
			"chr2\t0\t50\tf\t6\t+\n"+
			"chr2\t0\t40\tg\t7\t-\n")

	argSets := [][]string{
		{"-c", "4", "-o", "collapse"},
		{"-c", "4", "-o", "distinct"},
		{"-c", "5", "-o", "collapse"},
		{"-c", "4,5", "-o", "collapse,collapse"},
		{"-s", "-c", "4", "-o", "collapse"},
		{"-S", "-c", "4", "-o", "collapse"},
		{"-c", "5", "-o", "sum"},   // order-independent: must stay matching
		{"-c", "4", "-o", "count"}, // order-independent: must stay matching
	}
	for _, args := range argSets {
		args := args
		t.Run(joinArgs(args), func(t *testing.T) {
			want := runCapture(t, bt, append([]string{"map", "-a", a, "-b", b}, args...)...)
			got := runCapture(t, ours, append([]string{"-a", a, "-b", b}, args...)...)
			if !bytes.Equal(got, want) {
				t.Fatalf("map collapse-order mismatch for %v\nupstream:\n%s\nours:\n%s", args, want, got)
			}
		})
	}
}

func joinArgs(args []string) string {
	out := ""
	for i, a := range args {
		if i > 0 {
			out += "_"
		}
		out += a
	}
	return out
}

// TestUpstreamCompat_Map_CustomDelim verifies the multi-char `-delim` flag is
// accepted (single dash) and applied identically to upstream.
func TestUpstreamCompat_Map_CustomDelim(t *testing.T) {
	bt := upstreamBedtools(t)
	ours := buildOurs(t)
	dir := t.TempDir()
	a := writeFile(t, dir, "a.bed", "chr1\t10\t20\tg1\n")
	b := writeFile(t, dir, "b.bed", "chr1\t12\t18\tx\nchr1\t15\t19\ty\n")
	want := runCapture(t, bt, "map", "-a", a, "-b", b, "-c", "4", "-o", "collapse", "-delim", "|")
	got := runCapture(t, ours, "-a", a, "-b", b, "-c", "4", "-o", "collapse", "-delim", "|")
	if !bytes.Equal(got, want) {
		t.Fatalf("map -delim mismatch\nupstream:\n%s\nours:\n%s", want, got)
	}
}
