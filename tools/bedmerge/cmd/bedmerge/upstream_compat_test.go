package main

// CLI-level retro-compatibility tests for bedmerge against the live upstream
// bedtools binary. See the bedintersect equivalent for the full rationale:
// upstream uses multi-char single-dash flag NAMES (-c, -o, -d, -s, -S), which
// Go's `flag` package parses natively — no POSIX bundling, no cliflag.Parse.
//
// bedmerge implements `bedtools merge` with matching semantics for -i, -d,
// -c/-o aggregation, and -s, so these assert byte-for-byte equality with the
// live upstream binary. t.Fatalf (never t.Skip) on any gap.

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
	bin := filepath.Join(t.TempDir(), "bedmerge")
	cmd := exec.Command("go", "build", "-o", bin,
		"github.com/yassineS/bio_ai_experiment/tools/bedmerge/cmd/bedmerge")
	cmd.Dir = root
	if out, buildErr := cmd.CombinedOutput(); buildErr != nil {
		t.Fatalf("build bedmerge: %v\n%s", buildErr, out)
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

// TestUpstreamCompat_Merge_MultiCharSingleDash drives both binaries with the
// upstream multi-char single-dash command lines and asserts our output equals
// the live upstream binary's.
func TestUpstreamCompat_Merge_MultiCharSingleDash(t *testing.T) {
	bt := upstreamBedtools(t)
	ours := buildOurs(t)
	dir := t.TempDir()
	plain := writeFile(t, dir, "m.bed", "chr1\t10\t20\nchr1\t15\t25\nchr1\t40\t50\n")
	withCol := writeFile(t, dir, "mc.bed", "chr1\t10\t20\t5\nchr1\t15\t25\t3\nchr1\t40\t50\t9\n")

	cases := []struct {
		name string
		in   string
		args []string
	}{
		{"default", plain, nil},
		{"distance", plain, []string{"-d", "5"}},
		{"col-sum", withCol, []string{"-c", "4", "-o", "sum"}},
		{"col-collapse", withCol, []string{"-c", "4", "-o", "collapse"}},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			want := runCapture(t, bt, append([]string{"merge", "-i", tc.in}, tc.args...)...)
			got := runCapture(t, ours, append([]string{"-i", tc.in}, tc.args...)...)
			if !bytes.Equal(got, want) {
				t.Fatalf("merge mismatch for %v\nupstream:\n%s\nours:\n%s", tc.args, want, got)
			}
		})
	}
}

// TestUpstreamCompat_Merge_StdinDash verifies `-i -` reads from stdin.
func TestUpstreamCompat_Merge_StdinDash(t *testing.T) {
	bt := upstreamBedtools(t)
	ours := buildOurs(t)
	input := "chr1\t10\t20\nchr1\t15\t25\n"
	runStdin := func(bin string, args ...string) []byte {
		cmd := exec.Command(bin, args...)
		cmd.Stdin = bytes.NewReader([]byte(input))
		var out, errb bytes.Buffer
		cmd.Stdout = &out
		cmd.Stderr = &errb
		if err := cmd.Run(); err != nil {
			t.Fatalf("run %s %v: %v\nstderr:\n%s", bin, args, err, errb.String())
		}
		return out.Bytes()
	}
	want := runStdin(bt, "merge", "-i", "-")
	got := runStdin(ours, "-i", "-")
	if !bytes.Equal(got, want) {
		t.Fatalf("stdin merge mismatch\nupstream:\n%s\nours:\n%s", want, got)
	}
}
