package main

// CLI-level parity test for the default window. The bug fixed here: the -w
// default was 0, but upstream `bedtools window` uses a default window of
// 1000 bp. This drives the real run() entry point (where the default lives)
// and asserts byte-for-byte equality with the live upstream binary.

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
	upstreamOnce sync.Once
	upstreamPath string
	upstreamErr  error
)

func upstreamBedtools(t *testing.T) string {
	t.Helper()
	upstreamOnce.Do(func() {
		_, file, _, _ := runtime.Caller(0)
		dir := filepath.Dir(file)
		var root string
		for {
			if _, err := os.Stat(filepath.Join(dir, "go.mod")); err == nil {
				root = dir
				break
			}
			parent := filepath.Dir(dir)
			if parent == dir {
				upstreamErr = os.ErrNotExist
				return
			}
			dir = parent
		}
		btDir := filepath.Join(root, "reference_code", "bedtools")
		bin := filepath.Join(btDir, "bin", "bedtools")
		if _, statErr := os.Stat(bin); statErr != nil {
			cmd := exec.Command("make", "-j", "4")
			cmd.Dir = btDir
			if out, buildErr := cmd.CombinedOutput(); buildErr != nil {
				upstreamErr = &exec.ExitError{Stderr: out}
				return
			}
		}
		upstreamPath = bin
	})
	if upstreamErr != nil || upstreamPath == "" {
		t.Fatalf("upstream bedtools unavailable: %v", upstreamErr)
	}
	return upstreamPath
}

func writeTemp(t *testing.T, name, content string) string {
	t.Helper()
	p := filepath.Join(t.TempDir(), name)
	if err := os.WriteFile(p, []byte(content), 0o644); err != nil {
		t.Fatalf("write %s: %v", name, err)
	}
	return p
}

// TestCLI_DefaultWindowIs1000 builds our binary and compares its no-flag
// (default) output against upstream's, on a B record 800 bp upstream of A —
// in range only because the default window is 1000 bp.
func TestCLI_DefaultWindowIs1000(t *testing.T) {
	bin := upstreamBedtools(t)
	aContent := "chr1\t2000\t2100\ta1\t0\t+\n"
	bContent := "chr1\t1100\t1200\tb1\t0\t+\n"
	aFile := writeTemp(t, "a.bed", aContent)
	bFile := writeTemp(t, "b.bed", bContent)

	var up, upErr bytes.Buffer
	cmd := exec.Command(bin, "window", "-a", aFile, "-b", bFile)
	cmd.Stdout = &up
	cmd.Stderr = &upErr
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream window: %v\n%s", err, upErr.String())
	}
	if up.Len() == 0 {
		t.Fatal("expected upstream to report a hit at the default 1000bp window")
	}

	ours := buildOurs(t)
	var got, gotErr bytes.Buffer
	oc := exec.Command(ours, "-a", aFile, "-b", bFile)
	oc.Stdout = &got
	oc.Stderr = &gotErr
	if err := oc.Run(); err != nil {
		t.Fatalf("our bedwindow: %v\n%s", err, gotErr.String())
	}
	if !bytes.Equal(got.Bytes(), up.Bytes()) {
		t.Fatalf("default-window mismatch.\nupstream:\n%s\nours:\n%s", up.String(), got.String())
	}
}

func buildOurs(t *testing.T) string {
	t.Helper()
	_, file, _, _ := runtime.Caller(0)
	dir := filepath.Dir(file)
	var root string
	for {
		if _, err := os.Stat(filepath.Join(dir, "go.mod")); err == nil {
			root = dir
			break
		}
		dir = filepath.Dir(dir)
	}
	out := filepath.Join(t.TempDir(), "bedwindow")
	cmd := exec.Command("go", "build", "-o", out,
		"github.com/yassineS/bio_ai_experiment/tools/bedwindow/cmd/bedwindow")
	cmd.Dir = root
	if b, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("build bedwindow: %v\n%s", err, b)
	}
	return out
}
