package main

import (
	"bytes"
	"image"
	_ "image/png"
	"os"
	"os/exec"
	"path/filepath"
	"sort"
	"strings"
	"testing"
)

// TestGraphPNGCLI exercises the `graph_png` subcommand end-to-end:
// it renders the vendored example .gd into the full PNG set (plus an
// HTML index) and verifies the filenames and that every PNG decodes
// to a valid, sufficiently-large image.
func TestGraphPNGCLI(t *testing.T) {
	bin := buildBinary(t)
	gd := filepath.Join("..", "..", "testdata", "parity", "graphdata_example1.gd")
	dir := t.TempDir()
	prefix := filepath.Join(dir, "rep")

	cmd := exec.Command(bin, "graph_png", "-i", gd, "-o", prefix, "--png_all", "--html_all")
	var stdout, stderr bytes.Buffer
	cmd.Stdout = &stdout
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("run: %v; stderr=%s", err, stderr.String())
	}

	// The stdout list names each written file.
	lines := strings.Fields(stdout.String())
	if len(lines) != 17 { // 16 PNG + 1 HTML
		t.Fatalf("expected 17 output paths, got %d: %v", len(lines), lines)
	}

	wantSuffixes := []string{
		"_cd.png", "_ce.png", "_df.png", "_dl.png", "_dm.png", "_gc.png",
		"_ld.png", "_ns.png", "_or.png", "_pm.png", "_pv.png", "_qd.png",
		"_qd2.png", "_qd3.png", "_td3.png", "_td5.png",
	}
	var gotSuffixes []string
	for _, p := range lines {
		base := filepath.Base(p)
		if strings.HasSuffix(base, ".html") {
			continue
		}
		gotSuffixes = append(gotSuffixes, strings.TrimPrefix(base, "rep"))
	}
	sort.Strings(gotSuffixes)
	if strings.Join(gotSuffixes, ",") != strings.Join(wantSuffixes, ",") {
		t.Fatalf("PNG suffixes = %v, want %v", gotSuffixes, wantSuffixes)
	}

	// Every PNG decodes and is reasonably sized.
	for _, p := range lines {
		if strings.HasSuffix(p, ".html") {
			b, err := os.ReadFile(p)
			if err != nil {
				t.Fatalf("read html: %v", err)
			}
			if !bytes.Contains(b, []byte("<img")) {
				t.Fatalf("html missing <img> tags")
			}
			continue
		}
		f, err := os.Open(p)
		if err != nil {
			t.Fatalf("open %s: %v", p, err)
		}
		cfg, format, err := image.DecodeConfig(f)
		f.Close()
		if err != nil {
			t.Fatalf("decode %s: %v", p, err)
		}
		if format != "png" {
			t.Fatalf("%s is %q, want png", p, format)
		}
		if cfg.Width < 50 || cfg.Height < 50 {
			t.Fatalf("%s too small: %dx%d", p, cfg.Width, cfg.Height)
		}
	}
}

// TestGraphPNGCLIDefaultPrefix verifies the derived-prefix path (no
// -o) strips the .gd extension and writes alongside the input.
func TestGraphPNGCLIDefaultPrefix(t *testing.T) {
	bin := buildBinary(t)
	src := filepath.Join("..", "..", "testdata", "parity", "graphdata_example1.gd")
	b, err := os.ReadFile(src)
	if err != nil {
		t.Fatalf("read fixture: %v", err)
	}
	dir := t.TempDir()
	gd := filepath.Join(dir, "sample.gd")
	if err := os.WriteFile(gd, b, 0o644); err != nil {
		t.Fatalf("write copy: %v", err)
	}

	cmd := exec.Command(bin, "graph_png", "-i", gd)
	var stderr bytes.Buffer
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("run: %v; stderr=%s", err, stderr.String())
	}
	// Default prefix is "<dir>/sample" (no .gd), so the length plot is
	// "<dir>/sample_ld.png".
	if _, err := os.Stat(filepath.Join(dir, "sample_ld.png")); err != nil {
		t.Fatalf("expected sample_ld.png: %v", err)
	}
}

// TestGraphPNGCLIMissingInput asserts the CLI errors without -i.
func TestGraphPNGCLIMissingInput(t *testing.T) {
	bin := buildBinary(t)
	cmd := exec.Command(bin, "graph_png")
	if err := cmd.Run(); err == nil {
		t.Fatal("expected non-zero exit without -i")
	}
}
