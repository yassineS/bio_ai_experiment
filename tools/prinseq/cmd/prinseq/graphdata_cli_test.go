package main

import (
	"bytes"
	"encoding/json"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

// buildBinary compiles the prinseq CLI into a temp file so the tests
// can exercise the actual command-line wiring. Returns the binary
// path; the temp dir is cleaned up via t.Cleanup.
func buildBinary(t *testing.T) string {
	t.Helper()
	dir := t.TempDir()
	bin := filepath.Join(dir, "prinseq")
	cmd := exec.Command("go", "build", "-o", bin, ".")
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("go build: %v", err)
	}
	return bin
}

// TestGraphDataCLI exercises the new subcommand end-to-end on the
// upstream example fixture.
func TestGraphDataCLI(t *testing.T) {
	bin := buildBinary(t)
	fastq := filepath.Join("..", "..", "testdata", "parity", "graphdata_example1.fastq")
	outPath := filepath.Join(t.TempDir(), "out.gd")

	cmd := exec.Command(bin, "graph_data",
		"--fastq", fastq,
		"--graph_data", outPath,
		"--no-header",
	)
	var stderr bytes.Buffer
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("run: %v; stderr=%s", err, stderr.String())
	}

	raw, err := os.ReadFile(outPath)
	if err != nil {
		t.Fatalf("read output: %v", err)
	}
	// Body parses as valid JSON.
	var v map[string]any
	if err := json.Unmarshal(raw, &v); err != nil {
		t.Fatalf("output is not valid JSON: %v; head=%s", err, string(raw[:200]))
	}
	// Required top-level keys.
	for _, k := range []string{"numseqs", "numbases", "maxlength", "binval", "counts", "stats", "quals", "qualsmean", "dinucodds"} {
		if _, ok := v[k]; !ok {
			t.Errorf("missing top-level key %q", k)
		}
	}
	if got, want := v["numseqs"].(float64), 12.0; got != want {
		t.Errorf("numseqs=%v, want %v", got, want)
	}
}

// TestGraphDataCLIDefaultPath verifies the upstream-default filename
// convention (<input>__.gd) when --graph_data is omitted.
func TestGraphDataCLIDefaultPath(t *testing.T) {
	bin := buildBinary(t)

	// Copy the fixture into a temp dir so the .gd lands somewhere we
	// can clean up.
	tmp := t.TempDir()
	src := filepath.Join("..", "..", "testdata", "parity", "graphdata_example1.fastq")
	dst := filepath.Join(tmp, "in.fastq")
	srcBytes, err := os.ReadFile(src)
	if err != nil {
		t.Fatalf("read fixture: %v", err)
	}
	if err := os.WriteFile(dst, srcBytes, 0o644); err != nil {
		t.Fatalf("write fixture: %v", err)
	}

	cmd := exec.Command(bin, "graph_data", "--fastq", dst, "--no-header")
	cmd.Dir = tmp
	var stderr bytes.Buffer
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("run: %v; stderr=%s", err, stderr.String())
	}
	want := dst + "__.gd"
	info, err := os.Stat(want)
	if err != nil {
		t.Fatalf("expected output %s: %v", want, err)
	}
	if info.Size() < 100 {
		t.Fatalf("output too small: %d bytes", info.Size())
	}
}

// TestGraphDataCLIHeader confirms that the default invocation emits
// the upstream-shaped two-line #-comment header.
func TestGraphDataCLIHeader(t *testing.T) {
	bin := buildBinary(t)
	fastq := filepath.Join("..", "..", "testdata", "parity", "graphdata_example1.fastq")
	outPath := filepath.Join(t.TempDir(), "out.gd")

	cmd := exec.Command(bin, "graph_data", "--fastq", fastq, "--graph_data", outPath)
	var stderr bytes.Buffer
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("run: %v; stderr=%s", err, stderr.String())
	}
	raw, err := os.ReadFile(outPath)
	if err != nil {
		t.Fatalf("read: %v", err)
	}
	first := strings.SplitN(string(raw), "\n", 3)
	if len(first) < 3 {
		t.Fatalf("expected at least 3 lines, got %d", len(first))
	}
	if first[0] != "#Graph data" {
		t.Errorf("line 1 = %q, want %q", first[0], "#Graph data")
	}
	if !strings.HasPrefix(first[1], "#[prinseq-lite-0.20.4]") {
		t.Errorf("line 2 prefix unexpected: %q", first[1])
	}
}
