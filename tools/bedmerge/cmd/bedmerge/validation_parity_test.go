package main

// Differential-validation parity tests for bedmerge against the live upstream
// bedtools binary (gap A18). Upstream `bedtools merge` REJECTS malformed input
// that an earlier bedmerge accepted: unsorted (not coordinate-sorted) input and
// records with inconsistent field counts. These tests assert byte-for-byte
// parity of the exit code AND stderr for those rejections, plus that valid
// inputs are still accepted unchanged.
//
// Unlike upstream_compat_test.go (which Fatalf's when upstream is missing),
// these skip gracefully when the upstream binary cannot be located/built, so
// the suite stays green in environments without the bedtools submodule.

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"sync"
	"testing"
)

var (
	upstreamOptOnce sync.Once
	upstreamOptPath string
)

// upstreamBedtoolsOptional returns the path to a usable upstream bedtools
// binary, or "" (with t.Skip) when it cannot be located or built. It never
// fails the test, so CI without the submodule still passes.
func upstreamBedtoolsOptional(t *testing.T) string {
	t.Helper()
	if testing.Short() {
		t.Skip("skipping upstream-binary parity test in -short mode")
	}
	upstreamOptOnce.Do(func() {
		root, err := repoRoot()
		if err != nil {
			return
		}
		dir := filepath.Join(root, "reference_code", "bedtools")
		bin := filepath.Join(dir, "bin", "bedtools")
		if _, statErr := os.Stat(bin); statErr == nil {
			upstreamOptPath = bin
			return
		}
		// Try a build only if the source Makefile is present.
		if _, statErr := os.Stat(filepath.Join(dir, "Makefile")); statErr != nil {
			return
		}
		cmd := exec.Command("make", "-j", "4")
		cmd.Dir = dir
		if _, buildErr := cmd.CombinedOutput(); buildErr != nil {
			return
		}
		if _, statErr := os.Stat(bin); statErr == nil {
			upstreamOptPath = bin
		}
	})
	if upstreamOptPath == "" {
		t.Skip("upstream bedtools binary unavailable; skipping validation-parity test")
	}
	return upstreamOptPath
}

// runExit runs a command (with no stdin) and returns its exit code, stdout, and
// stderr. A non-zero exit is not a test failure here — the exit code is the
// value under test.
func runExit(t *testing.T, name string, args ...string) (code int, stdout, stderr []byte) {
	t.Helper()
	cmd := exec.Command(name, args...)
	var out, errb bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = &errb
	err := cmd.Run()
	code = 0
	if err != nil {
		if ee, ok := err.(*exec.ExitError); ok {
			code = ee.ExitCode()
		} else {
			t.Fatalf("run %s %v: %v", name, args, err)
		}
	}
	return code, out.Bytes(), errb.Bytes()
}

// TestValidationParity_RejectsMalformed asserts that bedmerge rejects unsorted
// and ragged-field input exactly as upstream bedtools merge does — identical
// exit code and identical stderr — and accepts valid input identically.
func TestValidationParity_RejectsMalformed(t *testing.T) {
	bt := upstreamBedtoolsOptional(t)
	ours := buildOurs(t)
	dir := t.TempDir()

	cases := []struct {
		name    string
		content string
	}{
		{
			// Out-of-order start within a chromosome: upstream errors with
			// "Sorted input specified, but the file ... out of order record".
			name:    "unsorted-start",
			content: "chr1\t100\t200\nchr1\t50\t150\n",
		},
		{
			// Same, with extra columns echoed in the offending-record print.
			name:    "unsorted-cols",
			content: "chr1\t100\t200\tA\t5\t+\nchr1\t50\t150\tB\t6\t-\n",
		},
		{
			// A chromosome reappears after a different one: out of order.
			name:    "chrom-revisit",
			content: "chr1\t100\t200\nchr2\t50\t150\nchr1\t10\t20\n",
		},
		{
			// Inconsistent field counts within the first four valid lines:
			// upstream's type checker rejects with the "wrong number of fields"
			// message.
			name:    "ragged-early",
			content: "chr1\t100\t200\tfoo\nchr1\t250\t350\n",
		},
		{
			// Inconsistent field count on a later line (5th valid): upstream's
			// per-line reader reports "line number N ... has X fields ...".
			name:    "ragged-late",
			content: "chr1\t1\t2\nchr1\t3\t4\nchr1\t5\t6\nchr1\t7\t8\nchr1\t9\t10\tx\n",
		},
		{
			// A zero-length out-of-order record: the printed line uses upstream's
			// adjusted (start-1, end+1) coordinates.
			name:    "zerolen-out-of-order",
			content: "chr1\t100\t200\nchr1\t99\t99\n",
		},
		{
			// Valid, sorted input: both accept it with identical output.
			name:    "valid-sorted",
			content: "chr1\t100\t200\nchr1\t250\t350\n",
		},
	}

	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			in := writeFile(t, dir, tc.name+".bed", tc.content)
			wantCode, wantOut, wantErr := runExit(t, bt, "merge", "-i", in)
			gotCode, gotOut, gotErr := runExit(t, ours, "-i", in)

			if gotCode != wantCode {
				t.Errorf("exit code: got %d, want %d", gotCode, wantCode)
			}
			if !bytes.Equal(gotErr, wantErr) {
				t.Errorf("stderr mismatch:\nupstream: %q\nours:     %q", wantErr, gotErr)
			}
			if !bytes.Equal(gotOut, wantOut) {
				t.Errorf("stdout mismatch:\nupstream: %q\nours:     %q", wantOut, gotOut)
			}
		})
	}
}

// TestValidationParity_Stdin asserts unsorted input on stdin (`-i -`) is
// rejected with the same exit code and stderr (filename reported as "-").
func TestValidationParity_Stdin(t *testing.T) {
	bt := upstreamBedtoolsOptional(t)
	ours := buildOurs(t)
	input := []byte("chr1\t100\t200\nchr1\t50\t150\n")

	runStdin := func(bin string, args ...string) (int, []byte, []byte) {
		cmd := exec.Command(bin, args...)
		cmd.Stdin = bytes.NewReader(input)
		var out, errb bytes.Buffer
		cmd.Stdout = &out
		cmd.Stderr = &errb
		code := 0
		if err := cmd.Run(); err != nil {
			if ee, ok := err.(*exec.ExitError); ok {
				code = ee.ExitCode()
			} else {
				t.Fatalf("run %s: %v", bin, err)
			}
		}
		return code, out.Bytes(), errb.Bytes()
	}

	wantCode, _, wantErr := runStdin(bt, "merge", "-i", "-")
	gotCode, _, gotErr := runStdin(ours, "-i", "-")
	if gotCode != wantCode {
		t.Errorf("stdin exit code: got %d, want %d", gotCode, wantCode)
	}
	if !bytes.Equal(gotErr, wantErr) {
		t.Errorf("stdin stderr mismatch:\nupstream: %q\nours:     %q", wantErr, gotErr)
	}
}
