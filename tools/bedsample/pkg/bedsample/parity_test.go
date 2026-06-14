package bedsample

// Parity tests against the upstream bedtools sample test suite.
//
// Cases mirror reference_code/bedtools/test/sample/test-sample.sh. The
// upstream suite generates its fixture at test time via `bedtools random`;
// we vendor an equivalent fixture (`mainFile.bed`, 1000 records with a
// `#header` line) under tools/bedsample/testdata/parity/.
//
// As of the std::mt19937_64 port (mt19937.go) this port also matches upstream
// BYTE-FOR-BYTE for a given seed; that is asserted directly against the live
// upstream binary in upstream_parity_test.go. The cases here cover the
// observable invariants without needing the binary:
//
//   - `-n N` yields exactly N records.
//   - Two runs with the same `-seed` are identical.
//   - `-header` forwards the header verbatim before the records.
//   - Requesting more records than the file has returns an error
//     (matching upstream's "Input file has fewer records ..." error).
//   - Each emitted record is a subset of the input (no fabrication).

import (
	"bytes"
	"errors"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func readParityFile(t *testing.T, name string) []byte {
	t.Helper()
	data, err := os.ReadFile(filepath.Join("..", "..", "testdata", "parity", name))
	if err != nil {
		t.Fatalf("read %s: %v", name, err)
	}
	return data
}

// sample.new.t07 — `-n 10` returns 10 records.
func TestParity_Sample_T07_TenRecords(t *testing.T) {
	in := readParityFile(t, "mainFile.bed")
	var buf bytes.Buffer
	n, err := Sample(bytes.NewReader(in), &buf, Options{N: 10, Seed: 1})
	if err != nil {
		t.Fatalf("Sample: %v", err)
	}
	if n != 10 {
		t.Errorf("count = %d, want 10", n)
	}
	lines := strings.Count(strings.TrimRight(buf.String(), "\n"), "\n") + 1
	if lines != 10 {
		t.Errorf("output line count = %d, want 10", lines)
	}
}

// sample.new.t08 — `-seed` is deterministic.
func TestParity_Sample_T08_DeterministicSeed(t *testing.T) {
	in := readParityFile(t, "mainFile.bed")
	var a, b bytes.Buffer
	if _, err := Sample(bytes.NewReader(in), &a, Options{N: 50, Seed: 4}); err != nil {
		t.Fatalf("sample a: %v", err)
	}
	if _, err := Sample(bytes.NewReader(in), &b, Options{N: 50, Seed: 4}); err != nil {
		t.Fatalf("sample b: %v", err)
	}
	if !bytes.Equal(a.Bytes(), b.Bytes()) {
		t.Fatalf("seeded runs disagree.\nA:\n%s\nB:\n%s", a.Bytes(), b.Bytes())
	}
}

// sample.new.t09 — `-header` keeps the header line first.
func TestParity_Sample_T09_HeaderForwarded(t *testing.T) {
	in := readParityFile(t, "mainFile.bed")
	var buf bytes.Buffer
	if _, err := Sample(bytes.NewReader(in), &buf, Options{N: 10, Seed: 1, Header: true}); err != nil {
		t.Fatalf("Sample: %v", err)
	}
	first := strings.SplitN(buf.String(), "\n", 2)[0]
	want := "#header line for the parity fixture"
	if first != want {
		t.Errorf("first line = %q, want %q", first, want)
	}
}

// sample.new.t06 — requesting N greater than total raises the upstream
// "Input file has fewer records ..." error.
func TestParity_Sample_T06_TooFewRecords(t *testing.T) {
	in := readParityFile(t, "mainFile.bed")
	var buf bytes.Buffer
	_, err := Sample(bytes.NewReader(in), &buf, Options{N: 2000, Seed: 1})
	if err == nil {
		t.Fatalf("expected error for N > total")
	}
	var tfr *ErrTooFewRecords
	if !errors.As(err, &tfr) {
		t.Errorf("err type = %T, want *ErrTooFewRecords", err)
	}
	if !strings.Contains(err.Error(), "fewer records") {
		t.Errorf("err text %q lacks 'fewer records'", err)
	}
}

// Additional invariant: sampled records are a SUBSET of input records.
func TestParity_Sample_SubsetOfInput(t *testing.T) {
	in := readParityFile(t, "mainFile.bed")
	var buf bytes.Buffer
	if _, err := Sample(bytes.NewReader(in), &buf, Options{N: 50, Seed: 4}); err != nil {
		t.Fatalf("Sample: %v", err)
	}
	have := make(map[string]bool, 1000)
	for _, line := range strings.Split(string(in), "\n") {
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		have[line] = true
	}
	for _, line := range strings.Split(strings.TrimRight(buf.String(), "\n"), "\n") {
		if !have[line] {
			t.Errorf("sampled line not present in input: %q", line)
		}
	}
}

// sample.t01 — "No input file given". This is intentionally retained as a
// CLI-only documented case, NOT a library gap: upstream emits the message from
// its argument layer, whereas the Sample library always takes an io.Reader and
// the bedsample CLI (cmd/bedsample/main.go) defaults a missing -i to stdin (a
// deliberate POSIX drop-in improvement). The argument-validation/exit-code
// behaviour therefore lives in main.go, not here.
//
// covered by: cmd/bedsample/main.go (-n validation + iohelper.OpenReader,
// which surfaces open errors and exits non-zero).
func TestParity_Sample_T01_NoArgs(t *testing.T) {
	t.Skip("CLI-only (intentional): input defaulting/validation lives in cmd/bedsample/main.go, not the Sample library")
}

// sample.new.t02 — "Unrecognized parameter". Also a CLI-only documented case:
// an unknown flag is rejected by the standard flag.FlagSet in
// cmd/bedsample/main.go (fs.Parse returns an error and main exits 2), so there
// is no library-level behaviour to assert here.
//
// covered by: cmd/bedsample/main.go (fs.Parse error path -> os.Exit(2)).
func TestParity_Sample_T02_UnrecognizedFlag(t *testing.T) {
	t.Skip("CLI-only (intentional): unknown-flag rejection lives in cmd/bedsample/main.go's flag.FlagSet, not the Sample library")
}
