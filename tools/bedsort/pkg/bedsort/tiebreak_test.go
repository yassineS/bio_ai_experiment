package bedsort

// Tie-break parity tests.
//
// The decisive bug these guard against: for records with an equal
// (chrom, start) key, upstream bedtools sort preserves INPUT order (its
// per-chromosome sort compares chromStart alone), whereas an earlier port
// ordered ties by chromEnd ascending. These tests assert byte-for-byte
// equality against the live upstream bedtools binary for every sort mode, on
// an input deliberately containing many equal-(chrom, start) records that are
// NOT in end-sorted order (plus equal sizes and equal scores to exercise the
// size/score-mode tie-breaks).
//
// When the binary is absent the tests t.Fatalf (never t.Skip), per the parity
// policy. A separate binary-free TestUnitSortTieBreak* set pins the comparator
// behaviour so the core logic is covered even without upstream.

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

// upstreamBedtools resolves the live upstream bedtools binary, building it from
// the submodule if necessary. It t.Fatalf's (never skips) when unavailable.
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
		t.Fatalf("upstream bedtools unavailable: %v\n"+
			"run: git submodule update --init reference_code/bedtools && "+
			"(cd reference_code/bedtools && make -j\"$(nproc)\")", upstreamErr)
	}
	return upstreamPath
}

// upstreamSort runs `bedtools sort <flags> -i <inputFile>` and returns stdout.
func upstreamSort(t *testing.T, bin, inputFile string, flags ...string) []byte {
	t.Helper()
	args := append([]string{"sort"}, flags...)
	args = append(args, "-i", inputFile)
	var out, errBuf bytes.Buffer
	cmd := exec.Command(bin, args...)
	cmd.Stdout = &out
	cmd.Stderr = &errBuf
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream bedtools sort %v failed: %v\nstderr: %s", args, err, errBuf.String())
	}
	return out.Bytes()
}

// TestParity_Sort_TieBreak_AllModes asserts byte-for-byte parity with the live
// upstream binary across every sort mode on the equal-key stressor input.
func TestParity_Sort_TieBreak_AllModes(t *testing.T) {
	bin := upstreamBedtools(t)
	const input = "ties.bed"
	cases := []struct {
		name  string
		flag  string
		mode  SortMode
		faidx []string
	}{
		{"default", "", ModeChrom, nil},
		{"sizeA", "-sizeA", ModeSizeA, nil},
		{"sizeD", "-sizeD", ModeSizeD, nil},
		{"chrThenSizeA", "-chrThenSizeA", ModeChrThenSizeA, nil},
		{"chrThenSizeD", "-chrThenSizeD", ModeChrThenSizeD, nil},
		{"chrThenScoreA", "-chrThenScoreA", ModeChrThenScoreA, nil},
		{"chrThenScoreD", "-chrThenScoreD", ModeChrThenScoreD, nil},
	}
	for _, tc := range cases {
		tc := tc
		t.Run(tc.name, func(t *testing.T) {
			var flags []string
			if tc.flag != "" {
				flags = append(flags, tc.flag)
			}
			want := upstreamSort(t, bin, filepath.Join("..", "..", "testdata", "parity", input), flags...)
			got := runParity(t, input, Options{Mode: tc.mode, ChromOrder: tc.faidx})
			if !bytes.Equal(got, want) {
				t.Fatalf("mode %s mismatch.\nwant (upstream):\n%s\ngot:\n%s", tc.name, want, got)
			}
		})
	}
}

// TestParity_Sort_TieBreak_Faidx asserts the (chrom, start) input-order
// tie-break is also honoured when a faidx fixes the chromosome order.
func TestParity_Sort_TieBreak_Faidx(t *testing.T) {
	bin := upstreamBedtools(t)
	// names.txt orders chr9, chr3, chr7 — none of which are in ties.bed — so
	// craft a faidx covering the ties.bed chromosomes.
	dir := t.TempDir()
	faidx := filepath.Join(dir, "names.txt")
	if err := os.WriteFile(faidx, []byte("chr10\nchr2\nchr1\n"), 0o644); err != nil {
		t.Fatalf("write faidx: %v", err)
	}
	want := upstreamSort(t, bin, filepath.Join("..", "..", "testdata", "parity", "ties.bed"), "-faidx", faidx)
	order := []string{"chr10", "chr2", "chr1"}
	got := runParity(t, "ties.bed", Options{ChromOrder: order})
	if !bytes.Equal(got, want) {
		t.Fatalf("faidx tie-break mismatch.\nwant (upstream):\n%s\ngot:\n%s", want, got)
	}
}

// TestUnitSortTieBreakDefault pins the input-order tie-break on equal
// (chrom, start) keys without invoking the upstream binary.
func TestUnitSortTieBreakDefault(t *testing.T) {
	recs, _, err := readAll(bytes.NewReader([]byte(
		"chr1\t10\t100\ta\n"+
			"chr1\t10\t50\tb\n"+
			"chr1\t10\t75\tc\n"+
			"chr1\t5\t200\td\n")), false)
	if err != nil {
		t.Fatalf("readAll: %v", err)
	}
	Sort(recs, Options{})
	var got bytes.Buffer
	if err := Write(&got, recs); err != nil {
		t.Fatalf("Write: %v", err)
	}
	// start=5 first, then the start=10 group in INPUT order a,b,c (NOT end-asc
	// b,c,a).
	want := "chr1\t5\t200\td\n" +
		"chr1\t10\t100\ta\n" +
		"chr1\t10\t50\tb\n" +
		"chr1\t10\t75\tc\n"
	if got.String() != want {
		t.Fatalf("default tie-break.\nwant:\n%q\ngot:\n%q", want, got.String())
	}
}

// TestUnitSortTieBreakSizeAndScoreDesc pins that the size/score-descending
// modes keep input order on key ties (no chromEnd tie-break).
func TestUnitSortTieBreakSizeAndScoreDesc(t *testing.T) {
	// Three equal-(chrom,start) records, two with equal SIZE (50): inputs a,b
	// share size 50; their relative order must follow input (a before b).
	input := "chr1\t10\t60\ta\t20\t+\n" + // size 50, score 20
		"chr1\t10\t60\tb\t20\t+\n" + // size 50, score 20
		"chr1\t10\t90\tc\t99\t+\n" // size 80, score 99
	recs, _, err := readAll(bytes.NewReader([]byte(input)), false)
	if err != nil {
		t.Fatalf("readAll: %v", err)
	}
	Sort(recs, Options{Mode: ModeChrThenSizeD})
	var got bytes.Buffer
	if err := Write(&got, recs); err != nil {
		t.Fatalf("Write: %v", err)
	}
	want := "chr1\t10\t90\tc\t99\t+\n" + // largest first
		"chr1\t10\t60\ta\t20\t+\n" + // ties keep input order a,b
		"chr1\t10\t60\tb\t20\t+\n"
	if got.String() != want {
		t.Fatalf("chrThenSizeD tie-break.\nwant:\n%q\ngot:\n%q", want, got.String())
	}
}
