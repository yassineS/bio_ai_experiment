package bedsplit

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func TestParseAlgorithm(t *testing.T) {
	cases := []struct {
		in   string
		want Algorithm
		err  bool
	}{
		{"simple", AlgSimple, false},
		{"size", AlgSize, false},
		{"SIMPLE", AlgSimple, false},
		{" Size ", AlgSize, false},
		{"foo", 0, true},
	}
	for _, c := range cases {
		got, err := ParseAlgorithm(c.in)
		if c.err {
			if err == nil {
				t.Errorf("ParseAlgorithm(%q) want err, got %v", c.in, got)
			}
			continue
		}
		if err != nil {
			t.Errorf("ParseAlgorithm(%q) err: %v", c.in, err)
		}
		if got != c.want {
			t.Errorf("ParseAlgorithm(%q) = %v, want %v", c.in, got, c.want)
		}
	}
}

func TestSimpleBins_RoundRobin(t *testing.T) {
	// 10 records, 5 bins, round-robin: bin0=[0,5], bin1=[1,6], ...
	bins := simpleBins(10, 5)
	if len(bins) != 5 {
		t.Fatalf("len = %d, want 5", len(bins))
	}
	for i, b := range bins {
		if len(b) != 2 {
			t.Errorf("bin %d size = %d, want 2", i, len(b))
		}
	}
	if bins[0][0] != 0 || bins[0][1] != 5 {
		t.Errorf("bin 0 = %v, want [0,5]", bins[0])
	}
	if bins[4][0] != 4 || bins[4][1] != 9 {
		t.Errorf("bin 4 = %v, want [4,9]", bins[4])
	}
}

func TestSimpleBins_UnevenSplit(t *testing.T) {
	// 7 records, 3 bins, round-robin -> bin0=[0,3,6], bin1=[1,4], bin2=[2,5]
	bins := simpleBins(7, 3)
	if len(bins[0]) != 3 || len(bins[1]) != 2 || len(bins[2]) != 2 {
		t.Errorf("uneven bin sizes = %d,%d,%d; want 3,2,2",
			len(bins[0]), len(bins[1]), len(bins[2]))
	}
	if bins[0][0] != 0 || bins[0][1] != 3 || bins[0][2] != 6 {
		t.Errorf("bin 0 = %v, want [0,3,6]", bins[0])
	}
}

func TestSizeBins_BalancesByLength(t *testing.T) {
	// Lengths: 100, 50, 50, 50 -> 2 bins. Whatever the algorithm, the totals
	// must sum to 250 and each bin must hold at least one record.
	records := []record{
		{length: 100, line: "a"},
		{length: 50, line: "b"},
		{length: 50, line: "c"},
		{length: 50, line: "d"},
	}
	bins := sizeBins(records, 2)
	totals := make([]int, len(bins))
	count := 0
	for i, b := range bins {
		for _, ix := range b {
			totals[i] += records[ix].length
		}
		count += len(b)
	}
	if totals[0]+totals[1] != 250 {
		t.Errorf("totals sum = %d, want 250", totals[0]+totals[1])
	}
	if count != 4 {
		t.Errorf("count = %d, want 4", count)
	}
	if len(bins[0]) == 0 || len(bins[1]) == 0 {
		t.Errorf("expected non-empty bins, got %v", bins)
	}
	// Diff between bins should be <= max-record-length (100); LPT/heuristic
	// will usually do better.
	diff := totals[0] - totals[1]
	if diff < 0 {
		diff = -diff
	}
	if diff > 100 {
		t.Errorf("bins very unbalanced: %v", totals)
	}
}

func TestSplit_EndToEnd_Simple(t *testing.T) {
	dir := t.TempDir()
	prefix := filepath.Join(dir, "out")
	in := strings.NewReader(
		"chr1\t0\t10\nchr1\t20\t30\nchr1\t40\t50\nchr1\t60\t70\n",
	)
	var manifest bytes.Buffer
	rows, err := Split(in, &manifest, Options{N: 2, Prefix: prefix, Algorithm: AlgSimple})
	if err != nil {
		t.Fatalf("Split: %v", err)
	}
	if len(rows) != 2 {
		t.Fatalf("rows = %d, want 2", len(rows))
	}
	for i, r := range rows {
		if r.NumRecords != 2 || r.TotalBP != 20 {
			t.Errorf("row[%d] = %+v, want NumRecords=2 TotalBP=20", i, r)
		}
		data, err := os.ReadFile(r.Filename)
		if err != nil {
			t.Errorf("reading %s: %v", r.Filename, err)
			continue
		}
		if len(data) == 0 {
			t.Errorf("file %s is empty", r.Filename)
		}
	}
	// First file must contain the first two records, in order.
	data, _ := os.ReadFile(rows[0].Filename)
	if !strings.HasPrefix(string(data), "chr1\t0\t10\n") {
		t.Errorf("simple-mode file 1 should start with first input record, got: %q", string(data))
	}
}

func TestSplit_EndToEnd_Size(t *testing.T) {
	dir := t.TempDir()
	prefix := filepath.Join(dir, "out")
	// Lengths: 100, 1, 1, 1 — size mode should put 100 alone in one shard.
	in := strings.NewReader(
		"chr1\t0\t100\nchr1\t200\t201\nchr1\t300\t301\nchr1\t400\t401\n",
	)
	var manifest bytes.Buffer
	rows, err := Split(in, &manifest, Options{N: 2, Prefix: prefix, Algorithm: AlgSize})
	if err != nil {
		t.Fatalf("Split: %v", err)
	}
	if len(rows) != 2 {
		t.Fatalf("rows = %d, want 2", len(rows))
	}
	// One shard has 1 record (length 100), the other has 3 records (length 3).
	totals := make([]int, 2)
	counts := make([]int, 2)
	for i, r := range rows {
		totals[i] = r.TotalBP
		counts[i] = r.NumRecords
	}
	if !((totals[0] == 100 && counts[0] == 1 && totals[1] == 3 && counts[1] == 3) ||
		(totals[0] == 3 && counts[0] == 3 && totals[1] == 100 && counts[1] == 1)) {
		t.Errorf("unexpected balance: totals=%v counts=%v", totals, counts)
	}
}

func TestSplit_FilenameFormat(t *testing.T) {
	dir := t.TempDir()
	prefix := filepath.Join(dir, "tmp")
	in := strings.NewReader("chr1\t0\t10\nchr1\t10\t20\n")
	var manifest bytes.Buffer
	rows, err := Split(in, &manifest, Options{N: 2, Prefix: prefix, Algorithm: AlgSimple})
	if err != nil {
		t.Fatalf("Split: %v", err)
	}
	if !strings.HasSuffix(rows[0].Filename, ".00001.bed") {
		t.Errorf("rows[0].Filename = %q, want suffix .00001.bed", rows[0].Filename)
	}
	if !strings.HasSuffix(rows[1].Filename, ".00002.bed") {
		t.Errorf("rows[1].Filename = %q, want suffix .00002.bed", rows[1].Filename)
	}
}

func TestSplit_NMoreThanRecords(t *testing.T) {
	dir := t.TempDir()
	prefix := filepath.Join(dir, "x")
	in := strings.NewReader("chr1\t0\t10\nchr1\t10\t20\n")
	var manifest bytes.Buffer
	rows, err := Split(in, &manifest, Options{N: 10, Prefix: prefix, Algorithm: AlgSimple})
	if err != nil {
		t.Fatalf("Split: %v", err)
	}
	// Should cap at len(records) = 2.
	if len(rows) != 2 {
		t.Errorf("rows = %d, want 2", len(rows))
	}
}

func TestSplit_EmptyInput(t *testing.T) {
	dir := t.TempDir()
	prefix := filepath.Join(dir, "x")
	rows, err := Split(strings.NewReader(""), nil, Options{N: 5, Prefix: prefix, Algorithm: AlgSimple})
	if err != nil {
		t.Fatalf("Split: %v", err)
	}
	if rows != nil {
		t.Errorf("rows = %v, want nil", rows)
	}
}

func TestSplit_BadOptions(t *testing.T) {
	if _, err := Split(strings.NewReader(""), nil, Options{N: 0, Prefix: "x"}); err == nil {
		t.Error("expected error for N=0")
	}
	if _, err := Split(strings.NewReader(""), nil, Options{N: 2, Prefix: ""}); err == nil {
		t.Error("expected error for empty prefix")
	}
}

func TestSplit_BadInput(t *testing.T) {
	dir := t.TempDir()
	prefix := filepath.Join(dir, "x")
	in := strings.NewReader("chr1\tnot-a-number\t100\n")
	if _, err := Split(in, nil, Options{N: 1, Prefix: prefix, Algorithm: AlgSimple}); err == nil {
		t.Error("expected error on bad input")
	}
}

func TestSplit_ManifestFormat(t *testing.T) {
	dir := t.TempDir()
	prefix := filepath.Join(dir, "x")
	in := strings.NewReader("chr1\t0\t10\nchr1\t20\t30\n")
	var manifest bytes.Buffer
	if _, err := Split(in, &manifest, Options{N: 2, Prefix: prefix, Algorithm: AlgSimple}); err != nil {
		t.Fatalf("Split: %v", err)
	}
	lines := strings.Split(strings.TrimRight(manifest.String(), "\n"), "\n")
	if len(lines) != 2 {
		t.Fatalf("manifest lines = %d, want 2", len(lines))
	}
	for _, ln := range lines {
		fields := strings.Split(ln, "\t")
		if len(fields) != 3 {
			t.Errorf("manifest line %q has %d fields, want 3", ln, len(fields))
		}
	}
}
