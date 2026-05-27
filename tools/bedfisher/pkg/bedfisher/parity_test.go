package bedfisher

// Parity tests against the upstream `bedtools fisher` test suite under
// reference_code/bedtools/test/fisher/. All five small cases (t1..t4, t6)
// match byte-for-byte; t5 depends on a long file path inside $TMPDIR which
// is a CLI/filesystem concern unrelated to the algorithm.

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"
)

func readFisherParity(t *testing.T, name string) []byte {
	t.Helper()
	data, err := os.ReadFile(filepath.Join("..", "..", "testdata", "parity", name))
	if err != nil {
		t.Fatalf("read %s: %v", name, err)
	}
	return data
}

func runFisherParity(t *testing.T, aFile, bFile, gFile string, opts Options) []byte {
	t.Helper()
	a := readFisherParity(t, aFile)
	b := readFisherParity(t, bFile)
	g := readFisherParity(t, gFile)
	var out bytes.Buffer
	if _, err := Run(bytes.NewReader(a), bytes.NewReader(b), bytes.NewReader(g), &out, opts); err != nil {
		t.Fatalf("Run failed: %v", err)
	}
	return out.Bytes()
}

// fisher.t1 — a.bed vs b.bed on a 500-bp genome.
func TestParity_Fisher_T1(t *testing.T) {
	got := runFisherParity(t, "a.bed", "b.bed", "t.500.genome", Options{})
	want := readFisherParity(t, "t1.expected.txt")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// fisher.t2 — same inputs as t1 but with a 60-bp genome.
func TestParity_Fisher_T2(t *testing.T) {
	got := runFisherParity(t, "a.bed", "b.bed", "t.60.genome", Options{})
	want := readFisherParity(t, "t2.expected.txt")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// fisher.t3 — a_merge.bed (4 raw records, 2 of which overlap) vs b.bed
// without -m. The 4 raw query intervals are preserved.
func TestParity_Fisher_T3(t *testing.T) {
	got := runFisherParity(t, "a_merge.bed", "b.bed", "t.60.genome", Options{})
	want := readFisherParity(t, "t3.expected.txt")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// fisher.t4 — same as t3 but with -m: the overlapping records in
// a_merge.bed are pre-merged, dropping the query count from 4 to 3.
func TestParity_Fisher_T4(t *testing.T) {
	got := runFisherParity(t, "a_merge.bed", "b.bed", "t.60.genome", Options{MergeA: true})
	want := readFisherParity(t, "t4.expected.txt")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// fisher.t5 — upstream verifies that `bedtools fisher` happily opens a
// deeply nested $TMPDIR path. Our port operates on io.Reader so the
// only filesystem concern is callers passing OS-opened files. This
// test exercises that path with a synthetic deep tempdir tree and
// asserts the same Fisher exact result as t1 (same fixtures, just
// addressed via a long absolute path). Intentional non-divergence
// pinned as a real assertion rather than a skip.
func TestParity_Fisher_T5_LongPath(t *testing.T) {
	dir := t.TempDir()
	// Build a deep nested path under t.TempDir. 32 segments of 8 bytes
	// each plus separators is comfortably under PATH_MAX on Linux/macOS
	// (4096 / 1024) but exercises any per-segment buffer caps.
	deep := dir
	for i := 0; i < 32; i++ {
		deep = filepath.Join(deep, "depth_seg")
	}
	if err := os.MkdirAll(deep, 0o755); err != nil {
		t.Fatalf("MkdirAll(%s): %v", deep, err)
	}
	for _, name := range []string{"a.bed", "b.bed", "t.500.genome", "t1.expected.txt"} {
		data := readFisherParity(t, name)
		if err := os.WriteFile(filepath.Join(deep, name), data, 0o644); err != nil {
			t.Fatalf("write %s: %v", name, err)
		}
	}
	openFromDeep := func(name string) *os.File {
		f, err := os.Open(filepath.Join(deep, name))
		if err != nil {
			t.Fatalf("open %s: %v", name, err)
		}
		return f
	}
	a := openFromDeep("a.bed")
	defer a.Close()
	b := openFromDeep("b.bed")
	defer b.Close()
	g := openFromDeep("t.500.genome")
	defer g.Close()
	var out bytes.Buffer
	if _, err := Run(a, b, g, &out, Options{}); err != nil {
		t.Fatalf("Run from deep path: %v", err)
	}
	want, err := os.ReadFile(filepath.Join(deep, "t1.expected.txt"))
	if err != nil {
		t.Fatalf("read expected: %v", err)
	}
	if !bytes.Equal(out.Bytes(), want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, out.Bytes())
	}
}

// fisher.t6 — issue 954 regression: 5 query intervals across 5 chromosomes,
// 9 db intervals, 2 overlaps; exercises a non-degenerate ratio (5.905).
func TestParity_Fisher_T6_Issue954(t *testing.T) {
	got := runFisherParity(t, "issue954_a.bed", "issue954_b.bed", "issue954.genome", Options{})
	want := readFisherParity(t, "t6.expected.txt")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}
