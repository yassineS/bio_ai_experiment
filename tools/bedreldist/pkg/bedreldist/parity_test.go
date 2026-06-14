package bedreldist

// Parity tests against the upstream `bedtools reldist` test suite under
// reference_code/bedtools/test/reldist/.
//
// All four upstream shell cases (t01..t04) are exercised here. t01..t03 use the
// upstream gzipped fixtures (refseq.chr1.exons.bed.gz, aluY.chr1.bed.gz,
// gerp.chr1.bed.gz), vendored under testdata/parity and decompressed at test
// time, and assert byte-for-byte against goldens generated from the upstream
// `bedtools reldist` binary (v2.31.1). t04 is the shipped issue_711 corner
// case, and a small hand-crafted self-intersect (small_self) is kept as a bonus.

import (
	"bytes"
	"compress/gzip"
	"os"
	"path/filepath"
	"testing"
)

// readReldistParityGz reads a gzipped fixture from testdata/parity and returns
// its decompressed bytes. The reldist fixtures are vendored compressed (the
// uncompressed BED would be several megabytes each).
func readReldistParityGz(t *testing.T, name string) []byte {
	t.Helper()
	f, err := os.Open(filepath.Join("..", "..", "testdata", "parity", name))
	if err != nil {
		t.Fatalf("open %s: %v", name, err)
	}
	defer f.Close()
	gr, err := gzip.NewReader(f)
	if err != nil {
		t.Fatalf("gzip %s: %v", name, err)
	}
	defer gr.Close()
	var buf bytes.Buffer
	if _, err := buf.ReadFrom(gr); err != nil {
		t.Fatalf("read %s: %v", name, err)
	}
	return buf.Bytes()
}

// runReldistParityGz runs reldist over two gzipped fixtures and returns the
// output bytes.
func runReldistParityGz(t *testing.T, aFile, bFile string, opts Options) []byte {
	t.Helper()
	a := readReldistParityGz(t, aFile)
	b := readReldistParityGz(t, bFile)
	var out bytes.Buffer
	if _, err := Run(bytes.NewReader(a), bytes.NewReader(b), &out, opts); err != nil {
		t.Fatalf("Run failed: %v", err)
	}
	return out.Bytes()
}

func readReldistParity(t *testing.T, name string) []byte {
	t.Helper()
	data, err := os.ReadFile(filepath.Join("..", "..", "testdata", "parity", name))
	if err != nil {
		t.Fatalf("read %s: %v", name, err)
	}
	return data
}

func runReldistParity(t *testing.T, aFile, bFile string, opts Options) []byte {
	t.Helper()
	a := readReldistParity(t, aFile)
	b := readReldistParity(t, bFile)
	var out bytes.Buffer
	if _, err := Run(bytes.NewReader(a), bytes.NewReader(b), &out, opts); err != nil {
		t.Fatalf("Run failed: %v", err)
	}
	return out.Bytes()
}

// reldist.t01 — self-intersect of refseq.chr1.exons (43k records); every
// relative distance is 0. Asserts byte-for-byte against the upstream golden.
func TestParity_Reldist_T01_SelfIntersect_LargeFixture(t *testing.T) {
	got := runReldistParityGz(t, "refseq.chr1.exons.bed.gz", "refseq.chr1.exons.bed.gz", Options{})
	want := readReldistParity(t, "t01_self.expected.tsv")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// reldist.t02 — refseq exons vs randomly distributed aluY elements; the
// relative distances are roughly uniform.
func TestParity_Reldist_T02_RandomDistribution(t *testing.T) {
	got := runReldistParityGz(t, "refseq.chr1.exons.bed.gz", "aluY.chr1.bed.gz", Options{})
	want := readReldistParity(t, "t02_random.expected.tsv")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// reldist.t03 — refseq exons vs gerp elements that cluster near the exons; the
// relative distances are biased towards 0.
func TestParity_Reldist_T03_BiasedToZero(t *testing.T) {
	got := runReldistParityGz(t, "refseq.chr1.exons.bed.gz", "gerp.chr1.bed.gz", Options{})
	want := readReldistParity(t, "t03_biased.expected.tsv")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// reldist.t04 — the shipped issue_711 corner case (a.bed has 2 records,
// b.bed has 3 partially-overlapping records). The expected output is
// byte-for-byte compared.
func TestParity_Reldist_T04_Issue711(t *testing.T) {
	got := runReldistParity(t, "issue_711.a.bed", "issue_711.b.bed", Options{})
	want := readReldistParity(t, "issue_711.expected.tsv")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// Bonus parity: self-intersect of a small fixture. Exercises the same
// algorithmic path as upstream's t01 (every query lands at distance 0) on a
// fixture we vendor in this repo.
func TestParity_Reldist_SmallSelf(t *testing.T) {
	got := runReldistParity(t, "small.bed", "small.bed", Options{})
	want := readReldistParity(t, "small_self.expected.tsv")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}
