package bedreldist

// Parity tests against the upstream `bedtools reldist` test suite under
// reference_code/bedtools/test/reldist/.
//
// Three of upstream's four shell cases (t01..t03) rely on multi-megabyte
// gzipped fixtures (refseq.chr1.exons.bed.gz, aluY.chr1.bed.gz,
// gerp.chr1.bed.gz) that we do not vendor — we exercise the same code path
// with a small hand-crafted self-intersect (small_self) instead, plus the
// shipped issue_711 corner case.

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"
)

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

// reldist.t01 — self-intersect. The full upstream case uses the
// refseq.chr1.exons.bed.gz fixture (~1MB gzipped, 43k records). We do not
// vendor that; the same code path is exercised by small_self below.
func TestParity_Reldist_T01_SelfIntersect_LargeFixture(t *testing.T) {
	t.Skip("not vendored: refseq.chr1.exons.bed.gz; covered by TestParity_Reldist_SmallSelf instead")
}

// reldist.t02 / t03 — same situation: external large fixtures.
func TestParity_Reldist_T02_RandomDistribution(t *testing.T) {
	t.Skip("not vendored: aluY.chr1.bed.gz / refseq.chr1.exons.bed.gz")
}
func TestParity_Reldist_T03_BiasedToZero(t *testing.T) {
	t.Skip("not vendored: gerp.chr1.bed.gz / refseq.chr1.exons.bed.gz")
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
