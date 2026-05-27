package bedreldist

// Parity tests against the upstream `bedtools reldist` test suite under
// reference_code/bedtools/test/reldist/.
//
// All four upstream cases (t01..t04) are now exercised byte-for-byte; the
// three large-fixture cases (t01..t03) read vendored gzipped corpora
// (refseq.chr1.exons.bed.gz, aluY.chr1.bed.gz, gerp.chr1.bed.gz, total
// ~1.7MB compressed). t04 is the shipped issue_711 corner case.

import (
	"bytes"
	"compress/gzip"
	"io"
	"os"
	"path/filepath"
	"strings"
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

func openReldistParity(t *testing.T, name string) io.Reader {
	t.Helper()
	data := readReldistParity(t, name)
	if strings.HasSuffix(name, ".gz") {
		gr, err := gzip.NewReader(bytes.NewReader(data))
		if err != nil {
			t.Fatalf("gzip %s: %v", name, err)
		}
		return gr
	}
	return bytes.NewReader(data)
}

func runReldistParity(t *testing.T, aFile, bFile string, opts Options) []byte {
	t.Helper()
	a := openReldistParity(t, aFile)
	b := openReldistParity(t, bFile)
	var out bytes.Buffer
	if _, err := Run(a, b, &out, opts); err != nil {
		t.Fatalf("Run failed: %v", err)
	}
	return out.Bytes()
}

// reldist.t01 — self-intersect of refseq.chr1.exons.bed.gz (43,424 records).
// All relative distances are 0 by construction.
func TestParity_Reldist_T01_SelfIntersect_LargeFixture(t *testing.T) {
	got := runReldistParity(t, "refseq.chr1.exons.bed.gz", "refseq.chr1.exons.bed.gz", Options{})
	want := readReldistParity(t, "t01_self.expected.tsv")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// reldist.t02 — refseq exons vs randomly distributed aluY repeats on chr1.
// Distances are uniformly spread (all bins ~0.02).
func TestParity_Reldist_T02_RandomDistribution(t *testing.T) {
	got := runReldistParity(t, "refseq.chr1.exons.bed.gz", "aluY.chr1.bed.gz", Options{})
	want := readReldistParity(t, "t02_random.expected.tsv")
	if !bytes.Equal(got, want) {
		t.Fatalf("mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

// reldist.t03 — refseq exons vs gerp constrained elements on chr1.
// Distances are biased toward 0 (the 0.00 bin holds ~48% of the mass).
func TestParity_Reldist_T03_BiasedToZero(t *testing.T) {
	got := runReldistParity(t, "refseq.chr1.exons.bed.gz", "gerp.chr1.bed.gz", Options{})
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
