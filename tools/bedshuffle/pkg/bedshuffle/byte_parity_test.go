package bedshuffle

// Byte-for-byte parity tests against the upstream `bedtools shuffle` binary.
//
// Each expected.*.bed golden under testdata/parity/ was produced by running
// reference_code/bedtools/bin/bedtools (v2.31.1, the default non-USE_RAND
// build that uses std::mt19937_64) with `shuffle -seed N` on the vendored
// fixtures. This port ports that exact MT19937-64 engine (mt19937.go), the
// rand_range rejection bound, the genome-file-order projection, and the
// per-mode draw/retry order from shuffleBed.cpp, so its output must match the
// golden byte-for-byte. Any divergence is a t.Fatalf — there is no skip path.

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"
)

// loadGenomeOrdered reads a genome fixture preserving chromosome file order
// (required for byte-exact parity in the default genome-projection mode).
func loadGenomeOrdered(t *testing.T, name string) (map[string]int, []string) {
	t.Helper()
	g, order, err := ParseGenomeOrdered(bytes.NewReader(readParityFixture(t, name)))
	if err != nil {
		t.Fatalf("parse genome %s: %v", name, err)
	}
	return g, order
}

// runShuffleBytes runs Shuffle on the named input fixture and returns the raw
// output bytes, failing on any error (over-large/unplaceable cases are not
// used by these byte tests).
func runShuffleBytes(t *testing.T, inputFixture string, opts Options) []byte {
	t.Helper()
	in := readParityFixture(t, inputFixture)
	var buf bytes.Buffer
	if _, err := Shuffle(bytes.NewReader(in), &buf, opts); err != nil {
		t.Fatalf("Shuffle(%s): %v", inputFixture, err)
	}
	return buf.Bytes()
}

// assertGolden compares got to the named golden file byte-for-byte, reporting
// the first differing line.
func assertGolden(t *testing.T, got []byte, goldenName string) {
	t.Helper()
	want, err := os.ReadFile(filepath.Join("..", "..", "testdata", "parity", goldenName))
	if err != nil {
		t.Fatalf("read golden %s: %v", goldenName, err)
	}
	if bytes.Equal(got, want) {
		return
	}
	gotLines := bytes.Split(got, []byte("\n"))
	wantLines := bytes.Split(want, []byte("\n"))
	n := len(gotLines)
	if len(wantLines) < n {
		n = len(wantLines)
	}
	for i := 0; i < n; i++ {
		if !bytes.Equal(gotLines[i], wantLines[i]) {
			t.Fatalf("byte mismatch vs %s at line %d:\n  got:  %q\n  want: %q",
				goldenName, i+1, gotLines[i], wantLines[i])
		}
	}
	t.Fatalf("byte mismatch vs %s: line counts differ got=%d want=%d",
		goldenName, len(gotLines), len(wantLines))
}

// TestByteParity_Default covers the default genome-projection mode on the real
// hg19 fixture across several seeds.
func TestByteParity_Default(t *testing.T) {
	g, order := loadGenomeOrdered(t, "human.hg19.genome")
	for _, seed := range []int64{1, 7, 42, 100} {
		seed := seed
		t.Run("seed_"+itoa(seed), func(t *testing.T) {
			got := runShuffleBytes(t, "simrep.bed", Options{Genome: g, GenomeOrder: order, Seed: seed})
			assertGolden(t, got, "expected.default.seed"+itoa(seed)+".bed")
		})
	}
}

// TestByteParity_Include covers -incl (size-weighted include selection).
func TestByteParity_Include(t *testing.T) {
	g, order := loadGenomeOrdered(t, "human.hg19.genome")
	incl := readParityBED(t, "incl.bed")
	got := runShuffleBytes(t, "simrep.bed", Options{Genome: g, GenomeOrder: order, Include: incl, Seed: 42})
	assertGolden(t, got, "expected.incl.seed42.bed")
}

// TestByteParity_Exclude covers -excl (redraw while overlapping an exclude
// region).
func TestByteParity_Exclude(t *testing.T) {
	g, order := loadGenomeOrdered(t, "human.hg19.genome")
	excl := readParityBED(t, "excl.bed")
	got := runShuffleBytes(t, "simrep.bed", Options{Genome: g, GenomeOrder: order, Exclude: excl, Seed: 42})
	assertGolden(t, got, "expected.excl.seed42.bed")
}

// TestByteParity_Chrom covers -chrom (keep each interval on its original
// chromosome).
func TestByteParity_Chrom(t *testing.T) {
	g, order := loadGenomeOrdered(t, "human.hg19.genome")
	got := runShuffleBytes(t, "simrep.bed", Options{Genome: g, GenomeOrder: order, Chrom: true, Seed: 42})
	assertGolden(t, got, "expected.chrom.seed42.bed")
}

// TestByteParity_ChromFirst covers -chromFirst (pick a chromosome uniformly,
// then a position) on the real fixture.
func TestByteParity_ChromFirst(t *testing.T) {
	g, order := loadGenomeOrdered(t, "human.hg19.genome")
	got := runShuffleBytes(t, "simrep.bed", Options{Genome: g, GenomeOrder: order, ChromFirst: true, Seed: 42})
	assertGolden(t, got, "expected.chromFirst.seed42.bed")
}

// TestByteParity_ChromFirstSmall covers -chromFirst on the small multi-chrom
// genome where placements spread across chromosomes.
func TestByteParity_ChromFirstSmall(t *testing.T) {
	g, order := loadGenomeOrdered(t, "small.genome")
	got := runShuffleBytes(t, "small.bed", Options{Genome: g, GenomeOrder: order, ChromFirst: true, Seed: 9})
	assertGolden(t, got, "expected.small_chromFirst.seed9.bed")
}

// TestByteParity_IncludeExclude covers the combined -incl -excl path on a
// non-degenerate small fixture.
func TestByteParity_IncludeExclude(t *testing.T) {
	g, order := loadGenomeOrdered(t, "small.genome")
	incl := readParityBED(t, "small.incl")
	excl := readParityBED(t, "small.excl")
	got := runShuffleBytes(t, "small.bed", Options{
		Genome: g, GenomeOrder: order, Include: incl, Exclude: excl, Seed: 7,
	})
	assertGolden(t, got, "expected.small_inclexcl.seed7.bed")
}

// TestByteParity_AllowBeyondChromEnd covers -allowBeyondChromEnd (clamp to the
// chromosome end instead of redrawing).
func TestByteParity_AllowBeyondChromEnd(t *testing.T) {
	g, order := loadGenomeOrdered(t, "small.genome")
	got := runShuffleBytes(t, "small.bed", Options{
		Genome: g, GenomeOrder: order, AllowBeyondChromEnd: true, Seed: 3,
	})
	assertGolden(t, got, "expected.small_allowbeyond.seed3.bed")
}

// itoa is a tiny int64→string helper (avoids importing strconv solely for
// test subtest names alongside the existing helpers in parity_test.go).
func itoa(n int64) string {
	if n == 0 {
		return "0"
	}
	neg := n < 0
	if neg {
		n = -n
	}
	var buf [20]byte
	i := len(buf)
	for n > 0 {
		i--
		buf[i] = byte('0' + n%10)
		n /= 10
	}
	if neg {
		i--
		buf[i] = '-'
	}
	return string(buf[i:])
}
