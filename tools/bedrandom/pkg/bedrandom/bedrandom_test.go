package bedrandom

import (
	"bytes"
	"strconv"
	"strings"
	"testing"
)

func mustGenome(t *testing.T, text string) *Genome {
	t.Helper()
	g, err := ParseGenome(strings.NewReader(text))
	if err != nil {
		t.Fatalf("ParseGenome: %v", err)
	}
	return g
}

// TestParseGenome_OrderAndOffsets verifies file order, cumulative offsets, and
// total size are computed exactly as upstream GenomeFile does.
func TestParseGenome_OrderAndOffsets(t *testing.T) {
	g := mustGenome(t, "chr1\t1000\nchr2\t2000\nchr3\t500\n")
	if got, want := g.NumChroms(), 3; got != want {
		t.Fatalf("NumChroms = %d, want %d", got, want)
	}
	if got, want := g.Size(), int64(3500); got != want {
		t.Fatalf("Size = %d, want %d", got, want)
	}
	wantOffsets := []int64{0, 1000, 3000}
	for i, w := range wantOffsets {
		if g.offsets[i] != w {
			t.Errorf("offsets[%d] = %d, want %d", i, g.offsets[i], w)
		}
	}
}

// TestParseGenome_SkipsCommentsBlanksAndBadSizes covers the parser's skip rules.
func TestParseGenome_SkipsCommentsBlanksAndBadSizes(t *testing.T) {
	g := mustGenome(t, "# comment\n\nchr1\t100\nchrX\tnotanumber\nchr2\t50\n")
	if g.NumChroms() != 2 {
		t.Fatalf("NumChroms = %d, want 2 (comment/blank/bad-size skipped)", g.NumChroms())
	}
	if g.names[0] != "chr1" || g.names[1] != "chr2" {
		t.Errorf("names = %v, want [chr1 chr2]", g.names)
	}
}

func TestParseGenome_TooFewFields(t *testing.T) {
	if _, err := ParseGenome(strings.NewReader("chr1\n")); err == nil {
		t.Fatal("expected error for line with one field")
	}
}

// TestProjectOnGenome reproduces upstream's lower_bound projection: a genome
// offset lands on the chromosome whose [offset, offset+size) range contains it.
func TestProjectOnGenome(t *testing.T) {
	g := mustGenome(t, "chr1\t1000\nchr2\t2000\nchr3\t500\n")
	cases := []struct {
		pos       int64
		wantChrom string
		wantStart int64
	}{
		{0, "chr1", 0},
		{999, "chr1", 999},
		{1000, "chr2", 0},
		{2999, "chr2", 1999},
		{3000, "chr3", 0},
		{3499, "chr3", 499},
	}
	for _, tc := range cases {
		chrom, start := g.projectOnGenome(tc.pos)
		if chrom != tc.wantChrom || start != tc.wantStart {
			t.Errorf("projectOnGenome(%d) = (%s,%d), want (%s,%d)",
				tc.pos, chrom, start, tc.wantChrom, tc.wantStart)
		}
	}
}

// TestGenerate_Geometry checks that every generated record is well-formed BED6:
// length-L interval inside a real chromosome, 1..n index, length in the score
// column, and a + or - strand.
func TestGenerate_Geometry(t *testing.T) {
	g := mustGenome(t, "chr1\t1000\nchr2\t2000\nchr3\t500\n")
	const n, l = 200, 100
	var buf bytes.Buffer
	written, err := Generate(g, &buf, Options{Length: l, Num: n, Seed: 42, HaveSeed: true})
	if err != nil {
		t.Fatalf("Generate: %v", err)
	}
	if written != n {
		t.Fatalf("written = %d, want %d", written, n)
	}
	lines := strings.Split(strings.TrimRight(buf.String(), "\n"), "\n")
	if len(lines) != n {
		t.Fatalf("line count = %d, want %d", len(lines), n)
	}
	sizes := map[string]int64{"chr1": 1000, "chr2": 2000, "chr3": 500}
	for i, ln := range lines {
		f := strings.Split(ln, "\t")
		if len(f) != 6 {
			t.Fatalf("row %d: %d columns, want 6: %q", i, len(f), ln)
		}
		size, ok := sizes[f[0]]
		if !ok {
			t.Errorf("row %d: chrom %q not in genome", i, f[0])
			continue
		}
		start, _ := strconv.ParseInt(f[1], 10, 64)
		end, _ := strconv.ParseInt(f[2], 10, 64)
		idx, _ := strconv.Atoi(f[3])
		score, _ := strconv.ParseInt(f[4], 10, 64)
		if end-start != l {
			t.Errorf("row %d: length = %d, want %d", i, end-start, l)
		}
		if start < 0 || end > size {
			t.Errorf("row %d: interval %d-%d outside chrom %q size %d", i, start, end, f[0], size)
		}
		if idx != i+1 {
			t.Errorf("row %d: index = %d, want %d", i, idx, i+1)
		}
		if score != l {
			t.Errorf("row %d: score = %d, want %d", i, score, l)
		}
		if f[5] != "+" && f[5] != "-" {
			t.Errorf("row %d: strand = %q, want + or -", i, f[5])
		}
	}
}

// TestGenerate_Deterministic confirms the same seed yields identical output and
// a different seed (very likely) differs.
func TestGenerate_Deterministic(t *testing.T) {
	g := mustGenome(t, "chr1\t1000\nchr2\t2000\n")
	run := func(seed int) string {
		var b bytes.Buffer
		if _, err := Generate(g, &b, Options{Length: 50, Num: 50, Seed: seed, HaveSeed: true}); err != nil {
			t.Fatalf("Generate: %v", err)
		}
		return b.String()
	}
	a, b := run(42), run(42)
	if a != b {
		t.Fatal("same seed produced different output")
	}
	if run(42) == run(99) {
		t.Fatal("different seeds produced identical output")
	}
}

// TestGenerate_KnownSequence pins the exact bytes for a small fixed seed so the
// engine + draw order can't silently drift even without the upstream binary.
// The golden was captured from upstream bedtools v2.31.1 (-seed 42).
func TestGenerate_KnownSequence(t *testing.T) {
	g := mustGenome(t, "chr1\t1000\nchr2\t2000\nchr3\t500\n")
	var buf bytes.Buffer
	if _, err := Generate(g, &buf, Options{Length: 100, Num: 5, Seed: 42, HaveSeed: true}); err != nil {
		t.Fatalf("Generate: %v", err)
	}
	want := "chr2\t406\t506\t1\t100\t-\n" +
		"chr1\t662\t762\t2\t100\t+\n" +
		"chr2\t1428\t1528\t3\t100\t-\n" +
		"chr1\t644\t744\t4\t100\t-\n" +
		"chr2\t1257\t1357\t5\t100\t+\n"
	if buf.String() != want {
		t.Fatalf("KnownSequence mismatch:\n got: %q\nwant: %q", buf.String(), want)
	}
}

func TestGenerate_Errors(t *testing.T) {
	g := mustGenome(t, "chr1\t1000\n")
	if _, err := Generate(g, &bytes.Buffer{}, Options{Length: 0, Num: 1}); err == nil {
		t.Error("expected error for non-positive length")
	}
	if _, err := Generate(g, &bytes.Buffer{}, Options{Length: 100, Num: -1}); err == nil {
		t.Error("expected error for negative num")
	}
}
