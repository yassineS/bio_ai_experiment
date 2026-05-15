package bedmultiinter

import (
	"bytes"
	"io"
	"strings"
	"testing"
)

func readers(in ...string) []io.Reader {
	out := make([]io.Reader, len(in))
	for i, s := range in {
		out[i] = strings.NewReader(s)
	}
	return out
}

// Hand-computed: 2-file overlap on a single chromosome.
//
//	A:  [0,10)  [20,30)
//	B:        [5,25)
//
// Event sweep:
//
//	0..5  -> {A}
//	5..10 -> {A,B}
//	10..20 -> {B}
//	20..25 -> {A,B}
//	25..30 -> {A}
func TestRun_TwoFileBasic(t *testing.T) {
	a := "chr1\t0\t10\nchr1\t20\t30\n"
	b := "chr1\t5\t25\n"
	var out bytes.Buffer
	if _, err := Run(readers(a, b), &out, Options{
		Filenames: []string{"a.bed", "b.bed"},
	}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	want := "chr1\t0\t5\t1\t1\t1\t0\n" +
		"chr1\t5\t10\t2\t1,2\t1\t1\n" +
		"chr1\t10\t20\t1\t2\t0\t1\n" +
		"chr1\t20\t25\t2\t1,2\t1\t1\n" +
		"chr1\t25\t30\t1\t1\t1\t0\n"
	if got := out.String(); got != want {
		t.Fatalf("output mismatch:\n got:\n%s\nwant:\n%s", got, want)
	}
}

// Hand-computed: the example from bedtools multiintersect_examples().
// Three files, eight expected rows. Mirror of the upstream docstring.
func TestRun_UpstreamExample(t *testing.T) {
	a := "chr1\t6\t12\nchr1\t10\t20\nchr1\t22\t27\nchr1\t24\t30\n"
	b := "chr1\t12\t32\nchr1\t14\t30\n"
	c := "chr1\t8\t15\nchr1\t10\t14\nchr1\t32\t34\n"
	var out bytes.Buffer
	if _, err := Run(readers(a, b, c), &out, Options{
		Filenames: []string{"a.bed", "b.bed", "c.bed"},
	}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	want := "chr1\t6\t8\t1\t1\t1\t0\t0\n" +
		"chr1\t8\t12\t2\t1,3\t1\t0\t1\n" +
		"chr1\t12\t15\t3\t1,2,3\t1\t1\t1\n" +
		"chr1\t15\t20\t2\t1,2\t1\t1\t0\n" +
		"chr1\t20\t22\t1\t2\t0\t1\t0\n" +
		"chr1\t22\t30\t2\t1,2\t1\t1\t0\n" +
		"chr1\t30\t32\t1\t2\t0\t1\t0\n" +
		"chr1\t32\t34\t1\t3\t0\t0\t1\n"
	if got := out.String(); got != want {
		t.Fatalf("output mismatch:\n got:\n%s\nwant:\n%s", got, want)
	}
}

// Hand-computed: -names and -header.
func TestRun_NamesAndHeader(t *testing.T) {
	a := "chr1\t0\t10\n"
	b := "chr1\t5\t15\n"
	var out bytes.Buffer
	if _, err := Run(readers(a, b), &out, Options{
		Filenames: []string{"a.bed", "b.bed"},
		Names:     []string{"alpha", "beta"},
		Header:    true,
	}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	want := "chrom\tstart\tend\tnum\tlist\talpha\tbeta\n" +
		"chr1\t0\t5\t1\talpha\t1\t0\n" +
		"chr1\t5\t10\t2\talpha,beta\t1\t1\n" +
		"chr1\t10\t15\t1\tbeta\t0\t1\n"
	if got := out.String(); got != want {
		t.Fatalf("output mismatch:\n got:\n%s\nwant:\n%s", got, want)
	}
}

// Hand-computed: -empty + -g emits the leading and trailing 0-count gaps,
// plus any internal gaps between covered spans.
func TestRun_EmptyWithGenome(t *testing.T) {
	a := "chr1\t10\t20\n"
	b := "chr1\t30\t40\n"
	sizes := map[string]int{"chr1": 100}
	var out bytes.Buffer
	if _, err := Run(readers(a, b), &out, Options{
		Filenames:  []string{"a.bed", "b.bed"},
		Empty:      true,
		ChromSizes: sizes,
	}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	want := "chr1\t0\t10\t0\tnone\t0\t0\n" +
		"chr1\t10\t20\t1\t1\t1\t0\n" +
		"chr1\t20\t30\t0\tnone\t0\t0\n" +
		"chr1\t30\t40\t1\t2\t0\t1\n" +
		"chr1\t40\t100\t0\tnone\t0\t0\n"
	if got := out.String(); got != want {
		t.Fatalf("output mismatch:\n got:\n%s\nwant:\n%s", got, want)
	}
}

// -cluster should collapse adjacent same-set rows.
func TestRun_Cluster(t *testing.T) {
	// Two A intervals that abut at coord 10 → without cluster, two rows;
	// with cluster, one row [0,20) with active={A}.
	a := "chr1\t0\t10\nchr1\t10\t20\n"
	b := "chr1\t50\t60\n"
	var out bytes.Buffer
	if _, err := Run(readers(a, b), &out, Options{
		Filenames: []string{"a.bed", "b.bed"},
		Cluster:   true,
	}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	want := "chr1\t0\t20\t1\t1\t1\t0\n" +
		"chr1\t50\t60\t1\t2\t0\t1\n"
	if got := out.String(); got != want {
		t.Fatalf("output mismatch:\n got:\n%s\nwant:\n%s", got, want)
	}
}

// -filler N/A swaps the "0" indicator for arbitrary text.
func TestRun_CustomFiller(t *testing.T) {
	a := "chr1\t0\t10\n"
	b := "chr1\t20\t30\n"
	var out bytes.Buffer
	if _, err := Run(readers(a, b), &out, Options{
		Filenames: []string{"a", "b"},
		Filler:    "N/A",
	}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	want := "chr1\t0\t10\t1\t1\t1\tN/A\n" +
		"chr1\t20\t30\t1\t2\tN/A\t1\n"
	if got := out.String(); got != want {
		t.Fatalf("output mismatch:\n got:\n%s\nwant:\n%s", got, want)
	}
}

// Within-file overlapping intervals must be merged before sweeping; they
// would otherwise inflate the depth count at the overlap.
func TestRun_WithinFileMerge(t *testing.T) {
	a := "chr1\t0\t10\nchr1\t5\t15\n" // overlaps -> merged to [0,15)
	b := "chr1\t10\t20\n"
	var out bytes.Buffer
	if _, err := Run(readers(a, b), &out, Options{
		Filenames: []string{"a", "b"},
	}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	want := "chr1\t0\t10\t1\t1\t1\t0\n" +
		"chr1\t10\t15\t2\t1,2\t1\t1\n" +
		"chr1\t15\t20\t1\t2\t0\t1\n"
	if got := out.String(); got != want {
		t.Fatalf("merge mismatch:\n got:\n%s\nwant:\n%s", got, want)
	}
}

// Multi-chrom: rows are emitted in sorted chrom order, with per-chrom
// sweeps independent.
func TestRun_MultiChrom(t *testing.T) {
	a := "chr1\t0\t10\nchr2\t0\t10\n"
	b := "chr1\t5\t15\nchr2\t5\t15\n"
	var out bytes.Buffer
	if _, err := Run(readers(a, b), &out, Options{
		Filenames: []string{"a", "b"},
	}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	want := "chr1\t0\t5\t1\t1\t1\t0\n" +
		"chr1\t5\t10\t2\t1,2\t1\t1\n" +
		"chr1\t10\t15\t1\t2\t0\t1\n" +
		"chr2\t0\t5\t1\t1\t1\t0\n" +
		"chr2\t5\t10\t2\t1,2\t1\t1\n" +
		"chr2\t10\t15\t1\t2\t0\t1\n"
	if got := out.String(); got != want {
		t.Fatalf("multichrom mismatch:\n got:\n%s\nwant:\n%s", got, want)
	}
}

// Argument validation.
func TestRun_Errors(t *testing.T) {
	var out bytes.Buffer
	if _, err := Run(readers("chr1\t0\t10\n"), &out, Options{Filenames: []string{"a"}}); err == nil {
		t.Fatal("expected error for <2 inputs")
	}
	if _, err := Run(readers("chr1\t0\t10\n", "chr1\t0\t10\n"), &out,
		Options{Filenames: []string{"a"}}); err == nil {
		t.Fatal("expected error for Filenames count mismatch")
	}
	if _, err := Run(readers("chr1\t0\t10\n", "chr1\t0\t10\n"), &out,
		Options{Filenames: []string{"a", "b"}, Names: []string{"A"}}); err == nil {
		t.Fatal("expected error for Names count mismatch")
	}
	if _, err := Run(readers("chr1\t0\t10\n", "chr1\t0\t10\n"), &out,
		Options{Filenames: []string{"a", "b"}, Empty: true}); err == nil {
		t.Fatal("expected error for -empty without sizes")
	}
}

// Header without -names should fall back to the supplied filenames.
func TestRun_HeaderWithoutNames(t *testing.T) {
	a := "chr1\t0\t10\n"
	b := "chr1\t20\t30\n"
	var out bytes.Buffer
	if _, err := Run(readers(a, b), &out, Options{
		Filenames: []string{"a.bed", "b.bed"},
		Header:    true,
	}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	want := "chrom\tstart\tend\tnum\tlist\ta.bed\tb.bed\n" +
		"chr1\t0\t10\t1\t1\t1\t0\n" +
		"chr1\t20\t30\t1\t2\t0\t1\n"
	if got := out.String(); got != want {
		t.Fatalf("default-header mismatch:\n got:\n%s\nwant:\n%s", got, want)
	}
}

// ReadGenomeSizes happy + error path.
func TestReadGenomeSizes(t *testing.T) {
	g, err := ReadGenomeSizes(strings.NewReader("# comment\nchr1\t1000\nchr2\t2000\n"))
	if err != nil {
		t.Fatalf("ReadGenomeSizes: %v", err)
	}
	if g["chr1"] != 1000 || g["chr2"] != 2000 {
		t.Fatalf("unexpected sizes: %v", g)
	}
	if _, err := ReadGenomeSizes(strings.NewReader("chr1\tNAN\n")); err == nil {
		t.Fatal("expected error on non-int size")
	}
	if _, err := ReadGenomeSizes(strings.NewReader("chr1_no_size\n")); err == nil {
		t.Fatal("expected error on short line")
	}
}

// DefaultNames strips the dir prefix and trailing extension, matching
// upstream `stl_basename`.
func TestDefaultNames(t *testing.T) {
	names := DefaultNames([]string{"a.bed", "/tmp/foo/b.bed.gz", "c"})
	want := []string{"a", "b.bed", "c"}
	for i := range names {
		if names[i] != want[i] {
			t.Fatalf("names[%d]=%q want %q", i, names[i], want[i])
		}
	}
}
