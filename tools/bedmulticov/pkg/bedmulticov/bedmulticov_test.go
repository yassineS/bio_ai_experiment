package bedmulticov

import (
	"bytes"
	"io"
	"strings"
	"testing"
)

// readers turns N input strings into N io.Readers — convenient for tests.
func readers(in ...string) []io.Reader {
	out := make([]io.Reader, len(in))
	for i, s := range in {
		out[i] = strings.NewReader(s)
	}
	return out
}

// Hand-computed: A has 2 intervals on chr1; two B files. B1 has 2 overlaps
// with A.row1 and 1 with A.row2; B2 has 0 with A.row1 and 2 with A.row2.
func TestRun_BasicTwoFiles(t *testing.T) {
	a := "chr1\t0\t100\tA1\n" +
		"chr1\t200\t300\tA2\n"
	b1 := "chr1\t10\t20\n" +
		"chr1\t50\t60\n" +
		"chr1\t250\t260\n"
	b2 := "chr1\t210\t220\n" +
		"chr1\t290\t295\n"
	var out bytes.Buffer
	n, err := Run(strings.NewReader(a), readers(b1, b2), &out, Options{})
	if err != nil {
		t.Fatalf("Run: %v", err)
	}
	if n != 2 {
		t.Fatalf("expected 2 A records, got %d", n)
	}
	want := "chr1\t0\t100\tA1\t2\t0\n" +
		"chr1\t200\t300\tA2\t1\t2\n"
	if got := out.String(); got != want {
		t.Fatalf("output mismatch:\n got:\n%s\nwant:\n%s", got, want)
	}
}

// Hand-computed: strand filter -s drops opposite-strand B records. A is
// BED6 ('+'). B1 has one + and one - overlap; the - should not count.
// With -S the result inverts.
func TestRun_StrandFilters(t *testing.T) {
	a := "chr1\t100\t200\tA1\t0\t+\n"
	b1 := "chr1\t110\t120\t.\t0\t+\n" +
		"chr1\t150\t160\t.\t0\t-\n"
	var out bytes.Buffer
	if _, err := Run(strings.NewReader(a), readers(b1), &out, Options{SameStrand: true}); err != nil {
		t.Fatalf("Run -s: %v", err)
	}
	if got, want := out.String(), "chr1\t100\t200\tA1\t0\t+\t1\n"; got != want {
		t.Fatalf("-s mismatch:\n got %q\nwant %q", got, want)
	}

	out.Reset()
	if _, err := Run(strings.NewReader(a), readers(b1), &out, Options{OppositeStrand: true}); err != nil {
		t.Fatalf("Run -S: %v", err)
	}
	if got, want := out.String(), "chr1\t100\t200\tA1\t0\t+\t1\n"; got != want {
		t.Fatalf("-S mismatch:\n got %q\nwant %q", got, want)
	}
}

// Hand-computed: -f 0.5 requires >= 50% of A covered by a single B record.
// A is 100bp; B1 has a 40bp hit (fails) and a 60bp hit (passes); B2 has
// two 20bp hits (each fails). Reciprocal: -r adds the same threshold on B.
func TestRun_FractionAndReciprocal(t *testing.T) {
	a := "chr1\t0\t100\n"
	b1 := "chr1\t0\t40\n" + // 40% of A
		"chr1\t10\t70\n" // 60% of A — passes -f 0.5
	b2 := "chr1\t0\t20\n" + // 20% of A
		"chr1\t30\t50\n" // 20% of A
	var out bytes.Buffer
	if _, err := Run(strings.NewReader(a), readers(b1, b2), &out, Options{FractionA: 0.5}); err != nil {
		t.Fatalf("Run -f: %v", err)
	}
	if got, want := out.String(), "chr1\t0\t100\t1\t0\n"; got != want {
		t.Fatalf("-f mismatch:\n got %q\nwant %q", got, want)
	}

	// -r with -f 0.5 also requires 50% of B covered: the 60bp B hit
	// extends from 10..70 and A overlap is 60bp = 100% of B, so it still
	// passes; but with -f 0.7 the 60% A coverage now fails entirely.
	out.Reset()
	if _, err := Run(strings.NewReader(a), readers(b1, b2), &out,
		Options{FractionA: 0.7, Reciprocal: true}); err != nil {
		t.Fatalf("Run -r: %v", err)
	}
	if got, want := out.String(), "chr1\t0\t100\t0\t0\n"; got != want {
		t.Fatalf("-r mismatch:\n got %q\nwant %q", got, want)
	}
}

// Conflicting strand flags should error out.
func TestRun_ConflictingStrandFlags(t *testing.T) {
	var out bytes.Buffer
	_, err := Run(strings.NewReader(""), nil, &out,
		Options{SameStrand: true, OppositeStrand: true})
	if err == nil {
		t.Fatal("expected error for -s + -S, got nil")
	}
}

// -r without -f is a user error.
func TestRun_ReciprocalWithoutFraction(t *testing.T) {
	var out bytes.Buffer
	_, err := Run(strings.NewReader(""), nil, &out, Options{Reciprocal: true})
	if err == nil {
		t.Fatal("expected error for -r without -f, got nil")
	}
}

// Range checks on -f / -F should surface user errors.
func TestRun_FractionRangeValidation(t *testing.T) {
	var out bytes.Buffer
	if _, err := Run(strings.NewReader(""), nil, &out, Options{FractionA: 1.5}); err == nil {
		t.Fatal("expected -f range error")
	}
	out.Reset()
	if _, err := Run(strings.NewReader(""), nil, &out, Options{FractionB: -0.1}); err == nil {
		t.Fatal("expected -F range error")
	}
}

// BED comments / track / browser lines and short rows should be handled
// gracefully — the former skipped, the latter surfaced as an error.
func TestRun_CommentsAndMalformed(t *testing.T) {
	a := "# header line\n" +
		"track name=foo\n" +
		"browser hide all\n" +
		"\n" +
		"chr1\t0\t10\n"
	var out bytes.Buffer
	if _, err := Run(strings.NewReader(a), nil, &out, Options{}); err != nil {
		t.Fatalf("Run on comments: %v", err)
	}
	if got, want := out.String(), "chr1\t0\t10\n"; got != want {
		t.Fatalf("comment handling mismatch:\n got %q\nwant %q", got, want)
	}

	// Malformed: 2 columns.
	out.Reset()
	if _, err := Run(strings.NewReader("chr1\t0\n"), nil, &out, Options{}); err == nil {
		t.Fatal("expected error on short record")
	}
	// Malformed: non-int start.
	out.Reset()
	if _, err := Run(strings.NewReader("chr1\tBAD\t10\n"), nil, &out, Options{}); err == nil {
		t.Fatal("expected error on non-int start")
	}
	// Malformed: non-int end.
	out.Reset()
	if _, err := Run(strings.NewReader("chr1\t0\tBAD\n"), nil, &out, Options{}); err == nil {
		t.Fatal("expected error on non-int end")
	}
}

// strand "" on either side under -s / -S should be treated as "no match".
func TestRun_MissingStrandUnderFilter(t *testing.T) {
	a := "chr1\t0\t100\tA1\t0\t+\n"
	b := "chr1\t10\t20\n" // no strand column
	var out bytes.Buffer
	if _, err := Run(strings.NewReader(a), readers(b), &out, Options{SameStrand: true}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	if got, want := out.String(), "chr1\t0\t100\tA1\t0\t+\t0\n"; got != want {
		t.Fatalf("missing-strand mismatch:\n got %q\nwant %q", got, want)
	}
	out.Reset()
	if _, err := Run(strings.NewReader(a), readers(b), &out, Options{OppositeStrand: true}); err != nil {
		t.Fatalf("Run -S: %v", err)
	}
	if got, want := out.String(), "chr1\t0\t100\tA1\t0\t+\t0\n"; got != want {
		t.Fatalf("missing-strand -S mismatch:\n got %q\nwant %q", got, want)
	}
}

// Multi-chrom inputs and intervals that span chrom gaps should still
// report a 0 for the file that has no records on a given chrom.
func TestRun_MultiChromAndMissingChrom(t *testing.T) {
	a := "chr1\t0\t100\nchr2\t0\t100\n"
	b1 := "chr1\t10\t20\n" // chr1 only
	b2 := "chr2\t10\t20\n" // chr2 only
	var out bytes.Buffer
	if _, err := Run(strings.NewReader(a), readers(b1, b2), &out, Options{}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	want := "chr1\t0\t100\t1\t0\n" +
		"chr2\t0\t100\t0\t1\n"
	if got := out.String(); got != want {
		t.Fatalf("output mismatch:\n got:\n%s\nwant:\n%s", got, want)
	}
}
