package bedshuffle

import (
	"bytes"
	"io"
	"strconv"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bed"
)

// run is a tiny helper that runs Shuffle on the given inputs and returns
// the output lines.
func run(t *testing.T, input string, opts Options) []string {
	t.Helper()
	var buf bytes.Buffer
	if _, err := Shuffle(strings.NewReader(input), &buf, opts); err != nil {
		t.Fatalf("Shuffle: %v", err)
	}
	out := strings.TrimRight(buf.String(), "\n")
	if out == "" {
		return nil
	}
	return strings.Split(out, "\n")
}

func TestShuffle_BasicDeterministic(t *testing.T) {
	in := "chr1\t100\t200\nchr2\t10\t20\n"
	opts := Options{
		Genome: map[string]int{"chr1": 1000, "chr2": 500, "chr3": 800},
		Seed:   42,
	}
	a := run(t, in, opts)
	b := run(t, in, opts)
	if len(a) != 2 || len(b) != 2 {
		t.Fatalf("expected 2 lines each, got %d / %d", len(a), len(b))
	}
	for i := range a {
		if a[i] != b[i] {
			t.Errorf("seed %d not deterministic: %q vs %q", opts.Seed, a[i], b[i])
		}
	}
}

func TestShuffle_DifferentSeedsDiffer(t *testing.T) {
	in := "chr1\t100\t200\n"
	g := map[string]int{"chr1": 1000, "chr2": 1000, "chr3": 1000}
	got1 := run(t, in, Options{Genome: g, Seed: 1})
	got2 := run(t, in, Options{Genome: g, Seed: 2})
	if got1[0] == got2[0] {
		t.Errorf("expected different output for different seeds, both %q", got1[0])
	}
}

func TestShuffle_LengthPreserved(t *testing.T) {
	in := "chr1\t100\t250\n"
	out := run(t, in, Options{Genome: map[string]int{"chr1": 100000}, Seed: 7})
	fields := strings.Split(out[0], "\t")
	if fields[0] != "chr1" {
		t.Errorf("chrom changed: %q", fields[0])
	}
	start := parseInt(t, fields[1])
	end := parseInt(t, fields[2])
	if end-start != 150 {
		t.Errorf("length not preserved: %d", end-start)
	}
}

func TestShuffle_ChromKeepsOriginal(t *testing.T) {
	in := "chr1\t100\t200\nchr2\t50\t150\nchr3\t10\t90\n"
	g := map[string]int{"chr1": 1000, "chr2": 1000, "chr3": 1000}
	out := run(t, in, Options{Genome: g, Seed: 11, Chrom: true})
	expectChrom := []string{"chr1", "chr2", "chr3"}
	for i, line := range out {
		if !strings.HasPrefix(line, expectChrom[i]+"\t") {
			t.Errorf("chrom not preserved at row %d: %q", i, line)
		}
	}
}

func TestShuffle_IncludeRestrictsToRegions(t *testing.T) {
	in := "chr1\t10\t30\n"
	g := map[string]int{"chr1": 1000}
	incl := []*bed.Record{{Chrom: "chr1", ChromStart: 100, ChromEnd: 200}}
	out := run(t, in, Options{Genome: g, Seed: 5, Include: incl})
	fields := strings.Split(out[0], "\t")
	s := parseInt(t, fields[1])
	e := parseInt(t, fields[2])
	if s < 100 || e > 200 {
		t.Errorf("placement %d-%d outside include [100,200)", s, e)
	}
}

func TestShuffle_ExcludeBlocksRegion(t *testing.T) {
	in := "chr1\t0\t10\n"
	g := map[string]int{"chr1": 1000}
	excl := []*bed.Record{{Chrom: "chr1", ChromStart: 0, ChromEnd: 990}}
	// Only [990, 1000) is allowed (10bp).
	out := run(t, in, Options{Genome: g, Seed: 5, Exclude: excl})
	fields := strings.Split(out[0], "\t")
	s := parseInt(t, fields[1])
	e := parseInt(t, fields[2])
	if s < 990 || e > 1000 {
		t.Errorf("placement %d-%d violates exclude", s, e)
	}
}

func TestShuffle_PreservesExtraColumns(t *testing.T) {
	in := "chr1\t100\t200\tname\t500\t+\textra1\textra2\n"
	g := map[string]int{"chr1": 10000}
	out := run(t, in, Options{Genome: g, Seed: 1})
	fields := strings.Split(out[0], "\t")
	if len(fields) != 8 {
		t.Errorf("expected 8 cols preserved, got %d: %v", len(fields), fields)
	}
	if fields[3] != "name" || fields[6] != "extra1" || fields[7] != "extra2" {
		t.Errorf("extras lost: %v", fields)
	}
}

func TestShuffle_HeadersAndCommentsSkipped(t *testing.T) {
	in := "#comment\ntrack name=foo\nchr1\t10\t20\n"
	out := run(t, in, Options{Genome: map[string]int{"chr1": 1000}, Seed: 0})
	if len(out) != 1 {
		t.Errorf("expected 1 data row, got %d: %v", len(out), out)
	}
}

func TestShuffle_NoGenomeError(t *testing.T) {
	_, err := Shuffle(strings.NewReader("chr1\t0\t10\n"), io.Discard, Options{})
	if err == nil {
		t.Errorf("expected error for missing genome")
	}
}

func TestShuffle_FailsAfterMaxRetries(t *testing.T) {
	// A 10bp interval, but the only chrom is 100bp AND everything is
	// excluded ⇒ no valid placement; should error after maxTries.
	in := "chr1\t0\t10\n"
	excl := []*bed.Record{{Chrom: "chr1", ChromStart: 0, ChromEnd: 100}}
	_, err := Shuffle(strings.NewReader(in), io.Discard, Options{
		Genome:     map[string]int{"chr1": 100},
		Exclude:    excl,
		Seed:       1,
		MaxRetries: 50,
	})
	if err == nil {
		t.Errorf("expected unplaceable error")
	} else if !strings.Contains(err.Error(), "tried 50 potential loci") {
		t.Errorf("unexpected error: %v", err)
	}
}

func TestShuffle_IntervalLargerThanChrom(t *testing.T) {
	in := "chr1\t0\t110\n"
	_, err := Shuffle(strings.NewReader(in), io.Discard, Options{
		Genome:     map[string]int{"chr1": 100},
		Seed:       0,
		MaxRetries: 100,
	})
	if err == nil {
		t.Errorf("expected error for over-large interval")
	}
}

func TestShuffle_InvalidBED(t *testing.T) {
	_, err := Shuffle(strings.NewReader("chr1\tbad\t10\n"), io.Discard, Options{
		Genome: map[string]int{"chr1": 1000},
	})
	if err == nil {
		t.Errorf("expected parse error")
	}
	_, err = Shuffle(strings.NewReader("chr1\t10\tbad\n"), io.Discard, Options{
		Genome: map[string]int{"chr1": 1000},
	})
	if err == nil {
		t.Errorf("expected parse error on end")
	}
	_, err = Shuffle(strings.NewReader("chr1\t10\t10\n"), io.Discard, Options{
		Genome: map[string]int{"chr1": 1000},
	})
	if err == nil {
		t.Errorf("expected non-positive length error")
	}
	_, err = Shuffle(strings.NewReader("chr1\t10\n"), io.Discard, Options{
		Genome: map[string]int{"chr1": 1000},
	})
	if err == nil {
		t.Errorf("expected too-few-fields error")
	}
}

func TestShuffle_ChromMissingFromGenome(t *testing.T) {
	// Input on chrX but genome only knows chr1. With -chrom we should fail.
	in := "chrX\t0\t10\n"
	_, err := Shuffle(strings.NewReader(in), io.Discard, Options{
		Genome:     map[string]int{"chr1": 1000},
		Chrom:      true,
		Seed:       0,
		MaxRetries: 50,
	})
	if err == nil {
		t.Errorf("expected unplaceable error")
	}
}

func TestShuffle_IncludeChromFiltersChroms(t *testing.T) {
	// Include only has chr2 regions; with default mode the sampler should
	// only ever pick chr2.
	in := "chr1\t10\t30\n"
	g := map[string]int{"chr1": 1000, "chr2": 1000}
	incl := []*bed.Record{{Chrom: "chr2", ChromStart: 100, ChromEnd: 500}}
	for i := 1; i <= 5; i++ {
		out := run(t, in, Options{Genome: g, Seed: int64(i), Include: incl})
		if !strings.HasPrefix(out[0], "chr2\t") {
			t.Errorf("seed %d landed off chr2: %q", i, out[0])
		}
	}
}

func TestParseGenome(t *testing.T) {
	g, err := ParseGenome(strings.NewReader("#comment\nchr1\t1000\nchr2\t2000\n\n"))
	if err != nil {
		t.Fatalf("ParseGenome: %v", err)
	}
	if g["chr1"] != 1000 || g["chr2"] != 2000 {
		t.Errorf("ParseGenome bad: %v", g)
	}
	_, err = ParseGenome(strings.NewReader("chr1\n"))
	if err == nil {
		t.Errorf("expected error for too few fields")
	}
	_, err = ParseGenome(strings.NewReader("chr1\tNaN\n"))
	if err == nil {
		t.Errorf("expected error for non-numeric size")
	}
}

func TestParseBED(t *testing.T) {
	recs, err := ParseBED(strings.NewReader("chr1\t10\t20\nchr2\t30\t40\n"))
	if err != nil {
		t.Fatalf("ParseBED: %v", err)
	}
	if len(recs) != 2 {
		t.Errorf("expected 2 records, got %d", len(recs))
	}
}

// parseInt is a tiny test helper.
func parseInt(t *testing.T, s string) int {
	t.Helper()
	n, err := strconv.Atoi(s)
	if err != nil {
		t.Fatalf("parseInt %q: %v", s, err)
	}
	return n
}
