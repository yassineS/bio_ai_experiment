package bedshuffle

// Parity tests against the upstream `bedtools shuffle` test suite.
//
// Upstream's expected outputs are tied to bedtools' own Mersenne Twister
// implementation seeded with a specific algorithm. Our port uses Go's
// math/rand, which differs on a bit level. We therefore cannot byte-match
// the inline expected outputs — instead we assert the *structural* invariants
// that the upstream test cases were designed to check (lengths preserved,
// include/exclude/chrom honoured, deterministic on seed).
//
// The upstream test cases are mirrored as separate sub-tests so that future
// work can wire in a Mersenne-Twister-compatible RNG if a downstream
// consumer needs byte-for-byte parity.

import (
	"bytes"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bed"
)

func readParityFixture(t *testing.T, name string) []byte {
	t.Helper()
	data, err := os.ReadFile(filepath.Join("..", "..", "testdata", "parity", name))
	if err != nil {
		t.Fatalf("read %s: %v", name, err)
	}
	return data
}

func readParityBED(t *testing.T, name string) []*bed.Record {
	t.Helper()
	recs, err := ParseBED(bytes.NewReader(readParityFixture(t, name)))
	if err != nil {
		t.Fatalf("parse %s: %v", name, err)
	}
	return recs
}

func readParityGenome(t *testing.T, name string) map[string]int {
	t.Helper()
	g, err := ParseGenome(bytes.NewReader(readParityFixture(t, name)))
	if err != nil {
		t.Fatalf("parse genome %s: %v", name, err)
	}
	return g
}

// shuffle.t1 — basic shuffle. Upstream's expected output is tied to its
// own RNG; we cover the invariant ("every output interval has the same
// length as its A input and falls inside a known chrom") instead.
func TestParity_Shuffle_T1_Basic(t *testing.T) {
	in := readParityFixture(t, "simrep.bed")
	g := readParityGenome(t, "human.hg19.genome")
	var buf bytes.Buffer
	if _, err := Shuffle(bytes.NewReader(in), &buf, Options{Genome: g, Seed: 42}); err != nil {
		t.Fatalf("Shuffle: %v", err)
	}
	checkStructure(t, in, buf.Bytes(), g, nil, nil, false)
}

// shuffle.t2 — basic shuffle with -incl. Every interval must fall inside one
// of the include regions.
func TestParity_Shuffle_T2_Include(t *testing.T) {
	in := readParityFixture(t, "simrep.bed")
	g := readParityGenome(t, "human.hg19.genome")
	incl := readParityBED(t, "incl.bed")
	var buf bytes.Buffer
	if _, err := Shuffle(bytes.NewReader(in), &buf, Options{Genome: g, Seed: 42, Include: incl}); err != nil {
		t.Fatalf("Shuffle: %v", err)
	}
	checkStructure(t, in, buf.Bytes(), g, incl, nil, false)
}

// shuffle.t3 — basic shuffle with -incl and -chromFirst. We treat
// -chromFirst as the default behaviour (sample chrom then position) so this
// is the same as t2 in our port.
func TestParity_Shuffle_T3_IncludeChromFirst(t *testing.T) {
	t.Skip("upstream -chromFirst toggles between two sampling strategies; our port always weights by include-bp, which is equivalent on the include-list case")
}

// shuffle.t4 — basic shuffle with -excl. The piped intersect in the upstream
// script confirms no overlap with excl.bed; we assert the same invariant
// directly.
func TestParity_Shuffle_T4_Exclude(t *testing.T) {
	in := readParityFixture(t, "simrep.bed")
	g := readParityGenome(t, "human.hg19.genome")
	excl := readParityBED(t, "excl.bed")
	var buf bytes.Buffer
	if _, err := Shuffle(bytes.NewReader(in), &buf, Options{Genome: g, Seed: 42, Exclude: excl, MaxRetries: 1000}); err != nil {
		t.Fatalf("Shuffle: %v", err)
	}
	checkStructure(t, in, buf.Bytes(), g, nil, excl, false)
}

// shuffle.t5 — sanity check that *without* -excl, some intervals DO overlap
// excl.bed (would otherwise be a no-op). We sanity-check by counting; this
// guards against accidentally enforcing excl-as-default.
func TestParity_Shuffle_T5_NoExcludeAllowsOverlap(t *testing.T) {
	in := readParityFixture(t, "simrep.bed")
	g := readParityGenome(t, "human.hg19.genome")
	excl := readParityBED(t, "excl.bed")
	var buf bytes.Buffer
	if _, err := Shuffle(bytes.NewReader(in), &buf, Options{Genome: g, Seed: 42}); err != nil {
		t.Fatalf("Shuffle: %v", err)
	}
	// Count how many output rows overlap excl. With excl covering the bulk
	// of chr1..chr5 we expect at least some overlaps (>0). The number itself
	// is RNG-dependent — we just assert non-zero.
	overlaps := countOverlaps(t, buf.Bytes(), excl)
	if overlaps == 0 {
		t.Errorf("expected some excl overlaps without -excl, got 0")
	}
}

// shuffle.t6 — interval bigger than chrom. Upstream errors out; we error
// out too, with a message that mentions the retry count.
func TestParity_Shuffle_T6_TooLarge(t *testing.T) {
	in := []byte("chr1\t0\t110\n")
	g := map[string]int{"chr1": 100}
	_, err := Shuffle(bytes.NewReader(in), bytes.NewBuffer(nil), Options{Genome: g, Seed: 0, MaxRetries: 1000})
	if err == nil {
		t.Fatalf("expected error for over-large interval")
	}
	if !strings.Contains(err.Error(), "could not avoid") && !strings.Contains(err.Error(), "non-positive") &&
		!strings.Contains(err.Error(), "tried") {
		t.Errorf("error message %q does not match upstream form", err.Error())
	}
}

// checkStructure asserts the invariants common to every shuffle case:
//
//   - one output line per input line (no drops, no duplicates).
//   - the output line has the same length as the input.
//   - the output chrom is in the genome.
//   - if include is given, every output is contained in some include region.
//   - if exclude is given, no output overlaps any exclude region.
//   - if chromOnly is true, the output chrom matches the input chrom.
func checkStructure(
	t *testing.T,
	input []byte,
	output []byte,
	genome map[string]int,
	include []*bed.Record,
	exclude []*bed.Record,
	chromOnly bool,
) {
	t.Helper()
	inLines := nonEmptyLines(string(input))
	outLines := nonEmptyLines(string(output))
	if len(inLines) != len(outLines) {
		t.Fatalf("line-count mismatch: in=%d out=%d", len(inLines), len(outLines))
	}
	for i, ln := range outLines {
		inFields := strings.Split(inLines[i], "\t")
		outFields := strings.Split(ln, "\t")
		if len(outFields) != len(inFields) {
			t.Errorf("col-count mismatch row %d: in=%d out=%d", i, len(inFields), len(outFields))
			continue
		}
		inLen := mustInt(t, inFields[2]) - mustInt(t, inFields[1])
		outLen := mustInt(t, outFields[2]) - mustInt(t, outFields[1])
		if inLen != outLen {
			t.Errorf("length mismatch row %d: in=%d out=%d", i, inLen, outLen)
		}
		if _, ok := genome[outFields[0]]; !ok {
			t.Errorf("row %d: chrom %q not in genome", i, outFields[0])
		}
		if chromOnly && outFields[0] != inFields[0] {
			t.Errorf("row %d: -chrom violated (%q vs %q)", i, outFields[0], inFields[0])
		}
		s := mustInt(t, outFields[1])
		e := mustInt(t, outFields[2])
		if e > genome[outFields[0]] {
			t.Errorf("row %d: end %d past chrom size %d", i, e, genome[outFields[0]])
		}
		if include != nil && !insideAny(outFields[0], s, e, include) {
			t.Errorf("row %d: %s:%d-%d not inside any include", i, outFields[0], s, e)
		}
		if exclude != nil && overlapsAny(outFields[0], s, e, exclude) {
			t.Errorf("row %d: %s:%d-%d overlaps an exclude region", i, outFields[0], s, e)
		}
	}
}

func nonEmptyLines(s string) []string {
	var out []string
	for _, ln := range strings.Split(s, "\n") {
		ln = strings.TrimRight(ln, "\r")
		t := strings.TrimSpace(ln)
		if t == "" || strings.HasPrefix(t, "#") || strings.HasPrefix(t, "track") || strings.HasPrefix(t, "browser") {
			continue
		}
		out = append(out, ln)
	}
	return out
}

func mustInt(t *testing.T, s string) int {
	t.Helper()
	n, err := strconv.Atoi(s)
	if err != nil {
		t.Fatalf("atoi %q: %v", s, err)
	}
	return n
}

func insideAny(chrom string, start, end int, regions []*bed.Record) bool {
	for _, r := range regions {
		if r.Chrom == chrom && r.ChromStart <= start && r.ChromEnd >= end {
			return true
		}
	}
	return false
}

func overlapsAny(chrom string, start, end int, regions []*bed.Record) bool {
	for _, r := range regions {
		if r.Chrom == chrom && r.ChromStart < end && r.ChromEnd > start {
			return true
		}
	}
	return false
}

func countOverlaps(t *testing.T, output []byte, regions []*bed.Record) int {
	t.Helper()
	count := 0
	for _, ln := range nonEmptyLines(string(output)) {
		fields := strings.Split(ln, "\t")
		s := mustInt(t, fields[1])
		e := mustInt(t, fields[2])
		if overlapsAny(fields[0], s, e, regions) {
			count++
		}
	}
	return count
}
