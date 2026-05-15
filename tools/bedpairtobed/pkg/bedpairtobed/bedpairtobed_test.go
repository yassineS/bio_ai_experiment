package bedpairtobed

import (
	"bytes"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/bed"
)

// helper to run with string inputs and return output.
func runStrings(t *testing.T, a, b string, opts Options) (string, *Stats) {
	t.Helper()
	var out bytes.Buffer
	s, err := Run(strings.NewReader(a), strings.NewReader(b), &out, opts)
	if err != nil {
		t.Fatalf("Run: %v", err)
	}
	return out.String(), s
}

func TestEither_TwoHitsOneOnEachEnd(t *testing.T) {
	a := "chr1\t10\t20\tchr2\t100\t200\tp1\t30\t+\t-\n"
	b := strings.Join([]string{
		"chr1\t5\t15\tb1\t0\t+",
		"chr2\t150\t250\tb2\t0\t-",
		"chr3\t0\t10\tb3\t0\t+", // no overlap
		"",
	}, "\n")
	got, stats := runStrings(t, a, b, Options{Type: TypeEither})
	if stats.OutputLines != 2 || stats.APairs != 1 || stats.BRecords != 3 || stats.EmittedPairs != 1 {
		t.Errorf("stats: %+v", stats)
	}
	want := "chr1\t10\t20\tchr2\t100\t200\tp1\t30\t+\t-\tchr1\t5\t15\tb1\t0\t+\n" +
		"chr1\t10\t20\tchr2\t100\t200\tp1\t30\t+\t-\tchr2\t150\t250\tb2\t0\t-\n"
	if got != want {
		t.Errorf("output mismatch:\n got=%q\nwant=%q", got, want)
	}
}

func TestBoth_RequiresHitOnEachEnd(t *testing.T) {
	a := strings.Join([]string{
		"chr1\t10\t20\tchr2\t100\t200\tboth\t.\t+\t-",
		"chr1\t10\t20\tchr2\t300\t400\tone-end-only\t.\t+\t-", // only end1 hits
		"",
	}, "\n")
	b := strings.Join([]string{
		"chr1\t0\t30\thit\t0\t+",
		"chr2\t100\t150\thit\t0\t-",
		"",
	}, "\n")
	got, stats := runStrings(t, a, b, Options{Type: TypeBoth})
	if stats.OutputLines != 2 || stats.EmittedPairs != 1 {
		t.Errorf("stats: %+v", stats)
	}
	if !strings.Contains(got, "both") || strings.Contains(got, "one-end-only") {
		t.Errorf("unexpected output: %q", got)
	}
}

func TestNeither_NoHits(t *testing.T) {
	a := "chr1\t10\t20\tchr2\t100\t200\tlonely\t.\t+\t-\n"
	b := "chr9\t0\t10\thit\t0\t+\n"
	got, _ := runStrings(t, a, b, Options{Type: TypeNeither})
	want := "chr1\t10\t20\tchr2\t100\t200\tlonely\t.\t+\t-\n"
	if got != want {
		t.Errorf("neither: %q want %q", got, want)
	}
}

func TestXor_OneEndOnly(t *testing.T) {
	a := strings.Join([]string{
		"chr1\t10\t20\tchr5\t100\t200\tone\t.\t+\t-",  // end2 on chr5 has no B
		"chr1\t10\t20\tchr2\t100\t200\tboth\t.\t+\t-", // both ends hit
		"",
	}, "\n")
	b := strings.Join([]string{
		"chr1\t0\t30\thit\t0\t+",
		"chr2\t150\t250\thit\t0\t-",
		"",
	}, "\n")
	got, stats := runStrings(t, a, b, Options{Type: TypeXor})
	if stats.OutputLines != 1 || stats.EmittedPairs != 1 {
		t.Errorf("stats: %+v want 1 output 1 emitted", stats)
	}
	if !strings.Contains(got, "\tone\t") || strings.Contains(got, "\tboth\t") {
		t.Errorf("xor output unexpected: %q", got)
	}
}

func TestNotboth_OnlyOneEndHits(t *testing.T) {
	a := strings.Join([]string{
		"chr1\t10\t20\tchr2\t100\t200\tboth\t.\t+\t-",
		"chr1\t10\t20\tchr5\t100\t200\tone\t.\t+\t-", // end2 on chr5 has no B
		"chr3\t10\t20\tchr4\t100\t200\tnone\t.\t+\t-",
		"",
	}, "\n")
	b := strings.Join([]string{
		"chr1\t0\t30\thit\t0\t+",
		"chr2\t150\t250\thit\t0\t-",
		"",
	}, "\n")
	got, stats := runStrings(t, a, b, Options{Type: TypeNotboth})
	// Expect: "one" with a single end1-hit line, and "none" emitted bare. "both" suppressed.
	if stats.EmittedPairs != 2 {
		t.Errorf("expected 2 emitted pairs, got %+v", stats)
	}
	if strings.Contains(got, "\tboth\t") {
		t.Errorf("notboth should not emit fully-paired record: %q", got)
	}
	if !strings.Contains(got, "\tone\t") || !strings.Contains(got, "\tnone\t") {
		t.Errorf("notboth missing expected pairs: %q", got)
	}
}

func TestFraction_Filters(t *testing.T) {
	// A end1 is 100bp. A B-record overlapping only 10bp -> fraction 0.10.
	a := "chr1\t0\t100\tchr1\t1000\t1100\tp\t.\t+\t+\n"
	b := "chr1\t90\t110\tb\t0\t+\n" // covers bp 90..99 -> 10bp -> 0.10 fraction
	// With default 1e-9 this should hit; with f=0.5 it should miss.
	_, st1 := runStrings(t, a, b, Options{Type: TypeEither})
	if st1.OutputLines != 1 {
		t.Errorf("default f: expected 1 hit, got %+v", st1)
	}
	_, st2 := runStrings(t, a, b, Options{Type: TypeEither, MinFractionA: 0.5})
	if st2.OutputLines != 0 {
		t.Errorf("f=0.5: expected 0 hits, got %+v", st2)
	}
}

func TestStrand_SameOpposite(t *testing.T) {
	a := "chr1\t10\t20\tchr2\t100\t200\tp\t.\t+\t-\n"
	bSamePos := "chr1\t0\t30\thit\t0\t+\n"
	bSameNeg := "chr1\t0\t30\thit\t0\t-\n"

	_, sSame := runStrings(t, a, bSamePos, Options{Type: TypeEither, SameStrand: true})
	if sSame.OutputLines != 1 {
		t.Errorf("same-strand + matching: %+v", sSame)
	}
	_, sMis := runStrings(t, a, bSameNeg, Options{Type: TypeEither, SameStrand: true})
	if sMis.OutputLines != 0 {
		t.Errorf("same-strand + mismatching: %+v", sMis)
	}
	_, sOpp := runStrings(t, a, bSameNeg, Options{Type: TypeEither, OppositeStrand: true})
	if sOpp.OutputLines != 1 {
		t.Errorf("opposite-strand + matching: %+v", sOpp)
	}
	// Ignore overrides both.
	_, sIg := runStrings(t, a, bSameNeg, Options{Type: TypeEither, SameStrand: true, IgnoreStrand: true})
	if sIg.OutputLines != 1 {
		t.Errorf("ignore-strand override: %+v", sIg)
	}
}

func TestInvalidType_Error(t *testing.T) {
	_, err := Run(strings.NewReader(""), strings.NewReader(""), &bytes.Buffer{}, Options{Type: "bogus"})
	if err == nil {
		t.Fatal("expected error for unknown type")
	}
}

func TestSameAndOppositeMutuallyExclusive(t *testing.T) {
	_, err := Run(strings.NewReader(""), strings.NewReader(""), &bytes.Buffer{}, Options{SameStrand: true, OppositeStrand: true})
	if err == nil {
		t.Fatal("expected -s/-S exclusivity error")
	}
}

func TestUnalignedEnd_NoHits(t *testing.T) {
	// Unaligned end (chrom ".") returns no hits.
	a := "chr1\t10\t20\t.\t-1\t-1\tp\t.\t+\t.\n"
	b := "chr1\t0\t30\thit\t0\t+\n"
	got, st := runStrings(t, a, b, Options{Type: TypeEither})
	if st.OutputLines != 1 {
		t.Errorf("unaligned: expected only end1 hit, got %d (%q)", st.OutputLines, got)
	}
}

func TestNotxor_BothOrNeither(t *testing.T) {
	a := strings.Join([]string{
		"chr1\t10\t20\tchr2\t100\t200\tboth\t.\t+\t-",
		"chr1\t10\t20\tchr5\t100\t200\tone\t.\t+\t-",  // end1 hits, end2 misses
		"chr9\t10\t20\tchr9\t100\t200\tnone\t.\t+\t-", // misses both
		"",
	}, "\n")
	b := strings.Join([]string{
		"chr1\t0\t30\thit\t0\t+",
		"chr2\t150\t250\thit\t0\t-",
		"",
	}, "\n")
	got, st := runStrings(t, a, b, Options{Type: TypeNotxor})
	if st.EmittedPairs != 2 {
		t.Errorf("expected 2 emitted (both + none), got %+v", st)
	}
	if !strings.Contains(got, "\tboth\t") || !strings.Contains(got, "\tnone\t") {
		t.Errorf("notxor missing expected pairs: %q", got)
	}
	if strings.Contains(got, "\tone\t") {
		t.Errorf("notxor should not emit xor pair: %q", got)
	}
}

func TestNotxor_NeitherHits_EmitsBare(t *testing.T) {
	a := "chr9\t10\t20\tchr9\t100\t200\tnone\t.\t+\t-\n"
	b := "chr1\t0\t30\thit\t0\t+\n"
	got, st := runStrings(t, a, b, Options{Type: TypeNotxor})
	if st.OutputLines != 1 || !strings.HasPrefix(got, "chr9\t10\t20") || strings.Count(got, "\t") != 9 {
		t.Errorf("notxor neither: %q (%+v)", got, st)
	}
}

func TestEither_NoHits_NoOutput(t *testing.T) {
	a := "chr9\t10\t20\tchr9\t100\t200\tn\t.\t+\t-\n"
	b := "chr1\t0\t30\thit\t0\t+\n"
	got, st := runStrings(t, a, b, Options{Type: TypeEither})
	if got != "" || st.OutputLines != 0 || st.EmittedPairs != 0 {
		t.Errorf("either no-hits: %q (%+v)", got, st)
	}
}

// errWriter fails after producing some bytes; used to exercise error paths.
type errWriter struct{ remaining int }

func (e *errWriter) Write(p []byte) (int, error) {
	if e.remaining <= 0 {
		return 0, errFakeWrite
	}
	n := len(p)
	if n > e.remaining {
		n = e.remaining
	}
	e.remaining -= n
	if e.remaining <= 0 {
		return n, errFakeWrite
	}
	return n, nil
}

var errFakeWrite = errFake{}

type errFake struct{}

func (errFake) Error() string { return "fake write error" }

func TestEmit_WriteError_Propagates(t *testing.T) {
	// Force a write error to exercise the error-return paths inside emit.
	a := "chr1\t10\t20\tchr2\t100\t200\tp\t.\t+\t-\n"
	b := "chr1\t0\t30\thit\t0\t+\nchr2\t150\t250\thit\t0\t-\n"
	_, err := Run(strings.NewReader(a), strings.NewReader(b), &errWriter{remaining: 0}, Options{Type: TypeEither})
	if err == nil {
		t.Fatal("expected error from failing writer")
	}
}

func TestStrand_DotAndEmpty(t *testing.T) {
	// Strand checks should not reject when either side has "" or "." strand.
	a := "chr1\t10\t20\tchr2\t100\t200\tp\t.\t.\t.\n"
	b := "chr1\t0\t30\thit\t0\t+\n"
	_, st := runStrings(t, a, b, Options{Type: TypeEither, SameStrand: true})
	if st.OutputLines != 1 {
		t.Errorf("dot-strand wildcard: %+v", st)
	}
	_, st = runStrings(t, a, b, Options{Type: TypeEither, OppositeStrand: true})
	if st.OutputLines != 1 {
		t.Errorf("dot-strand wildcard (opposite): %+v", st)
	}
}

func TestBadBEDPE_ReturnsError(t *testing.T) {
	// 9 columns -> reader returns error.
	bad := "chr1\t1\t2\tchr1\t3\t4\tn\t.\t+\n"
	_, err := Run(strings.NewReader(bad), strings.NewReader(""), &bytes.Buffer{}, Options{})
	if err == nil {
		t.Fatal("expected error from malformed BEDPE")
	}
}

func TestFormatBED_MinimalAndFull(t *testing.T) {
	bare := formatBED(&bed.Record{Chrom: "chr1", ChromStart: 1, ChromEnd: 2})
	if bare != "chr1\t1\t2" {
		t.Errorf("bare: %q", bare)
	}
	full := formatBED(&bed.Record{
		Chrom: "chr1", ChromStart: 1, ChromEnd: 2,
		Name: "n", Score: 7, Strand: "+",
		ExtraFields: []string{"x"},
	})
	if full != "chr1\t1\t2\tn\t7\t+\tx" {
		t.Errorf("full: %q", full)
	}
}
