package bedpairtopair

import (
	"bytes"
	"strings"
	"testing"
)

func runStrings(t *testing.T, a, b string, opts Options) (string, *Stats) {
	t.Helper()
	var out bytes.Buffer
	s, err := Run(strings.NewReader(a), strings.NewReader(b), &out, opts)
	if err != nil {
		t.Fatalf("Run: %v", err)
	}
	return out.String(), s
}

// Default mode = both: emits when SAME B pair has one end matched by A1 and
// the other matched by A2.
func TestBoth_BasicHit(t *testing.T) {
	a := "chr1\t10\t20\tchr1\t100\t200\tap\t.\t+\t-\n"
	b := "chr1\t5\t25\tchr1\t150\t250\tbp\t.\t+\t-\n"
	got, st := runStrings(t, a, b, Options{Type: TypeBoth})
	if st.OutputLines != 1 || st.EmittedPairs != 1 {
		t.Errorf("stats: %+v", st)
	}
	if !strings.Contains(got, "\tap\t") || !strings.Contains(got, "\tbp\t") {
		t.Errorf("missing pair labels: %q", got)
	}
}

// Both mode should also work in the swapped orientation: A1 hits B2 and A2 hits B1.
func TestBoth_SwappedOrientation(t *testing.T) {
	a := "chr1\t10\t20\tchr1\t100\t200\tap\t.\t+\t-\n"
	b := "chr1\t150\t250\tchr1\t5\t25\tbp\t.\t-\t+\n"
	got, st := runStrings(t, a, b, Options{Type: TypeBoth})
	if st.OutputLines != 1 {
		t.Errorf("swapped: expected 1 output, got %+v (%q)", st, got)
	}
}

// Both mode rejects when only one end matches.
func TestBoth_OnlyOneEndMatches(t *testing.T) {
	a := "chr1\t10\t20\tchr1\t100\t200\tap\t.\t+\t-\n"
	// B's end1 matches A's end1, but B's end2 is far away.
	b := "chr1\t5\t25\tchr1\t9000\t9100\tbp\t.\t+\t-\n"
	_, st := runStrings(t, a, b, Options{Type: TypeBoth})
	if st.OutputLines != 0 {
		t.Errorf("expected 0 emissions, got %+v", st)
	}
}

// Either mode: any end-vs-end overlap counts.
func TestEither_AnyEndMatches(t *testing.T) {
	a := "chr1\t10\t20\tchr1\t100\t200\tap\t.\t+\t-\n"
	b := "chr1\t5\t25\tchr1\t9000\t9100\tbp\t.\t+\t-\n"
	got, st := runStrings(t, a, b, Options{Type: TypeEither})
	if st.OutputLines != 1 || !strings.Contains(got, "\tap\t") {
		t.Errorf("either single-end hit: %+v %q", st, got)
	}
}

// Neither mode: emits when no end matches anything.
func TestNeither_NoMatches(t *testing.T) {
	a := "chr9\t10\t20\tchr9\t100\t200\tap\t.\t+\t-\n"
	b := "chr1\t5\t25\tchr1\t150\t250\tbp\t.\t+\t-\n"
	got, st := runStrings(t, a, b, Options{Type: TypeNeither})
	if st.OutputLines != 1 || !strings.HasPrefix(got, "chr9\t10\t20") {
		t.Errorf("neither: %+v %q", st, got)
	}
}

// Notboth mode: emits when no B pair matches both A ends.
func TestNotboth_OnlySingleEndMatches(t *testing.T) {
	a := strings.Join([]string{
		"chr1\t10\t20\tchr1\t100\t200\tboth\t.\t+\t-",
		"chr1\t10\t20\tchr1\t9000\t9100\tone\t.\t+\t-",
		"",
	}, "\n")
	b := "chr1\t5\t25\tchr1\t150\t250\tbp\t.\t+\t-\n"
	got, _ := runStrings(t, a, b, Options{Type: TypeNotboth})
	if strings.Contains(got, "\tboth\t") {
		t.Errorf("notboth: should not emit fully-matched pair: %q", got)
	}
	if !strings.Contains(got, "\tone\t") {
		t.Errorf("notboth: expected emission of single-end pair: %q", got)
	}
}

// Slop extends each A end before computing overlaps.
func TestSlop_ExtendsRange(t *testing.T) {
	// A end1 just barely misses B end1: A1 is chr1:90-100, B1 is chr1:120-130.
	// Slop of 30 makes A1 effectively chr1:60-130 -> overlap.
	a := "chr1\t90\t100\tchr1\t900\t1000\tap\t.\t+\t+\n"
	b := "chr1\t120\t130\tchr1\t950\t1050\tbp\t.\t+\t+\n"
	_, st := runStrings(t, a, b, Options{Type: TypeBoth, Slop: 0})
	if st.OutputLines != 0 {
		t.Errorf("no-slop: expected 0, got %+v", st)
	}
	_, st = runStrings(t, a, b, Options{Type: TypeBoth, Slop: 50})
	if st.OutputLines != 1 {
		t.Errorf("with slop 50: expected 1, got %+v", st)
	}
}

func TestStrandedSlop_Direction(t *testing.T) {
	// A1 is chr1:90-100 on +. B1 is chr1:120-130. Stranded slop +50 only extends end1.
	a := "chr1\t90\t100\tchr1\t900\t1000\tap\t.\t+\t+\n"
	b := "chr1\t120\t130\tchr1\t950\t1050\tbp\t.\t+\t+\n"
	_, st := runStrings(t, a, b, Options{Type: TypeBoth, Slop: 50, StrandedSlop: true})
	if st.OutputLines != 1 {
		t.Errorf("stranded slop +: expected 1, got %+v", st)
	}
	// Now A1 is on minus strand: stranded slop only extends start. B1 is to A's right -> miss.
	a = "chr1\t90\t100\tchr1\t900\t1000\tap\t.\t-\t+\n"
	_, st = runStrings(t, a, b, Options{Type: TypeBoth, Slop: 50, StrandedSlop: true, IgnoreStrand: true})
	if st.OutputLines != 0 {
		t.Errorf("stranded slop - on right-side B: expected 0, got %+v", st)
	}
}

func TestRDN_FilterSelfHits(t *testing.T) {
	a := "chr1\t10\t20\tchr1\t100\t200\tshared\t.\t+\t-\n"
	b := "chr1\t5\t25\tchr1\t150\t250\tshared\t.\t+\t-\n"
	_, st := runStrings(t, a, b, Options{Type: TypeBoth})
	if st.OutputLines != 1 {
		t.Errorf("default RDN off: expected 1, got %+v", st)
	}
	_, st = runStrings(t, a, b, Options{Type: TypeBoth, RequireDifferentNames: true})
	if st.OutputLines != 0 {
		t.Errorf("rdn on, same name: expected 0, got %+v", st)
	}
}

func TestIgnoreStrand_OverridesDefault(t *testing.T) {
	// Default behaviour enforces matching strands.
	a := "chr1\t10\t20\tchr1\t100\t200\tap\t.\t+\t+\n"
	b := "chr1\t5\t25\tchr1\t150\t250\tbp\t.\t-\t-\n"
	_, st := runStrings(t, a, b, Options{Type: TypeBoth})
	if st.OutputLines != 0 {
		t.Errorf("default: expected 0, got %+v", st)
	}
	_, st = runStrings(t, a, b, Options{Type: TypeBoth, IgnoreStrand: true})
	if st.OutputLines != 1 {
		t.Errorf("-is: expected 1, got %+v", st)
	}
}

func TestFraction_Filters(t *testing.T) {
	// A end1 = 100bp, B end1 covers only 5bp -> fraction 0.05.
	a := "chr1\t0\t100\tchr1\t1000\t1100\tap\t.\t+\t+\n"
	b := "chr1\t95\t105\tchr1\t1000\t1100\tbp\t.\t+\t+\n"
	_, st := runStrings(t, a, b, Options{Type: TypeBoth, MinFraction: 0.01})
	if st.OutputLines != 1 {
		t.Errorf("f=0.01: expected 1, got %+v", st)
	}
	_, st = runStrings(t, a, b, Options{Type: TypeBoth, MinFraction: 0.5})
	if st.OutputLines != 0 {
		t.Errorf("f=0.5: expected 0, got %+v", st)
	}
}

func TestInvalidType(t *testing.T) {
	_, err := Run(strings.NewReader(""), strings.NewReader(""), &bytes.Buffer{}, Options{Type: "bogus"})
	if err == nil {
		t.Fatal("expected invalid type error")
	}
}

func TestSameOpposite_MutuallyExclusive(t *testing.T) {
	_, err := Run(strings.NewReader(""), strings.NewReader(""), &bytes.Buffer{}, Options{SameStrand: true, OppositeStrand: true})
	if err == nil {
		t.Fatal("expected -s/-S exclusivity error")
	}
}

func TestSSWithoutSlop(t *testing.T) {
	_, err := Run(strings.NewReader(""), strings.NewReader(""), &bytes.Buffer{}, Options{StrandedSlop: true})
	if err == nil {
		t.Fatal("expected -ss-without-slop error")
	}
}

func TestBadBEDPE(t *testing.T) {
	_, err := Run(strings.NewReader("chr1\tbad\t2\tchr1\t3\t4\tn\t.\t+\t-\n"),
		strings.NewReader(""), &bytes.Buffer{}, Options{})
	if err == nil {
		t.Fatal("expected error from bad BEDPE A")
	}
	_, err = Run(strings.NewReader(""),
		strings.NewReader("chr1\tbad\t2\tchr1\t3\t4\tn\t.\t+\t-\n"), &bytes.Buffer{}, Options{})
	if err == nil {
		t.Fatal("expected error from bad BEDPE B")
	}
}

func TestEither_BothEndsMatchSameB_OnceEmitted(t *testing.T) {
	a := "chr1\t10\t20\tchr1\t100\t200\tap\t.\t+\t-\n"
	b := "chr1\t5\t25\tchr1\t150\t250\tbp\t.\t+\t-\n"
	_, st := runStrings(t, a, b, Options{Type: TypeEither})
	if st.OutputLines != 1 {
		t.Errorf("either both-ends-same-B: expected 1, got %+v", st)
	}
}

func TestSameStrandOption(t *testing.T) {
	// Default enforces same-strand; -s/--same-strand is functionally redundant
	// but should not break anything.
	a := "chr1\t10\t20\tchr1\t100\t200\tap\t.\t+\t-\n"
	b := "chr1\t5\t25\tchr1\t150\t250\tbp\t.\t+\t-\n"
	_, st := runStrings(t, a, b, Options{Type: TypeBoth, SameStrand: true})
	if st.OutputLines != 1 {
		t.Errorf("-s on matching: %+v", st)
	}
	_, st = runStrings(t, a, b, Options{Type: TypeBoth, OppositeStrand: true})
	if st.OutputLines != 0 {
		t.Errorf("-S on matching strands: expected 0 hits, got %+v", st)
	}
}

func TestUnalignedEnd_Skipped(t *testing.T) {
	// Unaligned ends never produce hits.
	a := "chr1\t10\t20\t.\t-1\t-1\tap\t.\t+\t.\n"
	b := "chr1\t5\t25\tchr1\t150\t250\tbp\t.\t+\t-\n"
	_, st := runStrings(t, a, b, Options{Type: TypeEither})
	if st.OutputLines != 1 {
		t.Errorf("either with unaligned end2: expected 1 (end1 only), got %+v", st)
	}
}
