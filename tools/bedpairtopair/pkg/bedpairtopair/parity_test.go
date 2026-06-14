package bedpairtopair

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// loadFixture reads a fixture file from tools/bedpairtopair/testdata/parity/.
func loadFixture(t *testing.T, name string) string {
	t.Helper()
	wd, err := os.Getwd()
	if err != nil {
		t.Fatalf("getwd: %v", err)
	}
	p := filepath.Join(wd, "..", "..", "testdata", "parity", name)
	b, err := os.ReadFile(p)
	if err != nil {
		t.Fatalf("read fixture %s: %v", p, err)
	}
	return string(b)
}

// Upstream `bedtools test/` ships no pairtopair subdirectory. The expected
// outputs below were computed by hand from the C++ source semantics and
// cross-checked against the upstream binary by hand. Each `Parity_*` test
// asserts one output mode.

func TestParity_Both(t *testing.T) {
	a := loadFixture(t, "a.bedpe")
	b := loadFixture(t, "b.bedpe")
	var out bytes.Buffer
	if _, err := Run(strings.NewReader(a), strings.NewReader(b), &out, Options{Type: TypeBoth}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	want := strings.Join([]string{
		"chr1\t10\t20\tchr1\t100\t200\ta1\t0\t+\t-\tchr1\t5\t25\tchr1\t150\t250\tb1\t0\t+\t-",
		"chr1\t500\t600\tchr2\t1000\t1100\ta2\t0\t+\t+\tchr1\t500\t600\tchr2\t1050\t1150\tb2\t0\t+\t+",
		"chr2\t2000\t2100\tchr2\t3000\t3100\ta3\t0\t-\t-\tchr2\t2050\t2150\tchr2\t3000\t3100\tb3\t0\t-\t-",
		"",
	}, "\n")
	if got := out.String(); got != want {
		t.Errorf("both:\n got=%q\nwant=%q", got, want)
	}
}

func TestParity_Either(t *testing.T) {
	a := loadFixture(t, "a.bedpe")
	b := loadFixture(t, "b.bedpe")
	var out bytes.Buffer
	if _, err := Run(strings.NewReader(a), strings.NewReader(b), &out, Options{Type: TypeEither}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	want := strings.Join([]string{
		"chr1\t10\t20\tchr1\t100\t200\ta1\t0\t+\t-\tchr1\t5\t25\tchr1\t150\t250\tb1\t0\t+\t-",
		"chr1\t10\t20\tchr1\t100\t200\ta1\t0\t+\t-\tchr1\t5\t25\tchr1\t900\t1000\tb4\t0\t+\t-",
		"chr1\t500\t600\tchr2\t1000\t1100\ta2\t0\t+\t+\tchr1\t500\t600\tchr2\t1050\t1150\tb2\t0\t+\t+",
		"chr2\t2000\t2100\tchr2\t3000\t3100\ta3\t0\t-\t-\tchr2\t2050\t2150\tchr2\t3000\t3100\tb3\t0\t-\t-",
		"",
	}, "\n")
	if got := out.String(); got != want {
		t.Errorf("either:\n got=%q\nwant=%q", got, want)
	}
}

func TestParity_Neither(t *testing.T) {
	a := loadFixture(t, "a.bedpe")
	b := loadFixture(t, "b.bedpe")
	var out bytes.Buffer
	if _, err := Run(strings.NewReader(a), strings.NewReader(b), &out, Options{Type: TypeNeither}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	want := "chr9\t0\t10\tchr9\t100\t200\tlonely\t0\t+\t+\n"
	if got := out.String(); got != want {
		t.Errorf("neither:\n got=%q\nwant=%q", got, want)
	}
}

func TestParity_Notboth(t *testing.T) {
	a := loadFixture(t, "a.bedpe")
	b := loadFixture(t, "b.bedpe")
	var out bytes.Buffer
	if _, err := Run(strings.NewReader(a), strings.NewReader(b), &out, Options{Type: TypeNotboth}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	// Pairs whose both ends both match SOME B pair on both ends are suppressed.
	// a1, a2, a3 each have a B with both-end match -> suppressed. lonely is emitted.
	want := "chr9\t0\t10\tchr9\t100\t200\tlonely\t0\t+\t+\n"
	if got := out.String(); got != want {
		t.Errorf("notboth:\n got=%q\nwant=%q", got, want)
	}
}

// TestParity_Slop documents that adding -slop widens the search window
// symmetrically. The fixture has a deliberate near-miss between a candidate
// B-pair and A; slop closes the gap.
func TestParity_Slop(t *testing.T) {
	a := "chr1\t100\t110\tchr1\t1000\t1010\tap\t0\t+\t+\n"
	b := "chr1\t150\t160\tchr1\t950\t960\tbp\t0\t+\t+\n"
	// Without slop: no overlap on either end (A1 vs B1: 100-110 vs 150-160 — gap; A2 vs B2: 1000-1010 vs 950-960 — gap).
	var out bytes.Buffer
	if _, err := Run(strings.NewReader(a), strings.NewReader(b), &out, Options{Type: TypeBoth}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	if out.Len() != 0 {
		t.Errorf("no-slop expected empty, got %q", out.String())
	}
	out.Reset()
	// With slop 100, A becomes effectively 0..210 / 900..1110 — overlaps both B ends.
	if _, err := Run(strings.NewReader(a), strings.NewReader(b), &out, Options{Type: TypeBoth, Slop: 100}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	if out.Len() == 0 {
		t.Errorf("slop=100 expected hit, got empty")
	}
}

// TestParity_StrandedSlop asserts -ss (strand-aware slop) byte-for-byte
// against the upstream `bedtools pairtopair` binary, using a fixture chosen so
// strand direction actually changes the result: A's second footprint is
// identical to B's (always overlaps), while B's first footprint sits just
// BEFORE A's start. With a "+" strand, -ss extends only the end, so the
// leading gap is never bridged and there is NO match; plain symmetric -slop
// extends the start too, bridging the gap and producing a match. The expected
// outputs (ss_slop10_ss.expected.bedpe = empty, ss_slop10_sym.expected.bedpe =
// the matched pair) were generated from the upstream binary.
func TestParity_StrandedSlop(t *testing.T) {
	a := loadFixture(t, "ss_a.bedpe")
	b := loadFixture(t, "ss_b.bedpe")

	var ssOut bytes.Buffer
	if _, err := Run(strings.NewReader(a), strings.NewReader(b), &ssOut,
		Options{Type: TypeBoth, Slop: 10, StrandedSlop: true}); err != nil {
		t.Fatalf("Run -ss: %v", err)
	}
	if got, want := ssOut.String(), loadFixture(t, "ss_slop10_ss.expected.bedpe"); got != want {
		t.Errorf("-ss:\n got=%q\nwant=%q", got, want)
	}

	var symOut bytes.Buffer
	if _, err := Run(strings.NewReader(a), strings.NewReader(b), &symOut,
		Options{Type: TypeBoth, Slop: 10}); err != nil {
		t.Fatalf("Run symmetric: %v", err)
	}
	if got, want := symOut.String(), loadFixture(t, "ss_slop10_sym.expected.bedpe"); got != want {
		t.Errorf("symmetric:\n got=%q\nwant=%q", got, want)
	}
}
