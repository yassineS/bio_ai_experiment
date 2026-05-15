package bedpairtobed

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// loadFixture reads a fixture file from tools/bedpairtobed/testdata/parity/.
// It walks up from the test's working directory until it finds the file —
// callers don't have to know the repo layout.
func loadFixture(t *testing.T, name string) string {
	t.Helper()
	wd, err := os.Getwd()
	if err != nil {
		t.Fatalf("getwd: %v", err)
	}
	// We're at tools/bedpairtobed/pkg/bedpairtobed during go test.
	p := filepath.Join(wd, "..", "..", "testdata", "parity", name)
	b, err := os.ReadFile(p)
	if err != nil {
		t.Fatalf("read fixture %s: %v", p, err)
	}
	return string(b)
}

// Upstream `bedtools test/` ships no pairtobed/pairtopair subdirectory, so we
// hand-rolled fixtures and compared against the upstream binary by hand. Each
// `Parity_*` test embeds the upstream output we observed for the same inputs.

func TestParity_Either(t *testing.T) {
	a := loadFixture(t, "a.bedpe")
	b := loadFixture(t, "b.bed")
	var out bytes.Buffer
	if _, err := Run(strings.NewReader(a), strings.NewReader(b), &out, Options{Type: TypeEither}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	want := strings.Join([]string{
		"chr1\t10\t20\tchr1\t30\t40\tpair1\t0\t+\t+\tchr1\t5\t25\thit_a\t0\t+",
		"chr1\t10\t20\tchr1\t30\t40\tpair1\t0\t+\t+\tchr1\t35\t45\thit_c\t0\t-",
		"chr1\t100\t200\tchr2\t300\t400\tpair2\t0\t+\t-\tchr1\t95\t195\thit_b\t0\t+",
		"chr1\t100\t200\tchr2\t300\t400\tpair2\t0\t+\t-\tchr2\t350\t450\thit_d\t0\t-",
		"chr2\t1000\t2000\tchr3\t3000\t4000\tpair3\t0\t-\t+\tchr3\t3500\t4500\thit_e\t0\t+",
		"",
	}, "\n")
	if got := out.String(); got != want {
		t.Errorf("either:\n got=%q\nwant=%q", got, want)
	}
}

func TestParity_Both(t *testing.T) {
	a := loadFixture(t, "a.bedpe")
	b := loadFixture(t, "b.bed")
	var out bytes.Buffer
	if _, err := Run(strings.NewReader(a), strings.NewReader(b), &out, Options{Type: TypeBoth}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	// pair1 (both ends hit) + pair2 (both ends hit); pair3 (only end2 hits) suppressed; lonely suppressed.
	got := out.String()
	if !strings.Contains(got, "\tpair1\t") || !strings.Contains(got, "\tpair2\t") {
		t.Errorf("both: missing expected pairs\n%s", got)
	}
	if strings.Contains(got, "\tpair3\t") || strings.Contains(got, "\tlonely\t") {
		t.Errorf("both: should suppress single-end-hitters\n%s", got)
	}
	// 4 output lines: 2 ends each for pair1 & pair2.
	if n := strings.Count(got, "\n"); n != 4 {
		t.Errorf("both: expected 4 output lines, got %d (%s)", n, got)
	}
}

func TestParity_Neither(t *testing.T) {
	a := loadFixture(t, "a.bedpe")
	b := loadFixture(t, "b.bed")
	var out bytes.Buffer
	if _, err := Run(strings.NewReader(a), strings.NewReader(b), &out, Options{Type: TypeNeither}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	// Only "lonely" has zero hits anywhere.
	want := "chr9\t0\t10\tchr9\t100\t200\tlonely\t0\t+\t+\n"
	if got := out.String(); got != want {
		t.Errorf("neither:\n got=%q\nwant=%q", got, want)
	}
}

func TestParity_Xor(t *testing.T) {
	a := loadFixture(t, "a.bedpe")
	b := loadFixture(t, "b.bed")
	var out bytes.Buffer
	if _, err := Run(strings.NewReader(a), strings.NewReader(b), &out, Options{Type: TypeXor}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	// pair3 is the only one with exactly-one end hit.
	want := "chr2\t1000\t2000\tchr3\t3000\t4000\tpair3\t0\t-\t+\tchr3\t3500\t4500\thit_e\t0\t+\n"
	if got := out.String(); got != want {
		t.Errorf("xor:\n got=%q\nwant=%q", got, want)
	}
}

func TestParity_Notboth(t *testing.T) {
	a := loadFixture(t, "a.bedpe")
	b := loadFixture(t, "b.bed")
	var out bytes.Buffer
	if _, err := Run(strings.NewReader(a), strings.NewReader(b), &out, Options{Type: TypeNotboth}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	got := out.String()
	// pair1 and pair2 both hit on both ends -> suppressed.
	if strings.Contains(got, "\tpair1\t") || strings.Contains(got, "\tpair2\t") {
		t.Errorf("notboth: should suppress both-end hitters\n%s", got)
	}
	// pair3: only end2 hits -> single line with end2 hit
	if !strings.Contains(got, "\tpair3\t") {
		t.Errorf("notboth: expected pair3 emission\n%s", got)
	}
	// lonely: no hits -> bare line
	if !strings.Contains(got, "\tlonely\t") {
		t.Errorf("notboth: expected lonely bare emission\n%s", got)
	}
}

// TestParity_Slop is a placeholder: upstream `pairtobed` does not currently
// expose -slop (only `pairtopair` does). When/if it ships we'll add a parity
// check here. Left as a documented skip so reviewers see the gap.
func TestParity_Slop(t *testing.T) {
	t.Skip("pairtobed does not accept -slop upstream; documented intentional gap.")
}
