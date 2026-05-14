package bedjaccard

import (
	"bytes"
	"strings"
	"testing"
)

func runOf(t *testing.T, a, b string, opts Options) (*Result, string) {
	t.Helper()
	var buf bytes.Buffer
	res, err := Run(strings.NewReader(a), strings.NewReader(b), &buf, opts)
	if err != nil {
		t.Fatalf("Run: %v", err)
	}
	return res, buf.String()
}

func TestRunBasic(t *testing.T) {
	// A: [0,10) + [20,30) = 20 bases
	// B: [5,25)           = 20 bases
	// Intersect: [5,10) + [20,25) = 5 + 5 = 10 bases. Pairs = 2.
	// Union = 20 + 20 - 10 = 30. Jaccard = 10/30 = 0.333...
	a := "chr1\t0\t10\nchr1\t20\t30\n"
	b := "chr1\t5\t25\n"
	res, out := runOf(t, a, b, Options{})
	if res.Intersection != 10 {
		t.Errorf("intersect=%d want 10", res.Intersection)
	}
	if res.Union != 30 {
		t.Errorf("union=%d want 30", res.Union)
	}
	if res.N != 2 {
		t.Errorf("n=%d want 2", res.N)
	}
	if res.Jaccard < 0.3333 || res.Jaccard > 0.3334 {
		t.Errorf("jaccard=%v want ~0.3333", res.Jaccard)
	}
	wantHeader := "intersection\tunion\tjaccard\tn_intersections\n"
	wantBody := "10\t30\t0.333333\t2\n"
	if out != wantHeader+wantBody {
		t.Errorf("output mismatch.\nwant:\n%s\ngot:\n%s", wantHeader+wantBody, out)
	}
}

func TestNoOverlap(t *testing.T) {
	a := "chr1\t0\t10\n"
	b := "chr1\t20\t30\n"
	res, _ := runOf(t, a, b, Options{})
	if res.Intersection != 0 || res.Union != 20 || res.Jaccard != 0 || res.N != 0 {
		t.Errorf("got %+v", res)
	}
}

func TestEmptyBoth(t *testing.T) {
	res, _ := runOf(t, "", "", Options{})
	if res.Intersection != 0 || res.Union != 0 || res.Jaccard != 0 || res.N != 0 {
		t.Errorf("got %+v", res)
	}
}

func TestEmptyA(t *testing.T) {
	res, _ := runOf(t, "", "chr1\t0\t10\n", Options{})
	if res.Intersection != 0 || res.Union != 10 || res.Jaccard != 0 || res.N != 0 {
		t.Errorf("got %+v", res)
	}
}

func TestEmptyB(t *testing.T) {
	res, _ := runOf(t, "chr1\t0\t10\n", "", Options{})
	if res.Intersection != 0 || res.Union != 10 || res.Jaccard != 0 || res.N != 0 {
		t.Errorf("got %+v", res)
	}
}

func TestMultipleChromosomes(t *testing.T) {
	a := "chr1\t0\t10\nchr2\t0\t10\n"
	b := "chr1\t5\t15\nchr2\t0\t5\n"
	// Intersect: chr1 [5,10)=5; chr2 [0,5)=5. Total=10. Pairs=2.
	// |A|=20, |B|=15. Union=20+15-10=25.
	res, _ := runOf(t, a, b, Options{})
	if res.Intersection != 10 || res.Union != 25 || res.N != 2 {
		t.Errorf("got %+v", res)
	}
}

func TestBOnlyOnLaterChroms(t *testing.T) {
	// A on chr1 only, B on chr2 only.
	a := "chr1\t0\t10\n"
	b := "chr2\t0\t10\n"
	res, _ := runOf(t, a, b, Options{})
	if res.Intersection != 0 || res.Union != 20 || res.N != 0 {
		t.Errorf("got %+v", res)
	}
}

func TestStrandSame(t *testing.T) {
	// A is +, B has a + and a - record covering the same region. Only + counts.
	a := "chr1\t0\t10\ta\t0\t+\n"
	b := "chr1\t0\t10\tx\t0\t+\nchr1\t0\t10\ty\t0\t-\n"
	res, _ := runOf(t, a, b, Options{SameStrand: true})
	// Pairs: one (only same-strand). intersect: 10. union: 10+20-10=20.
	if res.N != 1 || res.Intersection != 10 || res.Union != 20 {
		t.Errorf("got %+v", res)
	}
}

func TestStrandOpposite(t *testing.T) {
	a := "chr1\t0\t10\ta\t0\t+\n"
	b := "chr1\t0\t10\tx\t0\t+\nchr1\t0\t10\ty\t0\t-\n"
	res, _ := runOf(t, a, b, Options{OppositeStrand: true})
	if res.N != 1 || res.Intersection != 10 || res.Union != 20 {
		t.Errorf("got %+v", res)
	}
}

func TestStrandWithoutStrandField(t *testing.T) {
	// BED3: no strand. -s should match nothing.
	a := "chr1\t0\t10\n"
	b := "chr1\t0\t10\n"
	res, _ := runOf(t, a, b, Options{SameStrand: true})
	if res.N != 0 {
		t.Errorf("expected 0 pairs without strand, got %+v", res)
	}
}

func TestFractionA(t *testing.T) {
	// A=[0,100). B=[0,10). Overlap=10, fracA=0.1.
	a := "chr1\t0\t100\n"
	b := "chr1\t0\t10\n"
	// With -f 0.5, the pair should be dropped.
	res, _ := runOf(t, a, b, Options{FractionA: 0.5})
	if res.N != 0 || res.Intersection != 0 {
		t.Errorf("got %+v want N=0", res)
	}
	res, _ = runOf(t, a, b, Options{FractionA: 0.05})
	if res.N != 1 || res.Intersection != 10 {
		t.Errorf("got %+v", res)
	}
}

func TestFractionB(t *testing.T) {
	a := "chr1\t0\t10\n"
	b := "chr1\t0\t100\n"
	res, _ := runOf(t, a, b, Options{FractionB: 0.5})
	if res.N != 0 || res.Intersection != 0 {
		t.Errorf("got %+v", res)
	}
}

func TestFractionOutOfRange(t *testing.T) {
	var buf bytes.Buffer
	if _, err := Run(strings.NewReader(""), strings.NewReader(""), &buf,
		Options{FractionA: 2}); err == nil {
		t.Fatal("expected error for FractionA > 1")
	}
	if _, err := Run(strings.NewReader(""), strings.NewReader(""), &buf,
		Options{FractionB: -1}); err == nil {
		t.Fatal("expected error for FractionB < 0")
	}
}

func TestStrandFlagsMutuallyExclusive(t *testing.T) {
	var buf bytes.Buffer
	if _, err := Run(strings.NewReader(""), strings.NewReader(""), &buf,
		Options{SameStrand: true, OppositeStrand: true}); err == nil {
		t.Fatal("expected error for -s and -S together")
	}
}

func TestUnsortedAErrors(t *testing.T) {
	a := "chr1\t100\t200\nchr1\t0\t50\n"
	b := "chr1\t0\t300\n"
	var buf bytes.Buffer
	if _, err := Run(strings.NewReader(a), strings.NewReader(b), &buf, Options{}); err == nil {
		t.Fatal("expected unsorted A error")
	}
}

func TestUnsortedBErrors(t *testing.T) {
	a := "chr1\t0\t300\n"
	b := "chr1\t100\t200\nchr1\t0\t50\n"
	var buf bytes.Buffer
	if _, err := Run(strings.NewReader(a), strings.NewReader(b), &buf, Options{}); err == nil {
		t.Fatal("expected unsorted B error")
	}
}

func TestBadInputErrors(t *testing.T) {
	var buf bytes.Buffer
	if _, err := Run(strings.NewReader("chr1\tnope\t10\n"), strings.NewReader(""), &buf, Options{}); err == nil {
		t.Fatal("expected parse error on bad A")
	}
	if _, err := Run(strings.NewReader(""), strings.NewReader("chr1\t0\tnope\n"), &buf, Options{}); err == nil {
		t.Fatal("expected parse error on bad B")
	}
	// Bad B encountered mid-sweep too (after A consumed).
	if _, err := Run(strings.NewReader("chr1\t0\t10\n"), strings.NewReader("chr1\tabc\t10\n"), &buf, Options{}); err == nil {
		t.Fatal("expected parse error on bad B mid-sweep")
	}
}

func TestOverlappingPairsCounted(t *testing.T) {
	// A=[0,10). B=[0,5),[5,10). Overlap = 10 across two distinct pairs.
	a := "chr1\t0\t10\n"
	b := "chr1\t0\t5\nchr1\t5\t10\n"
	res, _ := runOf(t, a, b, Options{})
	if res.N != 2 || res.Intersection != 10 {
		t.Errorf("got %+v want N=2,intersect=10", res)
	}
}

func TestStreamingDoesNotKeepAllB(t *testing.T) {
	// Verify many distant B records don't all stay active (we can't observe
	// memory directly in a unit test, but exercising the path covers prune
	// branches).
	var a, b strings.Builder
	for i := 0; i < 100; i++ {
		fmtCoord(&a, "chr1", i*1000, i*1000+10)
	}
	for i := 0; i < 100; i++ {
		fmtCoord(&b, "chr1", i*1000+5, i*1000+15)
	}
	res, _ := runOf(t, a.String(), b.String(), Options{})
	// Each pair overlaps 5 bases.
	if res.N != 100 || res.Intersection != 500 {
		t.Errorf("got %+v", res)
	}
}

func TestFormatJaccard(t *testing.T) {
	if s := formatJaccard(0); s != "0" {
		t.Errorf("got %q", s)
	}
	if s := formatJaccard(1.0 / 3.0); s != "0.333333" {
		t.Errorf("got %q", s)
	}
}

// fmtCoord is a tiny helper for table-building in this test file.
func fmtCoord(sb *strings.Builder, chrom string, start, end int) {
	sb.WriteString(chrom)
	sb.WriteByte('\t')
	sb.WriteString(itoa(start))
	sb.WriteByte('\t')
	sb.WriteString(itoa(end))
	sb.WriteByte('\n')
}

func itoa(n int) string {
	// Keeps the test file independent of strconv.
	if n == 0 {
		return "0"
	}
	neg := n < 0
	if neg {
		n = -n
	}
	var b [20]byte
	i := len(b)
	for n > 0 {
		i--
		b[i] = byte('0' + n%10)
		n /= 10
	}
	if neg {
		i--
		b[i] = '-'
	}
	return string(b[i:])
}
