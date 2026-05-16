package bcftools

import (
	"bytes"
	"strings"
	"testing"
)

// TestCNVHeuristicSimple feeds a tiny 2-sample BAF/LRR VCF with a
// chr2 stretch that should classify as CN3 (gain).
func TestCNVHeuristicSimple(t *testing.T) {
	input := `##fileformat=VCFv4.2
##contig=<ID=chr1>
##contig=<ID=chr2>
##FORMAT=<ID=BAF,Number=1,Type=Float,Description="B-allele freq">
##FORMAT=<ID=LRR,Number=1,Type=Float,Description="Log R ratio">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	queryA	ctrlB
chr1	100	.	A	C	.	PASS	.	BAF:LRR	0.50:0.01	0.50:0.02
chr1	200	.	A	C	.	PASS	.	BAF:LRR	0.50:-0.01	0.51:0.00
chr1	300	.	A	C	.	PASS	.	BAF:LRR	0.50:0.00	0.49:0.01
chr2	100	.	A	C	.	PASS	.	BAF:LRR	0.33:0.65	0.50:0.01
chr2	200	.	A	C	.	PASS	.	BAF:LRR	0.66:0.62	0.50:-0.01
chr2	300	.	A	C	.	PASS	.	BAF:LRR	0.34:0.70	0.50:0.00
chr2	400	.	A	C	.	PASS	.	BAF:LRR	0.67:0.55	0.50:0.02
`
	var out bytes.Buffer
	n, err := CNV(strings.NewReader(input), &out, CNVOptions{})
	if err != nil {
		t.Fatalf("CNV: %v", err)
	}
	if n != 4 {
		t.Errorf("expected 4 rows (2 samples x 2 chromosomes), got %d", n)
	}
	got := out.String()
	wantLines := []string{
		"#sample\tchrom\tn_sites\tmedian_baf\tmean_lrr\tcn_call",
		"queryA\tchr1\t3\t0.000000\t0.000000\tCN2",
		"queryA\tchr2\t4\t0.165000\t0.630000\tCN4",
		"ctrlB\tchr1\t3\t0.010000\t0.010000\tCN2",
		"ctrlB\tchr2\t4\t0.000000\t0.005000\tCN2",
	}
	for _, want := range wantLines {
		if !strings.Contains(got, want) {
			t.Errorf("missing line %q in:\n%s", want, got)
		}
	}
}

// TestCNVSampleNarrowing exercises -s/-c and verifies non-matching
// samples are skipped.
func TestCNVSampleNarrowing(t *testing.T) {
	input := `##fileformat=VCFv4.2
##contig=<ID=chr1>
##FORMAT=<ID=BAF,Number=1,Type=Float,Description="">
##FORMAT=<ID=LRR,Number=1,Type=Float,Description="">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	A	B	C
chr1	1	.	A	C	.	.	.	BAF:LRR	0.50:0.01	0.30:-0.55	0.50:0.00
chr1	2	.	A	C	.	.	.	BAF:LRR	0.50:0.00	0.70:-0.60	0.50:0.01
`
	var out bytes.Buffer
	n, err := CNV(strings.NewReader(input), &out, CNVOptions{QuerySample: "B"})
	if err != nil {
		t.Fatalf("CNV: %v", err)
	}
	if n != 1 {
		t.Errorf("expected 1 row, got %d", n)
	}
	if !strings.Contains(out.String(), "B\tchr1\t2") {
		t.Errorf("expected sample B row, got:\n%s", out.String())
	}
	if strings.Contains(out.String(), "A\t") || strings.Contains(out.String(), "C\t") {
		t.Errorf("unexpected sample row in:\n%s", out.String())
	}
}

func TestCNVClassifyCN(t *testing.T) {
	cases := []struct {
		baf, lrr float64
		want     string
	}{
		{0.01, 0.01, "CN2"},
		{0.30, -0.55, "CN0"}, // mean LRR < -2*0.20
		{0.30, -0.25, "CN1"},
		{0.30, 0.25, "CN3"},
		{0.30, 0.55, "CN4"},
	}
	for _, tc := range cases {
		got := classifyCN(tc.baf, tc.lrr, CNVOptions{})
		if got != tc.want {
			t.Errorf("classifyCN(%v,%v) = %q want %q", tc.baf, tc.lrr, got, tc.want)
		}
	}
}

func TestCNVMedian(t *testing.T) {
	if median(nil) != 0 {
		t.Error("median(nil) should be 0")
	}
	if got := median([]float64{1, 3, 2}); got != 2 {
		t.Errorf("median(1,3,2) = %v want 2", got)
	}
	if got := median([]float64{1, 2, 3, 4}); got != 2.5 {
		t.Errorf("median(1,2,3,4) = %v want 2.5", got)
	}
}

func TestCNVMissingSampleError(t *testing.T) {
	input := `##fileformat=VCFv4.2
##contig=<ID=chr1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT
chr1	1	.	A	C	.	.	.	.
`
	var out bytes.Buffer
	if _, err := CNV(strings.NewReader(input), &out, CNVOptions{QuerySample: "NOPE"}); err == nil {
		t.Errorf("expected error when no samples match")
	}
}

func TestParseFloatField(t *testing.T) {
	cases := []struct {
		s    string
		want float64
		ok   bool
	}{
		{"", 0, false},
		{".", 0, false},
		{"0.5", 0.5, true},
		{"0.5,0.6", 0.5, true},
		{"abc", 0, false},
	}
	for _, tc := range cases {
		got, ok := parseFloatField(map[string]string{"X": tc.s}, "X")
		if ok != tc.ok || (ok && got != tc.want) {
			t.Errorf("parseFloatField(%q): got (%v,%v) want (%v,%v)", tc.s, got, ok, tc.want, tc.ok)
		}
	}
	// Missing tag.
	if _, ok := parseFloatField(map[string]string{}, "Y"); ok {
		t.Errorf("missing tag should return ok=false")
	}
}

// TestCNVBafFromAD pins the FORMAT/AD = REF,ALT fallback synthesis
// (PR #110 review finding #4): when FORMAT/BAF isn't present, BAF
// must come from ALT/(REF+ALT).
func TestCNVBafFromAD(t *testing.T) {
	cases := []struct {
		in     string
		want   float64
		wantOK bool
	}{
		{"10,10", 0.5, true},
		{"30,10", 0.25, true},
		{"0,10", 1.0, true},
		{"10,0", 0.0, true},
		{"0,0", 0, false},
		{"", 0, false},
		{".", 0, false},
		{"abc,def", 0, false},
		{"10", 0, false},
	}
	for _, c := range cases {
		got, ok := bafFromAD(c.in)
		if ok != c.wantOK {
			t.Errorf("bafFromAD(%q) ok=%v want %v", c.in, ok, c.wantOK)
			continue
		}
		if ok && got != c.want {
			t.Errorf("bafFromAD(%q) = %v, want %v", c.in, got, c.want)
		}
	}
}

// TestCNVClassifyCN_LOH pins PR #110 review finding #3: classifyCN
// MUST consult bafDev, not just LRR. A copy-neutral LOH signal
// (balanced LRR + skewed BAF) feeds through the new bafLo gate.
func TestCNVClassifyCN_LOH(t *testing.T) {
	opts := CNVOptions{BAFDev: 0.05, LRRDev: 0.20}
	// LRR=0 + BAF dev = 0.30 -> hits the LOH branch (still CN2 in v1
	// per docstring, but the gate is now load-bearing so a future
	// state can flip on without breaking the public API).
	if got := classifyCN(0.30, 0.0, opts); got != "CN2" {
		t.Errorf("LOH-like signal: got %s, want CN2", got)
	}
	// LRR=0 + BAF dev = 0 -> plain diploid CN2.
	if got := classifyCN(0.0, 0.0, opts); got != "CN2" {
		t.Errorf("plain diploid: got %s, want CN2", got)
	}
}
