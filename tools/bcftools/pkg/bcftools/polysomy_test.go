package bcftools

import (
	"bytes"
	"math"
	"strings"
	"testing"
)

// polysomyFixture: one sample, two chromosomes.
//
//   - chr1: 9 het sites with BAFs evenly around 0.5 -> CN2.
//   - chr2: 9 het sites with BAFs clustered around 0.33 -> CN3
//     (trisomy with RA bias).
//   - chr3: only homozygous sites -> n_het=0 -> CN1.
//
// We provide both FORMAT/BAF (preferred) and FORMAT/AD (fallback)
// so the same fixture exercises both code paths via mode toggles in
// individual tests.
func polysomyFixture() string {
	return `##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr1,length=1000>
##contig=<ID=chr2,length=1000>
##contig=<ID=chr3,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
##FORMAT=<ID=BAF,Number=1,Type=Float,Description="B-allele frequency">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	10	.	A	T	.	PASS	.	GT:BAF:AD	0/1:0.50:50,50
chr1	20	.	A	T	.	PASS	.	GT:BAF:AD	0/1:0.49:51,49
chr1	30	.	A	T	.	PASS	.	GT:BAF:AD	0/1:0.51:49,51
chr1	40	.	A	T	.	PASS	.	GT:BAF:AD	0/1:0.48:52,48
chr1	50	.	A	T	.	PASS	.	GT:BAF:AD	0/1:0.52:48,52
chr1	60	.	A	T	.	PASS	.	GT:BAF:AD	0/1:0.47:53,47
chr1	70	.	A	T	.	PASS	.	GT:BAF:AD	0/1:0.53:47,53
chr1	80	.	A	T	.	PASS	.	GT:BAF:AD	0/1:0.46:54,46
chr1	90	.	A	T	.	PASS	.	GT:BAF:AD	0/1:0.54:46,54
chr2	10	.	A	T	.	PASS	.	GT:BAF:AD	0/1:0.33:67,33
chr2	20	.	A	T	.	PASS	.	GT:BAF:AD	0/1:0.32:68,32
chr2	30	.	A	T	.	PASS	.	GT:BAF:AD	0/1:0.34:66,34
chr2	40	.	A	T	.	PASS	.	GT:BAF:AD	0/1:0.30:70,30
chr2	50	.	A	T	.	PASS	.	GT:BAF:AD	0/1:0.35:65,35
chr2	60	.	A	T	.	PASS	.	GT:BAF:AD	0/1:0.33:67,33
chr2	70	.	A	T	.	PASS	.	GT:BAF:AD	0/1:0.31:69,31
chr2	80	.	A	T	.	PASS	.	GT:BAF:AD	0/1:0.36:64,36
chr2	90	.	A	T	.	PASS	.	GT:BAF:AD	0/1:0.30:70,30
chr3	10	.	A	T	.	PASS	.	GT:BAF:AD	0/0:0.00:100,0
chr3	20	.	A	T	.	PASS	.	GT:BAF:AD	1/1:1.00:0,100
chr3	30	.	A	T	.	PASS	.	GT:BAF:AD	0/0:0.00:100,0
`
}

func TestPolysomyMean(t *testing.T) {
	cases := []struct {
		in   []float64
		want float64
	}{
		{nil, 0},
		{[]float64{0.5}, 0.5},
		{[]float64{0.1, 0.3, 0.5}, 0.3},
	}
	for _, c := range cases {
		got := polysomyMean(c.in)
		if math.Abs(got-c.want) > 1e-9 {
			t.Errorf("polysomyMean(%v) = %f want %f", c.in, got, c.want)
		}
	}
}

func TestPolysomyMedian(t *testing.T) {
	cases := []struct {
		in   []float64
		want float64
	}{
		{nil, 0},
		{[]float64{0.5}, 0.5},
		{[]float64{0.3, 0.1, 0.5}, 0.3},      // odd, middle after sort
		{[]float64{0.1, 0.3, 0.5, 0.7}, 0.4}, // even, mean of middles
	}
	for _, c := range cases {
		got := polysomyMedian(c.in)
		if math.Abs(got-c.want) > 1e-9 {
			t.Errorf("polysomyMedian(%v) = %f want %f", c.in, got, c.want)
		}
	}
}

func TestPolysomyIsHet(t *testing.T) {
	cases := []struct {
		in   string
		want bool
	}{
		{"0/1", true},
		{"1|0", true},
		{"0/0", false},
		{"1/1", false},
		{"./.", false},
		{"0/.", false},
		{"", false},
		{"1", false},
	}
	for _, c := range cases {
		if got := polysomyIsHet(c.in); got != c.want {
			t.Errorf("polysomyIsHet(%q) = %v want %v", c.in, got, c.want)
		}
	}
}

func TestPolysomyCNCall(t *testing.T) {
	opts := PolysomyOptions{CnPenalty: 0.7, MinBafDev: 0.1}
	cases := []struct {
		name string
		bafs []float64
		want float64
	}{
		{"no-hets", nil, 1.0},
		{"diploid-perfect", []float64{0.5, 0.5, 0.5}, 2.0},
		{"diploid-noisy", []float64{0.45, 0.5, 0.55}, 2.0},
		{"trisomy-low", []float64{0.33, 0.33, 0.34}, 3.0},
		{"trisomy-high", []float64{0.66, 0.67, 0.67}, 3.0},
		{"borderline", []float64{0.4, 0.4, 0.4}, 2.0}, // |0.4-0.5|=0.1 == threshold
		{"just-over", []float64{0.39, 0.39, 0.39}, 3.0},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			got := polysomyCNCall(c.bafs, opts)
			if got != c.want {
				t.Errorf("polysomyCNCall(%v) = %f want %f", c.bafs, got, c.want)
			}
		})
	}
}

func TestPolysomyCNCall_ForceCN(t *testing.T) {
	opts := PolysomyOptions{CnPenalty: 0.7, MinBafDev: 0.1, ForceCN: 4}
	if got := polysomyCNCall([]float64{0.5, 0.5, 0.5}, opts); got != 4.0 {
		t.Errorf("polysomyCNCall with ForceCN=4 returned %f, want 4.0", got)
	}
}

func TestPolysomyCNCall_PenaltyScaling(t *testing.T) {
	// With CnPenalty < 0.7 the effective threshold shrinks, so
	// borderline cases that were CN2 should flip to CN3.
	bafs := []float64{0.45, 0.45, 0.45} // dev=0.05
	relaxed := PolysomyOptions{CnPenalty: 0.7, MinBafDev: 0.1}
	if got := polysomyCNCall(bafs, relaxed); got != 2.0 {
		t.Errorf("relaxed: want CN2, got %f", got)
	}
	strict := PolysomyOptions{CnPenalty: 0.3, MinBafDev: 0.1} // threshold scales to ~0.043
	if got := polysomyCNCall(bafs, strict); got != 3.0 {
		t.Errorf("strict: want CN3, got %f", got)
	}
}

func TestPolysomy_BAFFromFormat(t *testing.T) {
	var out bytes.Buffer
	sum, err := Polysomy(strings.NewReader(polysomyFixture()), &out, PolysomyOptions{
		Sample: "S1",
	})
	if err != nil {
		t.Fatalf("Polysomy: %v", err)
	}
	if len(sum.Results) != 3 {
		t.Fatalf("want 3 result rows, got %d", len(sum.Results))
	}
	if sum.Results[0].Chrom != "chr1" || sum.Results[0].CN != 2.0 {
		t.Errorf("chr1: want CN=2.0, got chrom=%s CN=%v", sum.Results[0].Chrom, sum.Results[0].CN)
	}
	if sum.Results[1].Chrom != "chr2" || sum.Results[1].CN != 3.0 {
		t.Errorf("chr2: want CN=3.0, got chrom=%s CN=%v", sum.Results[1].Chrom, sum.Results[1].CN)
	}
	if sum.Results[2].Chrom != "chr3" || sum.Results[2].NHet != 0 || sum.Results[2].CN != 1.0 {
		t.Errorf("chr3: want CN=1.0 nhet=0, got chrom=%s nhet=%d CN=%v",
			sum.Results[2].Chrom, sum.Results[2].NHet, sum.Results[2].CN)
	}
	body := out.String()
	for _, want := range []string{
		"# sample\tchrom\tn_het\tmean_baf\tmedian_baf\tcn_call",
		"S1\tchr1\t9\t",
		"\t2.00",
		"S1\tchr2\t9\t",
		"\t3.00",
		"S1\tchr3\t0\t",
		"\t1.00",
	} {
		if !strings.Contains(body, want) {
			t.Errorf("output missing %q:\n%s", want, body)
		}
	}
}

func TestPolysomy_RegionsFilter(t *testing.T) {
	var out bytes.Buffer
	sum, err := Polysomy(strings.NewReader(polysomyFixture()), &out, PolysomyOptions{
		Sample:  "S1",
		Regions: []string{"chr2"},
	})
	if err != nil {
		t.Fatalf("Polysomy: %v", err)
	}
	if len(sum.Results) != 1 {
		t.Errorf("want 1 row (chr2 only), got %d", len(sum.Results))
	}
	if sum.Results[0].Chrom != "chr2" {
		t.Errorf("want chr2, got %s", sum.Results[0].Chrom)
	}
}

func TestPolysomy_UnknownSample(t *testing.T) {
	var out bytes.Buffer
	_, err := Polysomy(strings.NewReader(polysomyFixture()), &out, PolysomyOptions{Sample: "GHOST"})
	if err == nil {
		t.Error("expected error when sample is missing")
	}
}

func TestPolysomy_BAFFromAD(t *testing.T) {
	// Drop the FORMAT/BAF field so the AD fallback kicks in.
	src := strings.ReplaceAll(polysomyFixture(), ":BAF", "")
	// Remove BAF FORMAT meta header too, otherwise the synthetic VCF
	// would still parse it.
	src = strings.ReplaceAll(src,
		`##FORMAT=<ID=BAF,Number=1,Type=Float,Description="B-allele frequency">`+"\n",
		"")
	// In the GT:BAF:AD column we just deleted the BAF token from the
	// FORMAT header, but every per-sample value still has 3 fields.
	// Strip the middle field too.
	lines := strings.Split(src, "\n")
	for i, l := range lines {
		if strings.HasPrefix(l, "chr") {
			cols := strings.Split(l, "\t")
			// The sample column is the last one. Format is GT:AD now;
			// per-sample we need to drop the BAF middle slot.
			if len(cols) >= 10 {
				vals := strings.Split(cols[9], ":")
				if len(vals) == 3 {
					cols[9] = vals[0] + ":" + vals[2]
				}
			}
			lines[i] = strings.Join(cols, "\t")
		}
	}
	src = strings.Join(lines, "\n")

	var out bytes.Buffer
	sum, err := Polysomy(strings.NewReader(src), &out, PolysomyOptions{Sample: "S1"})
	if err != nil {
		t.Fatalf("Polysomy: %v", err)
	}
	// chr1 should still come out as CN2 (BAFs derived from AD are
	// the same 50:50 ratios as the explicit BAFs).
	var chr1, chr2 *PolysomyResult
	for i := range sum.Results {
		switch sum.Results[i].Chrom {
		case "chr1":
			chr1 = &sum.Results[i]
		case "chr2":
			chr2 = &sum.Results[i]
		}
	}
	if chr1 == nil || chr1.CN != 2.0 {
		t.Errorf("chr1 from AD fallback: want CN=2.0, got %+v", chr1)
	}
	if chr2 == nil || chr2.CN != 3.0 {
		t.Errorf("chr2 from AD fallback: want CN=3.0, got %+v", chr2)
	}
}
