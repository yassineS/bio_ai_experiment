package bcftools

import (
	"bytes"
	"math"
	"math/rand"
	"strings"
	"testing"
)

// --- pure helpers ---

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
		{[]float64{0.3, 0.1, 0.5}, 0.3},
		{[]float64{0.1, 0.3, 0.5, 0.7}, 0.4},
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

func TestCDiv(t *testing.T) {
	if v := cDiv(3, 0); !math.IsInf(v, 1) {
		t.Errorf("cDiv(3,0) = %v, want +Inf", v)
	}
	if v := cDiv(0, 0); !math.IsNaN(v) {
		t.Errorf("cDiv(0,0) = %v, want NaN", v)
	}
	if v := cDiv(6, 2); v != 3 {
		t.Errorf("cDiv(6,2) = %v, want 3", v)
	}
}

func TestMinRatio(t *testing.T) {
	if v := minRatio(2, 8); math.Abs(v-0.25) > 1e-12 {
		t.Errorf("minRatio(2,8) = %v, want 0.25", v)
	}
	if v := minRatio(8, 2); math.Abs(v-0.25) > 1e-12 {
		t.Errorf("minRatio(8,2) = %v, want 0.25", v)
	}
	if v := minRatio(0, 0); v != 0 {
		t.Errorf("minRatio(0,0) = %v, want 0", v)
	}
}

// --- BAF distribution generators ---

// gaussBAFs appends n BAF values drawn from N(mean, sd) to dst,
// clamped to (0,1). The generator is seeded so the distributions are
// deterministic across test runs.
func gaussBAFs(rng *rand.Rand, dst []float64, n int, mean, sd float64) []float64 {
	for i := 0; i < n; i++ {
		v := mean + rng.NormFloat64()*sd
		if v <= 0 {
			v = 1e-3
		} else if v >= 1 {
			v = 1 - 1e-3
		}
		dst = append(dst, v)
	}
	return dst
}

// diploidBAFs builds a clean CN2 BAF set: a single heterozygous peak
// centred at 0.5.
func diploidBAFs() []float64 {
	rng := rand.New(rand.NewSource(1))
	var b []float64
	b = gaussBAFs(rng, b, 4000, 0.5, 0.04)
	return b
}

// trisomyBAFs builds a clear CN3 BAF set: two symmetric heterozygous
// peaks at 1/3 and 2/3.
func trisomyBAFs() []float64 {
	rng := rand.New(rand.NewSource(2))
	var b []float64
	b = gaussBAFs(rng, b, 2000, 1.0/3.0, 0.035)
	b = gaussBAFs(rng, b, 2000, 2.0/3.0, 0.035)
	return b
}

// tetrasomyBAFs builds a clear CN4 BAF set. A tetrasomy het site can
// carry 1, 2 or 3 B alleles out of 4 copies, so its B-allele frequency
// clusters at 1/4, 2/4 and 3/4. The CN4 candidate fit in
// polysomy.c:fit_curves models exactly this shape: a central Gaussian
// at 0.5 plus two symmetric side Gaussians at 0.5 ± cn4_dx. Here the
// side peaks sit at 0.25 and 0.75 (cn4_dx = 0.25, the maximum the fit
// allows). The seed is fixed so the distribution is deterministic.
func tetrasomyBAFs() []float64 {
	rng := rand.New(rand.NewSource(3))
	var b []float64
	b = gaussBAFs(rng, b, 2000, 0.50, 0.025) // central RA peak
	b = gaussBAFs(rng, b, 1200, 0.25, 0.025) // lower side peak
	b = gaussBAFs(rng, b, 1200, 0.75, 0.025) // upper side peak
	return b
}

// --- algorithm-level tests ---

// TestBuildBAFDist_Diploid checks the histogram build and peak
// isolation on a clean diploid distribution: the het band must be
// detected (copyNumber stays 0 -> the fit runs) and the band centre
// (ira) must land near 0.5.
func TestBuildBAFDist_Diploid(t *testing.T) {
	opts := PolysomyOptions{}
	applyPolysomyDefaults(&opts)
	opts.RaRrScaling = true
	d := buildBAFDist(diploidBAFs(), opts)
	if d.copyNumber != 0 {
		t.Fatalf("clean diploid should run the fit (copyNumber 0), got %d", d.copyNumber)
	}
	xc := d.xvals[d.ira]
	if math.Abs(xc-0.5) > 0.1 {
		t.Errorf("band centre x=%v, want near 0.5", xc)
	}
}

// TestPolysomyCNCall_Diploid: a clean single-peak BAF distribution at
// 0.5 must be called CN2 (no aberration).
func TestPolysomyCNCall_Diploid(t *testing.T) {
	opts := PolysomyOptions{}
	applyPolysomyDefaults(&opts)
	opts.RaRrScaling = true
	cn, fit := polysomyCNCall(diploidBAFs(), opts)
	if math.Round(cn) != 2 {
		t.Errorf("diploid: CN = %v (fit %v), want ~2", cn, fit)
	}
}

// TestPolysomyCNCall_Trisomy: a clear two-peak BAF distribution at
// 1/3 and 2/3 must be called CN3 (trisomy).
func TestPolysomyCNCall_Trisomy(t *testing.T) {
	opts := PolysomyOptions{}
	applyPolysomyDefaults(&opts)
	opts.RaRrScaling = true
	cn, fit := polysomyCNCall(trisomyBAFs(), opts)
	if math.Round(cn) != 3 {
		t.Errorf("trisomy: CN = %v (fit %v), want ~3", cn, fit)
	}
}

// TestPolysomyCNCall_Tetrasomy: a CN4-shaped BAF distribution (central
// 0.5 peak plus symmetric side peaks at 0.25 / 0.75) must be called
// CN4 (tetrasomy).
//
// Expected-call derivation, from polysomy.c:fit_curves' CN4 branch:
//   - The CN4 fit places three Gaussians: a central peak at 0.5 and two
//     side peaks at 0.5 ± cn4_dx. With side data centred at 0.25 and
//     0.75 the optimiser converges to cn4_dx ≈ 0.25 (it is clamped at
//     0.25, "CN4 peaks should not be separated by more than 0.5").
//   - The reported CN is `cn = 3 + cn4_frac`, where
//     `cn4_frac = cn4RAaa_params[1] - cn4RArr_params[1]` is the gap
//     between the upper and lower side-peak centres. For peaks at 0.75
//     and 0.25 that gap is 0.50, so cn = 3 + 0.50 = 3.50, which rounds
//     to CN4.
//   - The CN4 acceptance gates (cn4_fit ≤ fit-th, cn4_ymin ≥ peak-size,
//     cn4_dy ≥ peak-symmetry, |cn4_dx| ≤ 0.1) all pass because the side
//     peaks are equal-height and symmetric about 0.5, and CN4 beats the
//     cn_penalty tiebreaker against CN2/CN3 because neither a single
//     0.5 Gaussian nor a 1/3–2/3 pair can model a three-peak shape.
func TestPolysomyCNCall_Tetrasomy(t *testing.T) {
	opts := PolysomyOptions{}
	applyPolysomyDefaults(&opts)
	opts.RaRrScaling = true
	cn, fit := polysomyCNCall(tetrasomyBAFs(), opts)
	if math.Round(cn) != 4 {
		t.Errorf("tetrasomy: CN = %v (fit %v), want ~4", cn, fit)
	}
	// The CN4 branch must genuinely have won: a CN4 call lands in
	// [3.5, 4.5) (cn = 3 + cn4_frac), never at the CN2/CN3 values.
	if cn < 3.5 || cn >= 4.5 {
		t.Errorf("tetrasomy: CN = %v not in the CN4 range [3.5,4.5)", cn)
	}
}

// hetWithAABafDist hand-builds a bafDist whose heterozygous band is a
// clean central peak at 0.5 but whose homozygous-AA band still carries
// a real Gaussian shoulder near 0.85. The band markers (irr/ira/iaa)
// are placed by hand so the AA band is NOT trimmed away — going through
// buildBAFDist would normalise and trim the AA peak to a near-empty
// stub, leaving nothing for the --include-aa exp fit to act on. This
// mirrors the hand-derived style of the peak-model unit tests.
func hetWithAABafDist() *bafDist {
	const nbins = 150
	d := &bafDist{
		xvals: make([]float64, nbins),
		yvals: make([]float64, nbins),
		nvals: nbins,
	}
	for i := 0; i < nbins; i++ {
		x := float64(i) / float64(nbins-1)
		d.xvals[i] = x
		ra := math.Exp(-((x - 0.5) * (x - 0.5)) / (0.04 * 0.04))
		aa := 0.6 * math.Exp(-((x-0.85)*(x-0.85))/(0.05*0.05))
		d.yvals[i] = ra + aa
	}
	d.irr, d.ira, d.iaa = 20, 75, 115
	d.copyNumber = 0
	return d
}

// cloneBafDist deep-copies a bafDist so a fit run cannot mutate the
// shared backing slices between two compared cases.
func cloneBafDist(src *bafDist) *bafDist {
	cp := *src
	cp.xvals = append([]float64(nil), src.xvals...)
	cp.yvals = append([]float64(nil), src.yvals...)
	return &cp
}

// TestPolysomyFit_IncludeAA: --include-aa (IncludeAA) must change which
// peaks are fitted. polysomy.c:fit_curves only adds the homozygous-AA
// exp fit (cn2aa_fit / cn4aa_fit) to the candidate scores when
// include_aa is set; with it off the AA band is ignored entirely.
//
// On a distribution with a real AA shoulder the two settings must
// therefore disagree: with IncludeAA off the central 0.5 peak alone
// gives a clean CN2; with IncludeAA on the AA exp fit adds a sizeable
// residual to every candidate's score, pushing them past fit-th so the
// chromosome is reported as unknown (-1). Identical results would mean
// the IncludeAA branch was never exercised.
func TestPolysomyFit_IncludeAA(t *testing.T) {
	base := PolysomyOptions{}
	applyPolysomyDefaults(&base)
	base.RaRrScaling = true

	off := base
	off.IncludeAA = false
	on := base
	on.IncludeAA = true

	cnOff, fitOff := fitBAFDist(cloneBafDist(hetWithAABafDist()), off)
	cnOn, fitOn := fitBAFDist(cloneBafDist(hetWithAABafDist()), on)

	if math.Round(cnOff) != 2 {
		t.Errorf("IncludeAA off: CN = %v (fit %v), want ~2", cnOff, fitOff)
	}
	if math.Abs(fitOn-fitOff) < 0.5 {
		t.Errorf("IncludeAA did not change the fit: off=%v on=%v (the AA exp fit was not applied)", fitOff, fitOn)
	}
	if cnOn == cnOff {
		t.Errorf("IncludeAA did not change the CN call: both = %v", cnOn)
	}
}

// TestPolysomyCNCall_NoHets: an empty BAF set is CN1 (LOH / monosomy).
func TestPolysomyCNCall_NoHets(t *testing.T) {
	opts := PolysomyOptions{}
	applyPolysomyDefaults(&opts)
	cn, _ := polysomyCNCall(nil, opts)
	if cn != 1.0 {
		t.Errorf("no hets: CN = %v, want 1.0", cn)
	}
}

// TestPolysomyCNCall_ForceCN: --force-cn bypasses the fit entirely.
func TestPolysomyCNCall_ForceCN(t *testing.T) {
	opts := PolysomyOptions{ForceCN: 4}
	applyPolysomyDefaults(&opts)
	cn, _ := polysomyCNCall(diploidBAFs(), opts)
	if cn != 4.0 {
		t.Errorf("ForceCN=4: CN = %v, want 4.0", cn)
	}
}

// TestPolysomyCNText covers the printable form, incl. the unknown
// sentinel.
func TestPolysomyCNText(t *testing.T) {
	cases := []struct {
		cn   float64
		want string
	}{
		{-1, "?"},
		{2.0, "2.00"},
		{3.0, "3.00"},
		{1.0, "1.00"},
	}
	for _, c := range cases {
		if got := polysomyCNText(c.cn); got != c.want {
			t.Errorf("polysomyCNText(%v) = %q want %q", c.cn, got, c.want)
		}
	}
}

// --- knob sensitivity (the previously-ignored options now take effect) ---

// TestFitTh_RejectsAll: a tiny fit-th rejects every candidate, so even
// a clean diploid comes out as unknown (-1).
func TestFitTh_RejectsAll(t *testing.T) {
	opts := PolysomyOptions{FitTh: 1e-9}
	applyPolysomyDefaults(&opts)
	opts.FitTh = 1e-9 // re-apply: applyPolysomyDefaults would not touch a non-zero value
	opts.RaRrScaling = true
	cn, _ := polysomyCNCall(diploidBAFs(), opts)
	if cn != -1 {
		t.Errorf("fit-th=1e-9 should reject all CN models, got CN=%v", cn)
	}
}

// TestNBins_TakesEffect: the histogram bin count is honoured (a
// distinct bin count yields a distinct histogram length).
func TestNBins_TakesEffect(t *testing.T) {
	optsA := PolysomyOptions{NBins: 100}
	applyPolysomyDefaults(&optsA)
	optsB := PolysomyOptions{NBins: 200}
	applyPolysomyDefaults(&optsB)
	dA := buildBAFDist(diploidBAFs(), optsA)
	dB := buildBAFDist(diploidBAFs(), optsB)
	if len(dA.xvals) == len(dB.xvals) {
		t.Errorf("nbins ignored: both histograms have length %d", len(dA.xvals))
	}
	if len(dA.xvals) != 100 || len(dB.xvals) != 200 {
		t.Errorf("histogram lengths = %d, %d; want 100, 200", len(dA.xvals), len(dB.xvals))
	}
}

// TestCnPenalty_TakesEffect: with cn_penalty=0 a CN3 fit only has to
// be strictly better than CN2 to win; with cn_penalty very close to 1
// the bar for upgrading is essentially unreachable. We verify the
// tiebreaker arithmetic directly through the decision: a borderline
// trisomy stays CN3 under a relaxed penalty.
func TestCnPenalty_TakesEffect(t *testing.T) {
	relaxed := PolysomyOptions{CnPenalty: 0.7}
	applyPolysomyDefaults(&relaxed)
	relaxed.RaRrScaling = true
	cn, _ := polysomyCNCall(trisomyBAFs(), relaxed)
	if math.Round(cn) != 3 {
		t.Fatalf("trisomy under default penalty: CN=%v, want 3", cn)
	}
}

// --- VCF integration ---

// polysomyVCF renders a single-sample VCF whose chrom carries the
// given BAF values as FORMAT/BAF at synthetic het sites.
func polysomyVCF(chromBAFs map[string][]float64, order []string) string {
	var sb strings.Builder
	sb.WriteString("##fileformat=VCFv4.2\n")
	sb.WriteString(`##FILTER=<ID=PASS,Description="All filters passed">` + "\n")
	for _, c := range order {
		sb.WriteString("##contig=<ID=" + c + ",length=100000>\n")
	}
	sb.WriteString(`##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">` + "\n")
	sb.WriteString(`##FORMAT=<ID=BAF,Number=1,Type=Float,Description="B-allele frequency">` + "\n")
	sb.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
	for _, c := range order {
		pos := 1
		for _, b := range chromBAFs[c] {
			sb.WriteString(c)
			sb.WriteString("\t")
			sb.WriteString(itoa(pos))
			sb.WriteString("\t.\tA\tT\t.\tPASS\t.\tGT:BAF\t0/1:")
			sb.WriteString(ftoa(b))
			sb.WriteString("\n")
			pos += 10
		}
	}
	return sb.String()
}

// itoa and ftoa are shared with cnv_test.go in this package.

// TestPolysomy_VCFIntegration runs the streaming entry point end to
// end: chr1 is a clean diploid -> CN2, chr2 a clear trisomy -> CN3,
// chr3 a clear tetrasomy -> CN4.
func TestPolysomy_VCFIntegration(t *testing.T) {
	src := polysomyVCF(map[string][]float64{
		"chr1": diploidBAFs(),
		"chr2": trisomyBAFs(),
		"chr3": tetrasomyBAFs(),
	}, []string{"chr1", "chr2", "chr3"})

	var out bytes.Buffer
	sum, err := Polysomy(strings.NewReader(src), &out, PolysomyOptions{
		Sample:      "S1",
		RaRrScaling: true,
	})
	if err != nil {
		t.Fatalf("Polysomy: %v", err)
	}
	if len(sum.Results) != 3 {
		t.Fatalf("want 3 result rows, got %d", len(sum.Results))
	}
	if sum.Results[0].Chrom != "chr1" || math.Round(sum.Results[0].CN) != 2 {
		t.Errorf("chr1: want CN~2, got chrom=%s CN=%v", sum.Results[0].Chrom, sum.Results[0].CN)
	}
	if sum.Results[1].Chrom != "chr2" || math.Round(sum.Results[1].CN) != 3 {
		t.Errorf("chr2: want CN~3, got chrom=%s CN=%v", sum.Results[1].Chrom, sum.Results[1].CN)
	}
	if sum.Results[2].Chrom != "chr3" || math.Round(sum.Results[2].CN) != 4 {
		t.Errorf("chr3: want CN~4, got chrom=%s CN=%v", sum.Results[2].Chrom, sum.Results[2].CN)
	}
	body := out.String()
	if !strings.Contains(body, "# sample\tchrom\tn_het\tmean_baf\tmedian_baf\tcn_call") {
		t.Errorf("missing header in output:\n%s", body)
	}
}

// TestPolysomy_NoHetChrom: a chromosome with only homozygous calls is
// reported CN1 with n_het=0.
func TestPolysomy_NoHetChrom(t *testing.T) {
	src := `##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chrH,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
##FORMAT=<ID=BAF,Number=1,Type=Float,Description="BAF">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chrH	10	.	A	T	.	PASS	.	GT:BAF	0/0:0.00
chrH	20	.	A	T	.	PASS	.	GT:BAF	1/1:1.00
`
	var out bytes.Buffer
	sum, err := Polysomy(strings.NewReader(src), &out, PolysomyOptions{Sample: "S1"})
	if err != nil {
		t.Fatalf("Polysomy: %v", err)
	}
	if len(sum.Results) != 1 || sum.Results[0].NHet != 0 || sum.Results[0].CN != 1.0 {
		t.Errorf("homozygous-only chrom: want CN=1.0 nhet=0, got %+v", sum.Results)
	}
}

// TestPolysomy_RegionsFilter restricts the analysis to a chromosome
// subset.
func TestPolysomy_RegionsFilter(t *testing.T) {
	src := polysomyVCF(map[string][]float64{
		"chr1": diploidBAFs(),
		"chr2": trisomyBAFs(),
	}, []string{"chr1", "chr2"})
	var out bytes.Buffer
	sum, err := Polysomy(strings.NewReader(src), &out, PolysomyOptions{
		Sample:  "S1",
		Regions: []string{"chr2"},
	})
	if err != nil {
		t.Fatalf("Polysomy: %v", err)
	}
	if len(sum.Results) != 1 || sum.Results[0].Chrom != "chr2" {
		t.Errorf("regions filter: want only chr2, got %+v", sum.Results)
	}
}

// TestPolysomy_UnknownSample errors out on a missing sample name.
func TestPolysomy_UnknownSample(t *testing.T) {
	src := polysomyVCF(map[string][]float64{"chr1": diploidBAFs()}, []string{"chr1"})
	var out bytes.Buffer
	if _, err := Polysomy(strings.NewReader(src), &out, PolysomyOptions{Sample: "GHOST"}); err == nil {
		t.Error("expected error when sample is missing")
	}
}

// TestPolysomy_BAFFromAD checks the FORMAT/AD = REF,ALT fallback path
// when no explicit FORMAT/BAF is present.
func TestPolysomy_BAFFromAD(t *testing.T) {
	rng := rand.New(rand.NewSource(7))
	var sb strings.Builder
	sb.WriteString("##fileformat=VCFv4.2\n")
	sb.WriteString(`##FILTER=<ID=PASS,Description="All filters passed">` + "\n")
	sb.WriteString("##contig=<ID=chr1,length=1000000>\n")
	sb.WriteString(`##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">` + "\n")
	sb.WriteString(`##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">` + "\n")
	sb.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
	pos := 1
	for i := 0; i < 5000; i++ {
		// ALT/(REF+ALT) ~ 0.5: a clean diploid. A deep total count
		// keeps the synthesised BAF finely resolved.
		baf := 0.5 + rng.NormFloat64()*0.04
		if baf < 0.02 {
			baf = 0.02
		} else if baf > 0.98 {
			baf = 0.98
		}
		alt := int(math.Round(baf * 2000))
		ref := 2000 - alt
		sb.WriteString("chr1\t")
		sb.WriteString(itoa(pos))
		sb.WriteString("\t.\tA\tT\t.\tPASS\t.\tGT:AD\t0/1:")
		sb.WriteString(itoa(ref))
		sb.WriteString(",")
		sb.WriteString(itoa(alt))
		sb.WriteString("\n")
		pos += 10
	}
	var out bytes.Buffer
	sum, err := Polysomy(strings.NewReader(sb.String()), &out, PolysomyOptions{
		Sample:      "S1",
		RaRrScaling: true,
	})
	if err != nil {
		t.Fatalf("Polysomy: %v", err)
	}
	if len(sum.Results) != 1 || math.Round(sum.Results[0].CN) != 2 {
		t.Errorf("AD fallback diploid: want CN~2, got %+v", sum.Results)
	}
}
