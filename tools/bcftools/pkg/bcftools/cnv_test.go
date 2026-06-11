package bcftools

import (
	"bytes"
	"math"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// vcfHeaderSingle is a minimal single-sample BAF/LRR VCF header.
const vcfHeaderSingle = `##fileformat=VCFv4.2
##contig=<ID=chr1>
##contig=<ID=chr2>
##FORMAT=<ID=BAF,Number=1,Type=Float,Description="B-allele freq">
##FORMAT=<ID=LRR,Number=1,Type=Float,Description="Log R ratio">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sampleA
`

// buildVCF appends per-site BAF:LRR records (one sample) to the header.
func buildVCF(header string, chrom string, baf, lrr []float64) string {
	var b strings.Builder
	b.WriteString(header)
	for i := range baf {
		b.WriteString(chrom)
		b.WriteString("\t")
		// 1-based positions, 1000 bp apart.
		b.WriteString(itoa((i + 1) * 1000))
		b.WriteString("\t.\tA\tC\t.\tPASS\t.\tBAF:LRR\t")
		b.WriteString(ftoa(baf[i]))
		b.WriteString(":")
		b.WriteString(ftoa(lrr[i]))
		b.WriteString("\n")
	}
	return b.String()
}

func itoa(n int) string {
	if n == 0 {
		return "0"
	}
	var d []byte
	for n > 0 {
		d = append([]byte{byte('0' + n%10)}, d...)
		n /= 10
	}
	return string(d)
}

func ftoa(f float64) string {
	return strings.TrimRight(strings.TrimRight(
		formatFixed(f, 4), "0"), ".")
}

func formatFixed(f float64, prec int) string {
	scale := 1.0
	for i := 0; i < prec; i++ {
		scale *= 10
	}
	neg := f < 0
	if neg {
		f = -f
	}
	n := int(math.Round(f * scale))
	intPart := n / int(scale)
	frac := n % int(scale)
	s := itoa(intPart) + "."
	fs := itoa(frac)
	for len(fs) < prec {
		fs = "0" + fs
	}
	s += fs
	if neg {
		s = "-" + s
	}
	return s
}

// cnvOnly returns the data ("RG") rows of a CNV summary, dropping the
// header comment line.
func cnvRows(out string) []string {
	var rows []string
	for _, ln := range strings.Split(strings.TrimRight(out, "\n"), "\n") {
		if strings.HasPrefix(ln, "RG\t") {
			rows = append(rows, ln)
		}
	}
	return rows
}

// TestCNVCleanDiploidIsAllCN2 is the internal consistency check: a
// clean CN=2 region (every site a balanced het BAF≈0.5 and LRR≈0)
// must decode to a single all-normal CN2 region.
func TestCNVCleanDiploidIsAllCN2(t *testing.T) {
	n := 40
	baf := make([]float64, n)
	lrr := make([]float64, n)
	for i := range baf {
		// Alternate the three CN2 genotype clusters around their
		// peaks (RR≈0, RA≈0.5, AA≈1) — all consistent with diploid.
		switch i % 3 {
		case 0:
			baf[i] = 0.02
		case 1:
			baf[i] = 0.50
		default:
			baf[i] = 0.98
		}
		lrr[i] = 0.0
	}
	input := buildVCF(vcfHeaderSingle, "chr1", baf, lrr)
	var out bytes.Buffer
	if _, err := CNV(strings.NewReader(input), &out, CNVOptions{}); err != nil {
		t.Fatalf("CNV: %v", err)
	}
	rows := cnvRows(out.String())
	if len(rows) != 1 {
		t.Fatalf("clean diploid: expected 1 region, got %d:\n%s", len(rows), out.String())
	}
	f := strings.Split(rows[0], "\t")
	if f[1] != "chr1" || f[4] != "CN2" {
		t.Errorf("clean diploid: expected chr1 CN2, got %v", f)
	}
	if f[2] != "1000" || f[3] != itoa(n*1000) {
		t.Errorf("clean diploid: expected span 1000..%d, got %s..%s", n*1000, f[2], f[3])
	}
}

// TestCNVHomozygousDeletionIsCN0 hand-derives a CN0 call. Upstream's
// emission model assigns CN0 a fixed 0.5 weight on a no-call (missing
// BAF) site and exactly 0 on a site with a real BAF, so CN0 only
// decodes for runs of missing BAF. A single-sample input drops no-call
// sites entirely (parse_lrr_baf returns 0), exactly as upstream does;
// CN0 is therefore only reachable in paired mode, where a record
// survives as long as either sample has a usable BAF. This test puts a
// run of query no-calls against a control with valid diploid BAF.
func TestCNVHomozygousDeletionIsCN0(t *testing.T) {
	header := `##fileformat=VCFv4.2
##contig=<ID=chr1>
##FORMAT=<ID=BAF,Number=1,Type=Float,Description="">
##FORMAT=<ID=LRR,Number=1,Type=Float,Description="">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	queryA	ctrlB
`
	var b strings.Builder
	b.WriteString(header)
	pos := 0
	write := func(q string) {
		pos += 1000
		b.WriteString("chr1\t" + itoa(pos) + "\t.\tA\tC\t.\tPASS\t.\tBAF:LRR\t" +
			q + "\t0.5:0.0\n")
	}
	for i := 0; i < 30; i++ {
		write("0.5:0.0")
	}
	for i := 0; i < 80; i++ {
		write(".:.")
	}
	for i := 0; i < 30; i++ {
		write("0.5:0.0")
	}
	var out bytes.Buffer
	// A looser transition prior lets the long no-call run flip to CN0;
	// the default xy-prob is too rigid for a region of this length.
	if _, err := CNV(strings.NewReader(b.String()), &out,
		CNVOptions{QuerySample: "queryA", ControlSample: "ctrlB", XYProb: 1e-3}); err != nil {
		t.Fatalf("CNV: %v", err)
	}
	rows := cnvRows(out.String())
	// The query path should be CN2, CN0, CN2; the control stays CN2.
	sawCN0 := false
	for _, r := range rows {
		f := strings.Split(r, "\t")
		if f[4] == "CN0" {
			sawCN0 = true
		}
		if f[5] != "CN2" {
			t.Errorf("control should stay CN2, got %s in %v", f[5], f)
		}
	}
	if !sawCN0 {
		t.Errorf("expected a query CN0 region, got:\n%s", out.String())
	}
}

// TestCNVSingleCopyGainCN3 hand-derives a CN3 call. In a single-copy
// gain the heterozygous BAF cluster splits into two bands at 1/3 and
// 2/3 (cell_frac=1: 1/(2+1) and 2/(2+1)) plus a positive LRR shift
// (~0.30). A run of such sites flanked by normal CN2 must produce a
// CN3 region in the middle.
func TestCNVSingleCopyGainCN3(t *testing.T) {
	var b strings.Builder
	b.WriteString(vcfHeaderSingle)
	pos := 0
	write := func(baf, lrr float64) {
		pos += 1000
		b.WriteString("chr1\t" + itoa(pos) + "\t.\tA\tC\t.\tPASS\t.\tBAF:LRR\t" +
			ftoa(baf) + ":" + ftoa(lrr) + "\n")
	}
	for i := 0; i < 15; i++ {
		// Normal het: BAF 0.5, LRR 0.
		write(0.5, 0.0)
	}
	for i := 0; i < 25; i++ {
		// CN3: BAF alternates between the 1/3 and 2/3 bands, LRR ~0.30.
		if i%2 == 0 {
			write(1.0/3.0, 0.30)
		} else {
			write(2.0/3.0, 0.30)
		}
	}
	for i := 0; i < 15; i++ {
		write(0.5, 0.0)
	}
	var out bytes.Buffer
	if _, err := CNV(strings.NewReader(b.String()), &out, CNVOptions{}); err != nil {
		t.Fatalf("CNV: %v", err)
	}
	rows := cnvRows(out.String())
	gain := false
	for _, r := range rows {
		if strings.Contains(r, "\tCN3\t") {
			gain = true
		}
	}
	if !gain {
		t.Errorf("expected a CN3 region, got:\n%s", out.String())
	}
}

// TestCNVErrProbKnobTakesEffect confirms that --err-prob is now a
// load-bearing HMM knob: a huge error floor flattens every non-CN0
// emission so the decode collapses to a single region.
func TestCNVErrProbKnobTakesEffect(t *testing.T) {
	var b strings.Builder
	b.WriteString(vcfHeaderSingle)
	pos := 0
	write := func(baf, lrr float64) {
		pos += 1000
		b.WriteString("chr1\t" + itoa(pos) + "\t.\tA\tC\t.\tPASS\t.\tBAF:LRR\t" +
			ftoa(baf) + ":" + ftoa(lrr) + "\n")
	}
	for i := 0; i < 15; i++ {
		write(0.5, 0.0)
	}
	for i := 0; i < 25; i++ {
		if i%2 == 0 {
			write(1.0/3.0, 0.30)
		} else {
			write(2.0/3.0, 0.30)
		}
	}
	for i := 0; i < 15; i++ {
		write(0.5, 0.0)
	}
	input := b.String()

	var defOut bytes.Buffer
	if _, err := CNV(strings.NewReader(input), &defOut, CNVOptions{}); err != nil {
		t.Fatalf("CNV default: %v", err)
	}
	var floorOut bytes.Buffer
	// A dominating error floor swamps the BAF/LRR signal.
	if _, err := CNV(strings.NewReader(input), &floorOut, CNVOptions{ErrProb: 1e6}); err != nil {
		t.Fatalf("CNV err-prob: %v", err)
	}
	if defOut.String() == floorOut.String() {
		t.Errorf("--err-prob did not change the decode; knob is not load-bearing")
	}
	if len(cnvRows(floorOut.String())) != 1 {
		t.Errorf("err-prob floor should collapse to one region, got:\n%s", floorOut.String())
	}
}

// TestCNVXYProbKnobTakesEffect confirms --xy-prob influences the
// decode: a near-1 transition probability makes state switching free,
// which fragments the path relative to the rigid default.
func TestCNVXYProbKnobTakesEffect(t *testing.T) {
	n := 30
	baf := make([]float64, n)
	lrr := make([]float64, n)
	for i := range baf {
		// Borderline-noisy signal so the transition prior matters.
		if i < n/2 {
			baf[i] = 0.5
			lrr[i] = 0.0
		} else {
			baf[i] = 1.0 / 3.0
			lrr[i] = 0.25
		}
	}
	input := buildVCF(vcfHeaderSingle, "chr1", baf, lrr)
	var rigid, loose bytes.Buffer
	if _, err := CNV(strings.NewReader(input), &rigid, CNVOptions{XYProb: 1e-9}); err != nil {
		t.Fatalf("rigid: %v", err)
	}
	if _, err := CNV(strings.NewReader(input), &loose, CNVOptions{XYProb: 0.2}); err != nil {
		t.Fatalf("loose: %v", err)
	}
	if rigid.String() == loose.String() {
		t.Errorf("--xy-prob did not change the decode; knob is not load-bearing")
	}
}

// TestCNVPairedMode exercises the 16-state paired HMM and confirms the
// per-sample CN columns are emitted.
func TestCNVPairedMode(t *testing.T) {
	input := `##fileformat=VCFv4.2
##contig=<ID=chr1>
##FORMAT=<ID=BAF,Number=1,Type=Float,Description="">
##FORMAT=<ID=LRR,Number=1,Type=Float,Description="">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	queryA	ctrlB
`
	var b strings.Builder
	b.WriteString(input)
	for i := 0; i < 20; i++ {
		b.WriteString("chr1\t" + itoa((i+1)*1000) + "\t.\tA\tC\t.\tPASS\t.\tBAF:LRR\t0.5:0.0\t0.5:0.0\n")
	}
	var out bytes.Buffer
	n, err := CNV(strings.NewReader(b.String()), &out, CNVOptions{
		QuerySample: "queryA", ControlSample: "ctrlB",
	})
	if err != nil {
		t.Fatalf("CNV paired: %v", err)
	}
	if n != 1 {
		t.Fatalf("paired clean diploid: expected 1 region, got %d:\n%s", n, out.String())
	}
	if !strings.Contains(out.String(), "Copy number:queryA\t[6]Copy number:ctrlB") {
		t.Errorf("paired header missing both sample columns:\n%s", out.String())
	}
	f := strings.Split(cnvRows(out.String())[0], "\t")
	// RG chrom start end queryCN ctrlCN qual nsites nhets
	if f[4] != "CN2" || f[5] != "CN2" {
		t.Errorf("paired clean diploid: expected CN2/CN2, got %s/%s", f[4], f[5])
	}
}

// TestCNVSampleNarrowing exercises -s and verifies the right sample
// drives the call.
func TestCNVSampleNarrowing(t *testing.T) {
	input := `##fileformat=VCFv4.2
##contig=<ID=chr1>
##FORMAT=<ID=BAF,Number=1,Type=Float,Description="">
##FORMAT=<ID=LRR,Number=1,Type=Float,Description="">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	A	B	C
`
	var b strings.Builder
	b.WriteString(input)
	for i := 0; i < 10; i++ {
		b.WriteString("chr1\t" + itoa((i+1)*1000) +
			"\t.\tA\tC\t.\t.\t.\tBAF:LRR\t0.5:0.0\t0.5:0.0\t0.5:0.0\n")
	}
	var out bytes.Buffer
	n, err := CNV(strings.NewReader(b.String()), &out, CNVOptions{QuerySample: "B"})
	if err != nil {
		t.Fatalf("CNV: %v", err)
	}
	if n != 1 {
		t.Errorf("expected 1 row, got %d", n)
	}
	if !strings.Contains(out.String(), "Copy number:B") {
		t.Errorf("expected sample B in header, got:\n%s", out.String())
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

func TestCNVAFFileTargetsFilter(t *testing.T) {
	// --AF-file acts as a targets filter: sites whose CHROM:POS is absent
	// from the file are dropped before being added to the per-contig
	// observation buffer (vcfcnv.c uses the AF-file as the targets index).
	input := `##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="gt">
##FORMAT=<ID=BAF,Number=1,Type=Float,Description="baf">
##contig=<ID=chr1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	100	.	A	C	.	.	.	GT:BAF	0/1:0.5
chr1	200	.	A	C	.	.	.	GT:BAF	0/1:0.5
chr1	300	.	A	C	.	.	.	GT:BAF	0/1:0.5
`
	// AF file lists only positions 100 and 300.
	afPath := writeCNVTempFile(t, "chr1\t100\tA,C\t0.2\nchr1\t300\tA,C\t0.3\n")

	af, err := loadCNVAFFile(afPath)
	if err != nil {
		t.Fatalf("loadCNVAFFile: %v", err)
	}
	// Position 200 is absent → not listed.
	if _, listed := af.lookup(&vcf.Variant{Chrom: "chr1", Pos: 200, Ref: "A", Alt: []string{"C"}}); listed {
		t.Errorf("pos 200 should not be listed")
	}
	// Position 100 with matching alleles → listed, AF 0.2.
	if got, listed := af.lookup(&vcf.Variant{Chrom: "chr1", Pos: 100, Ref: "A", Alt: []string{"C"}}); !listed || got != 0.2 {
		t.Errorf("pos 100: got af=%v listed=%v, want 0.2 true", got, listed)
	}
	// Position 300 with NON-matching alleles → listed (targets), default AF.
	if got, listed := af.lookup(&vcf.Variant{Chrom: "chr1", Pos: 300, Ref: "A", Alt: []string{"T"}}); !listed || got != cnvNonrefAFDflt {
		t.Errorf("pos 300 allele mismatch: got af=%v listed=%v, want %v true", got, listed, cnvNonrefAFDflt)
	}

	// Multiallelic match: the full ALT vector (joined) must match the
	// file's REF,ALT1,ALT2 entry, not just the first ALT.
	maPath := writeCNVTempFile(t, "chr2\t500\tA,C,G\t0.4\n")
	maAF, err := loadCNVAFFile(maPath)
	if err != nil {
		t.Fatalf("loadCNVAFFile (multiallelic): %v", err)
	}
	if got, listed := maAF.lookup(&vcf.Variant{Chrom: "chr2", Pos: 500, Ref: "A", Alt: []string{"C", "G"}}); !listed || got != 0.4 {
		t.Errorf("multiallelic match: got af=%v listed=%v, want 0.4 true", got, listed)
	}

	// End-to-end: the run must succeed (no rejection) and consume only
	// the two listed sites.
	var out bytes.Buffer
	if _, err := CNV(strings.NewReader(input), &out, CNVOptions{QuerySample: "S1", AFFile: afPath}); err != nil {
		t.Fatalf("CNV with --AF-file: %v", err)
	}
}

// writeCNVTempFile writes content to a temp file and returns its path.
func writeCNVTempFile(t *testing.T, content string) string {
	t.Helper()
	p := filepath.Join(t.TempDir(), "afs.tab")
	if err := os.WriteFile(p, []byte(content), 0o644); err != nil {
		t.Fatalf("write %s: %v", p, err)
	}
	return p
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
	if _, ok := parseFloatField(map[string]string{}, "Y"); ok {
		t.Errorf("missing tag should return ok=false")
	}
}

// TestCNVBafFromAD pins the FORMAT/AD = REF,ALT fallback synthesis:
// when FORMAT/BAF isn't present, BAF must come from ALT/(REF+ALT).
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

// --- Unit tests for the ported HMM internals ------------------------

// TestCNVNormCDF checks the truncated-Gaussian mass against hand
// values: a peak centred at 0.5 keeps almost all its mass inside [0,1].
func TestCNVNormCDF(t *testing.T) {
	got := cnvNormCDF(0.5, 0.04)
	if got <= 0.99 || got > 1.0 {
		t.Errorf("cnvNormCDF(0.5,0.04) = %v, expected ~1", got)
	}
	// A peak centred at 0 loses half its mass below 0.
	half := cnvNormCDF(0.0, 0.04)
	if math.Abs(half-0.5) > 0.02 {
		t.Errorf("cnvNormCDF(0,0.04) = %v, expected ~0.5", half)
	}
}

// TestCNVInitTprobSingle checks the 4-state transition matrix: every
// column sums to 1, the diagonal carries the staying probability, and
// off-diagonals carry the xy-prob.
func TestCNVInitTprobSingle(t *testing.T) {
	xy := 1e-3
	mat, err := cnvInitTprob(cnvNStates, xy, 0.5)
	if err != nil {
		t.Fatalf("cnvInitTprob: %v", err)
	}
	for j := 0; j < cnvNStates; j++ {
		sum := 0.0
		for i := 0; i < cnvNStates; i++ {
			sum += mat[i*cnvNStates+j]
		}
		if math.Abs(sum-1) > 1e-12 {
			t.Errorf("column %d sums to %v, want 1", j, sum)
		}
		if d := mat[j*cnvNStates+j]; math.Abs(d-(1-xy*3)) > 1e-12 {
			t.Errorf("diagonal[%d] = %v, want %v", j, d, 1-xy*3)
		}
	}
}

// TestCNVInitTprobRejectsBadXY checks the upstream guard: an xy-prob so
// large that P(x|x) < P(x|y) is rejected.
func TestCNVInitTprobRejectsBadXY(t *testing.T) {
	if _, err := cnvInitTprob(cnvNStates, 0.5, 0.5); err == nil {
		t.Errorf("expected error for xy-prob too high")
	}
}

// TestCNVInitTprobPaired checks the 16-state paired transition matrix
// is column-stochastic.
func TestCNVInitTprobPaired(t *testing.T) {
	mat, err := cnvInitTprob(cnvNStates*cnvNStates, 1e-3, 0.5)
	if err != nil {
		t.Fatalf("cnvInitTprob paired: %v", err)
	}
	nd := cnvNStates * cnvNStates
	for j := 0; j < nd; j++ {
		sum := 0.0
		for i := 0; i < nd; i++ {
			sum += mat[i*nd+j]
		}
		if math.Abs(sum-1) > 1e-9 {
			t.Errorf("paired column %d sums to %v, want 1", j, sum)
		}
	}
}

// TestCNVInitIProbs checks the initial-state vector favours CN2 and
// normalises to 1.
func TestCNVInitIProbs(t *testing.T) {
	p := cnvInitIProbs(cnvNStates, 0.5)
	sum := 0.0
	for _, x := range p {
		sum += x
	}
	if math.Abs(sum-1) > 1e-12 {
		t.Errorf("iprobs sum = %v, want 1", sum)
	}
	if p[cnvCN2] <= p[cnvCN0] {
		t.Errorf("CN2 prior %v should exceed CN0 prior %v", p[cnvCN2], p[cnvCN0])
	}
}

// TestCNVSmoothData checks the moving-average smoother flattens a step:
// every value stays within the [min,max] range of the input and the
// hard 0->1 jump becomes a gradient.
func TestCNVSmoothData(t *testing.T) {
	dat := []float64{0, 0, 0, 0, 1, 1, 1, 1}
	cnvSmoothData(dat, 4)
	for i, x := range dat {
		if x < 0 || x > 1 {
			t.Errorf("smoothed value dat[%d]=%v out of [0,1]", i, x)
		}
	}
	if dat[3] == 0 && dat[4] == 1 {
		t.Errorf("smoother left a hard step: %v", dat)
	}
	// A win<=1 smoother is a no-op.
	noop := []float64{3, 1, 4, 1, 5}
	cnvSmoothData(noop, 1)
	for i, x := range noop {
		want := []float64{3, 1, 4, 1, 5}[i]
		if x != want {
			t.Errorf("win=1 should be a no-op, dat[%d]=%v want %v", i, x, want)
		}
	}
}

// TestCNVObservedProbMissingBAF checks the no-call emission: CN0 gets a
// fixed 0.5 and the rest split the remainder.
func TestCNVObservedProbMissingBAF(t *testing.T) {
	s := &cnvSample{baf: []float64{cnvMissingBAF}, lrr: []float64{0}, bafDev2: 0.0016, lrrDev2: 0.04}
	s.cellFrac = 1
	s.setGaussParams()
	s.setObservedProb(0, 0.76, 0.14, 0.098, 1.0, 0.2, 1e-4)
	if s.pobs[cnvCN0] != 0.5 {
		t.Errorf("missing BAF: pobs[CN0] = %v, want 0.5", s.pobs[cnvCN0])
	}
	rest := s.pobs[cnvCN1] + s.pobs[cnvCN2] + s.pobs[cnvCN3]
	if math.Abs(rest-0.5) > 1e-9 {
		t.Errorf("missing BAF: non-CN0 sum = %v, want 0.5", rest)
	}
}

// TestCNVObservedProbDiploidPeak checks that a balanced het BAF makes
// CN2 the most likely state.
func TestCNVObservedProbDiploidPeak(t *testing.T) {
	s := &cnvSample{baf: []float64{0.5}, lrr: []float64{0.0}, bafDev2: 0.0016, lrrDev2: 0.04}
	s.cellFrac = 1
	s.setGaussParams()
	s.setObservedProb(0, 0.76, 0.14, 0.098, 1.0, 0.2, 1e-4)
	if s.pobs[cnvCN2] <= s.pobs[cnvCN1] || s.pobs[cnvCN2] <= s.pobs[cnvCN3] {
		t.Errorf("balanced het BAF should favour CN2: %v", s.pobs)
	}
}
