package vcftools

import (
	"bufio"
	"fmt"
	"io"
	"math"
	"strconv"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// ldSite stores the per-sample genotype information needed to compute LD
// between this site and any later site. genoCounts[i] is the number of ALT
// alleles for sample i at this site, or -1 if the sample is missing/multi-
// allelic and should be ignored. For --hap-r2 we additionally need the two
// haplotype alleles per sample (only meaningful when the genotype is phased);
// hapA[i] / hapB[i] hold those values, or -1 when unusable.
//
// We only consider the reference allele ("0") and the first alternate allele
// ("1") as in the upstream vcftools behaviour for these statistics; any sample
// carrying a higher-index alternate allele is treated as missing for the pair.
type ldSite struct {
	chrom string
	pos   int

	// chromIdx is the 1-based index of this site within its chromosome (after
	// filtering). Used to compute --ld-window / --ld-window-min SNP distances.
	chromIdx int

	// genoCounts[i] in {0,1,2} for diploid genotypes restricted to the first
	// two alleles, or -1 for missing/skip.
	genoCounts []int

	// hapA[i], hapB[i] in {0,1} for phased diploid genotypes restricted to the
	// first two alleles, or -1 when not usable for --hap-r2 (missing or
	// unphased).
	hapA []int
	hapB []int
}

// ldWriter wraps a per-output-file writer with a header.
type ldWriter struct {
	f   io.WriteCloser
	w   *bufio.Writer
	err error
}

func newLDWriter(path, header string) (*ldWriter, error) {
	f, err := iohelper.OpenWriter(path)
	if err != nil {
		return nil, fmt.Errorf("opening %s: %w", path, err)
	}
	w := bufio.NewWriter(f)
	if _, err := w.WriteString(header); err != nil {
		f.Close()
		return nil, fmt.Errorf("writing header to %s: %w", path, err)
	}
	return &ldWriter{f: f, w: w}, nil
}

func (l *ldWriter) writeLine(line string) {
	if l == nil || l.err != nil {
		return
	}
	if _, err := l.w.WriteString(line); err != nil {
		l.err = err
	}
}

func (l *ldWriter) close() error {
	if l == nil {
		return nil
	}
	if err := l.w.Flush(); err != nil {
		l.f.Close()
		return err
	}
	if err := l.f.Close(); err != nil {
		return err
	}
	return l.err
}

// extractLDSite turns a vcf.Variant into the compact ldSite representation. It
// returns (nil, false) for variants that cannot participate in any LD pair
// (no samples, or no alternate allele at all). Genotypes referencing any
// alternate beyond the first ALT are recorded as missing (-1) so the pair
// computation skips them, matching upstream "first ALT only" behaviour for
// multi-allelic sites.
func extractLDSite(v *vcf.Variant) (*ldSite, bool) {
	if len(v.Samples) == 0 || len(v.Alt) == 0 {
		return nil, false
	}
	s := &ldSite{
		chrom:      v.Chrom,
		pos:        v.Pos,
		genoCounts: make([]int, len(v.Samples)),
		hapA:       make([]int, len(v.Samples)),
		hapB:       make([]int, len(v.Samples)),
	}
	for i, sample := range v.Samples {
		gt, ok := sample.Data["GT"]
		if !ok {
			s.genoCounts[i] = -1
			s.hapA[i] = -1
			s.hapB[i] = -1
			continue
		}
		s.genoCounts[i], s.hapA[i], s.hapB[i] = parseGTForLD(gt)
	}
	return s, true
}

// parseGTForLD parses a VCF GT field for the LD statistics. It returns
// (genoCount, hapA, hapB) where:
//   - genoCount is the diploid ALT-allele count restricted to alleles {0,1}
//     (0, 1, or 2); -1 if the genotype is missing or references a higher-
//     index allele.
//   - hapA, hapB are the two haplotype alleles (0 or 1) when the genotype is
//     phased ("a|b") and both alleles are in {0,1}; -1 otherwise (missing or
//     unphased or out-of-range).
//
// This mirrors the upstream choice to use only the first ALT allele for
// pairwise LD on multi-allelic sites and to require phasing for --hap-r2.
func parseGTForLD(gt string) (genoCount, hapA, hapB int) {
	if gt == "" || gt == "." || gt == "./." || gt == ".|." {
		return -1, -1, -1
	}
	// Detect the separator: '|' means phased, '/' means unphased.
	phased := false
	sep := -1
	for i, r := range gt {
		if r == '|' {
			phased = true
			sep = i
			break
		}
		if r == '/' {
			sep = i
			break
		}
	}
	if sep < 0 {
		// Haploid call (e.g. "0" or "1"). Treat as a single chromosome:
		// genoCount unusable for diploid LD (skip), hap unusable.
		return -1, -1, -1
	}
	left := gt[:sep]
	right := gt[sep+1:]
	a, aOK := parseAlleleForLD(left)
	b, bOK := parseAlleleForLD(right)
	if !aOK || !bOK {
		return -1, -1, -1
	}
	if a < 0 || b < 0 || a > 1 || b > 1 {
		// Missing or allele beyond the first ALT: skip for this pair.
		return -1, -1, -1
	}
	gc := a + b
	if phased {
		return gc, a, b
	}
	return gc, -1, -1
}

// parseAlleleForLD converts a single allele token. Returns (-1, true) for
// missing, the integer index for numeric tokens, and (0, false) for tokens
// that don't parse.
func parseAlleleForLD(tok string) (int, bool) {
	if tok == "." {
		return -1, true
	}
	n, err := strconv.Atoi(tok)
	if err != nil {
		return 0, false
	}
	return n, true
}

// computeGenoR2 returns (nIndv, r2, ok) for two sites using diploid
// allele-count vectors. r2 = cov(g,h)^2 / (var(g)*var(h)) over samples that
// are non-missing at both sites. ok is false when n<2 or either variance is
// zero (monomorphic at one site within the shared sample set).
func computeGenoR2(a, b *ldSite) (int, float64, bool) {
	if len(a.genoCounts) != len(b.genoCounts) {
		return 0, 0, false
	}
	n := 0
	var sumA, sumB, sumAB, sumAA, sumBB float64
	for i := range a.genoCounts {
		ga := a.genoCounts[i]
		gb := b.genoCounts[i]
		if ga < 0 || gb < 0 {
			continue
		}
		fa := float64(ga)
		fb := float64(gb)
		sumA += fa
		sumB += fb
		sumAB += fa * fb
		sumAA += fa * fa
		sumBB += fb * fb
		n++
	}
	if n < 2 {
		return n, 0, false
	}
	fn := float64(n)
	meanA := sumA / fn
	meanB := sumB / fn
	cov := sumAB/fn - meanA*meanB
	varA := sumAA/fn - meanA*meanA
	varB := sumBB/fn - meanB*meanB
	if varA <= 0 || varB <= 0 {
		return n, 0, false
	}
	r2 := (cov * cov) / (varA * varB)
	return n, r2, true
}

// computeHapR2 returns (nChr, r2, D, Dprime, ok) for two phased diploid sites.
// nChr counts haplotypes (= 2 * number of individuals phased and non-missing
// at both sites). Skips if nChr<2 or either site is monomorphic in the shared
// haplotype set.
//
// Math: we pick the reference allele ("0") at each site as the haplotype of
// interest, so pA = freq("0") at site A, pB = freq("0") at site B, and
// pAB = freq(haplotype "0|0"). D = pAB - pA*pB. r² = D² / (pA(1-pA)pB(1-pB)).
// Dmax = min(pA*(1-pB), (1-pA)*pB) when D >= 0 and
//
//	min(pA*pB, (1-pA)*(1-pB)) when D < 0, then Dprime = D / Dmax.
//
// Using the "0" allele rather than the minor allele matches the documented
// formula; r² and |Dprime| are invariant to that choice anyway.
func computeHapR2(a, b *ldSite) (nChr int, r2, D, Dprime float64, ok bool) {
	if len(a.hapA) != len(b.hapA) {
		return 0, 0, 0, 0, false
	}
	var n00, n0, m0 int // n_AB, count of "0" at A, count of "0" at B
	for i := range a.hapA {
		for _, side := range []int{0, 1} {
			var ai, bi int
			if side == 0 {
				ai = a.hapA[i]
				bi = b.hapA[i]
			} else {
				ai = a.hapB[i]
				bi = b.hapB[i]
			}
			if ai < 0 || bi < 0 {
				continue
			}
			nChr++
			if ai == 0 {
				n0++
			}
			if bi == 0 {
				m0++
			}
			if ai == 0 && bi == 0 {
				n00++
			}
		}
	}
	if nChr < 2 {
		return nChr, 0, 0, 0, false
	}
	fn := float64(nChr)
	pA := float64(n0) / fn
	pB := float64(m0) / fn
	if pA <= 0 || pA >= 1 || pB <= 0 || pB >= 1 {
		return nChr, 0, 0, 0, false
	}
	pAB := float64(n00) / fn
	D = pAB - pA*pB
	var dMax float64
	if D >= 0 {
		dMax = pA * (1 - pB)
		if (1-pA)*pB < dMax {
			dMax = (1 - pA) * pB
		}
	} else {
		dMax = pA * pB
		if (1-pA)*(1-pB) < dMax {
			dMax = (1 - pA) * (1 - pB)
		}
	}
	if dMax == 0 {
		return nChr, 0, 0, 0, false
	}
	Dprime = D / dMax
	denom := pA * (1 - pA) * pB * (1 - pB)
	if denom <= 0 {
		return nChr, 0, 0, 0, false
	}
	r2 = (D * D) / denom
	return nChr, r2, D, Dprime, true
}

// withinLDWindow returns true when (a, b) satisfy all window constraints. a
// must come before b in the input stream; both must be on the same chromosome
// (vcftools does not emit LD across chromosomes for --geno-r2 / --hap-r2).
func withinLDWindow(a, b *ldSite, params *Params) bool {
	if a.chrom != b.chrom {
		return false
	}
	snpDist := b.chromIdx - a.chromIdx
	if snpDist < 0 {
		snpDist = -snpDist
	}
	bpDist := b.pos - a.pos
	if bpDist < 0 {
		bpDist = -bpDist
	}
	if params.LDWindow > 0 && snpDist > params.LDWindow {
		return false
	}
	if params.LDWindowBp > 0 && bpDist > params.LDWindowBp {
		return false
	}
	if params.LDWindowMin > 0 && snpDist < params.LDWindowMin {
		return false
	}
	if params.LDWindowBpMin > 0 && bpDist < params.LDWindowBpMin {
		return false
	}
	return true
}

// ldPositionAllowed implements --geno-r2-positions / --hap-r2-positions
// restriction: emit the pair only if at least one of the two sites is in the
// supplied positions set. When restricted is false, all pairs are allowed.
func ldPositionAllowed(a, b *ldSite, pos positionSet, restricted bool) bool {
	if !restricted {
		return true
	}
	if pos == nil {
		return false
	}
	if chromPos, ok := pos[a.chrom]; ok && chromPos[a.pos] {
		return true
	}
	if chromPos, ok := pos[b.chrom]; ok && chromPos[b.pos] {
		return true
	}
	return false
}

// ldRunner accumulates filtered sites and emits pairwise LD output for the
// currently-buffered chromosome. Because the main filtering pass is already
// streaming variants in order, we feed each kept variant into the runner; the
// runner emits all valid pairs (latest site × prior sites within the window)
// as soon as the latest site is appended, then prunes the window.
type ldRunner struct {
	params  *Params
	genoW   *ldWriter
	hapW    *ldWriter
	genoPos positionSet
	hapPos  positionSet

	wantGeno bool
	wantHap  bool

	window   []*ldSite
	chromIdx int
	chrom    string
}

// newLDRunner opens the output writers and returns a runner ready to consume
// per-site genotype data. Caller MUST call close() to flush output.
func newLDRunner(params *Params) (*ldRunner, error) {
	r := &ldRunner{params: params}
	r.wantGeno = params.GenoR2 || params.GenoR2Positions != ""
	r.wantHap = params.HapR2 || params.HapR2Positions != ""
	if params.GenoR2Positions != "" {
		p, err := loadPositions(params.GenoR2Positions)
		if err != nil {
			return nil, fmt.Errorf("loading --geno-r2-positions: %w", err)
		}
		r.genoPos = p
	}
	if params.HapR2Positions != "" {
		p, err := loadPositions(params.HapR2Positions)
		if err != nil {
			return nil, fmt.Errorf("loading --hap-r2-positions: %w", err)
		}
		r.hapPos = p
	}
	if r.wantGeno {
		w, err := newLDWriter(params.OutPrefix+".geno.ld", "CHR\tPOS1\tPOS2\tN_INDV\tR^2\n")
		if err != nil {
			return nil, err
		}
		r.genoW = w
	}
	if r.wantHap {
		w, err := newLDWriter(params.OutPrefix+".hap.ld", "CHR\tPOS1\tPOS2\tN_CHR\tR^2\tD\tDprime\n")
		if err != nil {
			// Best-effort close of the other writer to avoid leaking the fd.
			_ = r.genoW.close()
			return nil, err
		}
		r.hapW = w
	}
	return r, nil
}

// addVariant feeds a variant (already passed through the same filtering as the
// main pass) into the LD runner.
func (r *ldRunner) addVariant(v *vcf.Variant) {
	if r == nil || (!r.wantGeno && !r.wantHap) {
		return
	}
	site, ok := extractLDSite(v)
	if !ok {
		return
	}
	if site.chrom != r.chrom {
		r.window = r.window[:0]
		r.chromIdx = 0
		r.chrom = site.chrom
	}
	r.chromIdx++
	site.chromIdx = r.chromIdx

	for _, prev := range r.window {
		if !withinLDWindow(prev, site, r.params) {
			continue
		}
		if r.wantGeno {
			if ldPositionAllowed(prev, site, r.genoPos, r.params.GenoR2Positions != "") {
				n, r2, ok := computeGenoR2(prev, site)
				if ok && r2 >= r.params.MinR2 {
					r.genoW.writeLine(fmt.Sprintf("%s\t%d\t%d\t%d\t%s\n",
						prev.chrom, prev.pos, site.pos, n, formatCppDefault(r2)))
				}
			}
		}
		if r.wantHap {
			if ldPositionAllowed(prev, site, r.hapPos, r.params.HapR2Positions != "") {
				n, r2, D, Dp, ok := computeHapR2(prev, site)
				if ok && r2 >= r.params.MinR2 {
					r.hapW.writeLine(fmt.Sprintf("%s\t%d\t%d\t%d\t%s\t%s\t%s\n",
						prev.chrom, prev.pos, site.pos, n,
						formatCppDefault(r2), formatCppDefault(D), formatCppDefault(Dp)))
				}
			}
		}
	}

	r.window = append(r.window, site)
	r.pruneWindow(site)
}

// pruneWindow drops sites from the front of the window when *both* configured
// maximum-distance constraints have been exceeded relative to the latest
// site. Future sites can only be further along the chromosome so once a site
// drops out of every window it stays out.
func (r *ldRunner) pruneWindow(latest *ldSite) {
	cut := 0
	for cut < len(r.window) {
		prev := r.window[cut]
		snpDist := latest.chromIdx - prev.chromIdx
		bpDist := latest.pos - prev.pos
		if bpDist < 0 {
			bpDist = -bpDist
		}
		// A prev site is irrelevant for all future pairs iff at least one
		// configured maximum is already exceeded. Use the AND of "exceeded"
		// flags so that with neither maximum set we never prune (keep all,
		// matching unbounded default).
		exceededSNP := r.params.LDWindow > 0 && snpDist >= r.params.LDWindow
		exceededBp := r.params.LDWindowBp > 0 && bpDist >= r.params.LDWindowBp
		neitherSet := r.params.LDWindow == 0 && r.params.LDWindowBp == 0
		if neitherSet {
			break
		}
		if !exceededSNP && !exceededBp {
			break
		}
		// If only one window is configured, "exceeded" of that one is enough.
		// If both are configured, require the *less* restrictive (or rather:
		// either) to be exceeded to be safe — we only prune when the prev
		// can never re-enter via the other window either.
		if r.params.LDWindow > 0 && r.params.LDWindowBp > 0 {
			if !(exceededSNP && exceededBp) {
				break
			}
		}
		cut++
	}
	if cut == 0 {
		return
	}
	r.window = append(r.window[:0], r.window[cut:]...)
}

// close flushes and closes all open LD output files. Returns the first error
// encountered.
func (r *ldRunner) close() error {
	if r == nil {
		return nil
	}
	var firstErr error
	if err := r.genoW.close(); err != nil {
		firstErr = err
	}
	if err := r.hapW.close(); err != nil && firstErr == nil {
		firstErr = err
	}
	return firstErr
}

// computeGenoChiSq returns (nIndv, chi2, df, pValue, ok) for a 3×3 contingency
// table of diploid allele-count combinations (g_A in {0,1,2} × g_B in {0,1,2})
// across samples that are non-missing at both sites.
//
// We use Pearson's chi-square test of association on the 3×3 table. The
// degrees of freedom is (rows-1)*(cols-1) restricted to row/column categories
// that actually occur in the data: this matches upstream behaviour for the
// `--geno-chisq` output where monomorphic columns/rows aren't counted toward
// df. Expected counts are computed from the marginal totals.
//
// ok is false when n<2, or when df==0 (one of the sites is monomorphic across
// the shared samples, so no association can be tested).
//
// The p-value uses the regularised upper incomplete gamma function for the
// chi-square distribution with df degrees of freedom.
func computeGenoChiSq(a, b *ldSite) (n int, chi2 float64, df int, pValue float64, ok bool) {
	if len(a.genoCounts) != len(b.genoCounts) {
		return 0, 0, 0, 0, false
	}
	var table [3][3]int
	for i := range a.genoCounts {
		ga := a.genoCounts[i]
		gb := b.genoCounts[i]
		if ga < 0 || gb < 0 {
			continue
		}
		if ga > 2 || gb > 2 {
			// Defensive: extractLDSite caps at 2.
			continue
		}
		table[ga][gb]++
		n++
	}
	if n < 2 {
		return n, 0, 0, 0, false
	}
	var rowTot [3]int
	var colTot [3]int
	for i := 0; i < 3; i++ {
		for j := 0; j < 3; j++ {
			rowTot[i] += table[i][j]
			colTot[j] += table[i][j]
		}
	}
	// Count non-empty rows and columns: degrees of freedom = (R-1)*(C-1).
	rows, cols := 0, 0
	for i := 0; i < 3; i++ {
		if rowTot[i] > 0 {
			rows++
		}
		if colTot[i] > 0 {
			cols++
		}
	}
	if rows < 2 || cols < 2 {
		return n, 0, 0, 0, false
	}
	df = (rows - 1) * (cols - 1)
	fn := float64(n)
	for i := 0; i < 3; i++ {
		if rowTot[i] == 0 {
			continue
		}
		for j := 0; j < 3; j++ {
			if colTot[j] == 0 {
				continue
			}
			expected := float64(rowTot[i]) * float64(colTot[j]) / fn
			if expected <= 0 {
				continue
			}
			diff := float64(table[i][j]) - expected
			chi2 += diff * diff / expected
		}
	}
	pValue = chiSquareSurvival(chi2, df)
	return n, chi2, df, pValue, true
}

// chiSquareSurvival returns P(X >= x) for X ~ chi-square(df). df > 0, x >= 0.
// Implemented via the regularised upper incomplete gamma function
// Q(df/2, x/2). Returns 1.0 for x<=0 and 0.0 for df<=0 (defensive).
func chiSquareSurvival(x float64, df int) float64 {
	if df <= 0 {
		return 0
	}
	if x <= 0 {
		return 1
	}
	return regUpperGamma(float64(df)/2.0, x/2.0)
}

// regUpperGamma returns Q(s, x) = 1 - P(s, x) where P is the regularised
// lower incomplete gamma function. s>0, x>=0. Uses a series expansion for
// x < s+1 and a continued-fraction expansion otherwise, following Numerical
// Recipes §6.2.
func regUpperGamma(s, x float64) float64 {
	if x < 0 || s <= 0 {
		return math.NaN()
	}
	if x == 0 {
		return 1
	}
	if x < s+1 {
		return 1 - regLowerSeries(s, x)
	}
	return regUpperContFrac(s, x)
}

func regLowerSeries(s, x float64) float64 {
	// P(s,x) = e^(-x) * x^s / Γ(s+1) * Σ x^n / (s+1)(s+2)...(s+n)  (n=0..)
	gln, _ := math.Lgamma(s)
	ap := s
	sum := 1.0 / s
	del := sum
	for i := 1; i < 1000; i++ {
		ap++
		del *= x / ap
		sum += del
		if math.Abs(del) < math.Abs(sum)*1e-15 {
			break
		}
	}
	return sum * math.Exp(-x+s*math.Log(x)-gln)
}

func regUpperContFrac(s, x float64) float64 {
	// Q(s,x) via Lentz's algorithm on the continued fraction
	//   1/(x+1-s - 1·(1-s)/(x+3-s - 2·(2-s)/(x+5-s - ...)))
	gln, _ := math.Lgamma(s)
	fpmin := 1e-300
	b := x + 1 - s
	c := 1 / fpmin
	d := 1 / b
	h := d
	for i := 1; i < 1000; i++ {
		fi := float64(i)
		an := -fi * (fi - s)
		b += 2
		d = an*d + b
		if math.Abs(d) < fpmin {
			d = fpmin
		}
		c = b + an/c
		if math.Abs(c) < fpmin {
			c = fpmin
		}
		d = 1 / d
		del := d * c
		h *= del
		if math.Abs(del-1) < 1e-15 {
			break
		}
	}
	return math.Exp(-x+s*math.Log(x)-gln) * h
}

// interchromLDRunner buffers all per-chromosome sites for --interchrom-geno-r2
// / --interchrom-hap-r2. Because pairs cross chromosomes, we can't stream:
// upstream batches per-chrom and at end emits all (chrI, chrJ) pairs with
// I != J. We follow the same approach but limit the output to one row per
// distinct unordered pair.
type interchromLDRunner struct {
	params *Params

	wantGeno    bool
	wantHap     bool
	wantChiSq   bool
	genoW       *ldWriter
	hapW        *ldWriter
	chiSqW      *ldWriter
	sites       []*ldSite
	chromOrder  []string
	siteByChrom map[string][]*ldSite
}

// newInterchromLDRunner allocates output writers as requested.
func newInterchromLDRunner(params *Params) (*interchromLDRunner, error) {
	r := &interchromLDRunner{
		params:      params,
		wantGeno:    params.InterchromGenoR2,
		wantHap:     params.InterchromHapR2,
		wantChiSq:   params.GenoChiSq,
		siteByChrom: make(map[string][]*ldSite),
	}
	if r.wantGeno {
		w, err := newLDWriter(params.OutPrefix+".interchrom.geno.ld",
			"CHR1\tPOS1\tCHR2\tPOS2\tN_INDV\tR^2\n")
		if err != nil {
			return nil, err
		}
		r.genoW = w
	}
	if r.wantHap {
		w, err := newLDWriter(params.OutPrefix+".interchrom.hap.ld",
			"CHR1\tPOS1\tCHR2\tPOS2\tN_CHR\tR^2\tD\tDprime\n")
		if err != nil {
			_ = r.genoW.close()
			return nil, err
		}
		r.hapW = w
	}
	if r.wantChiSq {
		w, err := newLDWriter(params.OutPrefix+".geno.chisq",
			"CHR1\tPOS1\tCHR2\tPOS2\tN_INDV\tCHI^2\tDF\tP-VALUE\n")
		if err != nil {
			_ = r.genoW.close()
			_ = r.hapW.close()
			return nil, err
		}
		r.chiSqW = w
	}
	return r, nil
}

// addVariant buffers a variant for later pairwise computation.
func (r *interchromLDRunner) addVariant(v *vcf.Variant) {
	if r == nil || (!r.wantGeno && !r.wantHap && !r.wantChiSq) {
		return
	}
	site, ok := extractLDSite(v)
	if !ok {
		return
	}
	if _, seen := r.siteByChrom[site.chrom]; !seen {
		r.chromOrder = append(r.chromOrder, site.chrom)
	}
	r.siteByChrom[site.chrom] = append(r.siteByChrom[site.chrom], site)
}

// flush emits all interchromosomal pairs and, when --geno-chisq is set,
// emits chi-square for every (cross-chrom AND same-chrom) pair, matching
// upstream's all-pairs behaviour for --geno-chisq.
func (r *interchromLDRunner) flush() error {
	if r == nil || (!r.wantGeno && !r.wantHap && !r.wantChiSq) {
		return nil
	}
	// Inter-chromosomal pairs (cross-chrom).
	for i := 0; i < len(r.chromOrder); i++ {
		for j := i + 1; j < len(r.chromOrder); j++ {
			ci := r.chromOrder[i]
			cj := r.chromOrder[j]
			for _, a := range r.siteByChrom[ci] {
				for _, b := range r.siteByChrom[cj] {
					r.emitPair(a, b)
				}
			}
		}
	}
	// Same-chrom pairs for --geno-chisq only. For the interchrom-* outputs
	// upstream emits only cross-chromosome pairs.
	if r.wantChiSq {
		for _, chr := range r.chromOrder {
			sites := r.siteByChrom[chr]
			for i := 0; i < len(sites); i++ {
				for j := i + 1; j < len(sites); j++ {
					r.emitChiSqPair(sites[i], sites[j])
				}
			}
		}
	}
	return nil
}

func (r *interchromLDRunner) emitPair(a, b *ldSite) {
	if r.wantGeno {
		n, r2, ok := computeGenoR2(a, b)
		if ok && r2 >= r.params.MinR2 {
			r.genoW.writeLine(fmt.Sprintf("%s\t%d\t%s\t%d\t%d\t%s\n",
				a.chrom, a.pos, b.chrom, b.pos, n, formatCppDefault(r2)))
		}
	}
	if r.wantHap {
		n, r2, D, Dp, ok := computeHapR2(a, b)
		if ok && r2 >= r.params.MinR2 {
			r.hapW.writeLine(fmt.Sprintf("%s\t%d\t%s\t%d\t%d\t%s\t%s\t%s\n",
				a.chrom, a.pos, b.chrom, b.pos, n,
				formatCppDefault(r2), formatCppDefault(D), formatCppDefault(Dp)))
		}
	}
	if r.wantChiSq {
		r.emitChiSqPair(a, b)
	}
}

func (r *interchromLDRunner) emitChiSqPair(a, b *ldSite) {
	n, chi2, df, p, ok := computeGenoChiSq(a, b)
	if !ok {
		return
	}
	r.chiSqW.writeLine(fmt.Sprintf("%s\t%d\t%s\t%d\t%d\t%s\t%d\t%s\n",
		a.chrom, a.pos, b.chrom, b.pos, n,
		formatCppDefault(chi2), df, formatCppDefault(p)))
}

func (r *interchromLDRunner) close() error {
	if r == nil {
		return nil
	}
	var firstErr error
	if err := r.genoW.close(); err != nil {
		firstErr = err
	}
	if err := r.hapW.close(); err != nil && firstErr == nil {
		firstErr = err
	}
	if err := r.chiSqW.close(); err != nil && firstErr == nil {
		firstErr = err
	}
	return firstErr
}
