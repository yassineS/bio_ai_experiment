package bcftools

import (
	"bufio"
	"fmt"
	"io"
	"math"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bcf"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// StatsOptions controls the behaviour of Stats / StatsFile. It mirrors the
// `bcftools stats` upstream flags closely.
type StatsOptions struct {
	Samples         []string // -s
	SamplesGiven    bool     // whether -s/-S was supplied at all (gates PSC/PSI/HWE + DP genotype binning)
	SamplesFile     string   // -S
	Regions         []string // -r
	RegionsFile     string   // -R
	Targets         []string // -t
	TargetsFile     string   // -T
	IncludeExpr     string   // -i
	ExcludeExpr     string   // -e
	ApplyFilters    []string // -f
	DepthMin        int      // -d MIN,MAX,STEP
	DepthMax        int
	DepthStep       int
	AFBins          []float64 // -a; if nil we use the upstream default
	Collapse        string    // -c
	FirstAlleleOnly bool      // -1
	SplitByID       bool      // -I/--split-by-ID: per-section split into known vs novel
	AFTag           string    // --af-tag
	InputFile       string    // for the header line
}

// defaultAFBins keeps the legacy port's 11-bin layout for the internal
// `afSNPs`/`afTs`/`afTv`/`afNonS` accumulators that unit tests assert on.
// The textual `AF` section emitted by writeAF instead uses upstream's
// 101-bin scheme (mAFBins / computeUpstreamAFBins) so the output is
// byte-for-byte parity with `bcftools stats`.
var defaultAFBins = []float64{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99, 1.0}

// mAFBins mirrors upstream `args->m_af` (vcfstats.c:448) — 101 bins by
// default. Bin 0 is reserved for singletons (AC==1) and bins 1..mAFBins-1
// are populated via `floor(af*(mAFBins-2))+1`.
const mAFBins = 101

// naFHWE mirrors upstream `args->naf_hwe` (vcfstats.c:473) — 100
// het-fraction bins per allele-frequency bucket.
const naFHWE = 100

// statsResult holds all the accumulator state used while streaming records.
// It is exposed as a value-style struct so tests can assert on it directly
// without re-parsing the formatted output.
type statsResult struct {
	opts StatsOptions

	samples       []string
	sampleIndex   map[string]int
	headerSamples []string

	// SN counters
	numRecords int
	numNoALTs  int
	numSNPs    int
	numMNPs    int
	numIndels  int
	numOthers  int
	numMA      int
	numMASNP   int

	// AF accumulators (legacy 11-bin layout, exposed for unit tests)
	afBins []float64
	afSNPs []int
	afTs   []int
	afTv   []int
	afNonS []int

	// QUAL accumulators keyed by integer bucket index (floor(QUAL)).
	qualSNPs map[int]int
	qualTs   map[int]int
	qualTv   map[int]int
	qualNonS map[int]int

	// IDD: indel length distribution (length, count) — len < 0 = deletion.
	indelLen map[int]int

	// ST: substitution-type counts (REF>ALT for SNPs).
	subst map[string]int

	// DP: depth distribution.
	dpSites     map[int]int // bin -> sites
	dpGTs       map[int]int // bin -> per-sample GTs at that depth bucket
	dpTotalSite int         // total sites that contributed
	dpTotalGT   int

	// PSC / PSI: per-sample counters indexed by sampleIndex.
	pscNRefHom    []int
	pscNNonRefHom []int
	pscNHets      []int
	pscNTs        []int
	pscNTv        []int
	pscNIndels    []int
	pscDepthSum   []int
	pscDepthN     []int
	pscNSingleton []int // smpl_sngl in vcfstats.c
	pscNHapRef    []int // smpl_hapRef
	pscNHapAlt    []int // smpl_hapAlt
	pscNMissing   []int // smpl_missing

	// PSI per-sample indel het/hom ins/del counters
	// (vcfstats.c:1077-1083 / 1816-1817).
	psiInsHets []int
	psiDelHets []int
	psiInsHoms []int
	psiDelHoms []int

	// HWE: AF -> (n_obs, chi-square sum). Retained for the legacy
	// `TestStatsHWEChiSquare` unit test; the textual HWE output uses afHWE.
	hweObs    map[int]int
	hweChiSum map[int]float64

	// Upstream-parity accumulators feeding writeAF / writeQUAL / writeHWE.
	afSnpsUS  [mAFBins]int         // SNP count per upstream AF bin
	afTsUS    [mAFBins]int         // transition count per upstream AF bin
	afTvUS    [mAFBins]int         // transversion count per upstream AF bin
	afIndUS   [mAFBins]int         // indel count per upstream AF bin
	afRepNAUS [mAFBins]int         // upstream `af_repeats[2]` (na column)
	qualTsUS  map[int]int          // ts count keyed by iqual=1+int(qual*10)
	qualTvUS  map[int]int          // tv count keyed by iqual
	qualIndS  map[int]int          // indel count keyed by iqual
	tsAlt1    int                  // upstream `stats->ts_alt1` (1st-ALT transitions)
	tvAlt1    int                  // upstream `stats->tv_alt1` (1st-ALT transversions)
	afHWE     [mAFBins][naFHWE]int // [afBin][hetFreqBin] -> n records
}

// newStatsResult prepares accumulators sized to the requested sample set.
func newStatsResult(opts StatsOptions, headerSamples []string) *statsResult {
	r := &statsResult{
		opts:          opts,
		headerSamples: headerSamples,
		afBins:        opts.AFBins,
		qualSNPs:      make(map[int]int),
		qualTs:        make(map[int]int),
		qualTv:        make(map[int]int),
		qualNonS:      make(map[int]int),
		indelLen:      make(map[int]int),
		subst:         make(map[string]int),
		dpSites:       make(map[int]int),
		dpGTs:         make(map[int]int),
		hweObs:        make(map[int]int),
		hweChiSum:     make(map[int]float64),
		qualTsUS:      make(map[int]int),
		qualTvUS:      make(map[int]int),
		qualIndS:      make(map[int]int),
	}
	if len(r.afBins) == 0 {
		r.afBins = append([]float64{}, defaultAFBins...)
	}
	nBins := len(r.afBins) - 1
	if nBins < 1 {
		nBins = 1
	}
	r.afSNPs = make([]int, nBins)
	r.afTs = make([]int, nBins)
	r.afTv = make([]int, nBins)
	r.afNonS = make([]int, nBins)

	r.samples = filterSampleSet(headerSamples, opts.SamplesGiven, opts.Samples)
	r.sampleIndex = make(map[string]int, len(r.samples))
	for i, name := range r.samples {
		r.sampleIndex[name] = i
	}
	r.pscNRefHom = make([]int, len(r.samples))
	r.pscNNonRefHom = make([]int, len(r.samples))
	r.pscNHets = make([]int, len(r.samples))
	r.pscNTs = make([]int, len(r.samples))
	r.pscNTv = make([]int, len(r.samples))
	r.pscNIndels = make([]int, len(r.samples))
	r.pscDepthSum = make([]int, len(r.samples))
	r.pscDepthN = make([]int, len(r.samples))
	r.pscNSingleton = make([]int, len(r.samples))
	r.pscNHapRef = make([]int, len(r.samples))
	r.pscNHapAlt = make([]int, len(r.samples))
	r.pscNMissing = make([]int, len(r.samples))
	r.psiInsHets = make([]int, len(r.samples))
	r.psiDelHets = make([]int, len(r.samples))
	r.psiInsHoms = make([]int, len(r.samples))
	r.psiDelHoms = make([]int, len(r.samples))
	return r
}

// filterSampleSet returns the sample set used for per-sample stats. Upstream
// only populates per-sample counters (and gates the PSC/PSI/HWE sections and
// the DP genotype histogram) when `-s`/`-S` was given — i.e. when
// `args->files->n_smpl > 0`. When given is false we return an empty set so
// those sections are suppressed, matching `bcftools stats` without -s.
//
// When given is true, `-s -` (or an empty/`"-"` list) selects every header
// sample in header order; otherwise we return the intersection of
// headerSamples and want, preserving the requested order.
func filterSampleSet(headerSamples []string, given bool, want []string) []string {
	if !given {
		return nil
	}
	allSamples := len(want) == 0
	if len(want) == 1 && want[0] == "-" {
		allSamples = true
	}
	if allSamples {
		out := make([]string, len(headerSamples))
		copy(out, headerSamples)
		return out
	}
	seen := make(map[string]bool, len(headerSamples))
	for _, s := range headerSamples {
		seen[s] = true
	}
	var out []string
	for _, w := range want {
		if seen[w] {
			out = append(out, w)
		}
	}
	return out
}

// Stats consumes a VCF/BCF stream and writes the upstream-style tab-prefixed
// stats report to out. Returns the underlying accumulator for callers that
// need to assert structural values without re-parsing text output.
func Stats(in io.Reader, out io.Writer, opts StatsOptions) (*statsResult, error) {
	br := bufio.NewReader(in)
	head, err := br.Peek(5)
	if err != nil && err != io.EOF {
		return nil, err
	}
	if len(head) >= 5 && head[0] == 'B' && head[1] == 'C' && head[2] == 'F' {
		return statsFromBCF(br, out, opts)
	}
	return statsFromVCF(br, out, opts)
}

// StatsFile opens path through iohelper and emits stats.
func StatsFile(path string, out io.Writer, opts StatsOptions) (*statsResult, error) {
	in, err := iohelper.OpenReader(path)
	if err != nil {
		return nil, err
	}
	defer in.Close()
	if opts.InputFile == "" {
		opts.InputFile = path
	}
	return Stats(in, out, opts)
}

// statsFromVCF runs the streaming counter over a VCF text stream.
func statsFromVCF(in io.Reader, out io.Writer, opts StatsOptions) (*statsResult, error) {
	r := vcf.NewReader(in)
	hdr, err := r.ReadHeader()
	if err != nil {
		return nil, err
	}
	res := newStatsResult(opts, hdr.Samples)
	includeF, excludeF, err := compileStatsExpressions(opts)
	if err != nil {
		return nil, err
	}
	targets, err := parseStatsTargets(opts)
	if err != nil {
		return nil, err
	}
	for {
		v, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, err
		}
		if !keepStatsVariant(v, opts, includeF, excludeF, targets) {
			continue
		}
		accumulate(res, v)
	}
	if err := writeStats(out, res); err != nil {
		return res, err
	}
	return res, nil
}

// statsFromBCF runs the streaming counter over a BCF stream.
func statsFromBCF(in io.Reader, out io.Writer, opts StatsOptions) (*statsResult, error) {
	br, err := bcf.NewReader(in)
	if err != nil {
		return nil, err
	}
	hdr := br.Header().VCF
	res := newStatsResult(opts, hdr.Samples)
	includeF, excludeF, err := compileStatsExpressions(opts)
	if err != nil {
		return nil, err
	}
	targets, err := parseStatsTargets(opts)
	if err != nil {
		return nil, err
	}
	for {
		rec, err := br.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, err
		}
		v := rec.ToVariant(br.Header())
		if !keepStatsVariant(v, opts, includeF, excludeF, targets) {
			continue
		}
		accumulate(res, v)
	}
	if err := writeStats(out, res); err != nil {
		return res, err
	}
	return res, nil
}

// compileStatsExpressions reuses the view-side expression compiler.
func compileStatsExpressions(opts StatsOptions) (include, exclude *Filter, err error) {
	if opts.IncludeExpr != "" {
		include, err = CompileFilter(opts.IncludeExpr)
		if err != nil {
			return nil, nil, err
		}
	}
	if opts.ExcludeExpr != "" {
		exclude, err = CompileFilter(opts.ExcludeExpr)
		if err != nil {
			return nil, nil, err
		}
	}
	return include, exclude, nil
}

// parseStatsTargets merges -r/-R/-t/-T into a single region slice the
// streaming loop can apply as a post-filter. Without an index we always
// degrade to a post-filter — mirroring the bcftools view fallback.
func parseStatsTargets(opts StatsOptions) ([]region, error) {
	var specs []string
	specs = append(specs, opts.Regions...)
	specs = append(specs, opts.Targets...)
	if len(specs) == 0 {
		return nil, nil
	}
	return parseRegions(specs)
}

// keepStatsVariant runs every filter in opts. Returns true when the variant
// contributes to the counters.
func keepStatsVariant(v *vcf.Variant, opts StatsOptions, includeF, excludeF *Filter, targets []region) bool {
	if len(targets) > 0 && !overlapsAny(v, targets) {
		return false
	}
	if len(opts.ApplyFilters) > 0 {
		ok := false
		for _, want := range opts.ApplyFilters {
			for _, f := range v.Filter {
				if f == want {
					ok = true
					break
				}
			}
		}
		if !ok {
			return false
		}
	}
	if includeF != nil && !includeF.Eval(v) {
		return false
	}
	if excludeF != nil && excludeF.Eval(v) {
		return false
	}
	return true
}

// accumulate folds one variant into the result. Multi-allelic sites are
// counted once for the SN section; individual ALTs each contribute to AF,
// QUAL, ST, IDD bookkeeping unless --1st-allele-only restricts us.
func accumulate(r *statsResult, v *vcf.Variant) {
	r.numRecords++

	// Detect no-ALT records (gVCF-style spanning blocks, "." or no ALT).
	hasALT := false
	for _, a := range v.Alt {
		if a != "" && a != "." {
			hasALT = true
			break
		}
	}
	if !hasALT {
		r.numNoALTs++
		return
	}

	alts := v.Alt
	if r.opts.FirstAlleleOnly && len(alts) > 1 {
		alts = alts[:1]
	}
	if len(v.Alt) > 1 {
		r.numMA++
	}

	// SN / per-allele classification.
	sawSNP := false
	sawNonSNP := false
	for _, alt := range alts {
		if alt == "" || alt == "." || alt == "*" {
			continue
		}
		kind := classifyVariant(v.Ref, alt)
		switch kind {
		case "snp":
			sawSNP = true
		case "mnp":
			sawNonSNP = true
		case "indel":
			sawNonSNP = true
		default:
			sawNonSNP = true
		}
	}
	if sawSNP {
		r.numSNPs++
	}
	if len(v.Alt) > 1 && countSNPAlts(v) > 1 {
		r.numMASNP++
	}
	// Classify into MNP / indel / other for SN counters (per-record).
	if !sawSNP && sawNonSNP {
		nonSNPKind := dominantKind(v.Ref, alts)
		switch nonSNPKind {
		case "mnp":
			r.numMNPs++
		case "indel":
			r.numIndels++
		default:
			r.numOthers++
		}
	} else if sawSNP && sawNonSNP {
		// Mixed sites count as indel in upstream when an indel is present,
		// else MNP — but for SN-numbers upstream just bumps SNP and the
		// "indel" counter independently. We approximate that here.
		nonSNPKind := dominantKind(v.Ref, alts)
		switch nonSNPKind {
		case "indel":
			r.numIndels++
		case "mnp":
			r.numMNPs++
		default:
			r.numOthers++
		}
	}

	// Legacy QUAL binning (floor(qual)) feeds the internal qualSNPs map
	// used by `TestStatsQUALBinning`. The upstream `iqual` scheme below
	// feeds the textual QUAL output.
	qualBin := int(math.Floor(v.Qual))
	if v.Qual < 0 {
		qualBin = 0
	}
	iqual := 0
	if !math.IsNaN(v.Qual) && v.Qual >= 0 {
		iqual = 1 + int(v.Qual*10)
	}

	// Legacy AF binning (TestStatsAFBinning / TestStatsAFTag).
	af := computeAF(r, v)
	afIdx := afBinIndex(r.afBins, af)

	// Upstream per-ALT AF bin assignment (vcfstats.c:643).
	tmpUSBin := computeUpstreamAFBins(r, v)
	hweUSBin := 0
	if len(tmpUSBin) > 1 {
		hweUSBin = tmpUSBin[1]
	}

	for i, alt := range alts {
		if alt == "" || alt == "." || alt == "*" {
			continue
		}
		usBin := 0
		if i+1 < len(tmpUSBin) {
			usBin = tmpUSBin[i+1]
		}
		switch classifyVariant(v.Ref, alt) {
		case "snp":
			r.afSNPs[afIdx]++
			r.qualSNPs[qualBin]++
			r.afSnpsUS[usBin]++
			tsType := transitionType(v.Ref, alt)
			if tsType == "ts" {
				r.afTs[afIdx]++
				r.qualTs[qualBin]++
				r.afTsUS[usBin]++
				if i == 0 {
					r.tsAlt1++
					r.qualTsUS[iqual]++
				}
			} else {
				r.afTv[afIdx]++
				r.qualTv[qualBin]++
				r.afTvUS[usBin]++
				if i == 0 {
					r.tvAlt1++
					r.qualTvUS[iqual]++
				}
			}
			key := strings.ToUpper(v.Ref) + ">" + strings.ToUpper(alt)
			r.subst[key]++
		case "indel":
			r.afNonS[afIdx]++
			r.qualNonS[qualBin]++
			r.afIndUS[usBin]++
			r.afRepNAUS[usBin]++ // upstream's `af_repeats[2]++` when no indel_ctx (vcfstats.c:767-768)
			r.qualIndS[iqual]++
			delta := len(alt) - len(v.Ref)
			r.indelLen[delta]++
		case "mnp":
			r.afNonS[afIdx]++
			r.qualNonS[qualBin]++
		default:
			r.afNonS[afIdx]++
			r.qualNonS[qualBin]++
		}
	}

	// DP binning — INFO/DP first, falling back to per-sample DP sum.
	if dp, ok := siteDepth(v); ok {
		bin := dpBin(r.opts, dp)
		r.dpSites[bin]++
		r.dpTotalSite++
	}

	// Per-sample counters (PSC/PSI/HWE).
	nRefHom, nHet, nAlt := accumulateSamples(r, v, alts)

	// Legacy HWE chi-square accumulator (TestStatsHWEChiSquare).
	if sawSNP {
		chi, ok := hweChiSquare(v)
		if ok {
			bucket := int(math.Round(af * 1000))
			r.hweObs[bucket]++
			r.hweChiSum[bucket] += chi
		}
	}

	// Upstream HWE accumulator (vcfstats.c:1158-1173).
	total := nRefHom + nHet + nAlt
	if len(v.Alt) > 0 && total > 0 && hweUSBin >= 0 && hweUSBin < mAFBins {
		hetFrac := float64(nHet) / float64(total)
		iHet := int(hetFrac * float64(naFHWE-1))
		if iHet >= naFHWE {
			iHet = naFHWE - 1
		}
		if iHet < 0 {
			iHet = 0
		}
		r.afHWE[hweUSBin][iHet]++
	}
}

// computeUpstreamAFBins returns the upstream `tmp_iaf` per-allele AF bin
// assignment (vcfstats.c:643). The result has `len(v.Alt)+1` entries —
// index 0 reserved for REF (always 0), indexes 1..N for v.Alt[0..N-1].
//
// Bin 0 is the singleton (AC==1) bucket; bins 1..mAFBins-1 hold
// non-singletons via `floor(af*(mAFBins-2))+1`. When --af-tag is set the
// value comes straight from INFO; otherwise AC/AN are computed from
// FORMAT/GT first, then fall back to INFO/AC,AN. When no information is
// available every allele falls into bin 0 — matching upstream.
func computeUpstreamAFBins(r *statsResult, v *vcf.Variant) []int {
	out := make([]int, len(v.Alt)+1)
	if r.opts.AFTag != "" {
		raw, ok := v.Info[r.opts.AFTag]
		if !ok || raw == "" {
			return out
		}
		parts := strings.Split(raw, ",")
		if len(parts) != len(v.Alt) {
			return out
		}
		for i, p := range parts {
			f, err := strconv.ParseFloat(p, 32)
			if err != nil {
				continue
			}
			af := float32(f)
			if af < 0 {
				af = 0
			} else if af > 1 {
				af = 1
			}
			iaf := int(af * float32(mAFBins-2))
			if iaf >= mAFBins-1 {
				iaf = mAFBins - 2
			}
			out[i+1] = iaf + 1
		}
		return out
	}
	// Upstream `bcf_calc_ac` prefers INFO/AC,AN when both are present
	// (htslib vcfutils.c:39-92) and only falls back to FORMAT/GT when
	// they are not (vcfutils.c:94-131). Match that order here so the
	// AF bins agree with upstream byte-for-byte.
	ac, an := infoACAN(v)
	if an == 0 {
		ac, an = computeAlleleCountsAll(v)
	}
	if an == 0 {
		return out
	}
	for i := 1; i < len(out); i++ {
		count := 0
		if i < len(ac) {
			count = ac[i]
		}
		if count == 1 {
			out[i] = 0
			continue
		}
		// Upstream computes `af = (float)ac/an` (32-bit) and then
		// `iaf = af*(m_af-2)` truncated to int (vcfstats.c:692-695).
		// We mirror the same 32-bit precision so the truncation lands
		// on the same bin even at boundary fractions like 2/6.
		af := float32(count) / float32(an)
		if af < 0 {
			af = 0
		} else if af > 1 {
			af = 1
		}
		iaf := int(af * float32(mAFBins-2))
		if iaf >= mAFBins-1 {
			iaf = mAFBins - 2
		}
		out[i] = iaf + 1
	}
	return out
}

// computeAlleleCountsAll counts allele observations across FORMAT/GT,
// returning per-allele counts (REF at index 0, each ALT at 1..N) and the
// total AN. Returns (nil, 0) when no usable GT field is present.
func computeAlleleCountsAll(v *vcf.Variant) ([]int, int) {
	if len(v.Samples) == 0 {
		return nil, 0
	}
	ac := make([]int, len(v.Alt)+1)
	an := 0
	for _, s := range v.Samples {
		gt, ok := s.Data["GT"]
		if !ok || gt == "" || gt == "." {
			continue
		}
		for _, a := range splitGTAlleles(gt) {
			if a == "." {
				continue
			}
			n, err := strconv.Atoi(a)
			if err != nil || n < 0 || n >= len(ac) {
				continue
			}
			ac[n]++
			an++
		}
	}
	return ac, an
}

// splitGTAlleles splits a GT string ("0/1", "1|2|3") into per-allele tokens.
func splitGTAlleles(gt string) []string {
	return strings.FieldsFunc(gt, func(r rune) bool {
		return r == '/' || r == '|'
	})
}

// infoACAN reconstructs per-ALT AC and AN from INFO/AC and INFO/AN.
// When AN is known the REF count is filled in as `AN - sum(ALT AC)`.
func infoACAN(v *vcf.Variant) ([]int, int) {
	an := 0
	if raw, ok := v.Info["AN"]; ok && raw != "" {
		if n, err := strconv.Atoi(raw); err == nil {
			an = n
		}
	}
	ac := make([]int, len(v.Alt)+1)
	if raw, ok := v.Info["AC"]; ok && raw != "" {
		parts := strings.Split(raw, ",")
		for i, p := range parts {
			if i+1 >= len(ac) {
				break
			}
			if n, err := strconv.Atoi(p); err == nil {
				ac[i+1] = n
			}
		}
	}
	if an > 0 {
		ref := an
		for i := 1; i < len(ac); i++ {
			ref -= ac[i]
		}
		if ref < 0 {
			ref = 0
		}
		ac[0] = ref
	}
	return ac, an
}

// dominantKind returns whether a non-SNP record is MNP, indel, or "other".
func dominantKind(ref string, alts []string) string {
	indel := false
	mnp := false
	for _, alt := range alts {
		if alt == "" || alt == "." || alt == "*" {
			continue
		}
		switch classifyVariant(ref, alt) {
		case "indel":
			indel = true
		case "mnp":
			mnp = true
		}
	}
	switch {
	case indel:
		return "indel"
	case mnp:
		return "mnp"
	default:
		return "other"
	}
}

// countSNPAlts returns how many of v.Alt are SNPs (used for "multi-allelic
// SNP sites" in SN).
func countSNPAlts(v *vcf.Variant) int {
	n := 0
	for _, alt := range v.Alt {
		if classifyVariant(v.Ref, alt) == "snp" {
			n++
		}
	}
	return n
}

// classifyVariant returns one of "snp", "mnp", "indel", "other" for a single
// REF/ALT pair. We use the simple, structural rule: equal-length single-base
// = SNP; equal-length multi-base = MNP; different lengths = indel; otherwise
// other (symbolic alleles like <DEL>).
func classifyVariant(ref, alt string) string {
	if alt == "" || alt == "." || alt == "*" {
		return "other"
	}
	if strings.HasPrefix(alt, "<") || strings.ContainsAny(alt, "[]") {
		return "other"
	}
	if len(ref) == 1 && len(alt) == 1 {
		if isACGT(ref[0]) && isACGT(alt[0]) {
			return "snp"
		}
		return "other"
	}
	if len(ref) == len(alt) {
		return "mnp"
	}
	return "indel"
}

// isACGT reports whether b is one of A/C/G/T (case-insensitive).
func isACGT(b byte) bool {
	switch b {
	case 'A', 'C', 'G', 'T', 'a', 'c', 'g', 't':
		return true
	}
	return false
}

// transitionType returns "ts" (A<->G or C<->T) or "tv" for any other SNP.
func transitionType(ref, alt string) string {
	r := strings.ToUpper(ref)[0]
	a := strings.ToUpper(alt)[0]
	if (r == 'A' && a == 'G') || (r == 'G' && a == 'A') ||
		(r == 'C' && a == 'T') || (r == 'T' && a == 'C') {
		return "ts"
	}
	return "tv"
}

// computeAF returns the allele frequency of the first ALT allele. When
// opts.AFTag is set we read it from INFO; otherwise we compute it from
// genotypes. Returns 0.0 when no information is available.
func computeAF(r *statsResult, v *vcf.Variant) float64 {
	if r.opts.AFTag != "" {
		if raw, ok := v.Info[r.opts.AFTag]; ok && raw != "" {
			first := raw
			if i := strings.Index(raw, ","); i >= 0 {
				first = raw[:i]
			}
			if f, err := strconv.ParseFloat(first, 64); err == nil {
				return f
			}
		}
	}
	if raw, ok := v.Info["AF"]; ok && raw != "" {
		first := raw
		if i := strings.Index(raw, ","); i >= 0 {
			first = raw[:i]
		}
		if f, err := strconv.ParseFloat(first, 64); err == nil {
			return f
		}
	}
	ac, an := computeAC(v)
	if an == 0 {
		return 0
	}
	return float64(ac) / float64(an)
}

// afBinIndex returns the bin position for an allele frequency value.
func afBinIndex(bins []float64, af float64) int {
	if len(bins) < 2 {
		return 0
	}
	if af <= bins[0] {
		return 0
	}
	if af >= bins[len(bins)-1] {
		return len(bins) - 2
	}
	for i := 0; i < len(bins)-1; i++ {
		if af >= bins[i] && af < bins[i+1] {
			return i
		}
	}
	return len(bins) - 2
}

// siteDepth returns INFO/DP if present. Upstream's `dp_sites` histogram is
// driven solely by a scalar INFO/DP value (vcfstats.c:1307-1308 —
// `bcf_get_info_int32(...,"DP",...)==1`); when INFO/DP is absent the site does
// not contribute to the DP section at all. There is deliberately no FORMAT/DP
// fallback here: per-sample FORMAT/DP feeds the separate genotype histogram.
func siteDepth(v *vcf.Variant) (int, bool) {
	raw, ok := v.Info["DP"]
	if !ok || raw == "" || raw == "." {
		return 0, false
	}
	// Only a single scalar value counts (upstream requires the return of
	// bcf_get_info_int32 to be exactly 1).
	if strings.Contains(raw, ",") {
		return 0, false
	}
	n, err := strconv.Atoi(raw)
	if err != nil {
		return 0, false
	}
	return n, true
}

// dpBin returns the depth bin for value dp using opts.DepthMin/Max/Step. With
// defaults (0, 500, 1) each integer depth ends up in its own bucket — matching
// bcftools stats.
func dpBin(opts StatsOptions, dp int) int {
	min, max, step := opts.DepthMin, opts.DepthMax, opts.DepthStep
	if step <= 0 {
		step = 1
	}
	if max <= 0 {
		max = 500
	}
	if dp < min {
		return min
	}
	if dp > max {
		return max
	}
	return min + ((dp-min)/step)*step
}

// accumulateSamples implements upstream `do_sample_stats` (vcfstats.c:1094).
// Returns the (nRefHom, nHet, nAlt) tally for this record so the caller
// can populate the upstream HWE 2D accumulator.
//
// Key semantics that diverge from a naïve port and that we preserve for
// byte-for-byte parity:
//
//   - Missing genotypes (./.) increment pscNMissing regardless of the
//     site's variant type.
//   - Haploid REF/ALT GTs increment pscNHapRef / pscNHapAlt.
//   - The diploid SNP-only counters (homRR/homAA/hets/ts/tv) only count
//     SNP-typed alleles. An indel-only site (1/1 ATG>A) does NOT bump
//     pscNNonRefHom; pscNIndels covers it instead.
//   - For a 1/2 SNP het upstream scores ts/tv for both ial and jal.
//   - pscNIndels is bumped once per non-REF diploid GT at any INDEL site.
//   - pscNSingleton is bumped on the unique non-REF sample when a site
//     has exactly one non-REF diploid GT (vcfstats.c:1156).
func accumulateSamples(r *statsResult, v *vcf.Variant, alts []string) (nRefHom, nHet, nAlt int) {
	if len(r.samples) == 0 {
		return 0, 0, 0
	}
	altKinds := make([]string, len(alts))
	hasINDELALT := false
	for i, alt := range alts {
		altKinds[i] = classifyVariant(v.Ref, alt)
		if altKinds[i] == "indel" {
			hasINDELALT = true
		}
	}
	nNonRef, iNonRef := 0, -1

	for _, s := range v.Samples {
		idx, ok := r.sampleIndex[s.Name]
		if !ok {
			continue
		}
		// Depth is recorded regardless of GT (upstream calls
		// calc_sample_depth before sample_gt_stats — vcfstats.c:1115).
		if raw, ok := s.Data["DP"]; ok && raw != "" && raw != "." {
			// Upstream only records depths > 0 (vcfstats.c:1116) — a
			// zero or negative/missing depth is skipped entirely.
			if n, err := strconv.Atoi(raw); err == nil && n > 0 {
				r.pscDepthSum[idx] += n
				r.pscDepthN[idx]++
				bin := dpBin(r.opts, n)
				r.dpGTs[bin]++
				r.dpTotalGT++
			}
		}
		gtStr, hasGT := s.Data["GT"]
		if !hasGT || gtStr == "" {
			continue
		}
		kind, ial, jal, _ := classifySampleGT(gtStr)
		switch kind {
		case gtUnknown:
			r.pscNMissing[idx]++
			continue
		case gtHaploidRef:
			r.pscNHapRef[idx]++
			continue
		case gtHaploidAlt:
			r.pscNHapAlt[idx]++
			continue
		}

		if kind != gtHomRefRef {
			nNonRef++
			iNonRef = idx
		}

		ialKind := altKindAt(altKinds, ial)
		jalKind := altKindAt(altKinds, jal)
		hasSNPAllele := (ial > 0 && ialKind == "snp") || (jal > 0 && jalKind == "snp")
		// Upstream tallies HWE per-record using GT_HOM_RR -> nref,
		// GT_HET_RA -> nhet, and {GT_HET_AA, GT_HOM_AA} -> nalt
		// (vcfstats.c:1022-1029). The PSC group separately tracks
		// homAA/hets, but the HWE accumulator uses these three buckets.
		switch kind {
		case gtHomRefRef:
			nRefHom++
		case gtHetRefAlt:
			nHet++
		case gtHetAltAlt, gtHomAltAlt:
			nAlt++
		}
		if ial == 0 && jal == 0 {
			r.pscNRefHom[idx]++
		} else if hasSNPAllele {
			switch kind {
			case gtHetRefAlt, gtHetAltAlt:
				r.pscNHets[idx]++
			case gtHomAltAlt:
				r.pscNNonRefHom[idx]++
			}
			if ial > 0 && ialKind == "snp" {
				if transitionType(v.Ref, alts[ial-1]) == "ts" {
					r.pscNTs[idx]++
				} else {
					r.pscNTv[idx]++
				}
			}
			if jal > 0 && jalKind == "snp" && jal != ial {
				if transitionType(v.Ref, alts[jal-1]) == "ts" {
					r.pscNTs[idx]++
				} else {
					r.pscNTv[idx]++
				}
			}
		}
		// Non-SNP diploid GTs already contributed to nRefHom/nHet/nAlt
		// above via the upstream HWE accumulator block.

		// Indel counters (smpl_indels + PSI ins/del). Upstream bumps
		// these when the line has an INDEL allele AND the GT is not
		// HOM_RR (vcfstats.c:1058-1085).
		if hasINDELALT && kind != gtHomRefRef {
			r.pscNIndels[idx]++
			switch kind {
			case gtHetRefAlt, gtHetAltAlt:
				// Upstream (vcfstats.c:1063-1078) sets is_ins/is_del from
				// BOTH ial and jal, then bumps del_hets and/or ins_hets.
				// An alt-het with one ins and one del allele is counted in
				// both, hence the separate ifs (not else-if).
				isIns, isDel := classifyIndelAllele(v.Ref, alts, ial)
				ji, jd := classifyIndelAllele(v.Ref, alts, jal)
				if ji {
					isIns = true
				}
				if jd {
					isDel = true
				}
				if isDel {
					r.psiDelHets[idx]++
				}
				if isIns {
					r.psiInsHets[idx]++
				}
			case gtHomAltAlt:
				// HOM_AA: classify the single ALT allele (ial) as ins/del
				// (vcfstats.c:1080-1083).
				if isIns, isDel := classifyIndelAllele(v.Ref, alts, ial); isIns || isDel {
					if isDel {
						r.psiDelHoms[idx]++
					} else {
						r.psiInsHoms[idx]++
					}
				}
			}
		}
	}
	if nNonRef == 1 && iNonRef >= 0 {
		r.pscNSingleton[iNonRef]++
	}
	return nRefHom, nHet, nAlt
}

// classifyIndelAllele reports whether the allele at the 1-based ALT
// index `ial` is an insertion or deletion relative to REF. Returns
// (false,false) for non-indel or out-of-range indices.
func classifyIndelAllele(ref string, alts []string, ial int) (isIns, isDel bool) {
	if ial <= 0 || ial > len(alts) {
		return false, false
	}
	a := alts[ial-1]
	if classifyVariant(ref, a) != "indel" {
		return false, false
	}
	if len(a) > len(ref) {
		return true, false
	}
	return false, true
}

// altKindAt returns the variant kind for line allele index `idx` where
// 0 = REF (returned as the empty string) and 1..N indexes alts[0..N-1].
func altKindAt(altKinds []string, idx int) string {
	if idx <= 0 || idx > len(altKinds) {
		return ""
	}
	return altKinds[idx-1]
}

// Genotype-classification constants used by classifySampleGT.
const (
	gtUnknown    = iota // ./. or invalid
	gtHomRefRef         // 0/0 diploid
	gtHetRefAlt         // 0/N diploid
	gtHetAltAlt         // M/N diploid with both non-zero and M!=N
	gtHomAltAlt         // N/N diploid with N>0
	gtHaploidRef        // 0 (haploid)
	gtHaploidAlt        // N (haploid, N>0)
)

// classifySampleGT mirrors htslib `bcf_gt_type` (vcfutils.c:135) on its
// diploid / haploid / missing branches and additionally returns the
// smaller (`ial`) and larger (`jal`) non-zero allele indices and the
// ploidy. Indices are 1-based; 0 means the REF allele.
func classifySampleGT(gt string) (kind, ial, jal, ploidy int) {
	if gt == "" || gt == "." {
		return gtUnknown, 0, 0, 0
	}
	parts := splitGTAlleles(gt)
	if len(parts) == 0 {
		return gtUnknown, 0, 0, 0
	}
	hasRef, hasAlt := false, false
	for _, p := range parts {
		if p == "." {
			return gtUnknown, 0, 0, 0
		}
		n, err := strconv.Atoi(p)
		if err != nil || n < 0 {
			return gtUnknown, 0, 0, 0
		}
		if n == 0 {
			hasRef = true
			continue
		}
		hasAlt = true
		switch {
		case ial == 0:
			ial = n
		case n < ial:
			if jal == 0 || ial > jal {
				jal = ial
			}
			ial = n
		case n != ial:
			if jal == 0 || n > jal {
				jal = n
			}
		}
	}
	ploidy = len(parts)
	if ploidy == 1 {
		if hasRef {
			return gtHaploidRef, 0, 0, 1
		}
		return gtHaploidAlt, ial, 0, 1
	}
	if !hasAlt {
		return gtHomRefRef, 0, 0, ploidy
	}
	if !hasRef {
		if jal == 0 || jal == ial {
			return gtHomAltAlt, ial, ial, ploidy
		}
		return gtHetAltAlt, ial, jal, ploidy
	}
	return gtHetRefAlt, ial, 0, ploidy
}

// parseDiploidGT splits a "0/1" / "0|1" string into the integer allele
// indices. Retained for the legacy `TestStatsParseDiploidGT` unit test
// and the `hweChiSquare` helper.
func parseDiploidGT(gt string) (a, b int8, sep byte, ok bool) {
	if gt == "" || gt == "." {
		return 0, 0, '/', false
	}
	sepIdx := strings.IndexAny(gt, "/|")
	if sepIdx < 0 {
		// Haploid: treat as a/a.
		n, err := strconv.Atoi(gt)
		if err != nil {
			return 0, 0, '/', false
		}
		return int8(n), int8(n), '/', true
	}
	left := gt[:sepIdx]
	right := gt[sepIdx+1:]
	if left == "." || right == "." {
		return 0, 0, '/', false
	}
	la, err := strconv.Atoi(left)
	if err != nil {
		return 0, 0, '/', false
	}
	lb, err := strconv.Atoi(right)
	if err != nil {
		return 0, 0, '/', false
	}
	return int8(la), int8(lb), gt[sepIdx], true
}

// hweChiSquare returns Pearson's chi-square statistic for the genotypes at a
// SNP site. Returns false if there's not enough data.
func hweChiSquare(v *vcf.Variant) (float64, bool) {
	var nAA, nAa, naa int
	for _, s := range v.Samples {
		gt, ok := s.Data["GT"]
		if !ok {
			continue
		}
		a, b, _, parsed := parseDiploidGT(gt)
		if !parsed {
			continue
		}
		switch {
		case a == 0 && b == 0:
			nAA++
		case a != b:
			nAa++
		default:
			naa++
		}
	}
	n := nAA + nAa + naa
	if n == 0 {
		return 0, false
	}
	p := float64(2*nAA+nAa) / float64(2*n)
	q := 1.0 - p
	expAA := float64(n) * p * p
	expAa := 2.0 * float64(n) * p * q
	expaa := float64(n) * q * q
	chi := 0.0
	for _, pair := range [3][2]float64{{float64(nAA), expAA}, {float64(nAa), expAa}, {float64(naa), expaa}} {
		obs, exp := pair[0], pair[1]
		if exp <= 0 {
			continue
		}
		d := obs - exp
		chi += d * d / exp
	}
	return chi, true
}

// writeStats emits the upstream-style tab-prefixed report.
func writeStats(out io.Writer, r *statsResult) error {
	bw := bufio.NewWriter(out)

	// File header. The provenance lines below carry the version/command
	// and are stripped by the oracle comparison; the rest is byte-parity.
	fmt.Fprintln(bw, "# This file was produced by bcftools stats (pure-Go).")
	fmt.Fprintln(bw, "# The command line was:\tbcftools stats", r.opts.InputFile)
	fmt.Fprintln(bw, "#")

	// ID block — upstream lists each input file's ID; we only use one (0).
	fmt.Fprintln(bw, "# Definition of sets:")
	fmt.Fprintln(bw, "# ID\t[2]id\t[3]tab-separated file names")
	fmt.Fprintf(bw, "ID\t0\t%s\n", r.opts.InputFile)

	if err := writeSN(bw, r); err != nil {
		return err
	}
	if err := writeTSTV(bw, r); err != nil {
		return err
	}
	if err := writeSiS(bw, r); err != nil {
		return err
	}
	if err := writeAF(bw, r); err != nil {
		return err
	}
	if err := writeQUAL(bw, r); err != nil {
		return err
	}
	if err := writeIDD(bw, r); err != nil {
		return err
	}
	if err := writeST(bw, r); err != nil {
		return err
	}
	if err := writeDP(bw, r); err != nil {
		return err
	}
	// Upstream only emits PSC/PSI/HWE when -s/-S was given
	// (args->files->n_smpl>0). Without samples these sections are absent.
	if r.opts.SamplesGiven {
		if err := writePSC(bw, r); err != nil {
			return err
		}
		if err := writePSI(bw, r); err != nil {
			return err
		}
		if err := writeHWE(bw, r); err != nil {
			return err
		}
	}
	return bw.Flush()
}

func writeSN(bw *bufio.Writer, r *statsResult) error {
	fmt.Fprintln(bw, "# SN, Summary numbers:")
	fmt.Fprintln(bw, "#   number of records   .. number of data rows in the VCF")
	fmt.Fprintln(bw, "#   number of no-ALTs   .. reference-only sites, ALT is either \".\" or identical to REF")
	fmt.Fprintln(bw, "#   number of SNPs      .. number of rows with a SNP")
	fmt.Fprintln(bw, "#   number of MNPs      .. number of rows with a MNP, such as CC>TT")
	fmt.Fprintln(bw, "#   number of indels    .. number of rows with an indel")
	fmt.Fprintln(bw, "#   number of others    .. number of rows with other type, for example a symbolic allele or")
	fmt.Fprintln(bw, "#                          a complex substitution, such as ACT>TCGA")
	fmt.Fprintln(bw, "#   number of multiallelic sites     .. number of rows with multiple alternate alleles")
	fmt.Fprintln(bw, "#   number of multiallelic SNP sites .. number of rows with multiple alternate alleles, all SNPs")
	fmt.Fprintln(bw, "# ")
	fmt.Fprintln(bw, "#   Note that rows containing multiple types will be counted multiple times, in each")
	fmt.Fprintln(bw, "#   counter. For example, a row with a SNP and an indel increments both the SNP and")
	fmt.Fprintln(bw, "#   the indel counter.")
	fmt.Fprintln(bw, "# ")
	fmt.Fprintln(bw, "# SN\t[2]id\t[3]key\t[4]value")
	fmt.Fprintf(bw, "SN\t0\tnumber of samples:\t%d\n", len(r.headerSamples))
	fmt.Fprintf(bw, "SN\t0\tnumber of records:\t%d\n", r.numRecords)
	fmt.Fprintf(bw, "SN\t0\tnumber of no-ALTs:\t%d\n", r.numNoALTs)
	fmt.Fprintf(bw, "SN\t0\tnumber of SNPs:\t%d\n", r.numSNPs)
	fmt.Fprintf(bw, "SN\t0\tnumber of MNPs:\t%d\n", r.numMNPs)
	fmt.Fprintf(bw, "SN\t0\tnumber of indels:\t%d\n", r.numIndels)
	fmt.Fprintf(bw, "SN\t0\tnumber of others:\t%d\n", r.numOthers)
	fmt.Fprintf(bw, "SN\t0\tnumber of multiallelic sites:\t%d\n", r.numMA)
	fmt.Fprintf(bw, "SN\t0\tnumber of multiallelic SNP sites:\t%d\n", r.numMASNP)
	return nil
}

// writeTSTV emits the transitions/transversions section
// (vcfstats.c:1389-1399). Columns 3/4 are the total ts/tv summed over all
// AF bins; columns 6/7 are the 1st-ALT-only ts/tv. The ratios use %.2f and
// are 0 when the denominator is zero.
func writeTSTV(bw *bufio.Writer, r *statsResult) error {
	fmt.Fprintln(bw, "# TSTV, transitions/transversions")
	fmt.Fprintln(bw, "#   - transitions, see https://en.wikipedia.org/wiki/Transition_(genetics)")
	fmt.Fprintln(bw, "#   - transversions, see https://en.wikipedia.org/wiki/Transversion")
	fmt.Fprintln(bw, "# TSTV\t[2]id\t[3]ts\t[4]tv\t[5]ts/tv\t[6]ts (1st ALT)\t[7]tv (1st ALT)\t[8]ts/tv (1st ALT)")
	ts, tv := 0, 0
	for i := 0; i < mAFBins; i++ {
		ts += r.afTsUS[i]
		tv += r.afTvUS[i]
	}
	// Upstream computes the ratios in 32-bit float (vcfstats.c:1398).
	var ratio float32
	if tv != 0 {
		ratio = float32(ts) / float32(tv)
	}
	var ratio1 float32
	if r.tvAlt1 != 0 {
		ratio1 = float32(r.tsAlt1) / float32(r.tvAlt1)
	}
	fmt.Fprintf(bw, "TSTV\t0\t%d\t%d\t%.2f\t%d\t%d\t%.2f\n",
		ts, tv, ratio, r.tsAlt1, r.tvAlt1, ratio1)
	return nil
}

// writeSiS emits the singleton stats section (vcfstats.c:1439-1449). It uses
// the AF bin-0 (AC==1) accumulators, which is why it must run before writeAF
// folds bin 0 into bin 1. Columns 8/9 (repeat-consistent/inconsistent) come
// from the unimplemented --indel-context path and are always 0; column 10
// ("not applicable") is the indel count in bin 0.
func writeSiS(bw *bufio.Writer, r *statsResult) error {
	fmt.Fprintln(bw, "# SiS, Singleton stats:")
	fmt.Fprintln(bw, "#   - allele count, i.e. the number of singleton genotypes (AC=1)")
	fmt.Fprintln(bw, "#   - number of transitions, see above")
	fmt.Fprintln(bw, "#   - number of transversions, see above")
	fmt.Fprintln(bw, "#   - repeat-consistent, inconsistent and n/a: experimental and useless stats [DEPRECATED]")
	fmt.Fprintln(bw, "# SiS\t[2]id\t[3]allele count\t[4]number of SNPs\t[5]number of transitions\t[6]number of transversions\t[7]number of indels\t[8]repeat-consistent\t[9]repeat-inconsistent\t[10]not applicable")
	indels := r.afRepNAUS[0]
	fmt.Fprintf(bw, "SiS\t0\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
		1, r.afSnpsUS[0], r.afTsUS[0], r.afTvUS[0],
		indels, 0, 0, r.afRepNAUS[0])
	return nil
}

// writeAF emits the AF section in upstream format. Singletons (bin 0) are
// folded into bin 1 at print time (vcfstats.c:1451-1452 / 1826).
func writeAF(bw *bufio.Writer, r *statsResult) error {
	fmt.Fprintln(bw, "# AF, Stats by non-reference allele frequency:")
	fmt.Fprintln(bw, "# AF\t[2]id\t[3]allele frequency\t[4]number of SNPs\t[5]number of transitions\t[6]number of transversions\t[7]number of indels\t[8]repeat-consistent\t[9]repeat-inconsistent\t[10]not applicable")
	snps := r.afSnpsUS
	ts := r.afTsUS
	tv := r.afTvUS
	ind := r.afIndUS
	repNA := r.afRepNAUS
	snps[1] += snps[0]
	ts[1] += ts[0]
	tv[1] += tv[0]
	ind[1] += ind[0]
	repNA[1] += repNA[0]
	for i := 1; i < mAFBins; i++ {
		if snps[i]+ts[i]+tv[i]+ind[i]+repNA[i] == 0 {
			continue
		}
		af := float64(i-1) / float64(mAFBins-1)
		// Columns 8/9 are af_repeats[0]/[1] (repeat-consistent /
		// inconsistent). They are only populated by the
		// `--indel-context` path, which we do not implement, so we
		// always emit 0. Column 10 is af_repeats[2] ("not applicable"),
		// which upstream bumps unconditionally for every indel allele
		// when no indel_ctx is set (vcfstats.c:767-768).
		fmt.Fprintf(bw, "AF\t0\t%f\t%d\t%d\t%d\t%d\t0\t0\t%d\n",
			af, snps[i], ts[i], tv[i], ind[i], repNA[i])
	}
	return nil
}

// writeQUAL emits the QUAL section in upstream format
// (vcfstats.c:1489-1526). Bins are keyed by `iqual = 1 + int(qual*10)`
// with key 0 reserved for missing/negative QUAL ("."). The printed
// quality is `0.1 * (key - 1)` formatted as `%.1f`.
func writeQUAL(bw *bufio.Writer, r *statsResult) error {
	fmt.Fprintln(bw, "# QUAL, Stats by quality")
	fmt.Fprintln(bw, "# QUAL\t[2]id\t[3]Quality\t[4]number of SNPs\t[5]number of transitions (1st ALT)\t[6]number of transversions (1st ALT)\t[7]number of indels")
	keys := unionMapKeys(r.qualTsUS, r.qualTvUS, r.qualIndS)
	sort.Ints(keys)
	for _, k := range keys {
		nts := r.qualTsUS[k]
		ntv := r.qualTvUS[k]
		nin := r.qualIndS[k]
		if nts+ntv+nin == 0 {
			continue
		}
		fmt.Fprint(bw, "QUAL\t0\t")
		if k <= 0 {
			fmt.Fprint(bw, ".")
		} else {
			fmt.Fprintf(bw, "%.1f", 0.1*float64(k-1))
		}
		fmt.Fprintf(bw, "\t%d\t%d\t%d\t%d\n", nts+ntv, nts, ntv, nin)
	}
	return nil
}

// writeIDD emits the indel-length distribution. Upstream prints `0\t.`
// for the per-genotype-VAF columns when there is no VAF data
// (vcfstats.c:1558-1559).
func writeIDD(bw *bufio.Writer, r *statsResult) error {
	fmt.Fprintln(bw, "# IDD, InDel distribution:")
	fmt.Fprintln(bw, "# IDD\t[2]id\t[3]length (deletions negative)\t[4]number of sites\t[5]number of genotypes\t[6]mean VAF")
	keys := make([]int, 0, len(r.indelLen))
	for k := range r.indelLen {
		keys = append(keys, k)
	}
	sort.Ints(keys)
	for _, k := range keys {
		fmt.Fprintf(bw, "IDD\t0\t%d\t%d\t0\t.\n", k, r.indelLen[k])
	}
	return nil
}

func writeST(bw *bufio.Writer, r *statsResult) error {
	fmt.Fprintln(bw, "# ST, Substitution types:")
	fmt.Fprintln(bw, "# ST\t[2]id\t[3]type\t[4]count")
	canonical := []string{
		"A>C", "A>G", "A>T",
		"C>A", "C>G", "C>T",
		"G>A", "G>C", "G>T",
		"T>A", "T>C", "T>G",
	}
	for _, k := range canonical {
		fmt.Fprintf(bw, "ST\t0\t%s\t%d\n", k, r.subst[k])
	}
	return nil
}

func writeDP(bw *bufio.Writer, r *statsResult) error {
	fmt.Fprintln(bw, "# DP, depth:")
	fmt.Fprintln(bw, "#   - set id, see above")
	fmt.Fprintln(bw, "#   - the depth bin, corresponds to the depth (unless --depth was given)")
	fmt.Fprintln(bw, "#   - number of genotypes with this depth (zero unless -s/-S was given)")
	fmt.Fprintln(bw, "#   - fraction of genotypes with this depth (zero unless -s/-S was given)")
	fmt.Fprintln(bw, "#   - number of sites with this depth")
	fmt.Fprintln(bw, "#   - fraction of sites with this depth")
	fmt.Fprintln(bw, "# DP, Depth distribution")
	fmt.Fprintln(bw, "# DP\t[2]id\t[3]bin\t[4]number of genotypes\t[5]fraction of genotypes (%)\t[6]number of sites\t[7]fraction of sites (%)")
	bins := unionMapKeys(r.dpSites, r.dpGTs)
	sort.Ints(bins)
	for _, b := range bins {
		var fracGT, fracSite float64
		if r.dpTotalGT > 0 {
			fracGT = 100.0 * float64(r.dpGTs[b]) / float64(r.dpTotalGT)
		}
		if r.dpTotalSite > 0 {
			fracSite = 100.0 * float64(r.dpSites[b]) / float64(r.dpTotalSite)
		}
		fmt.Fprintf(bw, "DP\t0\t%d\t%d\t%.6f\t%d\t%.6f\n",
			b, r.dpGTs[b], fracGT, r.dpSites[b], fracSite)
	}
	return nil
}

func writePSC(bw *bufio.Writer, r *statsResult) error {
	fmt.Fprintln(bw, "# PSC, Per-sample counts. Note that the ref/het/hom counts include only SNPs, for indels see PSI. The rest include both SNPs and indels.")
	fmt.Fprintln(bw, "# PSC\t[2]id\t[3]sample\t[4]nRefHom\t[5]nNonRefHom\t[6]nHets\t[7]nTransitions\t[8]nTransversions\t[9]nIndels\t[10]average depth\t[11]nSingletons\t[12]nHapRef\t[13]nHapAlt\t[14]nMissing")
	for i, name := range r.samples {
		var avg float64
		if r.pscDepthN[i] > 0 {
			avg = float64(r.pscDepthSum[i]) / float64(r.pscDepthN[i])
		}
		// Upstream uses %.1f for the mean depth (vcfstats.c:1795).
		fmt.Fprintf(bw, "PSC\t0\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.1f\t%d\t%d\t%d\t%d\n",
			name,
			r.pscNRefHom[i], r.pscNNonRefHom[i], r.pscNHets[i],
			r.pscNTs[i], r.pscNTv[i], r.pscNIndels[i], avg,
			r.pscNSingleton[i], r.pscNHapRef[i], r.pscNHapAlt[i], r.pscNMissing[i])
	}
	return nil
}

func writePSI(bw *bufio.Writer, r *statsResult) error {
	fmt.Fprintln(bw, "# PSI, Per-Sample Indels. Note that alt-het genotypes with both ins and del allele are counted twice, in both nInsHets and nDelHets.")
	fmt.Fprintln(bw, "# PSI\t[2]id\t[3]sample\t[4]in-frame\t[5]out-frame\t[6]not applicable\t[7]out/(in+out) ratio\t[8]nInsHets\t[9]nDelHets\t[10]nInsAltHoms\t[11]nDelAltHoms")
	for i, name := range r.samples {
		fmt.Fprintf(bw, "PSI\t0\t%s\t0\t0\t0\t0.00\t%d\t%d\t%d\t%d\n",
			name, r.psiInsHets[i], r.psiDelHets[i], r.psiInsHoms[i], r.psiDelHoms[i])
	}
	return nil
}

// writeHWE emits the Hardy-Weinberg section in upstream byte-parity
// format (vcfstats.c:1821-1858). For each non-empty AF bucket it prints
// total observations and the 25th/median/75th percentile of the
// per-record heterozygous-fraction CDF, reported as `j/naFHWE`.
func writeHWE(bw *bufio.Writer, r *statsResult) error {
	fmt.Fprintln(bw, "# HWE")
	fmt.Fprintln(bw, "# HWE\t[2]id\t[3]1st ALT allele frequency\t[4]Number of observations\t[5]25th percentile\t[6]median\t[7]75th percentile")
	hwe := r.afHWE
	for j := 0; j < naFHWE; j++ {
		hwe[1][j] += hwe[0][j]
	}
	for i := 1; i < mAFBins; i++ {
		var sumTot int
		for j := 0; j < naFHWE; j++ {
			sumTot += hwe[i][j]
		}
		if sumTot == 0 {
			continue
		}
		af := float64(i-1) / float64(mAFBins-1)
		var p25, p50, p75 float64
		havep25, havep50, havep75 := false, false, false
		sumTmp := 0
		for j := 0; j < naFHWE; j++ {
			sumTmp += hwe[i][j]
			frac := float64(sumTmp) / float64(sumTot)
			val := float64(j) / float64(naFHWE)
			if !havep25 && frac >= 0.25 {
				p25 = val
				havep25 = true
			}
			if !havep50 && frac >= 0.5 {
				p50 = val
				havep50 = true
			}
			if !havep75 && frac >= 0.75 {
				p75 = val
				havep75 = true
				break
			}
		}
		// Upstream's loop guarantees all three are emitted even when the
		// CDF collapses to a single bin (vcfstats.c:1844). Fill in any
		// remaining quantile from the last reached value.
		if !havep25 {
			p25 = p50
			if !havep50 {
				p25 = p75
			}
		}
		if !havep50 {
			p50 = p75
			if !havep75 {
				p50 = p25
			}
		}
		if !havep75 {
			p75 = p50
		}
		fmt.Fprintf(bw, "HWE\t0\t%f\t%d\t%f\t%f\t%f\n", af, sumTot, p25, p50, p75)
	}
	return nil
}

// unionMapKeys returns the sorted set of int keys present in any of the maps.
func unionMapKeys(ms ...map[int]int) []int {
	set := make(map[int]struct{})
	for _, m := range ms {
		for k := range m {
			set[k] = struct{}{}
		}
	}
	out := make([]int, 0, len(set))
	for k := range set {
		out = append(out, k)
	}
	return out
}

// ParseDepthSpec parses "MIN,MAX,STEP" into its three components. Empty
// fields fall back to the defaults (0, 500, 1).
func ParseDepthSpec(s string) (min, max, step int, err error) {
	min, max, step = 0, 500, 1
	if s == "" {
		return min, max, step, nil
	}
	parts := strings.Split(s, ",")
	if len(parts) != 3 {
		return 0, 0, 0, fmt.Errorf("bcftools stats: -d must be MIN,MAX,STEP, got %q", s)
	}
	if parts[0] != "" {
		n, err := strconv.Atoi(parts[0])
		if err != nil {
			return 0, 0, 0, fmt.Errorf("bcftools stats: bad MIN in -d: %w", err)
		}
		min = n
	}
	if parts[1] != "" {
		n, err := strconv.Atoi(parts[1])
		if err != nil {
			return 0, 0, 0, fmt.Errorf("bcftools stats: bad MAX in -d: %w", err)
		}
		max = n
	}
	if parts[2] != "" {
		n, err := strconv.Atoi(parts[2])
		if err != nil {
			return 0, 0, 0, fmt.Errorf("bcftools stats: bad STEP in -d: %w", err)
		}
		step = n
	}
	return min, max, step, nil
}

// ParseAFBins turns "0,0.1,0.5,1" into a slice of bin edges.
func ParseAFBins(s string) ([]float64, error) {
	if s == "" {
		return nil, nil
	}
	parts := strings.Split(s, ",")
	out := make([]float64, 0, len(parts))
	for _, p := range parts {
		p = strings.TrimSpace(p)
		if p == "" {
			continue
		}
		f, err := strconv.ParseFloat(p, 64)
		if err != nil {
			return nil, fmt.Errorf("bcftools stats: bad AF bin %q: %w", p, err)
		}
		out = append(out, f)
	}
	if len(out) < 2 {
		return nil, fmt.Errorf("bcftools stats: need at least 2 AF bin edges, got %d", len(out))
	}
	return out, nil
}
