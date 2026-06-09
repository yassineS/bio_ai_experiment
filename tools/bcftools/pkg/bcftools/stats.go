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
	AFBins          []float64      // -a; if nil we use the upstream default
	Collapse        string         // -c
	FirstAlleleOnly bool           // -1
	AFTag           string         // --af-tag
	UserTSTV        []UserTSTVSpec // -u/--user-tstv
	InputFile       string         // for the header line
}

// UserTSTVSpec describes a single `-u/--user-tstv TAG[:min:max:n]` request:
// collect transition/transversion counts for 1st-ALT SNPs stratified by the
// numeric INFO tag Tag, binned into NBins buckets spanning [Min,Max]. Idx
// selects an element of a multi-valued tag (e.g. PV4[1]); it defaults to 0.
type UserTSTVSpec struct {
	Tag   string  // INFO tag name
	Idx   int     // value index for multi-valued tags (default 0)
	Min   float64 // lower edge of the binning range (default 0)
	Max   float64 // upper edge of the binning range (default 1)
	NBins int     // number of bins (default 100)
}

// userTSTVAcc holds the per-spec transition/transversion bin counters that
// accumulate over the streamed records. vals_ts/vals_tv mirror upstream's
// uint64 arrays sized to NBins. isFloat tracks whether the INFO tag is of
// Float type (it governs the bin-value print format in the USR section).
type userTSTVAcc struct {
	spec    UserTSTVSpec
	valsTs  []uint64
	valsTv  []uint64
	isFloat bool
}

// defaultAFBins matches the upstream `bcftools stats` default of 11 bins
// covering [0, 1] with a special last bucket for AF=1.0 ([0.99,1.0]).
var defaultAFBins = []float64{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99, 1.0}

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

	// AF accumulators (one entry per AF bin)
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

	psiNIns []int
	psiNDel []int
	psiNHet []int
	psiNAA  []int

	// HWE: AF -> (n_obs, accumulated chi-square sum).
	hweObs    map[int]int
	hweChiSum map[int]float64

	// USR: per-spec Ts/Tv-by-tag accumulators (-u/--user-tstv).
	userTSTV []*userTSTVAcc
}

// newStatsResult prepares accumulators sized to the requested sample set.
// metaInfo carries the header's `##` meta lines so that USR (-u/--user-tstv)
// can resolve the numeric type of each requested INFO tag; pass nil when no
// USR specs are configured.
func newStatsResult(opts StatsOptions, headerSamples []string, metaInfo []string) *statsResult {
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

	r.samples = filterSampleSet(headerSamples, opts.Samples)
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
	r.psiNIns = make([]int, len(r.samples))
	r.psiNDel = make([]int, len(r.samples))
	r.psiNHet = make([]int, len(r.samples))
	r.psiNAA = make([]int, len(r.samples))

	for _, spec := range opts.UserTSTV {
		acc := &userTSTVAcc{
			spec:    spec,
			valsTs:  make([]uint64, spec.NBins),
			valsTv:  make([]uint64, spec.NBins),
			isFloat: infoTagIsFloat(metaInfo, spec.Tag),
		}
		r.userTSTV = append(r.userTSTV, acc)
	}
	return r
}

// infoTagIsFloat reports whether the INFO tag named tag is declared with
// Type=Float in the header meta lines. Tags declared Type=Integer (or any
// non-Float numeric) return false. Mirrors upstream's BCF_HT_REAL test which
// chooses the `%e` vs `%.0f` bin-value print format in the USR section.
func infoTagIsFloat(metaInfo []string, tag string) bool {
	prefix := "##INFO=<"
	for _, line := range metaInfo {
		if !strings.HasPrefix(line, prefix) || !strings.HasSuffix(line, ">") {
			continue
		}
		fields := parseStructuredMeta(line[len(prefix) : len(line)-1])
		if fields["ID"] == tag {
			return strings.EqualFold(fields["Type"], "Float")
		}
	}
	return false
}

// parseStructuredMeta splits the body of an angle-bracket VCF meta line
// (e.g. `ID=DP,Number=1,Type=Integer,Description="..."`) into its key/value
// pairs, honouring double-quoted values that may themselves contain commas.
func parseStructuredMeta(body string) map[string]string {
	out := make(map[string]string)
	var key, val strings.Builder
	inKey := true
	inQuote := false
	flush := func() {
		if key.Len() > 0 {
			out[key.String()] = val.String()
		}
		key.Reset()
		val.Reset()
		inKey = true
	}
	for i := 0; i < len(body); i++ {
		c := body[i]
		switch {
		case c == '"':
			inQuote = !inQuote
		case c == '=' && inKey && !inQuote:
			inKey = false
		case c == ',' && !inQuote:
			flush()
		case inKey:
			key.WriteByte(c)
		default:
			val.WriteByte(c)
		}
	}
	flush()
	return out
}

// filterSampleSet returns the intersection of `headerSamples` and `want`
// preserving the requested order. When want is empty all header samples are
// returned in their original order.
func filterSampleSet(headerSamples, want []string) []string {
	if len(want) == 0 {
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
	res := newStatsResult(opts, hdr.Samples, hdr.MetaInfo)
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
	res := newStatsResult(opts, hdr.Samples, hdr.MetaInfo)
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

	// QUAL binning. Upstream rounds DOWN to the nearest integer.
	qualBin := int(math.Floor(v.Qual))
	if v.Qual < 0 {
		qualBin = 0
	}

	// AF binning. If --af-tag is set we look it up in INFO; otherwise we
	// compute AF from genotypes.
	af := computeAF(r, v)
	afIdx := afBinIndex(r.afBins, af)

	for _, alt := range alts {
		if alt == "" || alt == "." || alt == "*" {
			continue
		}
		switch classifyVariant(v.Ref, alt) {
		case "snp":
			r.afSNPs[afIdx]++
			r.qualSNPs[qualBin]++
			tsType := transitionType(v.Ref, alt)
			if tsType == "ts" {
				r.afTs[afIdx]++
				r.qualTs[qualBin]++
			} else {
				r.afTv[afIdx]++
				r.qualTv[qualBin]++
			}
			key := strings.ToUpper(v.Ref) + ">" + strings.ToUpper(alt)
			r.subst[key]++
		case "indel":
			r.afNonS[afIdx]++
			r.qualNonS[qualBin]++
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

	// USR (-u/--user-tstv): stratify the 1st-ALT Ts/Tv by a numeric INFO
	// tag. Upstream only counts the first ALT allele and only when it is a
	// SNP whose REF differs from ALT.
	if len(r.userTSTV) > 0 && len(alts) > 0 {
		first := alts[0]
		if classifyVariant(v.Ref, first) == "snp" &&
			strings.ToUpper(v.Ref) != strings.ToUpper(first) {
			isTs := transitionType(v.Ref, first) == "ts"
			accumulateUserTSTV(r, v, isTs)
		}
	}

	// DP binning — INFO/DP first, falling back to per-sample DP sum.
	if dp, ok := siteDepth(v); ok {
		bin := dpBin(r.opts, dp)
		r.dpSites[bin]++
		r.dpTotalSite++
	}

	// Per-sample counters (PSC/PSI/HWE).
	accumulateSamples(r, v, alts)

	// HWE: 1st-ALT AF -> (obs, chi^2). We aggregate by AF bucket only for
	// records where both REF and ALT are SNPs (the upstream default).
	if sawSNP {
		chi, ok := hweChiSquare(v)
		if ok {
			bucket := int(math.Round(af * 1000)) // bucket by AF * 1000
			r.hweObs[bucket]++
			r.hweChiSum[bucket] += chi
		}
	}
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

// accumulateUserTSTV folds a 1st-ALT SNP into every configured USR
// accumulator. For each spec it reads the numeric INFO tag at the requested
// index, maps the value to a bin, and bumps the ts or tv counter. Records
// where the tag is absent or the index is out of range are skipped — exactly
// as upstream does.
func accumulateUserTSTV(r *statsResult, v *vcf.Variant, isTs bool) {
	for _, acc := range r.userTSTV {
		val, ok := infoTagValue(v, acc.spec.Tag, acc.spec.Idx)
		if !ok {
			continue
		}
		idx := userTSTVBin(acc.spec, val)
		if isTs {
			acc.valsTs[idx]++
		} else {
			acc.valsTv[idx]++
		}
	}
}

// infoTagValue returns the idx-th comma-separated numeric value of INFO tag in
// v. It reports false when the tag is missing, empty, has fewer than idx+1
// values, or the selected value is not a number.
func infoTagValue(v *vcf.Variant, tag string, idx int) (float64, bool) {
	raw, ok := v.Info[tag]
	if !ok || raw == "" || raw == "." {
		return 0, false
	}
	parts := strings.Split(raw, ",")
	if idx < 0 || idx >= len(parts) {
		return 0, false
	}
	f, err := strconv.ParseFloat(parts[idx], 64)
	if err != nil {
		return 0, false
	}
	return f, true
}

// userTSTVBin maps a value to its bin index following upstream's rule:
// values at or below Min land in bin 0, values at or above Max land in the
// last bin, and intermediate values scale linearly across (NBins-1) steps.
func userTSTVBin(spec UserTSTVSpec, val float64) int {
	switch {
	case val <= spec.Min:
		return 0
	case val >= spec.Max:
		return spec.NBins - 1
	default:
		idx := int((val - spec.Min) / (spec.Max - spec.Min) * float64(spec.NBins-1))
		if idx < 0 {
			idx = 0
		}
		if idx > spec.NBins-1 {
			idx = spec.NBins - 1
		}
		return idx
	}
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

// siteDepth returns INFO/DP if present (else falls back to the sum of FORMAT/DP).
func siteDepth(v *vcf.Variant) (int, bool) {
	if raw, ok := v.Info["DP"]; ok && raw != "" && raw != "." {
		if n, err := strconv.Atoi(raw); err == nil {
			return n, true
		}
	}
	total := 0
	any := false
	for _, s := range v.Samples {
		if raw, ok := s.Data["DP"]; ok && raw != "" && raw != "." {
			if n, err := strconv.Atoi(raw); err == nil {
				total += n
				any = true
			}
		}
	}
	if any {
		return total, true
	}
	return 0, false
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

// accumulateSamples updates the PSC / PSI / HWE per-sample counters for v.
func accumulateSamples(r *statsResult, v *vcf.Variant, alts []string) {
	if len(r.samples) == 0 {
		return
	}
	// Quick lookup of REF/ALT kinds for this site.
	altKinds := make([]string, len(alts))
	for i, alt := range alts {
		altKinds[i] = classifyVariant(v.Ref, alt)
	}
	for _, s := range v.Samples {
		idx, ok := r.sampleIndex[s.Name]
		if !ok {
			continue
		}
		gt, ok := s.Data["GT"]
		if !ok {
			continue
		}
		a, b, sep, ok := parseDiploidGT(gt)
		if !ok {
			continue
		}
		// Track depth contributions.
		if raw, ok := s.Data["DP"]; ok && raw != "" && raw != "." {
			if n, err := strconv.Atoi(raw); err == nil {
				r.pscDepthSum[idx] += n
				r.pscDepthN[idx]++
				bin := dpBin(r.opts, n)
				r.dpGTs[bin]++
				r.dpTotalGT++
			}
		}
		_ = sep
		switch {
		case a == 0 && b == 0:
			r.pscNRefHom[idx]++
		case a == b && a > 0:
			r.pscNNonRefHom[idx]++
			if int(a) <= len(altKinds) {
				kind := altKinds[a-1]
				if kind == "snp" {
					if transitionType(v.Ref, alts[a-1]) == "ts" {
						r.pscNTs[idx]++
					} else {
						r.pscNTv[idx]++
					}
				} else if kind == "indel" {
					r.pscNIndels[idx]++
					psiBumpIndel(r, idx, v.Ref, alts[a-1])
					r.psiNAA[idx]++
				}
			}
		case a != b:
			r.pscNHets[idx]++
			// For hets we attribute SNP/indel using the non-zero allele.
			nonZero := a
			if a == 0 {
				nonZero = b
			}
			if int(nonZero) >= 1 && int(nonZero) <= len(altKinds) {
				kind := altKinds[nonZero-1]
				if kind == "snp" {
					if transitionType(v.Ref, alts[nonZero-1]) == "ts" {
						r.pscNTs[idx]++
					} else {
						r.pscNTv[idx]++
					}
				} else if kind == "indel" {
					r.pscNIndels[idx]++
					psiBumpIndel(r, idx, v.Ref, alts[nonZero-1])
					r.psiNHet[idx]++
				}
			}
		}
	}
}

// psiBumpIndel updates the PSI insertion/deletion totals for sample idx.
func psiBumpIndel(r *statsResult, idx int, ref, alt string) {
	if len(alt) > len(ref) {
		r.psiNIns[idx]++
	} else if len(alt) < len(ref) {
		r.psiNDel[idx]++
	}
}

// parseDiploidGT splits a "0/1" / "0|1" string into the integer allele indices.
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

	// File header.
	fmt.Fprintln(bw, "# This file was produced by bcftools stats (pure-Go).")
	fmt.Fprintln(bw, "# The command line was:\tbcftools stats", r.opts.InputFile)

	// ID block — upstream lists each input file's ID; we only use one (0).
	fmt.Fprintln(bw, "# ID, Definition of sets:")
	fmt.Fprintf(bw, "ID\t0\t%s\n", r.opts.InputFile)

	if err := writeSN(bw, r); err != nil {
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
	if err := writeUSR(bw, r); err != nil {
		return err
	}
	if err := writeDP(bw, r); err != nil {
		return err
	}
	if err := writePSC(bw, r); err != nil {
		return err
	}
	if err := writePSI(bw, r); err != nil {
		return err
	}
	if err := writeHWE(bw, r); err != nil {
		return err
	}
	return bw.Flush()
}

func writeSN(bw *bufio.Writer, r *statsResult) error {
	fmt.Fprintln(bw, "# SN, Summary numbers:")
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

func writeAF(bw *bufio.Writer, r *statsResult) error {
	fmt.Fprintln(bw, "# AF, Stats by non-reference allele frequency:")
	fmt.Fprintln(bw, "# AF\t[2]id\t[3]allele frequency\t[4]number of SNPs\t[5]number of transitions\t[6]number of transversions\t[7]number of indels\t[8]repeat-consistent\t[9]repeat-inconsistent\t[10]not applicable")
	for i := 0; i < len(r.afBins)-1; i++ {
		lo := r.afBins[i]
		// Upstream prints the LOWER bin edge with 6 decimals.
		fmt.Fprintf(bw, "AF\t0\t%.6f\t%d\t%d\t%d\t%d\t0\t0\t0\n",
			lo, r.afSNPs[i], r.afTs[i], r.afTv[i], r.afNonS[i])
	}
	return nil
}

func writeQUAL(bw *bufio.Writer, r *statsResult) error {
	fmt.Fprintln(bw, "# QUAL, Stats by quality:")
	fmt.Fprintln(bw, "# QUAL\t[2]id\t[3]Quality\t[4]number of SNPs\t[5]number of transitions (1st ALT)\t[6]number of transversions (1st ALT)\t[7]number of indels")
	keys := unionMapKeys(r.qualSNPs, r.qualTs, r.qualTv, r.qualNonS)
	sort.Ints(keys)
	for _, k := range keys {
		fmt.Fprintf(bw, "QUAL\t0\t%d\t%d\t%d\t%d\t%d\n",
			k, r.qualSNPs[k], r.qualTs[k], r.qualTv[k], r.qualNonS[k])
	}
	return nil
}

func writeIDD(bw *bufio.Writer, r *statsResult) error {
	fmt.Fprintln(bw, "# IDD, InDel distribution:")
	fmt.Fprintln(bw, "# IDD\t[2]id\t[3]length (deletions negative)\t[4]number of sites\t[5]number of genotypes\t[6]mean VAF")
	keys := make([]int, 0, len(r.indelLen))
	for k := range r.indelLen {
		keys = append(keys, k)
	}
	sort.Ints(keys)
	for _, k := range keys {
		fmt.Fprintf(bw, "IDD\t0\t%d\t%d\t0\t0.00\n", k, r.indelLen[k])
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

// writeUSR emits the USR (user-tstv) sections, one per -u/--user-tstv spec.
// Each row reports a non-empty bin's [value, total SNPs, transitions,
// transversions]; bins with no observations are omitted, mirroring upstream.
func writeUSR(bw *bufio.Writer, r *statsResult) error {
	for _, acc := range r.userTSTV {
		spec := acc.spec
		label := fmt.Sprintf("%s/%d", spec.Tag, spec.Idx)
		fmt.Fprintf(bw, "# USR:%s\t[2]id\t[3]%s\t[4]number of SNPs\t[5]number of transitions (1st ALT)\t[6]number of transversions (1st ALT)\n", label, label)
		for j := 0; j < spec.NBins; j++ {
			total := acc.valsTs[j] + acc.valsTv[j]
			if total == 0 {
				continue // skip empty bins
			}
			var val float64
			if spec.NBins > 1 {
				val = spec.Min + (spec.Max-spec.Min)*float64(j)/float64(spec.NBins-1)
			} else {
				val = spec.Min
			}
			if acc.isFloat {
				fmt.Fprintf(bw, "USR:%s\t0\t%e\t%d\t%d\t%d\n", label, val, total, acc.valsTs[j], acc.valsTv[j])
			} else {
				fmt.Fprintf(bw, "USR:%s\t0\t%.0f\t%d\t%d\t%d\n", label, val, total, acc.valsTs[j], acc.valsTv[j])
			}
		}
	}
	return nil
}

func writeDP(bw *bufio.Writer, r *statsResult) error {
	fmt.Fprintln(bw, "# DP, Depth distribution:")
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
	fmt.Fprintln(bw, "# PSC, Per-sample counts:")
	fmt.Fprintln(bw, "# PSC\t[2]id\t[3]sample\t[4]nRefHom\t[5]nNonRefHom\t[6]nHets\t[7]nTransitions\t[8]nTransversions\t[9]nIndels\t[10]average depth\t[11]nSingletons\t[12]nHapRef\t[13]nHapAlt\t[14]nMissing")
	for i, name := range r.samples {
		var avg float64
		if r.pscDepthN[i] > 0 {
			avg = float64(r.pscDepthSum[i]) / float64(r.pscDepthN[i])
		}
		fmt.Fprintf(bw, "PSC\t0\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.2f\t0\t0\t0\t0\n",
			name,
			r.pscNRefHom[i], r.pscNNonRefHom[i], r.pscNHets[i],
			r.pscNTs[i], r.pscNTv[i], r.pscNIndels[i], avg)
	}
	return nil
}

func writePSI(bw *bufio.Writer, r *statsResult) error {
	fmt.Fprintln(bw, "# PSI, Per-Sample Indels:")
	fmt.Fprintln(bw, "# PSI\t[2]id\t[3]sample\t[4]in-frame\t[5]out-frame\t[6]not applicable\t[7]out/(in+out) ratio\t[8]nInsHets\t[9]nDelHets\t[10]nInsAltHoms\t[11]nDelAltHoms")
	for i, name := range r.samples {
		fmt.Fprintf(bw, "PSI\t0\t%s\t0\t0\t0\t0.00\t%d\t%d\t%d\t%d\n",
			name, r.psiNIns[i], r.psiNDel[i], r.psiNHet[i], r.psiNAA[i])
	}
	return nil
}

func writeHWE(bw *bufio.Writer, r *statsResult) error {
	fmt.Fprintln(bw, "# HWE, Hardy-Weinberg equilibrium:")
	fmt.Fprintln(bw, "# HWE\t[2]id\t[3]1st ALT allele frequency\t[4]Number of observations\t[5]25th percentile\t[6]median\t[7]75th percentile")
	keys := make([]int, 0, len(r.hweObs))
	for k := range r.hweObs {
		keys = append(keys, k)
	}
	sort.Ints(keys)
	for _, k := range keys {
		af := float64(k) / 1000.0
		obs := r.hweObs[k]
		avg := 0.0
		if obs > 0 {
			avg = r.hweChiSum[k] / float64(obs)
		}
		fmt.Fprintf(bw, "HWE\t0\t%.6f\t%d\t%.6f\t%.6f\t%.6f\n", af, obs, avg, avg, avg)
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

// ParseUserTSTV parses a `-u/--user-tstv` argument of the form
// `TAG[:min:max:n]`. The TAG portion may carry a value index, e.g. `PV4[1]`.
// Omitted binning fields default to min=0, max=1, n=100 (matching upstream
// bcftools). Returns an error for an empty tag, a malformed index, an
// unparseable number, or a non-positive bin count.
func ParseUserTSTV(s string) (UserTSTVSpec, error) {
	spec := UserTSTVSpec{Min: 0, Max: 1, NBins: 100, Idx: 0}
	if s == "" {
		return spec, fmt.Errorf("bcftools stats: -u/--user-tstv requires a TAG")
	}
	tag := s
	if i := strings.IndexByte(s, ':'); i >= 0 {
		tag = s[:i]
		rest := strings.Split(s[i+1:], ":")
		if len(rest) > 3 {
			return spec, fmt.Errorf("bcftools stats: -u/--user-tstv: too many fields in %q", s)
		}
		if len(rest) >= 1 && rest[0] != "" {
			f, err := strconv.ParseFloat(rest[0], 64)
			if err != nil {
				return spec, fmt.Errorf("bcftools stats: -u/--user-tstv: bad min in %q: %w", s, err)
			}
			spec.Min = f
		}
		if len(rest) >= 2 && rest[1] != "" {
			f, err := strconv.ParseFloat(rest[1], 64)
			if err != nil {
				return spec, fmt.Errorf("bcftools stats: -u/--user-tstv: bad max in %q: %w", s, err)
			}
			spec.Max = f
		}
		if len(rest) >= 3 && rest[2] != "" {
			n, err := strconv.Atoi(rest[2])
			if err != nil {
				return spec, fmt.Errorf("bcftools stats: -u/--user-tstv: bad n in %q: %w", s, err)
			}
			if n <= 0 {
				return spec, fmt.Errorf("bcftools stats: -u/--user-tstv: number of bins must be positive in %q", s)
			}
			spec.NBins = n
		}
	}
	// Optional value index: TAG[idx].
	if strings.HasSuffix(tag, "]") {
		open := strings.LastIndexByte(tag, '[')
		if open < 0 {
			return spec, fmt.Errorf("bcftools stats: -u/--user-tstv: malformed index in %q", s)
		}
		idxStr := tag[open+1 : len(tag)-1]
		n, err := strconv.Atoi(idxStr)
		if err != nil || n < 0 {
			return spec, fmt.Errorf("bcftools stats: -u/--user-tstv: bad index in %q", s)
		}
		spec.Idx = n
		tag = tag[:open]
	}
	if tag == "" {
		return spec, fmt.Errorf("bcftools stats: -u/--user-tstv requires a TAG in %q", s)
	}
	spec.Tag = tag
	return spec, nil
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
