package bcftools

import (
	"bufio"
	"fmt"
	"io"
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

	// useSamples reports whether per-sample stats were requested (-s/-S),
	// which enables the PSC/PSI/HWE sections and per-genotype DP counts.
	useSamples bool

	// AF accumulators. Bins are indexed 0..mAF-1; bin 0 holds singletons and
	// is folded into bin 1 at print time, matching upstream. The printed AF
	// value for bin i (i>=1) is (i-1)/(mAF-1). afBins is non-nil only when a
	// custom --af-bins set was supplied.
	mAF    int
	afBins []float64
	afSNPs []int
	afTs   []int
	afTv   []int
	afNA   []int // indel/repeat "not applicable" counter (upstream af_repeats[2])

	// QUAL accumulators keyed by the upstream iqual bucket index
	// (1+int(QUAL*10), or 0 for missing/negative QUAL).
	qualTs     map[int]int
	qualTv     map[int]int
	qualIndels map[int]int

	// IDD: indel length distribution (length, count) — len < 0 = deletion.
	indelLen map[int]int

	// ST: substitution-type counts (REF>ALT for SNPs), indexed by ref<<2|alt.
	subst [16]int

	// DP: depth distribution. Bins follow upstream's idist layout: index 0 is
	// the underflow (<min) bucket and the final index is the overflow (>max)
	// bucket; the interior bin i maps to depth (i-1+min).
	dpMin, dpMax, dpStep int
	dpGTBins             []uint64 // per-genotype counts (needs -s/-S)
	dpSiteBins           []uint64 // per-site counts (INFO/DP)

	// PSC / PSI: per-sample counters indexed by sampleIndex.
	pscNRefHom    []int
	pscNNonRefHom []int
	pscNHets      []int
	pscNTs        []int
	pscNTv        []int
	pscNIndels    []int
	pscNSingleton []int
	pscNHapRef    []int
	pscNHapAlt    []int
	pscNMissing   []int
	pscDepthSum   []int
	pscDepthN     []int

	psiNInsHet []int
	psiNDelHet []int
	psiNInsHom []int
	psiNDelHom []int

	// HWE: het-fraction distribution keyed by [afBin][hetFreqBin], mirroring
	// upstream's af_hwe array. nafHWE is the number of het-fraction bins.
	nafHWE int
	afHWE  [][]int

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
		qualTs:        make(map[int]int),
		qualTv:        make(map[int]int),
		qualIndels:    make(map[int]int),
		indelLen:      make(map[int]int),
	}
	r.useSamples = len(opts.Samples) > 0 || opts.SamplesFile != ""

	// AF bins: m_af defaults to 101 (the [0,1] range split into 100 intervals),
	// growing to nsamples+1 when more samples are present so that low allele
	// frequencies still resolve. A user --af-bins set overrides the layout.
	if len(r.afBins) > 0 {
		r.mAF = len(r.afBins)
	} else {
		r.mAF = 101
		if r.useSamples && len(headerSamples)+1 > r.mAF {
			r.mAF = len(headerSamples) + 1
		}
	}
	r.afSNPs = make([]int, r.mAF)
	r.afTs = make([]int, r.mAF)
	r.afTv = make([]int, r.mAF)
	r.afNA = make([]int, r.mAF)

	// DP idist layout: index 0 = underflow, last index = overflow.
	r.dpMin, r.dpMax, r.dpStep = opts.DepthMin, opts.DepthMax, opts.DepthStep
	if r.dpMax <= 0 {
		r.dpMin, r.dpMax, r.dpStep = 0, 500, 1
	}
	if r.dpStep <= 0 {
		r.dpStep = 1
	}
	dpNBins := 4 + (r.dpMax-r.dpMin)/r.dpStep
	r.dpGTBins = make([]uint64, dpNBins)
	r.dpSiteBins = make([]uint64, dpNBins)

	if r.useSamples {
		r.samples = filterSampleSet(headerSamples, opts.Samples)
	}
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
	r.pscNSingleton = make([]int, len(r.samples))
	r.pscNHapRef = make([]int, len(r.samples))
	r.pscNHapAlt = make([]int, len(r.samples))
	r.pscNMissing = make([]int, len(r.samples))
	r.pscDepthSum = make([]int, len(r.samples))
	r.pscDepthN = make([]int, len(r.samples))
	r.psiNInsHet = make([]int, len(r.samples))
	r.psiNDelHet = make([]int, len(r.samples))
	r.psiNInsHom = make([]int, len(r.samples))
	r.psiNDelHom = make([]int, len(r.samples))

	r.nafHWE = 100
	if len(r.samples) > 0 {
		r.afHWE = make([][]int, r.mAF)
		for i := range r.afHWE {
			r.afHWE[i] = make([]int, r.nafHWE)
		}
	}

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
	// A bare "-" (as in `-s -`) selects every sample in the header, matching
	// upstream's convention; so does an empty request.
	if len(want) == 0 || (len(want) == 1 && want[0] == "-") {
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
	includeF, excludeF, err := compileStatsExpressions(opts, hdr)
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
	includeF, excludeF, err := compileStatsExpressions(opts, hdr)
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

// compileStatsExpressions reuses the view-side expression compiler, resolving
// bare identifiers against hdr (FORMAT-only tags become FORMAT; a tag declared
// as both INFO and FORMAT is rejected as ambiguous) exactly as upstream does.
func compileStatsExpressions(opts StatsOptions, hdr *vcf.Header) (include, exclude *Filter, err error) {
	if opts.IncludeExpr != "" {
		include, err = CompileFilterWithHeader(opts.IncludeExpr, hdr)
		if err != nil {
			return nil, nil, err
		}
	}
	if opts.ExcludeExpr != "" {
		exclude, err = CompileFilterWithHeader(opts.ExcludeExpr, hdr)
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

	alts := v.Alt
	if r.opts.FirstAlleleOnly && len(alts) > 1 {
		alts = alts[:1]
	}

	// SN counters: each variant type present in the line bumps its counter
	// independently, mirroring upstream's do_vcf_stats dispatch.
	lineType := lineVariantTypes(v)
	if lineType == vcfRef {
		r.numNoALTs++
	}
	if lineType&vcfSNP != 0 {
		r.numSNPs++
	}
	if lineType&vcfMNP != 0 {
		r.numMNPs++
	}
	if lineType&vcfIndel != 0 {
		r.numIndels++
	}
	if lineType&vcfOther != 0 {
		r.numOthers++
	}
	if len(v.Alt) > 1 {
		r.numMA++
		if lineType == vcfSNP {
			r.numMASNP++
		}
	}

	// QUAL bucket index, mirroring upstream's iqual = 1 + int(QUAL*10) with a
	// 0 bucket reserved for missing/negative QUAL.
	iqual := 0
	if v.Qual >= 0 {
		iqual = 1 + int(v.Qual*10)
	}

	// Per-allele AF bin indices (iaf), matching init_iaf: AC==1 → singleton
	// bin 0, otherwise int(af*(mAF-2))+1. iaf[0] (REF) is unused.
	iaf := r.allelesIAF(v)

	// do_snp_stats / do_indel_stats: iterate ALT alleles, bucketing each into
	// the AF/QUAL/ST/IDD accumulators by its precise variant type.
	for i, alt := range alts {
		t, n := statsVariantType(v.Ref, alt)
		ai := 0
		if i+1 < len(iaf) {
			ai = iaf[i+1]
		}
		if t&vcfSNP != 0 {
			refI := acgt2int(v.Ref[0])
			altI := acgt2int(alt[0])
			if refI < 0 || altI < 0 || refI == altI {
				continue
			}
			r.subst[refI<<2|altI]++
			r.afSNPs[ai]++
			isTs := abs(refI-altI) == 2
			if isTs {
				r.afTs[ai]++
			} else {
				r.afTv[ai]++
			}
			if i == 0 { // 1st ALT drives QUAL ts/tv and USR
				if isTs {
					r.qualTs[iqual]++
					if len(r.userTSTV) > 0 {
						accumulateUserTSTV(r, v, true)
					}
				} else {
					r.qualTv[iqual]++
					if len(r.userTSTV) > 0 {
						accumulateUserTSTV(r, v, false)
					}
				}
			}
		} else if t&vcfIndel != 0 {
			// IDD records counts the indel once per record (not per allele):
			// the qual_indels bucket is bumped in the per-record block below.
			r.afNA[ai]++
			length := n
			if length < 0 {
				length = -length
			}
			length--
			if length >= 60 {
				length = 59
			}
			if n < 0 {
				r.indelLen[-(length+1)]++
			} else {
				r.indelLen[length+1]++
			}
		}
	}
	// qual_indels is incremented once per indel record (do_indel_stats).
	if lineType&vcfIndel != 0 {
		r.qualIndels[iqual]++
	}

	// DP, sites: upstream uses INFO/DP only, and only a single integer value.
	if raw, ok := v.Info["DP"]; ok && raw != "" && raw != "." {
		if dp, err := strconv.Atoi(strings.TrimSpace(raw)); err == nil {
			r.dpSiteBins[idistIndex(r, dp)]++
		}
	}

	// Per-sample counters (PSC/PSI/HWE) and per-genotype DP.
	if len(r.samples) > 0 {
		accumulateSamples(r, v, iaf)
	}
}

// allelesIAF computes the per-allele AF bin index for a record, mirroring
// htslib's init_iaf. Index 0 (REF) is unused; ALT i lands in iaf[i].
func (r *statsResult) allelesIAF(v *vcf.Variant) []int {
	nAllele := len(v.Alt) + 1
	iaf := make([]int, nAllele)

	// --af-tag: bin by the INFO float tag value directly.
	if r.opts.AFTag != "" {
		if raw, ok := v.Info[r.opts.AFTag]; ok && raw != "" {
			parts := strings.Split(raw, ",")
			if len(parts) == nAllele-1 {
				for i, p := range parts {
					af, err := strconv.ParseFloat(strings.TrimSpace(p), 64)
					if err != nil {
						return make([]int, nAllele)
					}
					iaf[i+1] = r.afBinIndex(af) + 1
				}
				return iaf
			}
		}
		return iaf
	}

	ac, ok := statsAC(v, r.useSamples)
	if !ok {
		return iaf
	}
	an := 0
	for _, c := range ac {
		an += c
	}
	for i := 1; i < nAllele; i++ {
		switch {
		case ac[i] == 1:
			iaf[i] = 0 // singleton
		case an == 0:
			iaf[i] = 1
		default:
			af := float64(ac[i]) / float64(an)
			iaf[i] = r.afBinIndex(af) + 1
		}
	}
	return iaf
}

// afBinIndex maps an allele frequency in [0,1] to the upstream AF bin index
// (before the +1 singleton offset). It honours a custom --af-bins layout when
// one was supplied.
func (r *statsResult) afBinIndex(af float64) int {
	if af < 0 {
		af = 0
	} else if af > 1 {
		af = 1
	}
	if len(r.afBins) > 0 {
		return binGetIdx(r.afBins, af)
	}
	return int(af * float64(r.mAF-2))
}

// idistIndex returns the idist bin for a depth value: 0 underflow, last index
// overflow, interior 1+(val-min)/step.
func idistIndex(r *statsResult, val int) int {
	if val < r.dpMin {
		return 0
	}
	if val > r.dpMax {
		return len(r.dpSiteBins) - 1
	}
	return 1 + (val-r.dpMin)/r.dpStep
}

// binGetIdx returns the interval index in a sorted boundary slice that contains
// af, matching htslib's bin_get_idx for custom --af-bins.
func binGetIdx(bins []float64, af float64) int {
	if af <= bins[0] {
		return 0
	}
	if af >= bins[len(bins)-1] {
		return len(bins) - 2
	}
	for i := 1; i < len(bins); i++ {
		if af < bins[i] {
			return i - 1
		}
	}
	return len(bins) - 2
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

// accumulateSamples updates the PSC / PSI / HWE per-sample counters and the
// per-genotype DP distribution for v, mirroring htslib's do_sample_stats. iaf
// holds the per-allele AF bin indices computed by allelesIAF.
func accumulateSamples(r *statsResult, v *vcf.Variant, iaf []int) {
	// Per-allele variant types for this site.
	altType := make([]int, len(v.Alt)+1)
	altLen := make([]int, len(v.Alt)+1)
	for i, alt := range v.Alt {
		altType[i+1], altLen[i+1] = statsVariantType(v.Ref, alt)
	}

	nNonRef := 0
	iNonRef := -1
	var nRefTot, nHetTot, nAltTot int

	for _, s := range v.Samples {
		idx, ok := r.sampleIndex[s.Name]
		if !ok {
			continue
		}
		// Depth: upstream prefers FORMAT/DP (calc_sample_depth).
		if raw, ok := s.Data["DP"]; ok && raw != "" && raw != "." {
			if n, err := strconv.Atoi(strings.TrimSpace(raw)); err == nil && n > 0 {
				r.dpGTBins[idistIndex(r, n)]++
				r.pscDepthSum[idx] += n
				r.pscDepthN[idx]++
			}
		}

		gtRaw, ok := s.Data["GT"]
		if !ok {
			continue
		}
		gt, ial, jal := classifyGT(gtRaw)
		if gt == gtUnkn {
			r.pscNMissing[idx]++
			continue
		}

		varType := 0
		if ial > 0 && ial < len(altType) {
			varType |= altType[ial]
		}
		if jal > 0 && jal < len(altType) {
			varType |= altType[jal]
		}

		if gt == gtHaplR || gt == gtHaplA {
			if gt == gtHaplR {
				r.pscNHapRef[idx]++
			} else {
				r.pscNHapAlt[idx]++
			}
			continue
		}

		if gt != gtHomRR {
			nNonRef++
			iNonRef = idx
		}
		switch gt {
		case gtHomRR:
			nRefTot++
		case gtHetRA:
			nHetTot++
		case gtHetAA, gtHomAA:
			nAltTot++
		}

		// nRefHom/nNonRefHom/nHets and ts/tv only for SNP-typed (or REF) sites.
		if varType&vcfSNP != 0 || varType == vcfRef {
			switch gt {
			case gtHetRA, gtHetAA:
				r.pscNHets[idx]++
			case gtHomRR:
				r.pscNRefHom[idx]++
			case gtHomAA:
				r.pscNNonRefHom[idx]++
			}
			if gt != gtHomRR && ial < len(altType) && altType[ial]&vcfSNP != 0 {
				if tstv := snpTsTv(v.Ref, v.Alt, ial); tstv == 1 {
					r.pscNTs[idx]++
				} else if tstv == -1 {
					r.pscNTv[idx]++
				}
			}
			if gt != gtHomRR && ial != jal && jal < len(altType) && altType[jal]&vcfSNP != 0 {
				if tstv := snpTsTv(v.Ref, v.Alt, jal); tstv == 1 {
					r.pscNTs[idx]++
				} else if tstv == -1 {
					r.pscNTv[idx]++
				}
			}
		}

		// Indel-typed genotypes feed nIndels and the PSI ins/del counters.
		if varType&vcfIndel != 0 && gt != gtHomRR {
			r.pscNIndels[idx]++
			switch gt {
			case gtHetRA, gtHetAA:
				isIns, isDel := false, false
				if ial > 0 && ial < len(altType) && altType[ial]&vcfIndel != 0 {
					if altLen[ial] < 0 {
						isDel = true
					} else {
						isIns = true
					}
				}
				if jal > 0 && jal < len(altType) && altType[jal]&vcfIndel != 0 {
					if altLen[jal] < 0 {
						isDel = true
					} else {
						isIns = true
					}
				}
				if isDel {
					r.psiNDelHet[idx]++
				}
				if isIns {
					r.psiNInsHet[idx]++
				}
			case gtHomAA:
				if ial > 0 && ial < len(altType) {
					if altLen[ial] < 0 {
						r.psiNDelHom[idx]++
					} else {
						r.psiNInsHom[idx]++
					}
				}
			}
		}
	}

	if nNonRef == 1 && iNonRef >= 0 {
		r.pscNSingleton[iNonRef]++
	}

	// HWE: record the observed het fraction for the 1st ALT's AF bin.
	if r.afHWE != nil && len(v.Alt) >= 1 && (nRefTot != 0 || nHetTot != 0 || nAltTot != 0) {
		total := nRefTot + nHetTot + nAltTot
		hetFrac := float64(nHetTot) / float64(total)
		ihet := int(hetFrac * float64(r.nafHWE-1))
		ai := 0
		if len(iaf) > 1 {
			ai = iaf[1]
		}
		if ai >= 0 && ai < len(r.afHWE) && ihet >= 0 && ihet < r.nafHWE {
			r.afHWE[ai][ihet]++
		}
	}
}

// snpTsTv reports whether the SNP REF->ALT[ial] is a transition (1),
// transversion (-1), or undefined (0), using the htslib base ordering where a
// distance of 2 (A<->G, C<->T) marks a transition.
func snpTsTv(ref string, alts []string, ial int) int {
	if ial < 1 || ial > len(alts) {
		return 0
	}
	refI := acgt2int(ref[0])
	altI := acgt2int(alts[ial-1][0])
	if refI < 0 || altI < 0 {
		return 0
	}
	if abs(refI-altI) == 2 {
		return 1
	}
	return -1
}

// writeStats emits the upstream-style tab-prefixed report.
func writeStats(out io.Writer, r *statsResult) error {
	bw := bufio.NewWriter(out)

	// File header.
	fmt.Fprintln(bw, "# This file was produced by bcftools stats (pure-Go).")
	fmt.Fprintln(bw, "# The command line was:\tbcftools stats", r.opts.InputFile)
	fmt.Fprintln(bw, "#")

	// ID block — upstream lists each input file's ID; we only use one (0).
	// The input filename "-" is rendered as <STDIN>, matching upstream.
	fmt.Fprintln(bw, "# Definition of sets:")
	fmt.Fprintln(bw, "# ID\t[2]id\t[3]tab-separated file names")
	fname := r.opts.InputFile
	if fname == "-" {
		fname = "<STDIN>"
	}
	fmt.Fprintf(bw, "ID\t0\t%s\n", fname)

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
	if err := writeUSR(bw, r); err != nil {
		return err
	}
	if err := writeDP(bw, r); err != nil {
		return err
	}
	// PSC/PSI/HWE are only emitted when per-sample stats were requested
	// (-s/-S), mirroring upstream's `if (files->n_smpl)` guard.
	if len(r.samples) > 0 {
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

// writeTSTV emits the TSTV (transitions/transversions) section. The ts and tv
// totals are the sums of the AF-binned per-allele tallies over every bin
// (including the singleton bin 0), exactly as upstream's print_stats does
// before folding bin 0 into bin 1. The "1st ALT" columns are the first-ALT-only
// tallies — the same counters that drive the QUAL section — summed across all
// quality buckets. This must run before writeSiS/writeAF, which mutate the AF
// bins via the singleton fold.
func writeTSTV(bw *bufio.Writer, r *statsResult) error {
	fmt.Fprintln(bw, "# TSTV, transitions/transversions")
	fmt.Fprintln(bw, "#   - transitions, see https://en.wikipedia.org/wiki/Transition_(genetics)")
	fmt.Fprintln(bw, "#   - transversions, see https://en.wikipedia.org/wiki/Transversion")
	fmt.Fprintln(bw, "# TSTV\t[2]id\t[3]ts\t[4]tv\t[5]ts/tv\t[6]ts (1st ALT)\t[7]tv (1st ALT)\t[8]ts/tv (1st ALT)")
	ts, tv := 0, 0
	for i := 0; i < r.mAF; i++ {
		ts += r.afTs[i]
		tv += r.afTv[i]
	}
	tsAlt1, tvAlt1 := 0, 0
	for _, n := range r.qualTs {
		tsAlt1 += n
	}
	for _, n := range r.qualTv {
		tvAlt1 += n
	}
	var ratio, ratio1 float64
	if tv != 0 {
		ratio = float64(ts) / float64(tv)
	}
	if tvAlt1 != 0 {
		ratio1 = float64(tsAlt1) / float64(tvAlt1)
	}
	fmt.Fprintf(bw, "TSTV\t0\t%d\t%d\t%.2f\t%d\t%d\t%.2f\n", ts, tv, ratio, tsAlt1, tvAlt1, ratio1)
	return nil
}

// writeSiS emits the SiS (singleton stats) section. It reports the contents of
// the AF singleton bucket (bin 0) BEFORE that bucket is folded into bin 1 by
// writeAF. The allele count column is always 1 (the singleton AC). The repeat
// columns are emitted as zero: this port does not compute the experimental,
// deprecated indel-context repeat tallies, so the not-applicable bucket carries
// every singleton indel — matching upstream when --indel-context is off.
func writeSiS(bw *bufio.Writer, r *statsResult) error {
	fmt.Fprintln(bw, "# SiS, Singleton stats:")
	fmt.Fprintln(bw, "#   - allele count, i.e. the number of singleton genotypes (AC=1)")
	fmt.Fprintln(bw, "#   - number of transitions, see above")
	fmt.Fprintln(bw, "#   - number of transversions, see above")
	fmt.Fprintln(bw, "#   - repeat-consistent, inconsistent and n/a: experimental and useless stats [DEPRECATED]")
	fmt.Fprintln(bw, "# SiS\t[2]id\t[3]allele count\t[4]number of SNPs\t[5]number of transitions\t[6]number of transversions\t[7]number of indels\t[8]repeat-consistent\t[9]repeat-inconsistent\t[10]not applicable")
	// af_repeats[0][0]+af_repeats[1][0]+af_repeats[2][0] == total singleton
	// indels; here that total is r.afNA[0] and the consistent/inconsistent
	// buckets are zero, so the "not applicable" column equals r.afNA[0].
	fmt.Fprintf(bw, "SiS\t0\t1\t%d\t%d\t%d\t%d\t0\t0\t%d\n",
		r.afSNPs[0], r.afTs[0], r.afTv[0], r.afNA[0], r.afNA[0])
	return nil
}

func writeAF(bw *bufio.Writer, r *statsResult) error {
	fmt.Fprintln(bw, "# AF, Stats by non-reference allele frequency:")
	fmt.Fprintln(bw, "# AF\t[2]id\t[3]allele frequency\t[4]number of SNPs\t[5]number of transitions\t[6]number of transversions\t[7]number of indels\t[8]repeat-consistent\t[9]repeat-inconsistent\t[10]not applicable")
	// The singleton bin (0) is folded into bin 1, exactly as upstream does after
	// emitting the SiS line.
	if r.mAF > 1 {
		r.afSNPs[1] += r.afSNPs[0]
		r.afTs[1] += r.afTs[0]
		r.afTv[1] += r.afTv[0]
		r.afNA[1] += r.afNA[0]
	}
	for i := 1; i < r.mAF; i++ {
		if r.afSNPs[i]+r.afTs[i]+r.afTv[i]+r.afNA[i] == 0 {
			continue
		}
		fmt.Fprintf(bw, "AF\t0\t%f\t%d\t%d\t%d\t%d\t0\t0\t%d\n",
			r.afValue(i), r.afSNPs[i], r.afTs[i], r.afTv[i], r.afNA[i], r.afNA[i])
	}
	return nil
}

// afValue returns the printed allele-frequency label for AF bin i. With the
// default linear bins it is (i-1)/(mAF-1); with custom --af-bins it is the
// midpoint of the surrounding boundaries.
func (r *statsResult) afValue(i int) float64 {
	if len(r.afBins) > 0 {
		return (r.afBins[i] + r.afBins[i-1]) * 0.5
	}
	return float64(i-1) / float64(r.mAF-1)
}

func writeQUAL(bw *bufio.Writer, r *statsResult) error {
	fmt.Fprintln(bw, "# QUAL, Stats by quality")
	fmt.Fprintln(bw, "# QUAL\t[2]id\t[3]Quality\t[4]number of SNPs\t[5]number of transitions (1st ALT)\t[6]number of transversions (1st ALT)\t[7]number of indels")
	keys := unionMapKeys(r.qualTs, r.qualTv, r.qualIndels)
	sort.Ints(keys)
	for _, k := range keys {
		nts, ntv, nin := r.qualTs[k], r.qualTv[k], r.qualIndels[k]
		if nts+ntv+nin == 0 {
			continue
		}
		// Upstream prints "." for the missing/negative QUAL bucket (k==0) and
		// 0.1*(k-1) otherwise, always with one decimal place.
		if k == 0 {
			fmt.Fprintf(bw, "QUAL\t0\t.\t%d\t%d\t%d\t%d\n", nts+ntv, nts, ntv, nin)
		} else {
			fmt.Fprintf(bw, "QUAL\t0\t%.1f\t%d\t%d\t%d\t%d\n", 0.1*float64(k-1), nts+ntv, nts, ntv, nin)
		}
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
	// Upstream emits "." for the mean-VAF column when no AD-derived value is
	// available (the common case without FORMAT/AD).
	for _, k := range keys {
		fmt.Fprintf(bw, "IDD\t0\t%d\t%d\t0\t.\n", k, r.indelLen[k])
	}
	return nil
}

func writeST(bw *bufio.Writer, r *statsResult) error {
	fmt.Fprintln(bw, "# ST, Substitution types:")
	fmt.Fprintln(bw, "# ST\t[2]id\t[3]type\t[4]count")
	// Iterate ref<<2|alt the way upstream does, skipping the four ref==alt
	// codes (t>>2 == t&3).
	for t := 0; t < 16; t++ {
		if t>>2 == t&3 {
			continue
		}
		fmt.Fprintf(bw, "ST\t0\t%c>%c\t%d\n", "ACGT"[t>>2], "ACGT"[t&3], r.subst[t])
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
	fmt.Fprintln(bw, "# DP, depth:")
	fmt.Fprintln(bw, "#   - set id, see above")
	fmt.Fprintln(bw, "#   - the depth bin, corresponds to the depth (unless --depth was given)")
	fmt.Fprintln(bw, "#   - number of genotypes with this depth (zero unless -s/-S was given)")
	fmt.Fprintln(bw, "#   - fraction of genotypes with this depth (zero unless -s/-S was given)")
	fmt.Fprintln(bw, "#   - number of sites with this depth")
	fmt.Fprintln(bw, "#   - fraction of sites with this depth")
	fmt.Fprintln(bw, "# DP, Depth distribution")
	fmt.Fprintln(bw, "# DP\t[2]id\t[3]bin\t[4]number of genotypes\t[5]fraction of genotypes (%)\t[6]number of sites\t[7]fraction of sites (%)")
	var sumGT, sumSite uint64
	for i := range r.dpGTBins {
		sumGT += r.dpGTBins[i]
		sumSite += r.dpSiteBins[i]
	}
	last := len(r.dpGTBins) - 1
	for i := range r.dpGTBins {
		if r.dpGTBins[i] == 0 && r.dpSiteBins[i] == 0 {
			continue
		}
		var label string
		switch {
		case i == 0:
			label = fmt.Sprintf("<%d", r.dpMin)
		case i == last:
			label = fmt.Sprintf(">%d", r.dpMax)
		default:
			label = strconv.Itoa(i - 1 + r.dpMin)
		}
		var fracGT, fracSite float64
		if sumGT > 0 {
			fracGT = float64(r.dpGTBins[i]) * 100.0 / float64(sumGT)
		}
		if sumSite > 0 {
			fracSite = float64(r.dpSiteBins[i]) * 100.0 / float64(sumSite)
		}
		fmt.Fprintf(bw, "DP\t0\t%s\t%d\t%f\t%d\t%f\n",
			label, r.dpGTBins[i], fracGT, r.dpSiteBins[i], fracSite)
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
		// in-frame/out-frame/na require an exons file (--exons), which we do
		// not model; upstream emits zeros for them otherwise.
		fmt.Fprintf(bw, "PSI\t0\t%s\t0\t0\t0\t0.00\t%d\t%d\t%d\t%d\n",
			name, r.psiNInsHet[i], r.psiNDelHet[i], r.psiNInsHom[i], r.psiNDelHom[i])
	}
	return nil
}

func writeHWE(bw *bufio.Writer, r *statsResult) error {
	fmt.Fprintln(bw, "# HWE")
	fmt.Fprintln(bw, "# HWE\t[2]id\t[3]1st ALT allele frequency\t[4]Number of observations\t[5]25th percentile\t[6]median\t[7]75th percentile")
	if r.afHWE == nil {
		return nil
	}
	// Fold the singleton AF bin into bin 1, matching upstream.
	for j := 0; j < r.nafHWE; j++ {
		r.afHWE[1][j] += r.afHWE[0][j]
	}
	for i := 1; i < r.mAF; i++ {
		ptr := r.afHWE[i]
		sumTot := 0
		for _, c := range ptr {
			sumTot += c
		}
		if sumTot == 0 {
			continue
		}
		af := r.afValue(i)
		nprn := 3
		fmt.Fprintf(bw, "HWE\t0\t%f\t%d", af, sumTot)
		sumTmp := 0
		for j := 0; j < r.nafHWE; j++ {
			sumTmp += ptr[j]
			frac := float64(sumTmp) / float64(sumTot)
			q := float64(j) / float64(r.nafHWE)
			if frac >= 0.75 {
				for nprn > 0 {
					fmt.Fprintf(bw, "\t%f", q)
					nprn--
				}
				break
			}
			if frac >= 0.5 {
				for nprn > 1 {
					fmt.Fprintf(bw, "\t%f", q)
					nprn--
				}
				continue
			}
			if frac >= 0.25 {
				for nprn > 2 {
					fmt.Fprintf(bw, "\t%f", q)
					nprn--
				}
			}
		}
		fmt.Fprintln(bw)
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
