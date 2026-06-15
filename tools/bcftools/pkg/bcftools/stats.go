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

	// EnableSamples reports whether -s/-S was supplied. Upstream only
	// produces the per-sample sections (PSC, PSI, HWE, DP genotype counts)
	// when samples are explicitly requested; without it the per-sample
	// accumulators stay empty and those sections are suppressed.
	EnableSamples bool
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

	// AF accumulators. Upstream bins each ALT allele by its frequency into
	// one of afM buckets (the first bucket, index 0, is reserved for
	// singletons, i.e. AC==1). With the default binning afM = max(101,
	// nsamples+1) and the bin index for a non-singleton allele is
	// int(af32*(afM-2))+1; with --af-bins it is bin_get_idx(af)+1.
	afBins []float64 // user-supplied bin edges (nil => default linear binning)
	afM    int       // number of AF buckets (m_af in vcfstats.c)
	afSNPs []int
	afTs   []int
	afTv   []int
	// afIndel mirrors upstream's af_repeats[2] (the "not applicable" repeat
	// class): without -E exons every indel allele lands here, and it feeds
	// both the "number of indels" (column [7]) and "not applicable" (column
	// [10]) fields of the AF report.
	afIndel []int

	// QUAL accumulators keyed by upstream's iqual = 1 + int(qual*10); the
	// special value 0 denotes a missing/negative QUAL. qualSNPs counts only
	// the first ALT (matching upstream's ts/tv-by-quality semantics).
	qualSNPs map[int]int
	qualTs   map[int]int
	qualTv   map[int]int
	qualNonS map[int]int

	// IDD: indel length distribution (length, count) — len < 0 = deletion.
	indelLen map[int]int

	// ST: substitution-type counts (REF>ALT for SNPs).
	subst map[string]int

	// DP: depth distribution, binned with the idist scheme. dpMin/Max/Step
	// and dpM (= 4 + (max-min)/step) define the bucket layout; bucket 0 is
	// the underflow (<min) bin and bucket dpM-1 the overflow (>max) bin.
	dpMin, dpMax, dpStep, dpM int
	dpSites                   map[int]int // idist bin -> sites
	dpGTs                     map[int]int // idist bin -> per-sample GTs
	dpTotalSite               int         // total sites that contributed
	dpTotalGT                 int

	// PSC / PSI: per-sample counters indexed by sampleIndex.
	pscNRefHom    []int
	pscNNonRefHom []int
	pscNHets      []int
	pscNTs        []int
	pscNTv        []int
	pscNIndels    []int
	pscNSingleton []int // sites where this sample is the lone non-ref sample
	pscNHapRef    []int
	pscNHapAlt    []int
	pscNMissing   []int
	pscDepthSum   []int
	pscDepthN     []int

	psiNInsHet []int
	psiNDelHet []int
	psiNInsHom []int
	psiNDelHom []int

	// HWE: a 2-D histogram af_hwe[iaf*nafHWE + ihet] counting, per first-ALT
	// AF bucket, how many sites had a given observed heterozygous-genotype
	// fraction. nafHWE is the number of het-fraction bins (100 upstream).
	nafHWE int
	afHWE  []int

	// USR: per-spec Ts/Tv-by-tag accumulators (-u/--user-tstv).
	userTSTV []*userTSTVAcc

	// VAF (Variant Allele Frequency) accumulators, derived from FORMAT/AD over
	// the called depth. hasFmtAD records whether the header declares FORMAT/AD;
	// upstream only allocates the per-sample VAF distributions (and prints the
	// VAF section + the IDD mean-VAF columns) when it does and -s/-S is set.
	// vafSNV/vafIndel are per-sample 21-bucket histograms (bucket = round(vaf/
	// 0.05)). idVAFn/idVAFsum accumulate, per indel length, the genotype count
	// and the VAF sum used for the IDD "mean VAF" column.
	hasFmtAD bool
	vafSNV   [][21]int
	vafIndel [][21]int
	idVAFn   map[int]int
	idVAFsum map[int]float64
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
		nafHWE:        100,
		idVAFn:        make(map[int]int),
		idVAFsum:      make(map[int]float64),
		hasFmtAD:      formatTagDeclared(metaInfo, "AD"),
	}
	// Determine the number of AF buckets, mirroring init_stats() in
	// vcfstats.c. With --af-bins the bucket count equals the number of bin
	// boundaries (the last is unused, leaving room for the singleton bin);
	// otherwise it is max(101, nsamples+1).
	if len(r.afBins) >= 2 {
		r.afM = len(r.afBins)
	} else {
		r.afBins = nil
		r.afM = 101
		if len(headerSamples)+1 > r.afM {
			r.afM = len(headerSamples) + 1
		}
	}
	r.afSNPs = make([]int, r.afM)
	r.afTs = make([]int, r.afM)
	r.afTv = make([]int, r.afM)
	r.afIndel = make([]int, r.afM)

	// Per-sample sections are produced only when -s/-S was requested. A bare
	// "-" sample request (or no explicit list) selects every header sample.
	if opts.EnableSamples {
		want := opts.Samples
		if len(want) == 0 || (len(want) == 1 && want[0] == "-") {
			want = nil
		}
		r.samples = filterSampleSet(headerSamples, want)
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
	r.afHWE = make([]int, r.afM*r.nafHWE)
	if r.hasFmtAD {
		r.vafSNV = make([][21]int, len(r.samples))
		r.vafIndel = make([][21]int, len(r.samples))
	}

	// Depth-distribution bucket layout (idist), defaulting to 0..500 step 1.
	r.dpMin, r.dpMax, r.dpStep = opts.DepthMin, opts.DepthMax, opts.DepthStep
	if r.dpStep <= 0 {
		r.dpStep = 1
	}
	if r.dpMax <= 0 {
		r.dpMax = 500
	}
	r.dpM = 4 + (r.dpMax-r.dpMin)/r.dpStep

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

// formatTagDeclared reports whether the header meta lines declare a
// FORMAT field named tag. Mirrors upstream's bcf_hdr_id2int(..,BCF_DT_ID,..)
// lookup combined with the FORMAT-line presence test that gates the VAF
// section (has_fmt_ad in vcfstats.c).
func formatTagDeclared(metaInfo []string, tag string) bool {
	prefix := "##FORMAT=<"
	for _, line := range metaInfo {
		if !strings.HasPrefix(line, prefix) || !strings.HasSuffix(line, ">") {
			continue
		}
		fields := parseStructuredMeta(line[len(prefix) : len(line)-1])
		if fields["ID"] == tag {
			return true
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

// accumulate folds one variant into the result, mirroring do_vcf_stats() in
// vcfstats.c. The record's combined variant-type bitmask (the OR of every
// ALT's type) drives the SN counters — a row that contains both a SNP and an
// indel increments both — and each section uses the per-allele AF bucket
// array computed by computeIAF.
func accumulate(r *statsResult, v *vcf.Variant) {
	r.numRecords++

	// nAllele includes the reference, matching line->n_allele.
	nAllele := len(v.Alt) + 1
	iaf := computeIAF(r, v, nAllele)

	// Combined line type: OR of every ALT's variant-type bit (upstream's
	// bcf_get_variant_types). An empty mask means a reference-only row.
	lineType := 0
	for _, alt := range v.Alt {
		lineType |= variantTypeBit(v.Ref, alt)
	}

	if lineType == 0 {
		r.numNoALTs++
	}
	if lineType&vtSNP != 0 {
		r.doSNPStats(v, iaf)
	}
	if lineType&vtINDEL != 0 {
		r.doIndelStats(v, iaf)
	}
	if lineType&vtMNP != 0 {
		r.numMNPs++
	}
	if lineType&(vtOTHER|vtBND|vtOVERLAP) != 0 {
		r.numOthers++
	}

	if nAllele > 2 {
		r.numMA++
		if onlySNPAlts(v) {
			r.numMASNP++
		}
	}

	// DP binning — INFO/DP first, falling back to per-sample DP sum.
	if dp, ok := siteDepth(v); ok {
		r.dpSites[r.dpIdx(dp)]++
		r.dpTotalSite++
	}

	// Per-sample counters (PSC/PSI/HWE).
	accumulateSamples(r, v, iaf)
}

// doSNPStats mirrors do_snp_stats(): it bumps n_snps once, then for every SNP
// ALT (subject to --1st-allele-only) records the substitution, the AF-binned
// ts/tv split, and — for the first ALT only — the quality-binned ts/tv split
// and the user-tstv accumulators.
func (r *statsResult) doSNPStats(v *vcf.Variant, iaf []int) {
	r.numSNPs++
	ref := strings.ToUpper(v.Ref)
	if len(ref) == 0 || !isACGT(ref[0]) {
		return
	}
	iqual := qualIndex(v.Qual)
	for i := 1; i < len(iaf); i++ {
		if r.opts.FirstAlleleOnly && i > 1 {
			break
		}
		alt := v.Alt[i-1]
		if variantTypeBit(v.Ref, alt)&vtSNP == 0 {
			continue
		}
		altU := strings.ToUpper(alt)
		if len(altU) == 0 || !isACGT(altU[0]) || ref[0] == altU[0] {
			continue
		}
		r.subst[string(ref[0])+">"+string(altU[0])]++
		bin := iaf[i]
		r.afSNPs[bin]++
		if transitionType(v.Ref, alt) == "ts" {
			if i == 1 {
				r.qualTs[iqual]++
				accumulateUserTSTV(r, v, true)
			}
			r.afTs[bin]++
		} else {
			if i == 1 {
				r.qualTv[iqual]++
				accumulateUserTSTV(r, v, false)
			}
			r.afTv[bin]++
		}
	}
}

// doIndelStats mirrors do_indel_stats(): it bumps n_indels once, records the
// first-ALT quality bin, and for every indel ALT (subject to
// --1st-allele-only) records the length distribution and the AF-binned indel
// count (af_repeats[2], the "not applicable" repeat class without -E exons).
func (r *statsResult) doIndelStats(v *vcf.Variant, iaf []int) {
	r.numIndels++
	r.qualNonS[qualIndex(v.Qual)]++
	for i := 1; i < len(iaf); i++ {
		if r.opts.FirstAlleleOnly && i > 1 {
			break
		}
		alt := v.Alt[i-1]
		if variantTypeBit(v.Ref, alt)&vtINDEL == 0 {
			continue
		}
		r.indelLen[len(alt)-len(v.Ref)]++
		r.afIndel[iaf[i]]++
	}
}

// qualIndex maps a QUAL value to upstream's iqual bucket: 0 for missing or
// negative quality, otherwise 1 + int(qual*10) so the bucket width is 0.1.
func qualIndex(qual float64) int {
	if math.IsNaN(qual) || qual < 0 {
		return 0
	}
	return 1 + int(qual*10)
}

// onlySNPAlts reports whether every ALT of a multi-allelic record is a SNP,
// matching upstream's "multiallelic SNP sites" criterion.
func onlySNPAlts(v *vcf.Variant) bool {
	for _, alt := range v.Alt {
		if variantTypeBit(v.Ref, alt)&vtSNP == 0 {
			return false
		}
	}
	return len(v.Alt) > 0
}

// computeIAF mirrors init_iaf(): it returns the per-allele AF bucket index
// (length nAllele; element 0 is the reference and is always bucket 0). With
// --af-tag the frequency is read from the named INFO tag; otherwise it is
// derived from the allele counts (INFO/AC+AN or the genotypes, exactly like
// bcf_calc_ac). Singletons (AC==1) fall into bucket 0, the bin reserved for
// them; a non-singleton allele frequency af maps to bucket
// int(float32(af)*(afM-2))+1 (or bin_get_idx(af)+1 with --af-bins).
func computeIAF(r *statsResult, v *vcf.Variant, nAllele int) []int {
	iaf := make([]int, nAllele)
	if r.opts.AFTag != "" {
		if af, ok := infoFloatArray(v, r.opts.AFTag); ok && len(af) == nAllele-1 {
			for i := 1; i < nAllele; i++ {
				iaf[i] = r.afBucket(clamp01(af[i-1])) + 1
			}
		}
		return iaf
	}
	ac := computeACWithRef(v, nAllele)
	an := 0
	for _, c := range ac {
		an += c
	}
	if an == 0 {
		return iaf
	}
	for i := 1; i < nAllele; i++ {
		switch {
		case ac[i] == 1:
			iaf[i] = 0 // singletons into the first bucket
		default:
			af := clamp01(float64(float32(ac[i]) / float32(an)))
			iaf[i] = r.afBucket(af) + 1
		}
	}
	return iaf
}

// afBucket maps a clamped allele frequency to its zero-based bucket (before
// the +1 singleton offset). With user-supplied --af-bins it is the
// half-open-interval index; otherwise int(float32(af)*(afM-2)), matching
// vcfstats.c's `af*(m_af-2)` truncation in float32 arithmetic.
func (r *statsResult) afBucket(af float64) int {
	if len(r.afBins) >= 2 {
		return binGetIdx(r.afBins, af)
	}
	return int(float32(af) * float32(r.afM-2))
}

// clamp01 clamps x to the closed unit interval.
func clamp01(x float64) float64 {
	if x < 0 {
		return 0
	}
	if x > 1 {
		return 1
	}
	return x
}

// binGetIdx returns the index of the half-open [bins[i],bins[i+1]) interval
// containing af, mirroring htslib's bin_get_idx. Values at or below the first
// edge land in bin 0; values at or above the last edge land in the last bin.
func binGetIdx(bins []float64, af float64) int {
	if len(bins) < 2 {
		return 0
	}
	if af <= bins[0] {
		return 0
	}
	if af >= bins[len(bins)-1] {
		return len(bins) - 2
	}
	for i := 0; i+1 < len(bins); i++ {
		if af >= bins[i] && af < bins[i+1] {
			return i
		}
	}
	return len(bins) - 2
}

// infoFloatArray parses INFO tag as a comma-separated float array.
func infoFloatArray(v *vcf.Variant, tag string) ([]float64, bool) {
	raw, ok := v.Info[tag]
	if !ok || raw == "" || raw == "." {
		return nil, false
	}
	parts := strings.Split(raw, ",")
	out := make([]float64, len(parts))
	for i, p := range parts {
		f, err := strconv.ParseFloat(p, 64)
		if err != nil {
			return nil, false
		}
		out[i] = f
	}
	return out, true
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

// siteDepth returns INFO/DP if present. Upstream's dp_sites distribution is fed
// exclusively from a single-valued INFO/DP (vcfstats.c line ~1307); it does not
// fall back to summing FORMAT/DP, so neither do we.
func siteDepth(v *vcf.Variant) (int, bool) {
	if raw, ok := v.Info["DP"]; ok && raw != "" && raw != "." {
		if n, err := strconv.Atoi(raw); err == nil {
			return n, true
		}
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

// dpIdx returns the idist bucket index for depth dp: 0 for dp<min (the
// underflow bin), dpM-1 for dp>max (the overflow bin), and 1+(dp-min)/step
// otherwise. It mirrors htslib's idist().
func (r *statsResult) dpIdx(dp int) int {
	if dp < r.dpMin {
		return 0
	}
	if dp > r.dpMax {
		return r.dpM - 1
	}
	return 1 + (dp-r.dpMin)/r.dpStep
}

// dpLabel renders the idist bucket label for index i: "<min", ">max", or the
// bucket's depth value (i-1+min), mirroring print_stats() / idist_i2bin().
func (r *statsResult) dpLabel(i int) string {
	if i == 0 {
		return fmt.Sprintf("<%d", r.dpMin)
	}
	if i == r.dpM-1 {
		return fmt.Sprintf(">%d", r.dpMax)
	}
	return strconv.Itoa(i - 1 + r.dpMin)
}

// Genotype-type codes mirroring htslib's bcf_gt_type return values, as used
// by do_sample_stats() in vcfstats.c.
const (
	gtHomRR = iota // 0/0
	gtHetRA        // 0/x
	gtHetAA        // x/y, x!=y, both non-ref
	gtHomAA        // x/x, x non-ref
	gtHaplR        // single ref allele
	gtHaplA        // single non-ref allele
	gtUnkn         // missing
)

// gtType classifies a parsed genotype into one of the gt* codes and returns
// the two allele indices ial/jal (for a haploid genotype jal==ial). It mirrors
// bcf_gt_type: a missing allele in either position yields gtUnkn.
func gtType(als [2]int, kind genotypeKind) (typ, ial, jal int) {
	if kind == gtMissing {
		return gtUnkn, 0, 0
	}
	ial, jal = als[0], als[1]
	if kind == gtHemi {
		if ial == 0 {
			return gtHaplR, ial, jal
		}
		return gtHaplA, ial, jal
	}
	switch {
	case ial == 0 && jal == 0:
		return gtHomRR, ial, jal
	case ial == 0 || jal == 0:
		return gtHetRA, ial, jal
	case ial == jal:
		return gtHomAA, ial, jal
	default:
		return gtHetAA, ial, jal
	}
}

// accumulateSamples updates the PSC / PSI / DP-by-genotype / HWE per-sample
// counters for v, mirroring do_sample_stats() and sample_gt_stats().
func accumulateSamples(r *statsResult, v *vcf.Variant, iaf []int) {
	if len(r.samples) == 0 {
		return
	}
	var nrefTot, nhetTot, naltTot int
	nNonRef := 0
	iNonRef := 0
	for i := range v.Samples {
		idx, ok := r.sampleIndex[v.Samples[i].Name]
		if !ok {
			continue
		}
		// Track depth contributions (DP-by-genotype + PSC average depth).
		// Upstream's calc_sample_depth prefers FORMAT/DP and falls back to the
		// sum of FORMAT/AD when DP is absent or missing.
		dp, dpOK := calcSampleDepth(v, i)
		if dpOK && dp > 0 {
			r.pscDepthSum[idx] += dp
			r.pscDepthN[idx]++
			r.dpGTs[r.dpIdx(dp)]++
			r.dpTotalGT++
		}

		als, kind := parseGenotypeAlleles(v, i)
		typ, ial, jal := gtType(als, kind)

		// VAF distributions (FORMAT/AD over the called depth). Mirrors the
		// update_vaf / update_dvaf calls in do_sample_stats: only when DP>0 and
		// FORMAT/AD is present. For a missing GT, AD[ial] for the alt found by
		// get_ad is used; otherwise the alt allele(s) of the genotype are used.
		if r.hasFmtAD && dpOK && dp > 0 {
			r.accumulateVAF(v, i, idx, als, kind, typ, ial, jal, dp)
		}
		if typ == gtUnkn {
			r.pscNMissing[idx]++
			continue
		}

		varType := 0
		if ial > 0 {
			varType |= altVariantType(v, ial)
		}
		if jal > 0 {
			varType |= altVariantType(v, jal)
		}

		if typ == gtHaplR || typ == gtHaplA {
			if varType&vtINDEL != 0 {
				// frame-shift stats require -E exons; not tracked here.
			}
			if typ == gtHaplR {
				r.pscNHapRef[idx]++
			} else {
				r.pscNHapAlt[idx]++
			}
			continue
		}

		if typ != gtHomRR {
			nNonRef++
			iNonRef = idx
		}
		switch typ {
		case gtHomRR:
			nrefTot++
		case gtHetRA:
			nhetTot++
		case gtHetAA, gtHomAA:
			naltTot++
		}

		if varType&vtSNP != 0 || varType == 0 { // count ALT=. as SNP
			switch typ {
			case gtHetRA, gtHetAA:
				r.pscNHets[idx]++
			case gtHomRR:
				r.pscNRefHom[idx]++
			case gtHomAA:
				r.pscNNonRefHom[idx]++
			}
			if typ != gtHomRR && altVariantType(v, ial)&vtSNP != 0 {
				r.pscBumpTsTv(idx, v, ial)
			}
			if typ != gtHomRR && ial != jal && altVariantType(v, jal)&vtSNP != 0 {
				r.pscBumpTsTv(idx, v, jal)
			}
		}
		if varType&vtINDEL != 0 && typ != gtHomRR {
			r.pscNIndels[idx]++
			switch typ {
			case gtHetRA, gtHetAA:
				isIns, isDel := false, false
				if altVariantType(v, ial)&vtINDEL != 0 {
					if indelAlleleLen(v, ial) < 0 {
						isDel = true
					} else {
						isIns = true
					}
				}
				if altVariantType(v, jal)&vtINDEL != 0 {
					if indelAlleleLen(v, jal) < 0 {
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
				if indelAlleleLen(v, ial) < 0 {
					r.psiNDelHom[idx]++
				} else {
					r.psiNInsHom[idx]++
				}
			}
		}
	}
	if nNonRef == 1 {
		r.pscNSingleton[iNonRef]++
	}

	// HWE het-fraction histogram, keyed by the first-ALT AF bucket.
	if len(v.Alt) > 0 && (nrefTot != 0 || nhetTot != 0 || naltTot != 0) {
		hetFrac := float32(nhetTot) / float32(nrefTot+nhetTot+naltTot)
		ihet := int(hetFrac * float32(r.nafHWE-1))
		bucket := 0
		if len(iaf) > 1 {
			bucket = iaf[1]
		}
		r.afHWE[ihet+bucket*r.nafHWE]++
	}
}

// calcSampleDepth returns the called depth of sample i, mirroring
// calc_sample_depth(): FORMAT/DP if present and not missing, otherwise the sum
// of the non-missing FORMAT/AD values. Returns false when neither is available.
func calcSampleDepth(v *vcf.Variant, i int) (int, bool) {
	if raw, ok := v.Samples[i].Data["DP"]; ok && raw != "" && raw != "." {
		if n, err := strconv.Atoi(raw); err == nil {
			return n, true
		}
		return 0, false
	}
	ad, ok := sampleADInts(v, i)
	if !ok {
		return 0, false
	}
	sum, has := 0, false
	for _, a := range ad {
		if a < 0 {
			continue // missing entry
		}
		sum += a
		has = true
	}
	if !has {
		return 0, false
	}
	return sum, true
}

// sampleADInts parses sample i's FORMAT/AD into an integer slice. A "." entry
// (missing) is rendered as -1. Returns false when AD is absent or empty.
func sampleADInts(v *vcf.Variant, i int) ([]int, bool) {
	raw, ok := v.Samples[i].Data["AD"]
	if !ok || raw == "" || raw == "." {
		return nil, false
	}
	parts := strings.Split(raw, ",")
	out := make([]int, len(parts))
	for k, p := range parts {
		if p == "." || p == "" {
			out[k] = -1
			continue
		}
		n, err := strconv.Atoi(p)
		if err != nil {
			return nil, false
		}
		out[k] = n
	}
	return out, true
}

// accumulateVAF folds sample i's FORMAT/AD into the per-sample VAF histograms
// and the per-length indel mean-VAF accumulators, mirroring the update_vaf /
// update_dvaf calls in do_sample_stats. For a missing genotype the dominant ALT
// (largest AD) is used; otherwise the genotype's ALT allele(s) are used, with
// the second allele counted only when its AD differs from the first.
func (r *statsResult) accumulateVAF(v *vcf.Variant, i, idx int, als [2]int, kind genotypeKind, typ, ial, jal, dp int) {
	ad, ok := sampleADInts(v, i)
	if !ok {
		return
	}
	adAt := func(a int) int {
		if a < 0 || a >= len(ad) || ad[a] < 0 {
			return 0
		}
		return ad[a]
	}
	var iad, jad int
	var iAllele, jAllele int
	if typ == gtUnkn {
		// get_ad: pick the alt allele with the maximum AD.
		best, bestAl := 0, 0
		for a := 1; a < len(ad) && a < len(v.Alt)+1; a++ {
			if ad[a] < 0 {
				continue
			}
			if ad[a] > best {
				best, bestAl = ad[a], a
			}
		}
		iad, iAllele = best, bestAl
	} else {
		if ial != 0 {
			iad = adAt(ial)
		}
		iAllele = ial
		if jal != 0 {
			jad = adAt(jal)
		}
		jAllele = jal
	}
	if iad != 0 {
		r.recordVAF(v, idx, iAllele, float64(iad)/float64(dp))
	}
	if jad != 0 && iad != jad {
		r.recordVAF(v, idx, jAllele, float64(jad)/float64(dp))
	}
}

// recordVAF buckets one VAF observation: a SNP allele goes into the sample's
// SNV histogram, anything else into the indel histogram (vaf2bin =
// round(vaf/0.05), restricted to vaf in [0,1]); indel alleles additionally feed
// the per-length mean-VAF accumulators keyed by the allele's length change.
func (r *statsResult) recordVAF(v *vcf.Variant, idx, allele int, vaf float64) {
	if vaf >= 0 && vaf <= 1 {
		bin := int(roundHalfEven(vaf / 0.05))
		if bin >= 0 && bin < 21 {
			if altVariantType(v, allele)&vtSNP != 0 {
				r.vafSNV[idx][bin]++
			} else {
				r.vafIndel[idx][bin]++
			}
		}
	}
	if altVariantType(v, allele)&vtINDEL != 0 {
		d := indelAlleleLen(v, allele)
		if d < -60 {
			d = -60
		} else if d > 60 {
			d = 60
		}
		r.idVAFn[d]++
		r.idVAFsum[d] += vaf
	}
}

// roundHalfEven rounds x to the nearest integer using round-half-to-even,
// matching C's nearbyintf under the default (round-to-nearest) FP mode.
func roundHalfEven(x float64) float64 {
	return math.RoundToEven(x)
}

// pscBumpTsTv increments the per-sample transition or transversion counter for
// SNP allele a (1-based) at site v.
func (r *statsResult) pscBumpTsTv(idx int, v *vcf.Variant, a int) {
	if transitionType(v.Ref, v.Alt[a-1]) == "ts" {
		r.pscNTs[idx]++
	} else {
		r.pscNTv[idx]++
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
	if err := writeUSR(bw, r); err != nil {
		return err
	}
	if err := writeDP(bw, r); err != nil {
		return err
	}
	// The per-sample sections are emitted only when samples were requested,
	// matching upstream's `if (n_smpl)` guard.
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
		if err := writeVAF(bw, r); err != nil {
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

// writeTSTV emits the TSTV (transitions/transversions) summary row. The
// overall ts/tv counts are the sums of the AF-binned transition and
// transversion tallies across every bucket (including the singleton bucket 0,
// read before the SiS fold), mirroring print_stats()'s loop over m_af. The
// 1st-ALT counts reuse the quality-binned first-ALT ts/tv accumulators
// (ts_alt1/tv_alt1 upstream), which are incremented at the same place as the
// QUAL section's ts/tv. Both ratios use upstream's `tv ? ts/tv : 0` rule with
// %.2f formatting.
func writeTSTV(bw *bufio.Writer, r *statsResult) error {
	fmt.Fprintln(bw, "# TSTV, transitions/transversions")
	fmt.Fprintln(bw, "#   - transitions, see https://en.wikipedia.org/wiki/Transition_(genetics)")
	fmt.Fprintln(bw, "#   - transversions, see https://en.wikipedia.org/wiki/Transversion")
	fmt.Fprintln(bw, "# TSTV\t[2]id\t[3]ts\t[4]tv\t[5]ts/tv\t[6]ts (1st ALT)\t[7]tv (1st ALT)\t[8]ts/tv (1st ALT)")
	ts, tv := 0, 0
	for i := 0; i < r.afM; i++ {
		ts += r.afTs[i]
		tv += r.afTv[i]
	}
	tsAlt1, tvAlt1 := 0, 0
	for _, c := range r.qualTs {
		tsAlt1 += c
	}
	for _, c := range r.qualTv {
		tvAlt1 += c
	}
	fmt.Fprintf(bw, "TSTV\t0\t%d\t%d\t%.2f\t%d\t%d\t%.2f\n",
		ts, tv, tstvRatio(ts, tv), tsAlt1, tvAlt1, tstvRatio(tsAlt1, tvAlt1))
	return nil
}

// tstvRatio returns ts/tv as a float32-precision quotient, or 0 when tv==0,
// matching upstream's `tv ? (float)ts/tv : 0` expression.
func tstvRatio(ts, tv int) float64 {
	if tv == 0 {
		return 0
	}
	return float64(float32(ts) / float32(tv))
}

// writeSiS emits the deprecated SiS (Singleton stats) row. Upstream reports the
// contents of AF bucket 0 — the bucket reserved for singleton (AC==1) alleles —
// as: allele count (always 1), number of SNPs, transitions, transversions, and
// the indel count split into repeat-consistent/inconsistent/not-applicable
// classes. We only track the "not applicable" repeat class (af_repeats[2],
// stored in afIndel), so the number-of-indels and not-applicable columns both
// read afIndel[0] and the two repeat columns are 0. This must run before
// writeAF folds bucket 0 into bucket 1.
func writeSiS(bw *bufio.Writer, r *statsResult) error {
	fmt.Fprintln(bw, "# SiS, Singleton stats:")
	fmt.Fprintln(bw, "#   - allele count, i.e. the number of singleton genotypes (AC=1)")
	fmt.Fprintln(bw, "#   - number of transitions, see above")
	fmt.Fprintln(bw, "#   - number of transversions, see above")
	fmt.Fprintln(bw, "#   - repeat-consistent, inconsistent and n/a: experimental and useless stats [DEPRECATED]")
	fmt.Fprintln(bw, "# SiS\t[2]id\t[3]allele count\t[4]number of SNPs\t[5]number of transitions\t[6]number of transversions\t[7]number of indels\t[8]repeat-consistent\t[9]repeat-inconsistent\t[10]not applicable")
	snps, ts, tv, indel := 0, 0, 0, 0
	if r.afM > 0 {
		snps, ts, tv, indel = r.afSNPs[0], r.afTs[0], r.afTv[0], r.afIndel[0]
	}
	fmt.Fprintf(bw, "SiS\t0\t1\t%d\t%d\t%d\t%d\t0\t0\t%d\n", snps, ts, tv, indel, indel)
	return nil
}

// writeAF emits the per-AF-bucket statistics. Upstream folds the singleton
// bucket (index 0) into bucket 1, prints only non-empty buckets, and renders
// the bucket frequency as (i-1)/(afM-1) — or the midpoint of the user bin with
// --af-bins. The "number of indels" and "not applicable" columns both read the
// af_repeats[2] count (afIndel) because without -E exons every indel lands in
// the repeat-n/a class.
func writeAF(bw *bufio.Writer, r *statsResult) error {
	fmt.Fprintln(bw, "# AF, Stats by non-reference allele frequency:")
	fmt.Fprintln(bw, "# AF\t[2]id\t[3]allele frequency\t[4]number of SNPs\t[5]number of transitions\t[6]number of transversions\t[7]number of indels\t[8]repeat-consistent\t[9]repeat-inconsistent\t[10]not applicable")
	// Fold the singleton bucket into bucket 1, mirroring print_stats().
	if r.afM > 1 {
		r.afSNPs[1] += r.afSNPs[0]
		r.afTs[1] += r.afTs[0]
		r.afTv[1] += r.afTv[0]
		r.afIndel[1] += r.afIndel[0]
	}
	for i := 1; i < r.afM; i++ {
		if r.afSNPs[i]+r.afTs[i]+r.afTv[i]+r.afIndel[i] == 0 {
			continue
		}
		af := r.afFrequency(i)
		fmt.Fprintf(bw, "AF\t0\t%f\t%d\t%d\t%d\t%d\t0\t0\t%d\n",
			af, r.afSNPs[i], r.afTs[i], r.afTv[i], r.afIndel[i], r.afIndel[i])
	}
	return nil
}

// afFrequency returns the printed frequency for AF bucket i, matching
// print_stats(): the bin midpoint with --af-bins, else (i-1)/(afM-1).
func (r *statsResult) afFrequency(i int) float64 {
	if len(r.afBins) >= 2 {
		hi := binGetValue(r.afBins, i)
		lo := binGetValue(r.afBins, i-1)
		return (hi + lo) * 0.5
	}
	return float64(i-1) / float64(r.afM-1)
}

// binGetValue returns the i-th bin boundary, clamped to the available edges,
// mirroring htslib's bin_get_value.
func binGetValue(bins []float64, i int) float64 {
	if i < 0 {
		i = 0
	}
	if i >= len(bins) {
		i = len(bins) - 1
	}
	return bins[i]
}

// writeQUAL emits the quality distribution. Each bucket's quality is rendered
// as 0.1*(iqual-1) with one decimal — so integer QUAL=30 prints as "30.0" —
// and a missing/negative QUAL (iqual 0) prints as ".". The SNP column is the
// sum of first-ALT transitions and transversions; indels add their own column.
func writeQUAL(bw *bufio.Writer, r *statsResult) error {
	fmt.Fprintln(bw, "# QUAL, Stats by quality")
	fmt.Fprintln(bw, "# QUAL\t[2]id\t[3]Quality\t[4]number of SNPs\t[5]number of transitions (1st ALT)\t[6]number of transversions (1st ALT)\t[7]number of indels")
	keys := unionMapKeys(r.qualTs, r.qualTv, r.qualNonS)
	sort.Ints(keys)
	for _, k := range keys {
		nts, ntv, nin := r.qualTs[k], r.qualTv[k], r.qualNonS[k]
		if nts+ntv+nin == 0 {
			continue
		}
		fmt.Fprint(bw, "QUAL\t0\t")
		if k == 0 {
			fmt.Fprint(bw, ".")
		} else {
			fmt.Fprintf(bw, "%.1f", 0.1*float64(k-1))
		}
		fmt.Fprintf(bw, "\t%d\t%d\t%d\t%d\n", nts+ntv, nts, ntv, nin)
	}
	return nil
}

// writeIDD emits the indel-length distribution. Upstream prints deletions
// (negative lengths) first in descending magnitude, then insertions, and emits
// "0\t." for the unset genotype-count and mean-VAF trailing columns when no
// FORMAT/AD VAF data was collected.
func writeIDD(bw *bufio.Writer, r *statsResult) error {
	fmt.Fprintln(bw, "# IDD, InDel distribution:")
	fmt.Fprintln(bw, "# IDD\t[2]id\t[3]length (deletions negative)\t[4]number of sites\t[5]number of genotypes\t[6]mean VAF")
	keys := make([]int, 0, len(r.indelLen))
	for k := range r.indelLen {
		keys = append(keys, k)
	}
	sort.Ints(keys)
	for _, k := range keys {
		fmt.Fprintf(bw, "IDD\t0\t%d\t%d\t", k, r.indelLen[k])
		if r.hasFmtAD && r.idVAFn[k] > 0 {
			fmt.Fprintf(bw, "%d\t%.2f\n", r.idVAFn[k], r.idVAFsum[k]/float64(r.idVAFn[k]))
		} else {
			fmt.Fprint(bw, "0\t.\n")
		}
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

// writeDP emits the depth distribution over the idist buckets in bucket
// order, skipping buckets where both the genotype and site counts are zero.
// The first/last buckets carry the dynamic "<min"/">max" labels; fractions are
// percentages of the respective genotype and site totals.
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
	for i := 0; i < r.dpM; i++ {
		gt, site := r.dpGTs[i], r.dpSites[i]
		if gt == 0 && site == 0 {
			continue
		}
		var fracGT, fracSite float64
		if r.dpTotalGT > 0 {
			fracGT = 100.0 * float64(gt) / float64(r.dpTotalGT)
		}
		if r.dpTotalSite > 0 {
			fracSite = 100.0 * float64(site) / float64(r.dpTotalSite)
		}
		fmt.Fprintf(bw, "DP\t0\t%s\t%d\t%f\t%d\t%f\n",
			r.dpLabel(i), gt, fracGT, site, fracSite)
	}
	return nil
}

// writePSC emits the per-sample counts. The average depth uses upstream's
// %.1f formatting (float32 accumulation of the mean), and the trailing
// singleton/haploid/missing columns are now tracked rather than stubbed.
func writePSC(bw *bufio.Writer, r *statsResult) error {
	fmt.Fprintln(bw, "# PSC, Per-sample counts. Note that the ref/het/hom counts include only SNPs, for indels see PSI. The rest include both SNPs and indels.")
	fmt.Fprintln(bw, "# PSC\t[2]id\t[3]sample\t[4]nRefHom\t[5]nNonRefHom\t[6]nHets\t[7]nTransitions\t[8]nTransversions\t[9]nIndels\t[10]average depth\t[11]nSingletons\t[12]nHapRef\t[13]nHapAlt\t[14]nMissing")
	for i, name := range r.samples {
		var avg float32
		if r.pscDepthN[i] > 0 {
			avg = float32(r.pscDepthSum[i]) / float32(r.pscDepthN[i])
		}
		fmt.Fprintf(bw, "PSC\t0\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.1f\t%d\t%d\t%d\t%d\n",
			name,
			r.pscNRefHom[i], r.pscNNonRefHom[i], r.pscNHets[i],
			r.pscNTs[i], r.pscNTv[i], r.pscNIndels[i], avg,
			r.pscNSingleton[i], r.pscNHapRef[i], r.pscNHapAlt[i], r.pscNMissing[i])
	}
	return nil
}

// writePSI emits the per-sample indel statistics. The in-frame/out-frame
// frame-shift columns require -E exons (not supported here) and so stay zero;
// the nInsHets/nDelHets/nInsAltHoms/nDelAltHoms counts are tracked.
func writePSI(bw *bufio.Writer, r *statsResult) error {
	fmt.Fprintln(bw, "# PSI, Per-Sample Indels. Note that alt-het genotypes with both ins and del allele are counted twice, in both nInsHets and nDelHets.")
	fmt.Fprintln(bw, "# PSI\t[2]id\t[3]sample\t[4]in-frame\t[5]out-frame\t[6]not applicable\t[7]out/(in+out) ratio\t[8]nInsHets\t[9]nDelHets\t[10]nInsAltHoms\t[11]nDelAltHoms")
	for i, name := range r.samples {
		fmt.Fprintf(bw, "PSI\t0\t%s\t0\t0\t0\t0.00\t%d\t%d\t%d\t%d\n",
			name, r.psiNInsHet[i], r.psiNDelHet[i], r.psiNInsHom[i], r.psiNDelHom[i])
	}
	return nil
}

// writeHWE emits the Hardy-Weinberg section: for each first-ALT AF bucket with
// observations it reports the count and the 25th/median/75th percentiles of
// the observed heterozygous-genotype fraction, computed from the afHWE
// histogram exactly as print_stats() does. The singleton bucket is folded into
// bucket 1 first.
func writeHWE(bw *bufio.Writer, r *statsResult) error {
	fmt.Fprintln(bw, "# HWE")
	fmt.Fprintln(bw, "# HWE\t[2]id\t[3]1st ALT allele frequency\t[4]Number of observations\t[5]25th percentile\t[6]median\t[7]75th percentile")
	// Fold singletons (bucket 0) into bucket 1.
	for j := 0; j < r.nafHWE; j++ {
		r.afHWE[j+r.nafHWE] += r.afHWE[j]
	}
	for i := 1; i < r.afM; i++ {
		ptr := r.afHWE[i*r.nafHWE : (i+1)*r.nafHWE]
		sumTot := 0
		for _, c := range ptr {
			sumTot += c
		}
		if sumTot == 0 {
			continue
		}
		af := r.afFrequency(i)
		p25, p50, p75 := hwePercentiles(ptr, sumTot, r.nafHWE)
		fmt.Fprintf(bw, "HWE\t0\t%f\t%d\t%f\t%f\t%f\n", af, sumTot, p25, p50, p75)
	}
	return nil
}

// writeVAF emits the VAF section: per-sample SNV and indel VAF distributions
// over 21 buckets (round(vaf/0.05) for vaf in [0,1]). Upstream prints this only
// when the header declares FORMAT/AD and samples were requested; the buckets
// are comma-joined exactly as upstream does.
func writeVAF(bw *bufio.Writer, r *statsResult) error {
	if !r.hasFmtAD {
		return nil
	}
	fmt.Fprintln(bw, "# VAF, Variant Allele Frequency determined as fraction of alternate reads in FORMAT/AD")
	fmt.Fprintln(bw, "# VAF\t[2]id\t[3]sample\t[4]SNV VAF distribution\t[5]indel VAF distribution")
	for i, name := range r.samples {
		fmt.Fprintf(bw, "VAF\t0\t%s\t%s\t%s\n", name, joinInts(r.vafSNV[i][:]), joinInts(r.vafIndel[i][:]))
	}
	return nil
}

// joinInts renders a slice of ints as a comma-separated string.
func joinInts(xs []int) string {
	var b strings.Builder
	for i, x := range xs {
		if i > 0 {
			b.WriteByte(',')
		}
		b.WriteString(strconv.Itoa(x))
	}
	return b.String()
}

// hwePercentiles walks the cumulative het-fraction histogram and returns the
// 25th, 50th and 75th percentile fractions (j/nafHWE at the first bin whose
// cumulative share crosses each threshold), mirroring the percentile loop in
// print_stats().
func hwePercentiles(ptr []int, sumTot, nafHWE int) (p25, p50, p75 float64) {
	nprn := 3
	sumTmp := 0
	for j := 0; j < nafHWE; j++ {
		sumTmp += ptr[j]
		frac := float32(sumTmp) / float32(sumTot)
		val := float64(j) / float64(nafHWE)
		if frac >= 0.75 {
			for nprn > 0 {
				assignPercentile(&p25, &p50, &p75, nprn, val)
				nprn--
			}
			break
		}
		if frac >= 0.5 {
			for nprn > 1 {
				assignPercentile(&p25, &p50, &p75, nprn, val)
				nprn--
			}
			continue
		}
		if frac >= 0.25 {
			for nprn > 2 {
				assignPercentile(&p25, &p50, &p75, nprn, val)
				nprn--
			}
		}
	}
	return p25, p50, p75
}

// assignPercentile writes val into the percentile slot selected by the
// remaining-fields counter nprn (3 => 25th, 2 => median, 1 => 75th), matching
// the order in which print_stats() fills the three columns.
func assignPercentile(p25, p50, p75 *float64, nprn int, val float64) {
	switch nprn {
	case 3:
		*p25 = val
	case 2:
		*p50 = val
	case 1:
		*p75 = val
	}
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
