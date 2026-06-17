package vcftools

import (
	"fmt"
	"io"
	"math"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// statistics holds accumulated statistics for variants
type statistics struct {
	header *vcf.Header

	// Per-site statistics
	siteFrequencies []siteFreqStat
	siteDepths      []siteDepthStat
	siteQualities   []siteQualityStat
	siteMissing     []siteMissingStat
	siteHWE         []siteHWEStat
	sitePiValues    []sitePiStat
	windowPiSites   []windowPiSite

	// Per-individual statistics
	indvMissing map[string]*indvMissingStat
	indvHet     map[string]*indvHetStat
	indvDepth   map[string]*indvDepthStat

	// Ts/Tv statistics
	transitions   int
	transversions int
	// tsTvModelCounts holds per-substitution-model counts in the fixed
	// canonical order ("AC", "AG", "AT", "CG", "CT", "GT") — see upstream
	// output_TsTv_summary() in variant_file_output.cpp. AG and CT are
	// transitions, the rest are transversions.
	tsTvModelCounts [6]int
	// tsTvBins holds the per-chromosome --TsTv bin tallies. Upstream keys
	// Ts/Tv counts on (CHROM, bin index) and emits chromosomes in first-seen
	// order (output_TsTv, variant_file_output.cpp:2962). tsTvChromOrder
	// preserves that order.
	tsTvBins       map[string]map[int]*tsTvBinStat
	tsTvChromOrder []string
	tsTvByCount    map[int]*tsTvCountStat
	tsTvByQual     []tsTvQualStat

	// Phase 2: Population genetics statistics
	windowPiValues []windowPiStat
	tajimaDValues  []tajimaDStat
	// SNP-density bins per chromosome (variant counts indexed by bin). The
	// chromosome iteration order matches first-seen order. snpDenPrevPos /
	// snpDenPrevChrom dedup adjacent records at the same (CHROM, POS), as
	// upstream does.
	snpDensityBins  map[string][]int
	snpDensityChrs  []string
	snpDenPrevPos   int
	snpDenPrevChrom string
	fstValues       []fstStat
	// FILTER-summary accumulators keyed by the full FILTER-field string
	// (e.g. "PASS", ".", "q10;s50"): per-FILTER site count and Ts/Tv tallies.
	filterCounts   map[string]int
	filterTs       map[string]int
	filterTv       map[string]int
	singletonSites []singletonStat

	// Misc
	indelLenHist  map[int]int
	indelLenTotal int
	genoDepths    []genoDepthSite
	tajimaDSites  []tajimaDSite

	// Weir & Cockerham 1984 Fst accumulator (populated when --weir-fst-pop is
	// specified). Nil when Fst calculation was not requested.
	weirFst *weirFstAccumulator
}

// genoDepthSite holds the per-individual read depth at one site.
type genoDepthSite struct {
	chrom  string
	pos    int
	depths []int // one entry per sample, -1 if FORMAT/DP absent
}

// tajimaDSite holds the data needed to compute Tajima's D for one SNP.
type tajimaDSite struct {
	chrom string
	pos   int
	pi    float64 // per-site nucleotide diversity
	nChr  int     // number of non-missing chromosomes
}

type siteFreqStat struct {
	chrom    string
	pos      int
	nAlleles int
	nChr     int
	// alleles holds every allele in upstream index order (REF=0,
	// ALT[0]=1, ...). counts is the parallel per-allele chromosome count.
	// Multi-allelic and monomorphic (ALT=".") sites are represented here
	// in full, matching upstream's get_allele_counts over all N_alleles.
	alleles []string
	counts  []int
	// aaIdx is the ancestral-allele index when --derived was requested
	// (the matched allele is emitted first). 0 otherwise.
	aaIdx int
}

type siteDepthStat struct {
	chrom string
	pos   int
	// sum and sumsq are the sum and sum-of-squares of per-individual
	// read depths at the site; n is the number of non-missing depths.
	// Upstream computes MEAN_DEPTH = sum/n and VAR_DEPTH =
	// ((sumsq/n) - mean^2) * n/(n-1) (variant_file_output.cpp:3454-3458).
	sum   int
	sumsq int
	n     int
}

type siteQualityStat struct {
	chrom string
	pos   int
	qual  float64
}

type siteMissingStat struct {
	chrom string
	pos   int
	// nData is N_DATA (total non-filtered chromosomes at the site, 2 per
	// diploid call, 1 per haploid), nMiss is the number of missing
	// chromosomes. F_MISS = nMiss/nData. Mirrors upstream
	// output_site_missingness (variant_file_output.cpp:847-918).
	nData int
	nMiss int
}

type siteHWEStat struct {
	chrom   string
	pos     int
	obsHom1 int
	obsHet  int
	obsHom2 int
	expHom1 float64
	expHet  float64
	expHom2 float64
	chiSq   float64
	// pHWE is the two-sided exact-test p-value; pLo / pHi are the one-sided
	// het-deficit / het-excess p-values (upstream P_HET_DEFICIT /
	// P_HET_EXCESS).
	pHWE float64
	pLo  float64
	pHi  float64
}

type sitePiStat struct {
	chrom string
	pos   int
	pi    float64
}

// windowPiSite holds the per-site quantities --window-pi needs to bin: the
// number of actual pairwise mismatches at the site, the number of non-missing
// chromosomes (for the per-site pairwise-comparison count), and whether the
// site is polymorphic w.r.t. the reference (allele_counts[0] < N_chr). Only
// fully-diploid, non-fixed sites are recorded, mirroring upstream
// output_windowed_nucleotide_diversity (variant_file_output.cpp:4122-4178).
type windowPiSite struct {
	chrom         string
	pos           int
	nMismatches   uint64
	nChr          int
	isPolymorphic bool
}

type indvMissingStat struct {
	name string
	// nMissing counts genotypes whose call is missing (first allele "."),
	// among those that passed genotype-level filtering. Maps to upstream's
	// indv_N_missing (variant_file_output.cpp:830).
	nMissing int
	// nTotal counts genotypes that passed genotype-level filtering and were
	// parsed. Maps to upstream's indv_N_tot (variant_file_output.cpp:831)
	// and is emitted as the N_DATA column.
	nTotal int
	// nGenoFiltered counts genotypes excluded by a genotype-level filter
	// (--minDP/--maxDP/--minGQ/--remove-filtered-geno*). Maps to upstream's
	// indv_N_geno_filtered (variant_file_output.cpp:823) and is emitted as
	// the N_GENOTYPES_FILTERED column.
	nGenoFiltered int
}

type tsTvBinStat struct {
	ts int
	tv int
}

type tsTvCountStat struct {
	altCount int
	ts       int
	tv       int
}

// tsTvQualStat records the QUAL score of one biallelic SNP and whether it is a
// transition; used to build the --TsTv-by-qual report.
type tsTvQualStat struct {
	qual float64
	isTs bool
}

type indvDepthStat struct {
	name   string
	sum    int
	nSites int
}

// indvHetStat accumulates the PLINK-style per-individual inbreeding
// coefficient (--het) statistics over biallelic-SNP, fully-diploid sites
// where the non-reference allele is polymorphic. obsHom is the observed
// number of homozygous diploid genotypes; expHom is the expected number by
// chance (1 - 2pq * n/(n-1) per site); nSites is the count of included
// diploid genotypes. See variant_file_output.cpp:165-279.
type indvHetStat struct {
	name   string
	obsHom int
	expHom float64
	nSites int
}

type windowPiStat struct {
	chrom    string
	winStart int
	winEnd   int
	nSites   int
	pi       float64
}

type tajimaDStat struct {
	chrom    string
	binStart int
	binEnd   int
	nSites   int
	tajimaD  float64
}

type fstStat struct {
	chrom string
	pos   int
	fst   float64
}

type singletonStat struct {
	chrom string
	pos   int
	// kind is "S" (singleton: ALT count == 1) or "D" (private doubleton:
	// a single individual carries both ALT alleles). Matches upstream
	// vcftools .singletons SINGLETON/DOUBLETON column.
	kind   string
	allele string
	// indv is the sample id carrying the rare allele. Empty when we
	// cannot determine it.
	indv string
}

// newStatistics creates a new statistics collector
func newStatistics(header *vcf.Header) *statistics {
	return &statistics{
		header:         header,
		indvMissing:    make(map[string]*indvMissingStat),
		indvHet:        make(map[string]*indvHetStat),
		indvDepth:      make(map[string]*indvDepthStat),
		tsTvBins:       make(map[string]map[int]*tsTvBinStat),
		tsTvByCount:    make(map[int]*tsTvCountStat),
		snpDensityBins: make(map[string][]int),
		snpDenPrevPos:  -1,
		filterCounts:   make(map[string]int),
		filterTs:       make(map[string]int),
		filterTv:       make(map[string]int),
		indelLenHist:   make(map[int]int),
	}
}

// addVariant adds a variant to the statistics. genoFiltered names the samples
// whose genotype was dropped by a genotype-level filter at this site (nil when
// none were); it is consumed only by the per-individual missingness path.
func (s *statistics) addVariant(v *vcf.Variant, params *Params, genoFiltered map[string]bool) {
	// Allele frequency. The "2" (suppress-allele) variants drive the same
	// per-site collection as the plain --freq / --counts flags.
	if params.Freq || params.Counts || params.Freq2 || params.Counts2 {
		s.addFrequencyStat(v, params.Derived)
	}

	// Site depth
	if params.SiteDepth || params.SiteMeanDepth {
		s.addSiteDepthStat(v)
	}

	// Site quality
	if params.SiteQuality {
		s.addSiteQualityStat(v)
	}

	// Site missingness
	if params.MissingSite {
		s.addSiteMissingStat(v)
	}

	// Individual missingness
	if params.MissingIndv {
		s.addIndvMissingStat(v, genoFiltered)
	}

	// Hardy-Weinberg
	if params.Hardy {
		s.addHWEStat(v)
	}

	// Ts/Tv
	if params.TsTvSummary || params.TsTvBinSize > 0 {
		s.addTsTvStat(v, params.TsTvBinSize)
	}

	// Ts/Tv by alternate allele count
	if params.TsTvByCount {
		s.addTsTvByCountStat(v)
	}

	// Ts/Tv by quality score
	if params.TsTvByQual {
		s.addTsTvByQualStat(v)
	}

	// Per-individual mean depth
	if params.Depth {
		s.addIndvDepthStat(v)
	}

	// Per-genotype depth matrix
	if params.GenoDepth {
		s.addGenoDepthStat(v)
	}

	// Indel length histogram
	if params.HistIndelLen {
		s.addIndelLenStat(v)
	}

	// Site pi (nucleotide diversity)
	if params.SitePi {
		s.addSitePiStat(v)
	}
	// Windowed pi uses a distinct accumulation (mismatches/comparisons per
	// site, binned over windows) — see output_windowed_nucleotide_diversity.
	if params.WindowPi > 0 {
		s.addWindowPiSite(v)
	}

	// Tajima's D (collect per-SNP data; computed at output time)
	if params.TajimaD > 0 {
		s.addTajimaDStat(v)
	}

	// Phase 2: Population genetics statistics

	// Heterozygosity
	if params.Het {
		s.addHetStat(v)
	}

	// Singletons
	if params.Singletons {
		s.addSingletonStat(v)
	}

	// FILTER summary
	if params.FilterSummary {
		s.addFilterCount(v)
	}

	// SNP density
	if params.SNPDensity > 0 {
		s.addSNPDensityStat(v, params.SNPDensity)
	}

	// Weir & Cockerham 1984 Fst (per-site components, accumulated for the
	// per-site and optional windowed output).
	if s.weirFst != nil && len(params.WeirFstPop) >= 2 {
		s.weirFst.addVariant(v)
	}
}

// addFrequencyStat adds allele frequency statistics. When `derived` is
// true, the site's INFO/AA tag is consulted and the entry is dropped if
// (a) AA is missing/`.`/`?`, or (b) AA does not match REF or ALT (case-
// insensitive). When AA matches ALT the row is recorded with
// derivedSwap=true so outputFrequency emits ALT before REF, matching the
// reorder upstream applies via `aa_idx` (variant_file_output.cpp:67-101).
func (s *statistics) addFrequencyStat(v *vcf.Variant, derived bool) {
	if len(v.Samples) == 0 {
		return
	}

	// Upstream emits one row per passed site over ALL alleles (REF +
	// every ALT), including multi-allelic and monomorphic (ALT=".") sites.
	// See output_frequency (variant_file_output.cpp:10-163).
	alleles := siteAlleles(v)
	counts, nChr := siteAlleleCountsIndexed(v)

	stat := siteFreqStat{
		chrom:    v.Chrom,
		pos:      v.Pos,
		nAlleles: len(alleles),
		nChr:     nChr,
		alleles:  alleles,
		counts:   counts,
	}

	if derived {
		// Upstream uppercases INFO/AA before comparing
		// (variant_file_output.cpp:78). Sites with AA missing, ".",
		// "?", or AA that does not match any allele are skipped (the
		// `continue` branches at lines 81 and 97). The matched allele
		// is emitted first.
		aa, ok := v.Info["AA"]
		if !ok || aa == "" || aa == "." || aa == "?" {
			return
		}
		aaUp := strings.ToUpper(aa)
		found := -1
		for i, a := range alleles {
			if strings.ToUpper(a) == aaUp {
				found = i
				break
			}
		}
		if found < 0 {
			// AA does not match any allele: upstream emits a
			// one-off warning and drops the site.
			return
		}
		stat.aaIdx = found
	}

	s.siteFrequencies = append(s.siteFrequencies, stat)
}

// addSiteDepthStat adds site depth statistics
func (s *statistics) addSiteDepthStat(v *vcf.Variant) {
	if len(v.Samples) == 0 {
		return
	}

	// Upstream emits one row per passed site (no skip when no sample has
	// depth). It sums each individual's FORMAT/DP, ignoring missing/absent
	// depths (get_indv_DEPTH returns -1 for "." or a missing field).
	sum := 0
	sumsq := 0
	n := 0
	for _, sample := range v.Samples {
		dpStr, ok := sample.Data["DP"]
		if !ok || dpStr == "" || dpStr == "." {
			continue
		}
		dp, err := strconv.Atoi(dpStr)
		if err != nil || dp < 0 {
			continue
		}
		sum += dp
		sumsq += dp * dp
		n++
	}

	s.siteDepths = append(s.siteDepths, siteDepthStat{
		chrom: v.Chrom,
		pos:   v.Pos,
		sum:   sum,
		sumsq: sumsq,
		n:     n,
	})
}

// addSiteQualityStat adds site quality statistics
func (s *statistics) addSiteQualityStat(v *vcf.Variant) {
	stat := siteQualityStat{
		chrom: v.Chrom,
		pos:   v.Pos,
		qual:  v.Qual,
	}
	s.siteQualities = append(s.siteQualities, stat)
}

// addSiteMissingStat accumulates per-site missingness over chromosomes
// (alleles), mirroring upstream output_site_missingness
// (variant_file_output.cpp:847-918): each diploid genotype contributes 2 to
// N_DATA, each missing allele 1 to N_MISS, with phased-missing and haploid
// genotypes counting as a single chromosome.
func (s *statistics) addSiteMissingStat(v *vcf.Variant) {
	nMiss := 0
	nTot := 0
	for _, sample := range v.Samples {
		first, second, phased := parseGTAlleles(sample.Data["GT"])
		if first == -1 {
			nMiss++
		}
		if second == -1 {
			nMiss++
		}
		nTot += 2
		if second == -1 && phased {
			// Phased missing second slot ⇒ haploid genome: count one
			// chromosome (matches upstream's site_N_tot--/missing-- pair).
			nTot--
			nMiss--
		}
	}

	s.siteMissing = append(s.siteMissing, siteMissingStat{
		chrom: v.Chrom,
		pos:   v.Pos,
		nData: nTot,
		nMiss: nMiss,
	})
}

// parseGTAlleles parses a GT string into upstream's (first, second)
// allele-id pair plus whether the genotype is phased, mirroring
// vcf_entry::set_indv_GENOTYPE_and_PHASE (vcf_entry_setters.cpp). A missing
// allele "." maps to -1. A haploid genotype (a single allele slot, e.g. "0"
// or "." on male chrX) is stored exactly as upstream stores it: second = -1
// and phase = '|' (phased=true). An absent/empty GT is treated as a missing
// diploid so it contributes two missing chromosomes. Non-numeric tokens map
// to -1.
func parseGTAlleles(gt string) (first, second int, phased bool) {
	if gt == "" {
		return -1, -1, false
	}
	phased = strings.Contains(gt, "|")
	toks := strings.FieldsFunc(gt, func(r rune) bool { return r == '/' || r == '|' })
	parse := func(tok string) int {
		if tok == "." || tok == "" {
			return -1
		}
		n, err := strconv.Atoi(tok)
		if err != nil {
			return -1
		}
		return n
	}
	switch len(toks) {
	case 0:
		return -1, -1, false
	case 1:
		// Haploid: upstream sets the second allele to "." (-1) and the
		// phase to '|'.
		return parse(toks[0]), -1, true
	default:
		return parse(toks[0]), parse(toks[1]), phased
	}
}

// addIndvMissingStat adds individual missingness statistics. The
// genoFiltered set names the samples whose genotype was dropped by a
// genotype-level filter (--minDP/--maxDP/--minGQ/--remove-filtered-geno*)
// for this site; those are counted as N_GENOTYPES_FILTERED and excluded
// from N_DATA/N_MISS, mirroring upstream's include_genotype[ui]==false
// branch in variant_file_output.cpp:821-832.
func (s *statistics) addIndvMissingStat(v *vcf.Variant, genoFiltered map[string]bool) {
	for _, sample := range v.Samples {
		stat := s.indvMissing[sample.Name]
		if stat == nil {
			stat = &indvMissingStat{name: sample.Name}
			s.indvMissing[sample.Name] = stat
		}

		if genoFiltered[sample.Name] {
			stat.nGenoFiltered++
			continue
		}

		stat.nTotal++
		if genotypeIsMissing(sample.Data["GT"]) {
			stat.nMissing++
		}
	}
}

// genotypeIsMissing reports whether a GT string has a missing call, matching
// upstream's `alleles.first == -1` test in output_indv_missingness
// (variant_file_output.cpp:829). Only the FIRST allele decides: vcftools
// parses the substring up to the first '/' or '|' (or the whole string for a
// haploid call) and maps "." to -1 (vcf_entry_setters.cpp:69-92, 121-150).
// So "./1" is missing, "0/." is NOT, "1" (haploid) is NOT, and an absent or
// empty GT is missing.
func genotypeIsMissing(gt string) bool {
	if gt == "" {
		return true
	}
	first := gt
	if i := strings.IndexAny(gt, "/|"); i >= 0 {
		first = gt[:i]
	}
	return first == "" || first == "."
}

// addHWEStat adds Hardy-Weinberg equilibrium statistics for a biallelic,
// fully-diploid SNP. It mirrors upstream output_hwe
// (variant_file_output.cpp:278-365): the REF-allele frequency is taken over
// non-missing chromosomes, the expected genotype counts and chi-square follow
// the PLINK formulae, and the three p-values come from the exact SNPHWE test.
func (s *statistics) addHWEStat(v *vcf.Variant) {
	// Biallelic only (exactly one non-"." ALT).
	if len(siteAlleles(v)) != 2 {
		return
	}
	// Fully diploid only.
	if !siteIsDiploid(v) {
		return
	}

	counts, nChr := siteAlleleCountsIndexed(v)
	if nChr == 0 {
		return
	}
	hom1, het, hom2, _ := countDiploidGenotypes(v)
	n := hom1 + het + hom2
	if n == 0 {
		return
	}

	// REF-allele frequency over non-missing chromosomes (allele_counts[0]).
	freq := float64(counts[0]) / float64(nChr)
	tot := float64(n)
	expHom1 := freq * freq * tot
	expHet := 2.0 * freq * (1.0 - freq) * tot
	expHom2 := (1.0 - freq) * (1.0 - freq) * tot

	chiSq := math.Pow(float64(hom1)-expHom1, 2)/expHom1 +
		math.Pow(float64(het)-expHet, 2)/expHet +
		math.Pow(float64(hom2)-expHom2, 2)/expHom2

	pHWE, pLo, pHi := snpHWEFull(het, hom1, hom2)

	s.siteHWE = append(s.siteHWE, siteHWEStat{
		chrom:   v.Chrom,
		pos:     v.Pos,
		obsHom1: hom1,
		obsHet:  het,
		obsHom2: hom2,
		expHom1: expHom1,
		expHet:  expHet,
		expHom2: expHom2,
		chiSq:   chiSq,
		pHWE:    pHWE,
		pLo:     pLo,
		pHi:     pHi,
	})
}

// addTsTvStat adds transition/transversion statistics
func (s *statistics) addTsTvStat(v *vcf.Variant, binSize int) {
	// Only for biallelic SNPs
	if len(v.Alt) != 1 || isIndelVariant(v) {
		return
	}

	ref := strings.ToUpper(v.Ref)
	alt := strings.ToUpper(v.Alt[0])

	if len(ref) != 1 || len(alt) != 1 {
		return
	}

	// Upstream only counts substitutions whose REF+ALT is one of the 12
	// canonical nucleotide pairs (output_TsTv / output_TsTv_summary). A
	// monomorphic site (ALT ".") or a symbolic/N allele yields an unknown
	// model and is skipped, so we restrict to true nucleotide pairs here.
	if !isNucleotide(ref) || !isNucleotide(alt) {
		return
	}

	isTransition := isTransitionSNP(ref, alt)
	if isTransition {
		s.transitions++
	} else {
		s.transversions++
	}

	// Tally the substitution-model class (alphabetised pair) for the
	// --TsTv-summary output.
	if idx, ok := tstvModelIndex(ref, alt); ok {
		s.tsTvModelCounts[idx]++
	}

	// Add to the per-chromosome bin if --TsTv binning is active. Upstream
	// keys on (CHROM, bin index) where bin index = floor(POS / bin_size),
	// using the 1-based POS directly (variant_file_output.cpp:3007).
	if binSize > 0 {
		binIdx := v.Pos / binSize
		bins, ok := s.tsTvBins[v.Chrom]
		if !ok {
			bins = make(map[int]*tsTvBinStat)
			s.tsTvBins[v.Chrom] = bins
			s.tsTvChromOrder = append(s.tsTvChromOrder, v.Chrom)
		}
		bin := bins[binIdx]
		if bin == nil {
			bin = &tsTvBinStat{}
			bins[binIdx] = bin
		}
		if isTransition {
			bin.ts++
		} else {
			bin.tv++
		}
	}
}

// addSitePiStat records per-site nucleotide diversity.
func (s *statistics) addSitePiStat(v *vcf.Variant) {
	pi, ok := nucleotideDiversity(v)
	if !ok {
		return
	}
	s.sitePiValues = append(s.sitePiValues, sitePiStat{chrom: v.Chrom, pos: v.Pos, pi: pi})
}

// addWindowPiSite records the per-site mismatch/comparison quantities used by
// --window-pi. Only fully-diploid, non-fixed sites are recorded (upstream
// skips non-diploid sites and sites with zero pairwise mismatches).
func (s *statistics) addWindowPiSite(v *vcf.Variant) {
	if !siteIsDiploid(v) {
		return
	}
	counts, nChr := siteAlleleCountsIndexed(v)
	if nChr == 0 {
		return
	}
	var mismatches uint64
	for _, ac := range counts {
		mismatches += uint64(ac) * uint64(nChr-ac)
	}
	if mismatches == 0 {
		return // Site is fixed.
	}
	s.windowPiSites = append(s.windowPiSites, windowPiSite{
		chrom:         v.Chrom,
		pos:           v.Pos,
		nMismatches:   mismatches,
		nChr:          nChr,
		isPolymorphic: counts[0] < nChr,
	})
}

// siteAlleleCounts returns the count of each allele index ("0", "1", ...) across
// all non-missing chromosomes at the site, together with the total number of
// non-missing chromosomes.
func siteAlleleCounts(v *vcf.Variant) (map[string]int, int) {
	counts := make(map[string]int)
	total := 0
	for _, sample := range v.Samples {
		gt, ok := sample.Data["GT"]
		if !ok {
			continue
		}
		for _, a := range strings.FieldsFunc(gt, func(r rune) bool { return r == '/' || r == '|' }) {
			if a == "" || a == "." {
				continue
			}
			counts[a]++
			total++
		}
	}
	return counts, total
}

// siteAlleles returns the per-site allele list in upstream index order:
// REF is index 0, then each ALT in declaration order. Upstream's
// entry::add_ALT_allele drops the literal "." (monomorphic) ALT entirely, so
// a record with ALT="." yields just [REF] (N_alleles == 1). This mirrors
// entry::get_N_alleles == ALT.size()+1.
func siteAlleles(v *vcf.Variant) []string {
	alleles := make([]string, 0, len(v.Alt)+1)
	alleles = append(alleles, v.Ref)
	for _, a := range v.Alt {
		if a == "." {
			continue
		}
		alleles = append(alleles, a)
	}
	return alleles
}

// siteAlleleCountsIndexed counts alleles per allele-index (0=REF, 1=ALT[0],
// ...) the way upstream's entry::get_allele_counts does: a genotype slot is
// parsed as an integer allele index, slots outside the valid range or "."
// are skipped, and only valid slots contribute to N_non_missing_chr. The
// returned slice is parallel to siteAlleles(v). This faithfully reproduces
// upstream over multi-allelic and monomorphic sites alike.
func siteAlleleCountsIndexed(v *vcf.Variant) (counts []int, nChr int) {
	alleles := siteAlleles(v)
	counts = make([]int, len(alleles))
	for _, sample := range v.Samples {
		gt, ok := sample.Data["GT"]
		if !ok {
			continue
		}
		for _, a := range strings.FieldsFunc(gt, func(r rune) bool { return r == '/' || r == '|' }) {
			if a == "" || a == "." {
				continue
			}
			idx, err := strconv.Atoi(a)
			if err != nil || idx < 0 || idx >= len(counts) {
				continue
			}
			counts[idx]++
			nChr++
		}
	}
	return counts, nChr
}

// siteIsDiploid reports whether every non-empty-GT sample at the site is
// diploid (exactly two allele slots, as in "0/0", "1|2", or the fully-missing
// "./."). It mirrors upstream entry::is_diploid (entry_getters.cpp:94), which
// skips a record entirely for --site-pi when any included sample has ploidy
// != 2. A fully-missing call "./." still counts as ploidy 2; an empty GT field
// is ignored (upstream parses no genotype for it).
func siteIsDiploid(v *vcf.Variant) bool {
	split := func(r rune) bool { return r == '/' || r == '|' }
	for _, sample := range v.Samples {
		gt, ok := sample.Data["GT"]
		if !ok || gt == "" {
			continue
		}
		if len(strings.FieldsFunc(gt, split)) != 2 {
			return false
		}
	}
	return true
}

// nucleotideDiversity computes per-site nucleotide diversity as defined by
// vcftools --site-pi (output_per_site_nucleotide_diversity,
// variant_file_output.cpp:3870):
//
//	pi = mismatches / (n * (n - 1)),  mismatches = sum_a c_a * (n - c_a)
//	   = (n^2 - sum_a c_a^2) / (n * (n - 1))
//
// where c_a is the count of allele a across the n non-missing chromosomes at
// the site. The two forms are algebraically identical: upstream's per-allele
// pairwise-mismatch loop equals the textbook closed form, so our port matches
// upstream values exactly.
//
// Upstream only emits diploid sites (entry::is_diploid), so this returns
// ok=false for any site with a non-diploid included sample, and ok=false when
// fewer than two chromosomes have data.
func nucleotideDiversity(v *vcf.Variant) (pi float64, ok bool) {
	if !siteIsDiploid(v) {
		return 0, false
	}
	counts, n := siteAlleleCounts(v)
	if n < 2 {
		return 0, false
	}
	sumSq := 0
	for _, c := range counts {
		sumSq += c * c
	}
	return float64(n*n-sumSq) / float64(n*(n-1)), true
}

// isNucleotide reports whether s is a single canonical DNA base (A, C, G, or
// T, uppercase). Upstream's Ts/Tv model map only contains the 12 ordered pairs
// of these four bases, so non-nucleotide "alleles" (e.g. ".", "N", or symbolic
// ALTs) are excluded from the Ts/Tv tallies.
func isNucleotide(s string) bool {
	switch s {
	case "A", "C", "G", "T":
		return true
	default:
		return false
	}
}

// isTransitionSNP reports whether a single-nucleotide ref/alt pair is a
// transition (A<->G or C<->T); anything else is a transversion.
func isTransitionSNP(ref, alt string) bool {
	return (ref == "A" && alt == "G") || (ref == "G" && alt == "A") ||
		(ref == "C" && alt == "T") || (ref == "T" && alt == "C")
}

// tstvModelIndex maps an unordered (ref, alt) single-nucleotide
// substitution to its canonical index in the upstream `--TsTv-summary`
// table:
//
//	0: AC  1: AG  2: AT  3: CG  4: CT  5: GT
//
// The alphabetisation matches the upstream ordering exactly.
func tstvModelIndex(ref, alt string) (int, bool) {
	if len(ref) != 1 || len(alt) != 1 {
		return 0, false
	}
	r, a := ref[0], alt[0]
	if r > a {
		r, a = a, r
	}
	switch string([]byte{r, a}) {
	case "AC":
		return 0, true
	case "AG":
		return 1, true
	case "AT":
		return 2, true
	case "CG":
		return 3, true
	case "CT":
		return 4, true
	case "GT":
		return 5, true
	}
	return 0, false
}

// addTsTvByCountStat bins a biallelic SNP by its alternate-allele count and
// records whether it is a transition or transversion.
func (s *statistics) addTsTvByCountStat(v *vcf.Variant) {
	if len(v.Alt) != 1 || isIndelVariant(v) {
		return
	}
	ref := strings.ToUpper(v.Ref)
	alt := strings.ToUpper(v.Alt[0])
	if len(ref) != 1 || len(alt) != 1 {
		return
	}
	// Upstream only bins substitutions whose REF+ALT model is one of the 12
	// canonical nucleotide pairs (output_TsTv_by_count, the model-map lookup
	// at variant_file_output.cpp:3190). A "SNP" whose ALT is "." (a
	// monomorphic site) passes is_biallelic_SNP's size checks but yields an
	// unknown model and is logged-and-skipped upstream, so we skip it too.
	if !isNucleotide(ref) || !isNucleotide(alt) {
		return
	}
	counts, _ := siteAlleleCounts(v)
	ac := counts["1"]
	stat := s.tsTvByCount[ac]
	if stat == nil {
		stat = &tsTvCountStat{altCount: ac}
		s.tsTvByCount[ac] = stat
	}
	if isTransitionSNP(ref, alt) {
		stat.ts++
	} else {
		stat.tv++
	}
}

// addTsTvByQualStat records the QUAL score of a biallelic SNP together with
// whether it is a transition, for the --TsTv-by-qual report. SNPs with a
// missing QUAL (v.Qual < 0) are skipped.
func (s *statistics) addTsTvByQualStat(v *vcf.Variant) {
	if len(v.Alt) != 1 || isIndelVariant(v) {
		return
	}
	if v.Qual < 0 {
		return
	}
	ref := strings.ToUpper(v.Ref)
	alt := strings.ToUpper(v.Alt[0])
	if len(ref) != 1 || len(alt) != 1 {
		return
	}
	s.tsTvByQual = append(s.tsTvByQual, tsTvQualStat{qual: v.Qual, isTs: isTransitionSNP(ref, alt)})
}

// addIndvDepthStat accumulates per-individual depth from the FORMAT/DP field.
func (s *statistics) addIndvDepthStat(v *vcf.Variant) {
	for _, sample := range v.Samples {
		stat := s.indvDepth[sample.Name]
		if stat == nil {
			stat = &indvDepthStat{name: sample.Name}
			s.indvDepth[sample.Name] = stat
		}
		dpStr, ok := sample.Data["DP"]
		if !ok {
			continue
		}
		dp, err := strconv.Atoi(strings.TrimSpace(dpStr))
		if err != nil {
			continue
		}
		stat.sum += dp
		stat.nSites++
	}
}

// parseDP returns the FORMAT/DP value for a sample, or -1 if absent/unparseable.
func parseDP(sample vcf.Sample) int {
	dpStr, ok := sample.Data["DP"]
	if !ok {
		return -1
	}
	dp, err := strconv.Atoi(strings.TrimSpace(dpStr))
	if err != nil {
		return -1
	}
	return dp
}

// addGenoDepthStat records the per-individual depth at one site.
func (s *statistics) addGenoDepthStat(v *vcf.Variant) {
	depths := make([]int, len(v.Samples))
	for i, sample := range v.Samples {
		depths[i] = parseDP(sample)
	}
	s.genoDepths = append(s.genoDepths, genoDepthSite{chrom: v.Chrom, pos: v.Pos, depths: depths})
}

// addIndelLenStat records the length of each indel allele (positive for
// insertions, negative for deletions) relative to the reference.
func (s *statistics) addIndelLenStat(v *vcf.Variant) {
	refLen := len(v.Ref)
	for _, alt := range v.Alt {
		// Skip symbolic / structural alleles such as <DEL>.
		if strings.ContainsAny(alt, "<>[]*") {
			continue
		}
		d := len(alt) - refLen
		if d == 0 {
			continue
		}
		s.indelLenHist[d]++
		s.indelLenTotal++
	}
}

// addTajimaDStat collects the per-site data needed for Tajima's D over a
// biallelic SNP.
func (s *statistics) addTajimaDStat(v *vcf.Variant) {
	if len(v.Alt) != 1 || isIndelVariant(v) {
		return
	}
	pi, ok := nucleotideDiversity(v)
	if !ok {
		return
	}
	_, n := siteAlleleCounts(v)
	s.tajimaDSites = append(s.tajimaDSites, tajimaDSite{chrom: v.Chrom, pos: v.Pos, pi: pi, nChr: n})
}

// Phase 2 statistics methods

// addHetStat accumulates the PLINK-style per-individual inbreeding
// coefficient inputs (--het). It mirrors upstream output_het
// (variant_file_output.cpp:165-279): only biallelic SNPs that are fully
// diploid and polymorphic (the non-ref frequency strictly between 0 and 1)
// contribute. For each included diploid genotype it counts observed
// homozygotes and accumulates the per-site expected-homozygote term.
func (s *statistics) addHetStat(v *vcf.Variant) {
	// Ensure a stat slot exists for every sample (so an individual with no
	// included sites still has an entry — upstream tracks all N_indv).
	for _, sample := range v.Samples {
		if s.indvHet[sample.Name] == nil {
			s.indvHet[sample.Name] = &indvHetStat{name: sample.Name}
		}
	}

	// Biallelic only (N_alleles == 2, i.e. exactly one non-"." ALT).
	alleles := siteAlleles(v)
	if len(alleles) != 2 {
		return
	}
	// Fully diploid only.
	if !siteIsDiploid(v) {
		return
	}

	counts, nChr := siteAlleleCountsIndexed(v)
	if nChr == 0 {
		return
	}
	freq := float64(counts[1]) / float64(nChr)
	// Skip monomorphic sites (freq == 0 or 1, to within machine epsilon).
	const eps = 2.220446049250313e-16
	if freq <= eps || 1.0-freq <= eps {
		return
	}

	expTerm := 1.0 - (2.0 * freq * (1.0 - freq) * (float64(nChr) / (float64(nChr) - 1.0)))

	for _, sample := range v.Samples {
		gt, ok := sample.Data["GT"]
		if !ok || gt == "" {
			continue
		}
		gtAlleles := strings.FieldsFunc(gt, func(r rune) bool { return r == '/' || r == '|' })
		if len(gtAlleles) != 2 {
			continue
		}
		a0, a1 := gtAlleles[0], gtAlleles[1]
		// Both alleles must be non-missing (upstream: alleles.first > -1
		// && alleles.second > -1).
		if a0 == "." || a0 == "" || a1 == "." || a1 == "" {
			continue
		}
		stat := s.indvHet[sample.Name]
		stat.nSites++
		if a0 == a1 {
			stat.obsHom++
		}
		stat.expHom += expTerm
	}
}

// addSingletonStat identifies singleton sites (alleles present in only one sample)
func (s *statistics) addSingletonStat(v *vcf.Variant) {
	if len(v.Samples) == 0 {
		return
	}

	sampleNames := s.header.Samples
	alleles := siteAlleles(v)
	counts, _ := siteAlleleCountsIndexed(v)

	// Pre-parse each sample's diploid allele-id pair once.
	type genoPair struct{ first, second int }
	genos := make([]genoPair, len(v.Samples))
	for i, sample := range v.Samples {
		first, second, _ := parseGTAlleles(sample.Data["GT"])
		genos[i] = genoPair{first, second}
	}

	// Iterate allele index 0..N_alleles-1 in order (upstream a-loop). For a
	// singleton (count==1) emit the first individual carrying allele a in
	// either slot; for a private doubleton (count==2) emit the first
	// individual HOMOZYGOUS for a. The ALLELE column is the allele string.
	for a := 0; a < len(alleles); a++ {
		switch counts[a] {
		case 1:
			for i := range v.Samples {
				if genos[i].first == a || genos[i].second == a {
					name := ""
					if i < len(sampleNames) {
						name = sampleNames[i]
					}
					s.singletonSites = append(s.singletonSites, singletonStat{
						chrom: v.Chrom, pos: v.Pos, kind: "S", allele: alleles[a], indv: name,
					})
					break
				}
			}
		case 2:
			for i := range v.Samples {
				if genos[i].first == a && genos[i].second == a {
					name := ""
					if i < len(sampleNames) {
						name = sampleNames[i]
					}
					s.singletonSites = append(s.singletonSites, singletonStat{
						chrom: v.Chrom, pos: v.Pos, kind: "D", allele: alleles[a], indv: name,
					})
					break
				}
			}
		}
	}
}

// addFilterCount accumulates the --FILTER-summary tallies. The key is the
// WHOLE FILTER-field string (e.g. "PASS", ".", "q10;s50"), not individual
// tags, mirroring upstream output_FILTER_summary which keys on get_FILTER().
// Ts/Tv is derived from the REF+ALT[0] substitution model: A<->G / C<->T are
// transitions, the other four ordered pairs are transversions.
func (s *statistics) addFilterCount(v *vcf.Variant) {
	filter := filterFieldString(v)
	s.filterCounts[filter]++

	if len(v.Alt) == 0 {
		return
	}
	if idx, ok := tstvModelIndex(v.Ref, v.Alt[0]); ok {
		switch idx {
		case 1, 4: // AG, CT
			s.filterTs[filter]++
		case 0, 2, 3, 5: // AC, AT, CG, GT
			s.filterTv[filter]++
		}
	}
}

// filterFieldString reconstructs the original FILTER column string from the
// parsed slice (parseFilter joins on ';' and keeps "." / "PASS" as singletons).
func filterFieldString(v *vcf.Variant) string {
	if len(v.Filter) == 0 {
		return "."
	}
	return strings.Join(v.Filter, ";")
}

// addSNPDensityStat bins a variant for --SNPdensity, mirroring upstream
// output_SNP_density (variant_file_output.cpp): a record is counted only when
// its ALT is not "." (i.e. it is a real variant) and it is not a duplicate of
// the immediately preceding (CHROM, POS). Bins are per chromosome, indexed by
// POS/bin_size. Indels are counted (upstream counts all non-monomorphic
// variants here, despite the "SNP" name).
func (s *statistics) addSNPDensityStat(v *vcf.Variant, binSize int) {
	chrom := v.Chrom
	// Track first-seen chromosome order (recorded for every record, even
	// monomorphic ones, matching upstream's chrs.push_back placement).
	if _, seen := s.snpDensityBins[chrom]; !seen {
		s.snpDensityChrs = append(s.snpDensityChrs, chrom)
		s.snpDensityBins[chrom] = nil
	}

	altMonomorphic := len(v.Alt) == 0 || (len(v.Alt) == 1 && v.Alt[0] == ".")
	if !altMonomorphic && (v.Pos != s.snpDenPrevPos || chrom != s.snpDenPrevChrom) {
		idx := v.Pos / binSize
		bins := s.snpDensityBins[chrom]
		if idx >= len(bins) {
			grown := make([]int, idx+1)
			copy(grown, bins)
			bins = grown
			s.snpDensityBins[chrom] = bins
		}
		bins[idx]++
	}
	s.snpDenPrevPos = v.Pos
	s.snpDenPrevChrom = chrom
}

// Output functions

// formatFreq formats an allele frequency the way upstream's
// `ostream << double` does in output_frequency (variant_file_output.cpp:131).
// Upstream writes the value straight to a default-configured ostream, which
// uses C++ `defaultfloat` with precision 6: six *significant* digits with
// trailing zeros stripped (so 3/4 prints as "0.75", not "0.750000", and 1/3
// prints as "0.333333"). Go's `strconv.FormatFloat(v, 'g', 6, 64)` reproduces
// that byte-for-byte. Allele frequencies lie in [0, 1], so we never reach the
// exponent-notation threshold for ordinary inputs.
func formatFreq(v float64) string {
	return formatCppDefault(v)
}

// formatCppDefault renders a float the way a default-configured C++ ostream
// (defaultfloat, precision 6) would: up to 6 significant digits, trailing
// zeros stripped, switching to scientific notation only outside the
// [1e-4, 1e6) magnitude band. NaN prints as "-nan" and infinities as
// "inf"/"-inf", matching glibc's printf used by libstdc++. This is the
// formatting upstream relies on for nearly every floating-point statistic
// it writes straight to an ostream.
func formatCppDefault(v float64) string {
	switch {
	case math.IsNaN(v):
		return "-nan"
	case math.IsInf(v, 1):
		return "inf"
	case math.IsInf(v, -1):
		return "-inf"
	}
	return strconv.FormatFloat(v, 'g', 6, 64)
}

// outputFrequency outputs allele frequency statistics
func (s *statistics) outputFrequency(prefix string, counts bool, suppressAlleles bool) error {
	suffix := ".frq"
	if counts {
		suffix = ".frq.count"
	}

	f, err := iohelper.OpenWriter(prefix + suffix)
	if err != nil {
		return err
	}
	defer f.Close()

	// Header. Upstream vcftools emits a single literal `{ALLELE:FREQ}` /
	// `{ALLELE:COUNT}` column header — the curly braces are part of the
	// header text, not a per-allele wrapper. The data rows below have one
	// tab-separated `allele:value` entry per allele (no braces). See
	// reference_code/vcftools/src/cpp/variant_file_output.cpp around line 36.
	//
	// Under --freq2 / --counts2 (suppress_allele_output, suppressAlleles
	// here) the allele label is stripped from both header and rows: the
	// header collapses to `{FREQ}` / `{COUNT}` and each row prints the bare
	// tab-separated value for every allele with no `allele:` prefix
	// (variant_file_output.cpp:42-48, 118-127, 146-156). The output file is
	// the SAME `.frq` / `.frq.count` as --freq / --counts — upstream's
	// `--freq2` only toggles the suppress flag, it does not change the suffix.
	switch {
	case counts && suppressAlleles:
		fmt.Fprintln(f, "CHROM\tPOS\tN_ALLELES\tN_CHR\t{COUNT}")
	case counts:
		fmt.Fprintln(f, "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:COUNT}")
	case suppressAlleles:
		fmt.Fprintln(f, "CHROM\tPOS\tN_ALLELES\tN_CHR\t{FREQ}")
	default:
		fmt.Fprintln(f, "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:FREQ}")
	}

	// Data. Each site prints every allele (REF + all ALTs) in index order,
	// but the ancestral allele (aaIdx, set only under --derived) is emitted
	// first. This mirrors upstream's `aa_idx`-keyed loop in
	// variant_file_output.cpp:107-159: print allele[aaIdx] first, then the
	// remaining alleles ui != aaIdx in index order. Frequencies are
	// count/N_CHR with C++ defaultfloat formatting (6 significant digits).
	for _, stat := range s.siteFrequencies {
		// Build the emission order: aaIdx first, then the rest in order.
		order := make([]int, 0, len(stat.alleles))
		order = append(order, stat.aaIdx)
		for ui := range stat.alleles {
			if ui != stat.aaIdx {
				order = append(order, ui)
			}
		}

		var sb strings.Builder
		fmt.Fprintf(&sb, "%s\t%d\t%d\t%d", stat.chrom, stat.pos, stat.nAlleles, stat.nChr)
		for _, ui := range order {
			switch {
			case counts && suppressAlleles:
				fmt.Fprintf(&sb, "\t%d", stat.counts[ui])
			case counts:
				fmt.Fprintf(&sb, "\t%s:%d", stat.alleles[ui], stat.counts[ui])
			case suppressAlleles:
				fmt.Fprintf(&sb, "\t%s", formatFreq(siteFreq(stat.counts[ui], stat.nChr)))
			default:
				fmt.Fprintf(&sb, "\t%s:%s", stat.alleles[ui], formatFreq(siteFreq(stat.counts[ui], stat.nChr)))
			}
		}
		sb.WriteByte('\n')
		if _, err := io.WriteString(f, sb.String()); err != nil {
			return err
		}
	}

	return nil
}

// siteFreq returns count/nChr, matching upstream's division (which yields a
// NaN/Inf for nChr==0 sites; formatFreq renders those as upstream's C++
// ostream would).
func siteFreq(count, nChr int) float64 {
	return float64(count) / float64(nChr)
}

// outputSiteMeanDepth outputs mean depth per site
func (s *statistics) outputSiteMeanDepth(prefix string) error {
	f, err := iohelper.OpenWriter(prefix + ".ldepth.mean")
	if err != nil {
		return err
	}
	defer f.Close()

	fmt.Fprintln(f, "CHROM\tPOS\tMEAN_DEPTH\tVAR_DEPTH")

	for _, stat := range s.siteDepths {
		// Upstream: mean = sum/n, var = ((sumsq/n) - mean^2)*n/(n-1),
		// both written to a default ostream. For n==0 the division
		// yields NaN ("-nan"); for n==1 the variance denominator is 0
		// (also "-nan"). formatCppDefault renders those like upstream.
		mean := float64(stat.sum) / float64(stat.n)
		variance := ((float64(stat.sumsq) / float64(stat.n)) - mean*mean) * float64(stat.n) / float64(stat.n-1)
		fmt.Fprintf(f, "%s\t%d\t%s\t%s\n", stat.chrom, stat.pos,
			formatCppDefault(mean), formatCppDefault(variance))
	}

	return nil
}

// outputSiteDepth outputs sum depth per site
func (s *statistics) outputSiteDepth(prefix string) error {
	f, err := iohelper.OpenWriter(prefix + ".ldepth")
	if err != nil {
		return err
	}
	defer f.Close()

	// SUM_DEPTH is the summed per-individual depth; SUMSQ_DEPTH is the sum
	// of squared per-individual depths (variant_file_output.cpp:3460-3461).
	fmt.Fprintln(f, "CHROM\tPOS\tSUM_DEPTH\tSUMSQ_DEPTH")

	for _, stat := range s.siteDepths {
		fmt.Fprintf(f, "%s\t%d\t%d\t%d\n", stat.chrom, stat.pos, stat.sum, stat.sumsq)
	}

	return nil
}

// outputSiteQuality outputs quality scores per site
func (s *statistics) outputSiteQuality(prefix string) error {
	f, err := iohelper.OpenWriter(prefix + ".lqual")
	if err != nil {
		return err
	}
	defer f.Close()

	fmt.Fprintln(f, "CHROM\tPOS\tQUAL")

	for _, stat := range s.siteQualities {
		// Upstream writes get_QUAL() straight to a default ostream
		// (variant_file_output.cpp:32). Missing QUAL ("." in the VCF) is
		// stored as the sentinel -1 and prints literally as "-1". Present
		// values use C++ defaultfloat precision 6.
		fmt.Fprintf(f, "%s\t%d\t%s\n", stat.chrom, stat.pos, formatCppDefault(stat.qual))
	}

	return nil
}

// outputMissingSite outputs site missingness
func (s *statistics) outputMissingSite(prefix string) error {
	f, err := iohelper.OpenWriter(prefix + ".lmiss")
	if err != nil {
		return err
	}
	defer f.Close()

	// Upstream uses "CHR" (not "CHROM") in this header — same convention
	// as the upstream .hwe / .geno.ld / .hap.ld outputs (predating the
	// later "CHROM" convention used by --freq / --depth etc).
	fmt.Fprintln(f, "CHR\tPOS\tN_DATA\tN_GENOTYPE_FILTERED\tN_MISS\tF_MISS")

	for _, stat := range s.siteMissing {
		// N_GENOTYPE_FILTERED stays 0 (genotype-level filtering not wired
		// into this accumulator). F_MISS uses C++ defaultfloat precision 6.
		fMiss := float64(stat.nMiss) / float64(stat.nData)
		fmt.Fprintf(f, "%s\t%d\t%d\t0\t%d\t%s\n",
			stat.chrom, stat.pos, stat.nData, stat.nMiss, formatCppDefault(fMiss))
	}

	return nil
}

// outputMissingIndv outputs individual missingness
func (s *statistics) outputMissingIndv(prefix string) error {
	f, err := iohelper.OpenWriter(prefix + ".imiss")
	if err != nil {
		return err
	}
	defer f.Close()

	fmt.Fprintln(f, "INDV\tN_DATA\tN_GENOTYPES_FILTERED\tN_MISS\tF_MISS")

	// Sort by individual name
	var names []string
	for name := range s.indvMissing {
		names = append(names, name)
	}
	sort.Strings(names)

	for _, name := range names {
		stat := s.indvMissing[name]
		// F_MISS = N_MISS / N_DATA. Upstream prints it via the C++ ostream
		// default (six significant digits, %g style); N_DATA == 0 yields a
		// 0/0 division whose libstdc++ output is "-nan"
		// (variant_file_output.cpp:841). formatDiscordance reproduces both.
		fmt.Fprintf(f, "%s\t%d\t%d\t%d\t%s\n",
			stat.name, stat.nTotal, stat.nGenoFiltered, stat.nMissing,
			formatDiscordance(stat.nMissing, stat.nTotal))
	}

	return nil
}

// outputHWE outputs Hardy-Weinberg equilibrium statistics
func (s *statistics) outputHWE(prefix string) error {
	f, err := iohelper.OpenWriter(prefix + ".hwe")
	if err != nil {
		return err
	}
	defer f.Close()

	// Upstream uses "CHR" (not "CHROM"). The E(HOM1/HET/HOM2) triple is
	// fixed-notation precision 2; ChiSq_HWE and the three p-values are
	// scientific-notation precision 6 (variant_file_output.cpp:355-363).
	fmt.Fprintln(f, "CHR\tPOS\tOBS(HOM1/HET/HOM2)\tE(HOM1/HET/HOM2)\tChiSq_HWE\tP_HWE\tP_HET_DEFICIT\tP_HET_EXCESS")

	for _, stat := range s.siteHWE {
		// A monomorphic site has zero expected homozygotes, so the ChiSq sum is
		// 0/0 = NaN. Upstream is C++ and printf renders that NaN as glibc's
		// "-nan" (the quiet-NaN sign bit is set); Go's %e renders it as "NaN".
		// Emit "-nan" to match byte-for-byte. The p-values are computed by the
		// exact test and never go NaN here, so only ChiSq needs the guard.
		chiSqStr := fmt.Sprintf("%.6e", stat.chiSq)
		if math.IsNaN(stat.chiSq) {
			chiSqStr = "-nan"
		}
		fmt.Fprintf(f, "%s\t%d\t%d/%d/%d\t%.2f/%.2f/%.2f\t%s\t%.6e\t%.6e\t%.6e\n",
			stat.chrom, stat.pos,
			stat.obsHom1, stat.obsHet, stat.obsHom2,
			stat.expHom1, stat.expHet, stat.expHom2,
			chiSqStr, stat.pHWE, stat.pLo, stat.pHi)
	}

	return nil
}

// outputTsTvSummary outputs Ts/Tv summary in upstream vcftools format:
// one MODEL/COUNT table with one row per substitution class plus the two
// roll-up rows Ts and Tv. See
// reference_code/vcftools/src/cpp/variant_file_output.cpp
// output_TsTv_summary().
func (s *statistics) outputTsTvSummary(prefix string) error {
	f, err := iohelper.OpenWriter(prefix + ".TsTv.summary")
	if err != nil {
		return err
	}
	defer f.Close()

	fmt.Fprintln(f, "MODEL\tCOUNT")
	models := []string{"AC", "AG", "AT", "CG", "CT", "GT"}
	for i, m := range models {
		fmt.Fprintf(f, "%s\t%d\n", m, s.tsTvModelCounts[i])
	}
	// Ts = AG + CT, Tv = everything else.
	ts := s.tsTvModelCounts[1] + s.tsTvModelCounts[4]
	tv := s.tsTvModelCounts[0] + s.tsTvModelCounts[2] + s.tsTvModelCounts[3] + s.tsTvModelCounts[5]
	fmt.Fprintf(f, "Ts\t%d\n", ts)
	fmt.Fprintf(f, "Tv\t%d\n", tv)

	return nil
}

// outputTsTvByBin outputs Ts/Tv by genomic bins
// outputTsTvByBin writes the per-chromosome, binned Ts/Tv table produced by
// `--TsTv <bin_size>` (.TsTv). It matches upstream output_TsTv
// (variant_file_output.cpp:2962) byte-for-byte: the columns are
// CHROM/BinStart/SNP_count/Ts/Tv; chromosomes appear in first-seen order; and
// within each chromosome every bin from 0 up to the highest occupied bin is
// emitted (including empty interior bins, SNP_count 0). BinStart is the bin
// index times bin_size, SNP_count is Ts+Tv, and the ratio is Ts/Tv (0 when
// Tv == 0, matching upstream's `ratio = 0.0` initialiser).
func (s *statistics) outputTsTvByBin(prefix string, binSize int) error {
	f, err := iohelper.OpenWriter(prefix + ".TsTv")
	if err != nil {
		return err
	}
	defer f.Close()

	fmt.Fprintln(f, "CHROM\tBinStart\tSNP_count\tTs/Tv")

	for _, chrom := range s.tsTvChromOrder {
		bins := s.tsTvBins[chrom]
		maxBin := 0
		for idx := range bins {
			if idx > maxBin {
				maxBin = idx
			}
		}
		for idx := 0; idx <= maxBin; idx++ {
			ts, tv := 0, 0
			if bin := bins[idx]; bin != nil {
				ts, tv = bin.ts, bin.tv
			}
			ratio := 0.0
			if tv != 0 {
				ratio = float64(ts) / float64(tv)
			}
			fmt.Fprintf(f, "%s\t%d\t%d\t%s\n", chrom, idx*binSize, ts+tv, formatFreq(ratio))
		}
	}

	return nil
}

// outputSitePi outputs nucleotide diversity per site
func (s *statistics) outputSitePi(prefix string) error {
	f, err := iohelper.OpenWriter(prefix + ".sites.pi")
	if err != nil {
		return err
	}
	defer f.Close()

	fmt.Fprintln(f, "CHROM\tPOS\tPI")

	for _, stat := range s.sitePiValues {
		// Upstream writes pi straight to a default-configured ostream
		// (variant_file_output.cpp:3934), so "0.6" not "0.600000" and "0"
		// not "0.000000". formatFreq reproduces C++ defaultfloat precision 6.
		fmt.Fprintf(f, "%s\t%d\t%s\n", stat.chrom, stat.pos, formatFreq(stat.pi))
	}

	return nil
}

// Phase 2 output functions

// outputHet outputs heterozygosity statistics
func (s *statistics) outputHet(prefix string) error {
	f, err := iohelper.OpenWriter(prefix + ".het")
	if err != nil {
		return err
	}
	defer f.Close()

	fmt.Fprintln(f, "INDV\tO(HOM)\tE(HOM)\tN_SITES\tF")

	// Upstream iterates individuals in header order (variant_file_output.cpp
	// :263). Only those with N_SITES > 0 are emitted. E(HOM) is printed with
	// precision 1 (defaultfloat) and F with fixed precision 5; both follow
	// `out.setf(ios::fixed)` so all floats use fixed notation here.
	for _, name := range s.header.Samples {
		stat := s.indvHet[name]
		if stat == nil || stat.nSites == 0 {
			continue
		}
		fCoef := (float64(stat.obsHom) - stat.expHom) / (float64(stat.nSites) - stat.expHom)
		fmt.Fprintf(f, "%s\t%d\t%.1f\t%d\t%.5f\n",
			stat.name, stat.obsHom, stat.expHom, stat.nSites, fCoef)
	}

	return nil
}

// outputSingletons outputs singleton site analysis
func (s *statistics) outputSingletons(prefix string) error {
	f, err := iohelper.OpenWriter(prefix + ".singletons")
	if err != nil {
		return err
	}
	defer f.Close()

	// Upstream emits five columns: CHROM, POS, SINGLETON/DOUBLETON ("S" or
	// "D"), ALLELE, INDV. See reference_code/vcftools/src/cpp/
	// variant_file_output.cpp.
	fmt.Fprintln(f, "CHROM\tPOS\tSINGLETON/DOUBLETON\tALLELE\tINDV")

	for _, stat := range s.singletonSites {
		fmt.Fprintf(f, "%s\t%d\t%s\t%s\t%s\n",
			stat.chrom, stat.pos, stat.kind, stat.allele, stat.indv)
	}

	return nil
}

// outputFilterSummary writes the --FILTER-summary table
// (variant_file_output.cpp:output_FILTER_summary). Columns are
// FILTER/N_VARIANTS/N_Ts/N_Tv/Ts/Tv. Rows are sorted by ascending
// (N_VARIANTS, FILTER) and emitted in reverse, so the most frequent FILTER
// comes first (ties broken by descending FILTER string). Ts/Tv is the
// double ratio Ts/Tv via a default ostream (inf when Tv == 0).
func (s *statistics) outputFilterSummary(prefix string) error {
	f, err := iohelper.OpenWriter(prefix + ".FILTER.summary")
	if err != nil {
		return err
	}
	defer f.Close()

	fmt.Fprintln(f, "FILTER\tN_VARIANTS\tN_Ts\tN_Tv\tTs/Tv")

	type filterRow struct {
		filter string
		nsites int
	}
	rows := make([]filterRow, 0, len(s.filterCounts))
	for filter, n := range s.filterCounts {
		rows = append(rows, filterRow{filter, n})
	}
	// Ascending (nsites, filter); emitted in reverse below.
	sort.Slice(rows, func(i, j int) bool {
		if rows[i].nsites != rows[j].nsites {
			return rows[i].nsites < rows[j].nsites
		}
		return rows[i].filter < rows[j].filter
	})

	for i := len(rows) - 1; i >= 0; i-- {
		filter := rows[i].filter
		ts := s.filterTs[filter]
		tv := s.filterTv[filter]
		fmt.Fprintf(f, "%s\t%d\t%d\t%d\t%s\n",
			filter, rows[i].nsites, ts, tv,
			formatCppDefault(float64(ts)/float64(tv)))
	}

	return nil
}

// outputSNPDensity outputs SNP density in bins
func (s *statistics) outputSNPDensity(prefix string, binSize int) error {
	f, err := iohelper.OpenWriter(prefix + ".snpden")
	if err != nil {
		return err
	}
	defer f.Close()

	// Header and column layout match upstream output_SNP_density. Within a
	// chromosome, emission starts at the first non-empty bin and then prints
	// every bin (including empty interior bins) to the highest occupied bin.
	// BIN_START = bin_index*bin_size; VARIANTS/KB = count * 1000/bin_size
	// printed via a default ostream.
	fmt.Fprintln(f, "CHROM\tBIN_START\tSNP_COUNT\tVARIANTS/KB")

	perKb := 1000.0 / float64(binSize)
	for _, chrom := range s.snpDensityChrs {
		bins := s.snpDensityBins[chrom]
		output := false
		for sIdx, count := range bins {
			if count > 0 {
				output = true
			}
			if output {
				fmt.Fprintf(f, "%s\t%d\t%d\t%s\n",
					chrom, sIdx*binSize, count, formatCppDefault(float64(count)*perKb))
			}
		}
	}

	return nil
}

// outputTsTvByCount outputs Ts/Tv ratios grouped by alternate-allele count.
func (s *statistics) outputTsTvByCount(prefix string) error {
	f, err := iohelper.OpenWriter(prefix + ".TsTv.count")
	if err != nil {
		return err
	}
	defer f.Close()

	fmt.Fprintln(f, "ALT_ALLELE_COUNT\tN_Ts\tN_Tv\tTs/Tv")

	// Upstream sizes the Ts/Tv arrays at 2*N_kept_individuals and emits
	// every bin in [0, 2*N) — including empty bins whose Ts/Tv ratio is the
	// indeterminate 0.0/0.0 ("-nan"). See output_TsTv_by_count
	// (variant_file_output.cpp:3145). N_kept_individuals here is the full
	// sample set: --indv/--keep sample filtering has already been applied
	// upstream of statistics accumulation.
	nBins := 2 * len(s.header.Samples)
	for ac := 0; ac < nBins; ac++ {
		ts, tv := 0, 0
		if stat := s.tsTvByCount[ac]; stat != nil {
			ts, tv = stat.ts, stat.tv
		}
		fmt.Fprintf(f, "%d\t%d\t%d\t%s\n", ac, ts, tv, cppRatio(float64(ts), float64(tv)))
	}

	return nil
}

// cppRatio formats num/den exactly as a default-configured C++ ostream prints
// the resulting double, so vcftools' Ts/Tv columns are reproduced byte-for-byte:
//   - 0.0/0.0 prints "-nan" (the IEEE-754 quiet NaN with sign bit set that
//     glibc emits for 0.0/0.0);
//   - a positive value over zero prints "inf", a negative one "-inf";
//   - finite values use C++ defaultfloat with precision 6 (e.g. "0.5", "1").
func cppRatio(num, den float64) string {
	if den == 0 {
		switch {
		case num == 0:
			return "-nan"
		case num > 0:
			return "inf"
		default:
			return "-inf"
		}
	}
	return formatFreq(num / den)
}

// outputTsTvByQual writes transition/transversion counts and ratios bucketed by
// QUAL-score thresholds (.TsTv.qual). The thresholds are the sorted distinct
// QUAL values that appeared among biallelic SNPs. For each threshold q the
// "_LT_" columns count SNPs with qual < q and the "_GT_" columns count SNPs with
// qual >= q (both cumulative). Ts/Tv ratio columns are 0.0000 when the Tv count
// is zero.
func (s *statistics) outputTsTvByQual(prefix string) error {
	f, err := iohelper.OpenWriter(prefix + ".TsTv.qual")
	if err != nil {
		return err
	}
	defer f.Close()

	fmt.Fprintln(f, "QUAL_THRESHOLD\tN_Ts_LT_QUAL_THRESHOLD\tN_Tv_LT_QUAL_THRESHOLD\tTs/Tv_LT_QUAL_THRESHOLD\tN_Ts_GT_QUAL_THRESHOLD\tN_Tv_GT_QUAL_THRESHOLD\tTs/Tv_GT_QUAL_THRESHOLD")

	// Collect distinct sorted thresholds.
	seen := make(map[float64]bool)
	var thresholds []float64
	for _, st := range s.tsTvByQual {
		if !seen[st.qual] {
			seen[st.qual] = true
			thresholds = append(thresholds, st.qual)
		}
	}
	sort.Float64s(thresholds)

	ratio := func(ts, tv int) float64 {
		if tv == 0 {
			return 0.0
		}
		return float64(ts) / float64(tv)
	}

	for _, q := range thresholds {
		var tsLT, tvLT, tsGT, tvGT int
		for _, st := range s.tsTvByQual {
			if st.qual < q {
				if st.isTs {
					tsLT++
				} else {
					tvLT++
				}
			} else {
				if st.isTs {
					tsGT++
				} else {
					tvGT++
				}
			}
		}
		fmt.Fprintf(f, "%.4f\t%d\t%d\t%.4f\t%d\t%d\t%.4f\n",
			q, tsLT, tvLT, ratio(tsLT, tvLT), tsGT, tvGT, ratio(tsGT, tvGT))
	}

	return nil
}

// outputDepth outputs per-individual mean read depth (.idepth).
func (s *statistics) outputDepth(prefix string) error {
	f, err := iohelper.OpenWriter(prefix + ".idepth")
	if err != nil {
		return err
	}
	defer f.Close()

	fmt.Fprintln(f, "INDV\tN_SITES\tMEAN_DEPTH")

	var names []string
	for name := range s.indvDepth {
		names = append(names, name)
	}
	sort.Strings(names)

	for _, name := range names {
		stat := s.indvDepth[name]
		mean := 0.0
		if stat.nSites > 0 {
			mean = float64(stat.sum) / float64(stat.nSites)
		}
		// Upstream writes mean_depth straight to a default ostream
		// (variant_file_output.cpp:690), i.e. defaultfloat precision 6.
		fmt.Fprintf(f, "%s\t%d\t%s\n", stat.name, stat.nSites, formatCppDefault(mean))
	}

	return nil
}

// windowPiBin accumulates the four per-window quantities upstream tracks.
type windowPiBin struct {
	nVariantSites uint64 // sites with a VCF entry in the window
	nVariantPairs uint64 // pairwise comparisons at polymorphic sites
	nMismatches   uint64 // actual pairwise mismatches at polymorphic sites
	nPolymorphic  uint64 // sites polymorphic w.r.t. the reference
}

// outputWindowedPi writes nucleotide diversity in windows (.windowed.pi),
// mirroring upstream output_windowed_nucleotide_diversity
// (variant_file_output.cpp:4065-4283). Each site contributes to bins
// [first, last) where first = ceil((pos-window)/step) (clamped to 0) and
// last = ceil(pos/step). For each emitted bin, PI = N_mismatches / N_pairs,
// where N_pairs = polymorphic-site pairs + monomorphic-site pairs, and the
// reported N_VARIANTS column is the polymorphic-site count.
func (s *statistics) outputWindowedPi(prefix string, windowSize, stepSize int) error {
	if windowSize <= 0 {
		return nil
	}
	if stepSize <= 0 || stepSize > windowSize {
		stepSize = windowSize
	}

	f, err := iohelper.OpenWriter(prefix + ".windowed.pi")
	if err != nil {
		return err
	}
	defer f.Close()

	fmt.Fprintln(f, "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tN_MONOMORPHIC\tPI")

	var chromOrder []string
	bins := make(map[string][]windowPiBin)

	for _, st := range s.windowPiSites {
		if _, seen := bins[st.chrom]; !seen {
			chromOrder = append(chromOrder, st.chrom)
			bins[st.chrom] = make([]windowPiBin, 1)
		}
		first := int(math.Ceil(float64(st.pos-windowSize) / float64(stepSize)))
		if first < 0 {
			first = 0
		}
		last := int(math.Ceil(float64(st.pos) / float64(stepSize)))
		if last >= len(bins[st.chrom]) {
			grown := make([]windowPiBin, last+1)
			copy(grown, bins[st.chrom])
			bins[st.chrom] = grown
		}
		comparisons := uint64(st.nChr) * uint64(st.nChr-1)
		for idx := first; idx < last; idx++ {
			b := &bins[st.chrom][idx]
			b.nVariantSites++
			b.nVariantPairs += comparisons
			b.nMismatches += st.nMismatches
			if st.isPolymorphic {
				b.nPolymorphic++
			}
		}
	}

	nKeptChr := 2 * len(s.header.Samples)
	monoComparisons := uint64(nKeptChr) * uint64(nKeptChr-1)

	for _, chrom := range chromOrder {
		for sIdx, b := range bins[chrom] {
			if b.nPolymorphic == 0 && b.nMismatches == 0 {
				continue
			}
			// N_monomorphic = window_size - N_variant_sites (no mask).
			nMono := uint64(windowSize) - b.nVariantSites
			nPairs := b.nVariantPairs + nMono*monoComparisons
			pi := float64(b.nMismatches) / float64(nPairs)
			fmt.Fprintf(f, "%s\t%d\t%d\t%d\t%d\t%s\n",
				chrom, sIdx*stepSize+1, sIdx*stepSize+windowSize,
				b.nPolymorphic, nMono, formatCppDefault(pi))
		}
	}

	return nil
}

// outputGenoDepth writes the per-genotype read-depth matrix (.gdepth):
// CHROM, POS, then one column per individual (-1 where FORMAT/DP is absent).
func (s *statistics) outputGenoDepth(prefix string) error {
	f, err := iohelper.OpenWriter(prefix + ".gdepth")
	if err != nil {
		return err
	}
	defer f.Close()

	header := "CHROM\tPOS"
	for _, name := range s.header.Samples {
		header += "\t" + name
	}
	fmt.Fprintln(f, header)

	for _, site := range s.genoDepths {
		fmt.Fprintf(f, "%s\t%d", site.chrom, site.pos)
		for _, d := range site.depths {
			fmt.Fprintf(f, "\t%d", d)
		}
		fmt.Fprintln(f)
	}
	return nil
}

// outputIndelHist writes a histogram of indel lengths (.indel.hist):
// LENGTH (negative = deletion, positive = insertion), N_INDELS, PRCT.
func (s *statistics) outputIndelHist(prefix string) error {
	f, err := iohelper.OpenWriter(prefix + ".indel.hist")
	if err != nil {
		return err
	}
	defer f.Close()

	fmt.Fprintln(f, "LENGTH\tN_INDELS\tPRCT")

	var lengths []int
	for l := range s.indelLenHist {
		lengths = append(lengths, l)
	}
	sort.Ints(lengths)

	for _, l := range lengths {
		n := s.indelLenHist[l]
		pct := 0.0
		if s.indelLenTotal > 0 {
			pct = 100 * float64(n) / float64(s.indelLenTotal)
		}
		fmt.Fprintf(f, "%d\t%d\t%.4f\n", l, n, pct)
	}
	return nil
}

// outputTajimaD writes Tajima's D per non-overlapping window of binSize bases
// (.Tajima.D): CHROM, BIN_START, N_SNPS, TajimaD.
//
// D = (pi - thetaW) / sqrt(e1*S + e2*S*(S-1)), with pi the sum of per-site
// nucleotide diversity in the window, thetaW = S/a1, S the number of SNPs, and
// the a1/a2/e1/e2 constants derived from the number of sampled chromosomes n.
// n is taken from the SNPs in the window (the modal value); this matches the
// common case of complete data. Windows with fewer than two SNPs, or where the
// variance estimate is non-positive, are reported with TajimaD "nan".
func (s *statistics) outputTajimaD(prefix string, binSize int) error {
	if binSize <= 0 {
		return nil
	}
	f, err := iohelper.OpenWriter(prefix + ".Tajima.D")
	if err != nil {
		return err
	}
	defer f.Close()

	fmt.Fprintln(f, "CHROM\tBIN_START\tN_SNPS\tTajimaD")

	type binKey struct {
		chrom string
		start int
	}
	type binAcc struct {
		piSum  float64
		nSNPs  int
		nChrMC map[int]int // modal chromosome count
	}
	var order []binKey
	bins := make(map[binKey]*binAcc)

	for _, site := range s.tajimaDSites {
		start := ((site.pos - 1) / binSize) * binSize
		key := binKey{site.chrom, start}
		acc := bins[key]
		if acc == nil {
			acc = &binAcc{nChrMC: make(map[int]int)}
			bins[key] = acc
			order = append(order, key)
		}
		acc.piSum += site.pi
		acc.nSNPs++
		acc.nChrMC[site.nChr]++
	}

	sort.Slice(order, func(i, j int) bool {
		if order[i].chrom != order[j].chrom {
			return order[i].chrom < order[j].chrom
		}
		return order[i].start < order[j].start
	})

	for _, key := range order {
		acc := bins[key]
		d, ok := tajimasD(acc.piSum, acc.nSNPs, modalKey(acc.nChrMC))
		if ok {
			fmt.Fprintf(f, "%s\t%d\t%d\t%.5f\n", key.chrom, key.start, acc.nSNPs, d)
		} else {
			fmt.Fprintf(f, "%s\t%d\t%d\tnan\n", key.chrom, key.start, acc.nSNPs)
		}
	}
	return nil
}

// modalKey returns the most frequent key in counts (smallest on ties).
func modalKey(counts map[int]int) int {
	best, bestCount := 0, -1
	for k, c := range counts {
		if c > bestCount || (c == bestCount && k < best) {
			best, bestCount = k, c
		}
	}
	return best
}

// tajimasD computes Tajima's D from the summed per-site diversity (piSum), the
// number of segregating sites S, and the number of sampled chromosomes n.
func tajimasD(piSum float64, S, n int) (d float64, ok bool) {
	if S < 2 || n < 3 {
		return 0, false
	}
	a1, a2 := 0.0, 0.0
	for i := 1; i < n; i++ {
		a1 += 1.0 / float64(i)
		a2 += 1.0 / float64(i*i)
	}
	nf := float64(n)
	b1 := (nf + 1) / (3 * (nf - 1))
	b2 := 2 * (nf*nf + nf + 3) / (9 * nf * (nf - 1))
	c1 := b1 - 1.0/a1
	c2 := b2 - (nf+2)/(a1*nf) + a2/(a1*a1)
	e1 := c1 / a1
	e2 := c2 / (a1*a1 + a2)
	thetaW := float64(S) / a1
	variance := e1*float64(S) + e2*float64(S)*float64(S-1)
	if variance <= 0 {
		return 0, false
	}
	return (piSum - thetaW) / math.Sqrt(variance), true
}
