package vcftools

import (
	"fmt"
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
	// tsTvByBin holds per-chromosome bin counts for `--TsTv N`. Keyed
	// by chromosome, then dense slice indexed by bin (pos/binSize).
	// Mirrors upstream's `map<string, vector<int>>` Ts_counts /
	// Tv_counts at variant_file_output.cpp:2980-2981. tsTvBinChroms
	// records the first-seen order so output matches upstream's
	// declaration order.
	tsTvByBin     map[string][]tsTvBinStat
	tsTvBinChroms []string
	tsTvByCount   map[int]*tsTvCountStat
	tsTvByQual    []tsTvQualStat

	// Phase 2: Population genetics statistics
	windowPiValues []windowPiStat
	tajimaDValues  []tajimaDStat
	snpDensityBins map[int]*snpDensityStat
	fstValues      []fstStat
	filterCounts   map[string]int
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
	chrom     string
	pos       int
	nAlleles  int
	nChr      int
	refAllele string
	altAllele string
	refFreq   float64
	altFreq   float64
	refCount  int
	altCount  int
	// derivedSwap is set when --derived was requested AND the site's
	// ancestral allele (INFO/AA, case-insensitive) matched the ALT
	// allele. In that case the output puts ALT before REF so that the
	// first column is always the ancestral allele. Mirrors upstream
	// variant_file_output.cpp:67-101 (the `aa_idx` reorder loop).
	derivedSwap bool
}

type siteDepthStat struct {
	chrom      string
	pos        int
	sumDepth   int
	sumsqDepth int
	n          int
}

type siteQualityStat struct {
	chrom string
	pos   int
	qual  float64
}

type siteMissingStat struct {
	chrom string
	pos   int
	// nData is the number of alleles considered at the site (2 per
	// diploid sample, reduced for haploid/phased-missing genotypes),
	// matching upstream's site_N_tot. nMiss is the number of missing
	// alleles, nGenoFiltered the number of genotype-filtered samples.
	nData         int
	nMiss         int
	nGenoFiltered int
}

type siteHWEStat struct {
	chrom        string
	pos          int
	obsHom1      int
	obsHet       int
	obsHom2      int
	expHom1      float64
	expHet       float64
	expHom2      float64
	chiSq        float64
	pValue       float64
	pHetDeficit  float64
	pHetExcess   float64
}

type sitePiStat struct {
	chrom string
	pos   int
	pi    float64
}

type indvMissingStat struct {
	name     string
	nMissing int
	nTotal   int
	fMiss    float64
}

// tsTvBinStat holds Ts/Tv counts for one (chrom, bin) cell of the
// `--TsTv N` output. Bins are dense per-chromosome slices: a bin slot
// with zero Ts and zero Tv represents no biallelic SNPs falling into
// that bin (still emitted by upstream).
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

type indvHetStat struct {
	name    string
	nHet    int
	nHomAlt int
	nHomRef int
	nTotal  int
	hetRate float64
	// nSites is the number of biallelic, fully-diploid, non-missing sites
	// included for this individual (upstream N_sites_included).
	nSites int
	// nObsHom is the observed number of homozygous genotypes (upstream
	// N_obs_hom).
	nObsHom int
	// expHom is the accumulated expected number of homozygous genotypes,
	// summed per included site as 1 - 2*p*q*N/(N-1) (upstream
	// N_expected_hom).
	expHom float64
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

type snpDensityStat struct {
	binStart int
	binEnd   int
	nSNPs    int
	density  float64
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
		tsTvByBin:      make(map[string][]tsTvBinStat),
		tsTvByCount:    make(map[int]*tsTvCountStat),
		snpDensityBins: make(map[int]*snpDensityStat),
		filterCounts:   make(map[string]int),
		indelLenHist:   make(map[int]int),
	}
}

// addVariant adds a variant to the statistics
func (s *statistics) addVariant(v *vcf.Variant, params *Params) {
	// Allele frequency
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
		s.addIndvMissingStat(v)
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

	// Site pi (nucleotide diversity) - also required to build windowed pi
	if params.SitePi || params.WindowPi > 0 {
		s.addSitePiStat(v)
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

	// Count alleles (only for biallelic sites)
	if len(v.Alt) != 1 {
		return
	}

	refCount := 0
	altCount := 0
	totalAlleles := 0

	for _, sample := range v.Samples {
		gt, ok := sample.Data["GT"]
		if !ok || gt == "." || gt == "./." || gt == ".|." {
			continue
		}

		alleles := strings.FieldsFunc(gt, func(r rune) bool {
			return r == '/' || r == '|'
		})

		for _, allele := range alleles {
			if allele == "0" {
				refCount++
				totalAlleles++
			} else if allele == "1" {
				altCount++
				totalAlleles++
			}
		}
	}

	if totalAlleles == 0 {
		return
	}

	stat := siteFreqStat{
		chrom:     v.Chrom,
		pos:       v.Pos,
		nAlleles:  2,
		nChr:      totalAlleles,
		refAllele: v.Ref,
		altAllele: v.Alt[0],
		refCount:  refCount,
		altCount:  altCount,
		refFreq:   float64(refCount) / float64(totalAlleles),
		altFreq:   float64(altCount) / float64(totalAlleles),
	}

	if derived {
		// Upstream uppercases INFO/AA before comparing
		// (variant_file_output.cpp:78). Sites with AA missing, ".",
		// "?", or AA that does not match any allele are skipped (the
		// `continue` branches at lines 81 and 97). For our biallelic
		// fast path that means REF or ALT (upper-cased), otherwise
		// drop the row.
		aa, ok := v.Info["AA"]
		if !ok || aa == "" || aa == "." || aa == "?" {
			return
		}
		aaUp := strings.ToUpper(aa)
		refUp := strings.ToUpper(v.Ref)
		altUp := strings.ToUpper(v.Alt[0])
		switch aaUp {
		case refUp:
			// AA == REF: no reorder.
		case altUp:
			stat.derivedSwap = true
		default:
			// AA does not match any allele: upstream emits a
			// one-off warning and drops the site.
			return
		}
	}

	s.siteFrequencies = append(s.siteFrequencies, stat)
}

// addSiteDepthStat adds site depth statistics. Upstream
// (variant_file_output.cpp:3416-3452) emits one row per kept site even
// when no sample has a non-missing DP — sum/sumsq/n stay at zero, and
// the .ldepth.mean row formats as `-nan\t-nan` from the 0/0 mean and
// variance divisions. We match by always appending and only skipping
// per-sample DPs that are missing or negative.
func (s *statistics) addSiteDepthStat(v *vcf.Variant) {
	stat := siteDepthStat{chrom: v.Chrom, pos: v.Pos}

	for _, sample := range v.Samples {
		dpStr, ok := sample.Data["DP"]
		if !ok || dpStr == "." || dpStr == "" {
			continue
		}

		var dp int
		_, err := fmt.Sscanf(dpStr, "%d", &dp)
		if err != nil || dp < 0 {
			continue
		}

		stat.sumDepth += dp
		stat.sumsqDepth += dp * dp
		stat.n++
	}

	s.siteDepths = append(s.siteDepths, stat)
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

// addSiteMissingStat adds site missingness statistics
func (s *statistics) addSiteMissingStat(v *vcf.Variant) {
	if len(v.Samples) == 0 {
		return
	}

	// Upstream counts missingness per allele, not per sample:
	//   site_N_tot += 2 for each diploid sample, then both site_N_tot and
	//   site_N_missing are decremented by one for haploid genotypes (a
	//   single-allele entry, "-2") or phased-missing second alleles
	//   ("a|."). See variant_file_output.cpp:893-924.
	var nTot, nMiss, nGenoFiltered int
	for _, sample := range v.Samples {
		gt, ok := sample.Data["GT"]
		if !ok || gt == "" {
			// No GT field: treated as a fully-missing diploid genotype.
			nTot += 2
			nMiss += 2
			continue
		}

		phased := strings.ContainsRune(gt, '|')
		alleles := strings.FieldsFunc(gt, func(r rune) bool {
			return r == '/' || r == '|'
		})

		first := -1
		second := -2 // -2 marks "no second allele" (haploid)
		if len(alleles) >= 1 {
			first = parseHetAllele(alleles[0])
		}
		if len(alleles) >= 2 {
			second = parseHetAllele(alleles[1])
		}

		nTot += 2
		if first == -1 {
			nMiss++
		}
		if second == -1 {
			nMiss++
		}

		if second == -1 && phased {
			// Phased missing second allele indicates a haploid genome.
			nTot--
			nMiss--
		} else if second == -2 {
			// Haploid genotype (single allele).
			nTot--
		}
	}

	stat := siteMissingStat{
		chrom:         v.Chrom,
		pos:           v.Pos,
		nData:         nTot,
		nMiss:         nMiss,
		nGenoFiltered: nGenoFiltered,
	}

	s.siteMissing = append(s.siteMissing, stat)
}

// addIndvMissingStat adds individual missingness statistics
func (s *statistics) addIndvMissingStat(v *vcf.Variant) {
	for _, sample := range v.Samples {
		if s.indvMissing[sample.Name] == nil {
			s.indvMissing[sample.Name] = &indvMissingStat{
				name: sample.Name,
			}
		}

		stat := s.indvMissing[sample.Name]
		stat.nTotal++

		gt, ok := sample.Data["GT"]
		if !ok || gt == "." || gt == "./." || gt == ".|." {
			stat.nMissing++
		}
	}
}

// addHWEStat adds Hardy-Weinberg equilibrium statistics. Upstream
// (variant_file_output.cpp:340-351) only runs the test on biallelic SNPs
// where every kept individual reports a fully diploid call — any
// haploid genotype (e.g. "0" or "1" without a separator) disqualifies
// the site.
func (s *statistics) addHWEStat(v *vcf.Variant) {
	// Only for biallelic sites. A `.` ALT field means monomorphic
	// (REF-only) and isn't biallelic.
	if len(v.Alt) != 1 || v.Alt[0] == "." || v.Alt[0] == "" {
		return
	}

	// Count genotypes
	hom1 := 0 // 0/0
	het := 0  // 0/1 or 1/0
	hom2 := 0 // 1/1

	for _, sample := range v.Samples {
		gt, ok := sample.Data["GT"]
		if !ok || gt == "" || gt == "." {
			continue
		}

		// Reject the site outright if any included sample has a
		// non-diploid (haploid) genotype, matching upstream's
		// is_diploid() guard.
		if !strings.ContainsAny(gt, "/|") {
			return
		}

		alleles := strings.FieldsFunc(gt, func(r rune) bool {
			return r == '/' || r == '|'
		})

		if len(alleles) != 2 {
			return
		}
		if alleles[0] == "." || alleles[1] == "." {
			continue
		}

		if alleles[0] == "0" && alleles[1] == "0" {
			hom1++
		} else if alleles[0] == "1" && alleles[1] == "1" {
			hom2++
		} else if (alleles[0] == "0" && alleles[1] == "1") || (alleles[0] == "1" && alleles[1] == "0") {
			het++
		}
	}

	n := hom1 + het + hom2
	if n == 0 {
		return
	}

	// Calculate allele frequencies
	p := (float64(2*hom1) + float64(het)) / float64(2*n)
	q := 1 - p

	// Expected genotype counts under HWE
	expHom1 := p * p * float64(n)
	expHet := 2 * p * q * float64(n)
	expHom2 := q * q * float64(n)

	// Chi-square test
	chiSq := 0.0
	if expHom1 > 0 {
		chiSq += math.Pow(float64(hom1)-expHom1, 2) / expHom1
	}
	if expHet > 0 {
		chiSq += math.Pow(float64(het)-expHet, 2) / expHet
	}
	if expHom2 > 0 {
		chiSq += math.Pow(float64(hom2)-expHom2, 2) / expHom2
	}

	// Upstream uses the Wigginton/Cao/Abecasis exact test (SNPHWE)
	// for all three P-values, not the chi-square asymptotic; see
	// variant_file_output.cpp:365-372.
	pHWE, pLo, pHi := snpHWEAll(het, hom1, hom2)

	stat := siteHWEStat{
		chrom:       v.Chrom,
		pos:         v.Pos,
		obsHom1:     hom1,
		obsHet:      het,
		obsHom2:     hom2,
		expHom1:     expHom1,
		expHet:      expHet,
		expHom2:     expHom2,
		chiSq:       chiSq,
		pValue:      pHWE,
		pHetDeficit: pLo,
		pHetExcess:  pHi,
	}

	s.siteHWE = append(s.siteHWE, stat)
}

// addTsTvStat adds transition/transversion statistics. Mirrors
// upstream's biallelic-SNP gate (`if (!e->is_biallelic_SNP()) continue;`
// at variant_file_output.cpp:3001 / :3104 / :3183): only sites where
// REF and ALT are both single A/C/G/T bases participate. Monomorphic
// sites (ALT == ".") and multi-allelic sites are skipped.
func (s *statistics) addTsTvStat(v *vcf.Variant, binSize int) {
	if len(v.Alt) != 1 || isIndelVariant(v) {
		return
	}

	ref := strings.ToUpper(v.Ref)
	alt := strings.ToUpper(v.Alt[0])

	if len(ref) != 1 || len(alt) != 1 {
		return
	}

	// Only A/C/G/T substitutions count toward Ts/Tv. tstvModelIndex
	// returns ok=false for anything else (N, ".", etc.).
	idx, ok := tstvModelIndex(ref, alt)
	if !ok {
		return
	}
	s.tsTvModelCounts[idx]++
	isTransition := isTransitionSNP(ref, alt)
	if isTransition {
		s.transitions++
	} else {
		s.transversions++
	}

	// Add to per-chromosome bin if needed. Mirrors upstream
	// output_TsTv at variant_file_output.cpp:3007-3033 which keys
	// counts by (CHROM, idx=pos/binSize) and grows the per-chrom dense
	// slice on demand.
	if binSize > 0 {
		binIdx := v.Pos / binSize
		row, seen := s.tsTvByBin[v.Chrom]
		if !seen {
			s.tsTvBinChroms = append(s.tsTvBinChroms, v.Chrom)
		}
		if binIdx >= len(row) {
			grown := make([]tsTvBinStat, binIdx+1)
			copy(grown, row)
			row = grown
		}
		if isTransition {
			row[binIdx].ts++
		} else {
			row[binIdx].tv++
		}
		s.tsTvByBin[v.Chrom] = row
	}
}

// addSitePiStat records per-site nucleotide diversity.
//
// Mirrors upstream `output_per_site_nucleotide_diversity`
// (variant_file_output.cpp:3870) which gates output on
// `e->is_diploid()` and skips non-diploid sites with a one-off warning.
func (s *statistics) addSitePiStat(v *vcf.Variant) {
	if !isFullyDiploid(v) {
		return
	}
	pi, ok := nucleotideDiversity(v)
	if !ok {
		return
	}
	s.sitePiValues = append(s.sitePiValues, sitePiStat{chrom: v.Chrom, pos: v.Pos, pi: pi})
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

// nucleotideDiversity computes per-site nucleotide diversity as defined by
// vcftools --site-pi:
//
//	pi = (n^2 - sum_a c_a^2) / (n * (n - 1))
//
// where c_a is the count of allele a across the n non-missing chromosomes at
// the site. It returns ok=false when fewer than two chromosomes have data.
func nucleotideDiversity(v *vcf.Variant) (pi float64, ok bool) {
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

// addTsTvByCountStat bins a biallelic SNP by its alternate-allele count
// and records whether it is a transition or transversion. Mirrors
// upstream's biallelic-SNP gate at variant_file_output.cpp:3183 — only
// A/C/G/T substitutions participate; monomorphic ("."), multi-allelic,
// indel, and ambiguous-base sites are skipped.
func (s *statistics) addTsTvByCountStat(v *vcf.Variant) {
	if len(v.Alt) != 1 || isIndelVariant(v) {
		return
	}
	ref := strings.ToUpper(v.Ref)
	alt := strings.ToUpper(v.Alt[0])
	if len(ref) != 1 || len(alt) != 1 {
		return
	}
	if _, ok := tstvModelIndex(ref, alt); !ok {
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

// addHetStat adds heterozygosity statistics per individual
func (s *statistics) addHetStat(v *vcf.Variant) {
	// Ensure every sample has a stat entry so individuals absent from a
	// given site still appear in the output in a stable order.
	for _, sample := range v.Samples {
		if s.indvHet[sample.Name] == nil {
			s.indvHet[sample.Name] = &indvHetStat{name: sample.Name}
		}
	}

	// Upstream only uses biallelic SNPs for individual heterozygosity
	// (variant_file_output.cpp:219).
	if len(v.Alt) != 1 {
		return
	}

	// First pass: compute the site allele frequency of the non-reference
	// allele across all non-missing, diploid genotypes — and verify the
	// site is fully diploid (upstream is_diploid() check). Mixed ploidy
	// causes upstream to skip the whole site.
	type indvGT struct {
		idx     int
		a, b    int
		missing bool
	}
	gts := make([]indvGT, 0, len(v.Samples))
	altCount := 0
	nNonMissingChr := 0
	for i, sample := range v.Samples {
		gt, ok := sample.Data["GT"]
		if !ok {
			// Treated as fully-missing diploid (alleles -1/-1).
			gts = append(gts, indvGT{idx: i, a: -1, b: -1})
			continue
		}
		alleles := strings.FieldsFunc(gt, func(r rune) bool {
			return r == '/' || r == '|'
		})
		if len(alleles) != 2 {
			// Not diploid: upstream skips the entire site.
			return
		}
		ig := indvGT{idx: i}
		ig.a = parseHetAllele(alleles[0])
		ig.b = parseHetAllele(alleles[1])
		if ig.a == -1 || ig.b == -1 {
			ig.missing = true
		}
		if ig.a == 1 {
			altCount++
			nNonMissingChr++
		} else if ig.a == 0 {
			nNonMissingChr++
		}
		if ig.b == 1 {
			altCount++
			nNonMissingChr++
		} else if ig.b == 0 {
			nNonMissingChr++
		}
		gts = append(gts, ig)
	}

	if nNonMissingChr == 0 {
		return
	}
	freq := float64(altCount) / float64(nNonMissingChr)

	// Upstream skips monomorphic sites (freq at the numeric epsilon
	// boundaries), variant_file_output.cpp:240.
	const eps = 2.220446049250313e-16 // numeric_limits<double>::epsilon()
	if freq <= eps || 1.0-freq <= eps {
		return
	}

	// Per-site expected-homozygosity contribution, shared across all
	// included individuals: 1 - 2*p*q * N/(N-1).
	siteExpHom := 1.0 - (2.0 * freq * (1.0 - freq) * (float64(nNonMissingChr) / (float64(nNonMissingChr) - 1.0)))

	for _, ig := range gts {
		stat := s.indvHet[v.Samples[ig.idx].Name]
		if ig.missing {
			continue
		}
		stat.nSites++
		stat.nTotal++
		if ig.a == ig.b {
			stat.nObsHom++
			if ig.a == 0 {
				stat.nHomRef++
			} else {
				stat.nHomAlt++
			}
		} else {
			stat.nHet++
		}
		stat.expHom += siteExpHom
	}
}

// parseHetAllele parses a single genotype allele index, returning -1 for
// the missing allele ".".
func parseHetAllele(a string) int {
	if a == "." {
		return -1
	}
	n, err := strconv.Atoi(a)
	if err != nil {
		return -1
	}
	return n
}

// addSingletonStat identifies singleton sites (alleles present in only one sample)
func (s *statistics) addSingletonStat(v *vcf.Variant) {
	if len(v.Samples) == 0 {
		return
	}

	sampleNames := s.header.Samples

	// Per-allele totals + which sample(s) carry the allele.
	alleleCount := make(map[string]int)
	alleleSample := make(map[string]string)
	alleleDoubleton := make(map[string]bool)

	for i, sample := range v.Samples {
		gt, ok := sample.Data["GT"]
		if !ok || strings.Contains(gt, ".") {
			continue
		}
		alleles := strings.FieldsFunc(gt, func(r rune) bool {
			return r == '/' || r == '|'
		})
		perSample := make(map[string]int)
		for _, a := range alleles {
			perSample[a]++
		}
		for a, n := range perSample {
			alleleCount[a] += n
			if i < len(sampleNames) {
				if alleleSample[a] == "" {
					alleleSample[a] = sampleNames[i]
				} else if alleleSample[a] != sampleNames[i] {
					// More than one sample carries the allele —
					// disqualifies it from being a private doubleton.
					alleleSample[a] = ""
				}
			}
			if n == 2 {
				alleleDoubleton[a] = true
			}
		}
	}

	// Stable iteration order: sort the allele keys.
	keys := make([]string, 0, len(alleleCount))
	for a := range alleleCount {
		keys = append(keys, a)
	}
	sort.Strings(keys)

	for _, allele := range keys {
		count := alleleCount[allele]
		if allele == "0" {
			continue
		}
		var kind string
		switch {
		case count == 1:
			kind = "S"
		case count == 2 && alleleDoubleton[allele] && alleleSample[allele] != "":
			kind = "D"
		default:
			continue
		}
		s.singletonSites = append(s.singletonSites, singletonStat{
			chrom:  v.Chrom,
			pos:    v.Pos,
			kind:   kind,
			allele: allele,
			indv:   alleleSample[allele],
		})
	}
}

// addFilterCount tracks FILTER tag occurrences
func (s *statistics) addFilterCount(v *vcf.Variant) {
	if len(v.Filter) == 0 {
		s.filterCounts["PASS"]++
	} else {
		for _, filter := range v.Filter {
			s.filterCounts[filter]++
		}
	}
}

// addSNPDensityStat adds SNP density in bins
func (s *statistics) addSNPDensityStat(v *vcf.Variant, binSize int) {
	// Only count SNPs, not indels
	if isIndelVariant(v) {
		return
	}

	binIdx := v.Pos / binSize
	if s.snpDensityBins[binIdx] == nil {
		s.snpDensityBins[binIdx] = &snpDensityStat{
			binStart: binIdx * binSize,
			binEnd:   (binIdx + 1) * binSize,
		}
	}
	s.snpDensityBins[binIdx].nSNPs++
}

// Output functions

// outputFrequency outputs allele frequency statistics. When suppress is
// true (the --freq2 / --counts2 path), allele labels are omitted and only
// the numeric COUNT/FREQ columns are emitted, matching upstream's
// `params.suppress_allele_output` branches in variant_file_output.cpp:34-159.
func (s *statistics) outputFrequency(prefix string, counts, suppress bool) error {
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
	switch {
	case suppress && counts:
		fmt.Fprintln(f, "CHROM\tPOS\tN_ALLELES\tN_CHR\t{COUNT}")
	case suppress:
		fmt.Fprintln(f, "CHROM\tPOS\tN_ALLELES\tN_CHR\t{FREQ}")
	case counts:
		fmt.Fprintln(f, "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:COUNT}")
	default:
		fmt.Fprintln(f, "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:FREQ}")
	}

	// Data. When --derived was supplied and the site's ancestral allele
	// matched ALT (derivedSwap=true), the row prints ALT first so that
	// the leading column is always the ancestral (derived columns
	// follow), matching upstream's `aa_idx`-keyed loop in
	// variant_file_output.cpp:107-159.
	for _, stat := range s.siteFrequencies {
		firstAllele, secondAllele := stat.refAllele, stat.altAllele
		firstCount, secondCount := stat.refCount, stat.altCount
		firstFreq, secondFreq := stat.refFreq, stat.altFreq
		if stat.derivedSwap {
			firstAllele, secondAllele = stat.altAllele, stat.refAllele
			firstCount, secondCount = stat.altCount, stat.refCount
			firstFreq, secondFreq = stat.altFreq, stat.refFreq
		}
		switch {
		case suppress && counts:
			fmt.Fprintf(f, "%s\t%d\t%d\t%d\t%d\t%d\n",
				stat.chrom, stat.pos, stat.nAlleles, stat.nChr,
				firstCount, secondCount)
		case suppress:
			fmt.Fprintf(f, "%s\t%d\t%d\t%d\t%s\t%s\n",
				stat.chrom, stat.pos, stat.nAlleles, stat.nChr,
				formatCppDouble(firstFreq), formatCppDouble(secondFreq))
		case counts:
			fmt.Fprintf(f, "%s\t%d\t%d\t%d\t%s:%d\t%s:%d\n",
				stat.chrom, stat.pos, stat.nAlleles, stat.nChr,
				firstAllele, firstCount,
				secondAllele, secondCount)
		default:
			// Upstream emits frequencies via the default `ostream <<`
			// (six significant digits, trailing zeros stripped), not a
			// fixed %.6f — e.g. `A:0.5`, not `A:0.500000`. See
			// reference_code/vcftools/src/cpp/variant_file_output.cpp:135.
			fmt.Fprintf(f, "%s\t%d\t%d\t%d\t%s:%s\t%s:%s\n",
				stat.chrom, stat.pos, stat.nAlleles, stat.nChr,
				firstAllele, formatCppDouble(firstFreq),
				secondAllele, formatCppDouble(secondFreq))
		}
	}

	return nil
}

// outputSiteMeanDepth outputs mean depth per site (.ldepth.mean). Both
// MEAN_DEPTH and VAR_DEPTH follow upstream's bare `ostream <<` doubles
// (variant_file_output.cpp:3454-3458); when n==0 mean is 0/0 → -nan and
// variance is forced to -nan too (upstream divides by n-1 which on n==0
// yields -nan from the surrounding 0/0 mean²). When n==1 the n-1
// denominator in `(sumsq/n - mean²) * n/(n-1)` is 0 — that 0/0 also
// renders as -nan on glibc.
func (s *statistics) outputSiteMeanDepth(prefix string) error {
	f, err := iohelper.OpenWriter(prefix + ".ldepth.mean")
	if err != nil {
		return err
	}
	defer f.Close()

	fmt.Fprintln(f, "CHROM\tPOS\tMEAN_DEPTH\tVAR_DEPTH")

	for _, stat := range s.siteDepths {
		var mean, variance float64
		if stat.n == 0 {
			mean = math.NaN()
			variance = math.NaN()
		} else {
			mean = float64(stat.sumDepth) / float64(stat.n)
			if stat.n == 1 {
				variance = math.NaN()
			} else {
				variance = ((float64(stat.sumsqDepth) / float64(stat.n)) - mean*mean) * float64(stat.n) / float64(stat.n-1)
			}
		}
		fmt.Fprintf(f, "%s\t%d\t%s\t%s\n",
			stat.chrom, stat.pos,
			formatCppDouble(mean), formatCppDouble(variance))
	}

	return nil
}

// outputSiteDepth outputs sum depth per site (.ldepth). Both SUM_DEPTH
// and SUMSQ_DEPTH are unsigned ints in upstream (cpp:3434-3461); we
// emit them as decimal integers.
func (s *statistics) outputSiteDepth(prefix string) error {
	f, err := iohelper.OpenWriter(prefix + ".ldepth")
	if err != nil {
		return err
	}
	defer f.Close()

	fmt.Fprintln(f, "CHROM\tPOS\tSUM_DEPTH\tSUMSQ_DEPTH")

	for _, stat := range s.siteDepths {
		fmt.Fprintf(f, "%s\t%d\t%d\t%d\n", stat.chrom, stat.pos, stat.sumDepth, stat.sumsqDepth)
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
		if stat.qual < 0 {
			fmt.Fprintf(f, "%s\t%d\t.\n", stat.chrom, stat.pos)
		} else {
			// Upstream prints QUAL via the default ostream << (minimal
			// %g-style formatting): 30 not 30.0000. See
			// variant_file_output.cpp:1230-1234.
			fmt.Fprintf(f, "%s\t%d\t%s\n", stat.chrom, stat.pos, formatCppDouble(stat.qual))
		}
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
		// Upstream prints F_MISS = N_MISS / N_DATA via the default
		// ostream << (minimal %g-style formatting): 0 not 0.000000.
		// See variant_file_output.cpp:922-923.
		fMiss := 0.0
		if stat.nData != 0 {
			fMiss = float64(stat.nMiss) / float64(stat.nData)
		}
		fmt.Fprintf(f, "%s\t%d\t%d\t%d\t%d\t%s\n",
			stat.chrom, stat.pos, stat.nData, stat.nGenoFiltered,
			stat.nMiss, formatCppDouble(fMiss))
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
		stat.fMiss = float64(stat.nMissing) / float64(stat.nTotal)
		// Upstream prints F_MISS via the default ostream << (minimal
		// %g-style formatting): 0 not 0.000000. See
		// variant_file_output.cpp:1009-1013.
		fmt.Fprintf(f, "%s\t%d\t0\t%d\t%s\n",
			stat.name, stat.nTotal, stat.nMissing, formatCppDouble(stat.fMiss))
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

	// Upstream prints expected counts with `setprecision(2) << fixed`
	// (variant_file_output.cpp:368-369), then switches to default
	// precision + `scientific` for ChiSq_HWE/P_HWE/P_HET_DEFICIT/P_HET_EXCESS
	// (lines 370-372). Default `cout.precision()` is 6, so the scientific
	// output is %.6e.
	fmt.Fprintln(f, "CHR\tPOS\tOBS(HOM1/HET/HOM2)\tE(HOM1/HET/HOM2)\tChiSq_HWE\tP_HWE\tP_HET_DEFICIT\tP_HET_EXCESS")

	for _, stat := range s.siteHWE {
		fmt.Fprintf(f, "%s\t%d\t%d/%d/%d\t%.2f/%.2f/%.2f\t%.6e\t%.6e\t%.6e\t%.6e\n",
			stat.chrom, stat.pos,
			stat.obsHom1, stat.obsHet, stat.obsHom2,
			stat.expHom1, stat.expHet, stat.expHom2,
			stat.chiSq, stat.pValue, stat.pHetDeficit, stat.pHetExcess)
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
func (s *statistics) outputTsTvByBin(prefix string, binSize int) error {
	f, err := iohelper.OpenWriter(prefix + ".TsTv")
	if err != nil {
		return err
	}
	defer f.Close()

	// Header + per-(chrom, bin) row layout matches upstream
	// output_TsTv (variant_file_output.cpp:3057-3068): CHROM, BinStart
	// (= idx * binSize), SNP_count (Ts+Tv), and the C++-default-format
	// Ts/Tv ratio (0 when Tv==0). Bins are emitted dense per chromosome
	// in first-seen chromosome order.
	fmt.Fprintln(f, "CHROM\tBinStart\tSNP_count\tTs/Tv")

	for _, chrom := range s.tsTvBinChroms {
		row := s.tsTvByBin[chrom]
		for idx, stat := range row {
			ratio := 0.0
			if stat.tv != 0 {
				ratio = float64(stat.ts) / float64(stat.tv)
			}
			fmt.Fprintf(f, "%s\t%d\t%d\t%s\n",
				chrom, idx*binSize, stat.ts+stat.tv, formatCppDouble(ratio))
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
		fmt.Fprintf(f, "%s\t%d\t%s\n", stat.chrom, stat.pos, formatCppDouble(stat.pi))
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

	// Sort by individual name
	var names []string
	for name := range s.indvHet {
		names = append(names, name)
	}
	sort.Strings(names)

	for _, name := range names {
		stat := s.indvHet[name]
		// Upstream only prints individuals with at least one included
		// site (variant_file_output.cpp:268).
		if stat.nSites <= 0 {
			continue
		}

		// Method-of-moments inbreeding coefficient, following PLINK:
		//   F = (O - E) / (N - E)
		// where O = observed homozygotes, E = expected homozygotes
		// (summed per site as 1 - 2pq*N/(N-1)), N = sites included.
		// See variant_file_output.cpp:270.
		obsHom := stat.nObsHom
		expHom := stat.expHom
		fCoef := (float64(obsHom) - expHom) / (float64(stat.nSites) - expHom)

		// Upstream prints E(HOM) at fixed precision 1 and F at fixed
		// precision 5 (out.setf(ios::fixed); out.precision(1)/(5)).
		fmt.Fprintf(f, "%s\t%d\t%.1f\t%d\t%.5f\n",
			stat.name, obsHom, expHom, stat.nSites, fCoef)
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

// outputFilterSummary outputs FILTER tag summary
func (s *statistics) outputFilterSummary(prefix string) error {
	f, err := iohelper.OpenWriter(prefix + ".FILTER.summary")
	if err != nil {
		return err
	}
	defer f.Close()

	fmt.Fprintln(f, "FILTER\tN_SITES")

	// Sort by filter name
	var filters []string
	for filter := range s.filterCounts {
		filters = append(filters, filter)
	}
	sort.Strings(filters)

	for _, filter := range filters {
		count := s.filterCounts[filter]
		fmt.Fprintf(f, "%s\t%d\n", filter, count)
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

	fmt.Fprintln(f, "CHROM\tBIN_START\tBIN_END\tN_SNPs\tDENSITY")

	// Sort bins
	var bins []int
	for binIdx := range s.snpDensityBins {
		bins = append(bins, binIdx)
	}
	sort.Ints(bins)

	for _, binIdx := range bins {
		stat := s.snpDensityBins[binIdx]
		windowSize := float64(stat.binEnd - stat.binStart)
		if windowSize > 0 {
			stat.density = float64(stat.nSNPs) / windowSize
		}
		fmt.Fprintf(f, ".\t%d\t%d\t%d\t%.6f\n",
			stat.binStart, stat.binEnd, stat.nSNPs, stat.density)
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

	// Layout mirrors upstream output_TsTv_by_count
	// (variant_file_output.cpp:3220-3225): every count from 0 through
	// 2*N_kept_indv - 1 inclusive, with empty cells emitted as
	// "0\t0\t-nan" because upstream prints `double(Ts)/Tv` directly,
	// yielding glibc's signed-NaN literal for 0/0.
	fmt.Fprintln(f, "ALT_ALLELE_COUNT\tN_Ts\tN_Tv\tTs/Tv")

	nIndv := 0
	if s.header != nil {
		nIndv = len(s.header.Samples)
	}
	maxAC := 2 * nIndv
	for ac := 0; ac < maxAC; ac++ {
		stat := s.tsTvByCount[ac]
		ts, tv := 0, 0
		if stat != nil {
			ts, tv = stat.ts, stat.tv
		}
		ratio := math.NaN()
		if tv != 0 {
			ratio = float64(ts) / float64(tv)
		}
		fmt.Fprintf(f, "%d\t%d\t%d\t%s\n", ac, ts, tv, formatCppDouble(ratio))
	}

	return nil
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

	// Upstream's loop (variant_file_output.cpp:3260-3315) emits the
	// cumulative Ts/Tv sums STRICTLY BELOW and STRICTLY ABOVE each
	// threshold — values exactly equal to the threshold contribute to
	// neither side. The output formats every numeric column via the
	// default ostream << (`0`, not `0.0000`), with `0/0` yielding
	// glibc's signed-NaN literal (`-nan`).
	ratio := func(ts, tv int) float64 {
		if tv == 0 {
			return math.NaN()
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
			} else if st.qual > q {
				if st.isTs {
					tsGT++
				} else {
					tvGT++
				}
			}
		}
		fmt.Fprintf(f, "%s\t%d\t%d\t%s\t%d\t%d\t%s\n",
			formatCppDouble(q),
			tsLT, tvLT, formatCppDouble(ratio(tsLT, tvLT)),
			tsGT, tvGT, formatCppDouble(ratio(tsGT, tvGT)))
	}

	return nil
}

// outputDepth outputs per-individual mean read depth (.idepth).
// Upstream walks samples in VCF declaration order
// (variant_file_output.cpp:684-691) and emits MEAN_DEPTH with the
// default `ostream <<` formatter (%g-style 6-sig-figs). When an
// individual has no non-missing DP, mean = 0/0 → -nan.
func (s *statistics) outputDepth(prefix string) error {
	f, err := iohelper.OpenWriter(prefix + ".idepth")
	if err != nil {
		return err
	}
	defer f.Close()

	fmt.Fprintln(f, "INDV\tN_SITES\tMEAN_DEPTH")

	var names []string
	if s.header != nil {
		names = s.header.Samples
	} else {
		for name := range s.indvDepth {
			names = append(names, name)
		}
		sort.Strings(names)
	}

	for _, name := range names {
		stat := s.indvDepth[name]
		var mean float64
		var nSites int
		if stat == nil || stat.nSites == 0 {
			mean = math.NaN()
		} else {
			mean = float64(stat.sum) / float64(stat.nSites)
			nSites = stat.nSites
		}
		fmt.Fprintf(f, "%s\t%d\t%s\n", name, nSites, formatCppDouble(mean))
	}

	return nil
}

// windowPiAcc accumulates per-site nucleotide diversity within one window.
type windowPiAcc struct {
	winStart int
	nSites   int
	piSum    float64
}

// outputWindowedPi outputs nucleotide diversity summed over fixed-size windows
// (.windowed.pi). The PI column is the sum of per-site nucleotide diversity for
// variants falling in the window. If stepSize is zero or larger than windowSize
// the windows are non-overlapping.
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

	// Upstream emits an additional N_MONOMORPHIC column counting bases
	// inside the window that are monomorphic in the kept-individuals
	// subset. We don't track that yet — emit a literal 0 column so the
	// header and column count match upstream byte-for-byte. See
	// docs/PARITY_ROADMAP.md#vcftools.
	fmt.Fprintln(f, "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tN_MONOMORPHIC\tPI")

	var chromOrder []string
	windows := make(map[string]map[int]*windowPiAcc)

	for _, st := range s.sitePiValues {
		if _, seen := windows[st.chrom]; !seen {
			windows[st.chrom] = make(map[int]*windowPiAcc)
			chromOrder = append(chromOrder, st.chrom)
		}
		// Window starts are 1, 1+step, 1+2*step, ...; a variant at 1-based
		// position p belongs to every window [ws, ws+windowSize-1] containing p.
		p := st.pos
		kMax := (p - 1) / stepSize
		kMin := 0
		if p > windowSize {
			kMin = (p - windowSize + stepSize - 1) / stepSize
		}
		for k := kMin; k <= kMax; k++ {
			ws := 1 + k*stepSize
			acc := windows[st.chrom][ws]
			if acc == nil {
				acc = &windowPiAcc{winStart: ws}
				windows[st.chrom][ws] = acc
			}
			acc.nSites++
			acc.piSum += st.pi
		}
	}

	for _, chrom := range chromOrder {
		var starts []int
		for ws := range windows[chrom] {
			starts = append(starts, ws)
		}
		sort.Ints(starts)
		for _, ws := range starts {
			acc := windows[chrom][ws]
			fmt.Fprintf(f, "%s\t%d\t%d\t%d\t0\t%.6f\n", chrom, ws, ws+windowSize-1, acc.nSites, acc.piSum)
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

// chiSquareCDF approximates the chi-square CDF using gamma function
func chiSquareCDF(x, df float64) float64 {
	if x <= 0 {
		return 0
	}
	// For df=1, use a simple approximation
	// This is a simplified implementation
	// For production, use a proper statistical library
	return 1 - math.Exp(-x/2)
}

// formatCppDouble formats x the way upstream vcftools' C++ `ostream <<`
// does by default: six significant digits, fixed-or-scientific notation,
// trailing zeros stripped. Matches what Go's `strconv.FormatFloat(x, 'g',
// 6, 64)` produces for finite values.
//
// For NaN we emit the literal upstream produces on x86-64 glibc when an
// expression like `double(0)/0` is sent to an ostream: "-nan". Positive
// and negative infinity follow the C++ defaults ("inf" / "-inf").
func formatCppDouble(x float64) string {
	if math.IsNaN(x) {
		// Upstream's `double(Ts)/Tv` for Ts==Tv==0 renders as "-nan"
		// under x86-64 glibc (signed-NaN). Hardcode for byte parity.
		return "-nan"
	}
	if math.IsInf(x, 1) {
		return "inf"
	}
	if math.IsInf(x, -1) {
		return "-inf"
	}
	return strconv.FormatFloat(x, 'g', 6, 64)
}
