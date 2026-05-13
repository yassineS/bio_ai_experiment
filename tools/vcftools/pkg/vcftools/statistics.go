package vcftools

import (
	"fmt"
	"math"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/vcf"
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
	tsTvByBin     map[int]*tsTvBinStat
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
}

type siteDepthStat struct {
	chrom     string
	pos       int
	sumDepth  int
	meanDepth float64
}

type siteQualityStat struct {
	chrom string
	pos   int
	qual  float64
}

type siteMissingStat struct {
	chrom string
	pos   int
	fMiss float64
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
	pValue  float64
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

type tsTvBinStat struct {
	binStart int
	binEnd   int
	ts       int
	tv       int
	ratio    float64
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
	chrom  string
	pos    int
	allele string
}

// newStatistics creates a new statistics collector
func newStatistics(header *vcf.Header) *statistics {
	return &statistics{
		header:         header,
		indvMissing:    make(map[string]*indvMissingStat),
		indvHet:        make(map[string]*indvHetStat),
		indvDepth:      make(map[string]*indvDepthStat),
		tsTvByBin:      make(map[int]*tsTvBinStat),
		tsTvByCount:    make(map[int]*tsTvCountStat),
		snpDensityBins: make(map[int]*snpDensityStat),
		filterCounts:   make(map[string]int),
		indelLenHist:   make(map[int]int),
	}
}

// addVariant adds a variant to the statistics
func (s *statistics) addVariant(v *vcf.Variant, params *Params) {
	// Allele frequency
	if params.Freq || params.Counts {
		s.addFrequencyStat(v)
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

// addFrequencyStat adds allele frequency statistics
func (s *statistics) addFrequencyStat(v *vcf.Variant) {
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

	s.siteFrequencies = append(s.siteFrequencies, stat)
}

// addSiteDepthStat adds site depth statistics
func (s *statistics) addSiteDepthStat(v *vcf.Variant) {
	if len(v.Samples) == 0 {
		return
	}

	sumDepth := 0
	count := 0

	for _, sample := range v.Samples {
		dpStr, ok := sample.Data["DP"]
		if !ok {
			continue
		}

		var dp int
		_, err := fmt.Sscanf(dpStr, "%d", &dp)
		if err != nil {
			continue
		}

		sumDepth += dp
		count++
	}

	if count == 0 {
		return
	}

	stat := siteDepthStat{
		chrom:     v.Chrom,
		pos:       v.Pos,
		sumDepth:  sumDepth,
		meanDepth: float64(sumDepth) / float64(count),
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

	missing := 0
	for _, sample := range v.Samples {
		gt, ok := sample.Data["GT"]
		if !ok || gt == "." || gt == "./." || gt == ".|." {
			missing++
		}
	}

	stat := siteMissingStat{
		chrom: v.Chrom,
		pos:   v.Pos,
		fMiss: float64(missing) / float64(len(v.Samples)),
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

// addHWEStat adds Hardy-Weinberg equilibrium statistics
func (s *statistics) addHWEStat(v *vcf.Variant) {
	// Only for biallelic sites
	if len(v.Alt) != 1 {
		return
	}

	// Count genotypes
	hom1 := 0 // 0/0
	het := 0  // 0/1 or 1/0
	hom2 := 0 // 1/1

	for _, sample := range v.Samples {
		gt, ok := sample.Data["GT"]
		if !ok || strings.Contains(gt, ".") {
			continue
		}

		alleles := strings.FieldsFunc(gt, func(r rune) bool {
			return r == '/' || r == '|'
		})

		if len(alleles) != 2 {
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

	// P-value (1 degree of freedom)
	pValue := 1 - chiSquareCDF(chiSq, 1)

	stat := siteHWEStat{
		chrom:   v.Chrom,
		pos:     v.Pos,
		obsHom1: hom1,
		obsHet:  het,
		obsHom2: hom2,
		expHom1: expHom1,
		expHet:  expHet,
		expHom2: expHom2,
		chiSq:   chiSq,
		pValue:  pValue,
	}

	s.siteHWE = append(s.siteHWE, stat)
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

	isTransition := isTransitionSNP(ref, alt)
	if isTransition {
		s.transitions++
	} else {
		s.transversions++
	}

	// Add to bin if needed
	if binSize > 0 {
		binIdx := v.Pos / binSize
		if s.tsTvByBin[binIdx] == nil {
			s.tsTvByBin[binIdx] = &tsTvBinStat{
				binStart: binIdx * binSize,
				binEnd:   (binIdx + 1) * binSize,
			}
		}
		if isTransition {
			s.tsTvByBin[binIdx].ts++
		} else {
			s.tsTvByBin[binIdx].tv++
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
	for _, sample := range v.Samples {
		if s.indvHet[sample.Name] == nil {
			s.indvHet[sample.Name] = &indvHetStat{
				name: sample.Name,
			}
		}

		stat := s.indvHet[sample.Name]

		gt, ok := sample.Data["GT"]
		if !ok || strings.Contains(gt, ".") {
			continue
		}

		alleles := strings.FieldsFunc(gt, func(r rune) bool {
			return r == '/' || r == '|'
		})

		if len(alleles) != 2 {
			continue
		}

		stat.nTotal++
		if alleles[0] != alleles[1] {
			stat.nHet++
		} else if alleles[0] == "0" {
			stat.nHomRef++
		} else {
			stat.nHomAlt++
		}
	}
}

// addSingletonStat identifies singleton sites (alleles present in only one sample)
func (s *statistics) addSingletonStat(v *vcf.Variant) {
	if len(v.Samples) == 0 {
		return
	}

	// Count alleles
	alleleCounts := make(map[string]int)

	for _, sample := range v.Samples {
		gt, ok := sample.Data["GT"]
		if !ok || strings.Contains(gt, ".") {
			continue
		}

		alleles := strings.FieldsFunc(gt, func(r rune) bool {
			return r == '/' || r == '|'
		})

		for _, allele := range alleles {
			alleleCounts[allele]++
		}
	}

	// Check for singletons (alleles with count == 1)
	for allele, count := range alleleCounts {
		if count == 1 && allele != "0" { // Exclude reference allele
			s.singletonSites = append(s.singletonSites, singletonStat{
				chrom:  v.Chrom,
				pos:    v.Pos,
				allele: allele,
			})
		}
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

// outputFrequency outputs allele frequency statistics
func (s *statistics) outputFrequency(prefix string, counts bool) error {
	suffix := ".frq"
	if counts {
		suffix = ".frq.count"
	}

	f, err := iohelper.OpenWriter(prefix + suffix)
	if err != nil {
		return err
	}
	defer f.Close()

	// Header
	if counts {
		fmt.Fprintln(f, "CHROM\tPOS\tN_ALLELES\tN_CHR\t{REF:COUNT}\t{ALT:COUNT}")
	} else {
		fmt.Fprintln(f, "CHROM\tPOS\tN_ALLELES\tN_CHR\t{REF:FREQ}\t{ALT:FREQ}")
	}

	// Data
	for _, stat := range s.siteFrequencies {
		if counts {
			fmt.Fprintf(f, "%s\t%d\t%d\t%d\t{%s:%d}\t{%s:%d}\n",
				stat.chrom, stat.pos, stat.nAlleles, stat.nChr,
				stat.refAllele, stat.refCount,
				stat.altAllele, stat.altCount)
		} else {
			fmt.Fprintf(f, "%s\t%d\t%d\t%d\t{%s:%.6f}\t{%s:%.6f}\n",
				stat.chrom, stat.pos, stat.nAlleles, stat.nChr,
				stat.refAllele, stat.refFreq,
				stat.altAllele, stat.altFreq)
		}
	}

	return nil
}

// outputFrequency2 outputs alternative allele frequency format
func (s *statistics) outputFrequency2(prefix string) error {
	f, err := iohelper.OpenWriter(prefix + ".frq2")
	if err != nil {
		return err
	}
	defer f.Close()

	fmt.Fprintln(f, "CHROM\tPOS\tN_CHR\tREF_FREQ\tALT_FREQ")

	for _, stat := range s.siteFrequencies {
		fmt.Fprintf(f, "%s\t%d\t%d\t%.6f\t%.6f\n",
			stat.chrom, stat.pos, stat.nChr, stat.refFreq, stat.altFreq)
	}

	return nil
}

// outputCounts2 outputs alternative allele counts format
func (s *statistics) outputCounts2(prefix string) error {
	f, err := iohelper.OpenWriter(prefix + ".frq.count2")
	if err != nil {
		return err
	}
	defer f.Close()

	fmt.Fprintln(f, "CHROM\tPOS\tN_CHR\tREF_COUNT\tALT_COUNT")

	for _, stat := range s.siteFrequencies {
		fmt.Fprintf(f, "%s\t%d\t%d\t%d\t%d\n",
			stat.chrom, stat.pos, stat.nChr, stat.refCount, stat.altCount)
	}

	return nil
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
		// We don't calculate variance, so output 0
		fmt.Fprintf(f, "%s\t%d\t%.4f\t0\n", stat.chrom, stat.pos, stat.meanDepth)
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

	fmt.Fprintln(f, "CHROM\tPOS\tSUM_DEPTH")

	for _, stat := range s.siteDepths {
		fmt.Fprintf(f, "%s\t%d\t%d\n", stat.chrom, stat.pos, stat.sumDepth)
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
			fmt.Fprintf(f, "%s\t%d\t%.4f\n", stat.chrom, stat.pos, stat.qual)
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

	fmt.Fprintln(f, "CHROM\tPOS\tN_DATA\tN_GENOTYPE_FILTERED\tN_MISS\tF_MISS")

	for _, stat := range s.siteMissing {
		nData := len(s.header.Samples)
		nMiss := int(stat.fMiss * float64(nData))
		fmt.Fprintf(f, "%s\t%d\t%d\t0\t%d\t%.6f\n",
			stat.chrom, stat.pos, nData, nMiss, stat.fMiss)
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
		fmt.Fprintf(f, "%s\t%d\t0\t%d\t%.6f\n",
			stat.name, stat.nTotal, stat.nMissing, stat.fMiss)
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

	fmt.Fprintln(f, "CHROM\tPOS\tOBS(HOM1/HET/HOM2)\tE(HOM1/HET/HOM2)\tChiSq_HWE\tP_HWE")

	for _, stat := range s.siteHWE {
		fmt.Fprintf(f, "%s\t%d\t%d/%d/%d\t%.2f/%.2f/%.2f\t%.4f\t%.6g\n",
			stat.chrom, stat.pos,
			stat.obsHom1, stat.obsHet, stat.obsHom2,
			stat.expHom1, stat.expHet, stat.expHom2,
			stat.chiSq, stat.pValue)
	}

	return nil
}

// outputTsTvSummary outputs Ts/Tv summary
func (s *statistics) outputTsTvSummary(prefix string) error {
	f, err := iohelper.OpenWriter(prefix + ".TsTv.summary")
	if err != nil {
		return err
	}
	defer f.Close()

	ratio := 0.0
	if s.transversions > 0 {
		ratio = float64(s.transitions) / float64(s.transversions)
	}

	fmt.Fprintln(f, "Ts\tTv\tTs/Tv")
	fmt.Fprintf(f, "%d\t%d\t%.4f\n", s.transitions, s.transversions, ratio)

	return nil
}

// outputTsTvByBin outputs Ts/Tv by genomic bins
func (s *statistics) outputTsTvByBin(prefix string, binSize int) error {
	f, err := iohelper.OpenWriter(prefix + ".TsTv")
	if err != nil {
		return err
	}
	defer f.Close()

	fmt.Fprintln(f, "BIN_START\tBIN_END\tTs\tTv\tTs/Tv")

	// Sort bins
	var bins []int
	for binIdx := range s.tsTvByBin {
		bins = append(bins, binIdx)
	}
	sort.Ints(bins)

	for _, binIdx := range bins {
		stat := s.tsTvByBin[binIdx]
		ratio := 0.0
		if stat.tv > 0 {
			ratio = float64(stat.ts) / float64(stat.tv)
		}
		fmt.Fprintf(f, "%d\t%d\t%d\t%d\t%.4f\n",
			stat.binStart, stat.binEnd, stat.ts, stat.tv, ratio)
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
		fmt.Fprintf(f, "%s\t%d\t%.6f\n", stat.chrom, stat.pos, stat.pi)
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
		if stat.nTotal == 0 {
			continue
		}

		stat.hetRate = float64(stat.nHet) / float64(stat.nTotal)
		obsHom := stat.nHomRef + stat.nHomAlt

		// Expected homozygosity assuming HWE
		// This is a simplified calculation
		expHom := float64(stat.nTotal) * (1 - stat.hetRate)

		// Inbreeding coefficient F
		// F = (ExpectedHet - ObservedHet) / ExpectedHet
		expHet := float64(stat.nTotal) - expHom
		obsHet := float64(stat.nHet)
		f_coef := 0.0
		if expHet > 0 {
			f_coef = (expHet - obsHet) / expHet
		}

		fmt.Fprintf(f, "%s\t%d\t%.2f\t%d\t%.5f\n",
			stat.name, obsHom, expHom, stat.nTotal, f_coef)
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

	fmt.Fprintln(f, "CHROM\tPOS\tALLELE")

	for _, stat := range s.singletonSites {
		fmt.Fprintf(f, "%s\t%d\t%s\n", stat.chrom, stat.pos, stat.allele)
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

	fmt.Fprintln(f, "ALT_ALLELE_COUNT\tN_Ts\tN_Tv\tTs/Tv")

	var counts []int
	for ac := range s.tsTvByCount {
		counts = append(counts, ac)
	}
	sort.Ints(counts)

	for _, ac := range counts {
		stat := s.tsTvByCount[ac]
		ratio := 0.0
		if stat.tv > 0 {
			ratio = float64(stat.ts) / float64(stat.tv)
		}
		fmt.Fprintf(f, "%d\t%d\t%d\t%.4f\n", stat.altCount, stat.ts, stat.tv, ratio)
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
		fmt.Fprintf(f, "%s\t%d\t%.5f\n", stat.name, stat.nSites, mean)
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

	fmt.Fprintln(f, "CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tPI")

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
			fmt.Fprintf(f, "%s\t%d\t%d\t%d\t%.6f\n", chrom, ws, ws+windowSize-1, acc.nSites, acc.piSum)
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
