package vcftools

import (
	"fmt"
	"math"
	"sort"
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

	// Ts/Tv statistics
	transitions   int
	transversions int
	tsTvByBin     map[int]*tsTvBinStat

	// Phase 2: Population genetics statistics
	windowPiValues  []windowPiStat
	tajimaDValues   []tajimaDStat
	snpDensityBins  map[int]*snpDensityStat
	fstValues       []fstStat
	filterCounts    map[string]int
	singletonSites  []singletonStat
}

type siteFreqStat struct {
	chrom string
	pos   int
	nAlleles int
	nChr  int
	refAllele string
	altAllele string
	refFreq float64
	altFreq float64
	refCount int
	altCount int
}

type siteDepthStat struct {
	chrom string
	pos   int
	sumDepth int
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
	chrom string
	pos   int
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
	name        string
	nMissing    int
	nTotal      int
	fMiss       float64
}

type tsTvBinStat struct {
	binStart int
	binEnd   int
	ts       int
	tv       int
	ratio    float64
}

type indvHetStat struct {
	name        string
	nHet        int
	nHomAlt     int
	nHomRef     int
	nTotal      int
	hetRate     float64
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
	allele string
}

// newStatistics creates a new statistics collector
func newStatistics(header *vcf.Header) *statistics {
	return &statistics{
		header:         header,
		indvMissing:    make(map[string]*indvMissingStat),
		indvHet:        make(map[string]*indvHetStat),
		tsTvByBin:      make(map[int]*tsTvBinStat),
		snpDensityBins: make(map[int]*snpDensityStat),
		filterCounts:   make(map[string]int),
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

	// Site pi (nucleotide diversity)
	if params.SitePi {
		s.addSitePiStat(v)
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

	isTransition := false
	if (ref == "A" && alt == "G") || (ref == "G" && alt == "A") ||
		(ref == "C" && alt == "T") || (ref == "T" && alt == "C") {
		isTransition = true
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

// addSitePiStat adds nucleotide diversity statistics
func (s *statistics) addSitePiStat(v *vcf.Variant) {
	if len(v.Samples) == 0 {
		return
	}

	// Calculate pairwise differences
	var genotypes []string
	for _, sample := range v.Samples {
		gt, ok := sample.Data["GT"]
		if ok && !strings.Contains(gt, ".") {
			genotypes = append(genotypes, gt)
		}
	}

	if len(genotypes) < 2 {
		return
	}

	// Count pairwise differences
	differences := 0
	comparisons := 0

	for i := 0; i < len(genotypes); i++ {
		for j := i + 1; j < len(genotypes); j++ {
			diff := genotypeDistance(genotypes[i], genotypes[j])
			differences += diff
			comparisons++
		}
	}

	pi := 0.0
	if comparisons > 0 {
		pi = float64(differences) / float64(comparisons)
	}

	stat := sitePiStat{
		chrom: v.Chrom,
		pos:   v.Pos,
		pi:    pi,
	}

	s.sitePiValues = append(s.sitePiValues, stat)
}

// genotypeDistance calculates the number of differences between two genotypes
func genotypeDistance(gt1, gt2 string) int {
	alleles1 := strings.FieldsFunc(gt1, func(r rune) bool {
		return r == '/' || r == '|'
	})
	alleles2 := strings.FieldsFunc(gt2, func(r rune) bool {
		return r == '/' || r == '|'
	})

	if len(alleles1) != len(alleles2) {
		return 0
	}

	diff := 0
	for i := range alleles1 {
		if alleles1[i] != alleles2[i] {
			diff++
		}
	}
	return diff
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
