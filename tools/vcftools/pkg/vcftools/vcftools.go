// Package vcftools provides utilities for working with VCF files
package vcftools

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/vcf"
)

// Params holds all parameters for vcftools operations
type Params struct {
	// Output options
	OutPrefix     string
	UseStdout     bool
	Recode        bool
	RecodeInfoAll bool

	// Position filtering
	Chr                  string
	NotChr               string
	FromBp               int
	ToBp                 int
	PositionsFile        string
	ExcludePositionsFile string

	// Variant type filtering
	KeepOnlyIndels bool
	RemoveIndels   bool
	MinAlleles     int
	MaxAlleles     int

	// Quality filtering
	MinQ              float64
	RemoveFilteredAll bool

	// Allele frequency filtering
	Maf    float64
	MaxMaf float64
	Mac    int
	MaxMac int

	// Genotype filtering
	MaxMissing float64
	MinMeanDP  float64
	MaxMeanDP  float64

	// Statistics output
	Freq          bool
	Counts        bool
	Freq2         bool
	Counts2       bool
	Depth         bool
	SiteDepth     bool
	SiteMeanDepth bool
	SiteQuality   bool
	MissingIndv   bool
	MissingSite   bool
	Hardy         bool
	TsTvSummary   bool
	TsTvBinSize   int
	TsTvByCount   bool
	TsTvByQual    bool
	SitePi        bool
	Het           bool
	Singletons    bool
	HistIndelLen  bool
	GenoDepth     bool
	
	// Phase 2: Population genetics statistics
	WindowPi       int
	WindowPiStep   int
	TajimaD        int
	SNPDensity     int
	WeirFstPop     []string
	FstWindowSize  int
	FstWindowStep  int
	FilterSummary  bool
	
	// Phase 4: Format conversions
	Output012      bool
	OutputPlink    bool
	OutputPlinkTped bool
	ChromMap       string

	// Sample filtering
	IndvList       []string
	RemoveIndvList []string
	KeepFile       string
	RemoveFile     string
}

// positionSet represents a set of positions to include/exclude
type positionSet map[string]map[int]bool

// Run executes vcftools with the given parameters
func Run(input io.Reader, params *Params) error {
	// Read VCF
	reader := vcf.NewReader(input)
	header, err := reader.ReadHeader()
	if err != nil {
		return fmt.Errorf("reading VCF header: %w", err)
	}

	// Load position filters if needed
	var includePositions, excludePositions positionSet
	if params.PositionsFile != "" {
		includePositions, err = loadPositions(params.PositionsFile)
		if err != nil {
			return fmt.Errorf("loading positions file: %w", err)
		}
	}
	if params.ExcludePositionsFile != "" {
		excludePositions, err = loadPositions(params.ExcludePositionsFile)
		if err != nil {
			return fmt.Errorf("loading exclude positions file: %w", err)
		}
	}

	// Build sample filter set
	keepSamples, err := buildSampleFilter(header, params)
	if err != nil {
		return fmt.Errorf("building sample filter: %w", err)
	}

	// Filter header samples if needed
	filteredHeader := filterHeaderSamples(header, keepSamples)

	// Initialize statistics
	stats := newStatistics(filteredHeader)

	// Set up output writer for recode
	var recodeWriter *vcf.Writer
	if params.Recode {
		var w io.Writer
		if params.UseStdout {
			w = os.Stdout
		} else {
			outFile := params.OutPrefix + ".recode.vcf"
			f, err := iohelper.OpenWriter(outFile)
			if err != nil {
				return fmt.Errorf("opening output file: %w", err)
			}
			defer f.Close()
			w = f
		}
		recodeWriter = vcf.NewWriter(w, filteredHeader)
		if err := recodeWriter.WriteHeader(); err != nil {
			return fmt.Errorf("writing VCF header: %w", err)
		}
	}

	// Process variants
	keptSites := 0
	totalSites := 0
	var allVariants []*vcf.Variant // For format conversions that need all data

	for {
		variant, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return fmt.Errorf("reading variant: %w", err)
		}

		totalSites++

		// Apply filters
		if !passFilters(variant, params, includePositions, excludePositions) {
			continue
		}

		// Filter samples
		filteredVariant := filterVariantSamples(variant, keepSamples)

		// Update statistics
		stats.addVariant(filteredVariant, params)
		
		// Collect variants for format conversions
		if params.Output012 || params.OutputPlink || params.OutputPlinkTped {
			allVariants = append(allVariants, filteredVariant)
		}

		// Write to output if recoding
		if params.Recode {
			if params.RecodeInfoAll {
				// Keep all INFO fields
				if err := recodeWriter.Write(filteredVariant); err != nil {
					return fmt.Errorf("writing variant: %w", err)
				}
			} else {
				// Keep only essential INFO fields
				filtered := &vcf.Variant{
					Chrom:   filteredVariant.Chrom,
					Pos:     filteredVariant.Pos,
					ID:      filteredVariant.ID,
					Ref:     filteredVariant.Ref,
					Alt:     filteredVariant.Alt,
					Qual:    filteredVariant.Qual,
					Filter:  filteredVariant.Filter,
					Info:    make(map[string]string),
					Format:  filteredVariant.Format,
					Samples: filteredVariant.Samples,
				}
				if err := recodeWriter.Write(filtered); err != nil {
					return fmt.Errorf("writing variant: %w", err)
				}
			}
		}

		keptSites++
	}

	// Flush recode writer if needed
	if recodeWriter != nil {
		if err := recodeWriter.Flush(); err != nil {
			return fmt.Errorf("flushing output: %w", err)
		}
	}

	// Print summary to stderr
	fmt.Fprintf(os.Stderr, "\nAfter filtering, kept %d out of a possible %d Sites\n", keptSites, totalSites)

	// Output statistics
	if err := outputStatistics(stats, params); err != nil {
		return fmt.Errorf("outputting statistics: %w", err)
	}
	
	// Output format conversions
	if err := outputFormatConversions(allVariants, filteredHeader, params); err != nil {
		return fmt.Errorf("outputting format conversions: %w", err)
	}

	return nil
}

// passFilters checks if a variant passes all filters
func passFilters(v *vcf.Variant, params *Params, includePos, excludePos positionSet) bool {
	// Position filters
	if params.Chr != "" && v.Chrom != params.Chr {
		return false
	}
	if params.NotChr != "" && v.Chrom == params.NotChr {
		return false
	}
	if params.FromBp > 0 && v.Pos < params.FromBp {
		return false
	}
	if params.ToBp > 0 && v.Pos > params.ToBp {
		return false
	}

	// Position include/exclude
	if includePos != nil {
		if chromPos, ok := includePos[v.Chrom]; ok {
			if !chromPos[v.Pos] {
				return false
			}
		} else {
			return false
		}
	}
	if excludePos != nil {
		if chromPos, ok := excludePos[v.Chrom]; ok {
			if chromPos[v.Pos] {
				return false
			}
		}
	}

	// Variant type filters
	isIndel := isIndelVariant(v)
	if params.KeepOnlyIndels && !isIndel {
		return false
	}
	if params.RemoveIndels && isIndel {
		return false
	}

	// Allele count filter
	numAlleles := len(v.Alt) + 1 // +1 for reference
	if params.MinAlleles > 0 && numAlleles < params.MinAlleles {
		return false
	}
	if params.MaxAlleles > 0 && numAlleles > params.MaxAlleles {
		return false
	}

	// Quality filter
	if params.MinQ > 0 && (v.Qual < 0 || v.Qual < params.MinQ) {
		return false
	}

	// Filter flag
	if params.RemoveFilteredAll {
		if len(v.Filter) == 0 || (len(v.Filter) == 1 && v.Filter[0] != "PASS") {
			if len(v.Filter) == 0 || v.Filter[0] != "PASS" {
				return false
			}
		}
	}

	// Allele frequency filters
	if params.Maf > 0 || params.MaxMaf > 0 || params.Mac > 0 || params.MaxMac > 0 {
		maf, mac := calculateMAF(v)
		if params.Maf > 0 && maf < params.Maf {
			return false
		}
		if params.MaxMaf > 0 && maf > params.MaxMaf {
			return false
		}
		if params.Mac > 0 && mac < params.Mac {
			return false
		}
		if params.MaxMac > 0 && mac > params.MaxMac {
			return false
		}
	}

	// Genotype filters
	if params.MaxMissing < 1 {
		missingRate := calculateMissingRate(v)
		if missingRate > (1 - params.MaxMissing) {
			return false
		}
	}

	if params.MinMeanDP > 0 || params.MaxMeanDP > 0 {
		meanDP := calculateMeanDepth(v)
		if params.MinMeanDP > 0 && meanDP < params.MinMeanDP {
			return false
		}
		if params.MaxMeanDP > 0 && meanDP > params.MaxMeanDP {
			return false
		}
	}

	return true
}

// isIndelVariant checks if a variant is an indel
func isIndelVariant(v *vcf.Variant) bool {
	refLen := len(v.Ref)
	for _, alt := range v.Alt {
		if len(alt) != refLen {
			return true
		}
	}
	return false
}

// calculateMAF calculates minor allele frequency and count
func calculateMAF(v *vcf.Variant) (maf float64, mac int) {
	if len(v.Samples) == 0 {
		return 0, 0
	}

	// Count alleles
	alleleCounts := make(map[string]int)
	totalAlleles := 0

	for _, sample := range v.Samples {
		gt, ok := sample.Data["GT"]
		if !ok || gt == "." || gt == "./." || gt == ".|." {
			continue
		}

		// Parse genotype
		alleles := strings.FieldsFunc(gt, func(r rune) bool {
			return r == '/' || r == '|'
		})

		for _, allele := range alleles {
			if allele != "." {
				alleleCounts[allele]++
				totalAlleles++
			}
		}
	}

	if totalAlleles == 0 {
		return 0, 0
	}

	// Find minor allele (least frequent)
	minCount := totalAlleles
	for _, count := range alleleCounts {
		if count < minCount {
			minCount = count
		}
	}

	// If all alleles are the same, minor allele count is 0
	if len(alleleCounts) == 1 {
		return 0, 0
	}

	mac = minCount
	maf = float64(minCount) / float64(totalAlleles)
	return maf, mac
}

// calculateMissingRate calculates the proportion of missing genotypes
func calculateMissingRate(v *vcf.Variant) float64 {
	if len(v.Samples) == 0 {
		return 0
	}

	missing := 0
	for _, sample := range v.Samples {
		gt, ok := sample.Data["GT"]
		if !ok || gt == "." || gt == "./." || gt == ".|." {
			missing++
		}
	}

	return float64(missing) / float64(len(v.Samples))
}

// calculateMeanDepth calculates mean depth across samples
func calculateMeanDepth(v *vcf.Variant) float64 {
	if len(v.Samples) == 0 {
		return 0
	}

	totalDepth := 0
	count := 0

	for _, sample := range v.Samples {
		dpStr, ok := sample.Data["DP"]
		if !ok {
			continue
		}
		dp, err := strconv.Atoi(dpStr)
		if err != nil {
			continue
		}
		totalDepth += dp
		count++
	}

	if count == 0 {
		return 0
	}

	return float64(totalDepth) / float64(count)
}

// loadPositions loads a positions file
func loadPositions(filename string) (positionSet, error) {
	f, err := iohelper.OpenReader(filename)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	positions := make(positionSet)
	scanner := bufio.NewScanner(f)

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}

		fields := strings.Fields(line)
		if len(fields) < 2 {
			continue
		}

		chrom := fields[0]
		pos, err := strconv.Atoi(fields[1])
		if err != nil {
			continue
		}

		if positions[chrom] == nil {
			positions[chrom] = make(map[int]bool)
		}
		positions[chrom][pos] = true
	}

	return positions, scanner.Err()
}

// buildSampleFilter builds a set of samples to keep
func buildSampleFilter(header *vcf.Header, params *Params) (map[string]bool, error) {
	// If no sample filters, keep all
	if len(params.IndvList) == 0 && len(params.RemoveIndvList) == 0 &&
		params.KeepFile == "" && params.RemoveFile == "" {
		return nil, nil
	}

	// Start with all samples
	keep := make(map[string]bool)
	for _, sample := range header.Samples {
		keep[sample] = true
	}

	// Apply keep file
	if params.KeepFile != "" {
		samples, err := loadSampleFile(params.KeepFile)
		if err != nil {
			return nil, err
		}
		newKeep := make(map[string]bool)
		for _, sample := range samples {
			if keep[sample] {
				newKeep[sample] = true
			}
		}
		keep = newKeep
	}

	// Apply keep list
	if len(params.IndvList) > 0 {
		newKeep := make(map[string]bool)
		for _, sample := range params.IndvList {
			if keep[sample] {
				newKeep[sample] = true
			}
		}
		keep = newKeep
	}

	// Apply remove file
	if params.RemoveFile != "" {
		samples, err := loadSampleFile(params.RemoveFile)
		if err != nil {
			return nil, err
		}
		for _, sample := range samples {
			delete(keep, sample)
		}
	}

	// Apply remove list
	for _, sample := range params.RemoveIndvList {
		delete(keep, sample)
	}

	return keep, nil
}

// loadSampleFile loads a file with one sample name per line
func loadSampleFile(filename string) ([]string, error) {
	f, err := iohelper.OpenReader(filename)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	var samples []string
	scanner := bufio.NewScanner(f)

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line != "" && !strings.HasPrefix(line, "#") {
			samples = append(samples, line)
		}
	}

	return samples, scanner.Err()
}

// filterHeaderSamples filters header samples based on keep set
func filterHeaderSamples(header *vcf.Header, keepSamples map[string]bool) *vcf.Header {
	if keepSamples == nil {
		return header
	}

	filtered := &vcf.Header{
		MetaInfo: header.MetaInfo,
		Samples:  []string{},
	}

	for _, sample := range header.Samples {
		if keepSamples[sample] {
			filtered.Samples = append(filtered.Samples, sample)
		}
	}

	return filtered
}

// filterVariantSamples filters variant samples based on keep set
func filterVariantSamples(v *vcf.Variant, keepSamples map[string]bool) *vcf.Variant {
	if keepSamples == nil {
		return v
	}

	filtered := &vcf.Variant{
		Chrom:  v.Chrom,
		Pos:    v.Pos,
		ID:     v.ID,
		Ref:    v.Ref,
		Alt:    v.Alt,
		Qual:   v.Qual,
		Filter: v.Filter,
		Info:   v.Info,
		Format: v.Format,
	}

	for _, sample := range v.Samples {
		if keepSamples[sample.Name] {
			filtered.Samples = append(filtered.Samples, sample)
		}
	}

	return filtered
}

// outputStatistics outputs all requested statistics
func outputStatistics(stats *statistics, params *Params) error {
	if params.Freq {
		if err := stats.outputFrequency(params.OutPrefix, false); err != nil {
			return err
		}
	}

	if params.Counts {
		if err := stats.outputFrequency(params.OutPrefix, true); err != nil {
			return err
		}
	}

	if params.SiteMeanDepth {
		if err := stats.outputSiteMeanDepth(params.OutPrefix); err != nil {
			return err
		}
	}

	if params.SiteDepth {
		if err := stats.outputSiteDepth(params.OutPrefix); err != nil {
			return err
		}
	}

	if params.SiteQuality {
		if err := stats.outputSiteQuality(params.OutPrefix); err != nil {
			return err
		}
	}

	if params.MissingSite {
		if err := stats.outputMissingSite(params.OutPrefix); err != nil {
			return err
		}
	}

	if params.MissingIndv {
		if err := stats.outputMissingIndv(params.OutPrefix); err != nil {
			return err
		}
	}

	if params.Hardy {
		if err := stats.outputHWE(params.OutPrefix); err != nil {
			return err
		}
	}

	if params.TsTvSummary {
		if err := stats.outputTsTvSummary(params.OutPrefix); err != nil {
			return err
		}
	}

	if params.TsTvBinSize > 0 {
		if err := stats.outputTsTvByBin(params.OutPrefix, params.TsTvBinSize); err != nil {
			return err
		}
	}

	if params.SitePi {
		if err := stats.outputSitePi(params.OutPrefix); err != nil {
			return err
		}
	}
	
	// Phase 2: Population genetics statistics
	
	if params.Het {
		if err := stats.outputHet(params.OutPrefix); err != nil {
			return err
		}
	}
	
	if params.Singletons {
		if err := stats.outputSingletons(params.OutPrefix); err != nil {
			return err
		}
	}
	
	if params.FilterSummary {
		if err := stats.outputFilterSummary(params.OutPrefix); err != nil {
			return err
		}
	}
	
	if params.SNPDensity > 0 {
		if err := stats.outputSNPDensity(params.OutPrefix, params.SNPDensity); err != nil {
			return err
		}
	}

	return nil
}

// outputFormatConversions outputs requested format conversions
func outputFormatConversions(variants []*vcf.Variant, header *vcf.Header, params *Params) error {
	if params.Output012 {
		if err := output012Matrix(variants, header, params.OutPrefix); err != nil {
			return err
		}
	}
	
	if params.OutputPlink || params.OutputPlinkTped {
		// Load chromosome map if provided
		chromMap, err := loadChromMap(params.ChromMap)
		if err != nil {
			return fmt.Errorf("loading chromosome map: %w", err)
		}
		
		if params.OutputPlink {
			if err := outputPlink(variants, header, params.OutPrefix, chromMap); err != nil {
				return err
			}
		}
		
		if params.OutputPlinkTped {
			if err := outputPlinkTped(variants, header, params.OutPrefix, chromMap); err != nil {
				return err
			}
		}
	}
	
	return nil
}

// Helper function to get output file path
func getOutputPath(prefix, suffix string) string {
	return filepath.Join(".", prefix+suffix)
}
