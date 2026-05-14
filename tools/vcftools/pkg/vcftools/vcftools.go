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

	// SNP ID filtering
	SNP         string
	SNPs        string
	ExcludeSNP  string
	ExcludeSNPs string
	Thin        int

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
	MinDP      int
	MaxDP      int
	MinGQ      int

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
	WindowPi      int
	WindowPiStep  int
	TajimaD       int
	SNPDensity    int
	WeirFstPop    []string
	FstWindowSize int
	FstWindowStep int
	FilterSummary bool

	// Phase 4: Format conversions
	Output012       bool
	OutputPlink     bool
	OutputPlinkTped bool
	ChromMap        string

	// Sample filtering
	IndvList       []string
	RemoveIndvList []string
	KeepFile       string
	RemoveFile     string

	// Linkage disequilibrium analysis (--geno-r2 / --hap-r2 family).
	// GenoR2 enables --geno-r2 output (<prefix>.geno.ld). HapR2 enables
	// --hap-r2 output (<prefix>.hap.ld). GenoR2Positions / HapR2Positions
	// supply chrom/pos files that restrict pairs to those where at least one
	// endpoint is listed (analogous to upstream --geno-r2-positions /
	// --hap-r2-positions). The four window fields bound the pairwise SNP /
	// bp distance: zero means "no bound" for the maxima, zero for the minima
	// means "no minimum required". MinR2 thresholds the emitted r² (default
	// 0 = emit all pairs).
	GenoR2          bool
	HapR2           bool
	GenoR2Positions string
	HapR2Positions  string
	LDWindow        int
	LDWindowBp      int
	LDWindowMin     int
	LDWindowBpMin   int
	MinR2           float64

	// BED-based site filtering. Bed keeps only sites whose 1-based POS lies
	// inside any interval (0-based half-open) in the file. ExcludeBed is the
	// inverse. Both compose with other position/quality filters.
	Bed        string
	ExcludeBed string

	// VCF comparison (--diff family). Diff names the second VCF to compare
	// against; the boolean flags request individual output files. See
	// diff.go for the column layout of each output.
	Diff                string
	DiffSite            bool
	DiffIndv            bool
	DiffSiteDiscordance bool
	DiffIndvDiscordance bool

	// BEAGLE genotype-likelihood output. BEAGLEGL writes log10-scale GL
	// triplets derived from FORMAT/PL; BEAGLEPL writes the raw PL triplets.
	// Both are biallelic-SNP only.
	BEAGLEGL bool
	BEAGLEPL bool

	// Inter-chromosomal LD outputs. InterchromGenoR2 and InterchromHapR2
	// emit `<prefix>.interchrom.geno.ld` / `<prefix>.interchrom.hap.ld`
	// (only cross-chromosome pairs). GenoChiSq emits `<prefix>.geno.chisq`
	// (per-pair Pearson chi-square test across all pairs, same- and
	// cross-chromosome).
	InterchromGenoR2 bool
	InterchromHapR2  bool
	GenoChiSq        bool

	// Relatedness statistics. Relatedness enables <prefix>.relatedness
	// (Yang 2010 unadjusted A_jk). Relatedness2 enables
	// <prefix>.relatedness2 (KING-robust kinship; Manichaikul 2010).
	Relatedness  bool
	Relatedness2 bool

	// PhasedBlocks enables <prefix>.blocks reporting per-individual
	// contiguous runs of phased ("a|b") diploid genotypes.
	PhasedBlocks bool

	// Runs of homozygosity. LROH enables <prefix>.LROH. LROHMinVariants is
	// the minimum number of consecutive homozygous variants for a run to
	// be emitted (default 10 when LROH is true and the value is zero).
	LROH            bool
	LROHMinVariants int

	// Filter-name include/exclude (operate on the FILTER column). Both are
	// comma-separated lists. RemoveFiltered drops sites whose FILTER lists
	// any of the named tags; KeepFiltered keeps only sites that list at
	// least one of the named tags.
	RemoveFiltered string
	KeepFiltered   string

	// INFO tag selection for --recode output. Both are comma-separated lists
	// applied during recoding only. KeepINFO restricts the output INFO map
	// to the listed tags; RemoveINFO strips the listed tags. (--recode-INFO-all
	// already preserves everything; --keep-INFO and --remove-INFO compose
	// after it.)
	KeepINFO   string
	RemoveINFO string

	// --get-INFO TAG[,TAG]... extracts the named INFO tags as a TSV file
	// <prefix>.INFO with columns CHROM POS REF ALT <tags...>. The flag is
	// comma-separated; the upstream CLI accepts the same value-style.
	GetINFO string
}

// positionSet represents a set of positions to include/exclude
type positionSet map[string]map[int]bool

// Run executes vcftools with the given parameters
func Run(input io.Reader, params *Params) error {
	// Reject requested features that this port does not implement yet, instead
	// of silently producing no output.
	if err := checkUnsupported(params); err != nil {
		return err
	}

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

	// Load BED-based filters (--bed / --exclude-bed). Both are optional;
	// supplying both composes them (a site must pass include AND not be
	// excluded).
	var includeBed, excludeBed *bedRegions
	if params.Bed != "" {
		includeBed, err = loadBedRegions(params.Bed)
		if err != nil {
			return fmt.Errorf("loading --bed file: %w", err)
		}
	}
	if params.ExcludeBed != "" {
		excludeBed, err = loadBedRegions(params.ExcludeBed)
		if err != nil {
			return fmt.Errorf("loading --exclude-bed file: %w", err)
		}
	}

	// Load --weir-fst-pop population files if requested. We validate here so
	// that errors (missing file, sample appearing in multiple populations,
	// fewer than 2 populations) are surfaced before we start streaming the
	// VCF.
	var weirFstPops [][]string
	if len(params.WeirFstPop) > 0 {
		weirFstPops, err = loadPopulationFiles(params.WeirFstPop)
		if err != nil {
			return fmt.Errorf("loading --weir-fst-pop files: %w", err)
		}
	}

	// Load SNP ID filters
	var includeSNPs, excludeSNPs map[string]bool
	if params.SNP != "" {
		includeSNPs = make(map[string]bool)
		includeSNPs[params.SNP] = true
	}
	if params.SNPs != "" {
		includeSNPs, err = loadSNPIDs(params.SNPs)
		if err != nil {
			return fmt.Errorf("loading SNPs file: %w", err)
		}
	}
	if params.ExcludeSNP != "" {
		excludeSNPs = make(map[string]bool)
		excludeSNPs[params.ExcludeSNP] = true
	}
	if params.ExcludeSNPs != "" {
		excludeSNPs, err = loadSNPIDs(params.ExcludeSNPs)
		if err != nil {
			return fmt.Errorf("loading exclude SNPs file: %w", err)
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
	if len(weirFstPops) >= 2 {
		stats.weirFst = newWeirFstAccumulator(weirFstPops)
	}

	// Initialise LD runner (no-op when no LD flag is set).
	var ldRun *ldRunner
	if params.GenoR2 || params.HapR2 || params.GenoR2Positions != "" || params.HapR2Positions != "" {
		ldRun, err = newLDRunner(params)
		if err != nil {
			return fmt.Errorf("initialising LD analysis: %w", err)
		}
	}

	// Initialise --diff runner (no-op when --diff isn't set or no diff
	// sub-output is requested).
	diffRun, err := newDiffRunner(params, filteredHeader.Samples)
	if err != nil {
		return fmt.Errorf("initialising --diff analysis: %w", err)
	}

	// Initialise BEAGLE output writers. They're created lazily here so a
	// failure to open the file surfaces before we stream any variants.
	var beagleGL, beaglePL *beagleWriter
	if params.BEAGLEGL {
		beagleGL, err = newBEAGLEWriter(params.OutPrefix, beagleGLMode())
		if err != nil {
			return fmt.Errorf("initialising --BEAGLE-GL: %w", err)
		}
	}
	if params.BEAGLEPL {
		beaglePL, err = newBEAGLEWriter(params.OutPrefix, beaglePLMode())
		if err != nil {
			return fmt.Errorf("initialising --BEAGLE-PL: %w", err)
		}
	}

	// Inter-chromosomal LD / chi-square buffer (all-pairs after streaming).
	var interLD *interchromLDRunner
	if params.InterchromGenoR2 || params.InterchromHapR2 || params.GenoChiSq {
		interLD, err = newInterchromLDRunner(params)
		if err != nil {
			return fmt.Errorf("initialising interchrom LD: %w", err)
		}
	}

	// --relatedness accumulator.
	var rel *relatednessRunner
	if params.Relatedness {
		rel = newRelatednessRunner(filteredHeader.Samples)
	}

	// --relatedness2 accumulator.
	var rel2 *relatedness2Runner
	if params.Relatedness2 {
		rel2 = newRelatedness2Runner(filteredHeader.Samples)
	}

	// --LROH runner.
	var lroh *lrohRunner
	if params.LROH {
		lroh = newLROHRunner(filteredHeader.Samples, params.LROHMinVariants)
	}

	// --phased-blocks runner.
	var phasedBlocks *phasedBlockRunner
	if params.PhasedBlocks {
		phasedBlocks = newPhasedBlockRunner(filteredHeader.Samples)
	}

	// --get-INFO writer.
	var getInfo *getInfoRunner
	if params.GetINFO != "" {
		tags := splitCSV(params.GetINFO)
		getInfo, err = newGetInfoRunner(params.OutPrefix, tags)
		if err != nil {
			return fmt.Errorf("initialising --get-INFO: %w", err)
		}
	}

	// Pre-parse filter-name sets and INFO-tag sets so we don't re-tokenise
	// every line in the hot path.
	removeFilteredSet := parseFilterList(params.RemoveFiltered)
	keepFilteredSet := parseFilterList(params.KeepFiltered)
	keepInfoSet := parseInfoTagList(params.KeepINFO)
	removeInfoSet := parseInfoTagList(params.RemoveINFO)

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
	thinCounter := 0
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

		// Apply thinning filter
		if params.Thin > 0 {
			thinCounter++
			if thinCounter%params.Thin != 0 {
				continue
			}
		}

		// Apply filters
		if !passFilters(variant, params, includePositions, excludePositions, includeSNPs, excludeSNPs, includeBed, excludeBed) {
			continue
		}

		// --remove-filtered <names>: drop sites listing any of the named
		// FILTERs. --keep-filtered <names>: keep only sites listing at
		// least one of the named FILTERs. Both compose with
		// --remove-filtered-all (which is the union of all non-PASS sites).
		if !passRemoveFilteredNames(variant, removeFilteredSet) {
			continue
		}
		if !passKeepFilteredNames(variant, keepFilteredSet) {
			continue
		}

		// Filter samples
		filteredVariant := filterVariantSamples(variant, keepSamples)

		// Apply genotype-level filters
		filteredVariant = filterGenotypes(filteredVariant, params)

		// Update statistics
		stats.addVariant(filteredVariant, params)

		// Feed LD runner (writes pairwise output incrementally).
		if ldRun != nil {
			ldRun.addVariant(filteredVariant)
		}

		// Feed --diff runner.
		if diffRun != nil {
			if err := diffRun.addVariant(filteredVariant); err != nil {
				return fmt.Errorf("writing diff output: %w", err)
			}
		}

		// Emit BEAGLE rows (header is emitted lazily on the first call).
		if beagleGL != nil {
			if err := beagleGL.write(filteredVariant, filteredHeader.Samples); err != nil {
				return fmt.Errorf("writing BEAGLE-GL output: %w", err)
			}
		}
		if beaglePL != nil {
			if err := beaglePL.write(filteredVariant, filteredHeader.Samples); err != nil {
				return fmt.Errorf("writing BEAGLE-PL output: %w", err)
			}
		}

		// Inter-chromosomal LD: buffer for end-of-stream pair emission.
		if interLD != nil {
			interLD.addVariant(filteredVariant)
		}

		// Per-variant relatedness contribution.
		if rel != nil {
			rel.addVariant(filteredVariant)
		}
		if rel2 != nil {
			rel2.addVariant(filteredVariant)
		}

		// Per-variant LROH state update.
		if lroh != nil {
			lroh.addVariant(filteredVariant)
		}

		// Per-variant phased-block state update.
		if phasedBlocks != nil {
			phasedBlocks.addVariant(filteredVariant)
		}

		// Per-variant INFO extraction (--get-INFO).
		if getInfo != nil {
			if err := getInfo.addVariant(filteredVariant); err != nil {
				return fmt.Errorf("writing --get-INFO output: %w", err)
			}
		}

		// Collect variants for format conversions
		if params.Output012 || params.OutputPlink || params.OutputPlinkTped {
			allVariants = append(allVariants, filteredVariant)
		}

		// Write to output if recoding
		if params.Recode {
			var outInfo map[string]string
			switch {
			case len(keepInfoSet) > 0 || len(removeInfoSet) > 0:
				// --keep-INFO / --remove-INFO compose with --recode-INFO-all:
				// start from the full INFO map and project.
				outInfo = filterRecodeInfo(filteredVariant.Info, keepInfoSet, removeInfoSet)
			case params.RecodeInfoAll:
				outInfo = filteredVariant.Info
			default:
				outInfo = make(map[string]string)
			}
			outVariant := &vcf.Variant{
				Chrom:   filteredVariant.Chrom,
				Pos:     filteredVariant.Pos,
				ID:      filteredVariant.ID,
				Ref:     filteredVariant.Ref,
				Alt:     filteredVariant.Alt,
				Qual:    filteredVariant.Qual,
				Filter:  filteredVariant.Filter,
				Info:    outInfo,
				Format:  filteredVariant.Format,
				Samples: filteredVariant.Samples,
			}
			if err := recodeWriter.Write(outVariant); err != nil {
				return fmt.Errorf("writing variant: %w", err)
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

	// Flush LD outputs.
	if err := ldRun.close(); err != nil {
		return fmt.Errorf("closing LD output: %w", err)
	}

	// Flush --diff outputs (also emits file-2-only sites and per-individual
	// reports).
	if err := diffRun.close(); err != nil {
		return fmt.Errorf("closing --diff output: %w", err)
	}

	// Flush BEAGLE outputs.
	if err := beagleGL.close(); err != nil {
		return fmt.Errorf("closing --BEAGLE-GL output: %w", err)
	}
	if err := beaglePL.close(); err != nil {
		return fmt.Errorf("closing --BEAGLE-PL output: %w", err)
	}

	// Inter-chromosomal LD requires all sites in memory; emit pairs now.
	if interLD != nil {
		if err := interLD.flush(); err != nil {
			return fmt.Errorf("flushing interchrom LD: %w", err)
		}
		if err := interLD.close(); err != nil {
			return fmt.Errorf("closing interchrom LD output: %w", err)
		}
	}

	// Relatedness output.
	if rel != nil {
		if err := rel.writeOutput(params.OutPrefix); err != nil {
			return fmt.Errorf("writing --relatedness output: %w", err)
		}
	}
	if rel2 != nil {
		if err := rel2.writeOutput(params.OutPrefix); err != nil {
			return fmt.Errorf("writing --relatedness2 output: %w", err)
		}
	}

	// LROH output.
	if lroh != nil {
		if err := lroh.writeOutput(params.OutPrefix); err != nil {
			return fmt.Errorf("writing --LROH output: %w", err)
		}
	}

	// Phased blocks output.
	if phasedBlocks != nil {
		if err := phasedBlocks.writeOutput(params.OutPrefix); err != nil {
			return fmt.Errorf("writing --phased-blocks output: %w", err)
		}
	}

	// --get-INFO output.
	if err := getInfo.close(); err != nil {
		return fmt.Errorf("closing --get-INFO output: %w", err)
	}

	return nil
}

// splitCSV splits a comma-separated string into trimmed non-empty tokens, in
// order. Used for --get-INFO and similar list-valued flags.
func splitCSV(s string) []string {
	if s == "" {
		return nil
	}
	parts := strings.Split(s, ",")
	out := make([]string, 0, len(parts))
	for _, p := range parts {
		p = strings.TrimSpace(p)
		if p != "" {
			out = append(out, p)
		}
	}
	return out
}

// beagleGLMode and beaglePLMode are tiny helpers so we don't need to expose
// the unexported beagleMode constants outside the package.
func beagleGLMode() beagleMode { return beagleGL }
func beaglePLMode() beagleMode { return beaglePL }

// passFilters checks if a variant passes all filters
func passFilters(v *vcf.Variant, params *Params, includePos, excludePos positionSet, includeSNPs, excludeSNPs map[string]bool, includeBed, excludeBed *bedRegions) bool {
	// SNP ID filters
	if includeSNPs != nil && len(includeSNPs) > 0 {
		if !includeSNPs[v.ID] {
			return false
		}
	}
	if excludeSNPs != nil && len(excludeSNPs) > 0 {
		if excludeSNPs[v.ID] {
			return false
		}
	}

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

	// BED-based include/exclude. Sites must be inside any --bed interval and
	// must not be inside any --exclude-bed interval. Each is independent of
	// the other (no implicit subtraction).
	if includeBed != nil {
		if !includeBed.containsVCFPos(v.Chrom, v.Pos) {
			return false
		}
	}
	if excludeBed != nil {
		if excludeBed.containsVCFPos(v.Chrom, v.Pos) {
			return false
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

	// Genotype filters. Upstream's --max-missing is the MIN fraction of
	// non-missing genotypes (0.0 = allow all, 1.0 = require all
	// non-missing). The Params field is the same semantics: 0 means
	// "feature disabled" (no filter), >0 means apply. We guard against
	// the zero default explicitly so `--max-missing 1.0` (require all
	// non-missing) still applies — the previous `< 1` guard mistakenly
	// dropped that exact case.
	if params.MaxMissing > 0 {
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

// loadSNPIDs loads SNP IDs from a file
func loadSNPIDs(filename string) (map[string]bool, error) {
	f, err := iohelper.OpenReader(filename)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	snpIDs := make(map[string]bool)
	scanner := bufio.NewScanner(f)

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}

		// SNP ID is the first field (or whole line if no whitespace)
		fields := strings.Fields(line)
		if len(fields) > 0 {
			snpIDs[fields[0]] = true
		}
	}

	return snpIDs, scanner.Err()
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

// filterGenotypes applies genotype-level filters (sets genotypes to missing if they fail filters)
func filterGenotypes(v *vcf.Variant, params *Params) *vcf.Variant {
	// If no genotype filters specified, return as-is
	if params.MinDP == 0 && params.MaxDP == 0 && params.MinGQ == 0 {
		return v
	}

	// Create a copy to avoid modifying original
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
		// Check DP (depth) filter
		if params.MinDP > 0 || params.MaxDP > 0 {
			if dpStr, ok := sample.Data["DP"]; ok {
				if dp, err := strconv.Atoi(dpStr); err == nil {
					if (params.MinDP > 0 && dp < params.MinDP) || (params.MaxDP > 0 && dp > params.MaxDP) {
						// Set genotype to missing
						newSample := vcf.Sample{
							Name: sample.Name,
							Data: make(map[string]string),
						}
						for k, v := range sample.Data {
							if k == "GT" {
								newSample.Data[k] = "./."
							} else {
								newSample.Data[k] = v
							}
						}
						filtered.Samples = append(filtered.Samples, newSample)
						continue
					}
				}
			}
		}

		// Check GQ (genotype quality) filter
		if params.MinGQ > 0 {
			if gqStr, ok := sample.Data["GQ"]; ok {
				if gq, err := strconv.Atoi(gqStr); err == nil {
					if gq < params.MinGQ {
						// Set genotype to missing
						newSample := vcf.Sample{
							Name: sample.Name,
							Data: make(map[string]string),
						}
						for k, v := range sample.Data {
							if k == "GT" {
								newSample.Data[k] = "./."
							} else {
								newSample.Data[k] = v
							}
						}
						filtered.Samples = append(filtered.Samples, newSample)
						continue
					}
				}
			}
		}

		// Genotype passes filters, keep as-is
		filtered.Samples = append(filtered.Samples, sample)
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

	if params.Freq2 {
		if err := stats.outputFrequency2(params.OutPrefix); err != nil {
			return err
		}
	}

	if params.Counts2 {
		if err := stats.outputCounts2(params.OutPrefix); err != nil {
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

	if params.TsTvByCount {
		if err := stats.outputTsTvByCount(params.OutPrefix); err != nil {
			return err
		}
	}

	if params.TsTvByQual {
		if err := stats.outputTsTvByQual(params.OutPrefix); err != nil {
			return err
		}
	}

	if params.Depth {
		if err := stats.outputDepth(params.OutPrefix); err != nil {
			return err
		}
	}

	if params.SitePi {
		if err := stats.outputSitePi(params.OutPrefix); err != nil {
			return err
		}
	}

	if params.WindowPi > 0 {
		if err := stats.outputWindowedPi(params.OutPrefix, params.WindowPi, params.WindowPiStep); err != nil {
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

	if params.GenoDepth {
		if err := stats.outputGenoDepth(params.OutPrefix); err != nil {
			return err
		}
	}

	if params.HistIndelLen {
		if err := stats.outputIndelHist(params.OutPrefix); err != nil {
			return err
		}
	}

	if params.TajimaD > 0 {
		if err := stats.outputTajimaD(params.OutPrefix, params.TajimaD); err != nil {
			return err
		}
	}

	// Weir & Cockerham 1984 Fst (per-site, plus optional windowed output).
	if stats.weirFst != nil {
		if err := stats.weirFst.outputWeirFst(params.OutPrefix); err != nil {
			return err
		}
		if params.FstWindowSize > 0 {
			if err := stats.weirFst.outputWindowedWeirFst(params.OutPrefix, params.FstWindowSize, params.FstWindowStep); err != nil {
				return err
			}
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

// checkUnsupported returns an error if the parameters request a feature that
// this Go port does not implement yet. Previously these options were accepted
// and silently ignored, which produced no output and looked like success.
func checkUnsupported(params *Params) error {
	_ = params
	var missing []string
	if len(missing) > 0 {
		return fmt.Errorf("not implemented in this Go port yet: %s (see tools/vcftools/ROADMAP.md)", strings.Join(missing, ", "))
	}
	return nil
}
