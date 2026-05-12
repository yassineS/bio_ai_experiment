// vcftools - Utilities for working with VCF (Variant Call Format) files
package main

import (
	"flag"
	"fmt"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
	"github.com/yassineS/bio_ai_experiment/tools/vcftools/pkg/vcftools"
)

const usage = `vcftools - Utilities for VCF (Variant Call Format) files

Usage:
  vcftools [input options] [filtering options] [output options]

Description:
  A suite of functions for working with genetic variation data in VCF format.
  Tools are used to summarize data, filter data, and output in various formats.

Input Options:
  --vcf FILE            Input VCF file
  --gzvcf FILE          Input gzipped VCF file
  --stdin               Read VCF from stdin

Output Options:
  --out PREFIX          Output filename prefix (default: out)
  --stdout              Write output to stdout
  --recode              Output a new VCF file after filtering
  --recode-INFO-all     Include all INFO fields in recoded VCF

Position Filtering:
  --chr STRING          Include only this chromosome
  --not-chr STRING      Exclude this chromosome
  --from-bp INT         Include positions >= this value
  --to-bp INT           Include positions <= this value
  --positions FILE      Include positions listed in file
  --exclude-positions FILE  Exclude positions listed in file

Variant Type Filtering:
  --keep-only-indels    Keep only indel sites
  --remove-indels       Remove all indel sites
  --min-alleles INT     Minimum number of alleles (default: 2)
  --max-alleles INT     Maximum number of alleles

Quality Filtering:
  --minQ FLOAT          Minimum quality score
  --remove-filtered-all Remove sites with FILTER != PASS

Allele Frequency Filtering:
  --maf FLOAT           Minimum minor allele frequency
  --max-maf FLOAT       Maximum minor allele frequency
  --mac INT             Minimum minor allele count
  --max-mac INT         Maximum minor allele count

Genotype Filtering:
  --max-missing FLOAT   Maximum proportion of missing data (0-1)
  --min-meanDP FLOAT    Minimum mean depth across samples
  --max-meanDP FLOAT    Maximum mean depth across samples

Statistics Output:
  --freq                Output allele frequency
  --counts              Output allele counts
  --depth               Output mean depth per site
  --site-depth          Output depth for each site
  --site-mean-depth     Output mean depth per site
  --site-quality        Output quality scores per site
  --missing-indv        Output individual missingness
  --missing-site        Output site missingness
  --hardy               Test for Hardy-Weinberg equilibrium
  --TsTv-summary        Output Ts/Tv ratio summary
  --TsTv INT            Output Ts/Tv in bins of size INT
  --site-pi             Output nucleotide diversity per site

Sample Filtering:
  --indv STRING         Include only this individual (can use multiple times)
  --remove-indv STRING  Remove this individual (can use multiple times)
  --keep FILE           Keep only individuals listed in file
  --remove FILE         Remove individuals listed in file

Help:
  -h, --help            Show this help message

Examples:
  # Get allele frequency for chromosome 1
  vcftools --vcf input.vcf --chr 1 --freq --out chr1

  # Filter by quality and allele frequency
  vcftools --vcf input.vcf --minQ 30 --maf 0.05 --recode --out filtered

  # Remove indels and output new VCF
  vcftools --vcf input.vcf --remove-indels --recode --recode-INFO-all --out snps_only

  # Calculate Hardy-Weinberg equilibrium
  vcftools --vcf input.vcf --hardy --out hwe_test

  # Get statistics with compressed input
  vcftools --gzvcf input.vcf.gz --freq --depth --out stats

Notes:
  - Supports VCF v4.0, v4.1, and v4.2
  - Automatically handles gzipped files with --gzvcf
  - Use --stdin to read from standard input
  - Multiple filters can be combined
`

func main() {
	// Input options
	vcfFile := flag.String("vcf", "", "Input VCF file")
	gzvcfFile := flag.String("gzvcf", "", "Input gzipped VCF file")
	useStdin := flag.Bool("stdin", false, "Read from stdin")

	// Output options
	outPrefix := flag.String("out", "out", "Output filename prefix")
	useStdout := flag.Bool("stdout", false, "Write to stdout")
	recode := flag.Bool("recode", false, "Output a new VCF file")
	recodeInfoAll := flag.Bool("recode-INFO-all", false, "Include all INFO fields in recode")

	// Position filtering
	chr := flag.String("chr", "", "Include only this chromosome")
	notChr := flag.String("not-chr", "", "Exclude this chromosome")
	fromBp := flag.Int("from-bp", 0, "Include positions >= this value")
	toBp := flag.Int("to-bp", 0, "Include positions <= this value")
	positionsFile := flag.String("positions", "", "Include positions from file")
	excludePositionsFile := flag.String("exclude-positions", "", "Exclude positions from file")

	// SNP ID filtering
	snp := flag.String("snp", "", "Include only this SNP ID")
	snps := flag.String("snps", "", "Include SNP IDs from file")
	excludeSNP := flag.String("exclude", "", "Exclude this SNP ID")
	excludeSNPs := flag.String("exclude-snps", "", "Exclude SNP IDs from file")
	thin := flag.Int("thin", 0, "Thin sites by keeping every Nth site")

	// Variant type filtering
	keepOnlyIndels := flag.Bool("keep-only-indels", false, "Keep only indels")
	removeIndels := flag.Bool("remove-indels", false, "Remove indels")
	minAlleles := flag.Int("min-alleles", 2, "Minimum number of alleles")
	maxAlleles := flag.Int("max-alleles", 0, "Maximum number of alleles")

	// Quality filtering
	minQ := flag.Float64("minQ", 0, "Minimum quality score")
	removeFilteredAll := flag.Bool("remove-filtered-all", false, "Remove sites with FILTER != PASS")

	// Allele frequency filtering
	maf := flag.Float64("maf", 0, "Minimum minor allele frequency")
	maxMaf := flag.Float64("max-maf", 0, "Maximum minor allele frequency")
	mac := flag.Int("mac", 0, "Minimum minor allele count")
	maxMac := flag.Int("max-mac", 0, "Maximum minor allele count")

	// Genotype filtering
	maxMissing := flag.Float64("max-missing", 1, "Maximum proportion of missing data")
	minMeanDP := flag.Float64("min-meanDP", 0, "Minimum mean depth")
	maxMeanDP := flag.Float64("max-meanDP", 0, "Maximum mean depth")
	minDP := flag.Int("minDP", 0, "Minimum depth per genotype")
	maxDP := flag.Int("maxDP", 0, "Maximum depth per genotype")
	minGQ := flag.Int("minGQ", 0, "Minimum genotype quality")

	// Statistics output
	freq := flag.Bool("freq", false, "Output allele frequency")
	counts := flag.Bool("counts", false, "Output allele counts")
	freq2 := flag.Bool("freq2", false, "Alternative frequency output format")
	counts2 := flag.Bool("counts2", false, "Alternative counts output format")
	depth := flag.Bool("depth", false, "Output mean depth per site")
	siteDepth := flag.Bool("site-depth", false, "Output depth for each site")
	siteMeanDepth := flag.Bool("site-mean-depth", false, "Output mean depth per site")
	siteQuality := flag.Bool("site-quality", false, "Output quality per site")
	missingIndv := flag.Bool("missing-indv", false, "Output individual missingness")
	missingSite := flag.Bool("missing-site", false, "Output site missingness")
	hardy := flag.Bool("hardy", false, "Hardy-Weinberg equilibrium test")
	tsTvSummary := flag.Bool("TsTv-summary", false, "Ts/Tv ratio summary")
	tsTvBinSize := flag.Int("TsTv", 0, "Ts/Tv in bins of this size")
	tsTvByCount := flag.Bool("TsTv-by-count", false, "Ts/Tv by allele count")
	tsTvByQual := flag.Bool("TsTv-by-qual", false, "Ts/Tv by quality score")
	sitePi := flag.Bool("site-pi", false, "Nucleotide diversity per site")
	het := flag.Bool("het", false, "Heterozygosity statistics")
	singletons := flag.Bool("singletons", false, "Singleton site analysis")
	histIndelLen := flag.Bool("hist-indel-len", false, "Indel length histogram")
	genoDepth := flag.Bool("geno-depth", false, "Genotype depth distribution")

	// Phase 2: Population genetics statistics
	windowPi := flag.Int("window-pi", 0, "Nucleotide diversity in windows of this size")
	windowPiStep := flag.Int("window-pi-step", 0, "Step size for pi windows")
	tajimaD := flag.Int("TajimaD", 0, "Tajima's D in bins of this size")
	snpDensity := flag.Int("SNPdensity", 0, "SNP density in bins of this size")
	var weirFstPop []string
	flag.Func("weir-fst-pop", "Population file for Fst calculation (can use multiple times)", func(s string) error {
		weirFstPop = append(weirFstPop, s)
		return nil
	})
	fstWindowSize := flag.Int("fst-window-size", 0, "Window size for Fst")
	fstWindowStep := flag.Int("fst-window-step", 0, "Step size for Fst windows")
	filterSummary := flag.Bool("FILTER-summary", false, "FILTER tag summary")

	// Phase 4: Format conversions
	output012 := flag.Bool("012", false, "Output genotypes as 0/1/2 matrix")
	outputPlink := flag.Bool("plink", false, "Output PLINK PED/MAP format")
	outputPlinkTped := flag.Bool("plink-tped", false, "Output PLINK TPED/TFAM format")
	chromMap := flag.String("chrom-map", "", "Chromosome name to integer mapping file")

	// Sample filtering
	var indvList, removeIndvList []string
	flag.Func("indv", "Include individual (can use multiple times)", func(s string) error {
		indvList = append(indvList, s)
		return nil
	})
	flag.Func("remove-indv", "Remove individual (can use multiple times)", func(s string) error {
		removeIndvList = append(removeIndvList, s)
		return nil
	})
	keepFile := flag.String("keep", "", "Keep individuals from file")
	removeFile := flag.String("remove", "", "Remove individuals from file")

	// Help
	help := flag.Bool("h", false, "Show help message")
	flag.BoolVar(help, "help", false, "Show help message")

	flag.Parse()

	if *help {
		fmt.Fprint(os.Stderr, usage)
		os.Exit(0)
	}

	// Validate input
	inputCount := 0
	if *vcfFile != "" {
		inputCount++
	}
	if *gzvcfFile != "" {
		inputCount++
	}
	if *useStdin {
		inputCount++
	}

	if inputCount == 0 {
		fmt.Fprintln(os.Stderr, "Error: Must specify one of --vcf, --gzvcf, or --stdin")
		fmt.Fprintln(os.Stderr, "Use --help for usage information")
		os.Exit(1)
	}
	if inputCount > 1 {
		fmt.Fprintln(os.Stderr, "Error: Can only specify one input method")
		os.Exit(1)
	}

	// Determine input file
	var inputFile string
	if *vcfFile != "" {
		inputFile = *vcfFile
	} else if *gzvcfFile != "" {
		inputFile = *gzvcfFile
	} else {
		inputFile = "" // stdin
	}

	// Open input
	inputReader, err := iohelper.OpenReader(inputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error opening input: %v\n", err)
		os.Exit(1)
	}
	defer inputReader.Close()

	// Build params
	params := &vcftools.Params{
		OutPrefix:            *outPrefix,
		UseStdout:            *useStdout,
		Recode:               *recode,
		RecodeInfoAll:        *recodeInfoAll,
		Chr:                  *chr,
		NotChr:               *notChr,
		FromBp:               *fromBp,
		ToBp:                 *toBp,
		PositionsFile:        *positionsFile,
		ExcludePositionsFile: *excludePositionsFile,
		SNP:                  *snp,
		SNPs:                 *snps,
		ExcludeSNP:           *excludeSNP,
		ExcludeSNPs:          *excludeSNPs,
		Thin:                 *thin,
		KeepOnlyIndels:       *keepOnlyIndels,
		RemoveIndels:         *removeIndels,
		MinAlleles:           *minAlleles,
		MaxAlleles:           *maxAlleles,
		MinQ:                 *minQ,
		RemoveFilteredAll:    *removeFilteredAll,
		Maf:                  *maf,
		MaxMaf:               *maxMaf,
		Mac:                  *mac,
		MaxMac:               *maxMac,
		MaxMissing:           *maxMissing,
		MinMeanDP:            *minMeanDP,
		MaxMeanDP:            *maxMeanDP,
		MinDP:                *minDP,
		MaxDP:                *maxDP,
		MinGQ:                *minGQ,
		Freq:                 *freq,
		Counts:               *counts,
		Freq2:                *freq2,
		Counts2:              *counts2,
		Depth:                *depth,
		SiteDepth:            *siteDepth,
		SiteMeanDepth:        *siteMeanDepth,
		SiteQuality:          *siteQuality,
		MissingIndv:          *missingIndv,
		MissingSite:          *missingSite,
		Hardy:                *hardy,
		TsTvSummary:          *tsTvSummary,
		TsTvBinSize:          *tsTvBinSize,
		TsTvByCount:          *tsTvByCount,
		TsTvByQual:           *tsTvByQual,
		SitePi:               *sitePi,
		Het:                  *het,
		Singletons:           *singletons,
		HistIndelLen:         *histIndelLen,
		GenoDepth:            *genoDepth,
		WindowPi:             *windowPi,
		WindowPiStep:         *windowPiStep,
		TajimaD:              *tajimaD,
		SNPDensity:           *snpDensity,
		WeirFstPop:           weirFstPop,
		FstWindowSize:        *fstWindowSize,
		FstWindowStep:        *fstWindowStep,
		FilterSummary:        *filterSummary,
		Output012:            *output012,
		OutputPlink:          *outputPlink,
		OutputPlinkTped:      *outputPlinkTped,
		ChromMap:             *chromMap,
		IndvList:             indvList,
		RemoveIndvList:       removeIndvList,
		KeepFile:             *keepFile,
		RemoveFile:           *removeFile,
	}

	// Run vcftools
	if err := vcftools.Run(inputReader, params); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}
}
