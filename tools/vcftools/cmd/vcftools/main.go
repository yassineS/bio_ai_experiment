// vcftools - Utilities for working with VCF (Variant Call Format) files
package main

import (
	"flag"
	"fmt"
	"os"
	"strconv"
	"strings"

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
  --bed FILE            Keep only sites inside any interval in BED FILE
  --exclude-bed FILE    Remove sites inside any interval in BED FILE

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
  --non-ref-af FLOAT    Minimum non-reference allele frequency (per ALT)
  --non-ref-ac INT      Minimum non-reference allele count (per ALT)
  --max-non-ref-af FLOAT
                        Maximum non-reference allele frequency (per ALT)
  --max-non-ref-ac INT  Maximum non-reference allele count (per ALT)
  --non-ref-af-any FLOAT
                        Minimum non-reference allele frequency (any ALT
                        passes). NOTE: upstream parity quirk — see docs.
  --non-ref-ac-any INT  Minimum non-reference allele count (any ALT passes)
  --max-non-ref-af-any FLOAT
                        Maximum non-reference allele frequency (any ALT
                        passes). NOTE: upstream parity quirk — see docs.
  --max-non-ref-ac-any INT
                        Maximum non-reference allele count (any ALT passes)

Genotype Filtering:
  --max-missing FLOAT   Maximum proportion of missing data (0-1)
  --max-missing-count INT
                        Maximum number of missing chromosomes (haploid
                        alleles, NOT samples) tolerated per site. "0"
                        means "drop any site with any missing call".
  --hwe FLOAT           Minimum exact-test (Wigginton 2005) HWE p-value
                        per biallelic site. Setting --hwe also forces
                        --max-alleles 2 (matches upstream).
  --min-meanDP FLOAT    Minimum mean depth across samples
  --max-meanDP FLOAT    Maximum mean depth across samples
  --phased              Keep only sites where every kept-individual GT is
                        phased (separator '|' or haploid).

Statistics Output:
  --freq                Output allele frequency
  --counts              Output allele counts
  --freq2               Alternative allele frequency format
  --counts2             Alternative allele counts format
  --depth               Output mean read depth per individual (.idepth)
  --site-depth          Output summed depth for each site (.ldepth)
  --site-mean-depth     Output mean depth per site (.ldepth.mean)
  --site-quality        Output quality scores per site
  --missing-indv        Output individual missingness
  --missing-site        Output site missingness
  --hardy               Test for Hardy-Weinberg equilibrium
  --het                 Output per-individual heterozygosity / F
  --singletons          Output singleton/private sites
  --TsTv-summary        Output Ts/Tv ratio summary
  --TsTv INT            Output Ts/Tv in bins of size INT
  --TsTv-by-count       Output Ts/Tv grouped by alternate-allele count
  --TsTv-by-qual        Output Ts/Tv bucketed by quality-score thresholds (.TsTv.qual)
  --site-pi             Output nucleotide diversity per site (.sites.pi)
  --hist-indel-len      Output a histogram of indel lengths (.indel.hist)
  --geno-depth          Output a per-genotype read-depth matrix (.gdepth)
  --FILTER-summary      Output a summary of FILTER values
  --SNPdensity INT      Output SNP density in bins of size INT

Population Genetics:
  --window-pi INT       Nucleotide diversity summed over windows of size INT
  --window-pi-step INT  Step size for --window-pi windows (default: window size)
  --TajimaD INT         Tajima's D in non-overlapping windows of size INT
  --weir-fst-pop FILE   Population file (one sample per line) for Weir & Cockerham
                        1984 Fst; use the flag two or more times, once per pop.
                        Writes <prefix>.weir.fst.
  --fst-window-size INT Window size for Fst calculation; also writes
                        <prefix>.windowed.weir.fst
  --fst-window-step INT Step size for sliding Fst windows (default: window size)

Format Conversion:
  --012                 Output genotypes as a 0/1/2 matrix
  --plink               Output PLINK PED/MAP files
  --plink-tped          Output PLINK TPED/TFAM files
  --chrom-map FILE      Chromosome-name-to-integer map for PLINK output
  --BEAGLE-GL           Output <prefix>.BEAGLE.GL (log10 GL from PL)
  --BEAGLE-PL           Output <prefix>.BEAGLE.PL (raw PL triplets)
  --ldhat               Output phased LDhat format (<prefix>.ldhat.sites and
                        <prefix>.ldhat.locs). Requires --chr; implies --phased.
  --ldhat-geno          Output unphased LDhat format (same file names as
                        --ldhat). Requires --chr.
  --ldhelmet            Output LDhelmet format (<prefix>.ldhelmet.snps and
                        <prefix>.ldhelmet.pos). Requires --chr; implies
                        --phased and --remove-indels.
  --IMPUTE              Output IMPUTE reference-panel bundle
                        (<prefix>.impute.legend, <prefix>.impute.hap,
                        <prefix>.impute.hap.indv). Biallelic, phased SNPs
                        only; sites with missing data are dropped.

VCF Comparison (--diff family):
  --diff FILE                Compare against a second VCF file
  --diff-site                Emit <prefix>.diff.sites_in_files
  --diff-indv                Emit <prefix>.diff.indv_in_files
  --diff-site-discordance    Emit <prefix>.diff.sites
  --diff-indv-discordance    Emit <prefix>.diff.indv
  --diff-indv-map FILE       Two-column file renaming file-2 sample IDs
                             before matching against file-1
  --diff-discordance-matrix  Emit <prefix>.diff.discordance_matrix (4x4
                             genotype-by-genotype counts for biallelic loci)
  --diff-switch-error        Emit <prefix>.diff.switch (per-event log) and
                             <prefix>.diff.indv.switch (per-individual
                             phase-switch error rate) vs --diff file

Mendelian Inconsistency:
  --mendel FILE              PED file (four columns: family child father
                             mother) used to detect Mendelian errors in
                             trios. Emits <prefix>.mendel.

Linkage Disequilibrium:
  --geno-r2             Genotype-based LD r^2 within a window (.geno.ld)
  --hap-r2              Haplotype-based LD r^2 for phased data (.hap.ld)
  --geno-r2-positions FILE
                        Restrict --geno-r2 to pairs touching a position in FILE
  --hap-r2-positions FILE
                        Restrict --hap-r2 to pairs touching a position in FILE
  --ld-window INT       Maximum number of SNPs between LD pairs (default: unbounded)
  --ld-window-bp INT    Maximum bp distance between LD pairs (default: unbounded)
  --ld-window-min INT   Minimum number of SNPs between LD pairs
  --ld-window-bp-min INT  Minimum bp distance between LD pairs
  --min-r2 FLOAT        Only emit LD pairs with r^2 >= this threshold
  --interchrom-geno-r2  Inter-chromosomal genotype LD r^2 (.interchrom.geno.ld)
  --interchrom-hap-r2   Inter-chromosomal haplotype LD r^2 (.interchrom.hap.ld)
  --geno-chisq          Per-pair chi-square test of genotype association (.geno.chisq)

Advanced Analysis:
  --relatedness         Yang 2010 unadjusted A_jk relatedness (.relatedness)
  --relatedness2        KING-robust kinship coefficient (.relatedness2)
  --LROH                Runs of homozygosity per individual (.LROH)
  --LROH-min-variants INT
                        Minimum consecutive homozygous variants for a run
                        (default: 10)
  --phased-blocks       Per-individual contiguous phased-haplotype block
                        boundaries (.blocks)

FILTER / INFO Selection:
  --remove-filtered NAME[,NAME...]   Drop sites listing any of these FILTERs
  --keep-filtered NAME[,NAME...]     Keep only sites listing any of these
                                     FILTERs
  --keep-INFO TAG       Keep only this INFO tag in --recode output (repeatable)
  --remove-INFO TAG     Strip this INFO tag from --recode output (repeatable)
  --get-INFO TAG        Extract this INFO tag to <prefix>.INFO (repeatable)

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
	nonRefAF := flag.Float64("non-ref-af", 0, "Minimum non-reference allele frequency (per ALT)")
	nonRefAC := flag.Int("non-ref-ac", 0, "Minimum non-reference allele count (per ALT)")
	maxNonRefAF := flag.Float64("max-non-ref-af", 0, "Maximum non-reference allele frequency (per ALT)")
	maxNonRefAC := flag.Int("max-non-ref-ac", 0, "Maximum non-reference allele count (per ALT)")
	nonRefAFAny := flag.Float64("non-ref-af-any", 0, "Minimum non-reference allele frequency (any ALT passes); upstream parity quirk: no-op alone")
	nonRefACAny := flag.Int("non-ref-ac-any", 0, "Minimum non-reference allele count (any ALT passes)")
	maxNonRefAFAny := flag.Float64("max-non-ref-af-any", 0, "Maximum non-reference allele frequency (any ALT passes); upstream parity quirk: no-op alone")
	maxNonRefACAny := flag.Int("max-non-ref-ac-any", 0, "Maximum non-reference allele count (any ALT passes)")

	// Genotype filtering
	maxMissing := flag.Float64("max-missing", 1, "Maximum proportion of missing data")
	// --max-missing-count: use flag.Func so we can record whether the
	// flag was supplied at all (vs defaulted), since "0" is a meaningful
	// user-supplied value (drop any site with any missing call).
	var maxMissingCount int
	var maxMissingCountSet bool
	flag.Func("max-missing-count", "Maximum number of missing chromosomes (haploid alleles, NOT samples) tolerated per site", func(s string) error {
		v, err := strconv.Atoi(strings.TrimSpace(s))
		if err != nil {
			return fmt.Errorf("--max-missing-count: %w", err)
		}
		if v < 0 {
			return fmt.Errorf("--max-missing-count must be >= 0")
		}
		maxMissingCount = v
		maxMissingCountSet = true
		return nil
	})
	minMeanDP := flag.Float64("min-meanDP", 0, "Minimum mean depth")
	maxMeanDP := flag.Float64("max-meanDP", 0, "Maximum mean depth")
	minDP := flag.Int("minDP", 0, "Minimum depth per genotype")
	maxDP := flag.Int("maxDP", 0, "Maximum depth per genotype")
	minGQ := flag.Int("minGQ", 0, "Minimum genotype quality")
	// --hwe FLOAT: minimum exact-test HWE p-value per site (biallelic;
	// upstream also forces max_alleles=2 when this flag is set —
	// parameters.cpp:254). We apply max_alleles=2 below at param-build
	// time to mirror upstream's behaviour exactly.
	hwePvalue := flag.Float64("hwe", 0, "Minimum exact-test HWE p-value per biallelic site (Wigginton 2005)")

	// Statistics output
	freq := flag.Bool("freq", false, "Output allele frequency")
	counts := flag.Bool("counts", false, "Output allele counts")
	freq2 := flag.Bool("freq2", false, "Alternative frequency output format")
	counts2 := flag.Bool("counts2", false, "Alternative counts output format")
	depth := flag.Bool("depth", false, "Output mean read depth per individual (.idepth)")
	siteDepth := flag.Bool("site-depth", false, "Output summed depth for each site (.ldepth)")
	siteMeanDepth := flag.Bool("site-mean-depth", false, "Output mean depth per site (.ldepth.mean)")
	siteQuality := flag.Bool("site-quality", false, "Output quality per site")
	missingIndv := flag.Bool("missing-indv", false, "Output individual missingness")
	missingSite := flag.Bool("missing-site", false, "Output site missingness")
	hardy := flag.Bool("hardy", false, "Hardy-Weinberg equilibrium test")
	tsTvSummary := flag.Bool("TsTv-summary", false, "Ts/Tv ratio summary")
	tsTvBinSize := flag.Int("TsTv", 0, "Ts/Tv in bins of this size")
	tsTvByCount := flag.Bool("TsTv-by-count", false, "Ts/Tv grouped by alternate-allele count")
	tsTvByQual := flag.Bool("TsTv-by-qual", false, "Ts/Tv ratios bucketed by quality-score thresholds (.TsTv.qual)")
	sitePi := flag.Bool("site-pi", false, "Nucleotide diversity per site (.sites.pi)")
	het := flag.Bool("het", false, "Heterozygosity statistics")
	singletons := flag.Bool("singletons", false, "Singleton site analysis")
	histIndelLen := flag.Bool("hist-indel-len", false, "Histogram of indel lengths (.indel.hist)")
	genoDepth := flag.Bool("geno-depth", false, "Per-genotype read-depth matrix (.gdepth)")

	// Population genetics statistics
	windowPi := flag.Int("window-pi", 0, "Nucleotide diversity summed over windows of this size")
	windowPiStep := flag.Int("window-pi-step", 0, "Step size for --window-pi windows")
	tajimaD := flag.Int("TajimaD", 0, "Tajima's D in non-overlapping windows of this size")
	snpDensity := flag.Int("SNPdensity", 0, "SNP density in bins of this size")
	var weirFstPop []string
	flag.Func("weir-fst-pop", "Population file (one sample per line) for Weir & Cockerham 1984 Fst, repeatable", func(s string) error {
		weirFstPop = append(weirFstPop, s)
		return nil
	})
	fstWindowSize := flag.Int("fst-window-size", 0, "Window size for Fst calculation")
	fstWindowStep := flag.Int("fst-window-step", 0, "Step size for sliding Fst windows")
	filterSummary := flag.Bool("FILTER-summary", false, "FILTER tag summary")

	// Linkage disequilibrium analysis
	genoR2 := flag.Bool("geno-r2", false, "Genotype-based LD r^2 within a window (.geno.ld)")
	hapR2 := flag.Bool("hap-r2", false, "Haplotype-based LD r^2 for phased data (.hap.ld)")
	genoR2Positions := flag.String("geno-r2-positions", "", "Restrict --geno-r2 to pairs touching a position in this file (chrom pos)")
	hapR2Positions := flag.String("hap-r2-positions", "", "Restrict --hap-r2 to pairs touching a position in this file (chrom pos)")
	ldWindow := flag.Int("ld-window", 0, "Maximum number of SNPs between LD pairs (default: unbounded)")
	ldWindowBp := flag.Int("ld-window-bp", 0, "Maximum bp distance between LD pairs (default: unbounded)")
	ldWindowMin := flag.Int("ld-window-min", 0, "Minimum number of SNPs between LD pairs")
	ldWindowBpMin := flag.Int("ld-window-bp-min", 0, "Minimum bp distance between LD pairs")
	minR2 := flag.Float64("min-r2", 0, "Only emit LD pairs with r^2 >= this threshold")

	// Phase 4: Format conversions
	output012 := flag.Bool("012", false, "Output genotypes as 0/1/2 matrix")
	outputPlink := flag.Bool("plink", false, "Output PLINK PED/MAP format")
	outputPlinkTped := flag.Bool("plink-tped", false, "Output PLINK TPED/TFAM format")
	chromMap := flag.String("chrom-map", "", "Chromosome name to integer mapping file")

	// LDhat output. --ldhat is phased; --ldhat-geno is unphased.
	ldhat := flag.Bool("ldhat", false, "Output phased LDhat format (.ldhat.sites/.ldhat.locs); requires --chr")
	ldhatGeno := flag.Bool("ldhat-geno", false, "Output unphased LDhat format (.ldhat.sites/.ldhat.locs); requires --chr")

	// LDhelmet output. Implies --phased and --remove-indels (upstream
	// parameters.cpp:275). Requires --chr like the other LDhat-family
	// flags.
	ldhelmet := flag.Bool("ldhelmet", false, "Output LDhelmet format (.ldhelmet.snps/.ldhelmet.pos); requires --chr; implies --phased + --remove-indels")

	// IMPUTE output (case-sensitive flag name to match upstream). Implies
	// --phased, biallelic-only, and rejects any site with a missing GT.
	impute := flag.Bool("IMPUTE", false, "Output IMPUTE reference-panel format (.impute.legend/.impute.hap/.impute.hap.indv); phased biallelic SNPs with no missing data only")

	// --pca / --pca-no-norm / --pca-snp-loadings INT: registered for CLI
	// parity (so misuse is reported clearly) but not yet implemented.
	// Run() will reject these via checkUnsupported. See
	// docs/PARITY_ROADMAP.md#vcftools (wave 8 deferral note) for the
	// scope of what's required to land them.
	pca := flag.Bool("pca", false, "Principal component analysis (NOT IMPLEMENTED — see docs/PARITY_ROADMAP.md#vcftools)")
	pcaNoNorm := flag.Bool("pca-no-norm", false, "PCA without normalisation (NOT IMPLEMENTED — see docs/PARITY_ROADMAP.md#vcftools)")
	pcaSNPLoadings := flag.Int("pca-snp-loadings", 0, "Number of top PCs for SNP loadings (NOT IMPLEMENTED — see docs/PARITY_ROADMAP.md#vcftools)")

	// --phased: keep only sites where every kept-individual GT is phased.
	phased := flag.Bool("phased", false, "Keep only sites where every kept-individual GT is phased (separator '|' or haploid)")

	// BED-based filtering
	bedFile := flag.String("bed", "", "Keep only sites whose POS lies inside any interval in this BED file")
	excludeBedFile := flag.String("exclude-bed", "", "Remove sites whose POS lies inside any interval in this BED file")

	// VCF comparison (--diff family)
	diffFile := flag.String("diff", "", "Second VCF file to compare against")
	diffSite := flag.Bool("diff-site", false, "Emit <prefix>.diff.sites_in_files")
	diffIndv := flag.Bool("diff-indv", false, "Emit <prefix>.diff.indv_in_files")
	diffSiteDiscord := flag.Bool("diff-site-discordance", false, "Emit <prefix>.diff.sites with site-by-site discordance")
	diffIndvDiscord := flag.Bool("diff-indv-discordance", false, "Emit <prefix>.diff.indv with per-individual discordance")
	diffIndvMap := flag.String("diff-indv-map", "", "Two-column file mapping file-2 sample IDs to their file-1 equivalents")
	diffDiscMatrix := flag.Bool("diff-discordance-matrix", false, "Emit <prefix>.diff.discordance_matrix (4x4 genotype-by-genotype counts)")
	diffSwitchError := flag.Bool("diff-switch-error", false, "Emit <prefix>.diff.switch and <prefix>.diff.indv.switch (phase-switch error vs --diff file)")

	// --mendel takes a PED file path; emits <prefix>.mendel (Mendelian
	// inconsistencies across trios).
	mendelPed := flag.String("mendel", "", "PED file for Mendelian inconsistency check (emits <prefix>.mendel)")

	// BEAGLE genotype-likelihood output
	beagleGL := flag.Bool("BEAGLE-GL", false, "Emit <prefix>.BEAGLE.GL (log10 GL triplets from PL)")
	beaglePL := flag.Bool("BEAGLE-PL", false, "Emit <prefix>.BEAGLE.PL (raw PL triplets)")

	// Inter-chromosomal LD + chi-square LD
	interchromGenoR2 := flag.Bool("interchrom-geno-r2", false, "Inter-chromosomal genotype LD r^2 (.interchrom.geno.ld)")
	interchromHapR2 := flag.Bool("interchrom-hap-r2", false, "Inter-chromosomal haplotype LD r^2 (.interchrom.hap.ld)")
	genoChiSq := flag.Bool("geno-chisq", false, "Per-pair chi-square test of genotype association (.geno.chisq)")

	// Relatedness, LROH
	relatedness := flag.Bool("relatedness", false, "Yang 2010 unadjusted A_jk relatedness (.relatedness)")
	relatedness2 := flag.Bool("relatedness2", false, "KING-robust kinship coefficient (.relatedness2)")
	lroh := flag.Bool("LROH", false, "Runs of homozygosity per individual (.LROH)")
	lrohMin := flag.Int("LROH-min-variants", 0, "Minimum consecutive homozygous variants for an LROH run (default 10)")
	phasedBlocks := flag.Bool("phased-blocks", false, "Per-individual contiguous phased-haplotype block boundaries (.blocks)")

	// FILTER-name include/exclude
	removeFiltered := flag.String("remove-filtered", "", "Comma-separated FILTER names to drop")
	keepFiltered := flag.String("keep-filtered", "", "Comma-separated FILTER names to keep")

	// INFO tag handling. --keep-INFO / --remove-INFO are repeatable in
	// upstream; --get-INFO is upstream-repeatable too. We accept either
	// repeated single-tag invocations or one comma-separated value, joined
	// with commas in the same order seen on the command line.
	var keepINFOParts, removeINFOParts, getINFOParts []string
	flag.Func("keep-INFO", "INFO tag to keep in --recode output (repeatable)", func(s string) error {
		keepINFOParts = append(keepINFOParts, s)
		return nil
	})
	flag.Func("remove-INFO", "INFO tag to strip from --recode output (repeatable)", func(s string) error {
		removeINFOParts = append(removeINFOParts, s)
		return nil
	})
	flag.Func("get-INFO", "INFO tag to extract to <prefix>.INFO (repeatable)", func(s string) error {
		getINFOParts = append(getINFOParts, s)
		return nil
	})

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
		OutPrefix:             *outPrefix,
		UseStdout:             *useStdout,
		Recode:                *recode,
		RecodeInfoAll:         *recodeInfoAll,
		Chr:                   *chr,
		NotChr:                *notChr,
		FromBp:                *fromBp,
		ToBp:                  *toBp,
		PositionsFile:         *positionsFile,
		ExcludePositionsFile:  *excludePositionsFile,
		SNP:                   *snp,
		SNPs:                  *snps,
		ExcludeSNP:            *excludeSNP,
		ExcludeSNPs:           *excludeSNPs,
		Thin:                  *thin,
		KeepOnlyIndels:        *keepOnlyIndels,
		RemoveIndels:          *removeIndels,
		MinAlleles:            *minAlleles,
		MaxAlleles:            *maxAlleles,
		MinQ:                  *minQ,
		RemoveFilteredAll:     *removeFilteredAll,
		Maf:                   *maf,
		MaxMaf:                *maxMaf,
		Mac:                   *mac,
		MaxMac:                *maxMac,
		MinNonRefAF:           *nonRefAF,
		MinNonRefAC:           *nonRefAC,
		MaxNonRefAF:           *maxNonRefAF,
		MaxNonRefAC:           *maxNonRefAC,
		MinNonRefAFAny:        *nonRefAFAny,
		MinNonRefACAny:        *nonRefACAny,
		MaxNonRefAFAny:        *maxNonRefAFAny,
		MaxNonRefACAny:        *maxNonRefACAny,
		MaxMissing:            *maxMissing,
		MaxMissingCount:       maxMissingCount,
		MaxMissingCountSet:    maxMissingCountSet,
		MinHWEPvalue:          *hwePvalue,
		MinMeanDP:             *minMeanDP,
		MaxMeanDP:             *maxMeanDP,
		MinDP:                 *minDP,
		MaxDP:                 *maxDP,
		MinGQ:                 *minGQ,
		Freq:                  *freq,
		Counts:                *counts,
		Freq2:                 *freq2,
		Counts2:               *counts2,
		Depth:                 *depth,
		SiteDepth:             *siteDepth,
		SiteMeanDepth:         *siteMeanDepth,
		SiteQuality:           *siteQuality,
		MissingIndv:           *missingIndv,
		MissingSite:           *missingSite,
		Hardy:                 *hardy,
		TsTvSummary:           *tsTvSummary,
		TsTvBinSize:           *tsTvBinSize,
		TsTvByCount:           *tsTvByCount,
		TsTvByQual:            *tsTvByQual,
		SitePi:                *sitePi,
		Het:                   *het,
		Singletons:            *singletons,
		HistIndelLen:          *histIndelLen,
		GenoDepth:             *genoDepth,
		WindowPi:              *windowPi,
		WindowPiStep:          *windowPiStep,
		TajimaD:               *tajimaD,
		SNPDensity:            *snpDensity,
		WeirFstPop:            weirFstPop,
		FstWindowSize:         *fstWindowSize,
		FstWindowStep:         *fstWindowStep,
		FilterSummary:         *filterSummary,
		Output012:             *output012,
		OutputPlink:           *outputPlink,
		OutputPlinkTped:       *outputPlinkTped,
		ChromMap:              *chromMap,
		IndvList:              indvList,
		RemoveIndvList:        removeIndvList,
		KeepFile:              *keepFile,
		RemoveFile:            *removeFile,
		GenoR2:                *genoR2,
		HapR2:                 *hapR2,
		GenoR2Positions:       *genoR2Positions,
		HapR2Positions:        *hapR2Positions,
		LDWindow:              *ldWindow,
		LDWindowBp:            *ldWindowBp,
		LDWindowMin:           *ldWindowMin,
		LDWindowBpMin:         *ldWindowBpMin,
		MinR2:                 *minR2,
		Bed:                   *bedFile,
		ExcludeBed:            *excludeBedFile,
		Diff:                  *diffFile,
		DiffSite:              *diffSite,
		DiffIndv:              *diffIndv,
		DiffSiteDiscordance:   *diffSiteDiscord,
		DiffIndvDiscordance:   *diffIndvDiscord,
		DiffIndvMap:           *diffIndvMap,
		DiffDiscordanceMatrix: *diffDiscMatrix,
		DiffSwitchError:       *diffSwitchError,
		MendelPedFile:         *mendelPed,
		BEAGLEGL:              *beagleGL,
		BEAGLEPL:              *beaglePL,
		InterchromGenoR2:      *interchromGenoR2,
		InterchromHapR2:       *interchromHapR2,
		GenoChiSq:             *genoChiSq,
		Relatedness:           *relatedness,
		Relatedness2:          *relatedness2,
		PhasedBlocks:          *phasedBlocks,
		LROH:                  *lroh,
		LROHMinVariants:       *lrohMin,
		RemoveFiltered:        *removeFiltered,
		KeepFiltered:          *keepFiltered,
		KeepINFO:              strings.Join(keepINFOParts, ","),
		RemoveINFO:            strings.Join(removeINFOParts, ","),
		GetINFO:               strings.Join(getINFOParts, ","),
		Phased:                *phased,
		LDhat:                 *ldhat,
		LDhatGeno:             *ldhatGeno,
		LDhelmet:              *ldhelmet,
		IMPUTE:                *impute,
		PCA:                   *pca,
		PCANoNorm:             *pcaNoNorm,
		PCASNPLoadings:        *pcaSNPLoadings,
	}

	// --hwe implies max_alleles = 2 in upstream (parameters.cpp:254).
	// Apply the same coupling here so the user's CLI invocation behaves
	// identically (e.g. `--vcf x.vcf --hwe 0.05 --recode` drops multi-
	// allelic sites even when --max-alleles is at its default of 0).
	if *hwePvalue > 0 {
		if params.MaxAlleles == 0 || params.MaxAlleles > 2 {
			params.MaxAlleles = 2
		}
	}

	// Run vcftools
	if err := vcftools.Run(inputReader, params); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v\n", err)
		os.Exit(1)
	}
}
