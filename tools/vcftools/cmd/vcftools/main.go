// vcftools - Utilities for working with VCF (Variant Call Format) files
package main

import (
	"flag"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
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
  --stdout, -c          Write output to stdout
  --recode              Output a new VCF file after filtering
  --recode-INFO-all     Include all INFO fields in recoded VCF
  --recode-INFO TAG     Keep this INFO tag in --recode output (repeatable;
                        recode-column selector matching upstream
                        parameters.cpp:319)

Position Filtering:
  --chr STRING          Include only this chromosome
  --not-chr STRING      Exclude this chromosome
  --from-bp INT         Include positions >= this value
  --to-bp INT           Include positions <= this value
  --positions FILE      Include positions listed in file
  --exclude-positions FILE  Exclude positions listed in file
  --positions-overlap FILE  Include records overlapping any position in file
  --exclude-positions-overlap FILE  Exclude records overlapping any position in file
  --bed FILE            Keep only sites inside any interval in BED FILE
  --exclude-bed FILE    Remove sites inside any interval in BED FILE
  --mask FILE           FASTA-style positional mask: keep sites where mask
                        digit <= --mask-min (default 0). Streams forward-only;
                        VCF must be sorted to match the mask's chromosome
                        order (mirrors upstream).
  --invert-mask FILE    Like --mask but with keep/drop inverted.
  --mask-min INT        Maximum kept mask digit value (0-9, default 0).

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
  --freq2               Allele frequency without allele labels (.frq)
  --counts2             Allele counts without allele labels (.frq.count)
  --derived             With --freq/--counts, reorder so ancestral (INFO/AA)
                        allele is first; drops sites without a matching AA
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
  --kept-sites          Output CHROM/POS of sites that pass filters (.kept.sites)
  --removed-sites       Output CHROM/POS of sites that fail filters (.removed.sites)

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
  --keep-INFO TAG       SITE FILTER: keep sites where this INFO Flag is
                        present (repeatable; OR across tags). Matches
                        upstream parameters.cpp:266 + entry_filters.cpp:1033.
                        To restrict the INFO column in recoded output use
                        --recode-INFO TAG instead.
  --remove-INFO TAG     SITE FILTER: drop sites where this INFO Flag is
                        present (repeatable; OR-veto across tags). Matches
                        upstream parameters.cpp:328 + entry_filters.cpp:1068.
                        Composes with --keep-INFO (keep narrows first, then
                        remove vetoes the survivors).
  --get-INFO TAG        Extract this INFO tag to <prefix>.INFO (repeatable)
  --extract-FORMAT-info NAME
                        Extract per-genotype FORMAT NAME into
                        <prefix>.<NAME>.FORMAT (tab-separated; sites
                        whose FORMAT lacks NAME are skipped)
  --indv-burden         Per-individual diploid-burden counts (.iburden):
                        N_HOM_REF, N_HET, N_HOM_ALT, N_MISS. With
                        --derived, columns become N_HOM_ANC/HET/HOM_DER.
  --indv-freq-burden    Per-individual frequency-burden matrix
                        (.ifreqburden): rows are kept individuals,
                        columns are allele-count bins 0..2N.
  --indv-freq-burden2   Like --indv-freq-burden but hom-alt genotypes
                        contribute 1 (not 2) to the per-bin count.
  --hapcount BED        Per-bin haplotype-count summaries (.hapcount).
                        Auto-detects BED headers (lines starting with
                        '#', 'track', 'browser', or blank lines).
                        Bins must be non-overlapping. Implies --phased.

Diff:
  --diff FILE           Second VCF for --diff-* outputs.
  --gzdiff FILE         Alias for --diff (iohelper auto-sniffs gzip).

Misc:
  --temp DIR            Spill-file directory (accepted for parity;
                        this port doesn't spill to disk so the flag
                        has no observable effect).

Sample Filtering:
  --indv STRING         Include only this individual (can use multiple times)
  --remove-indv STRING  Remove this individual (can use multiple times)
  --keep FILE           Keep only individuals listed in file
  --remove FILE         Remove individuals listed in file
  --max-indv INT        Cap the number of kept individuals at INT. Upstream
                        picks randomly; this port keeps the first N in
                        input order (see docs/PARITY_ROADMAP.md#vcftools).

Per-Genotype FT Filtering:
  --remove-filtered-geno-all
                        Set GT to ./. for any genotype whose FORMAT FT
                        is not "PASS" or ".".
  --remove-filtered-geno NAME
                        Set GT to ./. for any genotype whose FORMAT FT
                        lists NAME (repeatable).

Help:
  -h, --help            Show this help message
  --version             Print VCFtools version and exit
  --keep-INFO-all       Deprecated synonym for --recode-INFO-all

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
	bcfFile := flag.String("bcf", "", "Input BCF file (BGZF-compressed BCF v2.2)")
	contigsFile := flag.String("contigs", "", "Supplemental ##contig= lines for BCF header construction (consulted only when input lacks contig declarations)")
	useStdin := flag.Bool("stdin", false, "Read from stdin")

	// Output options
	outPrefix := flag.String("out", "out", "Output filename prefix")
	useStdout := flag.Bool("stdout", false, "Write to stdout")
	// -c is upstream's short alias for --stdout (parameters.cpp:194).
	// Wired to the same boolean so either spelling toggles streaming
	// output. Upstream's logic is `stream_out = true` from either flag;
	// last-set wins which is the implicit behaviour of flag.BoolVar.
	flag.BoolVar(useStdout, "c", false, "Write to stdout (short alias for --stdout)")
	recode := flag.Bool("recode", false, "Output a new VCF file")
	recodeBCF := flag.Bool("recode-bcf", false, "Output a new BCF file (BGZF-compressed BCF v2.2)")
	recodeInfoAll := flag.Bool("recode-INFO-all", false, "Include all INFO fields in recode")

	// Position filtering
	chr := flag.String("chr", "", "Include only this chromosome")
	notChr := flag.String("not-chr", "", "Exclude this chromosome")
	fromBp := flag.Int("from-bp", 0, "Include positions >= this value")
	toBp := flag.Int("to-bp", 0, "Include positions <= this value")
	positionsFile := flag.String("positions", "", "Include positions from file")
	excludePositionsFile := flag.String("exclude-positions", "", "Exclude positions from file")
	positionsOverlapFile := flag.String("positions-overlap", "", "Include records that overlap any position in file (sweeps POS..POS+len(REF)-1)")
	excludePositionsOverlapFile := flag.String("exclude-positions-overlap", "", "Exclude records that overlap any position in file")

	// SNP ID filtering
	snp := flag.String("snp", "", "Include only this SNP ID")
	snps := flag.String("snps", "", "Include SNP IDs from file")
	excludeSNP := flag.String("exclude", "", "Exclude this SNP ID")
	excludeSNPs := flag.String("exclude-snps", "", "Exclude SNP IDs from file")
	thin := flag.Int("thin", 0, "Thin sites by keeping every Nth site")

	// Variant type filtering
	keepOnlyIndels := flag.Bool("keep-only-indels", false, "Keep only indels")
	removeIndels := flag.Bool("remove-indels", false, "Remove indels")
	minAlleles := flag.Int("min-alleles", 0, "Minimum number of alleles")
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
	maxMissing := flag.Float64("max-missing", 0, "Minimum site call rate (1 = no missing genotypes allowed; 0 = no filter, the upstream default)")
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
	// --derived: when combined with --freq or --counts, reorder the allele
	// columns so the ancestral allele (INFO/AA, case-insensitive) appears
	// first. Sites missing AA, with AA = "." / "?", or with AA that does
	// not match REF/ALT are dropped (upstream emits a one-off warning).
	// Ported from parameters.cpp:201 + variant_file_output.cpp:67-159.
	derived := flag.Bool("derived", false, "With --freq/--counts, reorder so ancestral (INFO/AA) is first; drops sites without a matching AA")
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
	// --kept-sites / --removed-sites emit a 2-column (CHROM, POS) TSV
	// listing the sites that pass / fail filtering. Mirrors upstream
	// parameters.cpp:268, 330 + variant_file_output.cpp:4285-4373.
	keptSites := flag.Bool("kept-sites", false, "Output CHROM/POS of sites that pass filters (.kept.sites)")
	removedSites := flag.Bool("removed-sites", false, "Output CHROM/POS of sites that fail filters (.removed.sites)")

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

	// --pca / --pca-no-norm / --pca-snp-loadings INT: build the N×N
	// Genomic Relatedness Matrix from centred (and optionally
	// variance-normalised) genotypes, eigendecompose via gonum, and
	// emit `<prefix>.pca` (and optionally `<prefix>.pca.loadings` for
	// the per-site projections). See pca.go for the algorithm and
	// docs/UPSTREAM_BUGS.md for the missing-data fix-on-port.
	pca := flag.Bool("pca", false, "Principal component analysis → <prefix>.pca (eigendecomposition of the N×N GRM; implies biallelic-only)")
	pcaNoNorm := flag.Bool("pca-no-norm", false, "PCA without per-SNP variance normalisation (implies --pca)")
	pcaSNPLoadings := flag.Int("pca-snp-loadings", 0, "K → <prefix>.pca.loadings: per-site projection onto the first K principal components")

	// --phased: keep only sites where every kept-individual GT is phased.
	phased := flag.Bool("phased", false, "Keep only sites where every kept-individual GT is phased (separator '|' or haploid)")

	// BED-based filtering
	bedFile := flag.String("bed", "", "Keep only sites whose POS lies inside any interval in this BED file")
	excludeBedFile := flag.String("exclude-bed", "", "Remove sites whose POS lies inside any interval in this BED file")

	// FASTA-like positional mask filtering. --mask and --invert-mask are
	// mutually exclusive (upstream's parameters.cpp:262 / :280 set the same
	// mask_file slot, the last one wins); we surface both flags and OR the
	// last-set into Params. --mask-min defaults to 0 (drop everything but
	// digit '0'); valid range 0-9 (upstream parameters.cpp:720).
	maskFile := flag.String("mask", "", "FASTA-style positional mask: keep sites with mask digit <= --mask-min")
	invertMaskFile := flag.String("invert-mask", "", "FASTA-style positional mask with inverted keep/drop semantics")
	maskMin := flag.Int("mask-min", 0, "Maximum kept mask digit value (0-9, default 0)")

	// VCF comparison (--diff family)
	diffFile := flag.String("diff", "", "Second VCF file to compare against")
	diffBCFFile := flag.String("diff-bcf", "", "Second BCF file to compare against (BGZF-compressed BCF v2.2)")
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

	// INFO tag handling. All four flags are repeatable upstream; this
	// port accepts repeated single-tag invocations or one
	// comma-separated value, joined with commas in the same order seen
	// on the command line.
	//
	// Wave 17 fixes a semantic divergence:
	//   - `--keep-INFO TAG` (parameters.cpp:266) is upstream a SITE
	//     FILTER (`site_INFO_flags_to_keep` → entry_filters.cpp:1033).
	//     A site passes only if at least one of the named INFO Flags is
	//     present. Upstream additionally errors out at runtime if a
	//     listed tag is not declared as Type=Flag in the header.
	//   - `--recode-INFO TAG` (parameters.cpp:319) is the recode-column
	//     selector (`recode_INFO_to_keep`). It restricts the INFO column
	//     in `.recode.vcf` output to the listed tags.
	//
	// Pre-wave-17 the port wired `--keep-INFO` to the recode-column
	// selector semantic and treated `--recode-INFO` as a synonym for
	// it. Wave 17 separated the two into distinct slices:
	// `keepINFOParts` drives `params.KeepINFO` (site filter) and
	// `recodeINFOParts` drives `params.RecodeINFO` (recode-column
	// selector).
	//
	// Wave 18 repoints `--remove-INFO` from the port-only recode-column
	// stripper at upstream's SITE FILTER semantic
	// (parameters.cpp:328 → site_INFO_flags_to_remove →
	// entry_filters.cpp:1068-1086): drop sites where any of the named
	// Flag-type tags is present (OR-veto), with the same Type=Flag
	// header check as --keep-INFO. See docs/UPSTREAM_BUGS.md
	// `Fix-on-port` section for both migration notes.
	var keepINFOParts, recodeINFOParts, removeINFOParts, getINFOParts []string
	flag.Func("keep-INFO", "INFO Flag-type tag to use as a SITE FILTER (upstream parameters.cpp:266; repeatable)", func(s string) error {
		keepINFOParts = append(keepINFOParts, s)
		return nil
	})
	flag.Func("recode-INFO", "INFO tag to keep in --recode output (recode-column selector; upstream parameters.cpp:319; repeatable)", func(s string) error {
		recodeINFOParts = append(recodeINFOParts, s)
		return nil
	})
	flag.Func("remove-INFO", "INFO Flag-type tag to use as a SITE FILTER, OR-veto (upstream parameters.cpp:328; repeatable)", func(s string) error {
		removeINFOParts = append(removeINFOParts, s)
		return nil
	})
	flag.Func("get-INFO", "INFO tag to extract to <prefix>.INFO (repeatable)", func(s string) error {
		getINFOParts = append(getINFOParts, s)
		return nil
	})

	// --extract-FORMAT-info NAME extracts the named per-genotype FORMAT
	// field across all kept samples into a tab-separated
	// <prefix>.<NAME>.FORMAT file. Single-valued upstream (the last value
	// wins if supplied multiple times — parameters.cpp:222 simply
	// overwrites). Ported from variant_file_format_convert.cpp:1204-1263.
	extractFormatInfo := flag.String("extract-FORMAT-info", "", "Per-genotype FORMAT tag to extract to <prefix>.<NAME>.FORMAT")

	// Per-individual burden flags (parameters.cpp:257-259).
	// --indv-burden writes <prefix>.iburden; --indv-freq-burden and
	// --indv-freq-burden2 both write <prefix>.ifreqburden (the latter
	// with doubleCountHomAlt=1; see burden.go). All three are
	// modifier-free and combine cleanly with --derived.
	indvBurden := flag.Bool("indv-burden", false, "Per-individual diploid-burden counts (.iburden)")
	indvFreqBurden := flag.Bool("indv-freq-burden", false, "Per-individual frequency-burden matrix (.ifreqburden)")
	indvFreqBurden2 := flag.Bool("indv-freq-burden2", false, "Same as --indv-freq-burden but hom-alt contributes 1 (not 2)")

	// --hapcount BED — per-BED-bin haplotype-count summaries
	// (.hapcount). Upstream parameters.cpp:248 sets `phased_only = true`
	// so unphased sites are dropped before this output is computed
	// (matched in Run()). Three upstream bugs in
	// `output_haplotype_count` are fixed-on-port (see hapcount.go and
	// docs/UPSTREAM_BUGS.md): prev_bin_idx shift on bin change,
	// end-of-stream read-after-free, and BED first-line silent skip.
	hapcountBED := flag.String("hapcount", "", "BED file of bins for per-bin haplotype-count summaries (.hapcount); implies --phased")

	// --temp DIR — upstream parameters.cpp:341 stores a directory used
	// for `mkstemp` spill files in LD / format-convert paths. This port
	// does not spill, so the flag is accepted for CLI parity but has
	// no observable effect (Run() logs to stderr that the value was
	// parsed-but-unused). Documented in docs/PARITY_ROADMAP.md.
	tempDir := flag.String("temp", "", "Directory for spill files (accepted for parity; this port does not spill)")

	// --gzdiff FILE — upstream parameters.cpp:237 is identical to
	// `--diff` plus a `diff_file_compressed=true` flag that selects a
	// gzip reader. This port's `iohelper.OpenReader` auto-sniffs gzip
	// from the magic bytes, so `--gzdiff` is wired as an alias for
	// `--diff` (last-set wins). Documented in
	// docs/PARITY_ROADMAP.md#vcftools.
	gzDiff := flag.String("gzdiff", "", "Alias for --diff (this port auto-sniffs gzip via iohelper)")

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

	// --max-indv N: cap kept-individual count at N. We use flag.Func so we
	// can record whether the flag was supplied (since N == 0 is meaningful
	// — drop every sample — and Go's default zero would otherwise look the
	// same as "no flag given"). See Params.MaxIndv docstring for the
	// upstream-parity note about deterministic input-order truncation.
	var maxIndv int
	var maxIndvSet bool
	flag.Func("max-indv", "Cap the number of kept individuals at N (input-order truncation, see ROADMAP)", func(s string) error {
		v, err := strconv.Atoi(strings.TrimSpace(s))
		if err != nil {
			return fmt.Errorf("--max-indv: %w", err)
		}
		maxIndv = v
		maxIndvSet = true
		return nil
	})

	// Per-genotype FT-based filtering. --remove-filtered-geno-all matches
	// upstream parameters.cpp:323; --remove-filtered-geno NAME (repeatable)
	// matches parameters.cpp:324. Both set the kept sample's GT to ./.
	// while leaving other FORMAT fields untouched in recoded output (see
	// vcf_entry.cpp:580-608).
	removeFilteredGenoAll := flag.Bool("remove-filtered-geno-all", false, "Set genotype to ./. when FT is anything other than PASS or .")
	var removeFilteredGenoList []string
	flag.Func("remove-filtered-geno", "FT flag name to drop genotypes by (repeatable)", func(s string) error {
		removeFilteredGenoList = append(removeFilteredGenoList, s)
		return nil
	})

	// --keep-INFO-all is the upstream-deprecated synonym for
	// --recode-INFO-all (parameters.cpp:267 — "Old command (soon to be
	// depreciated)"). Both set the same parameter flag in upstream
	// (recode_all_INFO = true); we OR them together below.
	keepINFOAll := flag.Bool("keep-INFO-all", false, "Synonym for --recode-INFO-all (deprecated upstream)")

	// --version prints "VCFtools (<version>)" and exits, matching
	// upstream parameters.cpp:648-652. Handled before flag.Parse'd flags
	// take effect (the boolean below is checked immediately after parse).
	versionFlag := flag.Bool("version", false, "Print VCFtools version and exit")

	// Help. Upstream vcftools (parameters.cpp:654) treats all of "-h",
	// "-?", "-help", "--?", "--help", "--h" as the help trigger. Go's
	// flag package matches a registered name under either a single or a
	// double dash, so registering the bare names "h", "help", and "?"
	// accepts every upstream spelling (e.g. flag name "h" matches both
	// "-h" and "--h"; "?" matches "-?" and "--?"; "help" matches "-help"
	// and "--help"). All three point at the same boolean.
	help := flag.Bool("h", false, "Show help message")
	flag.BoolVar(help, "help", false, "Show help message")
	flag.BoolVar(help, "?", false, "Show help message (upstream alias)")

	flag.Parse()

	if *versionFlag {
		// Match upstream's exact byte sequence (parameters.cpp:648-652).
		// The version string is hard-coded here to track upstream's
		// VCFTOOLS_VERSION at the time of port; bump alongside any
		// future re-baseline.
		fmt.Println("VCFtools (0.1.18)")
		os.Exit(0)
	}

	if *help {
		fmt.Fprint(os.Stderr, usage)
		os.Exit(0)
	}

	// --keep-INFO-all is the deprecated synonym for --recode-INFO-all.
	// Upstream parameters.cpp:267 / :318 both write to the same
	// recode_all_INFO bit; we OR them together so either flag (or both)
	// produces the same effect.
	recodeAllINFO := *recodeInfoAll || *keepINFOAll

	// Validate input
	inputCount := 0
	if *vcfFile != "" {
		inputCount++
	}
	if *gzvcfFile != "" {
		inputCount++
	}
	if *bcfFile != "" {
		inputCount++
	}
	if *useStdin {
		inputCount++
	}

	if inputCount == 0 {
		fmt.Fprintln(os.Stderr, "Error: Must specify one of --vcf, --gzvcf, --bcf, or --stdin")
		fmt.Fprintln(os.Stderr, "Use --help for usage information")
		os.Exit(1)
	}
	if inputCount > 1 {
		fmt.Fprintln(os.Stderr, "Error: Can only specify one input method")
		os.Exit(1)
	}

	// Determine input file. --bcf takes a distinct code path inside
	// Run (BGZF + BCF decoder); the io.Reader argument is unused in
	// that case, so we open /dev/null-equivalent (an empty
	// strings.Reader) below to satisfy the signature.
	var inputFile string
	if *vcfFile != "" {
		inputFile = *vcfFile
	} else if *gzvcfFile != "" {
		inputFile = *gzvcfFile
	} else {
		inputFile = "" // stdin or BCF (handled by Run via BCFInputFile)
	}

	// Open input. For --bcf we skip this entirely; Run opens the BCF
	// file directly when params.BCFInputFile is set.
	var inputReader io.ReadCloser
	if *bcfFile == "" {
		r, err := iohelper.OpenReader(inputFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "Error opening input: %v\n", err)
			os.Exit(1)
		}
		defer r.Close()
		inputReader = r
	} else {
		inputReader = io.NopCloser(strings.NewReader(""))
	}

	// Build params
	params := &vcftools.Params{
		OutPrefix:                   *outPrefix,
		UseStdout:                   *useStdout,
		Recode:                      *recode,
		RecodeBCF:                   *recodeBCF,
		BCFInputFile:                *bcfFile,
		ContigsFile:                 *contigsFile,
		RecodeInfoAll:               recodeAllINFO,
		Chr:                         *chr,
		NotChr:                      *notChr,
		FromBp:                      *fromBp,
		ToBp:                        *toBp,
		PositionsFile:               *positionsFile,
		ExcludePositionsFile:        *excludePositionsFile,
		PositionsOverlapFile:        *positionsOverlapFile,
		ExcludePositionsOverlapFile: *excludePositionsOverlapFile,
		SNP:                         *snp,
		SNPs:                        *snps,
		ExcludeSNP:                  *excludeSNP,
		ExcludeSNPs:                 *excludeSNPs,
		Thin:                        *thin,
		KeepOnlyIndels:              *keepOnlyIndels,
		RemoveIndels:                *removeIndels,
		MinAlleles:                  *minAlleles,
		MaxAlleles:                  *maxAlleles,
		MinQ:                        *minQ,
		RemoveFilteredAll:           *removeFilteredAll,
		Maf:                         *maf,
		MaxMaf:                      *maxMaf,
		Mac:                         *mac,
		MaxMac:                      *maxMac,
		MinNonRefAF:                 *nonRefAF,
		MinNonRefAC:                 *nonRefAC,
		MaxNonRefAF:                 *maxNonRefAF,
		MaxNonRefAC:                 *maxNonRefAC,
		MinNonRefAFAny:              *nonRefAFAny,
		MinNonRefACAny:              *nonRefACAny,
		MaxNonRefAFAny:              *maxNonRefAFAny,
		MaxNonRefACAny:              *maxNonRefACAny,
		MaxMissing:                  *maxMissing,
		MaxMissingCount:             maxMissingCount,
		MaxMissingCountSet:          maxMissingCountSet,
		MinHWEPvalue:                *hwePvalue,
		MinMeanDP:                   *minMeanDP,
		MaxMeanDP:                   *maxMeanDP,
		MinDP:                       *minDP,
		MaxDP:                       *maxDP,
		MinGQ:                       *minGQ,
		Freq:                        *freq,
		Counts:                      *counts,
		Freq2:                       *freq2,
		Counts2:                     *counts2,
		Depth:                       *depth,
		SiteDepth:                   *siteDepth,
		SiteMeanDepth:               *siteMeanDepth,
		SiteQuality:                 *siteQuality,
		MissingIndv:                 *missingIndv,
		MissingSite:                 *missingSite,
		Hardy:                       *hardy,
		TsTvSummary:                 *tsTvSummary,
		TsTvBinSize:                 *tsTvBinSize,
		TsTvByCount:                 *tsTvByCount,
		TsTvByQual:                  *tsTvByQual,
		SitePi:                      *sitePi,
		Het:                         *het,
		Singletons:                  *singletons,
		HistIndelLen:                *histIndelLen,
		GenoDepth:                   *genoDepth,
		WindowPi:                    *windowPi,
		WindowPiStep:                *windowPiStep,
		TajimaD:                     *tajimaD,
		SNPDensity:                  *snpDensity,
		WeirFstPop:                  weirFstPop,
		FstWindowSize:               *fstWindowSize,
		FstWindowStep:               *fstWindowStep,
		FilterSummary:               *filterSummary,
		Output012:                   *output012,
		OutputPlink:                 *outputPlink,
		OutputPlinkTped:             *outputPlinkTped,
		ChromMap:                    *chromMap,
		IndvList:                    indvList,
		RemoveIndvList:              removeIndvList,
		KeepFile:                    *keepFile,
		RemoveFile:                  *removeFile,
		GenoR2:                      *genoR2,
		HapR2:                       *hapR2,
		GenoR2Positions:             *genoR2Positions,
		HapR2Positions:              *hapR2Positions,
		LDWindow:                    *ldWindow,
		LDWindowBp:                  *ldWindowBp,
		LDWindowMin:                 *ldWindowMin,
		LDWindowBpMin:               *ldWindowBpMin,
		MinR2:                       *minR2,
		Bed:                         *bedFile,
		ExcludeBed:                  *excludeBedFile,
		Mask:                        *maskFile,
		InvertMask:                  false,
		MaskMin:                     *maskMin,
		Diff:                        *diffFile,
		DiffBCF:                     *diffBCFFile,
		DiffSite:                    *diffSite,
		DiffIndv:                    *diffIndv,
		DiffSiteDiscordance:         *diffSiteDiscord,
		DiffIndvDiscordance:         *diffIndvDiscord,
		DiffIndvMap:                 *diffIndvMap,
		DiffDiscordanceMatrix:       *diffDiscMatrix,
		DiffSwitchError:             *diffSwitchError,
		MendelPedFile:               *mendelPed,
		BEAGLEGL:                    *beagleGL,
		BEAGLEPL:                    *beaglePL,
		InterchromGenoR2:            *interchromGenoR2,
		InterchromHapR2:             *interchromHapR2,
		GenoChiSq:                   *genoChiSq,
		Relatedness:                 *relatedness,
		Relatedness2:                *relatedness2,
		PhasedBlocks:                *phasedBlocks,
		LROH:                        *lroh,
		LROHMinVariants:             *lrohMin,
		RemoveFiltered:              *removeFiltered,
		KeepFiltered:                *keepFiltered,
		KeepINFO:                    strings.Join(keepINFOParts, ","),
		RecodeINFO:                  strings.Join(recodeINFOParts, ","),
		RemoveINFO:                  strings.Join(removeINFOParts, ","),
		GetINFO:                     strings.Join(getINFOParts, ","),
		Phased:                      *phased,
		LDhat:                       *ldhat,
		LDhatGeno:                   *ldhatGeno,
		LDhelmet:                    *ldhelmet,
		IMPUTE:                      *impute,
		PCA:                         *pca,
		PCANoNorm:                   *pcaNoNorm,
		PCASNPLoadings:              *pcaSNPLoadings,
		KeptSites:                   *keptSites,
		RemovedSites:                *removedSites,
		MaxIndv:                     maxIndv,
		MaxIndvSet:                  maxIndvSet,
		RemoveFilteredGenoAll:       *removeFilteredGenoAll,
		RemoveFilteredGenoList:      removeFilteredGenoList,
		Derived:                     *derived,
		ExtractFormatInfo:           *extractFormatInfo,
		IndvBurden:                  *indvBurden,
		IndvFreqBurden:              *indvFreqBurden,
		IndvFreqBurden2:             *indvFreqBurden2,
		HapcountBED:                 *hapcountBED,
		TempDir:                     *tempDir,
	}

	// --gzdiff FILE: upstream-compatible alias for --diff. The
	// underlying iohelper.OpenReader auto-sniffs gzip from the magic
	// bytes, so the compressed/uncompressed distinction is irrelevant.
	// If both --diff and --gzdiff are supplied, --gzdiff wins (matches
	// upstream's "last-set" semantics for the shared `diff_file` slot —
	// parameters.cpp:209 vs :237 both write the same slot, last parsed
	// wins).
	if *gzDiff != "" {
		params.Diff = *gzDiff
	}

	// --diff-bcf is mutually exclusive with --diff / --gzdiff: upstream
	// writes them all to the same `diff_file` slot
	// (parameters.cpp:209,237) so only the last-parsed flag wins.
	//
	// DIVERGENCE: Go's flag package erases declaration order, so we
	// can't replicate upstream's true last-set-wins. Instead, when
	// both flags are present we prefer --diff-bcf unconditionally
	// (the BCF flag is the more specific intent) and emit a stderr
	// warning so the user knows which one we picked. To get true
	// upstream parity we'd need flag.Visit() to recover the
	// declaration order — not worth the complexity for this rarely
	// double-passed case.
	if params.DiffBCF != "" && params.Diff != "" {
		fmt.Fprintln(os.Stderr, "warning: both --diff and --diff-bcf set; --diff-bcf takes precedence")
		params.Diff = ""
	}

	// --invert-mask FILE shares the same mask_file slot upstream as --mask
	// (parameters.cpp:262 vs :280); whichever appears last wins. Go's
	// flag.String evaluation does not preserve last-set ordering when both
	// flags are present, so we just take --invert-mask as the override when
	// non-empty. Document the limitation in PARITY_ROADMAP.
	if *invertMaskFile != "" {
		params.Mask = *invertMaskFile
		params.InvertMask = true
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
