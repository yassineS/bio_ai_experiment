# vcftools Feature Comparison: Original vs Go Implementation

This document compares the original vcftools (C++/Perl) with the Go implementation to identify implemented and missing features.

## Summary

- **Original vcftools**: ~147 command-line options
- **Go implementation**: ~70 command-line options (after long-tail wave 1)
- **Coverage**: ~48% of total options, ~85% of commonly used features

## Implemented Features (✅)

### Input/Output

- ✅ `--vcf` - Input VCF file
- ✅ `--gzvcf` - Input gzipped VCF file
- ✅ `--stdin` - Read from stdin
- ✅ `--out` - Output prefix
- ✅ `--stdout` - Write to stdout
- ✅ `--recode` - Output filtered VCF
- ✅ `--recode-INFO-all` - Include all INFO fields in recode

### Position Filtering

- ✅ `--chr` - Filter by chromosome
- ✅ `--not-chr` - Exclude chromosome
- ✅ `--from-bp` - Minimum position
- ✅ `--to-bp` - Maximum position
- ✅ `--positions` - Include positions from file
- ✅ `--exclude-positions` - Exclude positions from file
- ✅ `--bed` - Include sites inside any BED interval
- ✅ `--exclude-bed` - Exclude sites inside any BED interval

### Variant Type Filtering

- ✅ `--keep-only-indels` - Keep only indels
- ✅ `--remove-indels` - Remove indels
- ✅ `--min-alleles` - Minimum number of alleles
- ✅ `--max-alleles` - Maximum number of alleles

### Quality Filtering

- ✅ `--minQ` - Minimum quality score
- ✅ `--remove-filtered-all` - Remove non-PASS sites

### Allele Frequency Filtering

- ✅ `--maf` - Minimum minor allele frequency
- ✅ `--max-maf` - Maximum minor allele frequency
- ✅ `--mac` - Minimum minor allele count
- ✅ `--max-mac` - Maximum minor allele count

### Genotype Filtering

- ✅ `--max-missing` - Maximum missing data proportion
- ✅ `--min-meanDP` - Minimum mean depth
- ✅ `--max-meanDP` - Maximum mean depth

### Sample Filtering

- ✅ `--indv` - Include individual(s)
- ✅ `--remove-indv` - Remove individual(s)
- ✅ `--keep` - Keep individuals from file
- ✅ `--remove` - Remove individuals from file

### Statistics Output

- ✅ `--freq` - Allele frequency
- ✅ `--counts` - Allele counts
- ✅ `--depth` - Mean depth per site
- ✅ `--site-depth` - Depth for each site
- ✅ `--site-mean-depth` - Mean depth per site
- ✅ `--site-quality` - Quality per site
- ✅ `--missing-indv` - Individual missingness
- ✅ `--missing-site` - Site missingness
- ✅ `--hardy` - Hardy-Weinberg equilibrium
- ✅ `--TsTv-summary` - Ts/Tv ratio summary
- ✅ `--TsTv` - Ts/Tv in bins
- ✅ `--site-pi` - Nucleotide diversity per site

### VCF Comparison (--diff family)

- ✅ `--diff FILE` - Compare against a second VCF (`.gz` auto-detected)
- ✅ `--diff-site` - `<prefix>.diff.sites_in_files`
- ✅ `--diff-indv` - `<prefix>.diff.indv_in_files`
- ✅ `--diff-site-discordance` - `<prefix>.diff.sites` (upstream 7-column
  layout: CHROM, POS, FILES, MATCHING_ALLELES, N_COMMON_CALLED, N_DISCORD,
  DISCORDANCE; file-1-only and file-2-only sites included)
- ✅ `--diff-indv-discordance` - `<prefix>.diff.indv` (4-column layout
  with DISCORDANCE over the union of both files' samples)

### BEAGLE Genotype-Likelihood Output

- ✅ `--BEAGLE-GL` - `<prefix>.BEAGLE.GL` (log10 GL triplets from FORMAT/PL)
- ✅ `--BEAGLE-PL` - `<prefix>.BEAGLE.PL` (raw PL triplets)

### Linkage Disequilibrium

- ✅ `--geno-r2` - Genotype-based LD (r²) → `<prefix>.geno.ld`
- ✅ `--hap-r2` - Haplotype-based LD (r²) for phased data → `<prefix>.hap.ld`
- ✅ `--geno-r2-positions FILE` - Restrict --geno-r2 to pairs touching a listed position
- ✅ `--hap-r2-positions FILE` - Restrict --hap-r2 to pairs touching a listed position
- ✅ `--ld-window INT` - Maximum number of SNPs between LD pairs
- ✅ `--ld-window-bp INT` - Maximum bp distance between LD pairs
- ✅ `--ld-window-min INT` - Minimum number of SNPs between LD pairs
- ✅ `--ld-window-bp-min INT` - Minimum bp distance between LD pairs
- ✅ `--min-r2 FLOAT` - Minimum r² threshold for output

### Inter-chromosomal LD + Chi-square

- ✅ `--interchrom-geno-r2` - Inter-chromosomal genotype LD → `<prefix>.interchrom.geno.ld`
- ✅ `--interchrom-hap-r2` - Inter-chromosomal haplotype LD → `<prefix>.interchrom.hap.ld`
- ✅ `--geno-chisq` - Per-pair Pearson chi-square test → `<prefix>.geno.chisq`

### Relatedness & Homozygosity

- ✅ `--relatedness` - Yang 2010 unadjusted A_jk → `<prefix>.relatedness`
- ✅ `--relatedness2` - KING-robust kinship → `<prefix>.relatedness2`
- ✅ `--LROH` - Runs of homozygosity per individual → `<prefix>.LROH`
- ✅ `--phased-blocks` - Per-individual phased-haplotype blocks → `<prefix>.blocks`

### FILTER / INFO Selection

- ✅ `--remove-filtered NAME[,NAME...]` - Drop sites listing any named FILTER
- ✅ `--keep-filtered NAME[,NAME...]` - Keep only sites listing any named FILTER
- ✅ `--keep-INFO TAG` (SITE FILTER, upstream parameters.cpp:266 +
  entry_filters.cpp:1033; repeatable, OR-composing)
- ✅ `--remove-INFO TAG` (SITE FILTER, upstream parameters.cpp:328 +
  entry_filters.cpp:1068; repeatable, OR-veto; composes with
  `--keep-INFO`)
- ✅ `--get-INFO TAG[,TAG...]` - Extract INFO tags to `<prefix>.INFO`

## Missing Features (❌)

### High Priority Features (commonly used)

#### Linkage Disequilibrium Analysis

(All `--geno-chisq`, `--interchrom-*-r2` now implemented; see above.)

#### Population Genetics Statistics

- ❌ `--weir-fst-pop` - Fst calculation (Weir & Cockerham)
- ❌ `--fst-window-size` - Window size for Fst
- ❌ `--fst-window-step` - Step size for Fst windows
- ❌ `--window-pi` - Nucleotide diversity in windows
- ❌ `--window-pi-step` - Step size for pi windows
- ❌ `--TajimaD` - Tajima's D statistic
- ❌ `--het` - Heterozygosity statistics
- ❌ `--singletons` - Singleton site analysis

#### Format Conversion

- ❌ `--plink` - Convert to PLINK format
- ❌ `--plink-tped` - Convert to PLINK transposed format
- ❌ `--chrom-map` - Chromosome mapping for PLINK
- ❌ `--012` - Output 0/1/2 matrix
- ❌ `--IMPUTE` - Output IMPUTE format
- ❌ `--ldhat` - Output LDhat format
- ❌ `--ldhat-geno` - Output LDhat genotype format
- ❌ `--ldhelmet` - Output LDhelmet format
- ✅ `--BEAGLE-GL` - Output BEAGLE genotype likelihoods (GL)
- ✅ `--BEAGLE-PL` - Output BEAGLE genotype likelihoods (PL)

#### VCF Comparison/Diff

- ✅ `--diff` - Compare with another VCF (supports `.gz` via auto-detect)
- ❌ `--gzdiff` - Explicit gzipped-diff alias (the regular `--diff` already
  decompresses `.gz`)
- ✅ `--diff-bcf` - Compare with BCF (routes through the shared
  variantSource adapter; composes with --diff-indv-map and every
  --diff-* output)
- ✅ `--diff-site` - Sites in common/unique (`.diff.sites_in_files`)
- ✅ `--diff-indv` - Individuals in common/unique (`.diff.indv_in_files`)
- ✅ `--diff-site-discordance` - Site-by-site discordance (`.diff.sites`)
- ✅ `--diff-indv-discordance` - Individual discordance (`.diff.indv`)
- ❌ `--diff-discordance-matrix` - Discordance matrix
- ❌ `--diff-switch-error` - Phasing switch errors
- ❌ `--diff-indv-map` - Map individuals between files

### Medium Priority Features

#### Additional Filtering

- ❌ `--snp` - Include specific SNP by ID
- ❌ `--snps` - Include SNPs from file
- ❌ `--exclude` - Exclude SNPs from file
- ✅ `--bed` - Include positions from BED file
- ✅ `--exclude-bed` - Exclude positions from BED file
- ❌ `--positions-overlap` - Overlap-based position inclusion
- ❌ `--exclude-positions-overlap` - Overlap-based exclusion
- ❌ `--mask` - Mask file filtering
- ❌ `--invert-mask` - Inverted mask filtering
- ❌ `--mask-min` - Mask minimum threshold
- ❌ `--thin` - Thin sites by distance
- ✅ `--keep-filtered NAME[,NAME...]` - Keep specific FILTER flags
- ✅ `--remove-filtered NAME[,NAME...]` - Remove specific FILTER flags
- ✅ `--keep-INFO TAG` - SITE FILTER: keep only sites where the
  named Flag-type INFO tag is present (upstream parameters.cpp:266)
- ✅ `--remove-INFO TAG` - SITE FILTER: drop sites where the named
  Flag-type INFO tag IS present (upstream parameters.cpp:328;
  polarity-inverted complement of `--keep-INFO`)
- ❌ `--keep-INFO-all` - Keep all INFO fields (use `--recode-INFO-all`)

#### Genotype-Level Filtering

- ❌ `--minDP` - Minimum depth per genotype
- ❌ `--maxDP` - Maximum depth per genotype
- ❌ `--minGQ` - Minimum genotype quality
- ❌ `--remove-filtered-geno` - Remove filtered genotypes
- ❌ `--remove-filtered-geno-all` - Remove all filtered genotypes

#### Advanced Allele Filters

- ❌ `--non-ref-ac` - Non-reference allele count
- ❌ `--max-non-ref-ac` - Maximum non-ref allele count
- ❌ `--non-ref-ac-any` - Any non-ref allele count
- ❌ `--max-non-ref-ac-any` - Maximum any non-ref count
- ❌ `--non-ref-af` - Non-reference allele frequency
- ❌ `--max-non-ref-af` - Maximum non-ref frequency
- ❌ `--non-ref-af-any` - Any non-ref frequency
- ❌ `--max-non-ref-af-any` - Maximum any non-ref frequency
- ❌ `--max-missing-count` - Maximum missing count

#### Additional Statistics

- ❌ `--freq2` - Alternative frequency output
- ❌ `--counts2` - Alternative counts output
- ❌ `--geno-depth` - Genotype depth distribution
- ❌ `--hist-indel-len` - Indel length histogram
- ❌ `--SNPdensity` - SNP density
- ❌ `--TsTv-by-count` - Ts/Tv by allele count
- ❌ `--TsTv-by-qual` - Ts/Tv by quality
- ❌ `--FILTER-summary` - FILTER tag summary
- ❌ `--hapcount` - Haplotype counts
- ❌ `--mendel` - Mendelian error check
- ❌ `--indv-burden` - Individual burden
- ❌ `--indv-freq-burden` - Individual frequency burden
- ❌ `--indv-freq-burden2` - Alternative frequency burden

#### Haplotype/Phase Analysis

- ❌ `--phased` - Work with phased data only
- ✅ `--LROH` - Runs of homozygosity → `<prefix>.LROH`
- ✅ `--phased-blocks` - Per-individual phased-haplotype blocks → `<prefix>.blocks`

#### Advanced Analysis

- ✅ `--relatedness` - Relatedness analysis (Yang 2010)
- ✅ `--relatedness2` - KING-robust kinship
- ✅ `--pca` - Principal component analysis → `<prefix>.pca`
  (eigendecomposition of the centred/normalised N×N GRM via gonum)
- ✅ `--pca-no-norm` - PCA without per-SNP variance normalisation
  (still mean-centres; implies `--pca`)
- ✅ `--pca-snp-loadings INT` - Per-site loadings on the first K
  principal components → `<prefix>.pca.loadings`

#### Format Info Extraction

- ✅ `--get-INFO TAG[,TAG...]` - Extract INFO field values → `<prefix>.INFO`
- ❌ `--extract-FORMAT-info` - Extract FORMAT field values
- ✅ `--recode-INFO TAG` - Recode-column selector (upstream
  parameters.cpp:319, distinct from `--keep-INFO` since wave 17)

#### Output Management

- ❌ `--kept-sites` - Output list of kept sites
- ❌ `--removed-sites` - Output list of removed sites
- ❌ `--temp` - Temporary directory

#### BCF Support

- ✅ `--bcf` - Input BCF file (BGZF-decompressed and decoded via the
  shared `pkg/bioformats/bcf` reader; composes with the full filter
  pipeline)
- ✅ `--recode-bcf` - Output BCF format (BGZF-compressed BCF v2.2,
  interop-tested against upstream's `--bcf` reader)

#### Miscellaneous

- ✅ `--contigs` - Supplemental `##contig=` lines for BCF header
  construction; only consulted when the source lacks contig
  declarations (matches upstream gating)
- ❌ `--derived` - Derived allele frequency
- ❌ `--hwe` - Alternative HWE test
- ❌ `--max-indv` - Maximum number of individuals
- ❌ `--version` - Show version

## Feature Coverage Analysis

### By Category

| Category | Implemented | Total | Coverage |
|----------|-------------|-------|----------|
| Input/Output | 7 | 10 | 70% |
| Position Filtering | 6 | 13 | 46% |
| Variant Type Filtering | 4 | 6 | 67% |
| Quality/Genotype Filtering | 5 | 12 | 42% |
| Allele Frequency Filtering | 4 | 12 | 33% |
| Sample Filtering | 4 | 6 | 67% |
| Statistics Output | 12 | 28 | 43% |
| Format Conversion | 2 | 12 | 17% |
| LD Analysis | 0 | 12 | 0% |
| Population Genetics | 1 | 8 | 13% |
| VCF Comparison | 0 | 10 | 0% |
| Advanced Analysis | 0 | 4 | 0% |

### Usage Priority Assessment

Based on typical vcftools usage patterns:

**High Priority (Commonly Used):**

- ✅ Basic filtering (position, quality, allele frequency) - **IMPLEMENTED**
- ✅ Basic statistics (frequency, depth, missingness, HWE, Ts/Tv) - **IMPLEMENTED**
- ✅ VCF recoding - **IMPLEMENTED**
- ❌ LD analysis (r²) - **NOT IMPLEMENTED**
- ❌ Fst statistics - **NOT IMPLEMENTED**
- ❌ PLINK conversion - **NOT IMPLEMENTED**

**Medium Priority (Occasionally Used):**

- ❌ VCF comparison/diff - **NOT IMPLEMENTED**
- ❌ Windowed statistics - **NOT IMPLEMENTED**
- ❌ Format conversion (BEAGLE, IMPUTE, LDhat) - **NOT IMPLEMENTED**
- ❌ Relatedness analysis - **NOT IMPLEMENTED**

**Low Priority (Rarely Used):**

- ✅ PCA analysis (`--pca`, `--pca-no-norm`, `--pca-snp-loadings INT`) — landed in wave 19 via gonum
- ❌ Mendelian error checking - **NOT IMPLEMENTED**
- ❌ Advanced haplotype analysis - **NOT IMPLEMENTED**

## Recommendations for Feature Parity

### Phase 1: Critical Missing Features (High Impact)

1. **Linkage Disequilibrium Analysis**
   - Implement `--geno-r2` (most commonly used)
   - Implement `--hap-r2` (for phased data)
   - Add LD window options
   - Priority: HIGH

2. **Population Genetics Statistics**
   - Implement `--weir-fst-pop` (Fst calculation)
   - Implement windowed pi (`--window-pi`)
   - Implement Tajima's D (`--TajimaD`)
   - Priority: HIGH

3. **PLINK Format Conversion**
   - Implement `--plink` and `--plink-tped`
   - Add chromosome mapping support
   - Priority: MEDIUM-HIGH

### Phase 2: Important Features

4. **VCF Comparison/Diff**
   - Implement basic diff operations
   - Site and individual concordance
   - Priority: MEDIUM

2. **Additional Format Conversions**
   - BEAGLE format (for imputation)
   - 012 matrix (for quick analysis)
   - Priority: MEDIUM

3. **Windowed Statistics**
   - Window-based nucleotide diversity
   - SNP density in windows
   - Priority: MEDIUM

### Phase 3: Nice-to-Have Features

7. **Advanced Filtering**
   - BED file filtering
   - Mask-based filtering
   - SNP ID filtering
   - Priority: LOW-MEDIUM

2. **Advanced Analysis**
   - Relatedness calculation
   - Heterozygosity statistics
   - Priority: LOW

3. **BCF Support**
   - Input/output BCF format
   - Priority: LOW

## Conclusion

The current Go implementation covers approximately **80% of commonly used features** but only **24% of total features**. The implementation successfully handles the core use cases:

✅ **Well Covered:**

- Basic filtering operations
- Essential statistics
- VCF recoding
- Sample management

❌ **Major Gaps:**

- Linkage Disequilibrium analysis (critical for population genetics)
- Fst and population structure statistics
- Format conversion (PLINK, BEAGLE, etc.)
- VCF comparison/diff operations
- Windowed statistics

For true feature parity with the original vcftools, the most important additions would be:

1. LD analysis (`--geno-r2`, `--hap-r2`)
2. Fst statistics (`--weir-fst-pop`)
3. PLINK conversion (`--plink`)
4. Windowed statistics (`--window-pi`, `--TajimaD`)
5. VCF comparison (`--diff` family)

These five feature sets would bring coverage up to approximately 90% of real-world use cases while still representing a manageable implementation scope.
