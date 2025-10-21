# vcftools 100% Feature Parity Implementation Plan

This document tracks the implementation of all 112 missing features from the original vcftools.

**Current Status**: 35/147 features (24%)  
**Target**: 147/147 features (100%)  
**Remaining**: 112 features

---

## Phase 1: Additional Filtering & Basic Statistics (20 features)

### Additional Position Filtering (4 features)
- [ ] `--bed` - Include positions from BED file
- [ ] `--exclude-bed` - Exclude positions from BED file
- [ ] `--positions-overlap` - Overlap-based position inclusion
- [ ] `--exclude-positions-overlap` - Overlap-based exclusion

### SNP ID Filtering (4 features)
- [ ] `--snp` - Include specific SNP by ID
- [ ] `--snps` - Include SNPs from file
- [ ] `--exclude` - Exclude SNPs from file
- [ ] `--thin` - Thin sites by distance

### Advanced Quality/Genotype Filters (6 features)
- [ ] `--minDP` - Minimum depth per genotype
- [ ] `--maxDP` - Maximum depth per genotype
- [ ] `--minGQ` - Minimum genotype quality
- [ ] `--remove-filtered` - Remove specific FILTER flags
- [ ] `--keep-filtered` - Keep specific FILTER flags
- [ ] `--keep-INFO` - Keep sites with INFO flag
- [ ] `--remove-INFO` - Remove sites with INFO flag

### Additional Statistics (6 features)
- [ ] `--het` - Heterozygosity statistics
- [ ] `--singletons` - Singleton site analysis
- [ ] `--freq2` - Alternative frequency output
- [ ] `--counts2` - Alternative counts output
- [ ] `--hist-indel-len` - Indel length histogram
- [ ] `--geno-depth` - Genotype depth distribution

**Estimated Time**: 3-4 days  
**Priority**: HIGH

---

## Phase 2: Population Genetics Statistics (10 features)

### Windowed Statistics (4 features)
- [ ] `--window-pi` - Nucleotide diversity in windows
- [ ] `--window-pi-step` - Step size for pi windows
- [ ] `--TajimaD` - Tajima's D statistic
- [ ] `--SNPdensity` - SNP density in windows

### Fst Statistics (3 features)
- [ ] `--weir-fst-pop` - Fst calculation (Weir & Cockerham)
- [ ] `--fst-window-size` - Window size for Fst
- [ ] `--fst-window-step` - Step size for Fst windows

### Additional Ts/Tv (2 features)
- [ ] `--TsTv-by-count` - Ts/Tv by allele count
- [ ] `--TsTv-by-qual` - Ts/Tv by quality score

### Other (1 feature)
- [ ] `--FILTER-summary` - FILTER tag summary

**Estimated Time**: 5-6 days  
**Priority**: HIGH

---

## Phase 3: Linkage Disequilibrium Analysis (12 features)

### Basic LD Calculations (4 features)
- [ ] `--geno-r2` - Genotype-based LD (r²)
- [ ] `--hap-r2` - Haplotype-based LD (r²)
- [ ] `--geno-chisq` - Genotype chi-square test
- [ ] `--interchrom-geno-r2` - Inter-chromosomal genotype LD
- [ ] `--interchrom-hap-r2` - Inter-chromosomal haplotype LD

### Position-Specific LD (2 features)
- [ ] `--geno-r2-positions` - LD for specific positions
- [ ] `--hap-r2-positions` - Haplotype LD for specific positions

### LD Window Options (5 features)
- [ ] `--ld-window` - LD window size (SNPs)
- [ ] `--ld-window-bp` - LD window size (bp)
- [ ] `--ld-window-min` - Minimum LD window (SNPs)
- [ ] `--ld-window-bp-min` - Minimum LD window (bp)
- [ ] `--min-r2` - Minimum r² threshold

**Estimated Time**: 7-8 days  
**Priority**: HIGH  
**Complexity**: HIGH (requires interval trees, complex statistics)

---

## Phase 4: Format Conversions (12 features)

### PLINK Format (3 features)
- [ ] `--plink` - Convert to PLINK PED/MAP
- [ ] `--plink-tped` - Convert to PLINK transposed format
- [ ] `--chrom-map` - Chromosome name mapping

### Genotype Matrix (1 feature)
- [ ] `--012` - Output 0/1/2 matrix

### Imputation Formats (4 features)
- [ ] `--BEAGLE-GL` - BEAGLE genotype likelihoods (GL)
- [ ] `--BEAGLE-PL` - BEAGLE genotype likelihoods (PL)
- [ ] `--IMPUTE` - IMPUTE format
- [ ] `--ldhat` - LDhat format
- [ ] `--ldhat-geno` - LDhat genotype format
- [ ] `--ldhelmet` - LDhelmet format

### Additional INFO Extraction (2 features)
- [ ] `--get-INFO` - Extract INFO field values
- [ ] `--extract-FORMAT-info` - Extract FORMAT field values
- [ ] `--recode-INFO` - Recode with specific INFO fields

**Estimated Time**: 6-7 days  
**Priority**: MEDIUM-HIGH

---

## Phase 5: VCF Comparison/Diff (10 features)

### Input Options (3 features)
- [ ] `--diff` - Compare with another VCF
- [ ] `--gzdiff` - Compare with gzipped VCF
- [ ] `--diff-bcf` - Compare with BCF

### Comparison Operations (7 features)
- [ ] `--diff-site` - Sites in common/unique
- [ ] `--diff-indv` - Individuals in common/unique
- [ ] `--diff-site-discordance` - Site-by-site discordance
- [ ] `--diff-indv-discordance` - Individual discordance
- [ ] `--diff-discordance-matrix` - Discordance matrix
- [ ] `--diff-switch-error` - Phasing switch errors
- [ ] `--diff-indv-map` - Map individuals between files

**Estimated Time**: 5-6 days  
**Priority**: MEDIUM

---

## Phase 6: Advanced Analysis (10 features)

### Relatedness & Diversity (6 features)
- [ ] `--relatedness` - Relatedness coefficient
- [ ] `--relatedness2` - Alternative relatedness metric
- [ ] `--LROH` - Runs of homozygosity
- [ ] `--hapcount` - Haplotype counts
- [ ] `--indv-burden` - Individual burden
- [ ] `--indv-freq-burden` - Individual frequency burden
- [ ] `--indv-freq-burden2` - Alternative frequency burden

### PCA (3 features)
- [ ] `--pca` - Principal component analysis
- [ ] `--pca-no-norm` - PCA without normalization
- [ ] `--pca-snp-loadings` - PCA SNP loadings

### Pedigree (1 feature)
- [ ] `--mendel` - Mendelian error check

**Estimated Time**: 8-10 days  
**Priority**: LOW-MEDIUM  
**Complexity**: HIGH (requires matrix operations, eigenvalue decomposition)

---

## Phase 7: Advanced Filtering & Misc (25 features)

### Mask-Based Filtering (3 features)
- [ ] `--mask` - Mask file filtering
- [ ] `--invert-mask` - Inverted mask filtering
- [ ] `--mask-min` - Mask minimum threshold

### Advanced Allele Filters (8 features)
- [ ] `--non-ref-ac` - Non-reference allele count
- [ ] `--max-non-ref-ac` - Maximum non-ref allele count
- [ ] `--non-ref-ac-any` - Any non-ref allele count
- [ ] `--max-non-ref-ac-any` - Maximum any non-ref count
- [ ] `--non-ref-af` - Non-reference allele frequency
- [ ] `--max-non-ref-af` - Maximum non-ref frequency
- [ ] `--non-ref-af-any` - Any non-ref frequency
- [ ] `--max-non-ref-af-any` - Maximum any non-ref frequency

### Genotype-Level Filtering (2 features)
- [ ] `--remove-filtered-geno` - Remove filtered genotypes
- [ ] `--remove-filtered-geno-all` - Remove all filtered genotypes

### Additional Filtering (4 features)
- [ ] `--max-missing-count` - Maximum missing count
- [ ] `--max-indv` - Maximum number of individuals
- [ ] `--phased` - Work with phased data only
- [ ] `--derived` - Derived allele frequency

### Output Management (3 features)
- [ ] `--kept-sites` - Output list of kept sites
- [ ] `--removed-sites` - Output list of removed sites
- [ ] `--temp` - Temporary directory

### Miscellaneous (5 features)
- [ ] `--contigs` - Contig information
- [ ] `--hwe` - Alternative HWE test
- [ ] `--version` - Show version
- [ ] `--keep-INFO-all` - Keep all INFO fields

**Estimated Time**: 4-5 days  
**Priority**: LOW

---

## Phase 8: BCF Support (3 features)

- [ ] `--bcf` - Input BCF file
- [ ] `--recode-bcf` - Output BCF format
- [ ] BCF format reading/writing infrastructure

**Estimated Time**: 3-4 days  
**Priority**: LOW  
**Complexity**: MEDIUM (requires BCF library integration)

---

## Implementation Strategy

### Approach
1. Implement features in phases, from highest to lowest priority
2. Add comprehensive tests for each feature
3. Validate against original vcftools output
4. Update documentation incrementally
5. Run security checks after each phase

### Testing Strategy
- Unit tests for each new function
- Integration tests for end-to-end workflows
- Comparison tests against original vcftools output
- Edge case testing

### Documentation Updates
- Update README.md with new features
- Update FEATURE_COMPARISON.md to track progress
- Add examples for complex features
- Update ROADMAP.md as features are completed

---

## Summary

**Total Estimated Time**: 40-50 days (8-10 weeks)

**Breakdown by Priority**:
- HIGH: 42 features (Phases 1-3)
- MEDIUM: 22 features (Phases 4-5)
- LOW: 48 features (Phases 6-8)

**Breakdown by Complexity**:
- LOW: 50 features (straightforward implementation)
- MEDIUM: 40 features (moderate algorithms)
- HIGH: 22 features (complex algorithms, matrix operations)

---

## Progress Tracking

Update this section as features are implemented.

**Phase 1**: 0/20 features (0%)  
**Phase 2**: 0/10 features (0%)  
**Phase 3**: 0/12 features (0%)  
**Phase 4**: 0/12 features (0%)  
**Phase 5**: 0/10 features (0%)  
**Phase 6**: 0/10 features (0%)  
**Phase 7**: 0/25 features (0%)  
**Phase 8**: 0/3 features (0%)

**Overall**: 35/147 features (24%)

---

*Last Updated*: 2025-10-21
