# vcftools Go Implementation Roadmap

This document outlines planned features and enhancements for the vcftools Go implementation.

## Version 1.0 (Current) ✅

**Status:** Complete

### Implemented Features

- Core filtering (position, quality, allele frequency, variant type)
- Sample management (keep/remove individuals)
- Basic statistics (frequency, depth, missingness, HWE, Ts/Tv, nucleotide diversity)
- VCF recoding with filtering
- Comprehensive test suite
- Complete documentation

### Statistics

- 35 command-line options
- ~80% coverage of commonly-used features
- 0 security vulnerabilities
- 100% test pass rate

## Version 1.1 (Planned)

**Focus:** Critical missing features for population genetics

### Linkage Disequilibrium Analysis

- [ ] `--geno-r2` - Genotype-based LD (r²) calculation
- [ ] `--geno-r2-positions` - LD for specific positions vs all others
- [ ] `--ld-window` - Window size in number of SNPs
- [ ] `--ld-window-bp` - Window size in base pairs
- [ ] `--min-r2` - Minimum r² threshold for output
- [ ] `--hap-r2` - Haplotype-based LD (for phased data)
- [ ] `--hap-r2-positions` - Haplotype LD for specific positions

**Priority:** HIGH  
**Estimated Effort:** 2-3 weeks  
**Impact:** Critical for GWAS and population genetics studies

### Population Structure Statistics

- [ ] `--weir-fst-pop` - Fst calculation (Weir & Cockerham 1984)
- [ ] `--fst-window-size` - Window size for Fst calculations
- [ ] `--fst-window-step` - Step size for sliding windows

**Priority:** HIGH  
**Estimated Effort:** 1-2 weeks  
**Impact:** Essential for population differentiation studies

## Version 1.2 (Planned)

**Focus:** Format conversion and interoperability

### PLINK Format Conversion

- [ ] `--plink` - Convert to PLINK PED/MAP format
- [ ] `--plink-tped` - Convert to PLINK transposed format
- [ ] `--chrom-map` - Chromosome name to integer mapping
- [ ] Handle multi-allelic sites appropriately

**Priority:** MEDIUM-HIGH  
**Estimated Effort:** 1-2 weeks  
**Impact:** Critical for integration with PLINK analyses

### Additional Format Conversions

- [ ] `--012` - Output 0/1/2 genotype matrix
- [ ] `--BEAGLE-GL` - BEAGLE genotype likelihoods (GL)
- [ ] `--BEAGLE-PL` - BEAGLE genotype likelihoods (PL)

**Priority:** MEDIUM  
**Estimated Effort:** 1 week  
**Impact:** Useful for imputation and downstream analysis

## Version 1.3 (Planned)

**Focus:** Windowed statistics and advanced analysis

### Windowed Statistics

- [ ] `--window-pi` - Nucleotide diversity in windows
- [ ] `--window-pi-step` - Step size for pi windows
- [ ] `--TajimaD` - Tajima's D statistic
- [ ] `--SNPdensity` - SNP density in windows

**Priority:** MEDIUM  
**Estimated Effort:** 2 weeks  
**Impact:** Important for selection scans and diversity analysis

### Additional Statistics

- [ ] `--het` - Individual heterozygosity
- [ ] `--singletons` - Singleton site analysis
- [ ] `--TsTv-by-count` - Ts/Tv by allele count
- [ ] `--TsTv-by-qual` - Ts/Tv by quality score

**Priority:** MEDIUM  
**Estimated Effort:** 1 week  
**Impact:** Useful for QC and population genetics

## Version 1.4 (Planned)

**Focus:** VCF comparison and validation

### VCF Comparison

- [ ] `--diff` - Compare two VCF files
- [ ] `--gzdiff` - Compare with gzipped VCF
- [ ] `--diff-site` - Sites in common/unique between files
- [ ] `--diff-indv` - Individuals in common/unique
- [ ] `--diff-site-discordance` - Site-by-site discordance
- [ ] `--diff-indv-discordance` - Per-individual discordance
- [ ] `--diff-indv-map` - Map individual names between files

**Priority:** MEDIUM  
**Estimated Effort:** 2 weeks  
**Impact:** Important for validation and QC

## Version 2.0 (Future)

**Focus:** Advanced analysis and optimization

### Advanced Filtering

- [ ] `--bed` - Filter by BED file regions
- [ ] `--exclude-bed` - Exclude BED file regions
- [ ] `--snp` / `--snps` - Filter by SNP IDs
- [ ] `--exclude` - Exclude SNPs by ID
- [ ] `--mask` - Mask-based filtering
- [ ] `--thin` - Thin sites by distance

**Priority:** LOW-MEDIUM  
**Estimated Effort:** 1-2 weeks  
**Impact:** Useful for specific use cases

### Advanced Analysis

- [ ] `--relatedness` - Relatedness coefficient
- [ ] `--relatedness2` - Alternative relatedness metric
- [ ] `--LROH` - Runs of homozygosity

**Priority:** LOW  
**Estimated Effort:** 2-3 weeks  
**Impact:** Specialized analyses

### Performance Optimization

- [ ] Parallel processing for statistics
- [ ] Memory optimization for large VCFs
- [ ] Indexed VCF support for random access
- [ ] BCF format support

**Priority:** LOW-MEDIUM  
**Estimated Effort:** Ongoing  
**Impact:** Improved performance for large datasets

## Not Planned

These features are unlikely to be implemented due to complexity or limited use:

- `--pca` / `--pca-no-norm` / `--pca-snp-loadings` - Better done with dedicated tools
- `--mendel` - Specialized pedigree analysis
- `--ldhat` / `--ldhat-geno` / `--ldhelmet` - Format-specific for specialized tools
- `--IMPUTE` - Format-specific for IMPUTE software
- Complex haplotype analysis beyond basic LD

## Contributing

We welcome contributions! If you're interested in implementing any of these features:

1. Open an issue to discuss the implementation approach
2. Review the existing code structure and patterns
3. Implement with comprehensive tests
4. Update documentation
5. Submit a pull request

### Priority Guidelines

**HIGH Priority:**

- Features used in >50% of typical analyses
- Critical for common workflows (GWAS, population genetics)
- High user demand

**MEDIUM Priority:**

- Features used in 20-50% of analyses
- Important for specific workflows
- Moderate user demand

**LOW Priority:**

- Features used in <20% of analyses
- Specialized use cases
- Low user demand or better alternatives exist

## Version History

- **v1.0** (2025-10-21) - Initial release with core features
  - 35 command-line options
  - ~80% coverage of common use cases
  - Comprehensive test suite
  - Full documentation

## Feedback

Please open an issue on GitHub to:

- Request new features
- Report bugs
- Suggest improvements
- Contribute code

---

*This roadmap is subject to change based on user feedback and community contributions.*
