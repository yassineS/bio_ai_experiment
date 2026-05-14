# vcftools Go Implementation Roadmap

This document outlines planned features and enhancements for the vcftools Go implementation.

## Current state

**Status:** Partial — roughly 45 of vcftools' ~147 options.

### Implemented

- Filtering: position, BED-interval (`--bed`/`--exclude-bed`), SNP ID + thinning, quality, allele frequency/count, variant type, genotype-level (`--minDP`/`--maxDP`/`--minGQ`)
- Sample management (`--indv`, `--remove-indv`, `--keep`, `--remove`)
- Per-site / per-individual statistics: `--freq`/`--counts`(+`2`), `--depth`,
  `--geno-depth`, `--site-depth`, `--site-mean-depth`, `--site-quality`,
  `--missing-site`, `--missing-indv`, `--hardy`, `--het`, `--singletons`,
  `--site-pi`, `--window-pi`(+`--window-pi-step`), `--TajimaD`,
  `--TsTv-summary`, `--TsTv`, `--TsTv-by-count`, `--TsTv-by-qual`
  (→ `.TsTv.qual`: Ts/Tv counts and ratios cumulative below and at-or-above
  each distinct QUAL threshold), `--hist-indel-len`,
  `--FILTER-summary`, `--SNPdensity`
- Population structure: `--weir-fst-pop` (Weir & Cockerham 1984 per-site Fst → `.weir.fst`),
  `--fst-window-size`/`--fst-window-step` (windowed `WEIGHTED_FST`/`MEAN_FST` → `.windowed.weir.fst`)
- VCF recoding (`--recode`, `--recode-INFO-all`)
- Format conversion: `--012`, `--plink`, `--plink-tped`, `--chrom-map`,
  `--BEAGLE-GL`, `--BEAGLE-PL`
- VCF comparison: `--diff` with `--diff-site`, `--diff-indv`,
  `--diff-site-discordance`, `--diff-indv-discordance`

### Recognised but **not implemented** (return an error)

- Inter-chromosomal LD (`--interchrom-geno-r2`, `--interchrom-hap-r2`),
  `--geno-chisq`, and many other upstream options

Test coverage of the `vcftools` package is ~51% of statements; the `cmd/`
entry point has no tests.

## Version 1.1 (Planned)

**Focus:** Critical missing features for population genetics

### Linkage Disequilibrium Analysis

- [x] `--geno-r2` - Genotype-based LD (r²) calculation
- [x] `--geno-r2-positions` - LD for specific positions vs all others
- [x] `--ld-window` - Window size in number of SNPs
- [x] `--ld-window-bp` - Window size in base pairs
- [x] `--ld-window-min` - Minimum SNP distance between LD pairs
- [x] `--ld-window-bp-min` - Minimum bp distance between LD pairs
- [x] `--min-r2` - Minimum r² threshold for output
- [x] `--hap-r2` - Haplotype-based LD (for phased data)
- [x] `--hap-r2-positions` - Haplotype LD for specific positions

*Implemented in this PR* (vcftools: linkage disequilibrium analysis). Writes
`<prefix>.geno.ld` (`CHR POS1 POS2 N_INDV R^2`) and `<prefix>.hap.ld`
(`CHR POS1 POS2 N_CHR R^2 D Dprime`). Pairs are restricted to the same
chromosome; multi-allelic sites use only the first ALT; `--hap-r2` requires
phased GTs. Upstream byte-for-byte parity hasn't been validated yet — see the
follow-up issue.

Still missing: `--interchrom-geno-r2`, `--interchrom-hap-r2`, `--geno-chisq`.

**Priority:** HIGH  
**Estimated Effort:** 2-3 weeks  
**Impact:** Critical for GWAS and population genetics studies

### Population Structure Statistics

- [x] `--weir-fst-pop` - Fst calculation (Weir & Cockerham 1984)
- [x] `--fst-window-size` - Window size for Fst calculations
- [x] `--fst-window-step` - Step size for sliding windows

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
- [x] `--BEAGLE-GL` - BEAGLE genotype likelihoods (GL)
- [x] `--BEAGLE-PL` - BEAGLE genotype likelihoods (PL)

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
- [x] `--TsTv-by-count` - Ts/Tv by allele count
- [x] `--TsTv-by-qual` - Ts/Tv by quality score (→ `.TsTv.qual`)

**Priority:** MEDIUM  
**Estimated Effort:** 1 week  
**Impact:** Useful for QC and population genetics

## Version 1.4 (Planned)

**Focus:** VCF comparison and validation

### VCF Comparison

- [x] `--diff` - Compare two VCF files
- [ ] `--gzdiff` - Compare with gzipped VCF (compose with `iohelper`; the
      same `--diff FILE` accepts `.gz` already)
- [x] `--diff-site` - Sites in common/unique between files (`.diff.sites_in_files`)
- [x] `--diff-indv` - Individuals in common/unique (`.diff.indv_in_files`)
- [x] `--diff-site-discordance` - Site-by-site discordance (`.diff.sites`)
- [x] `--diff-indv-discordance` - Per-individual discordance (`.diff.indv`)
- [ ] `--diff-indv-map` - Map individual names between files
- [ ] `--diff-discordance-matrix` - Full per-pair discordance matrices

**Priority:** MEDIUM  
**Estimated Effort:** 2 weeks  
**Impact:** Important for validation and QC

## Version 2.0 (Future)

**Focus:** Advanced analysis and optimization

### Advanced Filtering

- [x] `--bed` - Filter by BED file regions
- [x] `--exclude-bed` - Exclude BED file regions
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
