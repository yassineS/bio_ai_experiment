# vcftools Feature-Parity Plan (aspirational)

> **This is a plan, not a status report.** Despite the original title ("100%
> Feature Parity"), full parity has **not** been achieved or even closely
> approached. As of 2026-05-12 the Go port implements roughly 40 of vcftools'
> ~147 options. A few of the "Phase 1" statistics below were wired up
> (`--site-pi` fixed, `--window-pi`/`--window-pi-step`, `--TsTv-by-count`,
> `--depth`), and the ones that were never implemented (`--TsTv-by-qual`,
> `--hist-indel-len`, `--geno-depth`, `--TajimaD`, `--weir-fst-pop`,
> `--fst-window-*`) now return an error instead of being silently ignored.
> For the actual current state see
> [../PORTING_STATUS.md](../PORTING_STATUS.md) and [README.md](README.md).
> The phase breakdown below is kept as a roadmap.

**Original date**: 2025-10-21 · **Reality check**: 2026-05-12

---

## Executive Summary (planning estimate, not progress)

Reaching full feature parity with original vcftools would require roughly
**8-10 weeks** of dedicated development — adding ~110 missing features across
the phases below. None of that timeline has been executed; the items marked
"✅" in the phase notes refer to CLI flags being added, not necessarily to the
underlying logic working.

### Scope Assessment

**Total Features**: 147  
**Currently Implemented**: 35 (24%)  
**Remaining**: 112 (76%)

**Lines of Code Estimate**: 15,000-20,000 additional lines  
**Test Coverage Target**: Maintain >50% coverage  
**Documentation**: Comprehensive updates required for all features

---

## Implementation Phases

### Phase 1: Additional Filtering & Basic Statistics (20 features)

**Status**: In Progress  
**Duration**: 3-4 days  
**Started**: 2025-10-21  

#### New Statistics Flags Added

- ✅ `--freq2` - Alternative frequency output format
- ✅ `--counts2` - Alternative counts output format  
- ✅ `--TsTv-by-count` - Ts/Tv by allele count
- ✅ `--TsTv-by-qual` - Ts/Tv by quality score
- ✅ `--het` - Heterozygosity statistics
- ✅ `--singletons` - Singleton site analysis
- ✅ `--hist-indel-len` - Indel length histogram
- ✅ `--geno-depth` - Genotype depth distribution

#### Implementation Status

- [x] CLI flags added to main.go
- [x] Parameters added to Params struct
- [ ] Statistics collection logic
- [ ] Output formatting
- [ ] Unit tests
- [ ] Integration tests
- [ ] Documentation updates

#### Remaining Phase 1 Features

- [ ] BED file filtering (--bed, --exclude-bed)
- [ ] Position overlap filtering
- [ ] SNP ID filtering (--snp, --snps, --exclude, --thin)
- [ ] Genotype-level filters (--minDP, --maxDP, --minGQ)
- [ ] FILTER flag operations (--keep-filtered, --remove-filtered)
- [ ] INFO flag operations (--keep-INFO, --remove-INFO)

---

### Phase 2: Population Genetics Statistics (10 features)

**Status**: Planned  
**Duration**: 5-6 days  
**Dependencies**: Phase 1 complete

#### Features

- Windowed nucleotide diversity (--window-pi, --window-pi-step)
- Tajima's D statistic (--TajimaD)
- SNP density in windows (--SNPdensity)
- Fst calculations (--weir-fst-pop, --fst-window-size, --fst-window-step)
- Additional Ts/Tv metrics (--TsTv-by-count, --TsTv-by-qual)
- FILTER summary (--FILTER-summary)

#### Technical Requirements

- Window-based iteration over variants
- Population structure algorithms
- Fst calculation (Weir & Cockerham 1984)
- Tajima's D calculation
- Statistical summaries

---

### Phase 3: Linkage Disequilibrium Analysis (12 features)

**Status**: Planned  
**Duration**: 7-8 days  
**Complexity**: HIGH

#### Features

- Genotype-based LD (--geno-r2)
- Haplotype-based LD (--hap-r2)
- Chi-square test (--geno-chisq)
- Position-specific LD (--geno-r2-positions, --hap-r2-positions)
- Inter-chromosomal LD (--interchrom-geno-r2, --interchrom-hap-r2)
- LD window options (--ld-window, --ld-window-bp, --ld-window-min, --ld-window-bp-min, --min-r2)

#### Technical Requirements

- r² calculation for genotypes
- D and D' statistics for haplotypes
- Interval tree data structure for efficient window queries
- Chi-square independence tests
- Memory-efficient pairwise comparisons
- Phased genotype handling

#### Algorithm Complexity

- Pairwise comparisons: O(n²) for n variants
- With windowing: O(n × w) where w is window size
- Optimization strategies needed for large datasets

---

### Phase 4: Format Conversions (12 features)

**Status**: Planned  
**Duration**: 6-7 days

#### Features

- PLINK format (--plink, --plink-tped, --chrom-map)
- 012 matrix (--012)
- BEAGLE format (--BEAGLE-GL, --BEAGLE-PL)
- IMPUTE format (--IMPUTE)
- LDhat formats (--ldhat, --ldhat-geno, --ldhelmet)
- INFO/FORMAT extraction (--get-INFO, --extract-FORMAT-info, --recode-INFO)

#### Technical Requirements

- PLINK PED/MAP format writing
- PLINK TPED/TFAM format writing
- Genotype encoding (0/1/2 for PLINK)
- Genotype likelihood formatting
- Chromosome name to integer mapping
- Multi-allelic site handling

---

### Phase 5: VCF Comparison/Diff (10 features)

**Status**: Planned  
**Duration**: 5-6 days

#### Features

- Diff input (--diff, --gzdiff, --diff-bcf)
- Site comparison (--diff-site)
- Individual comparison (--diff-indv)
- Discordance analysis (--diff-site-discordance, --diff-indv-discordance)
- Discordance matrix (--diff-discordance-matrix)
- Switch error detection (--diff-switch-error)
- Individual mapping (--diff-indv-map)

#### Technical Requirements

- Synchronized iteration over two VCF files
- Site matching algorithms
- Genotype concordance calculations
- Phasing error detection
- Individual name mapping

---

### Phase 6: Advanced Analysis (10 features)

**Status**: Planned  
**Duration**: 8-10 days  
**Complexity**: HIGH

#### Features

- Relatedness (--relatedness, --relatedness2)
- Runs of homozygosity (--LROH)
- Haplotype counts (--hapcount)
- Individual burden (--indv-burden, --indv-freq-burden, --indv-freq-burden2)
- PCA (--pca, --pca-no-norm, --pca-snp-loadings)
- Mendelian errors (--mendel)

#### Technical Requirements

- Kinship coefficient calculations
- ROH detection algorithms
- PCA with eigenvalue decomposition
- Matrix operations (covariance, correlation)
- Pedigree file parsing
- Mendelian inheritance checks

#### External Dependencies

- May require matrix algebra library (gonum/mat)

---

### Phase 7: Advanced Filtering & Misc (25 features)

**Status**: Planned  
**Duration**: 4-5 days

#### Features

- Mask-based filtering (--mask, --invert-mask, --mask-min)
- Non-ref allele filters (8 options)
- Genotype-level filtering (--remove-filtered-geno, --remove-filtered-geno-all)
- Additional filtering (--max-missing-count, --max-indv, --phased, --derived)
- Output management (--kept-sites, --removed-sites, --temp)
- Miscellaneous (--contigs, --hwe, --version, --keep-INFO-all)

---

### Phase 8: BCF Support (3 features)

**Status**: Planned  
**Duration**: 3-4 days

#### Features

- BCF input (--bcf)
- BCF output (--recode-bcf)
- BCF format infrastructure

#### Technical Requirements

- BCF format parsing (binary VCF)
- BGZF compression handling
- Binary field encoding/decoding
- May require external library (htslib bindings or pure Go implementation)

---

## Technical Challenges

### High Complexity Features

1. **Linkage Disequilibrium**
   - Complex statistical calculations
   - O(n²) algorithm complexity
   - Memory management for large datasets
   - Requires interval trees for efficiency

2. **PCA Analysis**
   - Matrix operations
   - Eigenvalue decomposition
   - Requires robust linear algebra library
   - Normalization strategies

3. **BCF Format**
   - Binary format parsing
   - BGZF compression
   - Consider using external library vs pure Go implementation

4. **Fst Calculations**
   - Population structure algorithms
   - Weir & Cockerham 1984 method
   - Window-based calculations
   - Statistical variance estimates

### Testing Strategy

For each phase:

1. Unit tests for individual functions
2. Integration tests for end-to-end workflows
3. Comparison tests against original vcftools output
4. Edge case testing
5. Performance benchmarking

### Documentation Requirements

For each phase:

1. Update README.md with new features
2. Add usage examples
3. Update FEATURE_COMPARISON.md
4. Document algorithm details for complex features
5. Add migration notes for any behavioral differences

---

## Resource Estimates

### Development Time

- **Phase 1**: 3-4 days (20 features)
- **Phase 2**: 5-6 days (10 features)
- **Phase 3**: 7-8 days (12 features)
- **Phase 4**: 6-7 days (12 features)
- **Phase 5**: 5-6 days (10 features)
- **Phase 6**: 8-10 days (10 features)
- **Phase 7**: 4-5 days (25 features)
- **Phase 8**: 3-4 days (3 features)

**Total**: 41-50 days (8-10 weeks)

### Code Volume

- **New Go Code**: ~15,000-20,000 lines
- **Test Code**: ~8,000-10,000 lines
- **Documentation**: ~5,000-8,000 words

### External Dependencies

Potential needs:

- `gonum/mat` - Matrix operations for PCA
- `htslib` bindings or pure Go BCF parser
- Statistical libraries for advanced calculations

---

## Risk Assessment

### High Risk Items

1. **LD Analysis** - Complex algorithms, performance concerns
2. **PCA** - Requires robust matrix operations
3. **BCF Support** - Binary format complexity
4. **Fst Calculations** - Statistical algorithm complexity

### Mitigation Strategies

1. Incremental implementation with validation at each step
2. Use established algorithms from literature
3. Extensive testing against original vcftools
4. Performance benchmarking and optimization
5. Consider external libraries for complex operations

---

## Next Steps

### Immediate (This Week)

1. Complete Phase 1 statistics implementation
2. Add unit tests for new statistics
3. Add integration tests
4. Update documentation

### Short Term (Next 2 Weeks)

1. Complete Phase 1
2. Begin Phase 2 (Population genetics)
3. Implement windowed statistics
4. Implement Fst calculations

### Medium Term (Weeks 3-6)

1. Complete Phases 2-4
2. Implement LD analysis
3. Implement format conversions
4. Comprehensive testing

### Long Term (Weeks 7-10)

1. Complete Phases 5-8
2. VCF comparison
3. Advanced analysis (PCA, relatedness)
4. BCF support
5. Final testing and optimization

---

## Current Session Progress

### Completed Today

- ✅ Created IMPLEMENTATION_PLAN.md
- ✅ Added 8 new statistics flags to CLI
- ✅ Updated Params struct
- ✅ Created this progress report
- ✅ Committed Phase 1 initialization

### Next Actions

1. Implement statistics collection for new features
2. Add output formatters
3. Write unit tests
4. Write integration tests
5. Update documentation

---

## Recommendations

Given the 8-10 week scope, consider:

1. **Incremental Review**: Review and merge each phase separately
2. **Priority Reassessment**: Focus on most-used features first within each phase
3. **External Libraries**: Evaluate using established libraries for complex operations
4. **Performance**: Profile and optimize memory-intensive operations
5. **Validation**: Extensive comparison with original vcftools outputs

---

*Last Updated*: 2025-10-21T08:34:54Z  
*Next Update*: After Phase 1 completion
