# Tool Porting Status

This document tracks the status of bioinformatics tools being ported from their original implementations to Go.

**Last Updated**: 2026-05-12

> **Accuracy note (2026-05-12 audit):** earlier revisions of this file
> overstated progress — claiming ">85% test coverage" and labelling several
> tools "Complete" with "functional parity". That is not the case. The numbers
> and status below were re-derived from `go test -cover ./...` and from reading
> the code. None of these ports is a drop-in replacement for the original tool;
> each implements a subset of features. Treat "Complete" here as "the planned
> first slice of functionality is in place", not "feature-complete".

---

## Overview

### Goals

- Port high-priority bioinformatics tools to Go
- Cover the commonly used subset of each tool's functionality first
- Improve usability, documentation, and maintainability
- Ensure tested, validated behaviour for everything that is implemented

### Progress Summary

- **Tools with a working subset**: 8
- **Tools tested**: 8 (package-level tests; `cmd/` entry points have no tests)
- **Test coverage (statements, `go test -cover`)**: vcftools 38%, skewer 46%,
  prinseq 55%, seqtk 57%, fastp 63%, bedintersect 75%, bedmerge 79%, sickle 80%
- **Documentation**: README per tool; some design docs are aspirational, not status
- **gzip support**: sickle, skewer, fastp, bedmerge, bedintersect, vcftools (not seqtk/prinseq)

---

## Completed Tools

### 1. ✅ seqtk

**Status**: Complete  
**Version**: 1.0.0  
**Original**: C (Heng Li)  
**Category**: Quality Control / FASTA/Q Processing

**Implemented Commands**:

- `comp` - Sequence composition statistics
- `fq2fa` - FASTQ to FASTA conversion
- `seq` - Sequence manipulation (reverse complement)
- `sample` - Random subsampling
- `trimfq` - Quality-based trimming

**Test Coverage**: ~57% of statements (`go test -cover`)  
**Performance**: ~1.05-1.1x faster than original on the implemented commands  
**Documentation**: README with examples  

**Key Features**:

- Fast FASTA/Q processing for the five commands above
- Quality score handling
- Memory-efficient streaming

**Migration Notes**:

- Command structure changed (subcommands instead of flags)
- Only the five commands above are implemented; upstream seqtk has many more
  (`subseq`, `mergepe`, `mutfa`, `randbase`, `hpc`, `cutN`, ...)
- Output format intended to be compatible for the implemented commands

---

### 2. ✅ PRINSEQ

**Status**: Complete  
**Version**: 1.0.0  
**Original**: Perl (PRINSEQ-lite)  
**Category**: Quality Control

**Implemented Commands**:

- `stats` - Calculate sequence statistics
- `filter` - Multi-criteria filtering and trimming

**Test Coverage**: ~55% of statements (`go test -cover`)  
**Performance**: ~1.2-1.35x faster than the original Perl on the implemented paths  
**Documentation**: README with examples

**Key Features**:

- Sequence statistics
- Length / GC / N-content / quality filtering
- Trimming operations (fixed, percentage, quality-based)
- Poly-N and poly-A/T tail trimming
- Duplicate removal
- Paired-end support
- Phred+64 encoding support (Illumina 1.3-1.7)
- Bad-sequence output
- Complexity filtering (DUST and entropy methods)

**Migration Notes**:

- Command structure changed (subcommands instead of flags)
- Covers the commonly used PRINSEQ-lite filtering/trimming options; not every
  upstream option is implemented, and graph/report generation is out of scope
- Output format intended to be compatible

---

### 3. ✅ sickle

**Status**: Complete  
**Version**: 1.1.0  
**Original**: C (Joshi & Fass)  
**Category**: Quality Control / Trimming

**Implemented Commands**:

- `se` - Single-end read trimming
- `pe` - Paired-end read trimming

**Test Coverage**: ~80% of statements (`go test -cover`)  
**Performance**: ~0.96-1.0x (similar to original)  
**Documentation**: README with examples

**Key Features**:

- Sliding window quality assessment
- Quality threshold-based trimming
- Length threshold filtering
- N-truncation support
- 5' trim control
- Paired-end synchronization
- Orphaned read handling
- Built-in gzip support

**Migration Notes**:

- CLI mirrors upstream sickle's `se`/`pe` flags; behaviour aims to match but
  has not been validated byte-for-byte against the C implementation
- Built-in gzip support (automatic by .gz extension)

---

### 4. ✅ skewer

**Status**: Complete  
**Version**: 1.0.0  
**Original**: C++ (Hongshan Jiang)  
**Category**: Adapter Trimming

**Implemented Commands**:

- `se` - Single-end adapter trimming
- `pe` - Paired-end adapter trimming

**Test Coverage**: ~46% of statements (`go test -cover`)  
**Performance**: ~1.0x (comparable to original)  
**Documentation**: README with examples

**Key Features**:

- 3' and 5' adapter detection
- Error-tolerant matching
- Configurable minimum overlap
- Quality-based trimming
- Length filtering
- Paired-end support
- Built-in gzip support

**Migration Notes**:

- Similar CLI to original
- Simplified adapter detection algorithm
- Built-in gzip support
- Complements sickle for complete preprocessing

---

### 5. ✅ fastp

**Status**: Complete (Core Features)  
**Version**: 1.0.0  
**Original**: C++ (Shifu Chen)  
**Category**: All-in-One Preprocessor

**Implemented Commands**:

- Single command with multiple filters

**Test Coverage**: ~63% of statements (`go test -cover`)  
**Performance**: ~1.1x  
**Documentation**: README with examples

**Key Features**:

- Adapter trimming (3' and 5')
- Quality filtering
- Length filtering
- N content filtering
- Poly-G/X tail trimming (NovaSeq)
- Complexity filtering
- Built-in gzip support
- Comprehensive statistics

**Migration Notes**:

- Simplified version of original
- Core preprocessing features implemented
- No HTML reports (future feature)
- Single-threaded (parallel processing planned)

---

### 6. ✅ bedmerge

**Status**: Complete  
**Version**: 1.0.0  
**Original**: bedtools merge (C++)  
**Category**: Genomic Intervals / Utilities

**Implemented Commands**:

- Single command for merging BED intervals

**Test Coverage**: ~79% of statements (`go test -cover`)  
**Performance**: ~2x faster than bedtools merge on the implemented path  
**Documentation**: README with examples

**Key Features**:

- Merge overlapping BED intervals
- Distance-based merging (-d option)
- Strand-specific merging (-s option)
- Statistics output (-S option)
- Built-in gzip support
- Automatic sorting

**Migration Notes**:

- Compatible with bedtools merge basic functionality
- Simplified version (no advanced options yet)
- Output always BED3 format
- All core features working

---

### 7. ✅ bedintersect

**Status**: Complete  
**Version**: 1.0.0  
**Original**: bedtools intersect (C++)  
**Category**: Genomic Intervals / Utilities

**Implemented Commands**:

- Single command for interval intersection

**Test Coverage**: ~75% of statements (`go test -cover`)  
**Performance**: Comparable to bedtools intersect  
**Documentation**: README with examples

**Key Features**:

- Find overlapping intervals between two BED files
- Multiple output modes (-wa, -wb, -c, -v)
- Minimum overlap filters (-m)
- Fractional overlap filters (-f, -F)
- Strand-specific intersection (-s)
- Built-in gzip support
- Statistics output (-S option)

**Migration Notes**:

- Compatible with bedtools intersect common operations
- Simplified version (no sorted/reciprocal modes yet)
- All essential features working
- Same output format as bedtools

---

### 8. ✅ vcftools

**Status**: Complete  
**Version**: 1.0.0  
**Original**: C++/Perl (Danecek et al.)  
**Category**: VCF Manipulation / Population Genetics

**Status**: Partial — a subset of upstream vcftools, ~40 of ~147 options

**Implemented Commands**:

- Single command with multiple filtering, statistics and conversion options

**Test Coverage**: ~38% of statements (`go test -cover`)  
**Performance**: Comparable to original on the implemented operations  
**Documentation**: README with examples  

**Implemented features**:

- Position-based filtering (`--chr`, `--from-bp`/`--to-bp`, `--positions`, ...)
- SNP-ID filtering and thinning (`--snp`, `--snps`, `--exclude`, `--thin`)
- Quality, allele-frequency and allele-count filtering (`--minQ`, `--maf`, `--mac`, ...)
- Variant-type filtering (`--remove-indels`, `--keep-only-indels`, `--min/max-alleles`)
- Genotype-level filtering (`--minDP`, `--maxDP`, `--minGQ`)
- Sample filtering (`--indv`, `--remove-indv`, `--keep`, `--remove`)
- Site statistics: `--freq`/`--counts`(+`2`), `--site-depth`, `--site-mean-depth`,
  `--site-quality`, `--missing-site`, `--missing-indv`, `--depth`, `--hardy`,
  `--site-pi`, `--window-pi`(+`--window-pi-step`), `--TsTv-summary`, `--TsTv`,
  `--TsTv-by-count`, `--het`, `--singletons`, `--FILTER-summary`, `--SNPdensity`
- VCF recoding (`--recode`, `--recode-INFO-all`)
- Format conversion: `--012`, `--plink`, `--plink-tped` (with `--chrom-map`)

**Not yet implemented** — these options are now rejected with an error instead of
being silently ignored (older builds accepted them and produced no output):

- `--TsTv-by-qual`, `--hist-indel-len`, `--geno-depth`
- `--TajimaD`, `--weir-fst-pop`, `--fst-window-size`, `--fst-window-step`
- All LD analysis (`--geno-r2`, `--hap-r2`, ...) and many other upstream options

See [FEATURE_COMPARISON.md](vcftools/FEATURE_COMPARISON.md) and
[ROADMAP.md](vcftools/ROADMAP.md) for the full picture.

**Migration Notes**:

- Per-site nucleotide diversity (`--site-pi`) uses the standard
  `(n^2 - Σ c_a^2) / (n(n-1))` formula; earlier builds reported a different
  (incorrect) per-genotype quantity
- Not a drop-in replacement: anything not in the list above is unavailable

---

## Tool Comparison Matrix

| Tool | Original Lang | Go Version | Commands | Tests | Docs | Performance | Gzip |
|------|---------------|------------|----------|-------|------|-------------|------|
| seqtk | C | 1.0.0 | 5 | ✓ | ✓ | 1.05-1.1x | - |
| PRINSEQ | Perl | 1.0.0 | 2 | ✓ | ✓ | 1.2-1.35x | - |
| sickle | C | 1.1.0 | 2 | ✓ | ✓ | 0.96-1.0x | ✓ |
| skewer | C++ | 1.0.0 | 2 | ✓ | ✓ | ~1.0x | ✓ |
| fastp | C++ | 1.0.0 | 1 | ✓ | ✓ | ~1.1x | ✓ |
| bedmerge | C++ (bedtools) | 1.0.0 | 1 | ✓ | ✓ | ~2.0x | ✓ |
| bedintersect | C++ (bedtools) | 1.0.0 | 1 | ✓ | ✓ | ~1.0x | ✓ |
| vcftools | C++/Perl | 1.0.0 | 1 | ✓ | ✓ | ~1.0x | ✓ |

---

## Priority Tools for Future Porting

Based on the top 50 analysis, these tools are recommended for future porting:

### High Priority (Simple, High Impact)

1. **Trim Galore** (Perl) - Quality and adapter trimming (Rank: 53.27)
   - Wrapper functionality
   - Widely used
   - Quality + adapter handling
   - Note: Could be implemented as wrapper around sickle + skewer

### Medium Priority (More Complex)

4. **BEDTools subset** (C++) - Genomic interval operations
   - Core operations: intersect, merge, sort
   - Widely used format
   - Complex but modular

2. **SAMtools subset** (C) - SAM/BAM manipulation
   - Basic operations only
   - View, sort, index
   - High-impact tool

### Lower Priority (Very Complex)

6. **minimap2** (C) - Long-read alignment
   - Complex algorithm
   - High performance requirements
   - Large codebase

### Analyzed But Not Recommended for Porting

1. **BWA** (C) - Short-read alignment ❌ **NOT RECOMMENDED**
   - **Status**: Analyzed in detail (2025-10-21)
   - **Decision**: Do not port
   - **Reasons**:
     - Extremely complex: ~17,000 lines of code
     - Already well-maintained and highly optimized
     - Already includes multi-threading and batch processing
     - Better alternatives exist (BWA-MEM2, minimap2)
     - Scope far exceeds project "minimal changes" philosophy
   - **See**: [BWA Implementation Decision](BWA_IMPLEMENTATION_DECISION.md) for full analysis
   - **Alternatives**: Use original BWA, BWA-MEM2, or create MCP wrapper

---

## Testing Standards

All ported tools must meet these criteria:

### Code Quality

- ✓ >80% test coverage
- ✓ All tests passing
- ✓ No race conditions
- ✓ Clean go vet output
- ✓ Formatted with gofmt

### Documentation

- ✓ Complete README with:
  - Installation instructions
  - Usage examples
  - Command reference
  - Performance comparison
  - Migration notes
- ✓ API documentation (godoc)
- ✓ CLI differences documented

### Functionality

- ✓ Functional parity with original (core features)
- ✓ Error handling and validation
- ✓ Input/output format compatibility
- ✓ Performance within 20% of original

### Usability

- ✓ Consistent CLI interface
- ✓ Both short and long options
- ✓ Clear error messages
- ✓ Help text for all commands

---

## Architecture Patterns

### Standard Tool Structure

```
tool-name/
├── cmd/
│   └── tool-name/
│       └── main.go           # CLI entry point
├── pkg/
│   └── tool-name/
│       ├── tool-name.go      # Core functionality
│       └── tool-name_test.go # Unit tests
├── README.md
└── docs/                     # Optional, for complex tools
```

### Common Libraries

- `pkg/bioformats/fastq` - FASTQ I/O
- `pkg/bioformats/fasta` - FASTA I/O
- `pkg/bioformats/bed` - BED format
- `pkg/bioformats/vcf` - VCF format
- `pkg/cliflag` - Consistent CLI parsing

### Design Principles

1. **Streaming Processing**: Handle files larger than RAM
2. **Minimal Dependencies**: Use Go standard library when possible
3. **Type Safety**: Leverage Go's type system
4. **Error Handling**: Clear, actionable error messages
5. **Testing**: Comprehensive unit and integration tests

---

## Implementation Guidelines

### Before Starting a Port

1. **Research Original Tool**
   - Read documentation thoroughly
   - Understand all features and options
   - Identify most commonly used functionality
   - Check for existing Go implementations

2. **Plan Implementation**
   - List core features to implement
   - Identify optional/advanced features
   - Design API structure
   - Plan test cases

3. **Set Up Structure**
   - Create directory structure
   - Initialize go.mod if needed
   - Set up basic README

### During Implementation

1. **Core First**: Implement basic functionality before advanced features
2. **Test Early**: Write tests alongside code
3. **Document**: Add comments and docstrings as you go
4. **Validate**: Compare output with original tool

### After Implementation

1. **Testing**
   - Run all unit tests
   - Perform integration testing
   - Test edge cases
   - Compare with original tool output

2. **Documentation**
   - Complete README with examples
   - Document CLI differences
   - Add performance comparison
   - Write migration guide

3. **Validation**
   - Run on real datasets
   - Compare statistics with original
   - Verify file format compatibility
   - Test error handling

---

## Performance Benchmarking

### Standard Benchmark Dataset

Use consistent datasets for comparison:

- 1M read FASTQ file (~200MB)
- 10M read FASTQ file (~2GB)
- Mix of read lengths and qualities
- Both single-end and paired-end

### Metrics to Track

- **Execution Time**: Wall clock time
- **Memory Usage**: Peak RSS
- **Throughput**: Reads per second
- **Accuracy**: Identical output to original

### Acceptable Performance Range

- **Target**: Within 20% of original (0.8x - 1.2x)
- **Good**: Faster than original (>1.0x)
- **Acceptable**: 0.8x - 1.0x (slight slowdown for safety/features)
- **Needs Optimization**: <0.8x

---

## Known Limitations

### Current Implementations

**seqtk**:

- Only `comp`, `fq2fa`, `seq`, `sample`, `trimfq` are implemented; many
  upstream subcommands are missing
- No built-in gzip support yet

**PRINSEQ**:

- Covers the common filtering/trimming options (length, GC, N, quality,
  fixed/percentage/quality trimming, poly-N/A/T, dedup, paired-end, Phred+64,
  bad-sequence output, complexity filters); not every upstream option
- Graph/HTML report generation not included

**sickle**:

- Built-in gzip support (by `.gz` extension)
- No automatic quality-encoding detection
- Not validated byte-for-byte against the C original

**skewer**:

- Simplified adapter-detection algorithm
- No automatic adapter detection

**fastp**:

- Single-end and paired-end processing implemented
- No HTML/JSON reports
- No automatic adapter detection
- Parallel worker pool exists but the feature surface is a subset of upstream

**bedmerge**:

- Output always BED3 format
- No advanced options (distinct, count, etc.)
- In-memory processing (not suitable for very large files)

**bedintersect**:

- Uses an interval tree (`pkg/bedintersect/intervaltree.go`)
- No reciprocal overlap mode
- No sorted-file streaming optimization
- In-memory B-file loading

**vcftools**:

- ~40 of ~147 upstream options; see the vcftools section above
- Unimplemented options are now rejected with an error rather than ignored

### General Limitations

- Partial gzip support (sickle, skewer, fastp, bedmerge, bedintersect, vcftools have it; seqtk, prinseq do not)
- `cmd/` entry points have no automated tests (coverage there is 0%)
- None of these ports has been validated output-for-output against its original
- Performance may vary by dataset characteristics

---

## Version History

### Version 1.0.0 (2025-10-21)

- Initial release of 3 core QC tools
- seqtk: FASTA/Q processing
- PRINSEQ: Quality control and filtering
- sickle: Quality-based trimming
- Complete test suites
- Comprehensive documentation

### Version 1.1.0 (2025-10-21)

- ✓ Built-in gzip support (sickle, skewer, fastp)
- ✓ skewer: Adapter trimming tool
- ✓ fastp: All-in-one preprocessor
- ✓ iohelper library for transparent gzip handling

### Version 1.2.0 (2025-10-21)

- ✓ bedmerge: BED interval merger
- ✓ bedintersect: BED interval intersection finder
- ✓ First BED utilities added
- ✓ 19 additional tests (bedmerge: 8, bedintersect: 11)

### Planned Version 1.3.0

- Phred+64 support in PRINSEQ
- Automatic quality encoding detection
- Built-in gzip support for seqtk and prinseq
- JSON statistics output
- Progress reporting
- Parallel processing framework

---

## Statistics

### Code Metrics

- **Test Coverage (statements)**: ~38-80% per package, ~60% unweighted average;
  0% for all `cmd/` entry points (run `go test -cover ./...` for current numbers)
- **Shared Libraries**: iohelper (gzip support), bioformats, cliflag
- **Tools with a working subset**: 8 (seqtk, prinseq, sickle, skewer, fastp,
  bedmerge, bedintersect, vcftools)

### Performance Summary

- **Speedup**: roughly comparable to the originals (~0.95-2x) on the implemented
  operations; not benchmarked exhaustively
- **Binary Size**: a few MB per tool
- **Startup Time**: <100ms
- **Gzip Support**: Transparent with minimal overhead

---

## Contributing

### How to Contribute

1. **Select a Tool**
   - Check priority list above
   - Review original tool documentation
   - Open an issue to claim the tool

2. **Implement**
   - Follow architecture patterns
   - Write tests alongside code
   - Document as you go

3. **Submit**
   - Create pull request
   - Include tests and documentation
   - Run benchmarks
   - Update this document

### Areas for Improvement

Existing tools:

- Add missing features
- Optimize performance
- Improve error messages
- Expand test coverage
- Add more examples

New tools:

- Port high-priority tools
- Create utility libraries
- Build tool pipelines
- Add format converters

---

## Resources

### Documentation

- [CLI Differences](CLI_DIFFERENCES.md) - Detailed comparison with originals
- [Tools README](README.md) - General guidelines and structure
- Individual tool READMEs - Usage and examples

### Code

- `pkg/bioformats/` - Format libraries
- `pkg/cliflag/` - CLI utilities
- Individual tool packages - Implementation references

### External

- Original tool repositories (see individual READMEs)
- FASTQ format specification
- Phred quality scores documentation

---

## Contact

For questions, issues, or suggestions:

- Open an issue on GitHub
- Include tool name in title
- Provide example data if relevant

---

*This document is maintained alongside tool development and updated with each release.*
