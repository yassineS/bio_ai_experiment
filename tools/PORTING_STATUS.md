# Tool Porting Status

This document tracks the status of bioinformatics tools being ported from their original implementations to Go.

**Last Updated**: 2026-05-13

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

- **Tools with a working subset**: 16 (8 original + 8 bedtools subcommands)
- **Tools tested**: 16 (package-level tests; `cmd/` entry points have no tests)
- **Test coverage (statements, `go test -cover`)** — main tools: vcftools 58%,
  seqtk 66%, fastp 67%, bedintersect 75%, sickle 82%, prinseq 99.9%,
  skewer 100%, bedmerge 100%
- **Test coverage** — new bedtools tools: bedsort 92%, bedflank 92%,
  bedclosest 93%, bedsubtract 94%, bedgenomecov 94%, bedcomplement 95%,
  bedslop 95%, bedjaccard 96%
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
- `subseq` - Extract subsequences by name list or BED region
- `mergepe` - Interleave two paired FASTA/FASTQ files
- `cutN` - Cut sequences at runs of N

**Test Coverage**: ~66% of statements (`go test -cover`)  
**Performance**: ~1.05-1.1x faster than original on the implemented commands  
**Documentation**: README with examples  

**Key Features**:

- Fast FASTA/Q processing for the eight commands above
- Quality score handling
- Memory-efficient streaming

**Migration Notes**:

- Command structure changed (subcommands instead of flags)
- Only the eight commands above are implemented; upstream seqtk has more
  (`mutfa`, `randbase`, `hpc`, ...)
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

**Test Coverage**: ~99.9% of statements (`go test -cover`)  
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

**Test Coverage**: ~82% of statements (`go test -cover`)  
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

**Test Coverage**: **100%** of statements (`go test -cover`)  
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

**Test Coverage**: ~67% of statements (`go test -cover`)  
**Performance**: ~1.1x  
**Documentation**: README with examples

**Key Features**:

- Adapter trimming (3' and 5')
- Quality filtering
- Sliding-window quality trimming (`--cut_front`, `--cut_tail`, `--cut_right`)
- Length filtering
- N content filtering
- Poly-G/X tail trimming (NovaSeq)
- Complexity filtering
- Built-in gzip support
- Comprehensive statistics

**Migration Notes**:

- Simplified version of original
- Core preprocessing features and sliding-window quality trimming implemented
- No HTML/JSON reports (future feature)
- Parallel worker pool exists; many upstream knobs are still missing

---

### 6. ✅ bedmerge

**Status**: Complete  
**Version**: 1.0.0  
**Original**: bedtools merge (C++)  
**Category**: Genomic Intervals / Utilities

**Implemented Commands**:

- Single command for merging BED intervals

**Test Coverage**: **100%** of statements (`go test -cover`)  
**Performance**: ~2x faster than bedtools merge on the implemented path  
**Documentation**: README with examples

**Key Features**:

- Merge overlapping BED intervals
- Distance-based merging (`-d` option)
- Strand-specific merging (`-s` option)
- Column aggregation `bedtools merge`-style: `-c`/`--columns` and `-o`/`--operations`
  with `sum`, `min`, `max`, `mean`, `median`, `count`, `count_distinct`,
  `distinct`, `collapse`, `first`, `last`, `mode`, `antimode`
- Statistics output (`-S` option)
- Built-in gzip support
- Automatic sorting

**Migration Notes**:

- Compatible with `bedtools merge` for the documented common path
- Output is BED3 by default; `-c`/`-o` adds the requested aggregated columns
- CLI note: `-c` short form now means `--columns` (matches `bedtools merge`);
  `--count` is still available by its long name

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

**Status**: Partial — a subset of upstream vcftools, ~50 of ~147 options

**Implemented Commands**:

- Single command with multiple filtering, statistics and conversion options

**Test Coverage**: ~58% of statements (`go test -cover`)  
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
  `--site-quality`, `--missing-site`, `--missing-indv`, `--depth`, `--geno-depth`,
  `--hardy`, `--site-pi`, `--window-pi`(+`--window-pi-step`), `--TajimaD`,
  `--TsTv-summary`, `--TsTv`, `--TsTv-by-count`, `--TsTv-by-qual`, `--het`,
  `--singletons`, `--hist-indel-len`, `--FILTER-summary`, `--SNPdensity`
- Population genetics: Weir & Cockerham 1984 Fst (`--weir-fst-pop` ×2+) per site
  and over windows (`--fst-window-size`/`--fst-window-step`); mean and weighted
  summary printed to stderr
- VCF recoding (`--recode`, `--recode-INFO-all`)
- Format conversion: `--012`, `--plink`, `--plink-tped` (with `--chrom-map`)

Every `Params` field declared in the package is now wired to real logic;
`checkUnsupported` no longer rejects anything. The remaining gap vs upstream
vcftools is dominated by **LD analysis** (`--geno-r2`, `--hap-r2`, all the
`--ld-window-*` options) and a long tail of less-common options.

See [FEATURE_COMPARISON.md](vcftools/FEATURE_COMPARISON.md) and
[ROADMAP.md](vcftools/ROADMAP.md) for the full picture.

**Migration Notes**:

- Per-site nucleotide diversity (`--site-pi`) uses the standard
  `(n^2 - Σ c_a^2) / (n(n-1))` formula; earlier builds reported a different
  (incorrect) per-genotype quantity
- Not a drop-in replacement: anything not in the list above is unavailable

---

### 9. ✅ New bedtools subcommands (May 2026)

Eight more `bedtools` subcommands were ported in this round, alongside the
existing `bedmerge` and `bedintersect`. Every one of them has its own
`tools/bedX/cmd/bedX/main.go` + `pkg/bedX/*.go` + `README.md`, follows the
POSIX-compliant CLI conventions in [`../docs/CLI_CONVENTIONS.md`](../docs/CLI_CONVENTIONS.md),
and reuses `pkg/bioformats/bed` + `pkg/bioformats/iohelper`.

| Tool | Maps to | Coverage | Highlights |
|------|---------|---------:|------------|
| `bedsort` | `bedtools sort` | 91.6% | Lex / size / score sort modes; `-g`/`--faidx` for chrom order |
| `bedslop` | `bedtools slop` | 95.2% | `-b N` / `-l N -r N` (+ `--pct`), `-s` strand swap, clip to chrom |
| `bedcomplement` | `bedtools complement` | 94.6% | Gaps over chroms in `-g`; errors if input not sorted |
| `bedsubtract` | `bedtools subtract` | 93.7% | A − B with `-A` / `-N` / `-s` / `-S`; splits A around B |
| `bedflank` | `bedtools flank` | 92.2% | Flank-only `slop` variant |
| `bedclosest` | `bedtools closest` | 92.8% | Sweep on sorted input; `-D ref/a/b`, `-N`, `-t all/first/last` |
| `bedgenomecov` | `bedtools genomecov` | 94.0% | histogram / `-bg` / `-bga` / `-d` / `-dz`; `-strand`, `-max`, `-scale`, `-5`/`-3` |
| `bedjaccard` | `bedtools jaccard` | 96.3% | Streaming sweep; `-s`/`-S`, `-f`/`-F` |

Smoke tests for each are hand-verified against expected output (see the
respective PRs and READMEs). Validated parity against the upstream `bedtools`
test suite is **still outstanding** — coverage ≠ byte-for-byte output match.

---

## Tool Comparison Matrix

| Tool | Original Lang | Go Version | Commands | Tests | Docs | Performance | Gzip |
|------|---------------|------------|----------|-------|------|-------------|------|
| seqtk | C | 1.0.0 | 8 | ✓ | ✓ | 1.05-1.1x | - |
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

- Eight subcommands implemented (`comp`, `fq2fa`, `seq`, `sample`, `trimfq`,
  `subseq`, `mergepe`, `cutN`); upstream subcommands still missing include
  `mutfa`, `randbase`, `hpc`

**PRINSEQ**:

- Covers the common filtering/trimming options (length, GC, N, quality,
  fixed/percentage/quality trimming, poly-N/A/T, dedup, paired-end, Phred+64,
  bad-sequence output, complexity filters); not every upstream option
- Graph/HTML report generation not included
- Quirk: `trimQualityLeft`/`trimQualityRight` always assume Phred+33 regardless
  of `QualType`; documented but not yet fixed

**sickle**:

- Built-in gzip support (by `.gz` extension)
- Phred-encoding auto-detect via `bufio.Reader.Peek` (default `-t auto`); explicit
  `sanger`/`illumina`/`solexa` still accepted; one-line stderr notice on detection
- Not validated byte-for-byte against the C original

**skewer**:

- Simplified adapter-matching algorithm
- `--auto-detect` picks from a small built-in adapter list (deterministic, by
  declaration order on ties); ~90% statement test coverage

**fastp**:

- Single-end and paired-end processing implemented
- Sliding-window quality trimming (`--cut_front`/`--cut_tail`/`--cut_right` with
  `--cut_window_size`/`--cut_mean_quality`)
- No HTML/JSON reports
- No automatic adapter detection
- Parallel worker pool exists but the feature surface is a subset of upstream

**bedmerge**:

- Column aggregation via `-c`/`--columns` and `-o`/`--operations` (sum, min,
  max, mean, median, count, count_distinct, distinct, collapse, first, last,
  mode, antimode); BED3 by default
- In-memory processing (not suitable for very large files)

**bedintersect**:

- Uses an interval tree (`pkg/bedintersect/intervaltree.go`)
- No reciprocal overlap mode
- No sorted-file streaming optimization
- In-memory B-file loading

**vcftools**:

- ~50 of ~147 upstream options; every declared `Params` field is now wired
- Largest remaining gap: LD analysis (`--geno-r2`, `--hap-r2`, `--ld-window-*`)

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

- **Test Coverage (statements)**: ~58-90% per package, ~77% unweighted average;
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
