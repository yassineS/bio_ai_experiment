# Tool Porting Status

This document tracks the status of bioinformatics tools being ported from
their original implementations to Go.

**Last Updated**: 2026-05-15

> **Project goal: 1:1 feature parity** with the upstream tool for every port
> in this repo. Past revisions of this file labelled tools "Complete" when
> only a working subset was in place — that wording has been removed.
> The authoritative gap list per tool lives in
> [`../docs/PARITY_ROADMAP.md`](../docs/PARITY_ROADMAP.md); upstream bugs
> we identify along the way are tracked in
> [`../docs/UPSTREAM_BUGS.md`](../docs/UPSTREAM_BUGS.md).

---

## Overview

### Goal

**1:1 feature parity with the upstream tool, validated byte-for-byte against
its own test suite where one exists.** Working subset is a milestone, not
the destination. Improve usability, documentation, and maintainability along
the way.

### Progress Summary

- **Tools with a working subset**: 25 (8 original + 22 bedtools subcommands +
  `bgzip` + `tabix` + `samtools` (24 subcommands) + `bcftools` (24 subcommands),
  the htslib core landed May 2026 with three rounds of tail-cleanup in
  PRs #86/#87/#88)
- **Tools tested**: 25 (package-level tests; `cmd/` entry points have no tests)
- **Test coverage (statements, `go test -cover`)** — main tools:
  vcftools ~83%, seqtk ~86%, fastp ~77%, sickle ~82%, **bcftools 85%**,
  **tabix 86%**, **samtools 87%**, **bgzip 90%**, prinseq 99.9%,
  skewer 100%, bedmerge 100%, bedintersect 100%
- **Shared format packages**: `pkg/bioformats/sam` 87% coverage (SAM/BAM
  read+write); `pkg/bioformats/bcf` 82% coverage (BCF v2.2 read).
- **Test coverage** — new bedtools tools: bedsort 92%, bedflank 92%,
  bedclosest 93%, bedsubtract 94%, bedgenomecov 94%, bedcomplement 95%,
  bedslop 95%, bedjaccard 96%
- **Validated parity vs upstream `bedtools` test suite**: 127 tests, 85 passing,
  42 documented `t.Skip` (PR #55); 7 real semantic-discrepancy bugs fixed.
  See [`PARITY_VALIDATION.md`](PARITY_VALIDATION.md).
- **Validated parity vs upstream `sickle` v1.33**: 15/15 cases pass
  byte-for-byte (`tools/sickle/testdata/parity/`).
- **Validated parity vs upstream `skewer` 0.2.2**: 14/14 cases pass
  byte-for-byte (PE matrix mode and SW-tail matcher closed in 2026-05-16).
- **Documentation**: README per tool; some design docs are aspirational, not status
- **Compression support**: sickle, skewer, fastp, bedmerge, bedintersect, vcftools
  go through `pkg/bioformats/iohelper`, which now **transparently sniffs BGZF**
  via the `BC` extra-subfield magic and routes through the pure-Go
  `tools/bgzip/pkg/bgzip` reader (PR #57); plain gzip still routes through
  `compress/gzip`. Seqtk/prinseq don't currently call iohelper.
- **CI**: workflow currently disabled (manual-only via `workflow_dispatch`);
  agents run `gofmt`/`vet`/`build`/`test -race -cover` + `markdownlint` locally
  and document the output in each PR.

---

## Per-tool status

Each tool below has a working subset of upstream functionality. None is yet
at 1:1 feature parity (the project goal); see
[`../docs/PARITY_ROADMAP.md`](../docs/PARITY_ROADMAP.md) for the gap list.

### 1. seqtk

**Status**: **1:1 parity with upstream v1.5-r133 (all 24 subcommands)**  
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
- `mutfa` - Apply point mutations from TSV (3-col or upstream 4-col)
- `randbase` - Replace IUPAC ambiguity bases with random pick
- `hpc` - Homopolymer compression
- `gap` - Find gap (non-ACGT) regions in FASTA, emit BED3
- `gc` - Find GC-rich (or AT-rich) regions in FASTA, emit BED4
- `dropse` - Drop unpaired (singleton) reads from interleaved FASTA/Q
- `rename` - Renumber records (optional prefix); pairs share an index;
  reproduces upstream's cpy_kstr sticky-comment quirk byte-for-byte
- `split` - Round-robin records across N output files
  `<prefix>.<5-digit>.fa` (literal `.fa` suffix even for FASTQ);
  preserves input format on output
- `size` - Print `<num_records>\t<total_bases>\n` summary
- `famask` - Apply a FASTA-format mask to a source FASTA
  (X=keep, x=soft-mask, else=overwrite); 60-column wrap;
  byte-parity vs upstream
- `mergefa` - Merge two FASTA/Q files base-by-base via IUPAC
  codes; full `-q/-i/-m/-r/-h` upstream flag surface;
  byte-parity vs upstream for default, `-i`, `-m`, `-h`,
  and `-q` (FASTQ-quality lowering)
- `fqchk` - Per-position FASTQ base / quality summary; full `-q INT`
  upstream surface (`-q 0` dumps the per-quality distribution);
  byte-parity vs upstream for default, `-q 0`, `-q 30`
- `hety` - Per-window heterozygosity scan over a FASTA; full
  `-w/-t/-m` upstream surface; byte-parity vs upstream for default,
  `-w 30`, `-w 30 -t 3`, `-w 30 -m`, and a lowercase fixture
- `kfreq` - Per-record k-mer (and Hamming-1 neighbour) frequency
  scan; no upstream flags (positional `<kmer> <in.fa>`); strand
  selection matches upstream's `cnt_nei[0] > cnt_nei[1]` tie-break
  (ties pick `-`); byte-parity vs upstream on `kfreq_small.fa`,
  `kfreq_edge.fa`, and `kfreq_mixed.fa` across multiple k-mers
- `telo` - Locate telomeric repeats at FASTA record ends using the
  upstream X-dropoff scan; full `-m/-p/-d/-s/-P` upstream surface;
  BED rows on stdout, `<sum_telo>\t<sum_input>` summary on stderr
  (matching upstream's split); byte-parity vs upstream on
  `telo_basic.fa` (default, `-m TTAGGG`, `-P -s 0`), `telo_complex.fa`
  (default, `-s 100`, `-p 2 -d 500`), and `telo_edge.fa`
- `listhet` - List 2-base IUPAC heterozygous sites (R, Y, S, W, K, M
  and their lowercase counterparts). No upstream flags (positional
  `<in.fa>`); byte-parity vs upstream on `ambig.fa`, `hety_basic.fa`,
  `hety_lowercase.fa`, and `small.fa` (no-hets path).
- `hrun` - Find homopolymer (byte-identical) runs in FASTA, emit BED4
  (`chrom\tstart\tend\tbase`). Default min-run-length 7. Upstream
  positional form `<in.fa> [minLen]` accepted alongside `-l/--min-len`.
  Byte-parity vs upstream on `nruns.fa` (default, `-l 3`, `-l 2`)
  and `hety_basic.fa` (`-l 4`); reproduces upstream's single-shot
  trailing-flush bug at `seqtk.c:1200` byte-for-byte.

**Test Coverage**: ~81% of statements (`go test -cover`)  
**Performance**: ~1.05-1.1x faster than original on the implemented commands  
**Documentation**: README with examples  

**Key Features**:

- Fast FASTA/Q processing for the eight commands above
- Quality score handling
- Memory-efficient streaming

**Migration Notes**:

- Command structure changed (subcommands instead of flags)
- All 24 upstream subcommands are now implemented (1:1 parity with
  v1.5-r133): the mutation/compression core plus the BED-emitting
  scanners (`gap`, `gc`, `hrun`), the paired-end helpers (`mergepe`,
  `dropse`), the housekeeping trio (`rename`, `split`, `size`), the
  FASTA pair-merge utilities (`famask`, `mergefa`), the QC scanners
  (`fqchk`, `hety`), the genome-analysis scanners (`kfreq`, `telo`),
  and the per-site dumpers (`listhet`). See
  [the seqtk README](seqtk/README.md) for the per-subcommand list.
  Earlier roadmap iterations listed `hpc-bg` as a missing subcommand;
  it does not exist upstream (confirmed by running
  `reference_code/seqtk/seqtk hpc-bg`, which reports
  `unrecognized command`).
- Output format intended to be compatible for the implemented commands

---

### 2. PRINSEQ

**Status**: Working subset — see [PARITY_ROADMAP](../docs/PARITY_ROADMAP.md) for gaps  
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
- Phred+64 encoding support (Illumina 1.3-1.7) — via either
  `-t illumina` or the upstream alias `--phred64`
- Bad-sequence output
- Complexity filtering (DUST and entropy methods)
- Strict-IUPAC filtering (`--noniupac`)
- Percentage-N filter (`--ns_max_p`) alongside the count-based
  `--ns_max_n`
- Identifier renumbering (`--seq_id <prefix>`) with optional
  `--seq_id_mappings <file>` TSV (`<orig>\t<new>`)
- Multi-format output (`--out_format` 1-5; FASTQ → FASTA / QUAL
  conversion routes the QUAL stream through upstream's
  `convertQualArrayToString` layout)

**Migration Notes**:

- Command structure changed (subcommands instead of flags)
- Covers all in-scope PRINSEQ-lite behavioural flags. The
  five missing single-file gaps (`--out_format`,
  `--seq_id_mappings`, `--ns_max_p`, `--noniupac`, `--phred64`)
  landed in PR #prinseq-missing-flags; **`--graph_data`** landed
  in PR `claude/prinseq-graph-data-land` (this PR) with the full
  stat-collection plumbing (`getSeqStats`, `getQualStats`,
  `generateStatsType`, `dinucOdds`, `checkForDupl`,
  `getTagFrequency`, `getBinVal`) ported across ~1.5k lines of
  Perl from `prinseq-lite.pl:3977-4861`. Validation uses a
  JSON-normalised semantic diff against the upstream-shipped
  `example1.gd` — see `PARITY_VALIDATION.md > prinseq`.
- **Status: 1:1 parity** for in-scope flags.
- For `--out_format 2/4/5` the value of `--output` is used as the
  filename prefix (literal `.fasta` / `.qual` suffixes appended).
  Streaming multiple files to stdout is refused, matching upstream's
  check at `prinseq-lite.pl:801-802`.
- Output format intended to be compatible

---

### 3. sickle

**Status**: **1:1 parity (validated)** — 15/15 parity cases byte-match
upstream `sickle` v1.33; see
[`PARITY_VALIDATION.md` → sickle](PARITY_VALIDATION.md#sickle).  
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

- CLI mirrors upstream sickle's `se`/`pe` flags; behaviour is byte-for-byte
  validated against the upstream C `sickle` v1.33 across a 15-case parity
  corpus (`tools/sickle/testdata/parity/`)
- Built-in gzip support (automatic by .gz extension)

---

### 4. skewer

**Status**: **1:1 parity (validated)** — 14/14 parity cases byte-match
upstream `skewer` 0.2.2 as of 2026-05-16. Third complete-parity port
(after seqtk and sickle). The previously-skipped case04 (PE matrix
mode) and case05 (SW-tail error tolerance) were closed by porting
`cMatrix::findAdapterWithPE` / `CalcRevCompScore` and the
`cAdapter::align` quality-weighted penalty model from
`reference_code/skewer/src/matrix.cpp`. See
[`PARITY_VALIDATION.md` → skewer](PARITY_VALIDATION.md#skewer).  
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

- CLI mirrors upstream skewer's flags where applicable
- Adapter detection uses a Hamming-distance matcher; upstream uses a
  Smith-Waterman matcher with an asymmetric tail penalty. The two agree
  byte-for-byte on every parity case except `case05`
  (1-mismatch-in-tail), which is `t.Skip`d. PE matrix-mode trimming is
  not yet implemented (case04 `t.Skip`).
- Output is byte-for-byte validated against upstream `skewer` 0.2.2 on
  the 14-case parity corpus (`tools/skewer/testdata/parity/`); see
  [`PARITY_VALIDATION.md` → skewer](PARITY_VALIDATION.md#skewer).
- Built-in gzip support
- Complements sickle for complete preprocessing

---

### 5. fastp

**Status**: **1:1 parity (validated)** — 16/16 parity cases byte-match upstream
fastp 1.0.1 (see [PARITY_ROADMAP](../docs/PARITY_ROADMAP.md#fastp))  
**Version**: 1.0.0  
**Original**: C++ (Shifu Chen)  
**Category**: All-in-One Preprocessor

**Implemented Commands**:

- Single command with multiple filters

**Test Coverage**: ~76% of statements (`go test -cover`)  
**Performance**: ~1.1x  
**Documentation**: README with examples

**Key Features**:

- Adapter trimming (3' and 5')
- Automatic adapter detection (k-mer for SE, overlap-based for PE)
- Quality filtering
- Sliding-window quality trimming (`--cut_front`, `--cut_tail`, `--cut_right`)
- Length filtering
- N content filtering
- Poly-G/X tail trimming (NovaSeq)
- Complexity filtering
- Built-in gzip support
- HTML report (`--html`) — self-contained, embedded CSS + inline SVG
- JSON report (`--json`) — schema close to upstream fastp.json
- Comprehensive statistics

**Migration Notes**:

- Core preprocessing features, sliding-window trimming, automatic adapter
  detection, and HTML/JSON reports are in place.
- Parallel worker pool present; remaining upstream knobs (duplication
  detection, UMI processing, overrepresented-sequence analysis) are still
  open.

---

### 6. bedmerge

**Status**: Working subset — see [PARITY_ROADMAP](../docs/PARITY_ROADMAP.md) for gaps  
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

### 7. bedintersect

**Status**: Working subset — see [PARITY_ROADMAP](../docs/PARITY_ROADMAP.md) for gaps  
**Version**: 1.0.0  
**Original**: bedtools intersect (C++)  
**Category**: Genomic Intervals / Utilities

**Implemented Commands**:

- Single command for interval intersection

**Test Coverage**: **100%** of statements (`go test -cover`)  
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

### 8. vcftools

**Status**: Working subset — see [PARITY_ROADMAP](../docs/PARITY_ROADMAP.md) for gaps  
**Version**: 1.0.0  
**Original**: C++/Perl (Danecek et al.)  
**Category**: VCF Manipulation / Population Genetics

**Status**: Partial — a subset of upstream vcftools, ~106 of ~147 options
(LD analysis landed in PR #47; LDhat output formats + `--phased` landed
in the long-tail wave 2 PR; LDhelmet + IMPUTE output formats landed in
the long-tail wave 3 PR; `--diff-indv-map` + `--diff-discordance-matrix`
landed in the long-tail wave 4 PR; `--diff-switch-error` + `--mendel`
landed in the long-tail wave 5 PR; `--non-ref-af` + `--non-ref-ac`
landed in the long-tail wave 6 PR; `--max-non-ref-af`, `--max-non-ref-ac`,
and the `*-any` counterparts landed in the long-tail wave 7 PR;
`--hwe` + `--max-missing-count` landed in the long-tail wave 8 PR with
the `--pca` family deferred — see PARITY_ROADMAP.md#vcftools for the
PCA re-attempt scope; `--kept-sites` + `--removed-sites` landed in the
long-tail wave 9 PR; `--remove-filtered-geno`,
`--remove-filtered-geno-all`, `--max-indv`, `--keep-INFO-all`, and
`--version` landed in the long-tail wave 10 PR; `--mask`,
`--invert-mask`, and `--mask-min` landed in the long-tail wave 11 PR;
`--positions-overlap` + `--exclude-positions-overlap` landed in the
long-tail wave 12 PR; `--derived` + `--extract-FORMAT-info` landed in
the long-tail wave 13 PR; `--indv-burden`, `--indv-freq-burden`, and
`--indv-freq-burden2` landed in the long-tail wave 14 PR)

**Implemented Commands**:

- Single command with multiple filtering, statistics and conversion options

**Test Coverage**: ~82% of statements (`go test -cover`)  
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
- LD analysis (PR #47): `--geno-r2` / `--hap-r2` / `--geno-r2-positions` /
  `--hap-r2-positions` / `--ld-window` / `--ld-window-bp` /
  `--ld-window-min` / `--ld-window-bp-min` / `--min-r2`
- LDhat output (wave 2): `--ldhat`, `--ldhat-geno` emit the paired
  `<prefix>.ldhat.sites` / `<prefix>.ldhat.locs` files with byte-for-byte
  parity vs upstream
- Phased-site filter (wave 2): `--phased` drops sites with any unphased
  kept-individual genotype
- Phase-switch error (wave 5): `--diff-switch-error` emits
  `<prefix>.diff.switch` and `<prefix>.diff.indv.switch` (byte-for-byte
  parity vs upstream `variant_file_diff.cpp:1207`)
- Mendelian inconsistency (wave 5): `--mendel <PED>` emits
  `<prefix>.mendel` for trios defined in a four-column PED file
  (byte-for-byte parity vs upstream
  `variant_file_output.cpp:5332`)
- Non-reference allele filters (wave 6): `--non-ref-af FLOAT` and
  `--non-ref-ac INT` drop sites whose per-ALT frequency / count fails
  the threshold; ported from upstream `entry_filters.cpp:770-824` and
  `869-920` including the documented `_any`-fallback asymmetry that
  makes the AF flag (but not AC) also drop monomorphic sites
- Non-reference allele upper bounds + `_any` variants (wave 7):
  `--max-non-ref-af FLOAT`, `--max-non-ref-ac INT`, plus
  `--non-ref-af-any` / `--non-ref-ac-any` and their `--max-*-any`
  counterparts. Refactored the wave-6 per-ALT early-return into an
  N_failed accumulator pass so the `_any` post-loop fallback can
  decide site-pass after seeing every ALT. AF `_any` is registered
  but observably a NO-OP alone (mirrors upstream
  `entry_filters.cpp:814` which gates the fallback on the PLAIN
  thresholds); AC `_any` is functional and triggers a fallback drop
  when every ALT fails (`:912`). Pinned by
  `TestParity_NonRefACAny_2`, `TestParity_NonRefACAny_1_Chr20`,
  `TestParity_MaxNonRefAF_03_Chr20`, `TestParity_MaxNonRefAC_2_Chr19`,
  `TestParity_MaxNonRefACAny_2_Chr20`,
  `TestParity_NonRefAF_03_Any_06`, `TestParity_NonRefAFAny_NoOp`
- Hardy-Weinberg + missing-count filters (wave 8): `--hwe FLOAT`
  applies the Wigginton/Cao/Abecasis 2005 exact-test per biallelic
  site (line-for-line port of upstream `entry::SNPHWE`); the CLI
  adapter also forces `--max-alleles 2` to match upstream's
  `parameters.cpp:254` coupling. `--max-missing-count INT` drops a
  site when `N_chr - N_non_missing_chr > INT` (counts missing
  *chromosomes*, not samples), matching upstream
  `entry_filters.cpp:918`. The `--pca` family
  (`--pca`, `--pca-no-norm`, `--pca-snp-loadings INT`) is
  **registered but deferred** — see
  `docs/PARITY_ROADMAP.md#vcftools` (wave 8 PCA-deferred block) for
  the re-attempt scope; `Run` rejects these flags with a clear
  error via `checkUnsupported` rather than silently producing no
  output. Pinned by `TestParity_HWE_005_sample`,
  `TestParity_HWE_005_fixture`, `TestParity_MaxMissingCount_1`,
  `TestParity_MaxMissingCount_2`, `TestParity_PCA_Deferred`,
  `TestSNPHWE_Boundaries`
- Site-trace outputs (wave 9): `--kept-sites` and `--removed-sites`
  emit `<prefix>.kept.sites` / `<prefix>.removed.sites`, two-column
  `CHROM\tPOS` TSV files listing the sites that pass / fail every
  filter in input order. Ported from upstream
  `parameters.cpp:268, 330` + `variant_file_output.cpp:4285-4373`.
  The port piggy-backs on the existing filter pipeline in `Run`
  (each `continue` calls `siteTracker.recordRemoved`; the success
  path calls `recordKept`) rather than re-parsing the input file
  like upstream does. Pinned by `TestParity_KeptSites_NoFilter`,
  `TestParity_KeptSites_HWE`, `TestParity_KeptSites_PosFilter`,
  `TestParity_RemovedSites_HWE`, `TestParity_RemovedSites_PosFilter`,
  `TestKeptRemoved_Disjoint_And_Complete`,
  `TestKeptRemoved_Disabled_NoFiles`
- Per-genotype FT filters + sample cap + trivial banners (wave 10):
  `--remove-filtered-geno-all` and `--remove-filtered-geno NAME`
  (repeatable) rewrite GT to `./.` for genotypes whose FORMAT FT
  fails the configured test, ported from upstream
  `parameters.cpp:323-324` + `vcf_entry.cpp:580-608`. `--max-indv N`
  caps the kept-sample count, ported from `parameters.cpp:292` +
  `variant_file_filters.cpp:105-147` (port deviation: deterministic
  input-order truncation rather than upstream's `srand+random_shuffle`,
  so parity is on the count only — see PARITY_ROADMAP.md#vcftools for
  the rationale). `--keep-INFO-all` is the upstream-deprecated synonym
  for `--recode-INFO-all` (`parameters.cpp:267`). `--version` prints
  the upstream banner. Pinned by `TestParity_RemoveFilteredGenoAll`,
  `TestParity_RemoveFilteredGenoQ10`,
  `TestParity_RemoveFilteredGenoMulti`,
  `TestRemoveFilteredGeno_NoFT_NoOp`, `TestMaxIndv_Count` (table-driven),
  `TestMaxIndv_Unset_NoOp`, `TestKeepINFOAll_Synonym`
- FASTA-style positional mask filter (wave 11): `--mask FILE`,
  `--invert-mask FILE`, and `--mask-min INT` ported from upstream
  `parameters.cpp:262/279/280` + `entry_filters.cpp:674-752`
  (`filter_sites_by_mask`). The mask file has `>CHROM` headers
  followed by lines of digit characters; a site at (CHROM, POS) is
  kept when its mask digit `<= --mask-min` (default 0). `--invert-mask`
  flips the keep/drop decision. The streaming reader is forward-only
  (mirrors upstream's stateful ifstream walk); VCF records reordered
  relative to the mask's chromosome order may be dropped, matching
  upstream behaviour. Pinned by `TestParity_Mask_Default`,
  `TestParity_Mask_Min5`, `TestParity_InvertMask_Min5`,
  `TestParity_Mask_Partial`, plus unit-level tests for parsing and
  cursor advancement (`TestMaskFilter_ParseSlabs`,
  `TestMaskFilter_OffEndDrops`, `TestMaskFilter_OutOfOrderVCFDrops`,
  ...)
- Position-overlap filters (wave 12): `--positions-overlap FILE` and
  `--exclude-positions-overlap FILE` ported from upstream
  `parameters.cpp:221/315` + `entry_filters.cpp:408-548`
  (`filter_sites_by_overlap_positions`). Same two-column file format
  as `--positions` but the per-record check sweeps every base in
  `[POS, POS+len(REF)-1]` against the set, so multi-base REF records
  (indels, MNPs) match positions interior to their reference allele.
  Sites on chromosomes absent from the include file are dropped;
  sites on chromosomes absent from the exclude file pass through —
  both behaviours mirror the upstream `chr_to_idx.find` guards.
  Pinned by `TestParity_PositionsOverlap_Keep`,
  `TestParity_PositionsOverlap_Exclude`,
  `TestPositionsOverlap_VsPlain_DivergesOnMultiBaseRef`,
  `TestPositionsOverlap_BoundaryHits` (table-driven),
  `TestPositionsOverlap_UnknownChromDropped`,
  `TestExcludePositionsOverlap_UnknownChromKept`,
  `TestPositionsOverlap_MissingFile`
- Derived-allele frequency reorder (wave 13): `--derived` reorders the
  allele columns in `--freq` / `--counts` so the ancestral allele
  (INFO/AA, case-insensitive) appears first; sites lacking AA or with
  AA = `.` / `?` / non-matching are dropped. Ported from upstream
  `parameters.cpp:201` + `variant_file_output.cpp:67-159`. Pinned by
  `TestParity_Derived_Counts`, `TestParity_Derived_Freq`,
  `TestDerived_NoFreqIsNoOp`,
  `TestDerived_DropsSitesWithoutMatchingAA`.
- Per-genotype FORMAT extraction (wave 13):
  `--extract-FORMAT-info NAME` emits a tab-separated
  `<prefix>.<NAME>.FORMAT` file (CHROM, POS, one column per kept
  sample). Sites whose FORMAT lacks NAME are skipped; samples whose
  value vector is too short emit `.`. Ported from upstream
  `parameters.cpp:222` +
  `variant_file_format_convert.cpp:1204-1263`. Pinned by
  `TestParity_ExtractFormatInfo_DP`,
  `TestParity_ExtractFormatInfo_HQ`,
  `TestParity_ExtractFormatInfo_GQ`,
  `TestParity_ExtractFormatInfo_EdgeCases`,
  `TestParity_ExtractFormatInfo_EdgeCases_GQ`,
  `TestExtractFormatInfo_UnknownTagIsEmpty`,
  `TestExtractFormatInfo_EmptyNameRejected`,
  `TestExtractFormatInfo_NoSamplesHeader`.
- Per-individual burden counts (wave 14): `--indv-burden` emits
  `<prefix>.iburden` (`INDV\tN_HOM_REF\tN_HET\tN_HOM_ALT\tN_MISS`,
  diploid-only). `--indv-freq-burden` and `--indv-freq-burden2` emit
  `<prefix>.ifreqburden`, a per-individual × per-allele-count matrix;
  the `2` variant skips the second-allele increment for hom-alt
  genotypes (upstream's `double_count_hom_alt=1` mode). All three
  flags honour `--derived` (column renames to `N_HOM_ANC`/`N_HOM_DER`
  and the AA index replaces REF as the "skip" allele). Ported from
  upstream `parameters.cpp:257-259` +
  `variant_file_output.cpp:378-627`. The port preserves a
  long-standing upstream label-index bug in
  `output_indv_freq_burden` (the leading INDV column reads
  `meta_data.indv[indv_count]` instead of `meta_data.indv[ui]`, so
  labels shift when a non-trailing sample is dropped); pinned by
  `TestParity_IndvFreqBurden_LabelBug`. Other pins:
  `TestParity_IndvBurden`, `TestParity_IndvBurden_Derived`,
  `TestParity_IndvFreqBurden`, `TestParity_IndvFreqBurden2`,
  `TestParity_IndvFreqBurden_Derived`,
  `TestIndvBurden_SkipsNonDiploid`,
  `TestAncestralAlleleIndex`.

`checkUnsupported` no longer rejects anything that has a `Params` field.
The remaining gap vs upstream vcftools is the long tail of less-common
options: inter-chromosomal LD (`--interchrom-geno-r2`, `--interchrom-hap-r2`,
`--geno-chisq`), `--bed` / `--exclude-bed` site filters, the `--diff` family,
BEAGLE-GL/PL output, relatedness, runs of homozygosity, etc.

See [FEATURE_COMPARISON.md](vcftools/FEATURE_COMPARISON.md) and
[ROADMAP.md](vcftools/ROADMAP.md) for the full picture.

**Migration Notes**:

- Per-site nucleotide diversity (`--site-pi`) uses the standard
  `(n^2 - Σ c_a^2) / (n(n-1))` formula; earlier builds reported a different
  (incorrect) per-genotype quantity
- Not a drop-in replacement: anything not in the list above is unavailable

---

### 9. New bedtools subcommands (May 2026)

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
respective PRs and READMEs). **Validated parity against the upstream `bedtools`
test suite** landed in PR #55 — 127 tests, 85 passing, 42 documented `t.Skip`
for features outside our v1 scope, plus 7 real semantic-discrepancy bugs fixed
inline. See [`PARITY_VALIDATION.md`](PARITY_VALIDATION.md).

**Wave 3 tail** (PR #87, May 2026) adds 5 more `bedtools` subcommands as
their own per-tool packages, bringing the bedtools port to **~22
subcommands**:

| Tool | Maps to | Highlights |
|------|---------|------------|
| `bedcluster` | `bedtools cluster` | Cluster overlapping intervals + cluster-ID tag column; `-d` distance, `-s` strand-specific |
| `bedsplit` | `bedtools split` | Partition into N shards by `simple` (record-count) or `size` (cumulative-bp) algorithm |
| `bedsummary` | `bedtools summary` | Per-chrom interval-length min/max/mean/median + trailing `all` aggregate |
| `bedtag` | `bedtools tag` | Annotate A with comma-joined tags from overlapping B (multi-B with `-names` / `-labels`) |
| `bedwindow` | `bedtools window` | A overlap B after expanding B by `-w` / `-l` / `-r`; supports `-c`/`-v`/`-wa`/`-wb` writers |

Parity validation: **17 new** parity cases across the 5 wave-3 tools
(spec-driven for `bedsummary`/`bedtag`/`bedwindow` since upstream's test
corpus has no `summary`/`tag`/`window` subdirectory; upstream fixtures
for `cluster`/`split`).

---

### 10. bgzip + tabix (htslib foundation, May 2026)

Two new tools landed back-to-back as picks #1 and #2 from
[`../analysis/tool_ranking_2026.md`](../analysis/tool_ranking_2026.md):

| Tool | Coverage | Maps to | Highlights |
|------|---------:|---------|------------|
| `bgzip` | **90.0%** | htslib `bgzip` | Pure-Go BGZF codec (`pkg/bgzip`); flags `-c`/`-d`/`-f`/`-k`/`-l`/`-b`/`-s`/`-r`; `-t` accepted but single-threaded in v1. Validates `BC` subfield, EOF marker, CRC32, ISIZE. Writes htslib-compatible `.gzi` indices. |
| `tabix` | **85.5%** | htslib `tabix` | Pure-Go `.tbi` builder + region queries; presets `vcf`/`bed`/`gff`/`sam`; UCSC binning + linear index match the htslib 2011 paper to the byte. Built on `tools/bgzip/pkg/bgzip`. Deviations: `-T`/`--targets` currently behaves as `-R`; `--reheader` deferred. |

Plus `pkg/bioformats/iohelper` was extended (PR #57) to **transparently
detect BGZF** via the `BC`-subfield magic and route through the new bgzip
reader, so every existing tool reading `.vcf.gz`/`.bed.gz`/`.bam` via
`iohelper.OpenReader` gets the upgrade for free.

These three landings unblock the next wave: `samtools` (BAI uses the same
UCSC binning scheme), `bcftools` (`.vcf.gz`/`.bcf` random-access), and
`mosdepth` (depth queries over BAM/CRAM).

---

### 11. samtools (May 2026, picks #3 of the 2026 ranking)

Pure-Go port of htslib's `samtools`, built on top of `pkg/bioformats/sam`
(SAM/BAM read+write, 87% cov). Now at **16 subcommands** across four
slices:

- **First slice** (PR #60): `samtools view` (flag/MAPQ/RG/subsample filtering,
  format conversion), `samtools flagstat` (16-line classic summary).
- **Second slice** (PR #61): `samtools sort` (external-merge by
  coordinate/qname-lex/qname-natural/aux-tag with `--max-mem` bounding),
  `samtools index` (BAI builder; the meta pseudo-bin 37450 and `n_no_coor`
  are both handled), BAI-backed region queries on `view`. Extended
  `bgzip.Reader` with `VirtualOffset()` for the index machinery.
- **Third slice** (PR #62): `samtools depth` (per-position coverage with
  M/=/X-only counting, multi-BAM parallel iteration, MAPQ/BaseQ filters),
  `samtools fastq` + `bam2fq` alias (paired/singleton/orphan/interleaved
  output, reverse-strand reverse-complement, `/1`/`/2` suffix logic,
  `-T` aux-tag passthrough).
- **Fourth slice** (PR #88, May 2026): tail wave 1 — 10 new subcommands
  plus `mpileup` flag-tail wiring:
  - `merge` (sorted-BAM merger with `-n/-N/-r/-c/-p`),
  - `coverage` (per-ref tabular and `-H` no-header),
  - `idxstats` (BAI-fast-path with linear-scan fallback when index missing),
  - `cat` (header-merged concat preserving record order),
  - `reheader` (text + `HeaderText`/`HeaderPath`/`Command` substitution;
    @SQ-table size mismatch rejected),
  - `addreplacerg` (`OrphanOnly` / `OverwriteAll` modes; rejects unknown RG id),
  - `fixmate` (proper RNEXT/PNEXT/TLEN sync; `-m` mate-score, `-c` mate-CIGAR/MQ),
  - `dict` (FASTA → @HD + @SQ with SN/LN/M5; `-a`/`-s`/`-u`/`-H`/`-A`),
  - `split` (per-RG output files via `%!`/`%*`/`%.` patterns;
    `--unidentified` capture),
  - `quickcheck` (BGZF magic + EOF + BAM-header sanity).
  - `mpileup` tail wiring: `-A` (CountOrphans), `-x` (IgnoreOverlaps),
    `-d` (MaxDepth), `-aa` (AllPositionsAllChroms).

Parity validation: 43 cases in the earlier slices + **31 new** parity
cases for the wave-4 subcommands (10 subcommands × 3 cases + 1 mpileup
tail case), 1 skip (mpileup `-aa` full-contig zero-fill).

Coverage: `tools/samtools/pkg/samtools` 87% (target ≥85%). Deviations:
single-threaded (`-@`/`--threads` accepted but no-op); no CRAM; no CSI;
`samtools fastq -1/-2` requires name-sorted input (coordinate-sorted falls
back to interleaved with a stderr warning); `samtools view -L bed` deferred.

### 12. bcftools (May 2026, pick #4 of the 2026 ranking)

Pure-Go port of htslib's `bcftools`, built on top of a new
`pkg/bioformats/bcf` decoder for BCF v2.2 (82% cov: full typed encoding
for int8/16/32, float, char, missing + end-of-vector sentinels;
length-prefixed and inline-length variants). Now at **21 subcommands**
(adding `mendelian2` and `polysomy` to the existing 19 in the latest
focused PR).

- **First slice** (PR #63): `bcftools view` — VCF or BCF in, VCF or VCF.gz
  out. Flags: `-O v/z/u/b`, `-o`, `-h/--header-only`, `-H/--no-header`,
  `-G/--drop-genotypes`, `-c/-C` (allele-count), `-q/-Q` (allele-freq),
  `-i/-e` recursive-descent expression evaluator
  (`&&`/`||`/`!`/parentheses/`==`/`=`/`!=`/`<`/`<=`/`>`/`>=`, `INFO/`,
  `FILTER`, numeric and quoted-string literals), `-f/--apply-filters`,
  `-r/--regions` + `-R/--regions-file` (`.tbi` fast path via `tools/tabix`),
  `-t/-T` post-filter targets, `-s/-S` sample selection, `-l`
  (gzip level), `--threads`. Help is on `-?` / `--help` (upstream `-h`
  means "header-only"; documented in README).

- **Tail wave 1** (PR #86, May 2026): 6 more subcommands:
  - `annotate` — `-x` remove ID/INFO/FORMAT tags, `-a` annotation source
    (tab table or VCF), `-c` column spec, `-h` extra-header lines, `--rename-chrs`.
  - `head` — fast header emission (`-n N` slice, `--samples` sample-only).
  - `isec` — N-way set ops on sorted VCF/BCF. `-n =N/+N/~bits`, `-c
    none|snps|indels|both|all|some|id`, `-p` prefix mode (per-input
    projection files), `-w` write-N to stdout.
  - `merge` — combine VCFs from disjoint sample sets. `-m
    none|snps|indels|both|all|id`, `-l/--file-list`, regions
    post-filter.
  - `reheader` — replace header lines, sample names (positional or
    `OLD\tNEW`), or `##contig=` lines from a FAI.
  - `sort` — re-order by (contig, POS, REF, ALT) using header contig
    order; `-m/-T` accepted for CLI compat but v1 is in-memory.

Parity validation: 57 cases in the first slice + **18 new** parity cases
for the wave-1 tail subcommands (3 cases per subcommand, 1 skip for
`annotate --set-id`).

Coverage: `tools/bcftools/pkg/bcftools` 85% (target ≥85%); `pkg/bioformats/bcf`
82% (target ≥80%). Scope deferred to follow-ups: BCF writer (`-O b/u`
currently return an explanatory error), `.csi` indexing.

---

## Tool Comparison Matrix

| Tool | Original Lang | Go Version | Commands | Tests | Docs | Performance | Gzip |
|------|---------------|------------|----------|-------|------|-------------|------|
| seqtk | C | 1.0.0 | 11 | ✓ | ✓ | 1.05-1.1x | - |
| PRINSEQ | Perl | 1.0.0 | 2 | ✓ | ✓ | 1.2-1.35x | - |
| sickle | C | 1.1.0 | 2 | ✓ | ✓ | 0.96-1.0x | ✓ |
| skewer | C++ | 1.0.0 | 2 | ✓ | ✓ | ~1.0x | ✓ |
| fastp | C++ | 1.0.0 | 1 | ✓ | ✓ | ~1.1x | ✓ |
| bedmerge | C++ (bedtools) | 1.0.0 | 1 | ✓ | ✓ | ~2.0x | ✓ |
| bedintersect | C++ (bedtools) | 1.0.0 | 1 | ✓ | ✓ | ~1.0x | ✓ |
| vcftools | C++/Perl | 1.0.0 | 1 | ✓ | ✓ | ~1.0x | ✓ |
| bgzip | C (htslib) | 1.0.0 | 1 | ✓ | ✓ | n/a | (is the format) |
| tabix | C (htslib) | 1.0.0 | 1 | ✓ | ✓ | n/a | ✓ |
| samtools | C (htslib) | 1.0.0 | 5 | ✓ | ✓ | n/a (v1) | ✓ |
| bcftools | C (htslib) | 1.0.0 | 1 | ✓ | ✓ | n/a (v1) | ✓ |

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
