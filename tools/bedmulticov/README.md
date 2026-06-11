# bedmulticov - Per-interval overlap counts against N input files

A Go re-implementation of `bedtools multicov`. For each interval in a
primary BED file (`-bed`), the tool reports the count of overlapping
records from each of N input files (`-files` / `-bams`). Output is the
A row verbatim, followed by one integer count column per input file.

## Features

- N input files (variadic `-files`), independent per-file counts.
- Strand filters: `-s` (same strand) and `-S` (opposite strand).
- Fraction-of-A (`-f`), fraction-of-B (`-F`) and reciprocal (`-r`)
  overlap thresholds.
- Pure Go, no third-party dependencies; per-chrom interval-tree index
  on each input file.
- Transparent gzip/BGZF and `-` for stdin.

## Build

```bash
go build ./tools/bedmulticov/cmd/bedmulticov
```

## Usage

```bash
bedmulticov -bed <A.bed> -files <B1.bed> [<B2.bed> ...] [options]
```

### Options

- `-bed FILE`            A intervals (required; `-` for stdin; `.gz` ok).
- `-files FILE..`        One or more B files (BED). Variadic; everything
  between `-files` and the next dash-prefixed token is treated as a path.
- `-bams FILE..`         Alias for `-files` (kept for upstream
  compatibility). BAM files (`.bam`) are auto-detected by extension.
- `-q, --mapq N`         Minimum MAPQ for BAM inputs (default 0; ignored
  for BED inputs).
- `-D, --max-depth N`    Cap the reported count per A interval per BAM
  input (default 64000; 0 disables the cap; ignored for BED inputs).
- `-s, --strand`         Same-strand overlaps only.
- `-S, --opposite`       Opposite-strand overlaps only.
- `-f FRAC`              Minimum fraction of A overlapped (0,1].
- `-F FRAC`              Minimum fraction of B overlapped (0,1].
- `-r, --reciprocal`     Apply `-f` to BOTH A and B.
- `-o, --output FILE`    Output file (default: stdout).
- `-h, --help`           Show help.
- `-v, --version`        Show version.

### Output

```text
chrom  start  end  [<A's extra cols>...]  <count_1>  <count_2>  ... <count_N>
```

## Examples

```bash
# Default: report overlap counts for each B against A.
bedmulticov -bed a.bed -files b1.bed b2.bed

# Same-strand only, at least 50% of A covered:
bedmulticov -bed a.bed -files b1.bed b2.bed -s -f 0.5
```

## Parity with upstream

Behavioural parity covers BED-vs-BED and BAM-vs-BED multi-coverage with
the standard strand / fraction options, plus the `-q` MAPQ filter and
`-D` per-A-interval depth cap for BAM inputs.

The `reference_code/bedtools/test/multicov/` corpus is entirely BAM-based
(`one_block.sam`, `two_blocks.sam`, `test-multi.sam`, `test-multi.2.sam`,
all run through `htsutil samtoindexedbam`). All ten upstream cases now
pass byte-for-byte against real in-memory BAM streams — see
`pkg/bedmulticov/parity_test.go`:

- `t1` default overlap, `t2` `-s`, `t3` `-S`, `t4` two-block alignment
  without `-split`, `t10` multi-input.
- `t5` `-split` alone, `t6` `-split -s`, `t7` `-split -S`, `t8`
  `-split -f 0.01`, `t9` `-split -f 0.10` — all exercise BAM CIGAR
  `N`-op block-aware coverage on `15M10N15M` alignments. With
  `-split`, each alignment is decomposed into its contiguous
  reference-consuming op-runs (M/=/X, with D extending the current
  block to mirror upstream's `breakOnDeletionOps=false`); the
  alignment is counted at most once per A interval and `-f` is
  applied to `total_overlap / sum_of_BAM_block_lengths` using strict
  `>` (a bedtools 2.x quirk in `multiBamCov.cpp::FindBlockedOverlaps`
  preserved here for byte-for-byte parity).

BED, BAM, **and** CRAM inputs are all supported (BAM/CRAM are decoded via
`pkg/htsgo/alnio`, which routes CRAM to `pkg/htsgo/cram`); the `-q` MAPQ
filter and `-D` per-interval depth cap are honoured for both. A
reference-backed CRAM may be passed a FASTA for base reconstruction.

## Tests

```bash
go test -race -cover ./tools/bedmulticov/...
```

Coverage: ~88.8% on `pkg/bedmulticov` (race + cover, 2026-05-15).

## Limitations

- Stranded fraction semantics follow `bedintersect`'s convention: with
  `-s` / `-S`, records with an empty strand on either side do not match.
