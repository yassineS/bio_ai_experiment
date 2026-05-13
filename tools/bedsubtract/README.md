# bedsubtract - Subtract B intervals from A intervals

A Go re-implementation of `bedtools subtract`.

## Description

For each interval in A, subtract any overlap with intervals in B and emit the
remaining segments (in A's input order). When a B interval punches a hole in
the middle of an A interval, the A interval is split into multiple output rows
with all of A's non-coordinate columns preserved.

## Installation

```bash
go build ./tools/bedsubtract/cmd/bedsubtract
```

## Usage

```text
bedsubtract -a <fileA.bed> -b <fileB.bed> [options]
```

## Options

- `-a, --a FILE` - Input BED file A (use `-` for stdin)
- `-b, --b FILE` - Input BED file B (use `-` for stdin)
- `-o, --output FILE` - Output BED file (`-` for stdout, default: stdout)
- `-A` - If any part of A overlaps B, drop the entire A interval
- `-N, --min-fraction NUM` - Only subtract when overlap covers at least NUM
  (0..1) of A
- `-s, --strand` - Only subtract same-strand B intervals (BED6+)
- `-S` - Only subtract opposite-strand B intervals (BED6+)
- `-h, --help` - Show help message
- `-v, --version` - Show version (`1.0.0`)

## Examples

```bash
# Subtract peaks from genes
bedsubtract -a genes.bed -b peaks.bed > genes_minus_peaks.bed

# Drop A entries that touch B at all
bedsubtract -a genes.bed -b blacklist.bed -A > clean_genes.bed

# Strand-aware subtraction
bedsubtract -a a.bed -b b.bed -s > out.bed

# Stream A from stdin
cat a.bed | bedsubtract -a - -b b.bed > out.bed
```

## Format

- Input: BED (tab-delimited, minimum 3 columns); `.gz` files supported.
- Output: same column layout as A, with `start`/`end` rewritten as A is split.

## Notes

- Coordinates are 0-based half-open `[start, end)`.
- Only one of `-a` and `-b` may be `-` (stdin) per invocation.
- With `-s` or `-S`, B intervals lacking a strand are skipped.
