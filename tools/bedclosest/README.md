# bedclosest - Find the closest B interval for each A interval

A Go re-implementation of `bedtools closest`.

## Description

For each interval in A (sorted), find the closest interval in B (also sorted)
on the same chromosome and report A's columns, B's columns, and the signed
distance. Distance is 0 when A and B overlap. For tied distances the default
is to emit one row per tied B (`-t all`).

Both inputs MUST be sorted on `(chrom, start)`. bedclosest errors out clearly
when they are not.

Unlike upstream `bedtools closest`, the distance column is printed BY DEFAULT.
Use `--distance=false` to suppress it.

## Installation

```bash
go build ./tools/bedclosest/cmd/bedclosest
```

## Usage

```text
bedclosest -a <fileA.bed> -b <fileB.bed> [options]
```

## Options

- `-a, --a FILE` - Input BED file A (sorted; use `-` for stdin)
- `-b, --b FILE` - Input BED file B (sorted; use `-` for stdin)
- `-o, --output FILE` - Output BED file (`-` for stdout, default: stdout)
- `-d, --distance` - Print signed distance column (default: `true`; pass
  `--distance=false` to suppress)
- `-D MODE` - Strandedness of the distance sign:
  - `ref` (default): downstream is positive on the reference
  - `a`: relative to A's strand (BED6 col 6); flips on `-`-strand A
  - `b`: relative to B's strand
- `-N` - Require strict overlap; non-overlapping B intervals are treated as
  infinite (skipped)
- `-s` - Require the closest B to be on the SAME strand as A (BED6 col 6).
  Non-matching B intervals are skipped from candidate consideration.
- `-S` - Require the closest B to be on the OPPOSITE strand to A. Mutually
  exclusive with `-s`.
- `-t MODE` - Tie-break among equally-close B's:
  - `all` (default) - emit one row per tied B in B's input order
  - `first` - emit only the first tied B
  - `last` - emit only the last tied B
- `-h, --help` - Show help message
- `-v, --version` - Show version (`1.0.0`)

## Examples

```bash
# Closest peak for each gene
bedclosest -a genes.sorted.bed -b peaks.sorted.bed > out.bed

# Suppress the distance column
bedclosest -a a.bed -b b.bed --distance=false > out.bed

# Only report when A overlaps a B
bedclosest -a a.bed -b b.bed -N > out.bed

# Single hit per A (first in B input order on ties)
bedclosest -a a.bed -b b.bed -t first > out.bed

# Closest B on the same strand as A (skips opposite-strand B's)
bedclosest -a a.bed -b b.bed -s > out.bed

# Closest B on the opposite strand to A
bedclosest -a a.bed -b b.bed -S > out.bed
```

## Format

- Input: BED (tab-delimited, minimum 3 columns), sorted on `(chrom, start)`.
  `.gz` is supported.
- Output: A's columns, then B's columns, then signed distance (when `-d`).
- When A's chromosome has no B records, a sentinel B of
  `.\t-1\t-1` with distance `-1` is emitted (unless `-N` is set, in which case
  the A line is omitted).

## Algorithm

Records are read into memory and bucketed by chromosome. Per-chromosome,
bedclosest binary-searches A's start in B and then walks outward, pruning via
the running max-end of B (left side) and B.start (right side). This bounds the
work per A to O(log N + K) where K is the number of B's within the best
distance found so far.
