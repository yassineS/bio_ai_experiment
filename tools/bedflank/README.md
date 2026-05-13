# bedflank - Emit the flanking regions of each BED interval

A Go re-implementation of `bedtools flank`.

## Description

For each interval `[s, e)` on a chromosome of size `C`, bedflank emits up to
two new intervals:

- Left ("upstream") flank: `[max(0, s-l), s)`
- Right ("downstream") flank: `[e, min(C, e+r))`

The original interval itself is NOT emitted (use `bedslop` for that). Empty
flanks are skipped. With `-s/--strand`, the flanks are interpreted relative to
the transcribed strand (`l` and `r` swap on `-`-strand entries).

## Installation

```bash
go build ./tools/bedflank/cmd/bedflank
```

## Usage

```text
bedflank -i <input.bed> -g <genome.sizes> [-b N | -l L -r R] [options]
```

## Options

- `-i, --input FILE` - Input BED file (`-` for stdin, default: stdin)
- `-o, --output FILE` - Output BED file (`-` for stdout, default: stdout)
- `-g, --genome FILE` - Chromosome sizes file (`chrom<TAB>size` per line;
  samtools `.fai` also accepted). Required.
- `-b NUM` - Symmetric flank size (mutually exclusive with `-l`/`-r`)
- `-l NUM` - Left flank size (requires `-r`)
- `-r NUM` - Right flank size (requires `-l`)
- `-s, --strand` - Respect strand: swap `l`/`r` on `-`-strand records
- `-pct, --percentage` - Treat `-b`/`-l`/`-r` as fractions (0..1) of the
  interval length
- `-h, --help` - Show help message
- `-v, --version` - Show version (`1.0.0`)

## Examples

```bash
# 50bp flank on each side
bedflank -i input.bed -g hg38.sizes -b 50 > flanks.bed

# 100bp upstream, 25bp downstream
bedflank -i input.bed -g hg38.sizes -l 100 -r 25 > flanks.bed

# Strand-aware promoter-like flanks
bedflank -i tss.bed -g hg38.sizes -l 1000 -r 100 -s > regions.bed

# 10% of each interval as flank
bedflank -i input.bed -g hg38.sizes -b 0.1 --pct > flanks.bed
```

## Format

- Input: BED (tab-delimited, minimum 3 columns); `.gz` files supported.
- Output: same column layout as input, with `start`/`end` set to the flank
  coordinates. Each input interval produces 0, 1, or 2 output rows.

## Notes

- Coordinates are 0-based half-open `[start, end)`.
- Records on chromosomes missing from the genome file are dropped with a
  warning on stderr.
