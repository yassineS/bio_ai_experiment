# bedslop - Extend BED intervals

A Go re-implementation of `bedtools slop`. Reads a BED file, extends each
interval by a fixed number of bases (or a fraction of the interval length),
clips the result to chromosome boundaries, and writes the surviving records to
stdout.

## Features

- Symmetric extension (`-b N`) or asymmetric (`-l L -r R`).
- Fractional mode (`--pct`) treats the extension as a fraction of the interval
  length.
- Strand-aware mode (`-s`) swaps left/right semantics for `-` strand records,
  so "left" / "right" are interpreted relative to the transcribed strand.
- Clips results to `[0, chromSize]` using a chrom-sizes file.
- Drops intervals that shrink to non-positive length and prints a warning for
  each (to stderr) including the original input line.
- Preserves every input column (BED3, BED6, BED12, or extras).
- Built-in gzip support (`.gz` inputs/outputs).

## Build

```bash
go build ./tools/bedslop/cmd/bedslop
```

## Usage

```bash
bedslop [options] -g <genome.sizes>
```

## Options

| Option | Description |
|--------|-------------|
| `-i, --input FILE` | Input BED file (`-` for stdin, default stdin) |
| `-o, --output FILE` | Output BED file (`-` for stdout, default stdout) |
| `-g, --genome FILE` | Chrom-sizes file (`chrom<TAB>size` per line). Required. |
| `-b NUM` | Extend by NUM bases (or fraction) on both sides. Mutually exclusive with `-l`/`-r`. |
| `-l NUM` | Extend by NUM on the left side; must be paired with `-r`. |
| `-r NUM` | Extend by NUM on the right side; must be paired with `-l`. |
| `-s, --strand` | Respect strand: swap left/right on `-` strand records. |
| `-pct, --percentage` | Treat `-b`/`-l`/`-r` as fractions (0..1) of interval length. |
| `-h, --help` | Show help and exit. |
| `-v, --version` | Show version and exit. |

A negative value shrinks the interval. Intervals clipped to zero or negative
length are dropped with a warning to stderr.

## Examples

Extend every interval by 50bp on each side:

```bash
bedslop -i input.bed -g hg38.sizes -b 50 > extended.bed
```

Asymmetric extension:

```bash
bedslop -i input.bed -g hg38.sizes -l 100 -r 25 > extended.bed
```

Strand-aware (upstream = 5' of the gene):

```bash
bedslop -i input.bed -g hg38.sizes -l 100 -r 25 -s > extended.bed
```

Extend by 25% of the interval length on each side:

```bash
bedslop -i input.bed -g hg38.sizes -b 0.25 --pct > extended.bed
```

## Parity notes

- `-b`, `-l`, `-r` accept fractional values and integer values (the latter is
  the common case).
- Records on chromosomes not present in the genome file are dropped with a
  stderr warning. This matches upstream `bedtools slop`.
- Comment lines (`#...`), `track` lines, and `browser` lines are filtered out
  of the output.
- `--pct` rounds the extension to the nearest integer (round-half-away-from-
  zero), matching upstream's integer-output behaviour.

## Testing

```bash
go test ./tools/bedslop/...
go test -cover ./tools/bedslop/pkg/bedslop
```
