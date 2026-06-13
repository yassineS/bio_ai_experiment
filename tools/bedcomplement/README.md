# bedcomplement - Genomic regions NOT covered by a BED file

A Go re-implementation of `bedtools complement`. For each chromosome listed in
a chrom-sizes file, emits the half-open gaps not covered by the input BED
intervals.

## Features

- Streams a sorted BED file and emits the complementary intervals.
- Output is BED3 (`chrom<TAB>start<TAB>end`).
- Chromosomes with no intervals produce a single full-length record.
- Validates input ordering on the fly and emits a clear error if the input is
  not sorted on `(chrom, start)`.
- Intervals on chromosomes absent from the chrom-sizes file are skipped, with
  a single stderr warning per such chromosome.
- Overlapping/adjacent intervals on the same chromosome are merged before
  computing the complement.
- Built-in gzip support (`.gz` inputs/outputs).

## Build

```bash
go build ./tools/bedcomplement/cmd/bedcomplement
```

## Usage

```bash
bedcomplement [options] -g <genome.sizes>
```

## Options

| Option | Description |
|--------|-------------|
| `-i, --input FILE` | Input sorted BED file (`-` for stdin, default stdin) |
| `-o, --output FILE` | Output BED file (`-` for stdout, default stdout) |
| `-g, --genome FILE` | Chrom-sizes file (`chrom<TAB>size`). Required. |
| `-L, --limit` | Only emit chromosomes that had records in the input. |
| `-h, --help` | Show help and exit |
| `-v, --version` | Show version and exit |

## Examples

Print all regions not covered by `genes.sorted.bed`:

```bash
bedcomplement -i genes.sorted.bed -g hg38.sizes > intergenic.bed
```

Pipe sorted BED through bedcomplement:

```bash
bedsort -i genes.bed | bedcomplement -g hg38.sizes > intergenic.bed
```

## Parity notes

- Input must be sorted on `(chrom, start)`. If not, bedcomplement aborts with
  an error.
- Chromosomes are emitted in the order they appear in the chrom-sizes file;
  any chromosomes present in the file but absent from the input are still
  emitted as a single `chrom<TAB>0<TAB>size` record.
- Comment, `track`, and `browser` lines are ignored.
- Output is always BED3.

## Testing

```bash
go test ./tools/bedcomplement/...
go test -cover ./tools/bedcomplement/pkg/bedcomplement
```
