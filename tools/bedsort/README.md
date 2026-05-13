# bedsort - BED Interval Sorter

A small, fast re-implementation of `bedtools sort` in Go. Reads a BED file (or
stdin), sorts the records according to the requested mode, and writes them
back out. The full input is preserved column-for-column.

## Features

- Default sort by `(chrom, start, end)` (lexicographic chromosome order).
- Alternative sort modes: by interval size, by chromosome-then-size, by
  chromosome-then-score.
- Optional chromosome ordering from a `.fai` or chrom-sizes file
  (`--faidx`/`-g`).
- Preserves every input column (BED3, BED6, BED12, or extra fields).
- Built-in gzip support (`.gz` inputs/outputs).
- Stable sort: rows that tie under the chosen ordering keep their input order.
- Pure Go, standard library only.

## Build

```bash
go build ./tools/bedsort/cmd/bedsort
```

## Usage

```bash
bedsort [options] [<input.bed>]
```

If neither `-i` nor a positional argument is given, bedsort reads from stdin.
Output defaults to stdout. Use `-` explicitly for stdin/stdout.

## Options

| Option | Description |
|--------|-------------|
| `-i, --input FILE` | Input BED file (`-` for stdin, default stdin) |
| `-o, --output FILE` | Output BED file (`-` for stdout, default stdout) |
| `--sizeA` | Sort by interval size ascending |
| `--sizeD` | Sort by interval size descending |
| `--chrThenSizeA` | Sort by chromosome, then by interval size ascending |
| `--chrThenSizeD` | Sort by chromosome, then by interval size descending |
| `--chrThenScoreA` | Sort by chromosome, then by score (column 5) ascending |
| `--chrThenScoreD` | Sort by chromosome, then by score (column 5) descending |
| `--faidx FILE` | Order chromosomes by their appearance in FILE (`.fai` or chrom-sizes) |
| `-g, --genome FILE` | Alias for `--faidx` |
| `-h, --help` | Show help and exit |
| `-v, --version` | Show version and exit |

The size/score modes are mutually exclusive; combining two raises an error.

## Examples

Default sort:

```bash
bedsort -i input.bed -o sorted.bed
```

Stdin/stdout via pipe:

```bash
cat input.bed | bedsort > sorted.bed
```

Sort by interval size, descending:

```bash
bedsort --sizeD input.bed
```

Sort using a custom chromosome order:

```bash
bedsort --faidx hg38.fa.fai input.bed > sorted.bed
```

## Parity notes

- The default ordering matches upstream `bedtools sort`: chromosomes sort
  lexicographically as strings (so `chr10` precedes `chr2`). Use `--faidx` to
  get a natural / genome-ordered chromosome sort instead.
- Comment lines (`#...`), `track` lines, and `browser` lines are stripped from
  the output, matching upstream behaviour.
- All input columns round-trip through the sort, including any extra columns
  past BED12.

## Testing

```bash
go test ./tools/bedsort/...
go test -cover ./tools/bedsort/pkg/bedsort
```
