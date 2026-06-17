# bedsummary - Per-chromosome interval summary statistics

A Go re-implementation of `bedtools summary`. Reports summary statistics for
the intervals in a BED/GFF/VCF file, per chromosome, against a genome file.

## Description

For every chromosome listed in the genome (`-g`) file — in genome-file order —
`bedsummary` reports:

- the chromosome length (from `-g`),
- the number of intervals on it,
- the total interval base pairs,
- the chromosome's fraction of the whole genome,
- its fraction of all intervals and of all interval bp,
- the minimum, maximum, and mean interval length.

A final `all` row aggregates over the entire input. The output is a 10-column
TSV whose header, column set, genome-file ordering, fixed 9-decimal precision,
and per-data-row trailing tab match upstream `bedtools summary` byte-for-byte.

## Installation

```bash
go build ./tools/bedsummary/cmd/bedsummary
```

## Usage

```text
bedsummary -i FILE.bed -g GENOME [options]
```

## Options

- `-i, --input FILE` - Input BED/GFF/VCF file (required, `-` for stdin)
- `-g, --genome FILE` - Genome (chrom-sizes) file (required)
- `-o, --output FILE` - Output file (default: stdout)
- `--no-header` - Suppress the column-header line
- `-h, --help` - Show help message
- `-v, --version` - Show version (`1.0.0`)

## Output (TSV)

```text
chrom  chrom_length  num_ivls  total_ivl_bp  chrom_frac_genome
frac_all_ivls  frac_all_bp  min  max  mean
```

Example:

```text
chrom	chrom_length	num_ivls	total_ivl_bp	chrom_frac_genome	frac_all_ivls	frac_all_bp	min	max	mean
chr1	1000	3	260	0.263157895	0.600000000	0.339869281	10	150	86.666666667	
chr2	500	1	10	0.131578947	0.200000000	0.013071895	10	10	10.000000000	
chr3	2000	1	495	0.526315789	0.200000000	0.647058824	495	495	495.000000000	
chr4	300	0	0	0.078947368	0.000000000	0.000000000	-1	-1	-1
all	3800	5	765	1.0	1.0	1.0	10	495	153.000000000
```

## Notes

- A genome file is **required** (matching upstream). It is tab/space delimited:
  `<chromName><TAB><chromSize>`. The chromosome order in this file determines
  the output row order.
- Chromosomes present in the genome file but with no intervals are reported with
  `0` counts and `-1` for min/max/mean.
- An interval whose chromosome is absent from the genome file is a hard error
  (matching upstream).
- Coordinates are 0-based half-open `[start, end)`.

## Parity

`live_parity_test.go` runs the real vendored `bedtools summary` binary and
asserts byte-for-byte equality on a small multi-chromosome input (including an
interval-free chromosome) and a 4000-interval, 5-chromosome dataset. The
previous port emitted the wrong column set (no `chrom_length` / fraction
columns) and did not accept `-g`; the rewrite reproduces the exact upstream
report.

## Tests

```bash
go test ./tools/bedsummary/...
```
