# bedshift - Shift BED/GFF/VCF features

A Go re-implementation of `bedtools shift` (aka shiftBed). Reads a BED/GFF/VCF
file, shifts each feature by a requested number of base pairs, clamps the result
to the chromosome bounds, and writes the records to stdout.

## Features

- Uniform shift (`-s N`) or strand-aware shift (`-p` for `+`, `-m` for `-`).
- Fractional mode (`-pct`) treats the shift as a fraction of each feature's
  length.
- Clamps the start to `[0, chromSize-1]` and the end to `[1, chromSize]` using a
  chrom-sizes file, exactly like upstream (including upstream's truncation of the
  shifted coordinate toward zero).
- `-header` prints the input's leading header lines before the results.
- Preserves every input column (BED3, BED6, BED12, or extras).
- Built-in gzip support (`.gz` inputs/outputs).

## Build

```bash
go build ./tools/bedshift/cmd/bedshift
```

## Usage

```bash
bedshift [options] -i <bed/gff/vcf> -g <genome> [-s <int> | (-p and -m)]
```

## Options

| Option | Description |
|--------|-------------|
| `-i, --input FILE` | Input file (`-` for stdin, default stdin). |
| `-o, --output FILE` | Output file (`-` for stdout, default stdout). |
| `-g, --genome FILE` | Chrom-sizes file (`chrom<TAB>size` per line; also `.fai`). Required. |
| `-s NUM` | Shift every feature by NUM bp. Mutually exclusive with `-p`/`-m`. |
| `-p NUM` | Shift `+` strand features by NUM bp. Requires `-m`. |
| `-m NUM` | Shift `-` strand features by NUM bp. Requires `-p`. |
| `-pct, --percentage` | Treat `-s`/`-p`/`-m` as a fraction of the feature length. |
| `-header` | Print the input header before the results. |
| `-h, --help` | Show help and exit. |
| `-v, --version` | Show version and exit. |

Either `-s` alone, or `-p` and `-m` together, must be supplied.

## Examples

Shift every feature forward 5 bp:

```bash
bedshift -i a.bed -g genome.txt -s 5
```

Shift `+` and `-` strand features by different amounts:

```bash
bedshift -i a.bed -g genome.txt -p 100 -m 50
```

Shift each feature by half its own length:

```bash
bedshift -i a.bed -g genome.txt -s 0.5 -pct
```

## Parity notes

- Shift amounts are parsed with 32-bit float precision, matching upstream's C
  `float` fields, and the shifted coordinate is truncated toward zero on
  assignment to the integer coordinate (as the C++ `double`->`int` cast does).
- A record whose chromosome is absent from the genome file reproduces upstream's
  behaviour exactly: `getChromSize` returns `-1`, which yields negative output
  coordinates. This is a faithful re-creation of upstream's edge-case output.
- Validated byte-for-byte against the upstream `bedtools/test/shift/`
  `test-shift.sh` suite (cases t1–t11, including the issue-807 fractional case
  and the huge-genome and over-int8 shift cases).

## Testing

```bash
go test ./tools/bedshift/...
go test -cover ./tools/bedshift/pkg/bedshift
```
