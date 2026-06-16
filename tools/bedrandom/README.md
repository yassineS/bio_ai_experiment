# bedrandom - Generate random intervals across a genome

A Go re-implementation of `bedtools random` (the upstream `randomBed` tool).
Generates a requested number of fixed-length intervals placed uniformly at
random across a genome and writes them as BED6 records to stdout.

## Features

- Generates `-n` intervals, each of length `-l`, placed uniformly at random
  across the concatenated genome.
- Emits BED6: `chrom  start  end  index  length  strand`, where `index` counts
  `1..n`, the score column is the interval length, and `strand` is a random
  `+` or `-`.
- Rejects (and redraws) any placement whose end runs past the chromosome end,
  exactly as upstream does.
- With `-seed` the output is fully reproducible and **byte-for-byte identical**
  to upstream `bedtools random` (ports its `std::mt19937_64` engine and
  `rand_range` rejection-sampling bound).
- Built-in gzip support for the genome file and output (`.gz`).

## Build

```bash
go build ./tools/bedrandom/cmd/bedrandom
```

## Usage

```bash
bedrandom [options] -g <genome>
```

## Options

| Option | Description |
|--------|-------------|
| `-g, --genome FILE` | Genome (chrom-sizes) file: `chrom<TAB>size` per line. Required. |
| `-l, --length NUM` | Length of each interval (default: 100). |
| `-n, --number NUM` | Number of intervals to generate (default: 1000000). |
| `-seed NUM` | Seed for the RNG. If unset, a seed is derived from time+pid (output is then non-reproducible). |
| `-o, --output FILE` | Output file (`-` for stdout, default: stdout). |
| `-h, --help` | Show help and exit. |
| `-v, --version` | Show version and exit. |

## Examples

100 intervals of length 1000, reproducible:

```bash
bedrandom -g hg38.sizes -l 1000 -n 100 -seed 42 > random.bed
```

Default behaviour (1,000,000 intervals of length 100):

```bash
bedrandom -g hg38.sizes > random.bed
```

## Parity notes

- The placement algorithm mirrors upstream `randomBed.cpp` exactly: each
  interval draws `rand_range(genomeSize)` for a genome offset, projects it back
  to a `(chrom, start)` via the genome's file-order cumulative start offsets,
  redraws while `start + length` exceeds the chromosome size, then draws
  `rand_range(2)` for the strand. Reproducing this draw order and count is what
  makes the output byte-exact.
- The RNG is upstream's default (non-`USE_RAND`) `std::mt19937_64` engine. The
  seed is the integer passed to `-seed`; an unset seed becomes
  `time(0) + getpid()` (non-reproducible), matching upstream.
- Genome-file parsing matches upstream `GenomeFile`: file order is preserved,
  blank lines and lines whose first field begins with `#` are skipped, and a
  non-numeric size column is silently skipped.

## Testing

```bash
go test ./tools/bedrandom/...
go test -cover ./tools/bedrandom/pkg/bedrandom
```

The package includes a live-upstream parity test
(`upstream_parity_test.go`) that builds the real `bedtools` binary from the
vendored submodule and asserts byte-for-byte equality across several
`(genome, length, number, seed)` tuples.
