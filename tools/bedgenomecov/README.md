# bedgenomecov - Coverage histogram and depth tracks for BED intervals

A Go re-implementation of `bedtools genomecov`. Reads BED (or bedGraph)
intervals plus a chromosome-sizes file and produces a coverage histogram, a
bedGraph track, or per-base depth.

## Features

- BED input or BAM/SAM input (`--ibam`, genome taken from the alignment header)
- Histogram output (default), bedGraph (`-bg`/`-bga`) and per-base depth (`-d`/`-dz`)
- Strand-aware counting (`--strand +|-`) and 5'/3' end-only counting (`-5`/`-3`)
- Histogram depth cap (`--max`) and depth multiplier (`--scale`)
- Optional UCSC trackline header for bedGraph output
- Transparent gzip support and `-` for stdin/stdout
- Pure Go, no third-party dependencies

## Build

```bash
go build ./tools/bedgenomecov/cmd/bedgenomecov
```

## Usage

```bash
bedgenomecov -i <intervals.bed> -g <chrom.sizes> [options]
```

### Options

- `-i, --input FILE` Input BED file (default: stdin; `-` for stdin; `.gz` ok)
- `--ibam FILE` Input BAM/SAM file; the genome is taken from its `@SQ` header
  (no `-g` needed). Each alignment covers its reference span, or its CIGAR
  blocks under `--split`.
- `-pc, --pair-coverage` Coverage of paired-end fragments (BAM only): each
  proper pair contributes one fragment `[pos, pos+TLEN)`.
- `-fs, --fragment-size N` Force an N-base fragment per read instead of its
  alignment length (BAM only), anchored at the read's 5' end.
- `-g, --genome FILE` Chromosome sizes file, required unless `--ibam` (`chrom<TAB>size`)
- `--output FILE` Output file (default: stdout; `.gz` ok)
- `-bg, --bedGraph` Emit non-zero runs of constant depth as bedGraph
- `-bga` Emit every run of constant depth (includes zero)
- `-d, --per-base` Per-base depth (1-based positions)
- `-dz, --per-base-nonzero` Per-base depth, skip positions with depth 0
- `--strand +|-` Count only intervals on the given strand
- `--max N` (alias `--max-depth`) Cap histogram depth at N
- `--scale FLOAT` Multiply every depth by FLOAT (default 1.0)
- `-5, --five-prime` Count only the 5'-most base of each interval
- `-3, --three-prime` Count only the 3'-most base of each interval
- `--split` Treat BED12 records as their blocks (exon-aware coverage)
- `--trackline` Prepend a UCSC `track` line to `-bg`/`-bga` output
- `--trackopts STR` Extra trackline attributes appended after `track`
- `-h, --help` Show help
- `-v, --version` Show version

The output-mode flags (`-bg`, `-bga`, `-d`, `-dz`) are mutually exclusive; the
default is the histogram.

## Output formats

### Histogram (default)

Tab-separated: `chrom<TAB>depth<TAB>n_bases<TAB>chrom_size<TAB>frac`. After the
per-chromosome rows, a `genome` row aggregates across all chromosomes.

### bedGraph (`-bg`/`-bga`)

`chrom<TAB>start<TAB>end<TAB>depth`, one row per run of constant depth. `-bga`
also emits zero-depth runs.

### Per-base (`-d`/`-dz`)

`chrom<TAB>position<TAB>depth`. Positions are 1-based (bedtools convention).

## Examples

```bash
# Histogram
bedgenomecov -i reads.bed -g hg38.sizes

# bedGraph track, non-zero only
bedgenomecov -i reads.bed -g hg38.sizes -bg > coverage.bedgraph

# bedGraph track with a UCSC trackline
bedgenomecov -i reads.bed -g hg38.sizes -bga \
  --trackline --trackopts 'name=cov color=80,80,80' > coverage.bedgraph

# Per-base depth, skipping zero-depth bases
bedgenomecov -i reads.bed -g hg38.sizes -dz > depth.tsv

# Only count + strand intervals
bedgenomecov -i reads.bed -g hg38.sizes -bg --strand=+

# Cap histogram depth at 100
bedgenomecov -i reads.bed -g hg38.sizes --max 100
```

## Algorithm

For each chromosome listed in the genome file, the tool allocates an integer
depth array sized to the chromosome length. Every BED interval (after optional
strand filtering and 5'/3' end-only selection) increments the array. The
chosen output mode then sweeps each array. This linear-memory choice matches
bedtools' typical footprint and is the right tradeoff for vertebrate-scale
inputs; for genomes with very long single contigs you may prefer a streaming
event-based approach.

## Testing

```bash
go test ./tools/bedgenomecov/...
go test -cover ./tools/bedgenomecov/pkg/bedgenomecov
```

## License

Apache License 2.0.
