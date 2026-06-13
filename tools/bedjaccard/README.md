# bedjaccard - Jaccard similarity of two sorted BED files

A Go re-implementation of `bedtools jaccard`. Reads two BED files A and B,
both pre-sorted by `(chrom, start)`, and writes a one-line summary of the
intersection, union, Jaccard index and overlap count.

## Features

- Single streaming sweep over the inputs (memory is independent of file size)
- Strand-aware overlap (`-s` same-strand, `-S <+|->` single-strand filter)
- Overlap-fraction thresholds (`-f` for A, `-F` for B)
- Errors out on unsorted input rather than producing wrong answers
- Transparent gzip support and `-` for stdin/stdout
- Pure Go, no third-party dependencies

## Build

```bash
go build ./tools/bedjaccard/cmd/bedjaccard
```

## Usage

```bash
bedjaccard -a <fileA.bed> -b <fileB.bed> [options]
```

### Options

- `-a, --a FILE` First sorted BED file (required; `-` for stdin; `.gz` ok)
- `-b, --b FILE` Second sorted BED file (required)
- `--output FILE` Output file (default: stdout)
- `-s, --strand` Same-strand overlaps only (BED6 strand column required)
- `-S <+|->` Restrict both inputs to the given strand before overlap
- `-f FRACTION` Require >= FRACTION of A overlapped by B (0..1)
- `-F FRACTION` Require >= FRACTION of B overlapped by A (0..1)
- `-h, --help` Show help
- `-v, --version` Show version

`-s` and `-S` are mutually exclusive.

## Output

Tab-separated, two lines:

```text
intersection<TAB>union<TAB>jaccard<TAB>n_intersections
<int>           <int>      <float>     <int>
```

Where:

- `intersection` = total bases covered by both A and B
- `union` = `|A| + |B| - intersection`
- `jaccard` = `intersection / union` (`0` if union is 0)
- `n_intersections` = number of overlapping interval pairs

## Examples

```bash
# Plain Jaccard
bedjaccard -a peaksA.bed -b peaksB.bed

# Same strand only
bedjaccard -a a.bed -b b.bed -s

# Forward strand only (both inputs restricted to '+')
bedjaccard -a a.bed -b b.bed -S +

# Require >= 50% of A to be overlapped by B for the pair to count
bedjaccard -a a.bed -b b.bed -f 0.5

# Both files gzipped, output to a file
bedjaccard -a a.bed.gz -b b.bed.gz --output jaccard.tsv
```

## Algorithm

Both inputs are scanned once in `(chrom, start)` order. The tool keeps an
"active window" of B intervals that may still overlap subsequent A records,
appending newly read B records and dropping those whose end is before the
current A's start. Each new A is compared against the active window only.

## Testing

```bash
go test ./tools/bedjaccard/...
go test -cover ./tools/bedjaccard/pkg/bedjaccard
```

## License

Apache License 2.0.
