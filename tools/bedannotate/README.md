# bedannotate - Annotate A with overlap stats from N BED files

A Go re-implementation of `bedtools annotate`. For each interval in a
primary BED file (`-i`), it annotates with overlap statistics drawn from
one or more secondary BED files (`--files`).

## Features

- Default mode emits the fraction of A covered by each B (per-B float
  column, `%f`).
- `--counts` emits the count of overlapping records per B.
- `--both` interleaves count + fraction per B (2N columns total).
- Strand filters: `-s` (same-strand only), `-S` (opposite-strand only).
- Optional `--names` header (or basenames are used by default).
- Interval-tree based overlap (`pkg/bioformats/bed.IntervalTree`) — one
  tree per B file built once, then A is streamed.
- Pure Go, no third-party dependencies.
- Transparent gzip/BGZF input on every input and `-` for stdin on `-i`.

## Build

```bash
go build ./tools/bedannotate/cmd/bedannotate
```

## Usage

```bash
bedannotate -i <A.bed> --files <B1.bed> [<B2.bed> ...] [options]
```

### Options

- `-i, --input FILE`     A intervals (required, `-` for stdin)
- `--files FILE..`       One or more B files (variadic; stops at next flag)
- `--names N1,N2,..`     Comma-separated header labels (default: basenames)
- `--counts`             Emit per-B count of overlapping records
- `--both`               Emit count and coverage fraction per B (interleaved)
- `-s, --strand`         Same-strand overlaps only
- `-S, --opposite`       Opposite-strand overlaps only
- `-o, --output FILE`    Output file (default: stdout)
- `-h, --help`           Show help
- `-v, --version`        Show version

### Output

```text
chrom <TAB> start <TAB> end <TAB> ...A columns... <TAB> <col_per_B>...
```

With `--names`, a leading `#` header line precedes the data:

```text
# <TAB> exons <TAB> introns
```

With `--both`, the header splits each name into `<name>_cnt` and
`<name>_pct`.

## Deviations from upstream

- Upstream's leading-header padding aligns the `#` over A's bedType
  column; we emit a single `#` followed by the per-B labels. The
  parity test accepts the data rows as the source of truth and treats
  the header padding as a cosmetic deviation.
- Variadic `--files` is parsed by collecting positional values after
  `-files` / `--files` up to the next `-`-prefixed token. This mirrors
  upstream's parser.

## Tests

```bash
go test ./tools/bedannotate/...
```

The parity fixtures under `testdata/parity/` are hand-computed. Upstream
ships no `annotate/` test directory under
`reference_code/bedtools/test/`, so the expected outputs were derived
against the upstream algorithm in
`reference_code/bedtools/src/annotateBed/annotateBed.cpp`.

## Performance

`O(sum(|B_i|))` working set (per-B interval tree) and `O(|A| * log|B_i| + k)`
query time, where `k` is overlap count. Suitable for thousands-to-millions
of A intervals against megabyte-scale B files.
