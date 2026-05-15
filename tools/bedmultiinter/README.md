# bedmultiinter - Multi-way intersection across N BED files

A Go re-implementation of `bedtools multiinter` (a.k.a.
`multiIntersectBed`). Walks an event-driven sweep across N BED files and
emits one row per region where the active set of contributing files
changes.

## Features

- N input files (variadic `-i`), per-chromosome event sweep.
- `-empty` mode emits 0-count regions at chromosome heads, tails, and
  in any internal gaps (requires `-g` for chromosome sizes).
- `-cluster` collapses adjacent same-active-set rows.
- `-header` emits a column-header row.
- `-names` overrides the per-file labels (default: 1-based numeric
  indices in the `list` column, raw filenames in the header).
- `-filler` overrides the "absent" cell value (default `0`; `N/A` is
  another common choice).
- Within-file overlapping intervals are merged before sweeping, matching
  upstream's `GetNextMergedBed` behaviour.
- Pure Go, no third-party dependencies. Transparent gzip/BGZF input and
  `-` for stdin.

## Build

```bash
go build ./tools/bedmultiinter/cmd/bedmultiinter
```

## Usage

```bash
bedmultiinter -i <FILE1> <FILE2> [<FILE3> ...] [options]
```

### Options

- `-i FILE..`            Input BED files (variadic; >=2; `-` for stdin).
- `-names CSV`           Comma-separated labels (one per `-i` file).
- `-empty`               Emit 0-count regions; requires `-g`.
- `-g FILE`              Chrom-sizes genome file (required with `-empty`).
- `-cluster`             Collapse adjacent same-set rows.
- `-header`              Emit a column-header row.
- `-filler TEXT`         Indicator for absent cells (default `0`).
- `-o, --output FILE`    Output file (default: stdout).
- `-h, --help`           Show help.
- `-v, --version`        Show version.

### Output

```text
chrom  start  end  num  list  <per-file 0/1>...
```

- `num` is the count of files contributing to the region.
- `list` is a comma-separated list of file labels (or 1-based indices
  if `-names` is omitted). `none` is used when `num=0`.
- The trailing N columns are 0/1 indicators (or `-filler`/`1`).

## Examples

```bash
# Default mode (matches upstream's --examples heredoc).
bedmultiinter -i a.bed b.bed c.bed

# With header and labels.
bedmultiinter -i a.bed b.bed c.bed -header -names A,B,C

# Frame the output with leading/trailing/internal gap rows.
bedmultiinter -i a.bed b.bed c.bed -empty -g sizes.txt -header -names A,B,C
```

## Parity with upstream

Upstream ships no `multiinter/` subdirectory under
`reference_code/bedtools/test/`. The parity fixtures here are byte-for-byte
copies of the example documented in
`reference_code/bedtools/src/multiIntersectBed/multiIntersectBedMain.cpp::multiintersect_examples()`,
with the expected outputs lifted from the same heredoc. Three parity
tests cover the default mode, `-header -names`, and
`-empty -g -header -names`; all three pass byte-for-byte against the
upstream-documented example.

## Tests

```bash
go test -race -cover ./tools/bedmultiinter/...
```

Coverage: ~86.9% on `pkg/bedmultiinter` (race + cover, 2026-05-15).

## Limitations

- VCF / GFF input not implemented (upstream accepts these via its
  `BedFile` autodetect; this port reads BED only).
- Input files are assumed sorted by chrom/start; out-of-order records
  on a chromosome are tolerated only because each file's intervals are
  re-sorted and merged before the sweep. Strict upstream behaviour
  errors on unsorted input.
