# bedreldist - Relative-distance distribution between two BED files

A Go re-implementation of `bedtools reldist`. For each interval in A, finds
the two nearest B-interval midpoints that flank A's midpoint and reports
the relative distance:

```text
rel_dist = min(|m - left|, |m - right|) / (right - left)
```

where `m` is A's midpoint and `left`/`right` are the surrounding
B-midpoints. The value lies in `[0.0, 0.5]`.

## Features

- Histogram output with 0.01-wide bins (default), matching `bedtools reldist`.
- Per-A detail mode (`-detail` / `--detail`).
- Single linear pass over A; B is loaded into memory once and sorted by
  chromosome.
- Pure Go, no third-party dependencies.
- Transparent gzip/BGZF input and `-` for stdin.

## Build

```bash
go build ./tools/bedreldist/cmd/bedreldist
```

## Usage

```bash
bedreldist -a <fileA.bed> -b <fileB.bed> [options]
```

### Options

- `-a, --a FILE`  BED file A (queries; required; `-` for stdin; `.gz` ok)
- `-b, --b FILE`  BED file B (database; required)
- `--output FILE` Output file (default: stdout)
- `-detail, --detail` Emit per-A relative distances instead of the histogram
- `-h, --help`    Show help
- `-v, --version` Show version

### Output

Default histogram:

```text
reldist  count  total  fraction
0.00     <count>  <total>  <fraction>
...
0.50     <count>  <total>  <fraction>
```

Only non-empty bins are emitted, in ascending order. The fraction column
uses `%.3f` formatting; the reldist bin uses `%.2f`.

With `-detail`, each A interval is emitted as

```text
chrom  start  end[  other A fields]  <rel_dist>
```

where `<rel_dist>` is formatted with `%.3f`.

## Skipped / deviations from upstream

- The full upstream test suite uses multi-megabyte fixtures
  (`refseq.chr1.exons.bed.gz`, `aluY.chr1.bed.gz`, `gerp.chr1.bed.gz`) that
  this repo does not vendor; the algorithm is covered by smaller hand-built
  fixtures plus the upstream `issue_711` corner case.
- `-detail` reconstructs the BED prefix from parsed fields; for BED3-6
  inputs the output matches upstream verbatim, but tools that pass through
  custom trailing fields are not exercised in this port.

## Tests

```bash
go test ./tools/bedreldist/...
```
