# bedgroupby

Pure-Go reimplementation of `bedtools groupby`. Groups records by one or
more columns and applies column aggregation ops to the remaining columns.

## Usage

```bash
bedgroupby -g 1,2,3 -c 5 -o sum input.bed
```

## Flags

| Short | Long | Notes |
|---|---|---|
| `-g` | `--group-by` | Comma list of 1-based column indices to group by (default `1,2,3`). |
| `-c` | `--columns` | Comma list of 1-based column indices to apply ops to. |
| `-o` | `--operations` | Comma list of ops (parallel to `-c`). Supported ops: `sum`, `min`, `max`, `mean`, `median`, `count`, `count_distinct`, `distinct`, `collapse`, `first`, `last`, `mode`, `antimode`. Ops planned but not yet wired: `stdev`, `sstdev`, `absmin`, `absmax`, `cat`, `cat_uniq` — see `docs/PARITY_ROADMAP.md#bedtools`. |
| `-i` | `--input` | Input file (default stdin). `-` for stdin. |
| `-o` | `--output` | Output file (default stdout). `-` for stdout. |
| | `--in-header` | Treat first line as a header (passed through verbatim). |
| | `--full` | Emit every column of the first record per group, plus the aggregates (matches upstream `-full`). |
| `-h` | `--help` | Show help. |
| `-v` | `--version` | Show version. |

Reads transparent gzip / BGZF via `pkg/htsgo/iohelper`. Stdin / `-` is
supported.

## Deviations from upstream

- `stdev`/`sstdev`/`absmin`/`absmax`/`cat`/`cat_uniq` not yet supported
  (tracked in `docs/PARITY_ROADMAP.md`).
- `-ignorecase` accepts mixed-case grouping but always preserves the input
  record's case verbatim in the output (matches upstream's behaviour for
  the common path; minor edge case with multiple per-row case mappings
  documented in the parity tests).

See [`../PARITY_VALIDATION.md`](../PARITY_VALIDATION.md) for the validated
parity matrix and [`../../docs/PARITY_ROADMAP.md`](../../docs/PARITY_ROADMAP.md)
for the gap list.
