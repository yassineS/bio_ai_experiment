# bedmap

Pure-Go reimplementation of `bedtools map`. For each interval in A, find
overlapping intervals in B, take values from one or more B columns, apply an
aggregation op (sum, mean, collapse, ...), and append the result to A.

## Usage

```bash
bedmap -a A.bed -b B.bed                # default: sum of B's column 5
bedmap -a A.bed -b B.bed -c 4 -o collapse
bedmap -a A.bed -b B.bed -c 5 -o min,max
```

## Flags

| Short | Long | Notes |
|---|---|---|
| `-a` | `--input-a` | A intervals (required; `-` for stdin). |
| `-b` | `--input-b` | B intervals (required). |
|      | `--output` | Output file (default stdout). |
| `-c` | `--columns` | Comma-separated 1-based B column indices to aggregate (default `5`). |
| `-o` | `--operations` | Comma list of aggregation ops; one per column or one shared (default `sum`). Supported: `sum`, `min`, `max`, `mean`, `median`, `count`, `count_distinct`, `distinct`, `collapse`, `first`, `last`, `mode`, `antimode`. |
|      | `--null` | Placeholder for "no overlap" (default `.`). |
|      | `--delim` | Separator for `collapse` / `distinct` (default `,`). |
| `-s` | `--strand` | Same-strand only. |
| `-S` | `--opposite` | Opposite-strand only. |
| `-f` | `--fraction-a` | Min fraction of A that must overlap a B record. |
| `-F` | `--fraction-b` | Min fraction of B that must overlap A. |
| `-r` | `--reciprocal` | Require both `-f` and `-F`. |
| `-h` | `--help` | Show help. |
| `-v` | `--version` | Show version. |

Reads transparent gzip / BGZF via `pkg/htsgo/iohelper`. Stdin / `-` is
supported.

## Op semantics

`bedmap` reuses `bedmerge.ApplyOp` so the supported op vocabulary is exactly
the same as bedmerge / bedgroupby. See `tools/bedmerge/README.md` for the per-op
semantics.

When there are no overlapping B records:

- `count` and `count_distinct` always emit `0`.
- Every other op emits the `--null` placeholder (default `.`).

## Deviations from upstream

- GFF input is supported on the `-b` database (auto-detected: 1-based
  start/end in columns 4/5; `-c` extracts the literal GFF column). BAM/CRAM/VCF
  input is not yet supported. Tracked in `docs/PARITY_ROADMAP.md#bedtools`.
- `absmin`, `absmax`, `stdev`, `sstdev`, `cat`, `cat_uniq` not yet
  supported (same gap as `bedmerge` / `bedgroupby` — once they're added to
  `bedmerge.ApplyOp` they will be picked up here automatically).
- `-sorted` (sorted-stream fast path) is accepted as a no-op.

See [`../PARITY_VALIDATION.md`](../PARITY_VALIDATION.md) for the validated
parity matrix.
