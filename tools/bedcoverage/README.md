# bedcoverage

Pure-Go reimplementation of `bedtools coverage`. For each interval in A,
report how many features in B overlap, the total bp covered, the length of A,
and the fraction of A covered.

## Usage

```bash
bedcoverage -a A.bed -b B.bed
```

The default output appends four columns to each A line:

```text
<count> <bp_covered> <length_A> <fraction>
```

## Flags

| Short | Long | Notes |
|---|---|---|
| `-a` | `--input-a` | A intervals (required; `-` for stdin). |
| `-b` | `--input-b` | B intervals (required). |
|      | `--output`  | Output file (default stdout). |
|      | `--counts`  | Append only the overlap count. |
| `-d` | `--depth`   | Emit one line per base in A: A + 1-based position within A + depth. |
|      | `--hist`    | Per-A depth histogram + an "all" footer aggregated across all A records. |
|      | `--mean` / `--median` / `--min` / `--max` / `--sum` | Collapse the per-base depth vector with the requested op; append a single value. |
| `-s` | `--strand`  | Same-strand only. |
| `-S` | `--opposite`| Opposite-strand only. |
| `-f` | `--fraction-a` | Min fraction of A that must overlap a B record before B contributes. |
| `-F` | `--fraction-b` | Min fraction of B that must overlap A. |
| `-r` | `--reciprocal` | Require both `-f` and `-F` (default behaviour is already AND). |
| `-h` | `--help`    | Show help. |
| `-v` | `--version` | Show version. |

Reads transparent gzip / BGZF via `pkg/bioformats/iohelper`. Stdin / `-` is
supported.

## Output modes (mutually exclusive)

- **default** — count, bp covered, len A, fraction. Matches upstream's
  default block layout (7 decimal-place fraction, e.g. `0.7600000`).
- `--counts` — only the count.
- `--depth` — one row per base, A + (1-based) pos + depth.
- `--hist` — per-A histogram, depth bucket + bp at that depth + len + fraction.
  An `all` footer aggregates across all A records.
- `--mean` / `--median` / `--min` / `--max` / `--sum` — single-column summary
  of the per-base depth vector.

## Deviations from upstream

- BAM/SAM/CRAM input is not supported (BED only); upstream's `-abam` /
  `-b` BAM modes are out of scope. Tracked in
  `docs/PARITY_ROADMAP.md#bedtools`.
- `-split` (BED12 block-aware coverage) not yet supported.
- `-mean` prints with native float64 precision; upstream uses float32 and
  emits noise like `1.3200001`. We emit `1.32` instead. Semantic
  equivalence is covered by unit tests; the byte-parity test is `t.Skip`'d.
- `-sorted` (sorted-stream fast path) is accepted as a no-op since our
  default is already a single-pass interval-tree query.

See [`../PARITY_VALIDATION.md`](../PARITY_VALIDATION.md) for the validated
parity matrix and [`../../docs/PARITY_ROADMAP.md`](../../docs/PARITY_ROADMAP.md)
for the gap list.
