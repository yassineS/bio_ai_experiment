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
| `-a` | `--input-a` | A intervals (required; `-` for stdin). BED or BAM/SAM (auto-detected). |
| `-b` | `--input-b` | B intervals (required). BED or BAM/SAM (auto-detected). |
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
|      | `--split`   | Treat blocked `-b` records (BED12 lines or spliced BAM alignments) as their blocks (exon-aware). A blocked `-a` record under `--split` is rejected. |
| `-h` | `--help`    | Show help. |
| `-v` | `--version` | Show version. |

Reads transparent gzip / BGZF via `pkg/htsgo/iohelper`. Stdin / `-` is
supported. Both `-a` and `-b` accept BAM/SAM input: the stream is sniffed
(BAM/BGZF magic or a leading `@` SAM header) and each mapped alignment is
converted to a BED12 record via `pkg/htsgo/alnbed` (CIGAR blocks become BED12
blocks), so `--split` composes for free on the `-b` side. This mirrors
upstream's `-abam` / `-b <bam>` modes.

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

- BAM/SAM input is supported on both `-a` and `-b` (auto-detected); the
  `-b <bam>` `-split` cases (upstream coverage.t10–t13) pass byte-for-byte.
  CRAM input is not yet supported.
- A blocked query (`-a`) under `--split` — a BED12 line or a spliced BAM
  alignment — is rejected with a clear error rather than producing a wrong
  answer. Upstream splits the query into its blocks; that path is not yet
  ported (upstream coverage.t1's BAM `-a` echo also differs only by a
  cosmetic trailing comma in the BED12 block lists).
- `-mean` prints with native float64 precision; upstream uses float32 and
  emits noise like `1.3200001`. We emit `1.32` instead. Semantic
  equivalence is covered by unit tests; the byte-parity test is `t.Skip`'d.
- `-sorted` (sorted-stream fast path) is accepted as a no-op since our
  default is already a single-pass interval-tree query.

See [`../PARITY_VALIDATION.md`](../PARITY_VALIDATION.md) for the validated
parity matrix and [`../../docs/PARITY_ROADMAP.md`](../../docs/PARITY_ROADMAP.md)
for the gap list.
