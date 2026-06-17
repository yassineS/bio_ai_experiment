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
| `-f` | `--fraction-a` | Min fraction of A that must overlap a B record before B contributes. Ignored under `--split` (see below). |
| `-F` | `--fraction-b` | Min fraction of B that must overlap A. Ignored under `--split`. |
| `-r` | `--reciprocal` | Require the `-f` fraction reciprocally (A AND B; the B threshold equals `-f`). Ignored under `--split`. |
| `-e` | `--either`  | Require the minimum fraction for A **OR** B instead of the default AND across `-f`/`-F`. Ignored under `--split`. |
|      | `--split`   | Treat blocked records (BED12 lines or spliced BAM alignments) as their blocks (exon-aware), on both the `-b` database side and a blocked `-a` query. |
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

- BAM/SAM input is supported on both `-a` and `-b` (auto-detected). The
  BAM/BED12 `-a` echo is byte-for-byte with upstream (full 12 columns,
  trailing-comma block lists — coverage.t1 passes), and the `-b <bam>`
  `-split` cases (coverage.t10–t13) pass byte-for-byte. CRAM input is not yet
  supported.
- A blocked query (`-a`) under `--split` — a BED12 line or a spliced (`N`-CIGAR)
  BAM alignment — is fully supported and byte-validated against upstream across
  every mode (default / `--counts` / `--depth` / `--hist` / `--mean`, plus
  `-s`/`-S`). With `--split`, coverage is computed only over the query's
  sub-blocks (introns/gaps are excluded), while the reported length-of-A and the
  per-base depth vector still span the record's full `[start,end)` — intronic
  bases sit at depth 0 — matching upstream `coverageFile.cpp`
  (`_queryLen = endPos - startPos`). A single B feature that straddles an intron
  and overlaps two query blocks is counted once per block it touches, mirroring
  upstream's `findBlockedOverlaps`/`_hitCount`. See
  `TestUpstreamParity_SplitBlockedQuery{BED12,BAM}` (live upstream oracle) and
  `split_unit_test.go` (binary-free).
- Under `--split`, the overlap-fraction thresholds `-f` / `-F` (and therefore
  `-r` / `-e`) are **not applied** — every B feature overlapping any A block is
  counted, whatever the requested fractions. This matches upstream exactly: its
  blocked path (`coverageFile.cpp::checkSplits`) keeps the always-populated
  `BlockMgr` *overlapSet* rather than the fraction-filtered *resultSet* that the
  plain `intersect` path uses, so `-f`/`-F` have no effect on the count under
  `-split`. Verified against bedtools 2.31.1 (even `-f 1.0` / `-F 1.0` / `-r` /
  `-e` leave the count unchanged). The non-`--split` `-f`/`-F`/`-r`/`-e`
  behaviour is unaffected and still filters per record. See
  `TestUpstreamParity_FractionUnderSplitIgnored` and
  `TestUnitSplitSuppressesFractionFilter`.
- The BED12 `blockSizes` / `blockStarts` columns of an `-a` record are echoed
  **verbatim** — the optional trailing comma is preserved if present and omitted
  if absent, exactly as upstream re-emits them. (Earlier the port normalised the
  lists by always appending a trailing comma.) The `bed.Reader` now retains the
  raw column text (`RawBlockSizes`/`RawBlockStarts`); BAM-derived records carry
  no raw text and fall back to upstream's trailing-comma form (`50,50,`). See
  `TestUpstreamParity_VerbatimBED12BlockEcho` and `TestUnitVerbatimBlockEcho`.
- `-mean` reproduces upstream's float32-accumulated output (7 decimals,
  including float32 rounding noise such as `1.3200001`) — coverage.t6 passes.
- `-sorted` (sorted-stream fast path) is accepted as a no-op since our
  default is already a single-pass interval-tree query.

See [`../PARITY_VALIDATION.md`](../PARITY_VALIDATION.md) for the validated
parity matrix and [`../../docs/PARITY_ROADMAP.md`](../../docs/PARITY_ROADMAP.md)
for the gap list.
