# mosdepth — pure-Go reimplementation

`mosdepth` is a pure-Go port of brentp's
[mosdepth](https://github.com/brentp/mosdepth) — fast per-base and
per-region BAM depth-of-coverage.

The original is written in Nim, which is the install-pain story we keep
seeing in CI/conda environments. This port is a single static Go binary
sitting on top of our in-tree SAM/BAM, BGZF, and tabix libraries — no
htslib, no Nim toolchain.

`mosdepth` is pick #5 of the 2026 next-up list in
`analysis/tool_ranking_2026.md`. It is built on top of the SAM/BAM stack
that landed with `tools/samtools` (pick #3).

## Build

```bash
go build ./tools/mosdepth/cmd/mosdepth
go install ./tools/mosdepth/cmd/mosdepth
```

The binary has no third-party dependencies and no `cgo`.

## Usage

```bash
mosdepth [options] <prefix> <in.bam>
```

Outputs (always emitted unless a flag below suppresses one):

- `<prefix>.mosdepth.global.dist.txt` — cumulative coverage distribution
  per chromosome plus a `total` row. Format: `chrom\tdepth\tproportion`
  where `proportion` is the fraction of bases at depth ≥ `depth`.
- `<prefix>.mosdepth.summary.txt` — per-chromosome summary
  (`chrom\tlength\tbases\tmean\tmin\tmax`) plus a `total` row.
- `<prefix>.per-base.bed.gz` (with `.tbi`) — per-base depth as
  collapsed equal-depth BED runs. Omitted when `--by` is set or
  `--no-per-base` is passed.
- `<prefix>.regions.bed.gz` (with `.tbi`) — only when `--by` is set.
  Columns: `chrom\tstart\tend\t[region-name\t]mean-depth`.
- `<prefix>.thresholds.bed.gz` (with `.tbi`) — only when `-T/--thresholds`
  is set. Columns: `chrom\tstart\tend\tregion\tNXcount\t...` listing the
  number of bases at or above each integer threshold inside the region.

## Flags

| Short | Long | Description |
| --- | --- | --- |
| `-t` | `--threads INT` | Accepted; v1 is single-threaded. |
| `-b` | `--by FILE_OR_INT` | BED file of regions, or an integer window size in bases. |
| `-Q` | `--mapq INT` | Minimum MAPQ (default 0). |
| `-F` | `--flag INT` | Exclude reads with ANY of these flag bits (default `1796` = `0x704`: unmapped, secondary, QC-fail, duplicate, supplementary). |
| `-i` | `--include-flag INT` | Keep only reads with ALL of these flag bits. |
| `-x` | `--fast-mode` | Skip CIGAR walking; treat each read as covering POS..POS+ReferenceLength. ~3x faster, slightly inaccurate near indels. |
| `-n` | `--no-per-base` | Suppress the per-base output. |
| `-T` | `--thresholds LIST` | Comma list of integer thresholds (e.g. `1,5,10,30`). |
| `-c` | `--chrom STRING` | Restrict to one chromosome. |
| `-d` | `--d4` | Write the per-base depth track to `<prefix>.per-base.d4` in the dense D4 binary format instead of the bgzipped BED. |
| `-r` | `--read-groups LIST` | Comma list of allowed RG ids; prefix the first with `OPS:` to filter on the OPS aux tag instead. |
| `-l` | `--min-frag-len INT` | Minimum absolute TLEN to include. |
| `-u` | `--max-frag-len INT` | Maximum absolute TLEN to include. |
| `-h` | `--help` | Show help. |
| `-v` | `--version` | Show version. |

BED region input format (for `--by`): `chrom\tstart\tend\t[name]` — the
optional 4th column populates the per-region output as upstream does.

## Algorithm

mosdepth is single-pass. For each record we accumulate `+1/-1` events at
the start/end of every reference-consuming CIGAR run (`M`, `=`, `X`); in
fast mode the whole POS..POS+ReferenceLength is recorded as one run.
Events are sorted and swept across each chromosome to recover depth at
every base without materialising a per-base depth array. Each chromosome
is emitted before the next is processed, keeping memory bounded to one
chromosome's-worth of events.

## Deviations from upstream

mosdepth in Nim emits CSI indexes; our port emits TBI indexes built with
the in-tree `pkg/htsgo/tabix.Build`. Consumers that read either
format (e.g. `bcftools`, `tabix`) work transparently — the underlying
chunk/bin layout is identical. CSI emission is on the roadmap.

D4 output (`-d/--d4`) writes the per-base depth track to
`<prefix>.per-base.d4` in the real [D4](https://github.com/38/d4-format)
binary container, matching the upstream mosdepth binary **byte-for-byte**.
The file is a d4-framefile: the `d4\xdd\xdd` magic, a `.metadata` stream
(the JSON header with the chromosome list, the `SimpleRange{0,128}`
dictionary, and the denominator), a bit-width-packed `.ptab` primary
table (7 bits per base, all chromosomes concatenated in header order), a
`.stab` secondary-table sub-directory, and a `.index` sub-directory.
Depths ≥ 128 are clamped to the all-ones primary code (127), exactly as
upstream's d4 C-binding writer does for its per-base output (the
chromosome summary still reports the true maximum depth); the secondary
table is therefore empty. The on-disk size equals upstream's.

Parity is validated **against the real binary**, not a golden file:
`TestD4_UpstreamBinaryParity` downloads the official `mosdepth_d4`
release binary, runs it and our implementation on the same fixture BAM,
and asserts the two `.per-base.d4` files are byte-identical (and that our
reader decodes both to the same per-base depths). When `--d4` is set the
per-base BED is not written.

The runtime is single-threaded; the `-t/--threads` flag is accepted for
compatibility with existing pipelines. A future slice may parallelise
the per-chromosome event sweep.

**Overlap-pair detection** — upstream subtracts one copy of depth where
the two ends of a mate-paired fragment overlap on the reference. Our v1
engine doesn't implement this pairing pass, so our default-mode output
matches upstream's `--fast-mode` output rather than upstream's default.
Tracked at
[docs/UPSTREAM_BUGS.md#mosdepth-overlap-pair-detection](../../docs/UPSTREAM_BUGS.md#mosdepth-overlap-pair-detection).

## Testing

```bash
go test -race -cover ./tools/mosdepth/...
```

Coverage targets ≥85% on `pkg/mosdepth`. Tests cover:

- Per-base depth on a hand-computed mini-BAM (one CIGAR walk and one
  fast-mode walk).
- Per-region summary (BED-driven) — verifies the mean-depth column.
- Per-window summary (`--by 10`) over a small reference.
- Threshold proportions (`-T 1,2,5`).
- MAPQ filter (`-Q`), include-flag filter (`-i`), exclude-flag filter
  (`-F` default).
- `--no-per-base` suppression.
- `--chrom` restriction (only one chromosome appears in any output).
- Distribution file CDF — `total\t0` is always 1.00.
- Summary file — chr/total rows match hand-computed mean.
- `.tbi` indexes round-trip through `tabix.QueryBytes` so the produced
  files are tabix-readable end-to-end.
- `-d/--d4` D4 output is byte-identical to the upstream `mosdepth_d4`
  binary on the same BAM (`TestD4_UpstreamBinaryParity`), plus a
  writer/reader round-trip on a hand-built track.

## Status

- v1: per-base, per-region (BED), per-window, thresholds, distribution,
  summary, TBI indexes, byte-identical D4 per-base output.
- Roadmap: CSI output, multi-threaded chrom sweep, CRAM input (depends on
  the project's CRAM reader landing first).
