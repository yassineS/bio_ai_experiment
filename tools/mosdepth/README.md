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
- `<prefix>.per-base.bed.gz` (with `.csi`) — per-base depth as
  collapsed equal-depth BED runs. Omitted when `--by` is set or
  `--no-per-base` is passed.
- `<prefix>.regions.bed.gz` (with `.csi`) — only when `--by` is set.
  Columns: `chrom\tstart\tend\t[region-name\t]mean-depth`.
- `<prefix>.quantized.bed.gz` (with `.csi`) — only when `-q/--quantize`
  is set. Columns: `chrom\tstart\tend\tlabel`, one collapsed run per
  contiguous quantize bin.
- `<prefix>.thresholds.bed.gz` (with `.csi`) — only when `-T/--thresholds`
  is set. Columns: `chrom\tstart\tend\tregion\tNXcount\t...` listing the
  number of bases at or above each integer threshold inside the region.

## Flags

| Short | Long | Description |
| --- | --- | --- |
| `-t` | `--threads INT` | Number of BGZF/BAM decompression threads. Blocks are inflated across N goroutines and reassembled in order, so every output file is byte-identical for any thread count; it only affects throughput. `<2` runs single-threaded. |
| `-b` | `--by FILE_OR_INT` | BED file of regions, or an integer window size in bases. |
| `-Q` | `--mapq INT` | Minimum MAPQ (default 0). |
| `-F` | `--flag INT` | Exclude reads with ANY of these flag bits (default `1796` = `0x704`: unmapped, secondary, QC-fail, duplicate). Matches upstream; supplementary reads are NOT excluded by default. |
| `-i` | `--include-flag INT` | Keep only reads with ALL of these flag bits. |
| `-x` | `--fast-mode` | Skip CIGAR walking; treat each read as covering POS..POS+ReferenceLength. ~3x faster, slightly inaccurate near indels. |
| `-n` | `--no-per-base` | Suppress the per-base output. |
| `-T` | `--thresholds LIST` | Comma list of integer thresholds (e.g. `1,5,10,30`). |
| `-c` | `--chrom STRING` | Restrict to one chromosome. |
| — | `--d4` | Write the per-base depth track to `<prefix>.per-base.d4` in the dense D4 binary format instead of the bgzipped BED. Upstream has no short form; `-d` is a port-only alias. |
| `-R` | `--read-groups LIST` | Comma list of allowed RG ids; prefix the first with `OPS:` to filter on the OPS aux tag instead. (`-r` is a port-only lowercase alias.) |
| `-l` | `--min-frag-len INT` | Minimum absolute TLEN to include. |
| `-u` | `--max-frag-len INT` | Maximum absolute TLEN to include. |
| `-f` | `--fasta FILE` | FASTA reference for CRAM input. Accepted for parity; CRAM is not yet supported, so the value is ignored. |
| `-a` | `--fragment-mode` | Count coverage across the whole template (fragment) between properly-paired mates rather than the aligned reads only. Only read1 of a proper, non-supplementary pair contributes, covering `[min(read,mate) start, +\|TLEN\|)`. Byte-identical to upstream v0.3.14. Mutually exclusive with `-x/--fast-mode` (rejected, exit 2). |
| `-q` | `--quantize SEGS` | Bin per-base depth into the `:`-separated segments (e.g. `0:1:4:`) and write `<prefix>.quantized.bed.gz`. A leading `:` prepends `0`; a trailing `:` adds an open-ended top bin (`N:inf`). Labels default to `lo:hi` and can be overridden per bin with `MOSDEPTH_Q<i>` env vars. Depths outside the range leave a gap (no line). Byte-identical to upstream v0.3.14. |
| `-m` | `--use-median` | Upstream flag, not yet implemented in this port: supplying it is rejected (exit 2). |
| `-h` | `--help` | Show help. |
| `-v` | `--version` | Show version. |

Single-char short flags may be clustered docopt-style (`-nx` == `-n -x`)
and values may be concatenated (`-Q20` == `-Q 20`), matching upstream
mosdepth's docopt command-line parser. `-v/--version` is a port
convenience not present in upstream's docopt usage.

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

## Indexes

Like upstream mosdepth (which calls htslib's `tbx_index_build`), our port
emits a `.csi` index alongside each bgzipped BED output, built with the
in-tree `pkg/htsgo/tabix.BuildCSIFromDataFile`. The CSI uses `min_shift=14,
depth=5` — the htslib default — so it indexes any chromosome up to 1<<29
(~536 Mbp) and is read transparently by `tabix`, `bcftools`, and htslib.
Round-trip validation (build the index, query it back, confirm the right
records come out) lives in `TestRunCsiReadable` / `TestParity_IndexFiles_Csi`;
when a real `tabix` binary is on `PATH`, `TestRunCsiReadableByRealTabix`
additionally confirms it can read our `.bed.gz` + `.csi`.

## Performance: `--mapq 0` fast path

When no MAPQ filter is in effect (`--mapq 0`, the default), the record
filter binds a fast keep-predicate (`keepRecordNoMapq`) once and omits the
per-read MAPQ comparison from the hot loop, mirroring upstream mosdepth's
fast path. The output is byte-for-byte identical to the general path;
`TestMapqFastPathByteIdentical` proves this across every output file.

## Deviations from upstream

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

**Fragment mode** (`-a/--fragment-mode`) — counts coverage across the
whole template between properly-paired mates rather than the aligned
reads only. Following upstream, only read1 of a proper, non-supplementary
pair contributes; it covers `[min(read,mate) start, +|TLEN|)`. The
per-base output is byte-identical to the upstream v0.3.14 binary
(`TestUpstream_FragmentMode_Parity`). It is mutually exclusive with
`--fast-mode`.

**Quantize** (`-q/--quantize`) — bins per-base depth into the user's
`:`-separated segments and writes `<prefix>.quantized.bed.gz`, mirroring
upstream's segment parsing, `MOSDEPTH_Q*` labelling, and gap behaviour
for depths outside the range. Byte-identical to upstream v0.3.14
(`TestUpstream_Quantize_Parity`).

**Threads** (`-t/--threads`) — BGZF block decompression is spread across
N worker goroutines via the in-tree `pkg/htsgo/bgzf.MultiReader`, which
inflates blocks concurrently and reassembles the decoded byte stream in
order. The decoded bytes — and therefore every output file — are
byte-identical for any thread count; threading only affects decode
throughput. `<2` falls back to the sequential reader. Verified by
`TestThreads_OutputIdentical` and `bgzf.TestMultiReader_MatchesSequential`.

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
- `.csi` indexes round-trip through `CSI.QueryBytes` so the produced
  files are tabix/htslib-readable end-to-end.
- `--mapq 0` fast path is byte-identical to the general path.
- `-d/--d4` D4 output is byte-identical to the upstream `mosdepth_d4`
  binary on the same BAM (`TestD4_UpstreamBinaryParity`), plus a
  writer/reader round-trip on a hand-built track.
- `-a/--fragment-mode`, `-q/--quantize`, and `-t/--threads` are validated
  byte-for-byte against the upstream `mosdepth` v0.3.14 release binary
  (`TestUpstream_FragmentMode_Parity`, `TestUpstream_Quantize_Parity`,
  `TestThreads_OutputIdentical`); set `MOSDEPTH_BIN` to point at a local
  copy. Offline they fall back to internal-consistency checks and log the
  tier rather than skipping silently.

## Status

- v1: per-base, per-region (BED), per-window, thresholds, distribution,
  summary, CSI indexes, `--mapq 0` fast path, byte-identical D4
  per-base output, `--fragment-mode`, `--quantize`, and multi-threaded
  (`-t/--threads`) BGZF decode.
- Roadmap: `-m/--use-median`, default-mode overlap-pair correction, CRAM
  input (depends on the project's CRAM reader landing first).
