# bedmakewindows

Pure-Go reimplementation of `bedtools makewindows`. Partitions intervals
(either chromosomes from a genome-sizes file or BED records) into fixed-size
or fixed-count windows.

## Usage

```bash
# Fixed-width windows over a genome
bedmakewindows -g genome.txt -w 1000 -s 500 > windows.bed

# Fixed-count windows over a BED set
bedmakewindows -b regions.bed -n 5 -i srcwinnum > windows.bed
```

## Flags

| Short | Long | Notes |
|---|---|---|
| `-g` | `--genome` | Genome chrom-sizes file (`CHROM\tSIZE`). |
| `-b` | `--bed` | Source BED file (`-` = stdin). |
| `-w` | `--window-size` | Window width in bases. |
| `-s` | `--step-size` | Slide between consecutive windows (default = window-size). |
| `-n` | `--count` | Partition each interval into N equal windows. |
| `-i` | `--id-name` | Naming: `srcwinnum`, `winnum`, `src`, or `none` (default). |
| | `--reverse` | Reverse per-interval window numbering (last window = 1). |
| `-o` | `--output` | Output file (default stdout). |
| `-h` | `--help` | |
| `-v` | `--version` | |

Reads transparent gzip / BGZF via `pkg/htsgo/iohelper`. Stdin / `-`
supported for the BED input.

## Behaviour

Two partition strategies (mutually exclusive):

- **`-w`/`-s`** (fixed-width): each window is up to `-w` bases wide, sliding
  by `-s`. The final window per interval is clipped to the interval end.
- **`-n`** (fixed-count): split each interval into `-n` windows of as-equal
  length as possible. Intervals shorter than `-n` bases are skipped with the
  upstream warning to stderr.

`-i` naming (matching upstream `windowMaker.cpp`):

- default (no `-i`, internally `none`): emit **BED3** — no name column.
- `src`: append the source interval's name. For genome-file input the name
  is the chromosome (upstream builds each genome interval with `name ==
  chrom`).
- `winnum`: append the per-interval 1-based window number (reversed by
  `--reverse`).
- `srcwinnum`: append `<src_name>_<window_number>`.

Upstream accepts only `src`/`winnum`/`srcwinnum` and defaults to no name
column when `-i` is omitted; we additionally accept the explicit spelling
`-i none` as an alias for that default (a one-directional superset — see
[`../../docs/UPSTREAM_BUGS.md`](../../docs/UPSTREAM_BUGS.md#bedtools-makewindows-i-none)).

See [`../PARITY_VALIDATION.md`](../PARITY_VALIDATION.md) for the validated
parity matrix. `pkg/bedmakewindows/live_parity_test.go` asserts byte-for-byte
parity against the upstream `bedtools makewindows` binary for the default and
every `-i` value over both BED and genome input.
