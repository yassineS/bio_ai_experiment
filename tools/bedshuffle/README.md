# bedshuffle

Pure-Go reimplementation of `bedtools shuffle`. Randomly relocates BED
intervals across a genome, preserving interval length and any extra columns.

## Usage

```bash
bedshuffle -i input.bed -g hg19.genome > shuffled.bed
bedshuffle -i input.bed -g hg19.genome -incl incl.bed -seed 42
bedshuffle -i input.bed -g hg19.genome -excl centromeres.bed
bedshuffle -i input.bed -g hg19.genome -chrom        # keep on original chrom
```

## Flags

| Short | Long | Notes |
|---|---|---|
| `-i` | `--input` | Input BED (required; `-` for stdin). |
| `-g` | `--genome` | Chrom-size table (`chrom<TAB>size`). |
|      | `--output` | Output file (default stdout). |
| `-incl` | `--include` | Include-region BED: every placement must fall inside one of these regions. |
| `-excl` | `--exclude` | Exclude-region BED: no placement may overlap one of these regions. |
| `-chrom` | `--chromOnly` | Keep each interval on its original chromosome. |
| `-chromFirst` | | Pick the destination chromosome uniformly first, then a position within it (default: project a genome-wide position, weighting by chrom size). |
| `-allowBeyondChromEnd` | | Clamp to the chrom end instead of redrawing when an interval would exceed it. |
| `-seed` | `--seed` | Deterministic RNG seed (default 0). |
|         | `--maxTries` | Placement retries per interval (default 1000). |
| `-h` | `--help` | Show help. |
| `-v` | `--version` | Show version. |

Reads transparent gzip / BGZF via `pkg/htsgo/iohelper`. Stdin / `-` is
supported.

## Sampling strategy

Without `-incl`, each interval's new chromosome is drawn with probability
proportional to chromosome length, then the new start is drawn uniformly
inside that chromosome. With `-incl`, chromosomes are weighted by their
total bp listed in the include BED, then a region is picked weighted by its
length, then a uniform offset inside that region.

`-excl` is enforced post-draw: a candidate placement that overlaps any
exclude region is rejected and the draw retried, up to `--maxTries` times.
After that many failures the run aborts with a clear error message that
matches upstream's wording.

`-chrom` keeps the original chromosome and only randomises the start.

## Determinism and upstream parity

The RNG is a pure-Go port of `std::mt19937_64` — the exact 64-bit Mersenne
Twister that upstream bedtools' `Random.cpp` uses (the default, non-`USE_RAND`
build). Combined with the genome-file-order projection and the exact per-mode
draw/retry order from `shuffleBed.cpp`, **`bedshuffle -seed N` reproduces
`bedtools shuffle -seed N` byte-for-byte** for the default, `-incl`, `-excl`,
`-chrom`, `-chromFirst`, and `-allowBeyondChromEnd` modes. This is asserted
directly against the upstream binary in
`pkg/bedshuffle/byte_parity_test.go`. Identical inputs + identical seed always
produce identical output.

## Deviations from upstream

- `-noOverlapping` and `-f` (overlap fraction) upstream flags are not yet
  supported (low-priority; tracked in the roadmap). When a chromosome named in
  the input is absent from the genome under `-chrom`, this port reports the
  interval as unplaceable instead of emitting upstream's garbage (negative)
  coordinates — a fix-on-port.

See [`../PARITY_VALIDATION.md`](../PARITY_VALIDATION.md) for the validated
parity matrix.
