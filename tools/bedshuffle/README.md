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
| `-seed` | `--seed` | Deterministic RNG seed (default 0). |
|         | `--maxTries` | Placement retries per interval (default 1000). |
| `-h` | `--help` | Show help. |
| `-v` | `--version` | Show version. |

Reads transparent gzip / BGZF via `pkg/bioformats/iohelper`. Stdin / `-` is
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

## Determinism

The RNG is `math/rand` seeded by `-seed`. Identical inputs + identical seed
always produce identical output. Different seeds produce different output.

## Deviations from upstream

- **Byte-for-byte output is not parity** with upstream's bedtools shuffle
  for any seed: upstream uses its own Mersenne Twister implementation with
  a specific seeding regime that we do not reproduce. We instead validate
  the *structural* invariants the upstream tests were designed to check
  (length preserved, include / exclude / chrom honoured, error on
  unplaceable intervals). Tracked in `docs/PARITY_ROADMAP.md#bedtools`.
- `-chromFirst` (the alternative sampling order) is treated as the default
  in our port — the two are equivalent when `-incl` is present and the
  include list covers the same chrom subset.
- `-noOverlapping`, `-allowBeyondChromEnd`, `-f` upstream flags not yet
  supported (low-priority; tracked in the roadmap).

See [`../PARITY_VALIDATION.md`](../PARITY_VALIDATION.md) for the validated
parity matrix.
