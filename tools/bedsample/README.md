# bedsample

Pure-Go reimplementation of `bedtools sample`. Draws N random records from a
BED-like input without replacement, via single-pass reservoir sampling.

## Usage

```bash
bedsample -i regions.bed -n 100 -seed 42 > sampled.bed
sort -k1,1 -k2,2n regions.bed | bedsample -n 10
```

## Flags

| Short | Long | Notes |
|---|---|---|
| `-n` | `--number` | Number of records to draw. Required, must be > 0. |
| `-seed` | `--seed` | PRNG seed for deterministic output. `0` (default) = time-based. |
| `-header` | `--header` | Forward `#`, `track`, `browser` lines verbatim before the records. |
| `-i` | `--input` | Input BED (default stdin; `-` = stdin). Transparent gzip. |
| `-o` | `--output` | Output file (default stdout; `-` = stdout). |
| `-h` | `--help` | |
| `-v` | `--version` | |

## Behaviour

Single linear pass with reservoir sampling: O(N) extra memory, O(total)
time. Output preserves the input file's relative order of the sampled
records.

The PRNG is Go's `math/rand` seeded by `-seed`; this is reproducible across
runs of the same binary but is **not byte-identical** to upstream's C++
implementation, which uses a different sampler. With the same seed two
`bedsample` runs produce identical output (mirrors upstream's `-seed`
contract); cross-tool reproduction is not guaranteed.

Requesting more records than the file contains is an error and produces the
upstream-style message `Input file has fewer records than the requested
number of output records`.

See [`../PARITY_VALIDATION.md`](../PARITY_VALIDATION.md) for the validated
parity matrix.
