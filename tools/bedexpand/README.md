# bedexpand

Pure-Go reimplementation of `bedtools expand`. Replicates lines based on
columns whose values are comma-separated lists, emitting one output row per
list element.

## Usage

```bash
# Expand the 4th column
bedexpand -c 4 -i regions.bed > rows.bed

# Expand two columns in lock-step
bedexpand -c 4,5 -i regions.bed

# Swap two expanded columns (matches upstream expand.t3)
bedexpand -c 5,4 -i regions.bed
```

## Flags

| Short | Long | Notes |
|---|---|---|
| `-c` | `--columns` | 1-based comma-separated column list. Required. |
| `-i` | `--input` | Input file (default stdin; `-` = stdin). Transparent gzip. |
| `-o` | `--output` | Output file (default stdout; `-` = stdout). |
| `-h` | `--help` | |
| `-v` | `--version` | |

## Behaviour

Each input row is tab-split. For every column listed in `-c`, the value is
split on commas. Comma tokenization matches upstream's C++ `getline`
semantics exactly: a single **trailing** comma does not produce a spurious
empty final element (`10,20,30,` -> `10`, `20`, `30`, i.e. three rows, not
four), while leading and interior empty elements are preserved
(`,a,,b` -> `""`, `a`, `""`, `b`). All listed columns within a row must
produce lists of the same length; otherwise an error is raised (matching
upstream).

The output preserves the row's original column layout: when walking the
columns left to right, non-expanded columns are emitted verbatim and each
expanded column is replaced by the *k*-th list element, where *k* is the
position of that column inside `-c` (1-based). So `-c 5,4` substitutes
column-5 elements at position 4 and column-4 elements at position 5, which
mirrors `bedtools expand`'s output for `expand.t3`.

Lines starting with `#`, `track`, or `browser`, and empty lines, are
forwarded verbatim.

See [`../PARITY_VALIDATION.md`](../PARITY_VALIDATION.md) for the validated
parity matrix. `pkg/bedexpand/live_parity_test.go` additionally asserts
byte-for-byte parity against the upstream `bedtools expand` binary on
trailing-comma and leading/interior-empty inputs.
