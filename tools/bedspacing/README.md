# bedspacing

Pure-Go reimplementation of `bedtools spacing`. Reports the distance to the
previous interval on the same chromosome as a new column.

## Usage

```bash
bedspacing -i regions.bed > with-spacing.bed
sort -k1,1 -k2,2n regions.bed | bedspacing
bedspacing -i aln.bam > with-spacing.bed   # BAM/SAM auto-detected
```

## Flags

| Short | Long | Notes |
|---|---|---|
| `-i` | `--input` | Input BED (default stdin; `-` = stdin). Transparent gzip. |
| `-o` | `--output` | Output file (default stdout; `-` = stdout). |
| `-h` | `--help` | |
| `-v` | `--version` | |

## Output

For each input record, the original row is emitted verbatim plus one
appended tab-separated column:

| Value | Meaning |
|---|---|
| `.` | first interval on this chromosome |
| `-1` | overlaps the previous interval on this chromosome |
| `0` | exactly abuts (`this.start == prev.end`) |
| `N` (positive) | `this.start - prev.end`, the gap in bases |

## Behaviour notes

- bedspacing does **not** sort. Feed it a sorted BED (e.g. via `bedsort`).
- The "previous interval" is the most recent record seen on each
  chromosome. A contained record (e.g. `60-80` then `65-70`) advances the
  pointer; the following record's gap is measured from `70`, not `80`. This
  matches upstream `src/spacingFile/spacingFile.cpp`.
- Header lines (`#`, `track`, `browser`) and blank lines are passed
  through unchanged.
- SAM/BAM input is auto-detected from the leading bytes. Each mapped
  alignment is converted to its BED12 representation (the same conversion
  upstream applies under `bedtools spacing -i in.bam -bed`) before the
  spacing column is appended; spacing is measured on the whole reference
  span of the alignment.

See [`../PARITY_VALIDATION.md`](../PARITY_VALIDATION.md) for the validated
parity matrix.
