# bed12tobed6

Go re-implementation of `bedtools bed12tobed6`.

Splits each BED12 record into one BED6 record per block (`chrom`, `start`,
`end`, `name`, `score`, `strand`).

## Usage

```bash
bed12tobed6 [options] [< input.bed]
```

### Options

| Short | Long       | Description                                                                                          |
| ----- | ---------- | ---------------------------------------------------------------------------------------------------- |
| `-i`  | `--input`  | Input BED12 file (`-` for stdin, default: stdin).                                                    |
| `-o`  | `--output` | Output BED6 file (`-` for stdout, default: stdout).                                                  |
| `-n`  | `--number` | Number the blocks (1-based) into the score column. Reverses the numbering for any non-`+` strand.   |
| `-h`  | `--help`   | Show help.                                                                                           |
| `-v`  | `--version`| Show version.                                                                                        |

`-i -` and `-o -` use stdin / stdout. Gzipped input is auto-detected via
[`pkg/htsgo/iohelper`](../../pkg/htsgo/iohelper).

## Examples

```bash
# Split a BED12 file into BED6 blocks.
bed12tobed6 -i blocks.bed > out.bed

# Number the blocks per record.
bed12tobed6 -i blocks.bed -n > out.bed
```

## Behaviour notes

- BED12 input must have at least 12 columns. Records with fewer columns are
  passed through unchanged (matches upstream behaviour when given BED6/BED4).
- The output is always tab-separated. Each emitted BED6 block carries the
  parent record's score (column 5) unchanged, exactly like upstream
  (`GetBedBlocks` copies `bed.score` onto every block). When `-n` is set the
  score column is replaced by the 1-based block number instead.
- Under `-n`, blocks are numbered in reverse order for **any** strand that is
  not exactly `+` (i.e. `-`, `.`, or empty), so the first emitted block carries
  the highest index. This matches upstream's `strand == "+"` check in
  `bed12ToBed6.cpp` (covered by `t5` for the `-` strand).

## Parity

Validated against upstream's `test-bed12tobed6.sh` (cases `t1`-`t5`); see
[`PARITY_VALIDATION.md`](../PARITY_VALIDATION.md#bed12tobed6). A live-binary
parity test (`pkg/bed12tobed6/live_parity_test.go`) additionally proves
score propagation and `-n` numbering — including a `.`-strand record —
byte-for-byte against the upstream `bedtools bed12tobed6` binary.
