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
| `-n`  | `--number` | Number the blocks (1-based) into the score column. Reverses the numbering on `-` strand records.    |
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
- The output is always tab-separated. Score defaults to `0` (matches upstream)
  unless `-n` is set.
- On `-` strand records `-n` numbers blocks in reverse order so the first
  emitted block carries the highest index (this matches upstream's `t5`).

## Parity

Validated against upstream's `test-bed12tobed6.sh` (cases `t1`-`t5`); see
[`PARITY_VALIDATION.md`](../PARITY_VALIDATION.md#bed12tobed6).
