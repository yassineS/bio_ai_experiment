# bedmap

Pure-Go reimplementation of `bedtools map`. For each interval in A, find
overlapping intervals in B, take values from one or more B columns, apply an
aggregation op (sum, mean, collapse, ...), and append the result to A.

## Usage

```bash
bedmap -a A.bed -b B.bed                # default: sum of B's column 5
bedmap -a A.bed -b B.bed -c 4 -o collapse
bedmap -a A.bed -b B.bed -c 5 -o min,max
```

## Flags

| Short | Long | Notes |
|---|---|---|
| `-a` | `--input-a` | A intervals (required; `-` for stdin). |
| `-b` | `--input-b` | B intervals (required). |
|      | `--output` | Output file (default stdout). |
| `-c` | `--columns` | Comma-separated 1-based B column indices to aggregate (default `5`). |
| `-o` | `--operations` | Comma list of aggregation ops; one per column or one shared (default `sum`). Supported: `sum`, `min`, `max`, `mean`, `median`, `count`, `count_distinct`, `distinct`, `collapse`, `first`, `last`, `mode`, `antimode`. |
|      | `--null` | Placeholder for "no overlap" (default `.`). |
|      | `--delim` | Separator for `collapse` / `distinct` (default `,`). |
| `-s` | `--strand` | Same-strand only. |
| `-S` | `--opposite` | Opposite-strand only. |
| `-f` | `--fraction-a` | Min fraction of A that must overlap a B record. |
| `-F` | `--fraction-b` | Min fraction of B that must overlap A. |
| `-r` | `--reciprocal` | Require both `-f` and `-F`. |
| `-h` | `--help` | Show help. |
| `-v` | `--version` | Show version. |

Reads transparent gzip / BGZF via `pkg/htsgo/iohelper`. Stdin / `-` is
supported.

## Op semantics

`bedmap` reuses `bedmerge.ApplyOp` so the supported op vocabulary is exactly
the same as bedmerge / bedgroupby. See `tools/bedmerge/README.md` for the per-op
semantics.

When there are no overlapping B records:

- `count` and `count_distinct` always emit `0`.
- Every other op emits the `--null` placeholder (default `.`).

## Database (`-b`) input formats

The `-b` database format is auto-detected exactly as upstream
`BedFile::parseLine` (precedence: BED, then VCF, then GFF) plus a BAM magic
sniff:

- **BED / BED+** — columns 2,3 are the 0-based half-open span.
- **VCF** — POS-1 .. POS-1+len(REF); `-c` extracts the literal VCF column.
- **GFF** — 1-based start/end in columns 4/5; `-c` extracts the literal GFF column.
- **BAM** — mapped alignments are projected to a BED12 record (chrom, start,
  end, name, MAPQ, strand, thick, RGB, blockCount, sizes, starts), so e.g.
  `-c 5 -o mean` averages MAPQ.

## Additional flags

- `-header` echoes A's leading comment/`track`/`browser` header lines verbatim.
- `-g/--genome FILE` is accepted as a chromosome-order hint (the port does not
  enforce sort order, so it is effectively a no-op).
- `-split` treats BED12/BAM A and B records as their constituent blocks for
  overlap detection.
- `-prec INT` sets the significant digits used to format numeric op results
  (default 10, matching upstream).

## Error / warning parity

- Requesting a `-c` column outside the database's field range (including `0`
  and negatives) reproduces upstream's exact stderr block
  (`***** ERROR: Requested column N, but database file <name> only has fields 1 - M.`).
- A numeric op (`sum`, `mean`, `min`, `max`, `absmin`, `absmax`, `median`,
  `stdev`, `sstdev`) on a non-numeric value warns once per output row
  (`***** WARNING: Non numeric value <v> in <col>.`) and emits the null value
  when the result is NaN.

## Deviations from upstream

- `-sorted` (sorted-stream fast path) is accepted as a no-op.

See [`../PARITY_VALIDATION.md`](../PARITY_VALIDATION.md) for the validated
parity matrix.
