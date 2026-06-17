# bedannotate - Annotate A with overlap stats from N BED files

A Go re-implementation of `bedtools annotate`. For each interval in a
primary BED file (`-i`), it annotates with overlap statistics drawn from
one or more secondary BED files (`--files`).

## Features

- Default mode emits the fraction of A covered by each B (per-B float
  column, `%f`).
- `--counts` emits the count of overlapping records per B.
- `--both` interleaves count + fraction per B (2N columns total).
- Strand filters: `-s` (same-strand only), `-S` (opposite-strand only).
  Strand comparison matches upstream's raw-string semantics: two records with
  no strand (BED3, defaulting to `.`) count as same-strand under `-s`.
- A `#` header line is emitted **only** when `-names` is given (matching
  upstream — file basenames do NOT trigger a header).
- Records are reported grouped by chromosome (lexicographic), then by UCSC
  bin (matching upstream's `map<chrom, map<bin, ...>>` iteration order).
- Interval-tree based overlap (`pkg/htsgo/bed.IntervalTree`) — one
  tree per B file built once, then A is buffered, bin-sorted, and emitted.
- Pure Go, no third-party dependencies.
- Transparent gzip/BGZF input on every input and `-` for stdin on `-i`.

## Build

```bash
go build ./tools/bedannotate/cmd/bedannotate
```

## Usage

```bash
bedannotate -i <A.bed> --files <B1.bed> [<B2.bed> ...] [options]
```

### Options

- `-i, --input FILE`     A intervals (required, `-` for stdin)
- `--files FILE..`       One or more B files (variadic; stops at next flag)
- `--names N1 N2 ..`     Header labels (variadic, like upstream; a single
                         comma-separated token `b1,b2` is also accepted). A
                         header line is emitted only when this flag is given.
- `--counts`             Emit per-B count of overlapping records
- `--both`               Emit count and coverage fraction per B (interleaved)
- `-s, --strand`         Same-strand overlaps only
- `-S, --opposite`       Opposite-strand overlaps only
- `-o, --output FILE`    Output file (default: stdout)
- `-h, --help`           Show help
- `-v, --version`        Show version

### Output

```text
chrom <TAB> start <TAB> end <TAB> ...A columns... <TAB> <col_per_B>...
```

With `--names`, a leading `#` header line precedes the data. The `#` is padded
with `bedType-1` tabs (so the first label aligns under the first appended
column), exactly as upstream emits it — for a BED3 main file:

```text
#<TAB><TAB><TAB>exons<TAB>introns
```

With `--both`, the header splits each name into `<name>_cnt` and
`<name>_pct`.

## Parity

The fixtures under `testdata/parity/` have expected outputs generated directly
from the upstream `bedtools annotate` binary, and `live_parity_test.go` runs the
real vendored binary and asserts byte-for-byte equality over multi-file `-files`
inputs (default / `-counts` / `-both -names`), a 500-interval/4-chromosome
ordering stress, and the `-s` / `-S` strand cases. This locks in the two fixes
the parity pipeline found: the port no longer emits a spurious header without
`-names`, and it now reproduces upstream's per-chromosome/per-UCSC-bin record
ordering (previously it emitted records in input order).

## Tests

```bash
go test ./tools/bedannotate/...
```

## Performance

`O(sum(|B_i|))` working set (per-B interval tree) and `O(|A| * log|B_i| + k)`
query time, where `k` is overlap count. Suitable for thousands-to-millions
of A intervals against megabyte-scale B files.
