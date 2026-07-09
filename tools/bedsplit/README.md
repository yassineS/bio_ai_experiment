# bedsplit

Go re-implementation of `bedtools split`: split a single BED file into `N`
output shards, writing each shard to `<prefix>.NNNNN.bed` (1-based,
5-digit zero-padded) and emitting a manifest TSV
(`filename<TAB>total_bp<TAB>num_records`) to stdout.

## Usage

```bash
bedsplit -i input.bed -n 10 -p out -a size
```

## Flags

| Short | Long          | Notes                                                       |
| ----- | ------------- | ----------------------------------------------------------- |
| `-i`  | `--input`     | Input BED file (`-` / omitted = stdin).                     |
| `-p`  | `--prefix`    | Output filename prefix (required).                          |
| `-n`  | `--number`    | Number of output files (required, `>= 1`).                  |
| `-a`  | `--algorithm` | `size` (default) or `simple`.                               |
| `-h`  | `--help`      |                                                             |
| `-v`  | `--version`   |                                                             |

## Behaviour

Two partitioning algorithms (matching upstream `splitBed.cpp`):

- **`simple`**: route records round-robin so each file gets approximately
  equal record counts; record `i` lands in file `i % N` (like Unix
  `split`).
- **`size`** (default): a bin-packing heuristic that balances total bases
  across files. Records are loaded into memory and sorted by length
  **descending**; the first `N` records seed `N` bins, then each remaining
  record is placed into the bin that minimises the sum of absolute
  deviations of bin sizes from the expected mean (the first bin achieving
  the minimum wins). Within a bin, records are written in the order they
  were added — i.e. size-descending — exactly as upstream stores and
  writes its `items` vector (do **not** re-sort the bin back into input
  order, which was a prior parity bug).

When `N` exceeds the record count, only `min(N, records)` non-empty files
are produced (matching upstream).

## Parity

`pkg/bedsplit/parity_test.go` checks the manifest against the upstream
`test-split.sh` goldens (`-a simple` and `-a size`).
`pkg/bedsplit/live_parity_test.go` additionally runs the upstream
`bedtools split -a size` binary and asserts byte-for-byte parity of both
the manifest **and** every emitted shard file across several `-n` values
(including `-n` greater than the record count), on a fixture of
distinct-length records so the size-descending sort order is unambiguous.

### Equal-length tie order

When several records share the same length, which shard each lands in is
decided by the C++ standard library's `std::sort` **tie order for
equal-key elements** — a *stdlib-defined* detail, not a bedtools one. Our
`pkg/cppsort` is a libstdc++ introsort port, so it reproduces the
**libstdc++** upstream (the CI/container oracle) byte-for-byte; a
**libc++** oracle (e.g. a local arm64-macOS `bedtools`) may place a tied
record in a different shard. Only the per-shard membership of tied records
can differ — never the per-file bp totals, record counts, or the overall
partition. `live_parity_test.go`'s `TestLiveParity_SizeTie_*` cases assert
those order-**independent** invariants so they hold on either oracle.
Switching `pkg/cppsort` to a libc++ variant is a deferred, owner-gated
decision (it would flip parity from the CI oracle to a local one).
