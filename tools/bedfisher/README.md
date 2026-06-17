# bedfisher - Fisher's exact test of overlap enrichment

A Go re-implementation of `bedtools fisher`. Computes the 2x2 contingency
table of overlap counts between two BED files over a genome, then runs
Fisher's exact test (two-tailed) on the table.

## Features

- Same heuristic for the `n22` (genome-background) count as upstream
  (`bMean = (1 + queryUnion/queryCount) + (1 + dbUnion/dbCount)`).
- Two-tailed Fisher's exact test ported from htslib's `kfunc.cpp`
  (log-gamma + incremental hypergeometric accumulator).
- Pre-merge mode (`-m`, `--merge`) merges overlapping records in **both**
  inputs (A and B) before the test, matching upstream's `FileRecordMergeMgr`
  (which enables merging for every input file when `-m` is given — not A only).
- Strand filters (`-s`, `-S`) and overlap-fraction filters (`-f`, `-F`,
  `-r`).
- Pure Go, no third-party dependencies.
- Transparent gzip/BGZF input for `-a` and `-b` (the genome file is plain
  text).

## Build

```bash
go build ./tools/bedfisher/cmd/bedfisher
```

## Usage

```bash
bedfisher -a <A.bed> -b <B.bed> -g <genome> [options]
```

### Options

- `-a, --a FILE`  BED file A (queries; required; `-` for stdin; `.gz` ok)
- `-b, --b FILE`  BED file B (database; required)
- `-g, --g FILE`  Chrom-sizes / genome file (required)
- `--output FILE` Output file (default: stdout)
- `-f FRACTION`   Minimum fraction of A overlapped (0..1)
- `-F FRACTION`   Minimum fraction of B overlapped (0..1)
- `-r, --reciprocal`  Apply `-f` to both sides at the same threshold
- `-s, --strand`  Same-strand overlaps only (requires BED6)
- `-S, --opposite-strand`  Opposite-strand overlaps only
- `-m, --merge`   Pre-merge overlapping records in both A and B before the test
- `-h, --help`    Show help
- `-v, --version` Show version

### Output

```text
# Number of query intervals: <a_count>
# Number of db intervals: <b_count>
# Number of overlaps: <overlap_count>
# Number of possible intervals (estimated): <pseudo_count>
# phyper(<n11> - 1, <a_count>, <pseudo - a_count>, <b_count>, lower.tail=F)
# Contingency Table Of Counts
#_________________________________________
#           |  in -b       | not in -b    |
#     in -a | <n11>        | <n12>        |
# not in -a | <n21>        | <n22>        |
#_________________________________________
# p-values for fisher's exact test
left  right  two-tail  ratio
<left>  <right>  <two>  <ratio>
```

`ratio` is printed as `inf` when `n21 == 0`, `-nan` when `n12 == 0` or
`n22 == 0` (matching upstream's behaviour).

## Parity

All five small upstream parity cases (`fisher.t1`..`t4`, `t6`) pass
byte-for-byte. The sixth case (`t5`) only checks that the binary tolerates
a long $TMPDIR file path; it's a CLI / filesystem concern unrelated to the
algorithm and is skipped.

In addition, `live_parity_test.go` runs the real vendored `bedtools fisher`
binary against a heavy-overlap dataset (thousands of self-overlapping A and B
intervals) and asserts byte-for-byte equality, with and without `-m`. This
locks in the fix for the overlap-counting regression: the port previously
binary-searched on `ChromEnd` over a start-sorted B slice — invalid, because
`ChromEnd` is not monotonic there, so a long B that starts before A yet extends
past `A.Start` was skipped. The counter now scans every B with
`ChromStart < A.End` (an exact, monotonic upper bound) and tests the end
coordinate per record, matching upstream's chromsweep `intersects()` predicate.

## Tests

```bash
go test ./tools/bedfisher/...
```
