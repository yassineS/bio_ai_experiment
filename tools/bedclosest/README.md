# bedclosest - Find the closest B interval for each A interval

A Go re-implementation of `bedtools closest`. Output matches the upstream
binary (v2.31.1) byte-for-byte across the `closest` and `kclosest` test suites.

## Description

For each interval in A (sorted), find the closest interval(s) in B (also
sorted) on the same chromosome and report A's columns followed by the chosen
B's columns. A distance column is appended only when requested: `-d` appends
the unsigned distance, `-D <mode>` appends the signed distance. Distance is 0
when A and B overlap; touching intervals report 1 (upstream's `(gap + 1)`
convention). For tied distances the default is to emit one row per tied B
(`-t all`).

Both inputs MUST be sorted on `(chrom, start)`. bedclosest errors out clearly
when they are not.

## Installation

```bash
go build ./tools/bedclosest/cmd/bedclosest
```

## Usage

```text
bedclosest -a <fileA.bed> -b <fileB.bed> [options]
```

## Options

- `-a, --a FILE` - Input BED file A (sorted; use `-` for stdin)
- `-b, --b FILE...` - One or more sorted BED database files (use `-` for stdin).
  With multiple files a database-label column (the 1-based file index, or the
  `-names`/`-filenames` label) is inserted between A's and B's columns.
- `-names NAME...` - Labels for the `-b` databases (one per file, in order);
  replaces the numeric file-index column. Mutually exclusive with `-filenames`.
- `-filenames` - Use each `-b` file's name as its database-label column.
- `-mdb each|all` - Multi-database mode: `each` (default) reports the closest
  feature from every database on its own row; `all` reports the single overall
  closest across all databases.
- `-o, --output FILE` - Output BED file (`-` for stdout, default: stdout)
- `-d` - Append the unsigned distance to the closest B. Mutually exclusive with
  `-D`.
- `-D MODE` - Append the signed distance with the chosen sign convention:
  - `ref`: downstream of A (to the right) is positive, upstream negative
  - `a`: relative to A's strand (BED6 col 6); flips on a `-`-strand A
  - `b`: relative to B's strand; flips on a `+`-strand B
- `-k N` - Report the `N` closest hits per A interval (default 1).
- `-N` - Require the closest B to have a different name (BED column 4) than A;
  a B sharing A's name is skipped from candidate consideration.
- `-io` - Ignore B features that overlap A.
- `-iu` - Ignore B features upstream of A (requires `-D`).
- `-id` - Ignore B features downstream of A (requires `-D`).
- `-fu` - Report the closest upstream B in preference to overlaps/downstream
  (requires `-D`).
- `-fd` - The downstream-forcing counterpart of `-fu` (requires `-D`).
- `-s` - Require the closest B to be on the SAME strand as A (BED6 col 6).
  Non-matching B intervals are skipped from candidate consideration.
- `-S` - Require the closest B to be on the OPPOSITE strand to A. Mutually
  exclusive with `-s`.
- `-t MODE` - Tie-break among equally-close B's:
  - `all` (default) - emit one row per tied B in B's input order
  - `first` - emit only the first tied B
  - `last` - emit only the last tied B
- `-header` - Print A's header (comment/`track`/`browser`) lines before results.
- `-h, --help` - Show help message
- `-v, --version` - Show version (`1.0.0`)

## Examples

```bash
# Closest peak for each gene
bedclosest -a genes.sorted.bed -b peaks.sorted.bed > out.bed

# Append the unsigned distance column
bedclosest -a a.bed -b b.bed -d > out.bed

# Append the signed distance (upstream of A is negative)
bedclosest -a a.bed -b b.bed -D ref > out.bed

# Report the 5 closest features per A
bedclosest -a a.bed -b b.bed -k 5 > out.bed

# Only report a B with a different name than A
bedclosest -a a.bed -b b.bed -N > out.bed

# Ignore B's that overlap A
bedclosest -a a.bed -b b.bed -io > out.bed

# Closest B on the same strand as A (skips opposite-strand B's)
bedclosest -a a.bed -b b.bed -s > out.bed

# Closest feature from each of several databases (one row per database;
# the inserted column is the 1-based database index)
bedclosest -a a.bed -b db1.bed db2.bed db3.bed > out.bed

# Label the database column with names instead of indices
bedclosest -a a.bed -b db1.bed db2.bed db3.bed -names a b c > out.bed

# Single overall closest across all databases (still labelled by source DB)
bedclosest -a a.bed -b db1.bed db2.bed db3.bed -mdb all > out.bed
```

## Format

- Input: BED (tab-delimited, minimum 3 columns), sorted on `(chrom, start)`.
  `.gz` is supported.
- Output: A's columns, then B's columns, then the distance (only when `-d` or
  `-D` is given). With multiple `-b` databases, a database-label column (1-based
  index, or the `-names`/`-filenames` label) is inserted between A's and B's
  columns.
- When A's chromosome has no B records, a "null" placeholder B is emitted whose
  shape matches the B file's record type, exactly as upstream's
  `RecordOutputMgr::printNull`:
  - BED3: `.\t-1\t-1`
  - BED4: `.\t-1\t-1\t.`
  - BED5: `.\t-1\t-1\t.\t-1`
  - BED6: `.\t-1\t-1\t.\t-1\t.`
  - BED12: `.\t-1\t-1\t.\t-1\t.` followed by six more `.`
  - BedGraph: `.\t-1\t-1\t.`
  - BED4+/BED6+ (extra columns): `.\t-1\t-1` followed by a `.` per extra column

  With `-d`/`-D` a trailing `-1` distance is appended. With multiple databases
  the database column for the null row is a literal `.`, and exactly one null
  row is emitted only when *no* database yields any hit.

## Known parity gaps

The cross-file chromosome sort-order / naming-convention validation engine
(upstream's `testChromOrder` / `testNameConventions`, exercised by the
`sortAndNaming` sub-suite) is not implemented. bedclosest enforces only the
basic per-file `(chrom, start)` sort check. The exact lexicographic-vs-numeric
ordering detection, `chr`-prefix and leading-zero naming WARNINGs, and the
"contains chromosome X but query does not" ERROR messages are a distinct
subsystem tied to upstream's interleaved streaming merge.

## Algorithm

Records are read into memory and bucketed by chromosome. Per A interval, every
B on A's chromosome is classified as overlapping, upstream, or downstream and
collected into per-stream distance-bucketed lists capped at the `k` closest
distinct distances (mirroring upstream's `RecDistList`). The final hit list is
then assembled the way upstream's `CloseSweep::finalizeSelections` does: forced
streams first (`-fu`/`-fd`), then overlaps, then the closest of the remaining
upstream/downstream features merged outward by distance, honouring the tie mode
until `k` hits are produced.
