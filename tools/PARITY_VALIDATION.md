# bedtools parity validation

This document tracks the byte-for-byte output parity between the Go ports
under `tools/bed*/` and the upstream `bedtools` C++ test suite (vendored as a
git submodule at `reference_code/bedtools/`).

The methodology is golden-file: each upstream test's input and expected output
were copied into the corresponding tool's `testdata/parity/` directory, and a
`parity_test.go` file invokes the Go port's library function on the same input
and asserts `bytes.Equal` against the upstream expected output. Tests for
options we have not implemented are wrapped in `t.Skip` with a one-line
rationale rather than being deleted, so future work can revive them as the
features land.

## Summary

| Tool          | Tests added | Passed | Skipped | Notes |
| ------------- | -----------:| ------:| -------:| ----- |
| bedclosest    |          14 |      8 |       6 | Skips: `-N`, `-s`/`-S`, `-k`, multi-DB (`-names`/`-filenames`). |
| bedcomplement |          11 |      9 |       2 | Skips: `-L`, the out-of-range warning message form. |
| bedflank      |          11 |     11 |       0 | All upstream cases (`flank.t1..t11`) covered. |
| bedgenomecov  |           9 |      3 |       6 | Skips: BAM/SAM/CRAM input, `-pc`, `-fs`, `-split`. |
| bedintersect  |          13 |      8 |       5 | Skips: `-wo`/`-wao`/`-wa -wb` combined writer. |
| bedjaccard    |          14 |      4 |      10 | Skips: BAM input, `-split`, `-S <strand>`, large fixture t16, and cases that rely on upstream's auto-merging of B (real discrepancy). |
| bedmerge      |          16 |     10 |       6 | Skips: `-delim`, VCF/GFF input, `-S` strand filter, `-s` with mixed `.` strand fan-out (real discrepancy). |
| bedslop       |          15 |     13 |       2 | Skips: float-precision regression tests t13/t14 (require the full `human.hg19.genome` fixture). |
| bedsort       |          11 |      9 |       2 | Skips: `-header` (preserves leading `#` lines); one fixture-layout sanity check. |
| bedsubtract   |          13 |     10 |       3 | Skips: `-N` (union-coverage drop). |
| **TOTAL**     |     **127** | **85** | **42**  | |

(The discrepancy between this table and `go test`'s 87 passed / 42 skipped is
two helper / sanity sub-tests in `bedsort` and `bedintersect` that are not
direct mirrors of an upstream case.)

## What is validated

For each Go port we picked a representative subset (5-15 cases per tool, not
exhaustive) from the upstream `test-<tool>.sh` script, covering:

- the default happy path,
- edge cases (empty input, single record, overlapping inputs),
- strand-aware behaviour where applicable,
- fractional-overlap thresholds where applicable,
- common option combinations.

The inputs vendored under `tools/bed*/testdata/parity/` are byte-for-byte
copies of the upstream fixtures (e.g. `reference_code/bedtools/test/merge/a.bed`).
The `.expected.*` files are literal copies of the inline heredoc expected
strings from the upstream test scripts.

## Discrepancies found (and fixed in this PR)

The validation surfaced a handful of small but real semantic differences from
upstream that we fixed inline rather than masking with `t.Skip`:

- **bedclosest distance is now `(b.start - a.end) + 1`** instead of
  `b.start - a.end`. Upstream's `closest -d` counts touching intervals as 1bp
  apart (not 0); see `src/utils/NewChromsweep/CloseSweep.cpp` for the
  `+ 1` term. Existing bedclosest unit tests were updated to use the new
  convention.
- **bedclosest gained `DistanceAbsolute`**, a `DistanceMode` that emits the
  unsigned distance (matches upstream `-d`; `-D <mode>` is the signed flag).
  The parity tests use this mode.
- **bedsort now breaks size/score ties on the default (chrom asc, start asc,
  end asc) ordering** rather than relying on input order. Without this fix,
  `sort.t03` (`-sizeD`) produced a different ordering of equally-sized
  records than upstream.
- **bedslop now matches upstream's negative-slop semantics**: when slop
  inverts a record (newStart > newEnd) the two coordinates are swapped; when
  slop pushes a record entirely past a chromosome boundary the output is
  pinned to a 1bp slice at that boundary instead of being dropped. This was
  necessary for parity with `slop.t20`, `slop.t21`, and `slop.t22`.
- **bedintersect default mode now preserves A's full columns** (name, score,
  strand, extras) when clipping to the overlap range, instead of emitting
  bare BED3. Same for `-c` (`Count`): the count is appended as a trailing
  column rather than replacing the name. Upstream's `intersect.t01` and
  `intersect.t03` require this column preservation.
- **bedmerge's `ParseColumnOps` now accepts the single-column / many-ops
  form** (e.g. `-c 5 -o count,sum`): each op produces one output column
  computed from the same input column. Upstream `merge.t7` exercises this.
- **bedjaccard and bedgenomecov fraction output now uses %g precision 6**
  (matches C++ ostream's default), not 10. Existing unit tests were updated.

## Discrepancies found (NOT fixed in this PR)

These are real upstream/Go differences that the parity tests document with
`t.Skip("known discrepancy: ...")`. Each one is an open task for a later PR:

- **bedjaccard does not pre-merge B before computing intersection/union**.
  Upstream `bedtools jaccard` first merges overlapping records on both A and
  B (its `n_intersections` is counted against the merged sides). bedjaccard
  counts raw pairs and gets a slightly different intersection sum whenever B
  has overlapping records. Affected: `jaccard.t02/t03/t05/t06/t10/t11`.
- **bedmerge's `-s` does not implement upstream's `.` strand fan-out**.
  Upstream's `merge -s` treats a `.` strand record as belonging to BOTH `+`
  and `-` groups; bedmerge currently propagates the `.` strand through a
  single-pass group. Affected: `merge.t15`.

## What is NOT validated (skipped features)

The following upstream options are intentionally NOT covered by parity tests
because the Go ports do not implement them. Each skipped test names the
missing feature in its `t.Skip` call.

- **Non-BED input formats**: BAM, CRAM, SAM, VCF, GFF anywhere a tool accepts
  `-i`/`-a`/`-b`. The Go ports are BED-only.
- **`-split`** (BED12 block-aware overlap) on `bedintersect`, `bedjaccard`,
  `bedgenomecov`.
- **`-S <strand>`** (single-strand filter): on `bedmerge`, `bedjaccard`.
- **`-N`** (union-coverage drop) on `bedsubtract`. Upstream's `-N` differs
  from our `MinFraction` (which is `-f`).
- **`-N`** (force different names) on `bedclosest`.
- **`-s`/`-S`** strand filters on `bedclosest` (we do not expose these in
  the `Options` struct yet).
- **`-k`** (k-nearest) on `bedclosest`.
- **`-names`/`-filenames`** multi-database labelling on `bedclosest`.
- **`-pc`** (paired-end coverage) and **`-fs`** (fragment size) on
  `bedgenomecov`; both require BAM input.
- **`-wo`/`-wao`/combined `-wa -wb`** writer on `bedintersect`. We support
  `-wa` and `-wb` independently but not the side-by-side composite output.
- **`-delim`** on `bedmerge`. The collapse / distinct ops always join with
  `,`.
- **`-L`** (limit complement output to chromosomes seen in the input file)
  on `bedcomplement`.
- **`-header`** on `bedsort`. Upstream preserves the leading `#` line in the
  output; our reader strips comment lines unconditionally.

## How to run

```bash
# All parity tests, race detector on:
go test -race -run TestParity_ ./tools/bed.../...

# Just one tool:
go test -v -run TestParity_ ./tools/bedmerge/...

# View skipped tests with their rationale:
go test -v -run TestParity_ ./tools/bed.../... 2>&1 | grep SKIP
```

## How to add a new case

1. Find the upstream test under
   `reference_code/bedtools/test/<subcmd>/test-<subcmd>.sh`.
2. Copy the input fixture(s) it uses to
   `tools/bed<subcmd>/testdata/parity/`.
3. Copy the inline heredoc `exp` string to
   `tools/bed<subcmd>/testdata/parity/<case>.expected.*`.
4. Add a `TestParity_<Subcmd>_<TN>_<Name>` to the package's
   `parity_test.go`. The test should call the Go port's library function on
   the input and `bytes.Equal` the output against the expected file.
5. If the case exercises an option the Go port does not implement, replace
   the test body with `t.Skip("unimplemented: <option>")` instead of
   deleting the case.
