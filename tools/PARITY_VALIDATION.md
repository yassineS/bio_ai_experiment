# Parity validation

This document tracks the byte-for-byte output parity between the Go ports
in this repository and the upstream tools' regression test corpora
(vendored under `reference_code/`).

The initial corpus is bedtools (PR #55, 127 cases). Subsequent sections
extend the same methodology to the remaining ports.

## bedtools

Bytes-for-byte parity against the upstream `bedtools` C++ test suite
(vendored as a git submodule at `reference_code/bedtools/`).

The methodology is golden-file: each upstream test's input and expected output
were copied into the corresponding tool's `testdata/parity/` directory, and a
`parity_test.go` file invokes the Go port's library function on the same input
and asserts `bytes.Equal` against the upstream expected output. Tests for
options we have not implemented are wrapped in `t.Skip` with a one-line
rationale rather than being deleted, so future work can revive them as the
features land.

### Summary

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

The sickle and skewer ports each have a per-tool table in their respective
section below; the project-wide running total is 154 tests added, 107
passed, 44 skipped (47 if you count the sickle FixturesPresent + skewer
FixturesPresent + skewer PEHelperSmoke helper tests).

### What is validated

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

### Discrepancies found (and fixed in PR #55)

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

### Discrepancies found (NOT fixed)

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

### What is NOT validated (skipped features)

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

### How to run

```bash
# All parity tests, race detector on:
go test -race -run TestParity_ ./tools/bed.../...

# Just one tool:
go test -v -run TestParity_ ./tools/bedmerge/...

# View skipped tests with their rationale:
go test -v -run TestParity_ ./tools/bed.../... 2>&1 | grep SKIP
```

### How to add a new case

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

## sickle

Byte-for-byte parity against the upstream C `sickle` v1.33 binary built
from `reference_code/sickle`. Upstream compiles cleanly with `make`
against gcc 13 / glibc 2.38 on this environment (two
`-Wswitch-unreachable` warnings, no errors); the resulting binary is
~52 KB and self-contained. No extra setup is needed beyond
`git submodule update --init reference_code/sickle`.

Fixtures under `tools/sickle/testdata/parity/` were generated by running
the upstream binary with the exact flags noted in each test's docstring
and capturing stdout into `case*.expected.fq`. Go tests in
`tools/sickle/pkg/sickle/parity_test.go` drive the in-process library API
on the same inputs and assert `bytes.Equal` against the expected files.

### sickle summary

| Subcommand | Tests added | Passed | Skipped | Notes |
| ---------- | -----------:| ------:| -------:| ----- |
| se         |          10 |     10 |       0 | Basic / `-n` / `-x` / illumina (Phred+64) / empty / all-low / Q-threshold boundary / gzip / short-read filter. |
| pe         |           3 |      3 |       0 | Basic, `-s` singletons, synced pass/fail. |
| **TOTAL**  |       **13**|  **13**|    **0**| Includes one fixture-present smoke test. |

### sickle: discrepancies found in our port (fixed in this PR)

The validation surfaced one substantial divergence from upstream that we
fixed inline rather than masking with `t.Skip`:

- **Sliding-window trim rewritten to match upstream byte-for-byte.** The
  pre-existing `tools/sickle/pkg/sickle/sickle.go:trimRecord` implemented a
  separate-pass 5'/3' search with its own per-direction sliding window.
  That algorithm produced different cut sites from upstream's single-pass
  left-to-right scan whenever:
  - the window size differed (upstream uses `int(0.1 * len(seq))`; our
    port used a fixed configurable default of 10),
  - the integer-vs-float average semantics diverged (upstream compares a
    `double` window average to an `int` threshold; our port used integer
    division and dropped fractional information),
  - the read had **no** window whose average reached the threshold —
    upstream discards the whole read as "no 5' cut found", our port
    silently kept the untrimmed read.
  The new `trimRecord` follows `reference_code/sickle/src/sliding.c` step
  for step (see the comment block in the function). `opts.WindowSize` is
  now treated as a Go-port extension: `0` means "auto = int(0.1·L)"
  (upstream semantics, the default through the CLI), a positive value
  overrides. `opts.LengthThreshold` is also enforced as an up-front
  filter (upstream's `if (fqrec->seq.l < length_threshold)`).

- **`trimRecord` now takes a `fastq.QualityEncoding` parameter** so the
  Phred offset (33 vs 64) used to decode quality bytes matches upstream's
  per-encoding `quality_constants` table. Phred+64 reads were previously
  decoded with offset 33, producing nonsense decoded scores for illumina
  / solexa input.

- **The old standalone `trim5Prime` / `trim3Prime` helpers were removed**
  because the new single-pass algorithm has no equivalent split. The two
  unit tests that called them directly (`TestTrim5Prime`,
  `TestTrim3Prime`) were removed; equivalent coverage now comes from the
  parity tests and the `TestTrimRecord*` cases that exercise the public
  trim API.

- **`TestPairedAutoDetectFromR1` updated.** The pre-existing test
  asserted that a `det/1` Q0-quality "detection-hint" record would
  round-trip through the trimmer unchanged. With the upstream-faithful
  algorithm that record is correctly discarded (its window average never
  reaches the threshold), so the assertion was inverted to match —
  `det/1`/`det/2` must NOT appear in the paired output.

### sickle: discrepancies found (NOT fixed)

None. Every parity case passes byte-for-byte.

### sickle: what is NOT validated (skipped features)

The Go port exposes several conveniences upstream sickle does not
(`--json`, `--html`, `--progress`, `--auto-detect`, `--recalibrate`, the
`batch` subcommand). These are intentionally NOT covered by parity tests
— there is no upstream behaviour to validate against.

### sickle: bugs found in upstream

None. The behaviours that initially looked suspicious — discarding an
entire read when no window ever reaches the threshold; using a dynamic
`int(0.1·L)` window instead of a fixed default — are documented in
upstream's `src/sliding.c` and consistent with sickle's README. They are
features, not bugs.

## skewer

Byte-for-byte parity against the upstream C++ `skewer` 0.2.2 binary
built from `reference_code/skewer`.

**Build note.** Upstream's `src/matrix.h` declares an `ElementComparator`
whose `operator()` is not `const`-qualified. Modern libstdc++ (gcc 13+)
requires `set<>` comparators to be invocable as `const`, so the upstream
source fails to compile out of the box on this environment. We apply a
minimal patch (add `const` to `ElementComparator::operator()`) to the
submodule working tree before running `make`; the patch is not committed
back to the submodule and is documented in
[docs/UPSTREAM_BUGS.md](../docs/UPSTREAM_BUGS.md#skewer).

### skewer summary

| Subcommand | Tests added | Passed | Skipped | Notes |
| ---------- | -----------:| ------:| -------:| ----- |
| se         |          11 |      9 |       1 | Skip: case05 error-tolerance (algorithm difference, see below). |
| pe         |           1 |      0 |       1 | Skip: case04 PE matrix mode not implemented in Go port. |
| **TOTAL**  |       **14**|   **9**|    **2**| Plus 2 smoke tests (fixture presence + PE helper). |

### skewer: discrepancies found in our port (fixed in this PR)

None — the Go port already passed 9/12 cases byte-for-byte on first run
against upstream. The remaining 3 are documented divergences (one PE
algorithm gap, one matcher difference, one upstream non-issue).

### skewer: discrepancies found (NOT fixed)

These are documented divergences the parity tests record with `t.Skip`:

- **case04 — PE matrix mode (`-m pe`)**. With the default paired-end
  mode upstream's `Matrix::Detect` runs a paired-end overlap check
  between R1 and R2 and refuses to trim when the mates disagree on the
  insert size. Our Go port has no equivalent matrix logic —
  `TrimPairedEnd` just runs the per-read 3' adapter trimmer on each
  mate independently. For the test reads in this corpus the upstream
  behaviour is "don't trim" while ours is "trim each mate".
  Implementing the matrix path is tracked in
  [PARITY_ROADMAP.md](../docs/PARITY_ROADMAP.md#skewer).

- **case05 — error-tolerant matcher**. With `-r 0.1` over a 13 bp adapter
  the upstream Smith-Waterman-like matcher rejects a 1-mismatch match
  whose mismatch is in the last 4 bases of the adapter (asymmetric tail
  penalty). Our Go port uses a simpler Hamming-distance matcher that
  accepts the 1-mismatch alignment and over-trims one base. Bringing the
  matcher into byte parity needs a small Smith-Waterman implementation
  with the same tail-penalty curve — also tracked in
  [PARITY_ROADMAP.md](../docs/PARITY_ROADMAP.md#skewer).

### skewer: what is NOT validated (skipped features)

The Go port exposes adapter auto-detection, UMI extraction, the `batch`
subcommand, and `--json`/`--html-report` outputs. None of these exist
upstream — they are NOT covered by parity tests by design.

The following upstream knobs are not exercised because the Go port does
not implement them:

- **`-M, --matrix <file>`** valid-pairing matrix (PE).
- **`-j <str>`** junction adapter for mate-pair (Nextera Mate Pair).
- **`-b, --barcode`** demultiplexing.
- **`-c, --cut <int>,<int>`** hard-clip barcode length.
- **`-e, --cut3`** hard-clip 3' overhang past `-L`.
- **`-Q, --mean-quality`** drop reads below a mean quality.
- **`-L, --max`** maximum read length.
- **`-N, --fillNs`** replace trimmed bases with Ns.
- **`-A, --masked-output`** lowercase trimmed bases instead of cutting.
- **`-X, --excluded-output`** save discarded reads.
- **`--qiime`** QIIME barcode/mapping output.

### skewer: bugs found in upstream

One; recorded in
[docs/UPSTREAM_BUGS.md](../docs/UPSTREAM_BUGS.md#skewer).

- **Modern-libstdc++ compile failure**. `ElementComparator::operator()`
  is not `const`-qualified. **Fix-on-port** (build-side only): we patch
  the submodule working tree before building the parity binary. The Go
  port doesn't carry the underlying bug because the matrix-mode code
  path is not yet implemented in Go.
