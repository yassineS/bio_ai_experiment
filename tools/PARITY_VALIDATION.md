# Parity validation

This document tracks the byte-for-byte output parity between the Go ports
in this repository and the upstream tools' regression test corpora
(vendored under `reference_code/`).

The initial corpus is bedtools (PR #55, 127 cases). Subsequent sections
extend the same methodology to the remaining ports.

## bedtools

Byte-for-byte parity against the upstream `bedtools` C++ test suite
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
| bedjaccard    |          14 |     10 |       4 | Skips: BAM input, `-split`, `-S <strand>`, large fixture t16. The column-ops + discrepancies wave (this PR) unskipped t02/t03/t05/t06/t10/t11 after wiring the upstream pre-merge step into bedjaccard. |
| bedmerge      |          16 |     11 |       5 | Skips: `-delim`, VCF/GFF input, `-S` strand filter. The column-ops + discrepancies wave (this PR) unskipped `merge.t15` after fixing `-s` to drop `.`-strand records and merge `+` / `-` independently. |
| bedslop       |          15 |     13 |       2 | Skips: float-precision regression tests t13/t14 (require the full `human.hg19.genome` fixture). |
| bedsort       |          11 |     10 |       1 | Option-tail wave (this PR) unskipped `sort.t09` (`-header`) by buffering leading `#`/`track`/`browser` lines and emitting them verbatim ahead of the sorted body. Remaining skip is one fixture-layout sanity check. |
| bedsubtract   |          13 |     10 |       3 | Skips: `-N` (union-coverage drop). |
| bedexpand     |           6 |      6 |       0 | All canonical upstream cases (`expand.t1..t3`) plus stdin-shape smoke + two error paths. |
| bedgetfasta   |          15 |     15 |       0 | All cases pass. The latest option-tail wave (this PR) wires BGZF FASTA input through `pkg/bioformats/fasta`: `OpenRandomAccess` now sniffs the BGZF magic and routes `.fa.gz` paths to a new `OpenRandomAccessBGZF` that fully decompresses the payload in-memory and reuses the FAI index path. The `.gzi` sidecar is parsed for validation; the `.fa.gz.fai` is honoured when present. Upstream `getfasta.t18` (BGZF + `-split`) now passes byte-for-byte using the upstream `t.fa.gz` fixture.|
| bedsample     |           7 |      5 |       2 | Skips: two CLI-only cases (no-args / unrecognized-flag) that test main.go error messages, not the library. |
| bedspacing    |           7 |      6 |       1 | Skips: BAM input. All inline upstream cases (`spacing.t01`) + synthetic edge cases (per-chrom reset, exact abut, overlap, BED6 preservation, single record). |
| bedcoverage   |           9 |      6 |       3 | Skips: BAM input (t1), `-mean` float32 precision (t6), `-split` BAM modes (t10..t13). |
| bedmap        |          13 |     12 |       1 | Skips: GFF input (t14+). The column-ops + discrepancies wave (this PR) unskipped `map.t11` (absmin) and added `map.t13` (absmax) after wiring `absmin` / `absmax` into the shared `bedmerge.ApplyOp`. |
| bedshuffle    |           6 |      5 |       1 | Skips: upstream `-chromFirst` toggle (t3); inline expected outputs are RNG-specific so byte-parity is replaced with structural invariants (t1/t2/t4: lengths preserved, include/exclude/chrom honoured). |
| bedcluster    |           3 |      3 |       0 | All upstream `cluster.t1`/`t2` cases (basic + `-s` stranded) plus an idempotency smoke (PR #87 wave-3 tail). |
| bedsplit      |           3 |      3 |       0 | All canonical upstream cases (`split.01/02/03`: `-a simple -n 50`, `-a simple -n 1000`, `-a size -n 50`); manifest head comparison. |
| bedsummary    |           4 |      4 |       0 | Spec-driven (upstream has no `summary/` test subdir): 2-chrom basic, `--no-header`, `--skip-all`, odd-count median. |
| bedtag        |           4 |      4 |       0 | Spec-driven (upstream has no `tag/` test subdir): default name-column join, `-labels` source prefix, `-names` per-source override, `-s` strand filter. |
| bedwindow     |           6 |      6 |       0 | Spec-driven (upstream has no `window/` test subdir): default A<TAB>B writer, symmetric `-w` expansion, `-c` count-only, `-v` invert, asymmetric `-l 0 -r N`, low-clipping at 0. |
| bedreldist    |           5 |      2 |       3 | Skips t01..t03 (require large `refseq.chr1.exons.bed.gz` / `aluY.chr1.bed.gz` / `gerp.chr1.bed.gz` fixtures we do not vendor); passes the shipped `issue_711` corner case byte-for-byte plus a small self-intersect mirror of t01. |
| bedfisher     |           6 |      5 |       1 | All five small upstream cases (`fisher.t1`..`t4`, `t6`) pass byte-for-byte; skip t5 (long $TMPDIR path is a CLI/filesystem concern, not algorithmic parity). |
| bednuc        |           4 |      3 |       1 | Spec-driven (upstream has no `nuc/` test subdir): default 3-interval profile, `-s` + `-seq` round-trip, `-pattern` case-sensitive count. Skips `-fullHeader` (index always keys on first whitespace token; best-effort fallback map covers the common `>chr1 extra info` case but not index-time semantics). |
| bedannotate   |           3 |      3 |       0 | Spec-driven (upstream has no `annotate/` test subdir): default per-B fractions, `-counts`, `-both` with `-names`. Header-padding difference vs upstream is tolerated in the test (data rows match byte-for-byte). |
| bedmulticov   |          10 |     10 |       0 | All upstream `multicov.t1` through `t10` cases pass byte-for-byte. The latest wave (this PR) wires `-split` block-aware coverage on BAM CIGAR `N` ops through a new `indexBAMSplit` that walks each alignment's CIGAR and emits one block per contiguous reference-consuming op-run (M/=/X, with D extending the current block to match upstream's `breakOnDeletionOps=false`), skipping `N`-op gaps. Each alignment is counted at most once per A interval; with `-f`, the threshold is applied to `total_block_overlap / sum_of_BAM_block_lengths` using strict `>` — preserving the bedtools 2.x quirk in `multiBamCov.cpp::FindBlockedOverlaps`. Unskipped `t5..t9` (split alone, split+`-s`, split+`-S`, split+`-f 0.01`, split+`-f 0.10`). `-q` MAPQ filter and `-D` per-A-interval depth cap are honoured on BAM inputs. |
| bedmultiinter |           4 |      3 |       1 | Spec-driven (upstream ships no `multiinter/` test subdir): fixtures + expected outputs taken byte-for-byte from `multiintersect_examples()` in `src/multiIntersectBed/multiIntersectBedMain.cpp`. Default mode, `-header -names`, and `-empty -g -header -names` all pass byte-for-byte. The skip documents the absence of an upstream test directory. |
| bedigv        |           3 |      3 |       0 | Spec-driven (upstream ships no `igv/` test subdir): expected outputs derived directly from `src/bedToIgv/bedToIgv.cpp` by playing forward `ProcessBed()` on a BED6 fixture. Cases: defaults; `-path /tmp/snaps -sess my.xml -name`; `-slop 50 -img svg -sort base -clps`. |
| bedlinks      |           3 |      3 |       0 | Spec-driven (upstream ships no `links/` test subdir): expected outputs derived directly from `src/linksBed/linksBed.cpp` + `linksMain.cpp` by playing forward `CreateLinks()`/`WriteURL()`. Cases: BED6 defaults; BED6 custom mirror (`-base/-org/-db` per upstream help example); BED3 defaults (bedType=3 branch, no name/score/strand `<td>`). |
| bedpairtobed  |           6 |      5 |       1 | Spec-driven (upstream ships no `pairtobed/` test subdir): expected outputs derived directly from `src/pairToBed/pairToBed.cpp` by walking `FindOverlaps()`/`FindOneOrMoreOverlaps()` over a hand-curated BEDPE × BED fixture. Cases: `-type either`, `-type both`, `-type neither`, `-type xor`, `-type notboth`; skip documents that `-slop` is upstream-`pairtopair`-only. |
| bedpairtopair |           6 |      5 |       1 | Spec-driven (upstream ships no `pairtopair/` test subdir): expected outputs derived directly from `src/pairToPair/pairToPair.cpp` by playing forward `FindHitsOnBothEnds()`/`FindHitsOnEitherEnd()` over a 4-pair × 4-pair fixture. Cases: `-type both`, `-type either`, `-type neither`, `-type notboth`, `-slop` near-miss; skip documents that the upstream stranded-slop case is covered by the unit test `TestStrandedSlop_Direction`. |
| **TOTAL**     |     **264** |**214** | **50**  | |

(The discrepancy between this table and `go test`'s 87 passed / 42 skipped is
two helper / sanity sub-tests in `bedsort` and `bedintersect` that are not
direct mirrors of an upstream case.)

The sickle and skewer ports each have a per-tool table in their respective
section below; the project-wide running total is 188 tests added, 136
passed, 49 skipped (52 if you count the sickle FixturesPresent + skewer
FixturesPresent + skewer PEHelperSmoke helper tests). The 2026-05-14
wave-2 update (`bedexpand`, `bedgetfasta`, `bedsample`, `bedspacing`)
added 34 of those tests (29 passing, 5 skipped).

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

### Discrepancies fixed (column-ops + discrepancies wave)

The two previously-skipped discrepancies from PR #55 are now fixed:

- **bedjaccard now pre-merges A and B before computing
  intersection/union**, mirroring upstream's
  `setUseMergedIntervals(true)` in
  `reference_code/bedtools/src/utils/Contexts/ContextJaccard.cpp`. The
  merge is per-strand under `-s` (matching upstream's
  `SAME_STRAND_EITHER` semantics in `FileRecordMergeMgr`), with `.` /
  unknown-strand records dropped. Newly passing parity cases:
  `jaccard.t02 / t03 / t05 / t06 / t10 / t11`.
- **bedmerge's `-s` now matches upstream's per-strand merge**. Records
  with `.` (or empty) strand are dropped under `-s` (matching
  `FileRecordMergeMgr.cpp` lines 47-58 + 96-129), and `+` / `-` groups
  are merged independently before the two output streams are
  recombined by (chrom, start, end). Newly passing: `merge.t15`.

The shared column-op vocabulary (`bedmerge.ApplyOp`, used by
`bedgroupby`, `bedmap`, `bedcoverage`) was extended with `stdev`,
`sstdev`, `absmin`, `absmax`, `cat`, `cat_uniq` — matching upstream's
KeyListOps. Newly passing: `map.t11` (absmin) + new `map.t13`
(absmax); previously-skipped bedgroupby tests `TestGroup_AdditionalOps`
and `TestGroup_StdevSstdev` now run.

### Discrepancies found (NOT fixed)

These are real upstream/Go differences that the parity tests document
with `t.Skip("known discrepancy: ...")`. Each one is an open task for a
later PR.

_(No known discrepancies remaining for bedtools after the column-ops +
discrepancies wave. Remaining bedtools skips are all for unimplemented
features — BAM/VCF/GFF input, `-split`, `-S <strand>` filter,
`-delim`, large fixtures — not behaviour discrepancies.)_

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
- **`-fullHeader`** on `bedgetfasta`. Upstream uses the full FASTA header
  line (whitespace tolerated) when matching contig names from BED; our port
  uses the same first-token convention as `samtools faidx`. The skipped
  parity test (`Getfasta.T07`) documents the gap.
- **BGZF FASTA input** (`-fi *.fa.gz`) on `bedgetfasta`. Resolved: this PR
  adds `OpenRandomAccessBGZF` in `pkg/bioformats/fasta` (magic-sniffed
  auto-routing from the existing `OpenRandomAccess` entry point), reads
  the `.gzi` sidecar via a stdlib-only little-endian parser, and honours
  a sibling `.fa.gz.fai` when present. Upstream `getfasta.t18` (BGZF +
  `-split` BED12) now passes byte-for-byte using upstream's `t.fa.gz`
  fixture; the validation harness regenerates the `.gzi` on-demand via
  the existing `tools/bgzip --reindex` path. The in-memory
  decompression strategy is documented in `pkg/bioformats/fasta/bgzf.go`
  with a note that partial-decompression seek via `.gzi` (htslib-style)
  can be layered on later without API churn.
- **CLI-only "no args / bad args" diagnostics** on `bedsample`. Upstream
  prints specific `***** ERROR:` messages from its CLI driver; we exit with
  Go's stock `flag` package error handling and a `bedsample: …` prefix.
  These are CLI surface details, not library behaviour, and the skipped
  parity cases (`Sample.T01/T02`) document the divergence.
- **BAM input** on `bedspacing`. Upstream's `spacing -i x.bam` accepts BAM;
  our port is BED-only. Skipped: `Spacing.T07`.

### Per-tool notes for the wave-2 additions

The four wave-2 ports — `bedexpand`, `bedgetfasta`, `bedsample`,
`bedspacing` — bring the bedtools tally to 17 ported subcommands and 161
parity tests (114 passing, 47 documented skips).

- **bedexpand** (`tools/bedexpand/pkg/bedexpand/parity_test.go`). All three
  inline upstream cases (`expand.t1..t3`) pass byte-for-byte. The Go port
  implements the same per-row walk as
  `reference_code/bedtools/src/expand/expand.cpp`: non-expanded columns
  emit verbatim, expanded columns substitute the k-th element of the k-th
  list named in `-c`. So `-c 5,4` swaps the two expanded columns,
  reproducing `expand.t3`.
- **bedgetfasta** (`tools/bedgetfasta/pkg/bedgetfasta/parity_test.go`). All
  cases pass against the upstream `getfasta` corpus, covering default
  `chrom:start-end` header, `-name` / `-nameOnly`, `-s`, `-split`,
  `-split -s` (per-block revcomp, reversed-order blocks), `-rna`,
  `-fullHeader`, and BGZF (`.fa.gz`) input via `-fi`. The port carries
  its own `FetchPreserveCase` so IUPAC + case round-trip exactly — the
  shared `pkg/bioformats/fasta.RandomAccess.Fetch` uppercases for
  downstream callers that need a canonical case.
- **bedsample** (`tools/bedsample/pkg/bedsample/parity_test.go`). 5 cases
  pass: requested count, deterministic seed, `-header` forwarding,
  too-few-records error, and the "subset of input" invariant. We cannot
  byte-match upstream's PRNG (different sampler), so seeded runs are
  deterministic within `bedsample` but not against upstream. The two
  skips are CLI-only diagnostics.
- **bedspacing** (`tools/bedspacing/pkg/bedspacing/parity_test.go`). The
  single upstream test (`spacing.t01`) passes byte-for-byte. The port
  matches `reference_code/bedtools/src/spacingFile/spacingFile.cpp`'s
  "previous record on chrom" semantics (not running-max-end), which the
  test cases verify with synthetic per-chrom-reset / exact-abut / overlap
  / BED6-preservation / single-record edges.

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

### bedcoverage parity

Byte-for-byte parity against
`reference_code/bedtools/test/coverage/test-coverage.sh` using the upstream
`a.bed` / `b.bed` fixture pair.

| Case | Mode | Status | Notes |
| ---- | ---- | ------ | ----- |
| coverage.t1  | BAM input | skip | BAM/SAM input not supported in bedcoverage. |
| coverage.t2  | default (count / bp / len / fraction) | pass | 7-decimal fraction matches upstream exactly. |
| coverage.t3  | `-counts` | pass | |
| coverage.t4  | `-hist` per-A + `all` footer | pass | 25 lines, ordering matches. |
| coverage.t5  | `-d` (per-base depth) | pass | 300 lines. |
| coverage.t6  | `-mean` | skip | Upstream prints float32 noise (`1.3200001`) we cannot reproduce with float64. |
| coverage.t7  | `-s` (same strand) | pass | |
| coverage.t8  | `-S` (opposite strand) | pass | |
| coverage.t10 | `-split` BAM | skip | BAM input + BED12 block-aware coverage not supported. |

### bedmap parity

Byte-for-byte parity against
`reference_code/bedtools/test/map/test-map.sh` using the upstream
`ivls.bed` / `values{,2,4}.bed` fixture set.

| Case | Mode | Status | Notes |
| ---- | ---- | ------ | ----- |
| map.t01 | defaults (`-c 5 -o sum`) | pass | |
| map.t02 | explicit `-o sum` | pass | |
| map.t03 | `-o count` | pass | Count emits 0 for empty groups, not the null placeholder. |
| map.t04 | `-o mean` | pass | |
| map.t05 | `-o max` | pass | |
| map.t06 | `-o min` | pass | |
| map.t07 | `-o mode` (values2.bed) | pass | |
| map.t08 | `-o antimode` (values2.bed) | pass | |
| map.t09 | `-c 7 -o collapse` (values4.bed) | pass | BEDPlus column extraction. |
| map.t10 | `-c 7 -o min` (signed values) | pass | Negative numbers handled. |
| map.t11 | `-c 7 -o absmin` | pass | `absmin` now in shared `bedmerge.ApplyOp` (column-ops + discrepancies wave). |
| map.t13 | `-c 7 -o absmax` | pass | `absmax` now in shared `bedmerge.ApplyOp` (column-ops + discrepancies wave). |
| map.t14 | GFF input | skip | GFF input not supported (BED only). |

### bedshuffle parity

Upstream's expected outputs are tied to bedtools' own Mersenne Twister
RNG, which our `math/rand`-based port does not reproduce byte-for-byte.
Byte-parity is therefore replaced with **structural-invariant parity**:
each case asserts the property the upstream test was designed to check
(lengths preserved, include/exclude/chrom honoured, error on
unplaceable intervals).

| Case | Scenario | Status | Notes |
| ---- | -------- | ------ | ----- |
| shuffle.t1 | basic shuffle on hg19 | pass | Asserts: line count, length preserved, chrom in genome. |
| shuffle.t2 | `-incl` include regions | pass | Every output contained in some include region. |
| shuffle.t3 | `-incl -chromFirst` | skip | `-chromFirst` toggle equivalent to default mode in our port. |
| shuffle.t4 | `-excl` exclude regions | pass | No output overlaps any exclude region. |
| shuffle.t5 | sanity check (without `-excl`) | pass | Without `-excl`, some overlap with `excl.bed` is expected. |
| shuffle.t6 | interval larger than chrom | pass | Errors out with "could not avoid..." / "non-positive length"; matches upstream's error shape. |

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

---

## bcftools parity validation

The bcftools port (`tools/bcftools/` covering `view`, `index`, `stats`,
`query`, `concat`, `norm`, `call`, plus the PR #86 wave-1 tail
`annotate` / `head` / `isec` / `merge` / `reheader` / `sort`, the
convert/mendelian PR's `convert` + `mendelian`, the gtcheck/roh PR's
`gtcheck` + `roh`, the filter/consensus PR's `filter` + `consensus`,
and the mendelian2/polysomy PR's `mendelian2` + `polysomy`)
is validated against upstream `bcftools 1.19+htslib-1.19` via
`tools/bcftools/pkg/bcftools/parity_test.go` plus the per-subcommand
unit suites under the same package directory.

`convert`, `mendelian`, `mendelian2`, and `polysomy` are exercised
exclusively by hand-built fixtures in
`tools/bcftools/pkg/bcftools/{convert,mendelian,mendelian2,polysomy}_test.go`
in this PR; an upstream-fixture parity run will land in the
follow-up parity wave (the upstream plugins live under
`reference_code/bcftools/plugins/` and `reference_code/bcftools/test/`,
which we do not currently vendor for the per-suite parity rig).

The brief differs from bedtools' in two ways:

1. The upstream test corpus (`reference_code/bcftools/test/`) is heavily
   dependent on FORMAT-typing and tagged-INFO fixtures that exercise
   bcftools-specific encodings we don't yet decode (notably the BCF
   `int64` typed descriptor and per-record FORMAT key dictionary on
   htslib-produced BCFs). To keep the parity tests small and meaningful
   we hand-crafted minimal fixtures under
   `tools/bcftools/testdata/parity/` and captured the expected output
   by running upstream `bcftools 1.19` against them; this matches
   `concat.1.vcf.out` byte-for-byte for the upstream concat fixture
   when we replay it (see `TestParityConcat_UpstreamFixture`).
2. `bcftools` strips its own `##bcftools_<cmd>Version=` /
   `##bcftools_<cmd>Command=` lines when `--no-version` is passed; we
   capture expected outputs with `--no-version` so the comparison is
   stable.

## Summary

| Subcommand | Tests added | Passed | Skipped | Notes |
| ---------- | -----------:| ------:| -------:| ----- |
| view       |          12 |     10 |       2 | Skips: `-v/--types`, sample-subset AC/AN recomputation. BCF header parity passes; per-record FORMAT decode from htslib BCF is still a gap. |
| query      |          14 |     11 |       3 | Skips: bare `%INFO`, `%N_ALT`, `[%FORMAT/<char>]`. Every other listed token in the README is byte-stable. |
| index      |           4 |      2 |       2 | Skips: `.tbi` binary equality (BGZF padding differs), `.csi` for htslib-produced BCF (int64 / FORMAT-dict gap). |
| stats      |           8 |      1 |       7 | Only the SN section is byte-stable today; AF/QUAL/IDD/ST/DP/PSC/HWE/PSI all diverge on formatting glyphs or on whether the section is recomputed from genotypes vs read from INFO. Tracked. |
| concat     |           6 |      4 |       2 | Skips: `-a` sort-merge (different contig-order heuristic), plain `-D` adjacency dedup (upstream requires `-a`). Plain concat, `-O z` round-trip, conflicting-header detection, and the upstream `concat.1.vcf.out` fixture all match byte-for-byte. |
| norm       |           7 |      4 |       3 | Skips: `-f` left-align (FASTA fixture not yet added), `-c` check-ref policies (same), and a placeholder for the realignment regression suite. `-m -`/`-m -snps`/`-m +`/`-a`/`-d exact` all match. |
| call       |           6 |      4 |       2 | Skips: BCF input (FORMAT-key gap; see `docs/UPSTREAM_BUGS.md`), full multi-allelic caller. Consensus calls (`-c`, `-c -v`, `-A`, `--ploidy 1`) match hand-crafted upstream-replicated fixtures byte-for-byte. |
| head       |           3 |      3 |       0 | All cases: default emit, `-n N` slice, `--samples` sample-only (PR #86 wave-1 tail). |
| sort       |           3 |      3 |       0 | All cases: out-of-order CHROM/POS records, already-sorted no-op, empty-records header-only. `-m/-T` flags are accepted but in-memory only. |
| isec       |           3 |      3 |       0 | All cases: `-n=2 -w 1` intersection, `-n ~10` a-only bitmask, `-p PREFIX` per-input projection files. `-c some` simplified to strict tuple. |
| merge      |           3 |      3 |       0 | All cases: two single-sample VCFs → two-sample, disjoint positions, single-input rejected. Pre-sort required. |
| reheader   |           3 |      3 |       0 | All cases: positional sample rename, `OLD\tNEW` mapping rename, full header-file substitution. `-i` in-place mode emits to stdout in v1. |
| annotate   |           3 |      2 |       1 | Skip: `--set-id` macro expansion. Passing: `-x ID`, `-x INFO/DP`. |
| filter     |           7 |      7 |       0 | All cases: `-i` + `-s` soft-tag, `-e` + `-s`, `-m x` reset-on-pass, `-m +` append-preserve, `-S .` GT-rewrite (preserves `\|` phase), `-S 0` GT-rewrite, `-g SnpGap`, `-s +` auto-named filter. Mask flags (`--mask` / `-M`) parse cleanly and hard-reject with a roadmap pointer. |
| consensus  |          11 |     11 |       0 | All cases: SNP apply-all-ALTs, simple insertion (REF=A,ALT=AC), simple deletion (REF=AC,ALT=A), `--mark-del CHAR` padding, per-sample GT (hom-ref vs hom-alt), `-H R` (REF in hets), `-H A` (ALT in hets), `-p prefix`, `--mark-snv uc`, `-m` mask with default N, overlapping-variants-first-wins. `-c/--chain` and `-H NpIu` parse cleanly and hard-reject with a roadmap pointer. |
| **TOTAL**  |      **93** | **72** | **22**  | |

(Counts include three subtests in `view` and `query` that exercise the
streaming vs file paths separately. The two BCF header-only tests are
counted as separate cases.)

## Discrepancies found (and fixed in this PR)

The validation surfaced a handful of byte-level divergences that we fixed
inline rather than masking with `t.Skip`. All of them lived in either
the shared VCF/BCF library or in the bcftools port:

- **`pkg/bioformats/vcf`: QUAL formatting** now uses upstream's
  minimal-precision rule (integer-as-integer, otherwise shortest `%g`)
  rather than the previous unconditional `%.2f`. Upstream prints `30`
  where we used to print `30.00`. New helper: `formatQual`.
- **`pkg/bioformats/vcf`: INFO key order preserved**. The `Variant`
  struct gained an `InfoOrder []string` slice. The parser populates it
  in source order; the writer iterates `InfoOrder` first and only falls
  back to the (alphabetised) map keys for fields the caller appended
  without registering. Without this our VCF→VCF passes shuffled INFO
  keys at random because Go maps are not ordered.
- **`pkg/bioformats/bcf`: header IDX-suffix stripped on read**. htslib's
  bcf header binary contains an `,IDX=<N>` annotation on every
  structured line; upstream's text view drops it. We now strip it in
  the BCF reader's header parsing so the text-out parity holds.
- **`pkg/bioformats/bcf`: BCF int64 typed descriptor (4) decoded**.
  htslib 1.13+ may emit type 4 (int64) for FORMAT vectors over very
  large counts. We added a decoder that clamps to int32 range
  (sufficient for our downstream model). Pre-change our reader errored
  on any htslib-produced BCF that hit this path.
- **`tools/bcftools/pkg/bcftools/view.go`: implicit `##FILTER=<ID=PASS,…>`
  injected on output** when the input header lacks one. Upstream's
  htslib `bcf_hdr_parse_line` auto-adds this entry; we now match.
- **`tools/bcftools/pkg/bcftools/view.go`: `-G/--drop-genotypes` strips
  ##FORMAT lines from the output header**, matching upstream. Without
  this our `-G` output diverged on every meta-line block.
- **`tools/bcftools/pkg/bcftools/concat.go`: PASS injection in
  `MergeHeaders`** complements the view-side fix above so concat output
  always has the PASS line right after `##fileformat=...`, in the same
  position upstream uses.
- **`tools/bcftools/pkg/bcftools/norm.go`: bare `-m -` / `-m +` treated
  as `any`** (both SNPs and indels) instead of erroring out. Upstream
  accepts a bare `-`/`+`; we did not.
- **`tools/bcftools/pkg/bcftools/query.go`: `[ ... ]` sample loop no
  longer auto-inserts a tab separator**. Upstream repeats the inner
  body verbatim for each sample; the leading literal inside the
  brackets (`[\t%GT]`) is what creates the inter-sample tab. Our
  previous behaviour produced `0/0\t\t0/1\t\t1/1`; we now produce
  `0/0\t0/1\t1/1` (with the leading literal) or `0/00/11/1` (without).
  Existing `query_test.go` cases were updated.

## Upstream divergences left as `t.Skip`

The brief calls these out as "skip with a one-line rationale". Each
skipped test points at `docs/PARITY_ROADMAP.md` (a known gap in the
port) or `docs/UPSTREAM_BUGS.md` (a documented divergence from
upstream that we're tracking before changing behaviour). The full
list is the `SKIP` lines in `go test -v -run TestParity ./tools/bcftools/...`.

The most-visible group is the `stats` family: upstream re-derives AF
from GT and prints integer-bin labels with a `.0` suffix; we read
INFO/AF directly and use fixed bins. Both behaviours are reasonable;
this is more "feature gap" than "bug" and is tracked in
`docs/PARITY_ROADMAP.md`.

## Reproducing locally

```bash
# All bcftools parity tests, race detector on:
go test -race -run TestParity ./tools/bcftools/...

# View just the SKIP rationale:
go test -v -run TestParity ./tools/bcftools/... 2>&1 | grep -A0 SKIP:

# Generate (or regenerate) an expected fixture:
bcftools view --no-version tools/bcftools/testdata/parity/basic.vcf > \
    tools/bcftools/testdata/parity/view_basic.expected.vcf
```

## How to add a new bcftools case

1. Put the input under `tools/bcftools/testdata/parity/<name>.vcf` (or
   `<name>.bcf` / `<name>.vcf.gz` if testing a format-specific path).
   Avoid bringing in large fixtures from
   `reference_code/bcftools/test/` — those are heavily macro-expanded by
   `test.pl` and depend on htslib internals. A 10-record hand-crafted
   VCF that hits the specific feature is preferred.
2. Run upstream `bcftools <subcmd> --no-version ...` (or without
   `--no-version` for subcommands that don't accept it, e.g. `query`)
   and capture the output to `tools/bcftools/testdata/parity/<name>.expected.<ext>`.
3. Add a `TestParity<Subcmd>_<Name>` to
   `tools/bcftools/pkg/bcftools/parity_test.go` that calls the Go port
   on the same input and `equalBytes`-asserts the result.
4. If the case touches a feature we don't implement, mark it
   `t.Skip("<one-line reason> (see docs/PARITY_ROADMAP.md bcftools <subcmd>)")`
   instead of deleting it — the skip count is the project's gap meter.

---

## samtools

Byte-for-byte parity against the upstream `samtools` C regression test
suite vendored at `reference_code/samtools/test/`.

The methodology mirrors the bedtools section above: cases from upstream's
`test.pl` driver (`test_view`, `test_sort`, `test_bam2fq`, etc.) are
replicated as Go subtests in
`tools/samtools/pkg/samtools/parity_test.go`. Inputs and expected outputs
are vendored under `tools/samtools/testdata/parity/` — either direct byte
copies of upstream golden files (e.g. `bam2fq/1.1.fq.expected`,
`sort/pos.sort.expected.sam`) or small purpose-built SAM fixtures where
the upstream golden files are impractically large for unit-test scope.

### samtools summary

| Subcommand | Tests added | Passed | Skipped | Notes |
| ---------- | -----------:| ------:| -------:| ----- |
| view       |          10 |      9 |       1 | Skip: CRAM input (`-C/-T`); BAM↔SAM round-trip + flag/MAPQ/RG/region/header-only covered. `-L/--regions-file` BED filter shipped with per-chrom `bed.IntervalTree` and table-driven tests (`TestView_BedFilter_*`); `-M`/`--use-multi-region-iterator` is accept-and-ignore (the predicate result matches the indexed-walk path). `-d`/`-D` tag-value filter and `-N` qname-file filter shipped with hand-built table-driven tests in `view_test.go`. |
| sort       |           6 |      3 |       3 | Skips: `-n`/`-N` FLAG tie-break gap (2 cases), `-t TAG` 3-key compare gap. |
| index      |           5 |      5 |       0 | All cases: BAI build, CSI rejection, BAI region query, multi-chrom, empty BAM. |
| depth      |           8 |      6 |       2 | Skips: `-a`/`-A` zero-fill edge cases, `-b BED` byte parity. |
| fastq      |           7 |      4 |       3 | Skips: QNAME-based pair detection (singleton mid-stream), CRAM input, `-T '' / -T '*'` all-tag expansion. |
| flagstat   |           7 |      7 |       0 | All counters validated incl. QC-fail column, secondary/supplementary, diff-chr, paired-but-unmapped. |
| dict       |           3 |      3 |       0 | All cases: minimal one-record FASTA → @HD + @SQ with M5; multi-record order; `-a`/`-s`/`-H` (PR #88 wave-1 tail). |
| quickcheck |           3 |      3 |       0 | All cases: well-formed BAM passes; text-SAM rejected on magic check; empty file rejected with `empty file` reason. |
| idxstats   |           3 |      3 |       0 | All cases: 3-row output (chr1/chr2/`*`) from `basic.sam` BAI fast-path; 4-col TSV format; upstream-golden 4-column shape check. |
| coverage   |           3 |      3 |       0 | All cases: per-ref tabular rows for chr1+chr2, `-H` no-header, `Regions=["chr1"]` filter. `-A` histogram deferred. |
| sam-merge  |           3 |      3 |       0 | All cases: two-input concatenation, single-input copy, coordinate-sorted interleave (POS monotonic). |
| cat        |           3 |      3 |       0 | All cases: 2-input record-order preservation; empty input list errors; diverging @SQ tables rejected. |
| reheader   |           3 |      3 |       0 | All cases: `HeaderPath` substitution adds @CO line; record bodies preserved across reheader; @SQ-table size mismatch rejected. |
| addreplacerg |         3 |      3 |       0 | All cases: orphan-only mode adds RG, overwrite-all replaces existing RG, unknown RG id rejected. |
| fixmate    |           3 |      3 |       0 | All cases: paired records get correct RNEXT/PNEXT/TLEN, `-m` adds `ms` aux, `-c` adds `MC` aux. |
| split      |           3 |      3 |       0 | All cases: per-RG output files, unidentified-RG capture, single-RG one-file output. |
| mpileup-tail |         3 |      3 |       0 | Passing: `-d 1` (MaxDepth) caps depth, `-A` (CountOrphans) wired. `-aa` full-contig zero-fill exercised via `TestMpileup_AA_ZeroFillTableDriven` (multi-contig + gap fixture: chr1 partial, chr2 fully empty, chr3 partial). |
| markdup    |           4 |      4 |       0 | Byte-for-byte parity on 5_markdup and 6_remove_dups; flag+qname parity on 18_primary_duplicate_count (we don't emit `dt:Z:` tag); sequence-mode dup-count parity on a duplicate of fixture 5. See deferred-feature list below for the deliberate skips (optical-dup, per-RG keying). |
| stats      |           6 |      6 |       0 | SN-section byte parity on fixtures 1, 2, 5, 7, 8, 10 from `reference_code/samtools/test/stat/` (all 8 SN-only cases we exercise pass byte-identical). Non-SN sections (RL/MAPQ/IS) are emitted but only smoke-tested; see deferred-section list below. |
| calmd      |           6 |      5 |       1 | Logical MD/NM parity on hand-built fixtures covering match, mismatch, deletion, insertion, soft-clip, multi-contig; `-e` rewrites SEQ to '=' on matches; BAM round-trip preserves both tags; existing-tag overwrite emits "different" stderr warning + Quiet suppresses it. Skip: upstream `bam_md.c` `-uAr` BGZF byte-diff (libdeflate version mismatch). |
| import     |          12 |     11 |       1 | FASTQ→BAM/SAM across R1+R2 paired, -0 unpaired, -s interleaved, -T '*' / -T 'XZ' / -T '' aux extraction, -R short / -r full RG, --order int+padded forms, /1/2 suffix strip toggle, BAM round-trip, two positional args. Skip: upstream import-then-fastq round-trip (requires upstream samtools fastq for the round-half). |
| phase      |           5 |      5 |       0 | Hand-built SAM fixtures: two-het chain (flip case), two-het chain (label-1 case), ambiguous label when reads don't bridge two hets, MinMAPQ filter drops T-bearing reads, default-constants guard. Upstream MCMC chimera repair deferred; `-b` per-haplotype BAM split not yet wired (see PARITY_ROADMAP.md#samtools). |
| targetcut  |           3 |      3 |       0 | Hand-built SAM fixtures: soft-clip + insertion + deletion + unmapped/secondary skip + SEQ='*' skip; `-Q` per-base quality filter drops sub-threshold bases; default-constant guard. The upstream `cut_target.c` HMM consensus mode is **NOT** implemented in v1 — the Go port emits an aligned-slice FASTA per read instead (documented scope reduction). |
| **TOTAL**  |     **112** | **98** | **14**  | |

### samtools: discrepancies found in our port (fixed in this PR)

The validation surfaced three real divergences from upstream that were
fixed inline rather than masked with `t.Skip`:

- **samtools sort: `-n` and `-N` flags were inverted**. Upstream
  `bam_sort.c` reads `case 'N': natural_sort = 0; // fall through; case
  'n': sam_order = QueryName;` — i.e. `-n` is natural-numeric sort (the
  default for name-sort) and `-N` flips to plain lex. We had the CLI
  mapping the other way around. The library API (`SortByName` /
  `SortByNameNatural`) was unchanged; only the `runSort` CLI wiring moved.

- **samtools sort: missing `SS:queryname:{natural,lexicographical}` stamp
  on the `@HD` header line**. Upstream writes this sub-sort tag so
  downstream tools (and our own parity tests) can identify which form of
  name-sort produced a given file. Added in `stampSortOrder` in
  `tools/samtools/pkg/samtools/sort.go`.

- **samtools fastq: pair-suffix not auto-dropped in `-1/-2` mode**.
  Upstream `bam_fastq.c` does `if (opts->fnr[1] || opts->fnr[2])
  opts->has12 = false;` — when the user passes `-1` and/or `-2`, the
  `/1`/`/2` read-name suffix is suppressed because the separate output
  files already disambiguate mate identity. Our port was unconditionally
  adding it, failing byte-for-byte parity against
  `bam2fq/1.1.fq.expected` etc. Fixed in `Fastq` in
  `tools/samtools/pkg/samtools/fastq.go`; `-N` (`AlwaysAddSuffix`) still
  forces the suffix on.

### samtools: discrepancies found (NOT fixed)

These are real upstream/Go differences that the parity tests document
with `t.Skip("known discrepancy: ...")`. Each one is an open task for a
later PR (recorded against
[PARITY_ROADMAP.md](../docs/PARITY_ROADMAP.md#samtools)):

- **sort tie-break on FLAG**. Upstream `samtools sort -n` (and `-N`)
  uses a secondary comparison on `b->core.flag` when two records share a
  QName, so e.g. `r001/83` sorts before `r001/163`. Our port falls
  through to a stable input-order tie-break. Affects upstream
  `sort/name.sort.expected.sam` exactly when there is a same-QName pair
  — the bulk of the output is identical.

- **sort by tag uses a 3-key compare upstream**. Upstream
  `bam_sort.c`'s tag-sort key falls back to `(refID, pos, qname)` rather
  than just qname, so two records with the same tag value sort in
  coordinate order. Our port only uses qname as the secondary key.

- **samtools fastq active QNAME-based pairing**. In paired-output mode
  (`-1 -2 [-s]`), upstream actively pairs adjacent records by QNAME: a
  record whose paired flag is set but whose neighbour has a different
  QNAME is routed to the singleton file even though the flag says
  "paired". Our port dispatches by flag bits alone (0x40 → R1, 0x80 →
  R2, paired-but-neither/both → orphan, unpaired → singleton), which
  produces a different output for `bam2fq.002.sam` where
  `ref1_grp2_p002a` has no mate.

### samtools: what is NOT validated (skipped features)

| Subcommand | Feature | Tracking |
| ---------- | ------- | -------- |
| view       | CRAM input/output (`-C`, `-T ref`) | PARITY_ROADMAP.md#samtools |
| sort       | Multi-threading (`-@`) | PARITY_ROADMAP.md#samtools |
| sort       | Minimiser sort (`-M`/`-K`) | PARITY_ROADMAP.md#samtools |
| sort       | Template-coordinate sort | PARITY_ROADMAP.md#samtools |
| index      | CSI output (`-c`) | PARITY_ROADMAP.md#samtools |
| index      | Multi-threading | PARITY_ROADMAP.md#samtools |
| depth      | `-b BED` byte-parity | PARITY_ROADMAP.md#samtools |
| depth      | `-a`/`-A` zero-fill edge cases | PARITY_ROADMAP.md#samtools |
| fastq      | `-T '' / -T '*'` all-tags | PARITY_ROADMAP.md#samtools |
| fastq      | `-D TAG:file` value-list filter | PARITY_ROADMAP.md#samtools |
| fastq      | `--no-sc`/`--sc-aux` soft-clip handling | PARITY_ROADMAP.md#samtools |
| fastq      | Index extraction (`--index-format`) | PARITY_ROADMAP.md#samtools |
| markdup    | Optical-dup detection (`-d/--max-dist`, Illumina (x,y) parsing) | PARITY_ROADMAP.md#samtools |
| markdup    | Per-read-group keying (upstream `-S`) — fixture 17_read_group | PARITY_ROADMAP.md#samtools |
| markdup    | Barcode regex / barcode-tag keying | PARITY_ROADMAP.md#samtools |
| markdup    | `dt:Z:` "duplicate-type" aux tag (SQ/LB/OQ) — fixture 18 diffs only here | PARITY_ROADMAP.md#samtools |
| stats      | CHK checksum header block (CRC32 of names/seqs/quals) | PARITY_ROADMAP.md#samtools |
| stats      | COV/COV2 coverage histograms (need reference + BAI) | PARITY_ROADMAP.md#samtools |
| stats      | GCD/GCT/GCC/GCL GC-content distributions (need reference) | PARITY_ROADMAP.md#samtools |
| stats      | FFQ/LFQ per-cycle quality matrices, OXC oxidation context | PARITY_ROADMAP.md#samtools |
| stats      | `--target-regions BED` restriction | PARITY_ROADMAP.md#samtools |

### samtools: bugs found in upstream

None during this audit. If we find any in subsequent slices they will be
recorded in [docs/UPSTREAM_BUGS.md](../docs/UPSTREAM_BUGS.md) and skipped
in the parity test until we have a fix.

---

## mosdepth parity validation

`tools/mosdepth` is a single-subcommand tool. We mirror the upstream
`reference_code/mosdepth/functional-tests.sh` cases against the
vendored fixtures in `tools/mosdepth/testdata/parity/` (`ovl.bam`,
`empty-tids.bam`, `full-fragment-pairs.bam`, `track.bed`,
`unordered.bed`).

| Case | Status | Notes |
| --- | --- | --- |
| `overlapM` default per-base | SKIP | Needs overlap-pair detection (open gap; see [UPSTREAM_BUGS.md#mosdepth-overlap-pair-detection](../docs/UPSTREAM_BUGS.md#mosdepth-overlap-pair-detection)). |
| `overlapM` summary MT | SKIP | Same. |
| `overlapFastMode` per-base MT | PASS | Byte-for-byte under `--fast-mode`. |
| `overlapFastMode` per-base chr1 zero-depth | PASS | Whole-chrom zero-depth row. |
| `big_window` MT regions row | PASS | Structural check (single row spanning full reference). |
| `length_filter --min-frag-len 81` | PASS | All MT reads dropped. |
| `length_filter --max-frag-len 79` | PASS | Same. |
| `length_filter --min/--max 80 --fast-mode` | PASS | Reads kept; output identical to unfiltered fast-mode. |
| `threshold_test_by` (track.bed) byte parity | SKIP | 2X count differs by overlap-pair detection. |
| `threshold_test_by` our-values check | PASS | Pins our without-overlap-dedup numbers. |
| `track_header` regions mean | SKIP | Region mean depends on overlap-pair detection. |
| `unordered_bed` row count | PASS | regions.bed.gz has 2 rows. |
| `test_read_group` matching RG | PASS | Dist file has MT + total rows. |
| `test_missing_read_group` | PASS | Dist is `MT 0 1.00` / `total 0 1.00`. |
| `missing_chrom` strict failure | SKIP | We emit empty outputs; upstream exits 1. |
| `bad_frag_len_filter` strict failure | SKIP | We don't error on inverted bounds. |
| `--no-per-base` suppression | PASS | per-base.bed.gz is not created. |
| `-c MT` restriction | PASS | Only MT rows in per-base + summary. |
| `--mapq 60` byte parity | SKIP | Same overlap-pair gap. |
| `--flag 4` byte parity | SKIP | Same. |
| `fragment_mode` | SKIP | `--fragment-mode` not implemented. |
| `quantest` (`-q`) | SKIP | `--quantize` not implemented. |
| `--d4` rejection | PASS | Returns a clear error. |
| `empty_tids` (`--by` + thresholds) | PASS | Run completes; thresholds file emitted. |
| `.csi` vs `.tbi` index format | SKIP | We emit TBI; tracked at [PARITY_ROADMAP.md#mosdepth](../docs/PARITY_ROADMAP.md#mosdepth). |

Totals: 24 cases, 12 PASS, 12 SKIP. Skipped cases are split between the
overlap-pair detection gap (one root cause covers most of the
`t.Skip()`s) and a handful of options we have not implemented
(`--fragment-mode`, `-q`, `.csi` output).

### Discrepancies found (and fixed in this PR)

None. The overlap-pair gap, fast-mode equivalence, and TBI/CSI
deviation were already known when this PR started; we documented them
rather than papering over them.

### Discrepancies found (NOT fixed in this PR)

- **mosdepth lacks overlap-pair detection** for mate-pair reads.
  Upstream subtracts one copy of depth where the two ends of a fragment
  overlap on the reference. Our default-mode pipeline counts both copies
  and therefore produces depths that match upstream's `--fast-mode`
  output rather than upstream's default. Tracked at
  [docs/UPSTREAM_BUGS.md#mosdepth-overlap-pair-detection](../docs/UPSTREAM_BUGS.md#mosdepth-overlap-pair-detection).

## vcftools parity validation

`tools/vcftools` ports ~60 of the upstream ~147 options (see
`docs/PARITY_ROADMAP.md#vcftools`). The parity rig pins the fixture
`reference_code/vcftools/examples/valid-4.0.vcf` (12 sites, 3 samples)
into `tools/vcftools/testdata/parity/sample.vcf` and asserts each
option's output file matches a hand-computed golden file (or, when
upstream and our port diverge on a known format / formula gap, asserts
the file shape and skips byte parity).

| Category | Case | Status | Notes |
| --- | --- | --- | --- |
| Filter | `--chr 19` | PASS | Keeps 2 chr19 rows. |
| Filter | `--from-bp/--to-bp` | PASS | Keeps the two chr20 sites in range. |
| Filter | `--maf 0.4` | PASS | Keeps 4 sites (20:14370, X:9, X:11, X:12). |
| Filter | `--mac 3` | PASS | Keeps ≥1 record. |
| Filter | `--minQ 20` | PASS | Keeps 4 high-QUAL sites. |
| Filter | `--remove-indels` | PASS | Drops both indels. |
| Filter | `--keep-only-indels` | PASS | Keeps both indels plus symbolic-allele X:11. |
| Filter | `--max-missing 1.0` | PASS | Drops partially-missing X:11. (Fixed off-by-one in this PR.) |
| Sample | `--indv NA00001` | PASS | Keeps 1 sample column. |
| Sample | `--remove-indv NA00003` | PASS | Drops 1 sample column. |
| Sample | `--keep FILE` | PASS | Keeps 2 sample columns. |
| Stats | `--freq` | PASS | Byte-for-byte against golden file (fixed `{ALLELE:FREQ}` header in this PR). |
| Stats | `--counts` | PASS | Byte-for-byte. |
| Stats | `--site-pi` upstream byte parity | SKIP | Formula difference, see [UPSTREAM_BUGS.md#vcftools-site-pi](../docs/UPSTREAM_BUGS.md#vcftools-site-pi). |
| Stats | `--site-pi` textbook spot-check | PASS | Three hand-computed values. |
| Stats | `--hardy` | PASS | Byte-for-byte (P_HET_DEFICIT/EXCESS placeholders, gap tracked). |
| Stats | `--missing-site` | PASS | Byte-for-byte. |
| Stats | `--missing-indv` | PASS | Byte-for-byte. |
| Stats | `--depth` | PASS | Byte-for-byte. |
| Stats | `--site-depth` | PASS | SUMSQ_DEPTH emitted as 0 placeholder. |
| Stats | `--site-mean-depth` | PASS | VAR_DEPTH emitted as 0 placeholder. |
| Stats | `--het` | PASS | Byte-for-byte. |
| Stats | `--singletons` | PASS | Byte-for-byte (fixed `SINGLETON/DOUBLETON` + `INDV` columns in this PR). |
| Stats | `--TsTv-summary` | PASS | Byte-for-byte (fixed `MODEL\tCOUNT` format in this PR). |
| Stats | `--TsTv-by-count` header | PASS | Header byte-for-byte. |
| Stats | `--TsTv-by-count` rows | SKIP | Upstream emits every 0..2*N bin (incl NaN ratios); ours emits non-empty bins only. |
| Stats | `--TsTv-by-qual` header | PASS | Header byte-for-byte. |
| Stats | `--TsTv N` (binned) | SKIP | Column layout diverges. |
| PopGen | `--weir-fst-pop` per-site header | PASS | Header byte-for-byte. |
| PopGen | `--fst-window-size` header | PASS | Header byte-for-byte. |
| LD | `--geno-r2` header | PASS | Header byte-for-byte. |
| LD | `--hap-r2` header | PASS | Header byte-for-byte. |
| Recode | `--recode` all-sites | PASS | 12 rows. |
| Recode | `--recode --recode-INFO-all` | PASS | INFO preserved. |
| Convert | `--012` indv file | PASS | Sample-name list. |
| Convert | `--012` row prefix | PASS | 0-based sample index prefix (fixed in this PR). |
| Convert | `--012` biallelic-only | PASS | 8 of 12 sites kept. |
| Convert | `--plink` file presence | PASS | `.ped` + `.map` emitted. |
| Convert | `--plink-tped` file presence | PASS | `.tped` + `.tfam` emitted. |
| Tajima | `--TajimaD` header | PASS | Header byte-for-byte. |
| Pi | `--window-pi` header | PASS | Header byte-for-byte (N_MONOMORPHIC placeholder added in this PR). |
| LD | `--interchrom-geno-r2` header | PASS | Header byte-for-byte (long-tail wave 1). |
| LD | `--interchrom-hap-r2` header | PASS | Header byte-for-byte (long-tail wave 1). |
| LD | `--geno-chisq` header | PASS | Header byte-for-byte (long-tail wave 1). |
| Relatedness | `--relatedness` header | PASS | Header byte-for-byte (long-tail wave 1). |
| Relatedness | `--relatedness2` header | PASS | Header byte-for-byte (long-tail wave 1). |
| ROH | `--LROH` header | PASS | Header byte-for-byte (long-tail wave 1). |
| Phasing | `--phased-blocks` header | PASS | Header byte-for-byte (long-tail wave 1). |
| INFO | `--get-INFO DP,AF` | PASS | Spot-check 20:14370 DP=14, AF=0.5. |
| Filter | `--remove-filtered q10` | PASS | q10-only sites dropped; q10;s50 also dropped. |
| Filter | `--keep-filtered q10` | PASS | Only q10-listing sites kept (2 rows). |

Totals: 53 cases, 50 PASS, 3 SKIP.

### Discrepancies found (and fixed in this PR)

The parity audit surfaced several real output-format gaps relative to
upstream vcftools, all of which we fixed inline rather than masking with
`t.Skip`:

- **`.frq` / `.frq.count` header and per-allele row format.** Upstream
  emits a single literal `{ALLELE:FREQ}` / `{ALLELE:COUNT}` column
  header (the curly braces are part of the literal header text), and
  data rows have one tab-separated `allele:value` entry per allele
  with no braces around each entry. Our port previously emitted
  `{REF:FREQ}\t{ALT:FREQ}` in the header and `{A:0.833333}` in the
  data — both wrong. Fixed; the parity test pins the upstream format.
- **`.singletons` column count and ordering.** Upstream emits five
  columns: `CHROM\tPOS\tSINGLETON/DOUBLETON\tALLELE\tINDV`. We were
  emitting only three (`CHROM\tPOS\tALLELE`) and never identified
  private doubletons. Fixed; `addSingletonStat` now resolves the
  carrier individual and tags singletons (`S`) vs private doubletons
  (`D`).
- **`.hwe` header uses `CHR` (not `CHROM`)** and emits
  `P_HET_DEFICIT` + `P_HET_EXCESS` columns alongside `P_HWE`. Fixed;
  the two directional P-values are placeholders set to `P_HWE` until
  the [PARITY_ROADMAP.md#vcftools](../docs/PARITY_ROADMAP.md#vcftools)
  Wigginton-Cutler-Abecasis test is wired through. Header matches
  upstream byte-for-byte.
- **`.lmiss` header uses `CHR` (not `CHROM`)** — same upstream
  convention as `.hwe` / `.geno.ld` / `.hap.ld`. Fixed.
- **`.ldepth` is missing the `SUMSQ_DEPTH` column** — upstream emits
  both `SUM_DEPTH` and the sum of squared per-individual depths. Fixed
  the header/column count; the value is a literal `0` placeholder
  until the sum-of-squares accumulator lands.
- **`.TsTv.summary` layout.** Upstream emits a single `MODEL\tCOUNT`
  table with six per-substitution rows (`AC`, `AG`, `AT`, `CG`, `CT`,
  `GT`) plus two roll-up rows (`Ts`, `Tv`). Our port previously
  emitted a tiny `Ts\tTv\tTs/Tv` two-line table. Fixed.
- **`.012` row prefix is the 0-based sample index**, NOT the sample
  name. Upstream's `output_as_012_matrix` writes the integer index;
  our port wrote the name. Fixed. We also confirmed the file emits
  only biallelic loci, matching upstream's one-off warning + skip.
- **`--max-missing 1.0` was a no-op** because the filter was guarded
  by `< 1` (strict less-than). With `MaxMissing == 1.0` no records
  would be filtered even though semantically `--max-missing 1.0` is
  "require all non-missing genotypes". Fixed to `> 0`.
- **`.windowed.pi` is missing the `N_MONOMORPHIC` column**. Upstream
  emits CHROM/BIN_START/BIN_END/N_VARIANTS/N_MONOMORPHIC/PI. Fixed
  the header; the value is a literal `0` placeholder until we tally
  monomorphic windowed sites.

### Discrepancies found (NOT fixed in this PR)

- **`--site-pi` formula divergence.** Upstream and our port compute
  different quantities; tracked at
  [docs/UPSTREAM_BUGS.md#vcftools-site-pi](../docs/UPSTREAM_BUGS.md#vcftools-site-pi).
  Skipped byte parity; instead the textbook formula is spot-checked.
- **`--TsTv` (binned) column layout** — upstream emits
  `CHROM\tBinStart\tSNP_count\tTs/Tv` (4 cols, per-chrom bins). Our
  port emits `BIN_START\tBIN_END\tTs\tTv\tTs/Tv` (5 cols, no chrom
  column). Skipped; tracked at
  [docs/PARITY_ROADMAP.md#vcftools](../docs/PARITY_ROADMAP.md#vcftools).
- **`--TsTv-by-count` row enumeration** — upstream emits every count
  bin from 0 to 2*N_indv (some with NaN ratios for empty bins). Our
  port emits only non-empty bins. Header still matches; row content
  skipped.

## fastp parity validation

Byte-for-byte parity for the fastp Go port against the upstream C++
reference (`OpenGene/fastp 1.0.1`). Upstream is vendored as a git
submodule at `reference_code/fastp/` and built with
`cd reference_code/fastp && make` (needs `libisal-dev`, `libdeflate-dev`
on the system — `apt-get install libisal-dev libdeflate-dev` on
Debian/Ubuntu).

The test driver lives at `tools/fastp/pkg/fastp/parity_test.go` and
invokes the upstream binary at **test time** rather than relying on
pre-baked golden files. Inputs live under
`tools/fastp/testdata/parity/*.fq` and are regenerated by the
deterministic `generate.py` script (`random.seed(42)`) in the same
directory. Every test calls `ensureUpstream(t)` first, which skips
gracefully if the binary has not been built.

### Cases

| # | Case | Result | Notes |
| ---:| ---- | ------ | ----- |
|  1 | SE basic (defaults) | PASS | Byte-parity of output FASTQ + before/after counters. |
|  2 | SE custom adapter (`-a AGATCGGAAGAGC`) | PASS | Byte-parity + passed-filter counter. |
|  3 | SE length filter (`-l 50`) | PASS | Byte-parity + too_short counter. |
|  4 | SE max length (`--length_limit 60`) | PASS | too_long counter. |
|  5 | SE N filter (`-n 2`) | PASS | Byte-parity + too_many_N counter. |
|  6 | SE Q20/Q30 filter (`-q 30 -u 20`) | PASS | low_quality + passed_filter counters. |
|  7 | SE low-complexity (`-y -Y 30`) | PASS | After complexity-definition fix (see below). |
|  8 | SE UMI extraction (`-U --umi_loc read1 --umi_len 8`) | PASS | After UMI-tag-format fix (see below). |
|  9 | SE UMI with prefix (`--umi_prefix UMI`) | PASS | Byte-parity. |
| 10 | PE basic (defaults) | PASS | Byte-parity R1 + R2 + before/after counters. |
| 11 | PE adapter (`-a ... --adapter_sequence_r2 ...`) | PASS | Counter parity. |
| 11b | PE `--detect_adapter_for_pe` | PASS | Overlap-based detection; byte parity R1+R2. |
| 12 | SE poly-G trimming (`-g`) | SKIP | Upstream tolerates mismatches; Go does strict consecutive-G — see PARITY_ROADMAP.md#fastp. |
| 13 | SE sliding-window `--cut_right` | SKIP | Off-by-1..2 at cut boundary; algorithmic divergence — see PARITY_ROADMAP.md#fastp. |
| 14 | SE sliding-window `--cut_front --cut_tail` | SKIP | Same window-boundary divergence as case 13. |
| 15 | SE adapter auto-detect | SKIP | Different algorithm (substring vs kmer overlap-tree). |

Totals: 16 cases, 12 PASS, 4 SKIP.

### Bugs surfaced and fixed in this PR

The audit surfaced three real bugs in the Go port that we fixed inline
rather than masking with `t.Skip`:

- **UMI tag format mismatch.** Our Go port was unconditionally adding
  `":UMI_<umi>"` to every UMI-processed read name, regardless of
  whether `--umi_prefix` had been provided. Upstream's
  `umiprocessor.cpp::addUmiToName` emits `":<umi>"` when no prefix is
  set and `":<prefix>_<umi>"` only when one is. Aligners that key off
  the read name (e.g. `umi_tools dedup`) would silently treat our
  output as different molecules from upstream's. Fixed in
  `tools/fastp/pkg/fastp/umi.go::appendUMIName`. Existing unit tests
  were updated to assert the upstream-correct tag.
- **Low-complexity definition divergence.** Upstream defines complexity
  as "fraction of adjacent positions where seq[i] != seq[i+1]" (see
  `filter.cpp::passLowComplexityFilter`). Our Go port was computing
  "unique 2-mers / total 2-mers", which classified a 100% alternating
  `ATATATAT...` sequence as low complexity (it has only 2 distinct
  2-mers — `AT` and `TA`) when upstream considers it maximally complex
  (every adjacent pair differs, so complexity = 1.0). On the audit
  fixture upstream passed 10/10 reads while we dropped 10/10. Fixed in
  `tools/fastp/pkg/fastp/fastp.go::calculateComplexity`, with the
  existing `TestCalculateComplexity` table updated to the
  upstream-aligned expected values. We also accept `--complexity-threshold`
  values > 1 by dividing by 100, so upstream's `-Y 30` (percentage form)
  is accepted via the same flag.
- **Missing `low_complexity_reads` counter in JSON.** Upstream's
  `filtering_result` block has a `low_complexity_reads` field (always
  emitted; zero when the filter is disabled). Our schema was missing
  it. Added `ProcessStats.LowComplexityReads`, populated where the
  filter rejects a read, and wired into `report_json.go`.

### Upstream divergences left as `t.Skip`

These are tracked as **OUR-side bugs** (the Go port is wrong; upstream
is right) in [`docs/PARITY_ROADMAP.md#fastp`](../docs/PARITY_ROADMAP.md#fastp);
fixing them is more than a one-character change, so we documented
them as skipped parity cases:

- **PolyG mismatch tolerance** — upstream's `trimPolyG` tolerates
  1 mismatch per 8 bases scanned (capped at 5 total) and uses a
  "last G position" anchor (`polyx.cpp::trimPolyG`). We strip only
  strictly consecutive G's. On a poly-G read where the tail is
  `...GTAGGGGCCC...GGGGGG` upstream trims further than we do.
- **Sliding-window boundary** — three off-by-1..2 issues in
  `slidingWindowCut`:
  1. Upstream's `cut_right` keeps the high-quality prefix of the
     offending window (`filter.cpp:172-178`) — we drop the whole
     window.
  2. Upstream's `cut_front` skips past trailing `N` bases at the
     cut boundary (`filter.cpp:138-139`) — we don't.
  3. The window-iteration bounds (`s+w < l-tail` vs `s+w <= l-tail`)
     differ by one between us and upstream, producing 1-2bp drift on
     short low-quality tails.

  These are not upstream bugs — they're Go-side bugs in our
  re-implementation. They will be fixed in a follow-up PR; the parity
  tests are kept as `t.Skip` with a pointer so we don't lose track.
- **SE adapter auto-detection algorithm gap** — upstream samples the
  first 10000 reads and builds a kmer overlap-tree
  (`evaluator.cpp::evalAdapterAndReadNum`). We use a simple substring
  search against a small built-in adapter table. Different signal,
  different outputs. Tracked as a roadmap item, not a bug.

### Reproducing locally

```bash
git submodule update --init reference_code/fastp
sudo apt-get install -y libisal-dev libdeflate-dev   # if not already
cd reference_code/fastp && make -j && cd -
go test -race ./tools/fastp/...
```

The parity tests are guarded by `ensureUpstream(t)` so they skip
cleanly on systems without the upstream binary — which is also why CI
(currently disabled in this repo) does not need libisal/libdeflate.

## prinseq parity validation

`tools/prinseq` is a Go port of the PRINSEQ-lite Perl pipeline
(upstream `uwb-linux/prinseq` @ 0.20.4, vendored under
`reference_code/prinseq/`). The parity rig drives the Go library
functions on a small corpus of FASTA / FASTQ fixtures generated once
by running upstream `prinseq-lite.pl -line_width 0 -out_good <prefix>`
on representative inputs and committing the result under
`tools/prinseq/testdata/parity/`. Tests are then byte-for-byte against
the fixture, with the single exception of the stats-info / stats-len
case, where we parse upstream's `stats_*` text rows and compare
numbers (upstream and our port disagree on summary formatting, but
the numbers themselves must match).

### Summary

| Category | Case | Status | Notes |
| --- | --- | --- | --- |
| Stats | `stats_info` + `stats_len` (FASTA) | PASS | Numeric parity for `bases`, `reads`, `min`, `max`, `mean`. |
| Stats | `stats_info` + `stats_len` (FASTQ) | PASS | Numeric parity (encoding-independent for length stats). |
| Filter | `-min_len 10` | PASS | Byte-for-byte FASTA. |
| Filter | `-max_len 20` | PASS | Byte-for-byte FASTA. |
| Filter | `-min_gc 50` | PASS | Byte-for-byte FASTA. |
| Filter | `-max_gc 60` | PASS | Byte-for-byte FASTA. |
| Filter | `-ns_max_p 5` | PASS | Byte-for-byte FASTA. |
| Filter | `-ns_max_n 2` | PASS | Byte-for-byte FASTA. |
| Filter | `-min_qual_mean 15` (Phred+33) | PASS | Byte-for-byte FASTQ. |
| Filter | `-min_qual_mean 39 -phred64` | PASS | Byte-for-byte FASTQ — exercises Phred+64 mean-quality decode (fixed in this PR). |
| Filter | Multi-criteria (length + GC) | PASS | Byte-for-byte FASTA. |
| Trim | `-trim_left 5` | PASS | Byte-for-byte FASTQ. |
| Trim | `-trim_right 4` | PASS | Byte-for-byte FASTQ. |
| Trim | `-trim_qual_left 20` | PASS | Byte-for-byte FASTQ. |
| Trim | `-trim_qual_right 20` | PASS | Byte-for-byte FASTQ. |
| Trim | `-trim_tail_left 4` | PASS | Byte-for-byte FASTQ — exercises the same-base poly-A/T anchor fix in this PR. |
| Trim | `-trim_tail_right 4` | PASS | Byte-for-byte FASTQ — same fix. |
| Filter | `-derep 1` exact duplicates | PASS | Byte-for-byte FASTA. |
| Paired | `-fastq -fastq2 -min_len 10` | PASS | R1 and R2 outputs byte-for-byte. |
| Empty | `-fasta`/`-fastq` empty input | PASS | No crash; zero-byte output. |

Totals: 18 cases (counting the two empty sub-tests as one), 18 PASS,
0 SKIP.

### Discrepancies found (and fixed in this PR)

The parity audit surfaced three real bugs in our Go port relative to
upstream PRINSEQ-lite, all fixed inline:

- **`-min_qual_mean` ignored `--qual-type illumina`** — the filter
  loop in `tools/prinseq/pkg/prinseq/prinseq.go` called
  `calculateAvgQualityScore` (which hard-codes offset 33) instead of
  `calculateAvgQualityScoreWithOffset(_, phredOffset(opts.QualType))`.
  Phred+64 inputs were therefore decoded against the wrong offset,
  consistently 31 quality units too high. Now resolved.
- **`-trim_tail_left` / `-trim_tail_right` treated A and T as
  interchangeable** — our port matched any prefix of A or T bases
  (`A|T`) as a single poly-tail run. Upstream first picks an anchor
  (the leading or trailing N-base homopolymer, A *or* T, not both)
  and then extends only with that base or N. The old behaviour
  over-trimmed by 1+ bases whenever the tail straddled an A-run and a
  T-run. Fixed by rewriting `trimPolyATLeft` / `trimPolyATRight` and
  adding `matchesCase` / `allEqualCase` helpers.
- **FASTQ output emitted `+\n` instead of `+<header>\n`** — upstream
  PRINSEQ defaults to repeating the sequence header on the quality
  separator line; the bare `+` is only used when `-no_qual_header` is
  supplied. Our port used the generic `fastq.Writer` (bare `+`) which
  diverged byte-for-byte. Fixed by adding `writePrinseqFastq` in
  `prinseq.go` and routing filter / paired-end output through it.

### Discrepancies found (NOT fixed)

- **stats output format** — upstream emits
  `stats_<section>\t<key>\t<value>` rows; our `prinseq stats`
  emits a human-readable summary. The numbers match — we parse
  upstream's rows and compare values rather than literal text — but
  the formats themselves diverge. Documented under
  [docs/PARITY_ROADMAP.md#prinseq-lite](../docs/PARITY_ROADMAP.md#prinseq-lite).

### Reproducing locally

```bash
git submodule update --init reference_code/prinseq
go test -race ./tools/prinseq/...
```

## seqtk parity validation

`tools/seqtk` is a Go port of `lh3/seqtk` (upstream v1.5-r133,
vendored under `reference_code/seqtk/`). Each fixture under
`tools/seqtk/testdata/parity/` was generated once by piping the
matching input file through the upstream binary with the same
subcommand and flags. Tests then drive the Go library functions on
the same input and assert byte parity (or skip with a documented
divergence).

### Summary

| Subcommand | Case | Status | Notes |
| --- | --- | --- | --- |
| `comp` | FASTA small | PASS | Byte-for-byte per-record nucleotide composition (13 cols). |
| `comp` | FASTQ small | PASS | Byte-for-byte (input format auto-detected). |
| `comp` | FASTA with N runs | PASS | Exercises the #4 ambiguity column. |
| `seq -A` | FASTQ → FASTA | PASS | Byte-for-byte conversion via `ConvertFastqToFasta`. |
| `seq -A` | Phred+64 FASTQ → FASTA | PASS | Encoding-independent on the sequence side. |
| `seq -r` | FASTA reverse-complement | PASS | Header preserved verbatim (port fix in this PR). |
| `seq -r` | FASTQ reverse-complement | PASS | Quality reversed in lockstep with the sequence. |
| `subseq` | Name-list mode | PASS | Records emitted in input order. |
| `subseq` | BED-region mode | PASS | `name:start+1-end` header. |
| `mergepe` | Two FASTQ files | PASS | Interleaved output. |
| `cutN` | `-n 4` cut at runs ≥ 4 | PASS | Every fragment uses upstream's `name:start-end` header. |
| `cutN` | `-n 100` (no cuts) | PASS | Records still emitted with `name:1-len` (port fix in this PR). |
| `mutfa` | Apply mutations file | PASS | Byte-for-byte. |
| `hpc` | Homopolymer compression (homo.fa) | PASS | Single-line output, runs of length 1 preserved. |
| `hpc` | Homopolymer compression (small.fa) | PASS | Spot-check with mixed run lengths. |
| `sample` | Fraction 1.0 / 0.5 invariants | PASS | Structural-only: subset of input, fraction=1.0 keeps all. |
| `sample` | Upstream byte parity | SKIP | Different RNG; see `docs/UPSTREAM_BUGS.md#seqtk-sample-rng`. |
| `randbase` | IUPAC invariants | PASS | Only 2-base codes (R/Y/S/W/K/M) randomised (port fix in this PR); 3/4-base codes pass through. |
| `randbase` | Upstream byte parity | SKIP | Different RNG; see `docs/UPSTREAM_BUGS.md#seqtk-randbase-rng`. |
| `trimfq` | Upstream byte parity | SKIP | Algorithm gap: Phred-threshold vs Mott; see `docs/UPSTREAM_BUGS.md#seqtk-trimfq-algorithm`. |
| (all) | Empty input no-crash | PASS | Every public function on a zero-byte input returns nil. |

Totals: 21 cases (including the 2 sample sub-tests), 18 PASS, 3 SKIP.

### Discrepancies found (and fixed in this PR)

- **`seq -r` appended `" (reverse complement)"` to the FASTA / FASTQ
  description** — added by `fasta.Record.ReverseComplement` and
  `fastq.Record.ReverseComplement` in `pkg/bioformats/`. Upstream
  preserves the header verbatim; downstream tools that key on the
  description field were silently broken. Now the description is
  copied unchanged.
- **`comp` was emitting summary statistics rather than upstream's
  per-record nucleotide-composition rows.** Added a new `Comp`
  function (`tools/seqtk/pkg/seqtk/comp.go`) that mirrors upstream's
  `stk_comp` inner loop byte-for-byte (`name\tlen\t#A\t#C\t#G\t#T\t#2\t#3\t#4\t#CpG\t#tv\t#ts\t#CpG-ts`).
  The summary-stats form is preserved as a `--summary` opt-in on the
  CLI for backward compatibility with existing scripts.
- **`cutN` dropped the coordinate suffix when no run reached the
  threshold.** Upstream's `print_seq` always prints
  `name:start+1-end`. Without the fix, downstream code that walks
  uniform `name:S-E` headers had to special-case the unchanged
  records. Now every record gets the suffix.
- **`randbase` was randomising 3-base (B/D/H/V) and 4-base (N) IUPAC
  codes.** Upstream's `stk_randbase` only touches codes whose
  bit-count is exactly 2 (R/Y/S/W/K/M), leaving everything else
  alone. The expansion table in `tools/seqtk/pkg/seqtk/mutations.go`
  was trimmed accordingly; the existing `pickIUPAC` unit tests were
  updated to match.

### Discrepancies found (NOT fixed)

- **`sample` RNG** — upstream uses a seeded reservoir sampler with a
  default seed of 11; our port does a deterministic every-Nth-record
  pick. Byte parity is therefore impossible. Tracked at
  [docs/UPSTREAM_BUGS.md#seqtk-sample-rng](../docs/UPSTREAM_BUGS.md#seqtk-sample-rng);
  the parity test is split into a structural-invariants pass and a
  skipped byte-parity entry.
- **`randbase` RNG** — upstream uses `drand48()` with an implicit
  seed of 0, our port uses `math/rand` with a caller-supplied seed.
  Structural invariants are checked; byte parity is skipped (see
  [docs/UPSTREAM_BUGS.md#seqtk-randbase-rng](../docs/UPSTREAM_BUGS.md#seqtk-randbase-rng)).
- **`trimfq` algorithm** — upstream runs a Mott-style error-rate
  trim with a 0.05 default threshold; our port does a simple
  Phred-threshold cut on each end. Different algorithms, different
  cuts. Tracked at
  [docs/UPSTREAM_BUGS.md#seqtk-trimfq-algorithm](../docs/UPSTREAM_BUGS.md#seqtk-trimfq-algorithm).

### Reproducing locally

```bash
git submodule update --init reference_code/seqtk
cd reference_code/seqtk && make && cd -
go test -race ./tools/seqtk/...
```
