# Parity statistics: confidence intervals and max floating-point deviation

This document upgrades every reported parity rate to carry a **stated
denominator** and a **95% confidence interval** (both Wilson score and exact
Clopper-Pearson), and records the **maximum observed floating-point deviation**
for the tolerance-compared cells rather than a bare "within epsilon" pass/fail.
It backs manuscript claims **C2** (correctness / behavioral equivalence) and
**C3** rigor, and addresses the Tier-C item "CIs everywhere; max-FP-deviation
reporting" in
[`01_CLAIMS_AND_EXPERIMENTS.md`](../01_CLAIMS_AND_EXPERIMENTS.md).

All intervals are computed by the pure-Go, standard-library-only package
[`pipeline/stats`](../../../pipeline/stats/binomialci.go)
(`WilsonCI`, `ClopperPearsonCI`), validated against textbook / R `binom.test`
values in `pipeline/stats/binomialci_test.go`.

## Method

A parity rate is a binomial proportion: of `n` compared cells, `k` matched the
upstream oracle (byte-for-byte after provenance stripping, or within the
floating-point tolerance for similarity cells). The point estimate `k/n` is
reported with two 95% intervals:

- **Wilson score interval** — the recommended default; well-calibrated for
  small `n` and for proportions near 0 or 1, where the normal-approximation
  (Wald) interval is badly miscalibrated.
- **Clopper-Pearson (exact) interval** — conservative (coverage at least 95%);
  the honest choice for the all-pass cells (`k == n`), where Wilson can look
  optimistically narrow.

`z = 1.95996` (two-sided 95%) for Wilson; `alpha = 0.05` for Clopper-Pearson.

## Parity rates with confidence intervals

| Experiment | Successes / total | Point est. (%) | 95% Wilson CI (%) | 95% Clopper-Pearson CI (%) | Source |
|---|---|---:|---|---|---|
| Upstream-audit total | 223 / 314 | 71.02 | [65.77, 75.76] | [65.66, 75.98] | [`README.md`][r] audit table |
| — bedtools (10 subcmds) | 85 / 127 | 66.93 | [58.36, 74.51] | [58.03, 75.02] | [`README.md`][r], PR #55 |
| — sickle + skewer | 22 / 27 | 81.48 | [63.30, 91.82] | [61.92, 93.70] | [`README.md`][r], PR #73 |
| — samtools (6 subcmds) | 34 / 43 | 79.07 | [64.79, 88.58] | [63.96, 89.96] | [`README.md`][r], PR #75 |
| — bcftools (6 subcmds) | 32 / 52 | 61.54 | [47.96, 73.53] | [47.02, 74.70] | [`README.md`][r], PR #74 |
| — mosdepth + vcftools | 50 / 65 | 76.92 | [65.36, 85.49] | [64.81, 86.47] | [`README.md`][r], PR #76 |
| Spec conformance (htslib / htscodecs) | 89 / 89 | 100.00 | [95.86, 100.00] | [95.94, 100.00] | [`conformance_run.txt`](conformance_run.txt) |
| realbench (real GIAB HG002/GRCh38 chr20, all-tool sweep) | 129 / 131 | 98.47 | [94.60, 99.58] | [94.59, 99.81] | chr20 realbench (see notes) |
| Parity-pipeline (small scale) | 398 / 400 | 99.50 | [98.20, 99.86] | [98.21, 99.94] | [`01_CLAIMS_AND_EXPERIMENTS.md`][c1] C2 row |
| Parity-pipeline (medium scale) | DIVERGE = 0 (denominator not formally recorded) | — | — | — | [`01_CLAIMS_AND_EXPERIMENTS.md`][c1] C2 row |

[r]: ../../../README.md
[c1]: ../01_CLAIMS_AND_EXPERIMENTS.md

### Notes on each count

- **Upstream-audit total (223 / 314).** The denominator is the count of
  audit cases that were actually executed and compared; the `89` documented
  `t.Skip()`s are unimplemented features, not failures, and are **excluded
  from the denominator** (a skip is "not tested", not "tested and failed").
  Counting them would conflate coverage with correctness. Source: the
  "Validated parity against upstream test suites" table in
  [`README.md`](../../../README.md) (rows summing to 314 / 223 / 89) and
  `tools/PARITY_VALIDATION.md` for the
  per-subcommand breakdown.
- **Conformance (89 / 89).** Independent third-party corpora — the htscodecs
  rANS / arithmetic codec vectors and htslib SAM/BAM/CRAM `test/` fixtures —
  run through our binaries with **zero failures**. The denominator is the
  count of `PASS` lines in [`conformance_run.txt`](conformance_run.txt)
  (75 leaf subtests + 14 top-level test functions = 89, 0 `FAIL`). Because
  this is an all-pass cell, the Clopper-Pearson lower bound (95.94%) is the
  honest figure to quote: with 89 independent passes we can state with 95%
  confidence the true conformance rate is **at least ~96%**, not "100%".
- **realbench chr20 (129 / 131).** The real-data `realbench` sweep on the GIAB
  HG002/GRCh38 chr20 tier scored **PASS = 129, DIFF = 2, ERROR = 0, SKIP = 1**.
  The denominator is the **131 compared cells** (PASS + DIFF); the 1 SKIP
  (`bcftools csq`, an ours-only feature with no upstream oracle) is *not tested*,
  not *tested-and-failed*, so it is excluded — the same skip-accounting rule used
  for the upstream-audit total. The 2 DIFFs are the **accepted `samtools
  consensus` libm last-ULP `cq`-column residuals** (base/seq/qual bytes
  byte-exact; see [`max_fp_deviation.md`](max_fp_deviation.md)), not silent
  corruption. The **exome + wgs 60× whole-genome tiers are in progress** (Seqera
  run `mS3IH42QfGTWO`) and will be folded in by a follow-up.
- **Parity-pipeline small (398 / 400).** The `parity-pipeline` byte-compare
  (provenance-stripped) over the small-scale fixture matrix. Source: the C2
  status row in
  [`01_CLAIMS_AND_EXPERIMENTS.md`](../01_CLAIMS_AND_EXPERIMENTS.md)
  ("have 398/400 small, 0 DIVERGE medium"). The two non-matching cells are
  not silent divergences; see that claim row and
  [`bug_corpus.md`](../bug_corpus.md) for the tracked items.
- **Parity-pipeline medium.** Recorded only as "0 DIVERGE" — the explicit
  cell denominator at medium scale is **not formally recorded** in the source
  docs. A `0 / n` result with an unstated `n` carries no usable interval (with
  `n` unknown, the proportion is undefined); this is exactly the
  coverage-denominator gap called out below. To turn it into a real datum,
  the medium run must emit `cells_compared` alongside `DIVERGE`.

## The denominator / coverage caveat (cross-reference C1 / C2)

**A parity percentage means little without the cell-coverage denominator.**
Two distinct denominators are at play and must not be confused:

1. **Comparison denominator** (used above): of the cells we *did* compare, how
   many matched. This is what the CIs quantify.
2. **Coverage denominator** (`cells_exercised / cells_enumerable`, claim **C1**
   in [`01_CLAIMS_AND_EXPERIMENTS.md`](../01_CLAIMS_AND_EXPERIMENTS.md)): of the
   enumerable flag-combination space, what fraction did we exercise at all.
   Per C1 this is **NEED** — subcommand coverage is complete, but
   flag-combination coverage is **not yet formally measured**.

A high comparison-parity rate over a thin slice of the flag space is weaker
evidence than a moderate rate over a near-exhaustive slice. The CIs here bound
sampling error *given the cells we ran*; they say nothing about cells we never
enumerated. Until the C1 flag-space enumerator lands and reports a coverage
denominator per tool, every rate in this document should be read as
**"conditional on the exercised cell set"**. The skip accounting (89 documented
`t.Skip()`s, each cross-referenced to
[`docs/PARITY_ROADMAP.md`](../../../docs/PARITY_ROADMAP.md)) is the current best proxy
for the unexercised region, but it is a lower bound on what is untested, not the
full enumerable space.

## Maximum floating-point deviation

Some cells cannot be byte-exact because they emit floating-point fields whose
last digits depend on libm / formatting details (e.g. C++ `ostream` `%g` vs Go,
or glibc `pow`/`log` vs pure-Go `math`). These are compared with a **relative
tolerance** rather than byte equality.

### What the comparator records

The comparator is
[`pipeline/runner/compare.go`](../../../pipeline/runner/compare.go),
function `CompareSimilarity(ours, upstream []byte, eps float64)`. It splits both
streams (after `StripProvenance`) into lines and whitespace-delimited fields,
then for every field that parses as a number on both sides computes the
**relative deviation**

```text
relDev(a, b) = |a - b| / max(|a|, |b|)     (0 when a == b; 0 when both are 0)
```

It tracks the running maximum of this quantity across all numeric fields in the
cell and returns it in `CompareResult.MaxDeviation`. The cell passes when the
max stays at or below `eps`; the **first** field that exceeds `eps` fails the
cell and `MaxDeviation` / `Detail` pinpoint it
(`"line L field F numeric deviation D (a vs b)"`). Non-numeric tokens must match
exactly. `CompareOutputFiles` (used by the multi-file vcftools / mosdepth
matrices) propagates the per-file maximum the same way.

Key facts to report from this code:

- The recorded metric is the **maximum relative deviation** (`MaxDeviation`,
  a `float64`), not just a boolean pass.
- The **default tolerance is `eps = similarityEpsilon = 1e-6`** (relative); a
  matrix entry may widen it per-cell via `Entry.Tolerance` (see
  `resolveEpsilon`). So a passing similarity cell guarantees every numeric
  field agreed to **better than 1 part in 10^6** unless that cell declared a
  looser tolerance.
- The comparator records **relative** deviation. To also report an **absolute**
  max deviation (requested for the FP cells), read the failing field's `(a vs
  b)` pair from `Detail`, or extend `CompareResult` with an `AbsDeviation`
  accumulator (`|a-b|`) alongside the existing relative one — a one-field
  change in the same loop in `CompareSimilarity`.

### Observed maximum deviation

> **Now materialised — see [`max_fp_deviation.md`](max_fp_deviation.md).** The
> per-tool max abs/rel deviation is measured and tabulated: **byte-exact 0.0 for
> the vast majority of tools** (all of `vcftools`, `mosdepth`, most `samtools`/
> `bcftools`/`bed*`, `seqtk`/`fastp`/`prinseq`/`sickle`/`skewer`,
> `bgzip`/`tabix`/`htsfile`), with exactly **three documented last-ULP
> residuals**: `bcftools call -m` QUAL (max rel ~7.2e-6), `bedgenomecov`
> histogram fraction (max rel 9.7e-6, abs 1.0e-6), and `samtools consensus` gap5
> `cq` (discrete ±1 phred, base/seq/qual bytes byte-exact). `MaxAbsDeviation` was
> added to `CompareResult`/`Result` (JSON `max_abs_deviation`) alongside the
> relative one so both are recorded per cell. The paragraph below is retained for
> the bounding argument.

The largest *passing* relative deviation is, by construction, bounded by each
cell's tolerance: for every similarity cell that passes at the default setting,
`MaxDeviation <= 1e-6`; for a cell with a widened `Entry.Tolerance`, it is
bounded by that entry's value.

To obtain the concrete per-cell and aggregate maxima cheaply, run the
similarity cells and capture `MaxDeviation`:

```bash
# Small-scale parity over the float-bearing tools (vcftools/mosdepth carry the
# similarity cells); MaxDeviation is computed per cell in CompareSimilarity.
go run ./pipeline/cmd/parity-pipeline -scale=small -tools=vcftools,mosdepth
```

then read `CompareResult.MaxDeviation` (the report writer in
`pipeline/runner/report.go` can be extended to emit the per-cell and global
maxima; today the field is available to callers but not printed in the summary).
A full run requires the upstream binaries to be built (the suites are not
hermetic — they shell out to the real `vcftools` / `mosdepth`), so it is gated
behind the same upstream-build step as CI's live-parity job; for that reason the
headline number is reported here as the **mechanism and the guaranteed bound**
(`<= 1e-6` relative for default-tolerance passing cells) rather than a single
measured value, per the brief's "do not spend more than a few minutes on a slow
run" instruction.

## Reproducing the intervals

```bash
go test ./pipeline/stats/...        # validates the CI math against textbook values
```

The intervals in the table above were produced by calling
`stats.WilsonCI(k, n, 1.95996)` and `stats.ClopperPearsonCI(k, n, 0.05)` for
each `(k, n)` pair listed in the Source column.
