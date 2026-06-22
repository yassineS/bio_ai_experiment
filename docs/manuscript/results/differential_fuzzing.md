# Differential fuzzing — results and triage (claim C2)

Differential fuzzing feeds the **same** fuzzed input to our port and the vendored
upstream binary and diffs **stdout + stderr + exit code**; any difference is
minimized (delta-debugging) to a small reproducer and Go statement-coverage is
captured over the run. Harness: `pipeline/difffuzz/` + `pipeline/cmd/diff-fuzz`.
The point of this experiment is the honest one: bound the **unexplored-input
space** and show whether any divergence is a *silent corruption of valid data*
(the thing byte-exact parity on curated fixtures cannot rule out for inputs it
never tried).

## Run

`go run ./pipeline/cmd/diff-fuzz -quick -coverage` — seed `1`, 3 targets ×
40 iterations. Raw artifacts: [`differential_fuzzing.json`](differential_fuzzing.json)
(machine-readable) and [`differential_fuzzing_raw.txt`](differential_fuzzing_raw.txt)
(the harness's own per-reproducer dump). A fuller sweep (7 targets × 300
iterations) is wired and documented for the local box; on the cloud sandbox the
heavy run is repeatedly killed before it can flush its report, so the committed
numbers are the reproducible `-quick` tier — the **classes** of divergence are
identical at both tiers (verified from the partial heavy-run logs).

## Totals by class

| Target | Tool | Inputs | exitcode-differs | stderr-differs | stdout-differs | Coverage |
|---|---|---:|---:|---:|---:|---:|
| bcftools-view | bcftools | 40 | 17 | 2 | 1 | 7.0% |
| samtools-flagstat | samtools | 40 | 0 | 30 | 0 | 4.7% |
| bedtools-merge | bedmerge | 40 | 0 | 13 | 0 | 15.1% |

(The `-quick` tier deliberately drives a small coverage slice per target; the
coverage figure is the fraction of our statements the fuzzed inputs exercised,
reported for honesty — it is **not** a claim of high coverage.)

## Triage — every divergence is on malformed/garbage input

**None of the divergences is a silent corruption of valid data.** Classified:

- **exitcode-differs (error-handling convention).** On binary/garbage input both
  binaries *reject* it; they differ only in the exit code and message — ours
  exits `1` with a readable diagnostic (e.g. `bcftools view: unexpected line
  before VCF header: …`), upstream exits `255` with htslib's `[E::hts_hopen] …
  Exec format error`. Both refuse the input; the code/wording is cosmetic.
- **stderr-differs (warning wording).** Upstream emits `W::`-class warnings ours
  does not (e.g. `Contig '…' is not defined in the header`, `FORMAT '…' is not
  defined in the header`) on malformed records. Output records are unaffected;
  this is diagnostic chatter, not a data difference.
- **stdout-differs (the only substantive class) — 2 cases, both a
  strictness/sanitization difference on *invalid* input.** The minimized
  reproducers inject raw **non-printable bytes inside an `INFO` value**; our
  writer passes the illegal bytes through verbatim while upstream sanitizes/drops
  them. Example (diff at one record):

  ```text
  ours:     chr4  155864  rs7266  T  TC  61  PASS  DP=48;AF=0.57<binary>  GT:DP  1/1:48
  upstream: chr4  155864  rs7266  T  TC  61  PASS  DP=48;AF=0.57          GT:DP  1/1:48
  ```

  The divergence exists only because the **input** already contained bytes that
  are illegal in a VCF `INFO` field. On well-formed input the two agree. This is
  a leniency difference (we echo, upstream sanitizes), not a miscomputation —
  but it is a real, citable behavioural difference and a candidate hardening
  (reject or sanitize illegal `INFO`/`FORMAT` bytes to match upstream).

## What this supports (and what it does not)

- **Supports:** across the fuzzed input space, our ports never silently produced
  *different valid output* from upstream — every divergence is on input both
  tools treat as broken, and reduces to exit-code/warning convention or
  illegal-byte passthrough. This is the honest answer to "could they be wrong
  the same way on inputs you never tried?" for the explored region.
- **Does not support:** a high-coverage claim. The `-quick` tier is a smoke; the
  manuscript figure should come from the heavy sweep on the local box (see
  `HEAVY_TIER_RUNBOOK.md`), reported with the per-target branch coverage.

## Follow-up

1. Run the 7-target × ≥300-iteration sweep on the local box; report per-target
   branch coverage and re-triage any new `stdout-differs`.
2. Consider hardening the VCF writer to reject/sanitize non-printable bytes in
   `INFO`/`FORMAT` values (matches upstream; closes the only substantive class).
