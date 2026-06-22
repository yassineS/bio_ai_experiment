# Manuscript results

Concrete experimental artifacts for the manuscript claims (see
[`../01_CLAIMS_AND_EXPERIMENTS.md`](../01_CLAIMS_AND_EXPERIMENTS.md)). Everything
here was produced and validated in the cloud sandbox; the resource-heavy runs are
handed off to a local box via [`HEAVY_TIER_RUNBOOK.md`](HEAVY_TIER_RUNBOOK.md).

## Done in-sandbox (committed here)

| Artifact | Claim | Headline |
|---|---|---|
| [`conformance_run.txt`](conformance_run.txt) | C2 | htslib `test/` + htscodecs corpora + edge-case battery through our binaries: **89 PASS / 0 SKIP / 0 FAIL**. |
| [`parity_statistics.md`](parity_statistics.md) | C2/C3 | Wilson + Clopper-Pearson CIs on every parity rate (conformance all-pass → CP lower bound **≥ ~96%**, not a bare "100%"); max-FP-deviation mechanism. Backed by the stdlib-only `pipeline/stats` package. |
| [`differential_fuzzing.md`](differential_fuzzing.md) | C2 | stdout+stderr+exit diff over fuzzed inputs; **0 silent-corruption-on-valid-input** divergences — all are error-handling-convention cosmetics or illegal-byte passthrough on malformed input. (`differential_fuzzing.json`, `differential_fuzzing_raw.txt`.) |
| [`ablation.md`](ablation.md) | C7 | Validation-layer ablation — defects uniquely caught per layer + escape-rate-if-removed; completes the labeled bug corpus ([`../bug_corpus.md`](../bug_corpus.md)). |
| [`flag_compat.md`](flag_compat.md) | C4/C1 | Flag-compatibility **83.6%** weighted (2197/2628 upstream flag slots), per-tool with denominators. Backed by `pipeline/cmd/flagcompat`. |
| [`transpiler_baseline.md`](transpiler_baseline.md) | C5 ★★ | The counterfactual: a rule-based C→Go transpiler (c2go) cannot reach the bar — panics on trivial C, emits non-idiomatic `unsafe` + libc-shim code, and cannot ingest C++/Perl at all. (`../../../scripts/transpiler/`.) |

### Supporting code added for these results

- `pipeline/stats/` — Wilson + Clopper-Pearson binomial CIs (stdlib-only).
- `pipeline/cmd/flagcompat/` — flag-surface compatibility enumerator.
- `pipeline/cmd/realparity/` — real-data, multi-contig differential parity +
  performance runner (used by the local-box hand-off; reuses
  `runner.StripProvenance`/`CompareByteExact`).
- `pipeline/roundtrip/interop.go` — bidirectional container interop
  (ours-writes/upstream-reads **and** upstream-writes/ours-reads) for
  BGZF/BAM/CRAM/VCF.gz/BCF/FASTQ + `.bai`/`.csi`/`.tbi`, on multi-contig
  fixtures; wired into `full-validation`.
- `scripts/transpiler/run_c2go_baseline.sh` — runnable counterfactual harness.

## Handed off to the local box (`HEAVY_TIER_RUNBOOK.md`)

- **Real-data parity + interop + performance** on the GIAB files as
  whole-genome, **multi-contig** inputs (`realparity` + `full-validation`). No
  truth set — upstream is the oracle. (The GIAB *biological-concordance*
  experiment is intentionally **not** run.)
- **Large-tier** parity + round-trip + performance, with the perf-stats upgrade
  (median + IQR + ratio CI) and a pinned hardware spec.
- **K-run process-reproducibility** (needs LLM-API budget).

## Needs people, not a box

- Timed human-port anchor (C5 calibration) and a second independent bug-corpus
  labeler (C7 κ). See the runbook's "Not in scope" section.
