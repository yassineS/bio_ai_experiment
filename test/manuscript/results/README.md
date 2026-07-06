# Manuscript results

Concrete experimental artifacts for the manuscript claims (see
[`../01_CLAIMS_AND_EXPERIMENTS.md`](../01_CLAIMS_AND_EXPERIMENTS.md)). Everything
here was produced and validated in the cloud sandbox; the resource-heavy runs are
handed off to a local box via [`HEAVY_TIER_RUNBOOK.md`](HEAVY_TIER_RUNBOOK.md).

> **Note:** the performance harness has moved to `realbench`
> (`pipeline/realbench`, real GIAB data via the `nextflow/` Seqera pipeline); the
> synthetic-derived performance figures captured here are point-in-time records
> pending regeneration from a real-data run.

## Done in-sandbox (committed here)

| Artifact | Claim | Headline |
|---|---|---|
| [`conformance_run.txt`](conformance_run.txt) | C2 | htslib `test/` + htscodecs corpora + edge-case battery through our binaries: **89 PASS / 0 SKIP / 0 FAIL**. |
| [`parity_statistics.md`](parity_statistics.md) | C2/C3 | Wilson + Clopper-Pearson CIs on every parity rate (conformance all-pass → CP lower bound **≥ ~96%**, not a bare "100%"); max-FP-deviation mechanism. Backed by the stdlib-only `pipeline/stats` package. |
| [`differential_fuzzing.md`](differential_fuzzing.md) | C2 | stdout+stderr+exit diff over fuzzed inputs; **0 silent-corruption-on-valid-input** divergences — all are error-handling-convention cosmetics or illegal-byte passthrough on malformed input. (`differential_fuzzing.json`, `differential_fuzzing_raw.txt`.) |
| [`ablation.md`](ablation.md) | C7 | Validation-layer ablation — defects uniquely caught per layer + escape-rate-if-removed; completes the labeled bug corpus ([`../bug_corpus.md`](../bug_corpus.md)). |
| [`flag_compat.md`](flag_compat.md) | C4/C1 | Flag-compatibility **83.6%** weighted (2197/2628 upstream flag slots), per-tool with denominators. Backed by `pipeline/cmd/flagcompat`. |
| [`transpiler_baseline.md`](transpiler_baseline.md) | C5 ★★ | The counterfactual: a rule-based C→Go transpiler (c2go) cannot reach the bar — panics on trivial C, emits non-idiomatic `unsafe` + libc-shim code, and cannot ingest C++/Perl at all. (`../../../scripts/transpiler/`.) |
| [`branch_coverage.md`](branch_coverage.md) | C2 (G1) | Go **statement coverage of the parity-exercised code = 64.25 %**, per-tool breakdown — the sharpest answer to "which input regions are untested". |
| [`max_fp_deviation.md`](max_fp_deviation.md) | C2 (G2) | Per-tool **max abs/rel FP deviation**: byte-exact **0.0** for the vast majority; three documented last-ULP residuals (`bcftools call -m` QUAL, `bedgenomecov` fraction, `samtools consensus` `cq`). |
| [`view_speed.md`](view_speed.md) | C3 (G4) | `samtools view` region/BED→SAM: the historical ~12× is already ~2.5×; this pass cut hot-loop allocations 23–37× and peak RSS ~15 %; residual is the documented pure-Go inflate floor. |
| [`nfcore_dropin/`](nfcore_dropin/README.md) | C4 ★ (G3) | Our `samtools`/`flagstat` swapped **unchanged** into a real nf-core-style Nextflow module, run end-to-end — concrete drop-in usability evidence. |
| chr20 realbench | C2/C3 | Real GIAB HG002/GRCh38 chr20 all-tool sweep: **PASS = 129, DIFF = 2, ERROR = 0, SKIP = 1** (2 DIFF = accepted `consensus` libm `cq` residual; 1 SKIP = ours-only `bcftools csq`). Drove ~10 real-data bug fixes (bug corpus A26–A35). Exome + wgs 60× tiers **in progress** (Seqera run `mS3IH42QfGTWO`). |

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
  truth set — upstream is the oracle. GIAB *biological-concordance* (hap.py /
  vcfeval variant-calling accuracy) is **out of scope** for this project: it
  tests byte-exact parity against the upstream oracles, not variant-calling
  concordance. That experiment is **not** part of the manuscript — not deferred,
  not planned.
- **Large-tier** parity + round-trip + performance, with the perf-stats upgrade
  (median + IQR + ratio CI) and a pinned hardware spec.
- **K-run process-reproducibility** (needs LLM-API budget).

## Needs people, not a box

- Timed human-port anchor (C5 calibration) and a second independent bug-corpus
  labeler (C7 κ). See the runbook's "Not in scope" section.
