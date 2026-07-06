# 01 — Claims, experiments, metrics: have vs. need

The thesis decomposed into atomic, falsifiable claims. For each: the precise metric, the
experiment, the baseline/control, the statistic, and a **status** — `HAVE` (evidence exists),
`PARTIAL` (exists but must be strengthened), or `NEED` (new experiment required). The `NEED`
and `PARTIAL` rows **are the to-build list** for the manuscript.

Legend: ★ = load-bearing for acceptance (fix before submission).

---

## C1 — Feasibility (descriptive)

*Claim:* agents drove a bounded, named, interdependent legacy tool set to working,
parity-passing Go.

| Metric | Experiment | Baseline | Status |
|---|---|---|---|
| Per-tool parity-complete (Y/N) + subcommand coverage % | enumerate upstream subcommands/flags; report covered/total | upstream tool surface | **HAVE** — 53 CLIs / 13 families; per-tool tallies in `PROJECT_STATUS.md` |
| Flag-combo **cell** coverage (cells exercised / enumerable) | count matrix cells vs. enumerated flag space per tool | — | **NEED** — subcommand coverage is complete; **flag-combination coverage is not formally measured**. Build a coverage denominator. |

**To build:** a flag-space enumerator that reports, per tool, `cells_exercised / cells_enumerable`
so every parity rate has a stated denominator coverage.

---

## C2 — Correctness / behavioral equivalence ★ (the project's strength — formalize it)

*Claim:* output is byte-exact (or bounded-tolerance) vs the original across a defined input space.

| Metric | Experiment | Baseline (oracle) | Status |
|---|---|---|---|
| % parity cells byte-identical, **with Wilson/Clopper-Pearson CI** | `parity-pipeline` byte-compare (provenance-stripped) | upstream binary | **PARTIAL** — have 398/400 small, 0 DIVERGE medium; **add CIs and the coverage denominator** |
| FP-tolerance cells: **max abs/rel deviation per tool** | similarity comparator already emits deltas | upstream | **PARTIAL** — have ε-pass/fail; **report the max observed deviation, not just "within ε"** |
| Round-trip / metamorphic pass (oracle-free) | `pipeline/roundtrip` BGZF/BAM/CRAM/VCF↔BCF/FASTQ | self | **HAVE** — 0 FAIL medium; frame as metamorphic |
| Coverage over the input space (Go branch coverage during the parity run) | `-coverprofile` over the parity sweep | — | **NEED** — report what fraction of branches the oracle actually exercised |
| **Differential fuzzing**: feed both binaries the same fuzzed input, diff stdout+stderr+exit | extend the 10 existing `Fuzz` fns into a differential harness | upstream binary | **NEED ★** — strongest defense vs. "untested input regions" and the only honest answer to "both wrong the same way" for unexplored inputs |
| Metamorphic relations beyond round-trip (permutation invariance, compositionality, monotonicity, sort idempotence) | property tests over real + fuzzed inputs | self (no oracle) | **NEED** — addresses the "both implementations share a misconception" limit of differential testing |
| Conformance against the **spec's own corpus** | run htslib `test/` fixtures + htscodecs-corpus codec vectors through our binaries | htslib expected outputs | **NEED** — far more credible than only our seeded fixtures; ready-made oracle |

**To build (C2, priority order):**

1. ★ Differential fuzzing harness (both binaries, diff output+exit) + branch-coverage report.
2. Adopt htslib `test/` + htscodecs-corpus as conformance suites; report pass rates.
3. Explicit metamorphic relations (permutation/compositionality/monotonicity).
4. Upgrade every parity rate to carry a CI and a coverage denominator; report max FP deviation.

**The fix-on-port ledger as a first-class result:** maintain `docs/UPSTREAM_BUGS.md` as a
labeled table {tool, bug class, discovery method, disposition, upstream PR/issue, whether
upstream accepted the fix}. "LLM-agent differential porting surfaced N latent defects in
widely-used tools" is itself publishable and externally validates the method. **Status: PARTIAL**
(the ledger exists; complete it and pursue upstream confirmation).

---

## C3 — Performance (descriptive)

*Claim:* equal-or-better wall-clock, CPU, peak RSS, as a memory-safe cgo-free binary.

| Metric | Experiment | Baseline | Status |
|---|---|---|---|
| speedup = t_up/t_go; CPU×; RSS× per cell, **median + IQR over ≥5 reps** | `pipeline/bench` (wait4 rusage) at multiple scales | upstream, same machine | **PARTIAL** — have min-of-N at medium; **switch to median+IQR/Hodges-Lehmann + CI; report cold/warm; report the slower cells (mpileup, call, isec) plainly** |
| Hardware/scale anchoring | record CPU model/cores/RAM/kernel; run smoke→large | — | **NEED ★** — **no hardware spec is recorded**; the large tier was never run (disk-bound). Run large on a fat node and pin the spec |
| Memory-safety / deployment usability deltas | cgo-free static binary; single-binary deploy; CVE-class elimination argument | upstream C | **HAVE (qualitative)** — operationalize as a usability metric in C4 |

**To build:** statistical upgrade (median+IQR+CI, multi-scale, hardware spec pinned); a true
**large-tier** run on a ≥30 GB-scratch node; explicit per-cell honesty including regressions.

---

## C4 — Usability / drop-in (operationalize or demote)

*Claim:* the Go CLIs are drop-in and improve usability. "Usability" is the softest word in the
thesis — make it measurable.

| Metric | Experiment | Baseline | Status |
|---|---|---|---|
| **Flag-compat %**: upstream flags accepted with identical semantics | flag matrix vs upstream | upstream CLI surface | **PARTIAL** — measure and report as a number |
| **Pipeline drop-in**: substitute our binary in a real Nextflow/nf-core (or Snakemake/Galaxy) module and run an end-to-end workflow unchanged | run nf-core/rnaseq-style or a samtools/bcftools step with only a binary repoint | upstream tool in same workflow | **NEED ★** — concrete, credible usability evidence: every orchestrator addresses tools by literal CLI string + paths, so a clean swap is a real result |
| Memory-safety class elimination | enumerate CWE classes structurally impossible in safe Go (buffer overflow, UAF, etc.) vs C | upstream | **HAVE (qualitative)** — strengthen with the OSS-Fuzz htslib CVE history as motivation |
| Single static binary / no dependency-hell | build + run with zero shared-lib deps | upstream (htslib link chain) | **HAVE** |

**To build:** a flag-compat %; at least one **real-workflow drop-in demonstration**
(Nextflow/nf-core module repoint, end-to-end run, identical results). Demote remaining
"usability" prose to qualitative unless metricized.

---

## C5 — Effort / cost / feasibility counterfactual ★★ (the desk-reject risk)

*Claim:* agents make this feasible / cheaper than the human alternative.

| Metric | Experiment | Baseline / control | Status |
|---|---|---|---|
| Agent cost: wall-clock, tokens/$, interventions, retries per tool | instrument/replay the agent transcripts (raw logs exist) | — | **NEED ★** — quantify human-supervision cost; never report agent-cost alone |
| Non-LLM automated baseline reaches the bar? | run c2go/C2Rust-style transpiler on same inputs; measure parity %, compiles?, memory-safe?, idiomatic? | rule-based transpiler | **NEED ★★** — primary counterfactual; expected to fail the bar, isolating the LLM contribution |
| Human-effort anchor | 1–2 expert devs port a scoped slice (e.g. `sickle` + N bcftools subcommands) to the same parity bar, time-tracked | timed human | **NEED ★★** — calibration anchor; report as illustrative, CIs-not-possible |
| Human-effort bound | COCOMO-II on measured Go LOC; literature port-timelines | model + literature | **NEED** — bounded range with sensitivity table on day-rate and LOC-coefficient; triangulate vs the human anchor |
| Autonomy level + correction state | label the workflow on Human-Agency-Scale (H1–H5) / Levels-of-Autonomy; report correction rate | — | **NEED** — report autonomy + correction, not pass-rate-in-a-vacuum |

**To build (C5 is the make-or-break):** the non-LLM transpiler control + the timed human
anchor(s) + COCOMO/literature priors + supervision-cost accounting + autonomy-level labeling.
Report effort-saved only as a **bounded interval with explicit assumptions and a sensitivity
table** — never a headline multiplier. Subtract supervision cost from the denominator.

---

## C6 — Generality (descriptive, bounded)

*Claim:* this generalizes beyond one lucky tool.

| Metric | Experiment | Baseline | Status |
|---|---|---|---|
| Variance of C1–C5 across 13 families; success vs. tool size/language/format-count | treat each tool as a data point; correlate | internal | **PARTIAL** — have per-tool data; **add the cross-tool analysis and characterize where it stalled** |
| Scope of the claim | state the boundary explicitly | — | **HAVE (to write)** — scope to "tools with an executable reference + a checkable spec"; n=13 is small, do not over-generalize |
| Selection/survivorship transparency | report attempted-but-dropped tools, stall points, the "no new tools" pivot | — | **HAVE** — the evidence audit already surfaces the descoped MCP layer and the human-directed scope contraction; report them |

---

## C7 — Methodology contribution (D2) ★ (the novel transferable result — see `02`)

*Claim:* the differential-parity + round-trip + ablation methodology demonstrably caught bugs
that single-operation/test-based evaluation misses, and yields a quantified what-worked/what-
failed account.

| Metric | Experiment | Status |
|---|---|---|
| **Validation-layer ablation**: bugs uniquely caught per layer (unit / differential-parity / round-trip+metamorphic / fuzz / human) | label the bug corpus by detection layer; compute "escape rate if layer removed" | **NEED ★** — partially reconstructable from git history + `t.Skip`→`t.Fatalf` waves + `UPSTREAM_BUGS.md` |
| **Labeled bug corpus** {tool, class, caught-by, severity, introduced-by, fix-origin} | mine git history + PR bodies + the parity-completion commits | **PARTIAL/NEED ★** — the dataset is half-built; complete it |
| Agent-process metrics (first-pass-compile %, first-pass-parity %, retries) correlated with tool complexity | from transcripts | **NEED** |
| Process-reproducibility statistics (K repeated agent runs on a sample → J/K parity-reached, intervention distribution) | re-run the agent on a few tools K times | **NEED** — honest treatment of nondeterminism (see `04`) |

This is the table the field is missing and the manuscript's strongest methods novelty. See
`02_WHAT_WORKED_WHAT_FAILED.md` for the full design.

---

## Build order (recommended)

**Tier A — unblocks submission (do first):**

1. C5 counterfactual package (non-LLM transpiler + timed human anchor + COCOMO + supervision cost). ★★
2. C2 differential-fuzzing harness + branch-coverage report. ★
3. C7 validation-layer ablation + completed labeled bug corpus. ★
4. C3 large-tier run + pinned hardware spec + median/IQR/CI statistics. ★

**Tier B — strengthens:**
6. C4 real-workflow drop-in demo + flag-compat %.
7. C2 htslib/htscodecs conformance suites + metamorphic relations.
8. C6 cross-tool generality analysis.
9. C7 agent-process metrics + process-reproducibility (K-run) study.

**Tier C — polish:** CIs everywhere; max-FP-deviation reporting; usability metricization;
selection-bias transparency.

Each Tier-A item is small relative to the existing harness — the `pipeline/` infrastructure
already supports adding cells, comparators, and scales; most "NEED"s are new *experiments on
existing machinery*, not new machinery.
