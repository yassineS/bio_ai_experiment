# Manuscript plan — `bio_ai_experiment`

Working title (candidate):
**"Byte-exact modernization of legacy bioinformatics software with LLM agents:
a differential-parity methodology and a 53-tool case study."**

This directory is the **planning substrate** for a high-profile manuscript on the
feasibility of modernizing critical, aging bioinformatics command-line software
with LLM coding agents. It is the output of a structured research+brainstorm pass
(an internal repo-evidence audit + ~12 independent expert-research briefs + a
reviewer red-team). It is **not** the manuscript; it is the spec the manuscript and
its remaining experiments are built from.

## Read in this order

| File | What it gives you |
|---|---|
| [`00_STRATEGY.md`](00_STRATEGY.md) | Thesis, the defensible novelty wedge, target-venue analysis (honest), narrative arc, the one thing that gets it desk-rejected. |
| [`01_CLAIMS_AND_EXPERIMENTS.md`](01_CLAIMS_AND_EXPERIMENTS.md) | The atomic claims, and for each: metric → experiment → baseline → **have vs. need**. This is the **tests-and-metrics to build** list. |
| [`02_WHAT_WORKED_WHAT_FAILED.md`](02_WHAT_WORKED_WHAT_FAILED.md) | The methods contribution the project most wants: the validation-layer **ablation** + the labeled **bug corpus**, turned from anecdote into data. |
| [`03_RELATED_WORK.md`](03_RELATED_WORK.md) | Positioning against prior art + an annotated, verification-graded bibliography. |
| [`04_REPRODUCIBILITY_AND_THREATS.md`](04_REPRODUCIBILITY_AND_THREATS.md) | The artifact package, the nondeterminism problem, threats to validity (incl. memorization, self-graded-CI, oracle independence), the maintenance answer to Heng Li, and desk-reject pre-emptions. |
| [`05_REVIEWER_RED_TEAM.md`](05_REVIEWER_RED_TEAM.md) | Three independent internal-reviewer reports (skeptical methods, bioinformatics domain, editor), the new threats they surfaced, and the consolidated **P0/P1/P2 action list**. **Read this for the verdict.** |

## One-paragraph summary

Every prior high-profile result on LLM software modernization (Google's
int32→int64 migration, Amazon Q's "4,500 developer-years," IBM watsonx Code
Assistant for Z) is an **industry experience report on proprietary code with a
weak, self-estimated counterfactual**, and the dominant LLM-coding evaluations
(SWE-bench and kin) grade by **test-suite pass**, an oracle now shown to be
gameable and contaminated. This project occupies the empty quadrant: **public,
ubiquitous scientific tools anyone can re-verify, judged by byte-exact behavioral
equivalence against the original tool as an external oracle the agent cannot
edit.** As of mid-2026 no peer-reviewed work demonstrates byte-exact whole-tool
parity for a widely-used tool ecosystem; the closest cases (wedeo, RustQC,
MirrorCode, Meta's ProgramBench) are non-peer-reviewed and/or partial. The
manuscript's contribution is therefore twofold: (1) the **artifact** — 53
drop-in Go CLIs at near-perfect byte-exact parity, memory-safe and cgo-free, at
competitive-or-better performance — and (2) the **method** — a reusable
differential-parity + round-trip + ablation methodology that demonstrably caught
bugs single-operation testing misses, and a quantified account of *which
validation layers and agent practices actually worked and which failed.* The
single biggest gap to close before submission is a **counterfactual control** for
the causal "feasibility/effort" claim.

## Headline guidance (post-red-team)

- **Reframe the thesis from "feasibility" to "verifiability"** — the provable, novel claim is that an
  external reference-binary oracle makes LLM re-implementation *verifiable* (immune to the
  gaming/contamination undermining LLM-coding evals); "LLMs make it cheaper" is unsupported without a
  counterfactual and should not be the headline.
- **Primary venue: Nature Computational Science**; fallback Genome Biology; safety net GigaScience;
  **companion SE/ML paper for the method, sprinted for priority** (Seqera RustQC is the scoop risk).
- **The five make-or-break experiments** (all currently unrun): GIAB concordance (incl. CMRG +
  difficult regions), the named silent-corruption edge-case battery, differential fuzzing, the
  ablation-as-rigorous-data, and the C5 counterfactual (or drop the effort claim).
- **Fix first, cheaply:** re-enable CI + re-execute the suite independently; recompute every headline
  number once; reclassify the bug ledger to upstream-confirmed only; add the memorization threat.

## Status of this plan

- Built from: `docs/METRICS.md`, `docs/FULL_VALIDATION.md`, `docs/PARITY_ROADMAP.md`,
  `PROJECT_STATUS.md`, `pipeline/`, git history, ~12 external expert-research briefs, and a
  3-referee internal red-team (`05`).
- All external numbers cited here carry a verification grade in `03_RELATED_WORK.md`;
  many primary PDFs were egress-blocked during research and are flagged
  **[verify-PDF]** — confirm before they enter the manuscript.
- This is a living document. Update it as the "need" experiments in
  `01_CLAIMS_AND_EXPERIMENTS.md` are completed.
