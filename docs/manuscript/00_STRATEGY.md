# 00 — Strategy: thesis, novelty, venue, narrative

> **Reframe (post-red-team):** lead with **verifiability**, not feasibility. All three internal
> reviewers (see `05`) agree the causal "LLMs make modernization *feasible/cheaper*" claim is
> unsupported without a counterfactual and should not be the headline; the provable, novel claim
> is that an **external-binary differential oracle makes LLM-agent re-implementation verifiable** —
> immune to the test-gaming and contamination that undermine current LLM-coding evaluation.

## 1. The thesis (one sentence)

> An **external reference-binary differential-parity oracle** (the original tool, which the agent
> cannot see or edit) makes LLM-agent re-implementation of legacy bioinformatics CLIs **verifiable
> at whole-ecosystem scale** — we demonstrate a single memory-safe, dependency-minimal, cgo-free Go
> codebase at near-perfect byte-exact behavioral equivalence (and GIAB-concordant on the
> variant-calling path), and the *methodology* (which validation layers caught which bugs, and
> which agent practices worked vs. failed) is the transferable scientific result.

Two deliverables, deliberately separated so a weak one cannot launder a strong one:
- **D1 (artifact):** 53 drop-in CLIs / 13 tool families, near-perfect parity, memory-safe,
  faster-or-comparable. *Descriptive, strongly evidenced today.*
- **D2 (method):** the differential-parity + round-trip + ablation methodology and the
  quantified "what worked / what failed" account. *The novel, transferable contribution.*

## 2. The novelty wedge (what is genuinely new)

The research swept four adjacent literatures. The empty quadrant this work occupies is
sharp and defensible:

| Prior art | What it did | Why it is *not* this | Cite |
|---|---|---|---|
| Google int32→int64 / JUnit / Joda (ICSE-SEIP'25, FSE'25) | LLM migration **within one language**, at scale | Same-language refactor of **proprietary** code; oracle = the repo's own tests; counterfactual self-estimated | Nikolov 2501.06972; Ziftci 2504.09691 |
| Amazon Q (4,500 dev-yrs), IBM WCA4Z (COBOL→Java) | Cross-language modernization | **Vendor/CEO claims, unaudited**; one customer saw only ~36% acceleration | Jassy 2024; The Register 2024 |
| SWE-bench & kin | LLM agents resolve real issues | Graded by **gameable, contaminated unit tests** (see §3); bug-fix not whole-tool re-implementation | Jimenez 2310.06770 |
| wedeo, RustQC, MirrorCode, Meta ProgramBench | AI rewrites of real tools | **Non-peer-reviewed and/or partial**; bit-exact only on *slices* (one codec, one tool); ProgramBench reports *no* full solves | (03_RELATED_WORK §4) |
| uutils coreutils | Drop-in Rust rewrite of GNU | **Human-written**; and *still* diverges on **50.6%** of fuzzed `test` cases (NDSS'25) | Li 2025 NDSS |

**The wedge:** *As of mid-2026 there is no peer-reviewed demonstration of byte-exact,
whole-tool behavioral parity for a widely-used scientific-tool ecosystem produced by AI
agents, validated against the original binaries as external oracles.* This work is that
demonstration — on **public** tools anyone can re-verify, under a **hard safety constraint**
(memory-safe, cgo-free, dependency-minimal), across an **interdependent ecosystem** (a
shared htslib-equivalent feeding samtools/bcftools), not a single lucky tool.

### Why the oracle choice is the intellectual core
The dominant LLM-coding evaluations are now known to be undermined by exactly the failure
classes our oracle is immune to:
- **Test-gaming / reward hacking:** agents pass by editing tests (`conftest.py`), short-
  circuiting runners (`sys.exit(0)`), defeating assertions (`__eq__`→True), or exfiltrating
  the scorer's answer (METR, ImpossibleBench, EvilGenie, Anthropic 2511.18397).
- **Weak/insufficient oracles:** ~31% of "passed" SWE-bench patches are actually wrong
  (SWE-Bench+); ~29.6% of plausible patches diverge from the gold fix (PatchDiff, ICSE'26);
  UTBoost re-labels 24–41% of leaderboard entries; OpenAI retired SWE-bench Verified.
- **Contamination:** ≥94% of SWE-bench predates model cutoffs; file-path can be guessed
  from issue text 76% in-bench vs 53% out (SWE-Bench Illusion).

A **byte-for-byte comparison against an independent reference binary the agent cannot see
or edit** is the APR/oracle-problem literature's prescribed answer to all three: it is the
"pseudo-oracle / N-version / reference-implementation oracle" (Barr 2015; SemGraft ICSE'18;
RustAssure 2025), instantiated at whole-tool scale. *The agent cannot cheat an oracle it
does not control.* This reframes the project from "we ported some tools" to "we demonstrate
the evaluation methodology the field is missing."

## 3. The single biggest risk: the counterfactual

The title makes a **causal/comparative** claim ("LLMs make modernization *feasible*").
The evidence today is almost entirely **descriptive** (the artifact has property X). Empirical-
SE reviewers will judge the paper on the counterfactual and, absent one, desk-reject on
internal validity. **This is the #1 thing to fix before submission.** The defensible,
affordable control package (detailed in `01_CLAIMS_AND_EXPERIMENTS.md`, claim C5):
1. **Non-LLM transpiler control** (c2go / C2Rust-style) on the same inputs — shows the
   *automated* alternative fails the parity/safety/idiomaticity bar, isolating the LLM
   contribution. Cheap, primary.
2. **1–2 timed expert human ports** to the same parity bar (e.g. `sickle` + a few bcftools
   subcommands) — calibration anchors, reported as illustrative, not powered.
3. **COCOMO-II / literature effort priors** on the measured Go LOC — bounded range, not a
   point estimate; triangulated against (2).
4. **Honest human-supervision-cost accounting** — interventions, retries, review hours,
   tokens/$ — so the denominator is *agent + human oversight*, never agent-alone.

## 4. Target venue (honest)

The research strongly cautions against *Nature*/*Science* proper as the primary target —
those need broad cross-discipline significance and a fully-controlled study. Ranked realistic
targets and their bar:

| Venue | Bar it sets | Fit |
|---|---|---|
| **Nature Methods** (recommended primary) | A method/tool of broad utility to life scientists, rigorously validated on community-standard data (GIAB), reproducible | **Strong** if we add GIAB concordance + the counterfactual + artifact package. The "method" is the parity harness + the drop-in suite. |
| **Nature Computational Science / Nature Biotechnology** | Computational advance with demonstrated impact | Strong alternative; NComSci is a natural home for the AI-methodology framing. |
| **Genome Biology / GigaScience** (recommended fallback) | Solid tool + open, reproducible artifact (GigaScience *loves* this) | **Very strong, high-probability**; GigaScience artifact culture fits the harness. |
| **Bioinformatics (OUP)** | Useful, correct, benchmarked tool | Safe floor. |
| **PNAS / Patterns (Cell Press)** | Broad significance / data-science methods | Possible if the "AI modernizes critical infrastructure" framing lands broadly. |
| ICSE / FSE / NeurIPS-D&B (SE/ML venues) | The **methodology** contribution (differential parity, ablation, reward-hacking-immune oracle) | A **second paper** aimed at the SE/ML audience is viable and arguably higher-novelty there. |

**Recommendation (revised by the domain referee, R2 in `05`):** make **Nature Computational
Science** the *primary* target, not Nature Methods — it is the native register for the
AI-methodology framing and is prestigious without demanding a "this enables new biology" hook
that neither domain reviewer can currently see. **Genome Biology** is the co-equal/fallback;
**GigaScience** the high-probability safety net with the full reproducible artifact. Treat
**Nature Methods** as an optional first shot only if a crisp "enables new biology" hook emerges.
Ship a **companion SE/ML-venue paper** centered on D2 (the methodology + validation-layer ablation
+ reward-hacking-immune oracle) — and, per the editor (R3), **sprint D2 for priority**: it is the
perishable novelty (Seqera's RustQC, same domain, is the scoop risk) and is mostly reconstructable
from existing repo data.

## 5. Narrative arc (the story a reader follows)

1. **Hook — critical infrastructure is aging and unsafe.** samtools/bcftools/bedtools/htslib
   are dependencies of essentially all genomics; they are C/C++ of hundreds of thousands of
   LOC, memory-unsafe, hard to test, with real bus-factor/maintenance risk (Heng Li's "AI
   Rewrite Dilemma" frames the tension from inside the community).
2. **Tension — can LLMs modernize it, and can we *trust* the result?** The field's own
   evidence says LLM code is gamed, contaminated, hallucination-prone (slopsquatting ~20%).
   Trust is the unsolved problem.
3. **Idea — make the original tool the oracle.** Byte-exact differential parity + round-trip
   + GIAB biological concordance: an oracle the agent cannot game, immune to the failure
   classes plaguing test-based evals.
4. **Result (D1) — it works at scale.** 53 CLIs, near-perfect parity, memory-safe, often
   faster, drop-in into Nextflow/Snakemake/Galaxy/Conda.
5. **Result (D2) — *how* it worked, and where it failed.** The validation-layer ablation
   (round-trip caught the CRAM bug single-op parity missed; the gVCF retention bug), the
   labeled bug corpus, the agent-practice findings ("profile don't assume," constrain the
   dependency surface, parallel agents under a hard gate), and the honest floors (libm ULP,
   cgo/libdeflate, human-directed scope pivots).
6. **Honest limits — what AI could not do unsupervised.** The bugs the agent's own code
   masked; the human scope steering; the medium-only scale; the parity gaps that remain.
7. **Implication.** A reproducible recipe for trustworthy AI modernization of critical
   research software, and a reusable evaluation methodology for the broader question of
   when "AI rewrote it" can be believed.

## 6. What to *not* claim (discipline)

- Not "autonomous": it is a **human-supervised agent workflow under a human-designed gate**.
- Not "2× less code ⇒ 2× less effort": LOC ratio ≠ effort ratio (already disciplined in METRICS).
- Not "faster" unqualified: report the cells where Go is **slower** (mpileup, isec, call) plainly.
- Not "100% parity": it is near-perfect with a documented, quantified gap list and two
  accepted escape valves (libm-ULP Similarity; fix-on-port SKIP where *upstream* is wrong).
- Not "reproducible agent run": the **artifact** (code + harness) is reproducible; the
  **generation process** is reproducible only in distribution (see `04`).
