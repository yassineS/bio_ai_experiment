# 05 — Reviewer red-team (internal)

Three independent internal reviewers critiqued the plan (00–04) from different chairs:
**R1** a skeptical empirical-SE / ML-systems methodologist (top SE/ML venue); **R2** a
bioinformatics domain referee (Nature Methods / Genome Biology); **R3** a handling editor
(Nature Methods / Nature Computational Science triage). This file records their verdicts, the
**new threats they surfaced** (folded into `04`), and the consolidated action list.

## Verdicts

| Reviewer | Verdict | One-line |
|---|---|---|
| R1 (methods) | Reject-as-is → accept with `★` experiments | "Unusually honest plan; the path to acceptance is executing the experiments, not more red-teaming. Fix the CI-self-reporting hole first — it silently invalidates the independence story." |
| R2 (domain) | Major revision / not yet submittable | "Stop selling 'we cloned 53 tools' (a biologist shrugs); sell 'we solved the trust/validation problem that makes such clones believable.' GIAB + named silent-corruption edge cases are make-or-break." |
| R3 (editor) | Revise-then-review | "Reframe feasibility→verifiability (free, removes the desk-reject trigger); two papers, split by contribution; sprint D2 for priority — RustQC is the scoop threat; GIAB determines D1's venue tier." |

**Consensus:** the artifact and the oracle idea are strong and the self-critique is unusually
honest; but (a) the causal "feasibility" claim is unsupported, (b) the load-bearing experiments
(GIAB, counterfactual, ablation-as-data, fuzzing, large-tier) are unrun, and (c) two process
holes (CI disabled/self-reported; memorization not addressed) currently undercut the trust thesis.

## New threats they surfaced (now added to `04 §3`/`§6`)

1. **★ Memorization / solution contamination (R1-T4) — the most important missing threat.**
   samtools/bcftools source is in every training corpus; the agent may have **reproduced
   memorized upstream code**, and *byte-exact parity is exactly what memorization predicts*. The
   oracle-immunity argument defeats reward-*hacking* but not solution-*contamination*. Must be a
   named threat. Counter-evidence: the **genuine novel bugs the agent introduced** (gVCF
   retention, CRAM multi-ref) are anti-regurgitation evidence (a pure copier wouldn't introduce
   them); plus stratify parity by **spec-only vs upstream-source-derived** ports (only spec-only
   ports escape both the memorization and the shared-misconception critiques).

2. **★ CI disabled / self-reported validation (R1-T6, R2, R3).** `ci.yml` is
   `workflow_dispatch`-only and validation was run *by the agents* and self-reported in PR bodies
   — i.e., the gate the paper claims is "ungameable" was executed and reported by the agent, back
   in the SWE-bench self-grading regime. **Fix:** re-enable CI; re-run the *entire* parity +
   round-trip + bench suite in clean CI from committed code on specified hardware; report *those*
   as canonical; disclose the history. R2 notes this is *also* the literal embodiment of the
   maintenance answer to Heng Li (oracle-pinned continuous re-validation).

3. **Differential-testing independence is violated for source-derived ports (R1-T1).** Tools
   ported by reading the upstream C are **not independently derived**, so "both wrong the same
   way" risk is *elevated*, not negligible. Report parity stratified by derivation; lean on GIAB
   + htslib's own conformance corpus as genuinely independent oracles.

4. **The planned counterfactual is rigged (R1-T2).** Beating a rule-based transpiler (c2go) is a
   straw-man — everyone knows it produces non-idiomatic output. With n=1–2 unpowered human
   anchors, *no* comparative effort claim is supportable. **Resolution (all three agree): drop
   "feasibility/cheaper" as a headline; reframe to "verifiability"** and present human anchors as
   descriptive calibration only, with the counterfactual named as future work.

5. **The ablation is currently post-hoc storytelling (R1-§4).** Retrospective bug labeling by the
   authors invites confirmation bias; "escape rate if layer removed" is a counterfactual over a
   world that didn't run; n≈20 has no statistical power. **Make it rigorous:** pre-register the
   labeling schema; **independent second labeler + Cohen's κ**; make the ablation **prospective**
   (actually re-run the codebase-at-the-bug-commit through each validation layer in isolation and
   *measure* which catches it — the CRAM case proves this is doable); report as a **case series**,
   not a statistic. Mine transcripts for *actual* Axis-B gaming attempts — if zero exist, downgrade
   "our oracle prevents gaming" to "structurally cannot game" (don't claim measured prevention).

6. **n overstated (R1-T7).** 53 CLIs ≈ **1–2 independent technical cores** (htslib-core,
   bedtools-core) + thin CLI shells. Report results at the level of independent cores; "53 CLIs"
   stays as a deployment fact, not 53 independent successes.

7. **Survivorship at the bug level (R1-§6).** The caught-bug corpus excludes bugs that escaped
   *all* layers; the ablation measures relative layer value among *detected* bugs and is silent on
   the undetected tail. Differential fuzzing + GIAB are the only probes into that tail — state this.

8. **Comparator construct validity (R1-§6).** "Provenance-stripping" is a researcher degree of
   freedom; over-stripping manufactures parity. Specify and justify the normalization exactly.

9. **Oracle-version monoculture (R1-§6, R2).** Byte-exact is defined against *one* pinned upstream
   version, and samtools/bcftools output drifts across releases. Pin the SHA *and* show parity is
   stable across ≥2 upstream versions for a sample.

10. **The bug-finding claim is politically dangerous (R2-§3).** This community is small; the
    maintainers (Heng Li, Petr Daněček, …) will read/referee this. **Only report
    upstream-CONFIRMED defects** (filed PRs/issues with links — the "upstream accepted" column is
    the whole ballgame); reclassify the rest into **spec-gap** (collaborative framing) or
    **intended-but-surprising** (not bugs). Reverent tone; engage Heng Li's "AI Rewrite Dilemma"
    *with* him. Pin one audited number.

## Consensus reframes (apply to the whole plan)

+ **Thesis: feasibility → verifiability.** Lead with "an external-binary differential oracle makes
  LLM-agent re-implementation *verifiable* — immune to the gaming/contamination undermining
  current LLM-coding evals." Verifiability is provable; feasibility is not (yet). *(Applied to `00`.)*
+ **Venue (R2 reorders R3):** primary **Nature Computational Science** (native register for the
  AI-methodology framing; prestigious without demanding "new biology"), co-equal/fallback
  **Genome Biology**, safety net **GigaScience**; companion **SE/ML** paper for D2. Nature Methods
  only with a crisp "enables new biology" hook, which neither reviewer currently sees. *(Applied to `00 §4`.)*
+ **Two papers, split by contribution axis, co-cited (all three).** D1 (artifact → bio venue) for
  impact; **D2 (method) is the perishable novelty — sprint it for priority** (mostly
  reconstructable from existing repo data). Timing risk: **Seqera RustQC** (same domain, industrial
  backing) is the existential scoop.
+ **Maintenance answer = oracle-pinned continuous re-validation (R2-§5).** The reply to Heng Li:
  you never trust the Go code alone; you trust it because it is continuously diffed against the
  canonical C tool on every commit — which makes re-enabling CI thesis-critical, not housekeeping.
  Add a named stewardship subsection. *(Applied to `04 §5`.)*
+ **Memory-safety "why it matters" needs a constituency (R2-§1):** clinical/CLIA-CAP pipelines,
  biobank-scale, shared HPC — parsing *untrusted* BAM/CRAM/VCF in memory-unsafe C. Back it with an
  **OSS-Fuzz htslib CVE/crash table**. *(Strengthen C4.)*

## Consolidated, prioritized action list (supersedes `01` build-order where they conflict)

**P0 — removes desk-reject triggers (cheap/fast):**

+ Reframe title+abstract feasibility→verifiability; reorder venue to NComSci primary.
+ Re-enable CI; re-run full suite in clean CI; recompute *every* headline number from one scripted
  pass (kills metric drift); pin model ID+date + hardware spec.
+ Add the 4 missing threats (memorization, shared-derivation, version-monoculture, comparator) to `04`.
+ Reclassify the bug ledger; file upstream PRs; report only confirmed.

**P1 — the `★` experiments (the science):**

+ GIAB concordance to full GA4GH standard **incl. CMRG + difficult-region stratifications**,
  hap.py *and* vcfeval, vs *upstream* bcftools, ULP-floor proven to never flip a genotype/PASS.
+ Named silent-corruption edge-case battery as discrete results (CRAM M5/REF_CACHE/multi-ref;
  norm Number=A/R/G re-indexing; `.bai/.csi/.tbi` byte-identity; sort stability at `-@1`; MD/NM).
+ Differential fuzzer diffing **stdout+stderr+exit-code** + branch-coverage report.
+ Ablation as a **prospective, independently-labeled (κ) case series**; mine transcripts for gaming.
+ Counterfactual: keep only as descriptive human-anchor calibration + supervision-cost accounting;
  drop the comparative claim.
+ Large-tier perf run on a pinned node; median+IQR/Hodges-Lehmann + CIs; report regressions.

**P2 — strengthens:**

+ htslib `test/` + htscodecs-corpus conformance (**promote to Tier A** per R2).
+ Real end-to-end workflow drop-in (nf-core/Snakemake module repoint) + flag-compat %.
+ Process-reproducibility K-run study on a sample; cross-tool generality at the *core* level.
+ Wilson/Clopper-Pearson CIs + coverage denominators everywhere; max-FP-deviation per tool.

## What the reviewers agreed the plan got right

The self-awareness (most objections were already in the plan); the related-work positioning; the
oracle idea ("the agent cannot cheat an oracle it does not control"); the honest "what failed"
section; the two-paper instinct; the "what not to claim" discipline. R1: "the path to acceptance
is not more red-teaming; it's executing the `★` experiments and resisting the two over-claims you
keep almost making — that byte-exactness equals trustworthiness (run GIAB) and that the artifact
equals a feasibility result (drop the cost claim or power it)."
