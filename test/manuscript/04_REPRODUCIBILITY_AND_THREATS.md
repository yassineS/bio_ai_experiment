# 04 — Reproducibility package & threats to validity

## 1. The artifact package (target: ACM-style Available / Functional / Reproduced)

**Already exists (re-runnable by a reviewer):**

- `pipeline/cmd/parity-pipeline` and `pipeline/cmd/full-validation` (parity matrix + round-trip
  - bench, one gate, non-zero exit on any DIVERGE/ERROR/round-trip FAIL).
- Deterministic **seeded fixture generator** (`pipeline/fixtures`, default seed 1), built from
  vendored upstream tools; 4 scale tiers.
- Bench harness (`wait4` rusage, wall/CPU/peak-RSS).
- 4,665 tests + 52 benchmarks (`go test ./...`); 10 `Fuzz` functions.
- Pinned deps in `go.mod` (gonum v0.17.0, ulikunitz/xz v0.5.15, klauspost/compress v1.18.6).
- Per-PR byte-identity proofs in commit bodies.

**Must add for a credible package (these are `NEED`s):**

1. **Hardware spec** for every perf number (CPU model/cores/RAM/kernel) — currently *unrecorded*.
2. **Systematic upstream version pinning** alongside the metrics (vendored submodule SHAs +
   samtools/bcftools/bedtools versions used) — partially present in prose only.
3. **Model ID + date** for the agent runs (non-negotiable for the process claims).
4. **Container** (Docker/Apptainer) pinning Go toolchain + all upstream binaries + REF data.
5. **Agent transcripts/logs** (prompts, tool calls, retries) — the evidentiary basis for C5/C7.
6. **One-command repro** + a short artifact appendix (install → run → expected output), plus the
   **large-tier** run on a ≥30 GB-scratch node.
7. **Re-enable CI** (`.github/workflows/ci.yml` is currently `workflow_dispatch`-only) so
   validation is independently re-run, not self-reported; fix the stale CLAUDE.md claim.

## 2. The nondeterminism problem (address head-on)

LLM agents are **not deterministic**; re-running will not reproduce the same Go code. This does
**not** break artifact badging, but you must separate two reproducibility claims and never
conflate them:

- **Outcome reproducibility (deterministic, badge-eligible):** *given the committed Go code*, the
  parity/perf/round-trip results reproduce byte-for-byte. The artifact is **the code + harness**,
  not the agent. This is fully achievable and is the `Results Reproduced` target.
- **Process reproducibility (statistical, not deterministic):** the *generation process*
  reproduces only in distribution. Report it as a **process with outcome statistics over K
  repeated runs** (e.g., "across K independent agent runs on tool T, parity reached in J/K, with
  median N interventions"). State plainly: **the agent is the experimental *treatment*, not part
  of the reproducible *artifact*.**

Even LLM *inference* is nondeterministic at temperature 0 (batch-invariance; Thinking Machines
2025), and single-run agent pass rates vary 2.2–6.0 pp (arXiv:2602.07150) — so any agent-side
number is a distribution, not a point.

## 3. Threats to validity (construct / internal / external / conclusion)

### Construct (are we measuring what we claim?)

- **"Usability"/"feasible" unoperationalized** → metricize (flag-compat %, workflow drop-in,
  CWE-class elimination) or demote to qualitative (C4).
- **Byte-exact ≠ correct (both wrong the same way)** → mitigate with metamorphic relations +
  the fix-on-port bug ledger (positive evidence the method catches *real* upstream defects).
- **LOC ≠ effort** → already disciplined; keep it.

### Internal (is the causal attribution sound?) — **the #1 desk-reject risk**

- **No counterfactual** for "LLMs make this feasible" → the C5 control package (non-LLM
  transpiler + timed human anchor + bounded priors). Without it the central claim is unsupported.
- **Human-in-the-loop confound** → results are agent + expert operator; measure supervision cost;
  unit of analysis = the human-supervised workflow.
- **Selection/survivorship bias** → report attempted-but-dropped tools, stall points, the
  human-directed "no new tools" / MCP-descope pivots; state selection criteria.
- **Self-reported validation (CI disabled)** → re-enable CI; disclose that historical validation
  was agent-self-reported in PR bodies.

- **★ Memorization / solution contamination (added by red-team R1).** samtools/bcftools/bedtools
  source is in every training corpus; the agent may have **reproduced memorized upstream code**,
  and *byte-exact parity is precisely what memorization predicts*. The oracle-immunity argument
  defeats reward-*hacking* but not solution-*contamination*. This is the single most important
  missing threat for an ML-systems venue. Mitigate: (a) name it explicitly; (b) foreground the
  **genuine novel bugs the agent introduced** (gVCF retention, CRAM multi-ref) as anti-regurgitation
  evidence — a pure copier would not introduce them; (c) stratify parity by **spec-only vs
  upstream-source-derived** ports (only spec-only ports escape both memorization and the
  shared-misconception critique); (d) report Go-vs-C structural/n-gram overlap (you can translate
  memorized C, so this is necessary not sufficient).
- **★ Differential-testing independence is violated for source-derived ports (R1).** Tools ported
  by reading the upstream C are not independently derived → "both wrong the same way" risk is
  *elevated*. Report parity stratified by derivation; lean on GIAB + htslib's own conformance
  corpus as the genuinely independent oracles.
- **Comparator construct validity (R1).** "Provenance-stripping" before byte-comparison is a
  researcher degree of freedom — over-stripping can manufacture parity. Specify and justify exactly
  what is normalized (timestamps, `@PG`, version strings, hash-ordering) and show it hides no real
  divergence.
- **Oracle-version monoculture (R1, R2).** Byte-exact is defined against *one pinned* upstream
  version, and samtools/bcftools output drifts across releases. Pin the SHA *and* demonstrate parity
  is stable across ≥2 upstream versions for a sample, so a reviewer cannot suspect tuning to one
  convenient version.

### External (does it generalize?)

- **Bounded, format-heavy, well-specified tools with an executable reference** = a near-ideal
  oracle setting. Scope the claim to "tools with an executable reference and a checkable spec";
  do **not** generalize to "all legacy software." n=13 families is small — report per-tool.
- **Model/time dependence** → pin and disclose the model+date; frame as a capability snapshot.

### Conclusion (are the statistics sound?)

- **Perf: min-of-N on one unspecified box, medium-only** → median+IQR/Hodges-Lehmann + CI, ≥5
  reps, cold/warm disclosed, multi-scale incl. large, hardware pinned; **report the slower cells**.
- **Pass rates without CIs / coverage denominator** → Wilson/Clopper-Pearson CIs; always pair a
  rate with the coverage of the space it is over.
- **Self-reported metric drift** (e.g. non-test LOC 218k in docs vs 223k actual; call ×2.05 vs
  ×2.17 vs ×2.2 across prose) → recompute all headline numbers from a single scripted pass for
  the manuscript; a careful reviewer will notice hand-maintained drift.

## 4. Top desk-reject reasons → pre-emptions (checklist)

1. **No counterfactual for the comparative claim.** → C5 control package. *Highest priority.*
2. **Overclaiming reproducibility of a nondeterministic process.** → §2 split.
3. **Cherry-picked successes, no failure data.** → `02 §B` + the failure taxonomy + bug corpus.
4. **Vague "usability"/"better".** → C4 metrics or demotion.
5. **Differential-testing limits unacknowledged.** → metamorphic + fuzzing + coverage + GIAB.
6. **AI-productivity inflation** (ignoring supervision cost; LOC-as-effort). → C5 + existing discipline.
7. **Numbers not independently reproducible** (no hardware spec, CI disabled, medium-only). → §1.
8. **Memorization not addressed** (byte-exactness = what regurgitation predicts). → §3 + the
   novel-bugs-as-anti-regurgitation argument + spec-only-vs-source-derived stratification.
9. **Self-graded validation** (CI `workflow_dispatch`-only; agents self-reported parity in PR
   bodies — back in the SWE-bench self-grading regime the paper claims to escape). → re-enable CI,
   re-execute the full suite independently, report those as canonical, disclose the history.

## 4b. The maintenance / trust answer (the reply to Heng Li's "AI Rewrite Dilemma")

The objection that, unanswered, reduces the paper to a stunt: *a new AI-generated codebase is a
new maintenance liability; validation+maintenance, not generation, is the unsolved problem.* The
paper's genuine answer, and it should be a named subsection:

- **Oracle-pinned continuous re-validation *is* the maintenance model.** You never trust the Go
  code on its own; you trust it because it is **continuously diffed against the canonical C tool on
  every commit, forever, automatically.** This reframes the artifact from "unvalidatable AI output"
  to "a codebase pinned to a living oracle" — and makes **re-enabling CI thesis-critical, not
  housekeeping** (the embodiment of the answer; CI is currently disabled — fix it).
- **Concede the stewardship point honestly:** the transferable contribution is the *methodology*,
  which outlives any codebase; adoption decisions belong to the community; do not promise indefinite
  stewardship no one believes.
- **Engage Heng Li by name and on his terms** — cite "The AI Rewrite Dilemma" in the intro as the
  framing tension and structure a discussion thread as a respectful *response* (evidence on
  validation; honest position on maintenance), advancing his concern rather than refuting it.

## 4c. The bug-finding claim — handle as humble upstream collaboration (R2)

The "differential porting surfaced N latent defects in tools the field depends on" result is
powerful *and* politically dangerous (the maintainers will read/referee this). Requirements:
**report only upstream-CONFIRMED defects** (filed PR/issue links — the "upstream accepted" column
is the whole ballgame); **reclassify** the rest into *spec-gap* (collaborative: "we surfaced a spec
ambiguity") or *intended-but-surprising* (not bugs); reverent tone (these codecs are hard; 15 years
of correctness is a testament, not sloppiness); **pin one audited number**. Several current
fix-on-port SKIPs (`csq @pos`, `norm -m+` ID concat, `bedjaccard -s`) may be intended behavior or
spec ambiguity — do not call them "bugs" until upstream agrees.

## 5. Ethics / dual-use / disclosure notes

- **Disclose AI authorship** prominently (the code is agent-generated under human direction);
  this is itself part of the contribution, not something to hide.
- **Maintenance/sustainability honesty** (Heng Li's concern): a new codebase is a new
  maintenance burden; discuss the long-term stewardship plan, not just the one-time port.
- **Safety framing is a feature:** memory-safety + constrained dependency surface directly
  mitigate the CVE-class and slopsquatting risks the related-work section documents.
