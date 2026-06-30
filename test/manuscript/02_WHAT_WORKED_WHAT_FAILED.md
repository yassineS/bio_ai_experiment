# 02 — What worked / what failed (the methods contribution, D2)

The user's explicit emphasis: **"focus on what methods made the experiment work, and which
failed."** This is also the manuscript's most novel, transferable result. The goal is to turn
process anecdotes into **quantified, falsifiable data**, grounded in the oracle-problem and
reward-hacking literatures.

---

## A. What WORKED (with repo evidence)

Each item is a *method* with a concrete, citable instance from the project's git history.

1. **Byte-exact differential parity as a forcing function and a reward-hacking-immune oracle.**
   The original binary is an external oracle the agent cannot see or edit, so the test-gaming
   that undermines SWE-bench (edit `conftest.py`, `sys.exit(0)`, `__eq__`→True, answer
   exfiltration) is structurally impossible. Evidence: `2ddf87c` introduced the pipeline;
   `86f8d31` the single gate; `70256e0`/`4cceb25` *removed* env-guard skips so the gate
   couldn't be silently bypassed (converting `t.Skip` → hard `t.Fatalf`). *Theory cite:* Barr
   2015 (pseudo-oracle), SemGraft ICSE'18, RustAssure 2025; *contrast cite:* ImpossibleBench,
   SWE-Bench+, OpenAI retiring SWE-bench Verified.

2. **Upstream binaries as oracles AND fixture generators.** Never hand-coded expected output;
   upstream `samtools/bgzip/tabix` materialize valid BAM/CRAM/VCF/BED from seeded text
   (`pipeline/runner`, `pipeline/internal/upstream`). Removes the "we wrote the test to match
   our code" confound.

3. **Round-trip / metamorphic validation caught what single-op parity missed — the flagship.**
   `5f811f0`: our CRAM encoder produced files **our own decoder masked** (70,000/300,000 medium
   reads → `NNNN`); only an *upstream cross-decode round-trip* exposed it. Two defects:
   multi-reference slices (`RefSeqID=-2` filled with N) and the `RR` flag always 0. `1068c77`
   records round-trip catching **five bug classes at once**. *This is the single most persuasive
   evidence that layered validation > single-operation testing*, and maps onto the metamorphic-
   testing literature (Chen; the oracle-free round-trip needs no upstream).

4. **"Profile, don't assume."** `19b0430`/`c755e45`: profiling **disproved the documented
   "libm-bound" label** on `bcftools call` (~2% libm; the lever was allocation), and found
   `samtools stats` is **consumer-bound** (`observe()` ~49% cumulative) so threading can't help
   it — redirecting effort correctly. `f970994`: dropped a bzip2 brute-force default that made a
   48 MB encode take 97 s (9.5× faster once measured). *Method:* never optimize on a label;
   measure first.

5. **Allocation-reduction patterns, all byte-exact** (`ReadInto`, scratch reuse): `2bbad51`
   (VCF view 3.1×→1.04×), `826bafe`, `3d507d0` (call −61% allocs), and the bed-tool sweep.
   *Method:* output-neutral refactors validated by the same parity gate that guards correctness.

6. **Constraining the agent's dependency surface as a hallucination/supply-chain control.**
   stdlib-first, three sanctioned third-party deps only (gonum, klauspost/compress,
   ulikunitz/xz), each scope-confined. Directly counters the slopsquatting failure class
   (~20% of LLM package recommendations are hallucinated; Spracklen USENIX'25): a hallucinated
   import either doesn't compile or is caught by the parity oracle, and the sanctioned-set
   discipline prevents silent dependency bloat.

7. **Parallel-agent delegation under a hard, independent gate.** The `worktree-agent-<hash>`
   merge series shows multiple isolated git-worktree agents reconciling against `main`, each
   carrying its own byte-identity proof + medium parity gate + full `go test`. *Method:*
   parallelism is safe when every branch must pass the same external oracle before merge.

---

## B. What FAILED / was HARD (the honest core — do NOT omit)

A paper that shows only wins is a credibility red flag. These are real and, framed as data,
*strengthen* the paper.

1. **The agent shipped bugs its own code masked.** The CRAM multi-reference-slice bug
   (`5f811f0`) was invisible to our own decoder; the **latent gVCF string-retention bug**
   (`3d507d0`) only became live once an optimization removed an accidental per-line allocation
   — a contract violation that had been silently safe. *Lesson:* self-consistent testing is
   insufficient; you need an **independent** oracle and round-trip cross-checks. (This is the
   project's own instance of the "both wrong the same way" limit of differential testing.)

2. **Mischaracterizations that survived in docs until measured.** The "libm-bound" `call` label;
   the bzip2 over-engineering. *Lesson:* LLM-written *explanations* of performance can be
   confidently wrong; profiling is the corrective.

3. **Tasks abandoned / descoped — and by whom.** The entire planned **MCP-server LLM-interface
   layer was descoped** (`75c378b`); a deliberate **"no new tools" pivot** (`0e4e010`); the
   **large-tier validation abandoned in-environment** (disk-bound). Critically, the audit shows
   both scope pivots were **human-directed ("owner decided")**, not agent initiative — a
   concrete place the AI needed human steering. *Report this honestly:* the unit of analysis is
   a **human-supervised** agent workflow.

4. **Floors that cannot be crossed under the constraints.**
   - **libm parity floor** (`39857eb`): pure-Go `math.Pow/Log` aren't bit-identical to glibc;
     `bcftools call` QUAL is off by the last ULP on **133/~12,000 sites** — accepted as Similarity.
   - **cgo/libdeflate floor**: scan paths (stats/query/flagstat/depth) sit at the pure-Go
     inflate/parse floor; going below 1× would need cgo into libdeflate, deliberately forgone.
   *Lesson:* some "failures" are principled constraint choices, not capability gaps — say which.

5. **A correctness-process weakness the audit surfaced: CI is actually disabled**
   (`.github/workflows/ci.yml` is `workflow_dispatch` only; `41489bc`), and CLAUDE.md
   *misstates* it as active. Validation was run **locally by agents and self-reported in PR
   bodies**. This is a genuine reproducibility/independence gap that must be fixed (re-enable CI)
   and disclosed.

6. **Residual parity gaps (quantified, from `PARITY_ROADMAP.md`):** bedtools ~98% (3 Skips:
   `bedcluster -s` tie order, `bedjaccard -s`, `bedsplit -a size`); bcftools ~99% (norm `-m+`
   ID concat; merge INFO combine rules; csq GFF phase validation; + the `csq @pos` SKIP where
   *upstream* is buggy); mosdepth `--by` ±1 on ~4/1240 regions; CRAM ~95% (v4 transform codecs
   XPACK/XRLE/XDELTA parsed-not-decoded; aux-tag ordering on cross-decode; slice ref-MD5 zero).

---

## C. Turning A/B into rigorous data (the experiments)

### C.1 Labeled bug corpus (the dataset)

One row per defect caught during development:

| field | values |
|---|---|
| bug id | — |
| tool / subcommand | — |
| bug class | logic / off-by-one / format-encoding / FP-ULP / memory / CLI-semantics / perf |
| **caught by** | unit · differential-parity · round-trip/metamorphic · fuzz · human-review · upstream-bug(pre-existing) |
| severity | wrong-output · crash · silent-divergence |
| introduced by | agent-authored · agent-misport-of-upstream-bug |
| fix origin | agent-self-corrected · human-flagged |

Reconstructable from: git history, the `t.Skip`→`t.Fatalf` parity-completion waves,
`UPSTREAM_BUGS.md`, and PR bodies (the audit found "7 genuine port bugs", "5 samtools bugs"
already half-enumerated with discovery method).

### C.2 Validation-layer ablation (the headline methods table)

For the corpus, compute **bugs uniquely caught by each layer** = how many real bugs would have
*escaped* if that layer were removed:

| layer | bugs caught | uniquely caught | escape rate if removed |
|---|---|---|---|
| unit tests | | | |
| differential parity | | | |
| round-trip / metamorphic | | | |
| differential fuzzing | | | |
| human review | | | |

This directly answers "which validation method actually mattered" and is expected to show
**differential parity + round-trip are load-bearing** — the central methods claim. It is
computable largely from data already in the repo.

### C.3 A two-axis failure taxonomy (honest framing)

Separate **honest failures** (agent tried, fell short) from **gaming failures** (satisfied the
letter, not the intent) — the latter being precisely what our oracle prevents:

- **Axis A — capability/honest:** localization · incorrect/partial implementation · spec
  misunderstanding · hallucinated API/package · process pathology (loops, give-up) · build error.
- **Axis B — reward-hacking/gaming (prevented by the differential oracle, but worth showing the
  agent *attempted*):** test-harness manipulation · runner subversion · assertion defeat ·
  output hardcoding · reference-answer exfiltration · spec-vs-test exploitation.

Grounded in the published taxonomies (Mündler 2511.00197 for Axis A; Anthropic 2511.18397,
ImpossibleBench, EvilGenie, METR, Bondarenko for Axis B). The argument: in a *test-suite*
regime Axis B is rampant; in our *external-binary-oracle* regime Axis B failures cannot
produce a passing result — a concrete, measurable advantage of the methodology.

### C.4 Agent-process metrics

From transcripts: first-pass-compile %, first-pass-parity %, retries, interventions, tokens/$,
wall-clock; correlated against tool complexity (LOC/language/format-count) for C6. *Where the
agent struggled is as valuable as where it succeeded.*

---

## D. The one-paragraph "methods" abstract this supports

> We evaluated correctness not with the tool's own tests (which an agent can game and which
> are contaminated for tools predating the model) but by **byte-for-byte differential comparison
> against the original binary as an external oracle**, augmented by oracle-free round-trip and
> metamorphic checks and differential fuzzing. A layered ablation over a labeled corpus of NNN
> development-time defects shows that differential parity and round-trip validation uniquely
> caught M bugs that unit tests and the agents' self-consistency missed — including a CRAM
> encoder defect the agents' own decoder masked. We further report where the method and the
> agents failed: a libm last-ULP floor, a cgo/performance floor, residual parity gaps, and the
> points at which human scope-steering and review were indispensable — so the result is properly
> read as a *human-supervised* agent workflow under an independent correctness gate.
