# Manuscript gap analysis (validation / parity / performance workstream)

Status of the claim-backing experiments **as of the real-data + perf campaign**
(2026-06-25), focused on the validation/parity/performance side. Cross-reference
[`../01_CLAIMS_AND_EXPERIMENTS.md`](../01_CLAIMS_AND_EXPERIMENTS.md) — several of
its `NEED` markers are now **stale** (the experiment has since been built); this
page reconciles them.

## Closed this campaign

- **C2 parity with CIs + real data** — `realparity` on the real GIAB exome:
  **15/15 byte-exact, 0 DIVERGE** ([`realdata/`](realdata/README.md)); parity
  rates now carry Wilson + Clopper-Pearson 95% CIs
  ([`performance_tables.md`](performance_tables.md), `parity_statistics.md`).
- **C3 perf: median + IQR + bootstrap CI** — H1a upgrade landed
  (`pipeline/stats`, `pipeline/bench`); reported per cell across smoke→large +
  real data, slow cells (`mpileup`/`call`/`isec`) flagged plainly; figures in
  [`figures/`](figures/).
- **C3 hardware/scale anchoring** — environment pinned
  ([`hardware.md`](hardware.md)); large tier run per-format-group
  ([`large_tier/`](large_tier/README.md)) with the heavy >12 GB cells OOM-noted.
- **12 correctness + 5 scalability/perf bugs** found via real data and fixed
  byte-exact (CRAM reference-free decode, M5/UR, aux ordering; stats CIGAR;
  depth OOM→indexed; faidx/fqidx; bcftools norm 436×→stream). Strengthens the
  C6 "agent finds + fixes real bugs" narrative and the bug corpus.

## Already covered (claims-doc `NEED` markers that are now stale)

- **Differential fuzzing (C2 ★)** — harness `pipeline/difffuzz` + `cmd/diff-fuzz`
  and results [`differential_fuzzing.md`](differential_fuzzing.md). Present.
- **Flag-compat % (C4)** — [`flag_compat.md`](flag_compat.md) reports per-tool
  flag-acceptance rates. Present (could refresh the numbers).
- **Conformance corpus (C2)** — htscodecs/htslib vectors via
  `pipeline/conformance` ([`conformance_run.txt`](conformance_run.txt)). Present.
- **Round-trip / metamorphic (C2)** — `pipeline/roundtrip`, 14/14 interop at
  every tier. Present.

## Genuine remaining gaps (validation/parity/perf)

| # | Gap | Claim | Effort | Notes |
|---|---|---|---|---|
| G1 | **Go branch coverage of the parity sweep** — what fraction of port branches the oracle actually exercises | C2 | low | `-coverprofile` over the parity matrix; report a single coverage %. The strongest answer to "untested input regions". |
| G2 | **Max abs/rel FP deviation per tool** (not just within-ε) | C2 | low | the similarity comparator already emits per-cell deltas; aggregate the max per tool into a table. |
| G3 | **Pipeline drop-in demo** — swap our binary into a real nf-core/Nextflow (or Snakemake) step, run end-to-end unchanged | C4 ★ | medium | concrete usability evidence; an nf-core samtools/bcftools step is the cleanest. |
| G4 | **`samtools view` region→SAM speed** (~12×) and the remaining ~11× RSS cells | C3 | medium | correct today; perf optimisation (next task). |
| G5 | **GIAB biological concordance** (hap.py/vcfeval F1, GA4GH-stratified) | C2 ★ | high | **deliberately deferred** by the owner ("skip the official GIAB variant-calling test"); the byte-exact differential parity stands in for now. |

## Out of this workstream (separate effort — agent-process / meta)

These are real manuscript gaps but belong to the agent-process/effort side, not
validation: **C5** transpiler counterfactual + human-effort anchor + agent-cost
instrumentation; **C6** labelled bug corpus + validation-layer ablation +
K-run process-reproducibility; and the **CI self-reporting** credibility fix
(R1-T6) — make the upstream-parity job the authoritative, non-self-reported
gate. Flagged here for completeness; tracked in `01`/`05`.

## Recommendation (next, in priority order)

1. **G4 `samtools view` speed** — already queued; closes the last C3 perf
   outlier.
2. **G1 branch coverage** + **G2 max-FP table** — both low-effort, high-rigour;
   directly answer the two sharpest C2 reviewer questions.
3. **G3 Nextflow drop-in** — the single most credible C4 usability result.
4. **G5 GIAB concordance** — when the owner re-scopes it in (highest domain
   credibility, highest cost).
