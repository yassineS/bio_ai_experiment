# Manuscript gap analysis (validation / parity / performance workstream)

Status of the claim-backing experiments **as of the real-data + perf campaign**
(updated 2026-07-06), focused on the validation/parity/performance side.
Cross-reference [`../01_CLAIMS_AND_EXPERIMENTS.md`](../01_CLAIMS_AND_EXPERIMENTS.md)
— several of its `NEED` markers are now **stale** (the experiment has since been
built); this page reconciles them.

## chr20 realbench headline (real GIAB HG002/GRCh38, cloud tier)

The real-data `realbench` sweep on the GIAB HG002/GRCh38 **chr20** tier
(`pipeline/realbench`, via the `test/nextflow/` Seqera pipeline on AWS Batch)
scored **PASS = 129, DIFF = 2, ERROR = 0, SKIP = 1** across the samtools /
bcftools / bedtools / seqtk cell matrix. The two DIFFs are the **accepted
`consensus` libm last-ULP** residuals (documented in
[`max_fp_deviation.md`](max_fp_deviation.md) — base/seq/qual bytes are
byte-exact, only the derived `cq` score column differs); the single SKIP is
`bcftools csq` (an ours-only feature with no upstream oracle cell). Zero
byte-exact ERRORs. This drove ~10 real-data bug fixes this campaign (see the
[bug corpus](../bug_corpus.md) and the "Closed this campaign" list below).

The **exome + wgs 60× tiers** (whole-genome, the heavy tier) are **in progress**
(Seqera run `mS3IH42QfGTWO`, `realbench-exome-wgs-60x`); their PASS/DIFF/ERROR/SKIP
will be folded into these tables by a follow-up once the run SUCCEEDs. They do
not block the chr20 headline above.

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
  ([`large_tier/`](large_tier/README.md)). The heavy `mpileup`/`call`/`isec`
  cells now have real large numbers: the earlier "OOM" was **scratch-disk
  exhaustion** (a ~17 GB temp output on a 5 GB-free overlay), not memory — peak
  RSS is bounded (106/18/464 MB ours vs 42/10/11 MB upstream). Re-run with a
  big-disk `TMPDIR`, all three complete (see [`figures/`](figures/README.md)).
- **12 correctness + 5 scalability/perf bugs** found via real data and fixed
  byte-exact (CRAM reference-free decode, M5/UR, aux ordering; stats CIGAR;
  depth OOM→indexed; faidx/fqidx; bcftools norm 436×→stream). Strengthens the
  C6 "agent finds + fixes real bugs" narrative and the bug corpus.

### Closed in the chr20 realbench + parity/perf follow-up (2026-07 campaign)

- **~10 more real-data bugs fixed byte-exact** on the chr20 realbench sweep:
  `samtools fixmate` singleton+TLEN, `markdup` dup-selection + peak RSS
  **1.4 GiB → 18.5 MiB**, CRAM partial-reference `M5`, `prinseq` stable-sort,
  `bcftools gtcheck` multi-allelic; plus several `bed*` cells un-SKIPped into real
  upstream parity (pairtopair/pairtobed/tobam/split, and the RNG cells
  bedrandom/bedshuffle/bedsample under `-seed`), the `skewer` oracle build, and
  `mosdepth`. See the [bug corpus](../bug_corpus.md) (rows A26–A35).
- **G1 — Go statement coverage of the parity sweep = 64.25 %**
  ([`branch_coverage.md`](branch_coverage.md)); per-tool breakdown, the strongest
  answer to "which input regions are untested".
- **G2 — max abs/rel FP deviation per tool**
  ([`max_fp_deviation.md`](max_fp_deviation.md)): byte-exact **0.0** for the vast
  majority; three documented last-ULP residuals (`bcftools call -m` QUAL,
  `bedgenomecov` histogram fraction, `samtools consensus` gap5 `cq`).
- **G4 — `samtools view` region/BED→SAM** ([`view_speed.md`](view_speed.md)):
  the historical ~12× outlier is already ~2.5×; this pass cut hot-loop allocations
  23–37× and peak RSS ~15 %, with the residual wall gap documented as the
  irreducible pure-Go inflate floor.
- **G3 — nf-core drop-in demo** ([`nfcore_dropin/`](nfcore_dropin/README.md)): our
  `samtools`/`flagstat` swapped unchanged into a real Nextflow module, end-to-end.
- **G7 — streaming ByteExact-stdout parity compare** (PR #461): the runner's
  stdout compare now streams both children through a bounded digester
  (256 MiB → 16 MiB heap); genome-scale cells no longer OOM. (This is the gap the
  table below tracks as **G6**.)

## Already covered (claims-doc `NEED` markers that are now stale)

- **Differential fuzzing (C2 ★)** — harness `pipeline/difffuzz` + `cmd/diff-fuzz`
  and results [`differential_fuzzing.md`](differential_fuzzing.md). Present.
- **Flag-compat % (C4)** — [`flag_compat.md`](flag_compat.md) reports per-tool
  flag-acceptance rates. Present (could refresh the numbers).
- **Conformance corpus (C2)** — htscodecs/htslib vectors via
  `pipeline/conformance` ([`conformance_run.txt`](conformance_run.txt)). Present.
- **Round-trip / metamorphic (C2)** — `pipeline/roundtrip`, 14/14 interop at
  every tier. Present.
- **Parity coverage % + max-FP-deviation (C2)** — now materialised as concrete
  artifacts: [`branch_coverage.md`](branch_coverage.md) (statement coverage
  64.25 %, per tool) and [`max_fp_deviation.md`](max_fp_deviation.md) (per-tool
  max abs/rel deviation; byte-exact 0.0 for most). These retire the
  `parity_statistics.md` "mechanism-and-bound only" caveat with measured numbers.

## Genuine remaining gaps (validation/parity/perf)

**All six validation/parity/perf gaps below are now closed** (see the "Closed
this campaign" lists above); the table is retained as the record. (Note: the
campaign labelled the G6 streaming-compare work "G7"; it is the same gap.)
**GIAB biological concordance (hap.py / vcfeval variant-calling accuracy) is not
a gap here and never was** — it is **out of scope** for this project, which
validates byte-exact parity against the upstream oracles, not variant-calling
concordance. It is not deferred, not planned, and not part of the manuscript.

| # | Gap | Claim | Effort | Notes |
|---|---|---|---|---|
| G1 | **Go statement coverage of the parity sweep — DONE** | C2 | done | **64.25 %** across the parity-exercised code, per-tool breakdown → [`branch_coverage.md`](branch_coverage.md). |
| G2 | **Max abs/rel FP deviation per tool — DONE** | C2 | done | byte-exact **0.0** for the vast majority; three documented last-ULP residuals → [`max_fp_deviation.md`](max_fp_deviation.md). |
| G3 | **Pipeline drop-in demo — DONE** | C4 ★ | done | our `samtools`/`flagstat` swapped unchanged into a real nf-core-style Nextflow module, end-to-end → [`nfcore_dropin/`](nfcore_dropin/README.md). |
| G4 | **`samtools view` region/BED→SAM speed — DONE (documented floor)** | C3 | done | the ~12× was already ~2.5×; this pass cut hot-loop allocations 23–37× and peak RSS ~15 %; residual wall gap is the irreducible pure-Go inflate floor → [`view_speed.md`](view_speed.md). |
| G5 | **`bcftools isec` multi-sample FORMAT over-decode — FIXED** | C3 | done | memory fixed (streaming k-way merge + byte-bounded batch: RSS 39×→~6.5×); the apparent ~15× wall was a slow-bind-mount + under-buffered-write artifact (fixed, 256 KiB buffers); and the residual FORMAT over-decode is now fixed too — the vcf reader's `KeepRawSamples` mode keeps the FORMAT+sample columns verbatim (`RawTail`) instead of parsing every sample into per-sample maps and rebuilding them on write, which is byte-identical for well-formed records. isec wall **2.0× → 1.33× upstream** (0.56 s → 0.37 s on the 16-sample fixture), byte-exact vs upstream incl. PL. |
| G6 | **Parity-matrix streaming compare — FIXED** | C2 | done | the `pipeline/runner` ByteExact-stdout compare (the path the ~17 GB heavy cells `sam_mpileup`/`bcf_call`/`bcf_isec` take) now streams **both** children's stdout straight through a `StreamDigester` (provenance-strip + running md5 + 64 KiB head) and compares digests — memory is O(64 KiB)/side, never O(output), so it no longer OOMs. This reuses the exact streaming normaliser `realparity`/`realbench` already use (`timedRunStreaming` + `CompareDigests` in `pipeline/runner`). Verdicts are byte-exact-identical to the old buffered `CompareByteExact` (digest equality over the stripped stream = stream equality); proven by `TestCompareDigestsMatchesCompareByteExact` (every provenance fixture pair, streaming verdict == buffered verdict) and `TestStreamDigestBoundedMemory` (256 MiB stream, heap grows <16 MiB). The decode/similarity/output-file modes (which inherently need full bytes for a subprocess decoder or field-by-field compare) keep the buffered path; they are not the genome-scale OOM concern. |

## Out of this workstream (separate effort — agent-process / meta)

These are real manuscript gaps but belong to the agent-process/effort side, not
validation: **C5** transpiler counterfactual + human-effort anchor + agent-cost
instrumentation; **C6** labelled bug corpus + validation-layer ablation +
K-run process-reproducibility; and the **CI self-reporting** credibility fix
(R1-T6) — make the upstream-parity job the authoritative, non-self-reported
gate. Flagged here for completeness; tracked in `01`/`05`.

## Recommendation (next, in priority order)

All six validation/parity/perf gaps (G1–G6) are closed. The remaining
manuscript-facing work is:

1. **Fold in the exome + wgs 60× realbench tiers** once Seqera run
   `mS3IH42QfGTWO` (`realbench-exome-wgs-60x`) SUCCEEDs — pull
   `realbench.exome.*` / `realbench.wgs.*` from
   `s3://realbench-209855136877-work/work/results` and add their
   PASS/DIFF/ERROR/SKIP alongside the chr20 headline above.
2. **Agent-process / meta gaps** (out of this workstream): C5 human-effort
   anchor + agent-cost instrumentation; C6 second-rater κ on the bug corpus;
   the CI self-reporting credibility fix. Tracked in `01`/`05`.
