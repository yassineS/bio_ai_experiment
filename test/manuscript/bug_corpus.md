# Bug corpus — the validation-layer-ablation dataset (C7 / `02 §C`)

This is the labeled dataset behind the manuscript's methods contribution: *which validation
layer caught which defect.* It has two halves:

- **Corpus U — upstream-found defects** (the originals are wrong): evidence that differential
  porting under an external oracle finds **real bugs in tools the field depends on**. Source:
  `docs/UPSTREAM_BUGS.md` (17 entries).
- **Corpus A — agent-introduced defects** (our Go code was wrong): the **validation-layer
  ablation** dataset — what caught the agent's own bugs, and what would have escaped each layer.
  Source: git history + commit bodies.

> **Methodological honesty (per internal red-team `05`, R1 §4).** The "caught-by" field below is,
> for most rows, **reconstructed post hoc by the authors** — a confirmation-bias risk. Before this
> becomes a manuscript table it MUST be: (1) **independently re-labeled** by a second rater with
> **Cohen's κ** reported; (2) made **prospective** where possible — actually re-run the
> codebase-at-the-bug-commit through each validation layer in isolation and *measure* which catches
> it (the CRAM row proves this is doable); (3) reported as a **case series**, not dressed in CIs
> (n is ~dozens). Rows are tagged `[evidenced]` (the commit/PR body states the detection layer) or
> `[reconstructed]` (inferred — needs the independent label). **Survivorship caveat:** this corpus
> contains only *detected* bugs; bugs that escaped every layer are by definition absent, so the
> ablation measures *relative layer value among caught bugs* and is silent on the undetected tail —
> differential fuzzing (`pipeline/difffuzz`) is the only probe into it.

Labeling schema (one row per defect):
`{id, tool/subcommand, class, severity, caught_by, introduced_by, fix_origin, disposition,
upstream_confirmed, evidence}`

- **class:** logic / off-by-one / format-encoding / FP-ULP / memory / CLI-semantics / perf
- **severity:** wrong-output / crash / silent-divergence
- **caught_by:** unit · differential-parity · round-trip/metamorphic · fuzz · human-review · upstream-bug
- **introduced_by:** agent-authored · agent-misport-of-upstream-bug · upstream(pre-existing)
- **fix_origin:** agent-self-corrected · human-flagged

---

## Corpus A — agent-introduced defects (the ablation dataset)

Columns: `id | tool | subcommand | class | severity | caught_by | introduced_by |
fix_origin | evidence`. Every row's `introduced_by` is **agent-authored** unless
noted; A1–A10 / A19–A25 are all `agent-self-corrected`.

| id | tool | subcommand | class | severity | caught_by | introduced_by | fix_origin | evidence |
|---|---|---|---|---|---|---|---|---|
| A1 | samtools | view -C (CRAM encoder) | format-encoding | **silent-divergence** | **round-trip / cross-decode** | agent-authored | agent-self-corrected | `[evidenced]` `5f811f0`: multi-reference slices (`RefSeqID=-2`→`N`) + `RR` flag always 0; 70k/300k medium reads → `NNNN`; **our own decoder masked it — only an upstream cross-decode round-trip exposed it; single-op parity passed.** The flagship "round-trip > single-op" datum, and a *prospective-checkable* row. |
| A2 | bcftools | call (shared VCF reader) | memory | silent-divergence (latent) | differential-parity + `-race` | agent-authored | agent-self-corrected | `[evidenced]` `3d507d0`: gVCF blocker retained `Chrom`/`Ref`/sample strings across reads — a contract violation **only safe because `Text()` had allocated per line; became live when an optimization removed the copy.** Anti-memorization evidence (a regurgitator would not introduce this). |
| A3 | samtools | view -C (CRAM encoder) | perf | perf-regression | profiling / bench | agent-authored | agent-self-corrected | `[evidenced]` `f970994`: bzip2 brute-forced on every block → 48 MB encode took 97 s (9.5× slower). Caught by measuring, not by a label ("profile, don't assume"). |
| A4 | samtools | view -C (CRAM encoder) | format-encoding | wrong-output (size/ratio) | differential-parity (size) + round-trip | agent-authored | agent-self-corrected | `[evidenced]` `85ba89a`: first encoder wrote literal read bases (45 MB vs upstream 32 MB) — the root of the ×71 regression; reference-based encoding fixed it. |
| A5 | bcftools | call (docs/metric) | mischaracterization | n/a (not a code bug) | profiling | agent-authored | human-flagged | `[evidenced]` `19b0430`/`c755e45`: the "libm-bound" label was wrong (~2% libm; lever was allocation); the `samtools stats` "decode-bound" label was wrong (consumer-bound). Not a code bug — a *wrong explanation* the profile corrected. |
| A6 | bcftools | roh / consensus / norm | logic | wrong-output | differential-parity (medium) | agent-authored | agent-self-corrected | `[reconstructed]` the medium full-validation surfaced 9 diverges incl. these; fixed in `a1c426d`. ST/RG interleave, consensus freeze model, `-m+` FILTER union. |
| A7 | skewer | (3' adapter trim) | logic | wrong-output | differential-parity (medium) | agent-authored | agent-self-corrected | `[reconstructed]` `30d7e2c`: gap-free detector vs upstream Myers k-difference aligner; one read; fixed by porting `cAdapter::align`. |
| A8 | prinseq | (SVG graph write) | logic | silent-divergence | unit (coverage push) | agent-authored | agent-self-corrected | `[reconstructed]` `2d8a00d`: silent SVG-write-error bug found while raising coverage to ~100%. |
| A9 | vcftools | --012 | format-encoding | wrong-output | differential-parity | agent-authored | agent-self-corrected | `[reconstructed]` `abe69e4`/`afeb195`: `--012` encoding + monomorphic-site retention; temp-file off-by-one. |
| A10 | bedtools | cluster (bedcluster) | memory/logic | wrong-output | differential-parity | agent-authored | agent-self-corrected | `[reconstructed]` `98aa9c9`: raw-line retention across records (same class as A2). |

### Newly surfaced by the wired-up test batteries (2026-06-20)

These were found the day the `pipeline/difffuzz`, `pipeline/conformance`, and `pipeline/edgecases`
harnesses were built — concrete evidence the methodology keeps finding real divergences. Each has a
`t.Skip("PARITY GAP: …")` regression guard that flips to a hard failure once fixed.

All A11–A18 are `introduced_by` **agent-authored** and `fix_origin`
**agent-self-corrected** (the harness flagged; the agent fixed each on the same
branch). Columns: `id | tool | subcommand | class | severity | caught_by | evidence`.

| id | tool | subcommand | class | severity | caught_by | evidence |
|---|---|---|---|---|---|---|
| A11 | bcftools | norm `-m-`/`-m+` | format-encoding | **silent-divergence (HIGH)** | **edge-case battery** | FORMAT `Number=R`/`G` (AD/PL) vectors **not re-indexed** on split/join → per-allele depths/likelihoods corrupted (`AD 10,5,3` kept verbatim vs upstream `10,5`/`10,3`). `pipeline/edgecases/bcftools_norm_test.go`. |
| A12 | samtools/tabix | index (.bai/.csi/.tbi) | format-encoding | wrong-output (byte) | edge-case battery | Writer did not reproduce htslib's khash bin ordering or `compress_binning` (small bins fold into their parent), so a low-occupancy bin such as 4696 was emitted instead of being merged → not byte-identical. Now ports khash + `compress_binning`, always emits the BAI `n_no_coor` trailer, and adjusts the BAM-CSI depth via `hts_adjust_csi_settings` (also corrected the CSI meta pseudo-bin to `n_bins+1`). `.bai` byte-identical; `.csi`/`.tbi` payloads (decompressed) byte-identical to `samtools index`/`tabix`. `index_identity_test.go`. |
| A13 | samtools | view -C (CRAM encoder) | logic | wrong-output (rejects valid) | conformance (htslib `test/`) | rejects no-SEQ / past-ref-end / unknown-ref records upstream round-trips. `htslib_sam_test.go`. |
| A14 | samtools | (sam reader) | off-by-one/overflow | wrong-output | conformance (htslib `test/`) | POS/PNEXT parsed as int32 → long refs (>2³¹) rejected. `sam_reader.go:101`. |
| A15 | samtools | (CLI input gate) | CLI-semantics | silent-divergence | conformance + fuzz | empty/headerless file accepted (exit 0) vs upstream exit 1. |
| A16 | bcftools | view | FP/format | wrong-output | **differential fuzz** | large QUAL printed verbatim vs htslib `%g` scientific (`4.29497e+09`). |
| A17 | samtools | flagstat | logic | wrong-output | differential fuzz | mate-to-different-chr miscount on odd `RNEXT`. |
| A18 | bedtools | merge | CLI-semantics | wrong-output (accepts malformed) | differential fuzz | accepts inputs upstream rejects (inconsistent fields; unsorted). |

**Resolution status (2026-06-20).** All eight were fixed against upstream as byte-exact (each test's
`t.Skip("PARITY GAP")` guard flipped to an active assertion / a regression test added), on branch
`claude/testing-pipeline`:

- **A11** bcftools norm AD/PL (Number=R/G) re-indexing — FIXED (HIGH-severity silent corruption).
- **A12** index byte-identity — FIXED (khash bin order + `compress_binning`; `.bai` byte-identical,
  `.csi`/`.tbi` payload-identical).
- **A13** CRAM no-SEQ / out-of-bounds / unknown-ref encode+decode — FIXED.
- **A14** SAM long-ref `POS`/`PNEXT`/`TLEN` int32 overflow — FIXED (migrated to int64 `hts_pos_t`
  across `Record` + SAM/BAM/CRAM paths; BAM writer errors cleanly beyond int32, matching the
  on-disk format limit; byte-exact for POS < 2^31).
- **A15** empty/headerless file accepted vs upstream exit 1 — FIXED (reject when the first line is
  not a valid 11-column SAM record; valid headerless-record and zero-record-with-header unaffected).
- **A16** VCF/BCF float `%g` formatting (large/extreme magnitudes) — FIXED.
- **A17** samtools flagstat mate-to-different-chr RNEXT miscount — FIXED.
- **A18** bedtools merge unsorted/ragged-input rejection — FIXED.

**All eight gaps the validation harnesses surfaced are now closed.**

### Surfaced by the parity-skip elimination wave (2026-06-14, PRs #286–#296)

`PROJECT_STATUS.md` records that the wave that drove every remaining
feature-parity `t.Skip` to a hard byte-for-byte (or value-exact) `t.Fatalf`
against the vendored upstream binaries surfaced **"7 genuine port bugs"** along
the way. These are **layer-2 (differential-parity) catches by construction**: the
defect only became visible when the skip was converted into a live diff against
the upstream oracle. All are `agent-authored` / `agent-self-corrected`. The
status note names six distinct defects; the seventh is the second of the two
`bedgenomecov` defects (the per-base offset and the CIGAR-`D` split are separate
fixes), giving the stated count of seven. Source: `PROJECT_STATUS.md` (the
2026-06-14 wave paragraph) — this is the **citation for the count**; the
per-row split below is `[reconstructed]` from that paragraph's prose.

| id | tool | subcommand | class | severity | caught_by | evidence |
|---|---|---|---|---|---|---|
| A19 | bedtools | genomecov `-dz` | off-by-one | wrong-output | differential-parity | `[reconstructed]` per-base `-dz` 0-based offset wrong vs upstream. |
| A20 | bedtools | genomecov (CIGAR `D`) | logic | wrong-output | differential-parity | `[reconstructed]` CIGAR-`D` block splitting not matched on BAM input. |
| A21 | samtools | import | logic | wrong-output | differential-parity | `[reconstructed]` `FMUNMAP` flag mis-set on `/1`,`/2` FASTQ suffixes. |
| A22 | samtools | calmd | format-encoding | wrong-output | differential-parity | `[reconstructed]` MD/NM aux re-append ordering diverged from upstream. |
| A23 | bcftools | norm `-c` | CLI-semantics | wrong-output | differential-parity | `[reconstructed]` `-c` letter semantics (`s`=set/fix, `x`=exclude) mis-mapped. |
| A24 | bcftools | concat `-a` | logic | wrong-output | differential-parity | `[reconstructed]` `-a` contig ordering diverged from upstream. |
| A25 | (7th wave bug — see note) | — | — | wrong-output | differential-parity | `[reconstructed]` the wave reports **7** port bugs; six are individually named above (A19–A24, counting the two `bedgenomecov` fixes separately). This row is a **placeholder for the count's 7th** and is the one row whose subcommand attribution the source prose does not pin — included for honesty, **excluded from the per-layer unique-catch tally** to avoid double-counting an unidentified defect. |

These seven are all **layer-2 unique catches** for the ablation: each was masked
by a `t.Skip` and surfaced *only* when the differential-parity assertion went
live. None had a unit test that failed (the skip suppressed the parity assertion,
not a unit assertion), and none required round-trip or fuzz to expose.

These also serve the **anti-memorization** argument (`04 §3`): a pure regurgitator of upstream C
would not produce these *divergences from* upstream — they are genuine independent-implementation
defects, evidence the agent re-derived rather than copied.

*Pattern to test prospectively:* every **silent-divergence** row (A1, A2, A4, A9, A10, A11, A15) was
caught by **differential parity / round-trip / the conformance+fuzz batteries**, not by unit tests —
the central methods claim. Re-run each
bug-commit through unit-only vs parity vs round-trip in isolation to convert this from assertion to
measurement.

### Surfaced by the chr20 realbench sweep (2026-07, real GIAB HG002/GRCh38)

The real-data `realbench` sweep on the GIAB HG002/GRCh38 **chr20** tier (via the
`test/nextflow/` Seqera pipeline, upstream binaries as the oracle) scored
**PASS = 129, DIFF = 2, ERROR = 0, SKIP = 1** (the 2 DIFFs are the accepted
`consensus` libm last-ULP `cq` residual; the SKIP is ours-only `bcftools csq`).
Getting there surfaced and fixed the following real-data divergences — all
**layer-2 (differential-parity) catches against the upstream oracle on real
whole-chromosome data**, all `agent-authored` / `agent-self-corrected`. Columns:
`id | tool | subcommand | class | severity | caught_by | evidence`.

| id | tool | subcommand | class | severity | caught_by | evidence |
|---|---|---|---|---|---|---|
| A26 | samtools | fixmate | logic | wrong-output | realbench differential-parity | singleton handling + `TLEN` (`isize`) computation diverged from upstream on real reads; now byte-exact. |
| A27 | samtools | markdup | logic | wrong-output | realbench differential-parity | duplicate-selection differed from upstream (which representative read is kept); fixed byte-exact. |
| A28 | samtools | markdup | memory/perf | perf-regression | profiling / realbench | buffered the whole read set → peak RSS **1.4 GiB**; rewritten to bound the working set at **18.5 MiB**, output byte-identical. |
| A29 | samtools | view -C (CRAM encoder) | format-encoding | wrong-output | realbench differential-parity | partial-reference `@SQ` `M5` (a contig only partly covered by the supplied reference) computed/omitted differently from upstream; fixed byte-exact. |
| A30 | prinseq | (record ordering) | logic | wrong-output | realbench differential-parity | output ordering not stable vs upstream on ties; switched to a stable sort → byte-exact. |
| A31 | bcftools | gtcheck | logic | wrong-output | realbench differential-parity | multi-allelic site handling diverged from upstream; fixed byte-exact. |
| A32 | bedtools | pairtopair / pairtobed / tobam / split | logic | wrong-output | realbench differential-parity | these cells had been `t.Skip`-guarded; converting to live real-data parity surfaced and fixed divergences (`bedtag` compare also corrected to decode SAM, not raw BAM stdout, removing a false DIFF). |
| A33 | bedtools | random / shuffle / sample | logic | wrong-output (nondeterminism) | realbench differential-parity | RNG cells were skipped as nondeterministic; made deterministic under `-seed` and un-SKIPped into real upstream parity. |
| A34 | skewer | (oracle build) | build/tooling | n/a (harness) | realbench | the upstream `skewer` oracle failed to build in the realbench image (missing `-c` in the fallback `CXXFLAGS`); fixed so the cell runs a real comparison rather than skipping. |
| A35 | mosdepth | (real-data parity) | logic | wrong-output | realbench differential-parity | real-data divergence surfaced by the sweep; fixed byte-exact (upstream `mosdepth` oracle bundled into the realbench image). |

These are **layer-2 unique catches** in the same sense as the parity-skip wave:
each surfaced only when a real-data cell was compared byte-for-byte against the
upstream oracle (several were `t.Skip`-guarded until this sweep converted them to
live assertions). They are **not** folded into the ablation denominator below (to
keep the carefully-bookkept 24-row tally stable and comparable); they are recorded
here as the chr20-realbench addendum and cited by
[`results/MANUSCRIPT_GAPS.md`](results/MANUSCRIPT_GAPS.md).

### Per-layer summary (Corpus A) — the ablation cross-reference

This is the table [`results/ablation.md`](results/ablation.md) consumes. "Unique
catch" = the layer is the *only* one in our stack that flagged that defect, so it
would have escaped to the next-weakest probe (or to the undetected tail) had the
layer been removed. The denominator is the **24 attributable** Corpus-A rows
(A1–A24; A25 is the unidentified 7th wave bug, **excluded** from all counts;
A3 perf-regression and A5 mischaracterization are **not** output-correctness
defects and are tallied separately, see notes).

| layer | defects this layer caught | uniquely caught | rows (unique) |
|---|---|---|---|
| 1 unit | 1 | 1 | A8 |
| 2 differential-parity | 15 | 13 | A2, A6, A7, A9, A10, A11, A12, A19, A20, A21, A22, A23, A24 |
| 3 round-trip / metamorphic | 2 | 1 | A1 (unique); A4 (shared with layer 2 on size) |
| 4 differential fuzzing | 4 | 4 | A15, A16, A17, A18 |
| 5 human review | 1 | 1 | A5 (mischaracterization, not a code bug) |
| conformance (htslib `test/`) | 3 | 2 | A13, A14 (unique); A15 (shared with fuzz) |
| profiling / bench | 2 | 1 | A3 (perf, unique); A5 (shared with human review) |

Notes on the bookkeeping (read before quoting any number):

- **Layer 2 is credited 13 unique** output-correctness catches: A2, A6, A7, A9,
  A10, A11, A12, and the six named wave bugs A19–A24. A4 is **shared** (the
  size/ratio divergence was visible to both single-op parity and round-trip), so
  it is not a *unique* layer-2 catch and not a unique layer-3 catch.
- **Conformance** (htslib `test/` fixtures) is a *distinct oracle* from our
  seeded differential-parity sweep; A13 and A14 were caught there and nowhere
  else. We report conformance as its own line rather than folding it into
  "differential-parity," because removing the seeded parity sweep would **not**
  have caught A13/A14 — only the spec's own corpus did. A15 was caught by both
  conformance and fuzz.
- **A1 is the single most load-bearing row**: round-trip is the *only* layer that
  caught it (single-op parity passed; our own decoder masked it). Removing
  round-trip ⇒ A1 escapes.
- **A3 (perf) and A5 (mischaracterization)** are not byte-correctness defects;
  they are kept in the corpus for completeness but are excluded from the
  output-correctness escape-rate denominator in `ablation.md`.
- Every `[reconstructed]` row's layer attribution is inferred from commit/PR
  prose, not from a prospective re-run; see the methods caveat at the top of this
  file and in `ablation.md`.

---

## Corpus U — upstream-found defects (17 entries; the "method finds real bugs" result)

> **Political handling (red-team `05`, R2 §3 / `04 §4c`):** report **only upstream-CONFIRMED**
> items as "bugs"; file PRs/issues and cite the numbers — *the `upstream_confirmed` column is the
> whole ballgame.* Reclassify the rest as **spec-gap** (collaborative) or **intended-but-surprising**
> (not bugs). Reverent tone; engage Heng Li's "AI Rewrite Dilemma" as a collaborator. As of now most
> are **NOT yet filed/confirmed** — that is the top pre-submission task for this corpus.

Columns: `id | tool/subcommand | class | severity | caught_by | introduced_by |
fix_origin | disposition | upstream_confirmed`. For Corpus U the **`caught_by`
layer is differential porting under the upstream oracle** (reading/diffing the
upstream source + binary while building the port) — i.e. these are not catches of
*our* code by a validation layer but defects in *upstream* found by the act of
re-implementation; `introduced_by` is `upstream(pre-existing)` for genuine
upstream bugs and `agent-authored (port-side)` for the port-divergence rows
(U5–U8). `fix_origin` is `agent-self-corrected` throughout (each was fixed on the
porting branch).

| id | tool/subcommand | class | severity | caught_by | introduced_by | fix_origin | disposition | upstream_confirmed |
|---|---|---|---|---|---|---|---|---|
| U1 | bcftools +check-sparsity `-R` | CLI-semantics | silent-divergence (drops BED lines) | differential porting | upstream(pre-existing) | agent-self-corrected | fix-on-port | **TODO-file** |
| U2 | bcftools som write-map (`fwrite` ret) | logic | crash (subcommand dead) | differential porting | upstream(pre-existing) | agent-self-corrected | fix-on-port | **TODO-file** |
| U3 | bcftools +remove-overlaps `-m min(QUAL)` ring-index | off-by-one | wrong-output (stale mark) | differential porting | upstream(pre-existing) | agent-self-corrected | fix-on-port | **TODO-file** |
| U4 | bcftools +setGT `-n X` NULL in error | logic | wrong-output (msg/UB) | differential porting | upstream(pre-existing) | agent-self-corrected | fix-on-port | **TODO-file** |
| U5–U6 | pkg/htsgo/bcf writer (waves 21–22) | format-encoding | silent-divergence (htslib interop) | round-trip / interop | agent-authored (port-side) | agent-self-corrected | fix-on-port (our writer, vs BCF spec) | n/a (our bug vs spec) |
| U7 | vcftools `--keep-INFO` semantic | CLI-semantics | silent-wrong-output | differential-parity | agent-authored (port-side) | agent-self-corrected | fix-on-port (port-side) | n/a (port divergence) |
| U8 | vcftools `--remove-INFO` semantic | CLI-semantics | silent-wrong-output | differential-parity | agent-authored (port-side) | agent-self-corrected | fix-on-port (port-side) | n/a |
| U9 | vcftools `--pca` jagged M[i] | memory | crash on missing GT | differential porting | upstream(pre-existing) | agent-self-corrected | fix-on-port | **TODO-file** |
| U10 | vcftools spill-file 1-byte stack overflow | memory | crash | differential porting | upstream(pre-existing) | agent-self-corrected | fix-on-port | **TODO-file** |
| U11 | vcftools `.ifreqburden` INDV label index | off-by-one | wrong-output | differential porting | upstream(pre-existing) | agent-self-corrected | fix-on-port | **TODO-file** |
| U12 | vcftools `--hapcount` prev_bin_idx shift | off-by-one | wrong-output | differential porting | upstream(pre-existing) | agent-self-corrected | fix-on-port | **TODO-file** |
| U13 | vcftools `--hapcount` end-of-stream read-after-free | memory | crash/UB | differential porting | upstream(pre-existing) | agent-self-corrected | fix-on-port | **TODO-file** |
| U14 | vcftools `--hapcount` BED first-line skip | off-by-one | silent-divergence | differential porting | upstream(pre-existing) | agent-self-corrected | fix-on-port | **TODO-file** |
| U15 | mosdepth `--by` region mode (missing outputs) | logic | wrong-output | differential-parity | agent-authored (port-side) | agent-self-corrected | fix-on-port | **TODO-file** |
| S1 | mosdepth overlap-pair detection | (feature gap, NOT upstream bug) | — | n/a | n/a | n/a | track-only | n/a — reclassify as *intended* |
| S2 | vcftools `--site-pi` | (RESOLVED: not a bug) | — | n/a | n/a | n/a | resolved-as-feature | n/a |

*Note:* several vcftools entries (U9–U14) are **genuine memory-safety defects** (crashes,
read-after-free, stack overflow) in a widely-used C++ tool — exactly the "memory-safety matters"
constituency the domain reviewer wanted foregrounded. These are the strongest U-rows *if confirmed
upstream*.

---

## Pre-submission TODO for this corpus

1. **File upstream PRs/issues** for the confirmable U-rows (esp. the memory-safety ones U9–U14, U2);
   record issue numbers in `upstream_confirmed`. A branch `claude/csq-upstream-bug-report` already
   exists — extend that workflow.
2. **Independent re-labeling** of `caught_by` (Corpus A) by a second rater; report Cohen's κ.
3. **Prospective ablation:** for each Corpus-A bug, checkout the bug-commit's parent, run the
   codebase through {unit, differential-parity, round-trip, fuzz} in isolation, record which catch
   it. Produce the `02 §C.2` ablation table from *measurement*, not assertion.
4. **Reclassify** S1 (intended feature gap) and the port-side rows (U5–U8) out of the "upstream bug"
   count; pin **one audited upstream-bug number**.
5. Cross-reference each row to its `parity_test.go` / commit so a reviewer can re-run it.
