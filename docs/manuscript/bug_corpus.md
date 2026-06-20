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
> differential fuzzing (`pipeline/difffuzz`) and GIAB concordance are the only probes into it.

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

| id | tool | class | severity | caught_by | fix_origin | evidence |
|---|---|---|---|---|---|---|
| A1 | CRAM encoder (samtools view -C) | format-encoding | **silent-divergence** | **round-trip / cross-decode** | agent-self-corrected | `[evidenced]` `5f811f0`: multi-reference slices (`RefSeqID=-2`→`N`) + `RR` flag always 0; 70k/300k medium reads → `NNNN`; **our own decoder masked it — only an upstream cross-decode round-trip exposed it; single-op parity passed.** The flagship "round-trip > single-op" datum, and a *prospective-checkable* row. |
| A2 | bcftools call (shared VCF reader) | memory | silent-divergence (latent) | differential-parity + `-race` | agent-self-corrected | `[evidenced]` `3d507d0`: gVCF blocker retained `Chrom`/`Ref`/sample strings across reads — a contract violation **only safe because `Text()` had allocated per line; became live when an optimization removed the copy.** Anti-memorization evidence (a regurgitator would not introduce this). |
| A3 | CRAM encoder | perf | (perf regression) | profiling / bench | agent-self-corrected | `[evidenced]` `f970994`: bzip2 brute-forced on every block → 48 MB encode took 97 s (9.5× slower). Caught by measuring, not by a label ("profile, don't assume"). |
| A4 | CRAM encoder | format-encoding | wrong-output (size/ratio) | differential-parity (size) + round-trip | agent-self-corrected | `[evidenced]` `85ba89a`: first encoder wrote literal read bases (45 MB vs upstream 32 MB) — the root of the ×71 regression; reference-based encoding fixed it. |
| A5 | bcftools call (docs/metric) | (mischaracterization) | n/a | profiling | human-flagged | `[evidenced]` `19b0430`/`c755e45`: the "libm-bound" label was wrong (~2% libm; lever was allocation); the `samtools stats` "decode-bound" label was wrong (consumer-bound). Not a code bug — a *wrong explanation* the profile corrected. |
| A6 | bcftools roh/consensus/norm | logic | wrong-output | differential-parity (medium) | agent-self-corrected | `[reconstructed]` the medium full-validation surfaced 9 diverges incl. these; fixed in `a1c426d`. ST/RG interleave, consensus freeze model, `-m+` FILTER union. |
| A7 | skewer 3' adapter | logic | wrong-output | differential-parity (medium) | agent-self-corrected | `[reconstructed]` gap-free detector vs upstream Myers k-difference aligner; one read; fixed by porting `cAdapter::align`. |
| A8 | prinseq | logic | silent-divergence | unit (coverage push) | agent-self-corrected | `[reconstructed]` `2d8a00d`: silent SVG-write-error bug found while raising coverage to ~100%. |
| A9 | vcftools --012 | format-encoding | wrong-output | differential-parity | agent-self-corrected | `[reconstructed]` `abe69e4`: `--012` encoding + monomorphic-site retention. |
| A10 | bedcluster | memory/logic | wrong-output | differential-parity | agent-self-corrected | `[reconstructed]` `98aa9c9`: raw-line retention across records (same class as A2). |

*Pattern to test prospectively:* every **silent-divergence** row (A1, A2, A4, A9, A10) was caught
by **differential parity or round-trip**, not by unit tests — the central methods claim. Re-run each
bug-commit through unit-only vs parity vs round-trip in isolation to convert this from assertion to
measurement.

---

## Corpus U — upstream-found defects (17 entries; the "method finds real bugs" result)

> **Political handling (red-team `05`, R2 §3 / `04 §4c`):** report **only upstream-CONFIRMED**
> items as "bugs"; file PRs/issues and cite the numbers — *the `upstream_confirmed` column is the
> whole ballgame.* Reclassify the rest as **spec-gap** (collaborative) or **intended-but-surprising**
> (not bugs). Reverent tone; engage Heng Li's "AI Rewrite Dilemma" as a collaborator. As of now most
> are **NOT yet filed/confirmed** — that is the top pre-submission task for this corpus.

| id | tool | class | severity | disposition | upstream_confirmed |
|---|---|---|---|---|---|
| U1 | bcftools +check-sparsity `-R` | CLI-semantics | silent-divergence (drops BED lines) | fix-on-port | **TODO-file** |
| U2 | bcftools som write-map (`fwrite` ret) | logic | crash (subcommand dead) | fix-on-port | **TODO-file** |
| U3 | bcftools +remove-overlaps `-m min(QUAL)` ring-index | off-by-one | wrong-output (stale mark) | fix-on-port | **TODO-file** |
| U4 | bcftools +setGT `-n X` NULL in error | logic | wrong-output (msg/UB) | fix-on-port | **TODO-file** |
| U5–U6 | pkg/htsgo/bcf writer (waves 21–22) | format-encoding | silent-divergence (htslib interop) | fix-on-port (our writer, vs BCF spec) | n/a (our bug vs spec) |
| U7 | vcftools `--keep-INFO` semantic | CLI-semantics | silent-wrong-output | fix-on-port (port-side) | n/a (port divergence) |
| U8 | vcftools `--remove-INFO` semantic | CLI-semantics | silent-wrong-output | fix-on-port (port-side) | n/a |
| U9 | vcftools `--pca` jagged M[i] | memory | crash on missing GT | fix-on-port | **TODO-file** |
| U10 | vcftools spill-file 1-byte stack overflow | memory | crash | fix-on-port | **TODO-file** |
| U11 | vcftools `.ifreqburden` INDV label index | off-by-one | wrong-output | fix-on-port | **TODO-file** |
| U12 | vcftools `--hapcount` prev_bin_idx shift | off-by-one | wrong-output | fix-on-port | **TODO-file** |
| U13 | vcftools `--hapcount` end-of-stream read-after-free | memory | crash/UB | fix-on-port | **TODO-file** |
| U14 | vcftools `--hapcount` BED first-line skip | off-by-one | silent-divergence | fix-on-port | **TODO-file** |
| U15 | mosdepth `--by` region mode (missing outputs) | logic | wrong-output | fix-on-port | **TODO-file** |
| S1 | mosdepth overlap-pair detection | (feature gap, NOT upstream bug) | — | track-only | n/a — reclassify as *intended* |
| S2 | vcftools `--site-pi` | (RESOLVED: not a bug) | — | resolved-as-feature | n/a |

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
