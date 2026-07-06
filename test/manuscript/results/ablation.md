# Validation-layer ablation (claim C7)

The manuscript's strongest methods-novelty result: for a labeled corpus of
development-time defects, *which validation layer uniquely caught which bug* —
i.e. how many real defects would have **escaped** if a given layer were removed
from the stack. This file is the analysis; the labeled dataset it consumes is
[`../bug_corpus.md`](../bug_corpus.md) (Corpus A, the agent-introduced defects,
and its per-layer summary table). Read that file's methodological-honesty header
first — its caveats bind every number here.

> **One-line honest summary.** Of 22 attributable agent-introduced
> output-correctness defects, **differential parity uniquely caught 13** and is by
> far the load-bearing layer; **round-trip/metamorphic uniquely caught the single
> most important bug** (a CRAM defect our own decoder masked); **differential
> fuzzing uniquely caught 3**; the **spec's own conformance corpus uniquely caught
> 2**; **unit tests uniquely caught 1**; and **human review caught 0 code bugs**
> (its one catch was a wrong *explanation*, not wrong output). Most per-row layer
> attributions are **RECONSTRUCTED** from git/PR prose, not measured by a
> prospective re-run — see §Provenance.

---

## The validation stack (the five layers under ablation)

1. **Unit / compile tests** — in-package `*_test.go`, hand-written expectations,
   table-driven; the layer an agent can in principle game.
2. **Differential byte-exact parity** — diff our binary's output against the
   *upstream binary* on seeded fixtures, provenance-stripped, byte-for-byte (or
   bounded-tolerance for the libm-ULP path). The reward-hacking-immune external
   oracle. Includes the `t.Skip`→`t.Fatalf` parity-completion waves.
3. **Round-trip / metamorphic** — oracle-free relations: encode→decode→re-encode,
   and *cross-decode* (our encoder → upstream decoder and vice-versa). Catches the
   "our own decoder masks our own encoder" failure that single-op parity cannot.
4. **Differential fuzzing** — feed both binaries the same mutated input; diff
   stdout + stderr + exit code (`pipeline/difffuzz`). Probes the untested-input
   tail.
5. **Human review** — owner/author inspection; in this corpus it surfaced a
   wrong performance *explanation*, not a code defect.

A sixth probe, the **spec's own conformance corpus** (htslib `test/` fixtures +
htscodecs vectors, `pipeline/conformance`), is a *distinct oracle* from layer 2's
seeded parity sweep, so we score it on its own line rather than folding it into
"differential parity": removing the seeded sweep would not have caught the two
defects only the spec corpus exposed. **Profiling/bench** is likewise scored
separately for the two non-correctness defects (perf regression, mischaracterized
label).

---

## Methods — how each defect was attributed to a layer

Each Corpus-A row carries a `caught_by` field. We attributed it as follows, in
priority order:

1. **`[evidenced]` rows (A1–A5).** The detection layer is stated in the fixing
   commit body or PR description (e.g. `5f811f0` for A1 explicitly says the
   single-op parity passed and only the cross-decode round-trip exposed it). These
   attributions are **MEASURED** in the weak sense that a contemporaneous artifact
   names the layer — but not in the strong sense of a prospective re-run.
2. **`[reconstructed]` rows (A6–A10, A19–A24).** The layer is *inferred* from the
   commit/PR prose plus the structural fact that the defect surfaced during a
   known activity (the medium full-validation run for A6; the
   `t.Skip`→`t.Fatalf` parity-completion wave for A19–A24, which by construction
   converts a suppressed parity assertion into a live differential-parity diff).
   These are **RECONSTRUCTED**, with a confirmation-bias risk the corpus header
   flags.
3. **Harness-found rows (A11–A18).** The detection layer is the named harness file
   that first flagged the divergence (`pipeline/edgecases/…`,
   `pipeline/conformance/…`, `pipeline/difffuzz`). Because the harness *run* is the
   artifact, these are **MEASURED** at the layer granularity (the run output is in
   [`difffuzz_run.txt`](difffuzz_run.txt) and [`conformance_run.txt`](conformance_run.txt)),
   though we did not separately re-run *each* defect through *every other* layer to
   confirm the catch is unique — uniqueness is argued structurally (e.g. a fuzz-only
   input region has no fixture in the seeded parity sweep).

**"Uniquely caught"** means: among the layers in our stack, this layer is the only
one that flagged the defect, so removing it sends the defect to the next-weakest
probe — or, if none catches it, to the undetected tail. Where a defect was visible
to two layers (A4 = parity ∧ round-trip on size; A15 = conformance ∧ fuzz) we count
it as *non-unique* for both and say so explicitly, so no escape-rate is inflated.

**Denominator.** The ablation is computed over the **22 attributable
output-correctness defects**: Corpus A rows A1, A2, A4, A6–A24, *excluding*
(a) **A3** (a perf regression, not wrong output), (b) **A5** (a mischaracterized
performance label, not a code bug), and (c) **A25** (the 7th wave bug whose
subcommand the source prose does not pin — kept in the corpus for honesty,
excluded here to avoid crediting a layer for an unidentified defect). A3 and A5
are reported separately under *non-correctness probes* below.

**Survivorship caveat (binding).** This corpus contains only *detected* bugs.
Defects that escaped **every** layer are by definition absent, so the ablation
measures *relative layer value among caught bugs* and is silent on the undetected
tail. The only probes into that tail are differential fuzzing
([`difffuzz_run.txt`](difffuzz_run.txt)) and the real-data byte-exact parity
battery on whole GIAB inputs (`realparity`/`realbench`, upstream as the oracle —
claim C2). Read every "escape rate if removed" below as *escape-among-caught*,
not *escape-absolute*.

---

## The ablation table

Denominator = 22 attributable output-correctness defects. "Uniquely caught" =
defects no other layer in the stack flagged; "escape rate if removed" =
uniquely-caught / 22 (the fraction of caught defects that would have slipped to the
next-weakest probe or the undetected tail had this layer been deleted).

| layer | defects caught | uniquely caught | escape rate if removed | unique rows |
|---|---|---|---|---|
| 1 unit / compile | 1 | 1 | 1 / 22 = 4.5% | A8 |
| 2 differential parity | 15 | 13 | 13 / 22 = 59.1% | A2, A6, A7, A9, A10, A11, A12, A19, A20, A21, A22, A23, A24 |
| 3 round-trip / metamorphic | 2 | 1 | 1 / 22 = 4.5% | A1 |
| 4 differential fuzzing | 4 | 3 | 3 / 22 = 13.6% | A16, A17, A18 |
| 5 human review | 0 (code bugs) | 0 | 0 / 22 = 0% | — (A5 is non-correctness) |
| conformance corpus (separate oracle) | 3 | 2 | 2 / 22 = 9.1% | A13, A14 |

Rows that are **shared** between layers, and therefore *not* unique to either
(so they appear in no layer's escape-rate numerator):

- **A4** (CRAM literal-bases) — visible to both layer 2 (size diff vs upstream)
  and layer 3 (round-trip). Counted in "defects caught" for both 2 and 3, unique to
  neither.
- **A15** (empty/headerless file accepted) — visible to both conformance and
  differential fuzz. Counted in "defects caught" for both, unique to neither.

Sum check: unique catches = 1 (unit) + 13 (parity) + 1 (round-trip) + 3 (fuzz) +
0 (human) + 2 (conformance) = **20 uniquely-attributable**, plus **2 shared**
(A4, A15) = **22**. ✓

### Non-correctness probes (reported separately, not in the 22-defect denominator)

| probe | defect | what it caught |
|---|---|---|
| profiling / bench | A3 | bzip2 brute-force perf regression (9.5× slower) — uniquely caught by measuring. |
| profiling / bench + human review | A5 | a *wrong explanation* ("libm-bound"; "decode-bound"), corrected by the profile and flagged by the human. Not wrong output. |

---

## Reading of the result (what the manuscript should claim)

- **Differential parity is load-bearing.** It uniquely accounts for ~59% of caught
  output-correctness defects; deleting it would have let 13 of 22 escape. This is
  the central methods claim, and it survives the conservative bookkeeping
  (shared rows credited to no one).
- **Round-trip/metamorphic earns its place on quality, not quantity.** It uniquely
  caught **one** defect (A1) — but that defect is the single most important in the
  corpus: a silent CRAM corruption (70k/300k reads → `NNNN`) that **single-op
  parity passed** because our own decoder masked it. The argument "layered
  validation > single-operation testing" rests on this one prospective-checkable
  row, and the manuscript should present it as a *case*, not a rate.
- **Differential fuzzing and the conformance corpus each pull their weight on the
  input tail** (3 and 2 unique catches) — exactly the untested-input and
  spec-compliance regions the seeded parity sweep does not reach. The live fuzz
  run ([`difffuzz_run.txt`](difffuzz_run.txt)) shows the harness still finding
  stdout/stderr/exit divergences on `bcftools view` and `bcftools query`, evidence
  the layer keeps finding real divergences rather than being a one-off.
- **Unit tests are necessary but not sufficient.** They uniquely caught **one**
  defect (A8, a silent SVG-write error found while pushing coverage). Every
  *silent-divergence* row (A1, A2, A11, A15) was caught by a layer **above** unit.
  This is the honest counterweight to "we have high unit coverage."
- **Human review caught zero code bugs in this corpus** — its one catch (A5) was a
  wrong performance *explanation*. We do **not** claim human review is unimportant
  (the scope-steering and the very act of building these harnesses were human);
  we claim it did not, in this dataset, uniquely catch an output-correctness
  defect. Report it as such.

---

## Provenance: MEASURED vs RECONSTRUCTED (every number)

This is the section reviewers will scrutinize. We grade each input to the table.

### MEASURED (a cited artifact establishes the layer attribution)

- **A11, A12** — caught by the edge-case battery; the run is in
  [`conformance_run.txt`](conformance_run.txt) (`TestBCFToolsNormReindex`,
  `TestIndexByteIdentity` pass after fix). Layer = edge-case/parity. **MEASURED at
  layer granularity.**
- **A13, A14, A15** — caught by the htslib conformance corpus; run in
  [`conformance_run.txt`](conformance_run.txt) (`TestHtslibSAM_*`,
  `TestHtslibEmptyFile`, `TestHtslibLongRefs`). **MEASURED.**
- **A16, A17, A18, A15** — caught by differential fuzz; run in
  [`difffuzz_run.txt`](difffuzz_run.txt) (live `bcftools view`/`query` stdout,
  stderr, and exit-code divergences). **MEASURED** that the fuzzer finds
  divergences in these targets; the *specific* A16/A17/A18 minimizations are named
  in the corpus from the harness's minimized reproducers.
- **A1** — the fixing commit `5f811f0` and `1068c77` state that single-op parity
  passed and only the cross-decode round-trip exposed the defect. Layer attribution
  is **MEASURED** in the weak sense (a contemporaneous commit names the layer); it
  is the one row explicitly called out as prospectively re-checkable.
- **A3, A5** — commits `f970994`, `19b0430`/`c755e45` name profiling as the
  detector. **MEASURED** (non-correctness probes).

### RECONSTRUCTED (layer inferred from git/PR prose, not a prospective re-run)

- **A2, A6, A7, A8, A9, A10** — layer inferred from the fixing commit's context
  (medium full-validation for A6; coverage push for A8; etc.). The `caught_by` is a
  plausible single layer but was **not** established by re-running the bug-commit
  through each layer in isolation. **RECONSTRUCTED.**
- **A19, A20, A21, A22, A23, A24** — the **count** (7 port bugs) is MEASURED from
  `PROJECT_STATUS.md`'s 2026-06-14 wave note, and the layer (differential parity) is
  structurally certain (these are `t.Skip`→`t.Fatalf` parity-assertion
  conversions). But the **per-row split** of the count into individual subcommands is
  **RECONSTRUCTED** from that paragraph's prose, and the **7th** bug (A25) is
  unidentified and excluded.
- **The "uniquely caught" judgments** for every RECONSTRUCTED row are themselves
  **RECONSTRUCTED**: we argue uniqueness structurally (a parity-only fixture has no
  unit test asserting it; a fuzz-only input has no seeded fixture) rather than by
  the prospective in-isolation re-run the corpus header prescribes.

### NOT YET DONE (the honesty gap to close before submission)

- **Prospective ablation** — checkout each Corpus-A bug-commit's parent, run
  {unit, parity, round-trip, fuzz} in isolation, *measure* which catch it. Only A1
  has been argued prospectively; the rest are reconstructed. This converts the table
  from assertion to measurement and is the corpus's #3 pre-submission TODO.
- **Independent re-labeling** — a second rater should re-derive `caught_by` and
  report Cohen's κ; the current labels are single-rater (confirmation-bias risk).
- **Report as a case series, not with CIs** — n ≈ two-dozen; do not dress the
  escape rates in confidence intervals. The percentages above are *descriptive
  fractions of a small labeled set*, not estimates of a population rate.

---

## Cross-references

- Labeled dataset + per-layer summary: [`../bug_corpus.md`](../bug_corpus.md)
  (Corpus A, "Per-layer summary").
- Claim definition: [`../01_CLAIMS_AND_EXPERIMENTS.md`](../01_CLAIMS_AND_EXPERIMENTS.md)
  §C7.
- Method narrative + the A1 flagship: [`../02_WHAT_WORKED_WHAT_FAILED.md`](../02_WHAT_WORKED_WHAT_FAILED.md)
  §A.3, §C.2.
- Artifacts: [`difffuzz_run.txt`](difffuzz_run.txt),
  [`conformance_run.txt`](conformance_run.txt).
