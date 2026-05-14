# Upstream bug tracker

When we hit a real bug in an upstream tool — not a feature gap, not a design
choice we disagree with — we record it here. The Go port shouldn't carry the
bug over; we either fix it on the way (the better option when the fix is
small and obvious) or accept the discrepancy and add a `parity_test.go` skip
with a pointer to this file.

## What counts as a bug

- Off-by-one in a clearly-bounded operation (e.g. interval arithmetic
  that produces wrong output for an exactly-on-boundary input).
- Integer overflow, undefined behaviour, signed/unsigned confusion that
  produces wrong output (not just a crash).
- Specification non-compliance: upstream emits output the format spec
  doesn't allow (e.g. invalid VCF line, malformed BGZF block).
- Documented contract violated: the man page or `--help` promises X and
  the binary does Y.
- Crash on input that the tool is meant to accept (e.g. empty file).

## What is **not** a bug

- Different default values than we'd pick.
- Awkward but documented CLI semantics (`bcftools view -h` = "header-only").
- Slow algorithms.
- Missing features.
- Output that's surprising but consistent with the spec or a long-standing
  convention.

If unsure, file it here under "to investigate" rather than silently
diverge — it's easier to fix a wrong call than to chase a silent
discrepancy later.

## Disposition policy

For each entry:

- **Fix-on-port** — small, obvious correction; do it as part of the
  porting work, link the PR.
- **Track-only** — the right behaviour is non-trivial to determine, or
  the fix is large; leave a `t.Skip("upstream bug, see UPSTREAM_BUGS.md#<anchor>")`
  in the parity test until we have time to do it properly.
- **Resolved-as-feature** — turned out to be a documented behaviour; demote
  the entry to the "non-bugs we considered" section at the bottom.

Order of work is judgement-based: fix small obvious ones as we port the
relevant tool; track-only the ones that need real thought.

## Open entries

_Empty as of 2026-05-14. Will be populated by the Phase 2 validated-parity
audits against upstream test suites (samtools, bcftools, mosdepth, vcftools)
and by the bedtools long-tail closure. Existing notes from prior work that
warrant a closer look:_

### To investigate

- **vcftools `--site-pi` formula** — upstream computes a per-genotype
  pairwise-distance quantity rather than the textbook
  `(n² − Σ cₐ²) / (n(n-1))`. Our Go port uses the textbook formula (#24
  flagged this as a likely bug; we then fixed it on our side). Open
  question: is upstream's behaviour a real bug, or a deliberate
  per-genotype variant of nucleotide diversity? Need to read the
  vcftools paper and compare to PopGenome / scikit-allel. If it's a
  bug, this entry is **fix-on-port (done)**; if intentional, we should
  add a `--site-pi-vcftools-compat` flag.

- **bedtools `groupby` empty-group handling** (when we get to porting
  it) — Aaron Quinlan has acknowledged upstream emits a blank line on
  empty groups in some flag combinations; we should not.

- **htslib BGZF EOF block tolerance** — htslib historically warned but
  continued on missing EOF blocks for years before tightening to an
  error in recent versions. Our `bgzip.NewReader` errors out by
  default; should we add a `--ignore-truncation` flag matching newer
  htslib?

### Fix-on-port (resolved)

_None yet. PR #55 (bedtools parity) fixed 7 discrepancies but those were
in our Go code, not upstream._

The sickle/skewer parity audit added one build-side fix-on-port:

- **skewer `src/matrix.h` is not `const`-correct** under modern
  libstdc++.
  `ElementComparator::operator()` lacks the trailing `const` that gcc 13's
  `<bits/stl_tree.h>` `static_assert` now requires for `std::set<>`
  comparators. We patch the submodule working tree
  (`bool operator()(...) const`) before building the upstream binary used
  to generate parity fixtures. The patch is not pushed back to upstream
  (the repo is dormant) and is not part of our Go code. See the [skewer
  section below](#skewer) for the diff.

### Track-only (parity skipped, fix later)

- **skewer adapter matcher: Hamming vs Smith-Waterman.** Upstream's
  matcher rejects a 1-mismatch alignment whose mismatch is in the last 4
  bases of the adapter, even when the overall error rate is within `-r`.
  Our Go port's `improvedFindAdapter` is plain Hamming distance and
  accepts the alignment, over-trimming one base. `t.Skip` set on the
  case05 parity test with a pointer to the [skewer
  section](#skewer-case05). This is a Go-port limitation, not an
  upstream bug.

### Non-bugs we considered (closed)

The sickle audit found three behaviours that looked like upstream bugs
but turned out to be documented features:

- **`int(0.1 * read_length)` window sizing.** Looks like a quirky
  default; turns out to be the only window sizing rule upstream supports
  (sickle has no `-w` flag) and is described in the project README.

- **Discarding reads when no window reaches threshold.** Looks like an
  off-by-one in the discard path. Reading `src/sliding.c` confirms it
  is intentional: `if (found_five_prime == 0 && !no_fiveprime)` flags
  the read for discard via `three_prime_cut = -1`. Documented in
  upstream's source comments.

- **Phred-offset per encoding.** Both 33 (sanger) and 64 (illumina /
  solexa) decoded against their respective offsets via the
  `quality_constants` table in `src/sickle.h`. Not a bug — our Go port
  was just decoding all qualities with offset 33 regardless of `-t`. Fixed
  inline in this PR (see `tools/PARITY_VALIDATION.md > sickle`).

## skewer <a id="skewer"></a>

### skewer compile failure on modern libstdc++

- **Symptom.** Building `reference_code/skewer` with gcc 13 + libstdc++ 13
  fails with `static_assert failed: 'comparison object must be invocable
  as const'` in `<bits/stl_tree.h>`, traced back to `src/matrix.h`'s
  `class ElementComparator` whose `bool operator()(...)` is not
  `const`-qualified.

- **Root cause.** Modern libstdc++ tightened `std::set<>` to require its
  `_Compare` template parameter be invocable as `const`. Upstream's
  declaration is from the older relaxed era (last upstream commit 2017).
  Upstream is dormant.

- **Disposition.** Fix-on-port (build-side only). We apply a minimal
  patch to the submodule working tree before building the parity binary:

```diff
--- a/src/matrix.h
+++ b/src/matrix.h
@@ -52,3 +52,3 @@
-    bool operator()(const ELEMENT &elem1, const ELEMENT &elem2){
+    bool operator()(const ELEMENT &elem1, const ELEMENT &elem2) const {
         return elem1.idx.pos < elem2.idx.pos;
     }
```

  The patch is not committed back to the submodule (the submodule pointer
  stays at `978e8e4`). The Go port doesn't carry the underlying bug
  because the matrix-mode code path is not yet implemented in Go.

### skewer case05 SE error-tolerance matcher <a id="skewer-case05"></a>

- **Symptom.** Upstream rejects a 1-mismatch adapter alignment when the
  mismatch falls in the tail 4 bp of the adapter, even when
  `mismatches/len(adapter) <= -r`. Our Go port's `improvedFindAdapter`
  accepts it and over-trims one base.

- **Root cause.** Upstream uses a Smith-Waterman alignment with a
  position-dependent tail penalty; we use plain Hamming distance.

- **Disposition.** Track-only — Go-port limitation, not an upstream bug.
  `tools/skewer/pkg/skewer/parity_test.go > case05` has `t.Skip` with a
  pointer to
  [tools/PARITY_VALIDATION.md > "skewer"](../tools/PARITY_VALIDATION.md#skewer).

## sickle

No upstream bugs surfaced in the sickle audit. The behaviours that
initially looked suspicious are documented in
[the "Non-bugs we considered" section](#non-bugs-we-considered-closed).
