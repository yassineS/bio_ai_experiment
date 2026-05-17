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

#### vcftools-site-pi

- **vcftools `--site-pi` formula** — upstream computes a per-genotype
  pairwise-distance quantity rather than the textbook
  `(n² − Σ cₐ²) / (n(n-1))`. Our Go port uses the textbook formula (#24
  flagged this as a likely bug; we then fixed it on our side). Open
  question: is upstream's behaviour a real bug, or a deliberate
  per-genotype variant of nucleotide diversity? Need to read the
  vcftools paper and compare to PopGenome / scikit-allel. If it's a
  bug, this entry is **fix-on-port (done)**; if intentional, we should
  add a `--site-pi-vcftools-compat` flag.

  Parity test status: `TestParity_SitePi` is `t.Skip("known deviation,
  see docs/UPSTREAM_BUGS.md#vcftools-site-pi")`. A separate
  `TestParity_SitePi_TextbookFormula` spot-checks three hand-computed
  values against our textbook implementation.

#### mosdepth-overlap-pair-detection

- **mosdepth overlap-pair detection** — upstream subtracts one copy of
  depth where the two ends of a mate-paired fragment overlap on the
  reference, so a 80bp read pair with a 100bp insert contributes depth
  1 to the overlapped region (not depth 2). Our v1 engine doesn't
  implement this pairing: every aligned base of every read contributes
  to depth. The net effect is that our default-mode output matches
  upstream's `--fast-mode` output rather than upstream's default
  output. This is **NOT an upstream bug** — it's a feature gap in our
  port — but it lives here because every affected parity test cites
  this entry from a `t.Skip("known deviation, see
  docs/UPSTREAM_BUGS.md#mosdepth-overlap-pair-detection")`.

  Disposition: **track-only** until we add a read-name-keyed pairing
  pass. Five mosdepth parity tests reference this anchor.

- **bedtools `groupby` empty-group handling** (when we get to porting
  it) — Aaron Quinlan has acknowledged upstream emits a blank line on
  empty groups in some flag combinations; we should not.

- **htslib BGZF EOF block tolerance** — htslib historically warned but
  continued on missing EOF blocks for years before tightening to an
  error in recent versions. Our `bgzip.NewReader` errors out by
  default; should we add a `--ignore-truncation` flag matching newer
  htslib?

### Fix-on-port (resolved)

#### vcftools `.ifreqburden` INDV label-index bug

Upstream `variant_file_output.cpp:621` emits
`out << meta_data.indv[indv_count] << ...` for the leading INDV column
of the `.ifreqburden` output, where `indv_count` is the local kept-row
index. The correct sibling write at `output_indv_burden:494` uses
`meta_data.indv[ui]` (the original-index). When `--remove-indv` (or
any other sample filter) drops a non-trailing sample from
`[S1,S2,S3,S4]`, upstream emits `S1 S2 S3` next to the S1/S3/S4 burden
values — wrong-labelled rows that downstream analyses will silently
misread.

**Fixed in port** (PR #138, wave 14). We use the kept-sample list for
both the row order AND the INDV label. The sibling `.iburden` output
was already correct upstream and is unaffected. Pinned by
`TestIndvFreqBurden_RemoveIndv_FixesLabelBug` in
`tools/vcftools/pkg/vcftools/burden_test.go`.



PR #55 (bedtools parity) fixed 7 discrepancies in our Go code; see
`tools/PARITY_VALIDATION.md` for the bedtools list.

The column-ops + discrepancies wave (this PR) fixed the last two
PR-#55-era bedtools discrepancies on our Go side:

- **bedjaccard did not pre-merge A or B before computing
  intersection / union.** Upstream `bedtools jaccard` sets
  `setUseMergedIntervals(true)` on its context
  (`reference_code/bedtools/src/utils/Contexts/ContextJaccard.cpp`),
  which makes its FileRecordMergeMgr stream-merge both inputs before
  the sweep. Our bedjaccard now wraps each input reader in a
  `mergingReader` that does the same; the merge runs per-strand
  under `-s`. Newly passing parity cases: `jaccard.t02 / t03 / t05 /
  t06 / t10 / t11`.
- **bedmerge's `-s` did not match upstream's per-strand merge.**
  Upstream's `FileRecordMergeMgr`
  (`reference_code/bedtools/src/utils/FileRecordTools/FileRecordMergeMgr.cpp`
  lines 47-58 + 96-129) drops UNKNOWN-strand (`.` / empty) records
  under `-s` and merges `+` / `-` independently. Our `mergeIntervals`
  and `mergeWithColumnOps` now split into per-strand buckets under
  `StrandSpec` (dropping `.` records) and re-combine the two streams
  in `(chrom, start, end)` order. Newly passing: `merge.t15`.

The samtools parity audit (PR #75) surfaced three bugs **in our Go
port** (not upstream):

- **samtools sort `-n` / `-N` CLI mapping inverted.** Upstream's `-n` is
  natural numeric name sort (the default for name-sort) and `-N` is plain
  lexicographic; we had it reversed. Fixed by swapping the CLI binding in
  `tools/samtools/cmd/samtools/main.go` (the library API
  `SortByName` / `SortByNameNatural` was already correct). Surfaced by
  comparing record ordering against
  `reference_code/samtools/test/sort/name.sort.expected.sam`.

- **samtools sort missing `SS:queryname:*` sub-sort tag.** Upstream
  writes `SS:queryname:natural` or `SS:queryname:lexicographical` on the
  `@HD` line so downstream tooling can recognise the sub-form. Our
  `stampSortOrder` in `tools/samtools/pkg/samtools/sort.go` only stamped
  the `SO` field; added the `SS` field too.

- **samtools fastq pair-suffix not auto-dropped in `-1/-2` mode.**
  Upstream `bam_fastq.c` sets `has12 = false` whenever `-1` or `-2` is
  given because the separate file names already disambiguate mate
  identity; we were unconditionally appending `/1`/`/2`. Fixed in
  `Fastq` in `tools/samtools/pkg/samtools/fastq.go`; `-N` still forces
  the suffix.

The sickle/skewer parity audit (PR #73) added one build-side fix-on-port:

- **skewer `src/matrix.h` is not `const`-correct** under modern
  libstdc++.
  `ElementComparator::operator()` lacks the trailing `const` that gcc 13's
  `<bits/stl_tree.h>` `static_assert` now requires for `std::set<>`
  comparators. We patch the submodule working tree
  (`bool operator()(...) const`) before building the upstream binary used
  to generate parity fixtures. The patch is not pushed back to upstream
  (the repo is dormant) and is not part of our Go code.

The bcftools parity audit (PR #74) fixed 9 discrepancies, all on our
side rather than upstream's. They're recorded in
`tools/PARITY_VALIDATION.md#discrepancies-found-and-fixed-in-this-pr`.

The vcftools + mosdepth validated-parity audit (this PR) found 9 small
discrepancies in our Go code (not upstream), all fixed inline:

- **vcftools `.frq` / `.frq.count` header + row format** — upstream
  emits a single literal `{ALLELE:FREQ}` / `{ALLELE:COUNT}` column
  header with one tab-separated `allele:value` entry per allele in
  the data rows; we were emitting `{REF:FREQ}\t{ALT:FREQ}` headers
  and `{A:0.833333}`-wrapped data cells. Fixed.
- **vcftools `.singletons` column count** — upstream emits five
  columns (`CHROM`, `POS`, `SINGLETON/DOUBLETON`, `ALLELE`, `INDV`);
  we were emitting three. Fixed; we now resolve which individual
  carries the rare allele and tag singletons (S) vs private doubletons
  (D).
- **vcftools `.hwe` header** uses `CHR` (not `CHROM`) and includes
  `P_HET_DEFICIT` + `P_HET_EXCESS` columns. Fixed header; directional
  P-values are placeholders pending Wigginton-Cutler-Abecasis impl.
- **vcftools `.lmiss` header** uses `CHR` (not `CHROM`). Fixed.
- **vcftools `.ldepth`** needs a `SUMSQ_DEPTH` column. Fixed header +
  data; value is a literal 0 placeholder pending the sum-of-squares
  accumulator.
- **vcftools `.TsTv.summary` layout** — upstream emits a `MODEL\tCOUNT`
  table with six per-substitution rows plus two roll-up rows; we were
  emitting a `Ts\tTv\tTs/Tv` two-line table. Fixed.
- **vcftools `.012` row prefix** — upstream prefixes each row with the
  0-based sample index; we were prefixing with the sample name. Fixed.
  Also confirmed `.012` is biallelic-only (matches upstream's one-off
  warning + skip on multi-allelic).
- **vcftools `--max-missing 1.0` was a no-op** because the filter was
  guarded by `< 1`. Fixed to `> 0` so `--max-missing 1.0` (require all
  non-missing) actually filters.
- **vcftools `.windowed.pi` header** needs `N_MONOMORPHIC` column.
  Fixed header + data; value is a 0 placeholder.

### Track-only (parity skipped, fix later)

#### vcftools `--hapcount` prev_bin_idx shift on bin change

**Severity:** numerical (mis-attributed per-bin counts and merged
multiplicity histograms in `.hapcount`).

**Status:** PRESERVED in the Go port for byte-for-byte parity. Pinned
by `TestParity_Hapcount` in `tools/vcftools/pkg/vcftools/hapcount_test.go`.

Upstream `variant_file_output.cpp:1314-1315` unconditionally assigns
`prev_bin_idx = bin_idx; bin_idx = ui;` at the start of every per-site
search loop's successful BED-bin match. The flush-trigger predicate at
line 1322 is then `if ((found == false) || (prev_bin_idx != bin_idx))`
— so after a within-chromosome bin-transition flush has fired, the
next per-site iteration leaves `prev_bin_idx` pointing at the OLD bin
index even though the data has moved on. The next time a flush fires
(at the next bin change or at end-of-chromosome), `SNP_count[prev_bin_idx]`
and `haplotype_count[prev_bin_idx]` get OVERWRITTEN with the new bin's
values, AND `haplotype_frequencies[prev_bin_idx]` gets the new bin's
histogram ADDED to the old bin's histogram. The N_GROUPS column
reflects the union and the N_SNP / N_UNIQ_HAPS columns reflect the
latest write.

We mirror this verbatim because the project's parity bar is
byte-for-byte vs upstream output. Disposition: **track-only** until
either upstream is patched (unlikely; project is dormant) or the port
gains an opt-in `--hapcount-correct-binning` flag.

#### vcftools `--hapcount` end-of-stream read-after-free

**Severity:** crash-on-input / truncated output.

**Status:** PRESERVED in the Go port for byte-for-byte parity. Pinned
indirectly by `TestParity_Hapcount` (fixture ends with a sentinel
chromosome to force the last real-data chromosome's emission via the
chrom-transition path).

Upstream `variant_file_output.cpp:1370-1400` reads `e->include_indv[ui]`
inside the final-flush block AFTER `delete e;` at line 1370. Observed
behaviour on the upstream binary in this repo:

- When `have_data == true` at EOF, the read-after-free corrupts the
  final write path and the last chromosome's bins are SILENTLY DROPPED
  from the output file (the process exits 0 but the buffered final
  rows are never flushed).
- When `have_data == false` at EOF, the freed-pointer access happens
  to skip every iteration of the inner per-individual loop (likely
  because the freed memory reads as "false" for the
  `include_indv[ui]` check), so the final chromosome's bins are
  emitted with all-zero values.

The Go port replicates both branches in `hapcount.close()`. Test
fixtures should end with a sentinel chromosome (one row whose only job
is to trigger the chrom transition for the last real-data chromosome)
to avoid relying on the EOF code path.

- **skewer adapter matcher: Hamming vs Smith-Waterman.** Upstream's
  matcher rejects a 1-mismatch alignment whose mismatch is in the last 4
  bases of the adapter, even when the overall error rate is within `-r`.
  Our Go port's `improvedFindAdapter` is plain Hamming distance and
  accepts the alignment, over-trimming one base. `t.Skip` set on the
  case05 parity test with a pointer to the [skewer
  section](#skewer-case05). This is a Go-port limitation, not an
  upstream bug.

<a id="bcf-int64"></a>

- **`htslib` BCF type-4 (int64) emission for fields that fit in int32**
  — htslib 1.13+ optionally writes some FORMAT counters as the BCF
  `int64` typed descriptor even when the values are well within int32
  range. The BCF 2.2 spec allows it; in practice it just enlarges
  records and forces every consumer to add a 64-bit path. Not a bug per
  se but worth raising with the htslib maintainers about whether the
  encoder should clamp to int32 when it can. Tracked here so the next
  port can decide whether to round-trip 64-bit verbatim or downcast on
  read (we currently downcast — see `pkg/bioformats/bcf/typed.go`).

<a id="bcf-fmt-keys-missing"></a>

- **Our BCF reader drops per-record FORMAT keys on htslib-produced
  input** — after the int64 and IDX-strip fixes in this PR, the header
  parses fine but per-record `FmtKeys` come back as
  `[<resolved>, -1, -1, ...]`: only the first key resolves correctly,
  the rest decode as `MissingInt32`. This is almost certainly a bug in
  our `decodeIndiv` (probably mis-counting the dictionary index width
  for n_sample > 0 when the key entry is itself a typed-int vector,
  but we haven't bottomed it out). Workaround: use VCF / VCF.gz input.

<a id="bcf-info-order"></a>

- **Our BCF writer does not preserve `InfoOrder` on encode** — the
  reader-side fix in this PR populates `Variant.InfoOrder` so the VCF
  writer can preserve key order, but the BCF writer still iterates the
  map directly. Consequently a VCF→BCF→VCF round-trip shuffles INFO
  keys. The fix is to teach `bcf.NewWriterFromVCFHeader` /
  `bcf.Writer.Write` to consult `InfoOrder` in addition to the existing
  dict-order pass.

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

## seqtk

The seqtk parity audit (PR for `prinseq-seqtk-parity-validation`)
fixed four discrepancies on **our side** (not upstream); they are
listed under
[tools/PARITY_VALIDATION.md > seqtk parity validation](../tools/PARITY_VALIDATION.md).
The audit also surfaced three behavioural divergences where byte
parity is impractical without porting upstream-specific RNG / algorithm
machinery; we track them here so the parity tests can point at a
stable anchor.

### seqtk-sample-rng <a id="seqtk-sample-rng"></a>

- **Symptom.** Upstream `seqtk sample` uses a seeded reservoir
  sampler (`drand48()`-based, default seed 11) and produces a
  deterministic, fraction-correct subset for any input. Our Go port
  implements a "deterministic every-Nth-record" sampler that does
  not match upstream's selection regardless of seed.

- **Root cause.** Upstream's algorithm is the streaming reservoir
  pass over `n` records with `keep = (drand48() < frac)` per record
  in one-pass mode, and a two-pass random subset in `-2` mode. Our
  port short-circuits to `written/count < fraction` per record.

- **Disposition.** Track-only — Go-port limitation, not an upstream
  bug. The parity test for `sample` is split into a structural
  invariants pass (`TestParity_Seqtk_Sample_StructuralInvariants`,
  passing) and a byte-parity case
  (`TestParity_Seqtk_Sample_UpstreamByteParity`, skipped).
  Fixing this is straightforward (port `drand48` against the same
  seed); deferred because no caller currently relies on the
  byte-for-byte output.

### seqtk-randbase-rng <a id="seqtk-randbase-rng"></a>

- **Symptom.** Upstream `seqtk randbase` uses `drand48()` (with
  implicit seed 0) and is therefore deterministic across runs but
  not seed-controllable. Our Go port uses `math/rand` with a
  caller-supplied seed.

- **Root cause.** Different RNGs.

- **Disposition.** Track-only — Go-port limitation, not an upstream
  bug. The structural invariants
  (`TestParity_Seqtk_Randbase_StructuralInvariants`) verify that
  upstream's rules (only 2-base IUPAC codes are randomised; 3-base
  and 4-base codes pass through; case is preserved) are honoured
  on our side. Byte parity is skipped.

### seqtk-trimfq-algorithm <a id="seqtk-trimfq-algorithm"></a>

- **Symptom.** Upstream `seqtk trimfq` runs a modified Mott
  algorithm with an error-rate threshold (default `-q 0.05`) and a
  `-l 30` minimum-length floor. Our port's `TrimQuality` does a
  simple Phred-quality threshold trim. The two algorithms produce
  different cuts on every non-trivial input.

- **Disposition.** Track-only — feature-gap on our side, not an
  upstream bug. Tracked in `docs/PARITY_ROADMAP.md#seqtk` under the
  `trimfq` option-tail gaps; will be closed once we re-implement
  the Mott trim. The parity test
  (`TestParity_Seqtk_Trimfq_UpstreamByteParity`) is skipped with a
  pointer here.

## prinseq

The prinseq parity audit (PR for `prinseq-seqtk-parity-validation`)
fixed three discrepancies on **our side** (not upstream); they are
listed under
[tools/PARITY_VALIDATION.md > prinseq parity validation](../tools/PARITY_VALIDATION.md).
No upstream bugs surfaced — PRINSEQ-lite's documented behaviour
agreed with the corpus we tested for every option we exercised
(see the prinseq table for the 18 cases).

### `--graph_data` JSON key order is non-deterministic (`prinseq-lite.pl:2050-2287`)

**Severity:** behavioural / reproducibility (not a numerical bug).

**Status:** worked around in our Go port (lexicographic key order).

Upstream's `.gd` emitter walks `keys %hash` and `each %hash` directly
when writing the JSON-shaped string. Since Perl 5.18 (released
mid-2013) every interpreter start re-seeds the hash function, so
running the same `prinseq-lite.pl ... -graph_data <out>` twice produces
two `.gd` files with the same numerical content but different key
order. This is fine for downstream renderers that parse with the
`JSON` module, but it makes byte-for-byte regression tests impossible
and complicates diff-based review of stat output.

The Go port at `tools/prinseq/pkg/prinseq/graphdata.go` sorts every
map key (string-comparison order on the printed form, matching how
upstream renders the JSON) before emission. The parity test
(`tools/prinseq/pkg/prinseq/graphdata_test.go:TestGraphDataParityExample1`)
compares the upstream-shipped `example1.gd` and our emit by parsing
both as JSON and recursively diffing the structures with a 1e-3
absolute tolerance — which removes the key-order issue from the
comparison while keeping the numerical content authoritative.
`TestGraphDataDeterminism` separately asserts that two consecutive
runs of our emitter produce byte-identical output, the property
upstream lacks.

We classify this as an upstream behavioural defect rather than an
intentional design choice because:

1. The emission code predates Perl 5.18's randomised hashing — the
   author wrote a portable JSON-string assembler that happened to
   inherit hash randomness only after Perl was updated under their
   feet.
2. There is no downstream component that benefits from random key
   order; the companion `prinseq-graphs.pl` parses the file and
   reorders fields itself.

Per the disposition policy at the top of this file, fixing on our
side and documenting the deviation is the right call.
