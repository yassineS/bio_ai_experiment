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

- **vcftools `--site-pi` formula** — **resolved-as-feature**. On
  closer reading of `variant_file_output.cpp:3870-3936`, upstream's
  per-site Pi computes `mismatches / (total_alleles *
  (total_alleles - 1))` where
  `mismatches = Σ_a c_a × (total - c_a) = total² - Σ c_a²`. This is
  algebraically identical to the textbook
  `(n² − Σ cₐ²) / (n(n-1))` formula. The earlier "deviation" was
  a byte-format divergence: our port emitted `%.6f` ("0.600000")
  while upstream uses the C++ ostream default (six significant
  digits, shortest representation, "0.6"). After the
  deferred-test-closure wave the port emits with `formatCppDouble`
  (`strconv.FormatFloat(x, 'g', 6, 64)`) and gates Pi on
  `is_diploid()` to skip haploid sites the way upstream does. Byte
  parity holds on `sample.vcf` (`TestParity_SitePi`); the textbook
  spot-check `TestParity_SitePi_TextbookFormula` was updated to
  assert the C++ default literal form ("0.6" rather than "0.600000").

#### mosdepth-overlap-pair-detection

- **mosdepth overlap-pair detection** — upstream subtracts one copy of
  depth where the two ends of a mate-paired fragment overlap on the
  reference, so a 80bp read pair with a 100bp insert contributes depth
  1 to the overlapped region (not depth 2).

  **Resolved.** The default-mode engine now buffers a QName -> reference
  interval map per chromosome and emits a sign-flipped event pair over
  the per-fragment overlap interval the second time a QName appears.
  `--fast-mode` keeps the no-pairing fast path (matching upstream).
  The six previously-skipped parity tests
  (`TestParity_OverlapM_DefaultPerBase`, `_OverlapM_SummaryMT`,
  `_ThresholdByBED`, `_TrackHeader`, `_MAPQFilter`, `_FlagExclude`)
  now assert upstream byte-for-byte values.

- **bedtools `groupby` empty-group handling** (when we get to porting
  it) — Aaron Quinlan has acknowledged upstream emits a blank line on
  empty groups in some flag combinations; we should not.

- **htslib BGZF EOF block tolerance** — htslib historically warned but
  continued on missing EOF blocks for years before tightening to an
  error in recent versions. Our `bgzip.NewReader` errors out by
  default; should we add a `--ignore-truncation` flag matching newer
  htslib?

### Fix-on-port (resolved)

#### BCF writer correctness fixes (wave 21)

While wiring `--recode-bcf` (the wave-21 vcftools flag) we discovered
three latent bugs in `pkg/htsgo/bcf`'s writer that produced
self-consistent output (our own reader could roundtrip it) but
diverged from the BCF v2.2 spec and broke htslib interop:

1. **Missing-ID encoding.** The writer emitted missing string IDs as
   the descriptor byte `0x00` (type 0, "missing scalar"). Per the BCF
   spec a missing string is a zero-length typed-char vector,
   descriptor byte `0x07` (type 7, length 0). htslib's reader rejects
   the type-0 form with `Expected type 7 for string. Found type 0.`
   and then mis-aligns on subsequent FORMAT-block parsing.

2. **Unified dictionary numbering / missing IDX annotations.** BCF
   records refer to FILTER, INFO, and FORMAT tags by a single shared
   integer index, not a slice-local one. htslib annotates each
   `##INFO/##FILTER/##FORMAT` text-header line with a `,IDX=N`
   suffix so consumers can build the dictionary without re-deriving
   the numbering. The wave-21 writer now: (a) parses `,IDX=N` if
   present and otherwise auto-assigns IDX in declaration order across
   all three groups, (b) stores the IDX on each `DictEntry`, (c) maps
   tag-name → IDX (not name → slice position) for wire emission, and
   (d) round-trips `,IDX=N` annotations on every output header.

3. **FORMAT descriptor size = total flat length.** The writer encoded
   FORMAT fields with `size = nSample × per-sample-dim` in the
   descriptor's high nibble (or overflow). The spec says `size` is
   the **per-sample dimension**; htslib reads it that way and then
   multiplies by `n_sample` from the record header to walk the
   payload. The wave-21 writer (and reader) now interpret `size`
   per-spec.

The reader was updated symmetrically — `DecodeFormatTyped(buf, off,
nSample)` reads `nSample × per-sample-dim` elements; `splitPerSample`
uses `tv.Length` as the per-sample dim directly.

Pinned by `TestRun_RecodeBCF_Roundtrip`,
`TestRun_RecodeBCF_HeaderHasIDXAnnotations`, the existing
`TestWriterPerSampleRoundTrip`, and an out-of-band interop check
against upstream `vcftools --bcf <ours.recode.bcf>` (decoded VCF
matches the source).

#### BCF header name-dedup across INFO/FILTER/FORMAT (wave 22)

The wave-21 BCF correctness work assigned each `##INFO/##FILTER/
##FORMAT` line a fresh monotonically-increasing unified IDX. htslib's
actual policy (vcf.c:1500-1530) is to **deduplicate by name across
all three groups**: if a `##FORMAT=<ID=DP>` line follows a
`##INFO=<ID=DP>` line, the FORMAT entry **reuses** the INFO entry's
IDX rather than getting a new value. Without IDX annotations on the
text-header lines (htslib's writer commonly omits them on simple
inputs), our auto-assignment diverged from upstream's wire numbering
and `FmtTag(idx)` returned nil for upstream-produced FORMAT/DP
records — silently dropping the field at decode time.

Reproduced via `vcftools --vcf foo.vcf --recode-bcf` (writes BCF
without `,IDX=` annotations) then our `vcftools --bcf foo.bcf
--recode`: the output had FORMAT field columns mis-aligned by one
slot.

**Fixed in port** (wave 22). `parseTextHeader` now maintains a
`nameIDX` map and on each `##INFO/##FILTER/##FORMAT` line: (a)
honours an explicit `,IDX=N` annotation if present, (b) reuses the
name's existing IDX if seen earlier in either group, (c) falls back
to a fresh auto-IDX otherwise. Pinned by
`TestParseTextHeader_NameDedupAcrossINFOAndFORMAT` and an end-to-end
upstream-BCF → port roundtrip check.

#### vcftools `--keep-INFO TAG` semantic divergence (port → upstream)

This is a port-side divergence (not an upstream bug): pre-wave-17 the
Go port mapped `--keep-INFO TAG` onto upstream's recode-column
selector semantic, when upstream actually defines it as a SITE FILTER.

Upstream `parameters.cpp:266`:

```cpp
else if (in_str == "--keep-INFO") { site_INFO_flags_to_keep.insert(get_arg(i+1)); i++; }
```

The filter routine `entry::filter_sites_by_INFO` in
`entry_filters.cpp:1033-1063` requires every named tag to be declared
`Type=Flag` in the header (LOG.error otherwise) and then DROPS the
site unless at least one of the named tags has value "1" in the INFO
column (OR semantics across multiple tags). It is invoked from
`entry::apply_filters` at `entry_filters.cpp:44`. The semantic is
"keep sites where any of these Flag tags is present", not "restrict
the INFO column in the recoded output".

The recode-column selector is a separate flag — `--recode-INFO TAG`
(`parameters.cpp:319`, `recode_INFO_to_keep`). Pre-wave-16 the port
exposed the recode-column-selector semantic only via the misnamed
`--keep-INFO` flag; wave 16 (PR #141) added `--recode-INFO` as a
synonym pointing at the same internal slice and documented the
residual divergence as a follow-up.

**Severity:** silently-wrong output. A user invoking
`vcftools --vcf X.vcf --keep-INFO FLAG_A --recode` against the port
got every site in their VCF emitted with the INFO column restricted
to FLAG_A; against upstream they got just the sites where FLAG_A is
present, with INFO stripped to ".". Both row count and INFO content
differed.

**Fixed in port** (wave 17 — this PR). `--keep-INFO TAG` is now wired
to a new `Params.KeepINFO` site-filter codepath that mirrors
`entry_filters.cpp:1033-1063`: it errors at runtime if any named tag
is not declared `Type=Flag`, then drops sites where none of the named
tags is present in INFO. `--recode-INFO TAG` is now the sole flag
driving the recode-column selector (new `Params.RecodeINFO` field).

The two flags are independent and may be combined (e.g.
`--keep-INFO FLAG_A --recode-INFO DP` filters sites by FLAG_A
presence then restricts the recoded INFO column to DP).

Pinned by `TestPassKeepINFOSite`, `TestLookupInfoMeta`,
`TestRun_KeepINFO_SiteFilter_Integration`,
`TestRun_KeepINFO_SiteFilter_OR`, and
`TestRun_KeepINFO_SiteFilter_NonFlagType` in
`tools/vcftools/pkg/vcftools/info_filters_test.go`, plus the
byte-for-byte `TestParity_KeepINFO_SingleFlag` and
`TestParity_KeepINFO_OR` in
`tools/vcftools/pkg/vcftools/parity_test.go` against goldens generated
by the upstream binary (built with the FORTIFY_SOURCE workaround,
documented in `tools/vcftools/testdata/parity/keep_info_flags.vcf`
header comment).

**Sibling `--remove-INFO TAG` divergence resolved in wave 18** (see
below).

#### vcftools `--remove-INFO TAG` semantic divergence (port → upstream)

The polarity-inverted sibling of the wave-17 `--keep-INFO` fix. Same
shape of divergence: pre-wave-18 the Go port wired `--remove-INFO TAG`
as a recode-column stripper (a port-only invention), when upstream
defines it as a SITE FILTER.

Upstream `parameters.cpp:328`:

```cpp
else if (in_str == "--remove-INFO") { site_INFO_flags_to_remove.insert(get_arg(i+1)); i++; }
```

The remove path lives in the same `filter_sites_by_INFO` routine that
houses keep — `entry_filters.cpp:1068-1086`. Like the keep path it
requires every named tag to be declared `Type=Flag` in the header and
errors via `LOG.error` otherwise. The semantic is "drop the site if
any of the named tags has `value == "1"` in the INFO column"
(OR-veto). Composition with `--keep-INFO` is keep-then-remove: keep
narrows first (sites without any keep flag are dropped), then remove
vetoes the survivors. Upstream short-circuits at line 1066 if keep
already dropped the site.

The recode-column stripping behaviour the port previously implemented
had no upstream equivalent: upstream's only recode-column control is
`--recode-INFO TAG` (the keep-set selector, parameters.cpp:319), and
that landed as a distinct port flag in wave 17.

**Severity:** silently-wrong output. A user invoking
`vcftools --vcf X.vcf --remove-INFO FLAG_A --recode` against the port
got every site in their VCF emitted with FLAG_A stripped from the
INFO column; against upstream they got just the sites where FLAG_A
is absent, with INFO stripped to "." (unless `--recode-INFO-all`
preserves it). Both row count and INFO content differed.

**Fixed in port** (wave 18 — this PR). `--remove-INFO TAG` is now
wired to a new `passRemoveINFOSite` codepath in
`tools/vcftools/pkg/vcftools/info_filters.go` that mirrors
`entry_filters.cpp:1068-1086`. Header validation is shared with the
keep path via `validateFlagTypeINFO` (single per-tag check at Run
start; result is header-invariant).

The recode-column stripper code path (`filterRecodeInfo`'s
`removeInfo` parameter) was deleted as dead code: no CLI flag drove
it after the rewire, and no upstream-equivalent exists. The helper
now takes only the `--recode-INFO` keep set.

Pinned by `TestPassRemoveINFOSite`,
`TestRun_RemoveINFO_SiteFilter_Integration`,
`TestRun_RemoveINFO_SiteFilter_OR`,
`TestRun_RemoveINFO_SiteFilter_NonFlagType`, and
`TestRun_KeepAndRemoveINFO_Compose` in
`tools/vcftools/pkg/vcftools/info_filters_test.go`, plus the
byte-for-byte `TestParity_RemoveINFO_SingleFlag`,
`TestParity_RemoveINFO_OR`, and `TestParity_KeepAndRemoveINFO_Compose`
in `tools/vcftools/pkg/vcftools/parity_test.go` against goldens
generated by the same upstream binary used for wave 17.

The compose parity test deliberately omits `--recode-INFO-all` so the
INFO column collapses to `.` on both sides — this side-steps the
known port-vs-upstream INFO-key-ordering quirk
(`tools/vcftools/cmd/vcftools/aliases_cli_test.go:138-144`) and
keeps the parity check focused on the SITE SET, which is the
behaviour-of-record for this fix.

#### vcftools `--pca` jagged-M[i] crash on any missing genotype

Upstream's `output_PCA` (`variant_file_output.cpp:4954-4972`) appends
the centred/normalised genotype value to `M[ui_prime]` only when the
individual has a non-missing call at that site:

```cpp
e->get_indv_GENOTYPE_ids(ui, geno_id);
x = geno_id.first + geno_id.second;
if (x > -1) {
    if (use_normalisation == true)
        M[ui_prime].push_back((x - mu) * div);
    else
        M[ui_prime].push_back((x - mu));
}
ui_prime++;
```

…but the GRM accumulation at `variant_file_output.cpp:4988-4991`
iterates a SINGLE `s in 0..N_sites` index across every individual:

```cpp
for (unsigned int ui=0; ui<N_indvs; ui++)
    for (unsigned int uj=ui; uj<N_indvs; uj++)
        for (unsigned int s=0; s<N_sites; s++)
            X[ui][uj] += M[ui][s] * M[uj][s];
```

With any missing data, the per-individual `M[i]` vectors are
**jagged** (different lengths). The triple loop reads past the end of
the shortest vector — undefined behaviour that may segfault, return
garbage GRM entries, or read leftover heap memory depending on the
allocator. The bug is well-hidden because the Patterson/Price/Reich
2006 paper this code claims to implement explicitly **mean-imputes**
missing genotypes; upstream's "skip the push" implementation is a
silent divergence from that reference.

**Severity:** memory unsafety / silently-wrong output. With any
missing data the user gets garbage eigenvalues or a crash, depending
on the C++ runtime.

**Fixed in port** (wave 19 — this PR). The Go port drops any site
where at least one kept individual has a missing genotype, keeping
M rectangular. This is a deliberate deviation from upstream's
broken loop and from the Patterson method (which would mean-impute);
it matches upstream's effective intent ("compute PCA on
complete-data sites") while avoiding the indexing bug. The decision
is documented in `tools/vcftools/pkg/vcftools/pca.go` and pinned by
`TestPCA_MissingDataSkipsSite`.

Parity tests use no-missing-data fixtures so the divergence is
invisible to the byte-level comparison.

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

#### vcftools `--hapcount` prev_bin_idx shift on bin change

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
histogram ADDED to the old bin's histogram. Concretely: chr-1 bin 0
of a four-site fixture should report `N_SNP=4` but the buggy upstream
binary reports `N_SNP=1` (the single-site count from bin 1 leaking
backward).

**Severity:** numerical (silently-incorrect per-bin SNP counts and
merged multiplicity histograms in `.hapcount`).

**Fixed in port** (PR #140, wave 15). On every successful bin match
the runner FIRST flushes the OLD bin's data into its OWN slot (when
the bin index actually changed), then reassigns `binIdx`. This means
each bin's slot only ever receives that bin's own counts. Pinned by
`TestHapcount_CorrectBinTransitions` in
`tools/vcftools/pkg/vcftools/hapcount_test.go`. The `.expected.hapcount`
fixture was hand-traced from the corrected semantics, not generated
from the upstream binary.

#### vcftools `--hapcount` end-of-stream read-after-free

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

Either way the user silently loses real data for whatever chromosome
appears last in the VCF.

**Severity:** crash-on-input / silently truncated output (UB).

**Fixed in port** (PR #140, wave 15). The `hapcount.close()` handler
flushes any in-flight bin into its OWN slot and then emits the last
seen chromosome's rows unconditionally, using the kept-sample list we
already own (no freed-pointer access, no `have_data` discriminator).
Pinned by `TestHapcount_EndOfStreamFlush` in
`tools/vcftools/pkg/vcftools/hapcount_test.go`, which uses a VCF that
naturally ends mid-chromosome (no sentinel chrom needed).

#### vcftools `--hapcount` BED first-line silent skip

Upstream `variant_file_output.cpp:1183` runs
`BED.ignore(numeric_limits<streamsize>::max(), '\n');` BEFORE the BED
parse loop, unconditionally discarding the BED file's first line. A
user with a header-less BED therefore silently loses one bin (the
chr-1 bin 0 of a five-bin BED would silently become a four-bin BED
under upstream's reader).

**Severity:** silent data loss (one BED bin dropped per invocation
with a header-less BED).

**Fixed in port** (PR #140, wave 15). The runner inspects the first
BED line and skips it ONLY if it looks like a header (blank, or
starts with `#`, `track`, or `browser`). Otherwise the line is parsed
as data. Pinned by `TestHapcount_BEDFirstLineWithData` and
`TestShouldSkipBEDHeader` in
`tools/vcftools/pkg/vcftools/hapcount_test.go`. The auto-detection
covers the three header conventions commonly seen in BED files in the
wild.



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
  range. **Fixed in port** (this PR). The BCF reader at
  `pkg/htsgo/bcf/typed.go` decodes int64 typed payloads and downcasts
  values that fit (clamping the ±2³¹ sentinels onto the int32
  Missing/EndOfVector markers) so upstream-produced BCFs no longer
  trip the index/decode pipeline. CSI parity for upstream BCFs is now
  pinned by `TestParityIndex_CSIForBCF` in
  `tools/bcftools/pkg/bcftools/parity_test.go`. The encoder still
  emits int32 because no in-tree path produces values that need int64;
  if a future caller does, encoding will need a parallel int64 path.

<a id="bcf-fmt-keys-missing"></a>

- **Our BCF reader drops per-record FORMAT keys on htslib-produced
  input.** **Resolved** (this PR). Originally surfaced as
  `FmtKeys = [<resolved>, -1, -1, ...]` after the int64 + IDX-strip
  landings; the residual was actually masked by the int64 typed
  descriptor not being decoded — once `decodeTypedInternal` learned
  the int64 path (see `bcf-int64`) the FORMAT keys all resolve.
  Pinned by `TestParityView_BCFInput` in
  `tools/bcftools/pkg/bcftools/parity_test.go`, which asserts
  byte-for-byte equality against `view_basic.expected.vcf` on the
  upstream-produced `basic.bcf` fixture.

<a id="bcf-info-order"></a>

- **Our BCF writer does not preserve `InfoOrder` on encode.**
  **Fixed in port** (this PR). `Writer.encodeRecord` in
  `pkg/htsgo/bcf/typed_write.go` now walks `v.InfoOrder` first (any
  map-only keys are appended in sorted order for determinism) so a
  VCF→BCF→VCF cycle no longer shuffles INFO keys. The same change
  also surfaced a latent encoder bug where MissingInt32 was getting
  bitcast to a literal 0 at int8/int16 widths; the
  `encodeInts` / `encodeFormatTypedInts` switches now translate the
  sentinel to the width-appropriate `MissingInt8`/`MissingInt16` (and
  `EndOfVector*`) marker so downstream readers see the missing
  marker rather than a literal value. Pinned by
  `TestWriterPreservesInfoOrder` and
  `TestWriterPreservesMissingFormatInts` in
  `pkg/htsgo/bcf/writer_test.go`, plus
  `TestParityView_RoundTrip_OurBCF` in
  `tools/bcftools/pkg/bcftools/parity_test.go`.

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
  sampler (kr_drand-based MT19937-64, default seed 11). Our Go port
  used to do a "deterministic every-Nth-record" sampler that did
  not match upstream's selection regardless of seed.

- **Resolution.** Closed in the deferred-test wave. Ported upstream's
  MT19937-64 PRNG in-tree as
  `tools/seqtk/pkg/seqtk/krand.go` (matches `kr_srand` / `kr_rand`
  / `kr_drand` at seqtk.c:300-357 byte-for-byte). `Sample` now
  defaults to seed `SampleSeed = 11` and a streaming
  `kr.drand() < fraction` decision; `SampleSeeded(..., seed)` is
  available for explicit seed control. The parity test
  `TestParity_Seqtk_Sample_UpstreamByteParity` is unskipped and
  passes on `sample20.fq` at fraction 0.3.

### seqtk-randbase-rng <a id="seqtk-randbase-rng"></a>

- **Symptom.** Upstream `seqtk randbase` uses `drand48()` (with
  glibc's default implicit state X0 = 0) and is therefore
  deterministic across runs but not seed-controllable. Our Go port
  used `math/rand` with a caller-supplied seed.

- **Resolution.** Closed in the deferred-test wave. Ported glibc's
  drand48 in-tree as `tools/seqtk/pkg/seqtk/drand48.go` (48-bit
  LCG, default state 0, matches glibc byte-for-byte; verified
  against `man drand48`'s sample sequence). `Randbase` now uses
  `drand48State` and re-formats output to upstream's 60-char wrap
  rule (`if (i%60 == 0) putchar('\n')` at seqtk.c:557). The `seed`
  parameter is accepted for API stability but intentionally
  ignored. Parity test
  `TestParity_Seqtk_Randbase_UpstreamByteParity` passes on
  `ambig.fa`.

### seqtk-trimfq-algorithm <a id="seqtk-trimfq-algorithm"></a>

- **Symptom.** Upstream `seqtk trimfq` runs a modified Mott
  algorithm with an error-rate threshold (default `-q 0.05`) and a
  `-l 30` minimum-length floor. The port's `TrimQuality` does a
  simple Phred-quality threshold trim — a different feature.

- **Resolution.** Closed in the deferred-test wave. Added a new
  `TrimfqMott` function in `tools/seqtk/pkg/seqtk/seqtk.go` that
  ports upstream's main Mott pass plus the `imax / sliding-window`
  fallback (seqtk.c:397-426). The `seqtk trimfq` CLI now wires
  through `TrimfqMott` with `-q 0.05 / -l 30` defaults to mirror
  upstream's option-tail. The legacy `TrimQuality` library API is
  retained as a separate Phred-threshold helper. Parity test
  `TestParity_Seqtk_Trimfq_UpstreamByteParity` has two cases: an
  all-high-quality pass-through (`p64.fq`) and a low-quality-border
  case (`trimfq_borders.fq`) that actually exercises the trim — both
  pass byte-for-byte.

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
