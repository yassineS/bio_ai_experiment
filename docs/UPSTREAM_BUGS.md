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

#### bcftools +check-sparsity: `-R FILE` silently drops BED/TSV lines <a id="bcftools-check-sparsity-regions-file"></a>

- **Disposition:** Fix-on-port. Ours keeps the verbatim colon region-list parse
  byte-identical to upstream, but no longer silently drops BED/TSV lines — it
  parses them the synced-reader way (a strict superset).
- **What happens:** Unlike every other in-tree plugin (which loads `-R`/`-T`
  files through htslib's synced reader `bcf_sr_regions_init`, a TSV/BED parser),
  `check-sparsity` reads its `-R` file with `hts_readlist()` and hands each line
  verbatim to `tbx_itr_querys()`. `tbx_itr_querys` only understands the colon
  region-list syntax (`chr`, `chr:beg-end`). A tab-separated / BED line such as
  `chr1<TAB>0<TAB>10000` therefore fails to parse and the region is **silently
  skipped — no error, no output** — even though a normal BED file is exactly what
  a user would reach for.
- **Reproducer** (upstream 1.23.1, `BCFTOOLS_PLUGINS` pointed at the vendored
  `.so` dir; the input must be bgzipped + tabix-indexed because `-R` uses the
  index):

  ```sh
  printf 'chr1\t0\t10000\nchr2\t0\t10000\n' > all.bed
  bcftools +check-sparsity -n 1 -R all.bed gt_plugins.vcf.gz   # prints NOTHING
  printf 'chr2:1-10000\n' > colon.txt
  bcftools +check-sparsity -n 1 -R colon.txt gt_plugins.vcf.gz # prints "chr2:1-10000\tS3"
  ```

  The BED form produces no output; the colon-syntax line works. Compare with the
  synced-reader plugins (e.g. `+smpl-stats -R all.bed`) which accept the BED.
- **Our port:** `loadCheckSparsityRegionFile` keeps single-token lines (`chr`,
  `chr:beg-end`) verbatim — byte-identical to upstream, label and all — but
  instead of silently dropping a multi-column TSV/BED line it parses it the way
  htslib's synced reader / regidx does (`.bed`/`.bed.gz` => 0-based half-open;
  otherwise 1-based, two columns = a single position, three+ = `beg..end`),
  converting it to the equivalent 1-based `chr:beg-end` token (which also becomes
  the report label). The colon cases stay byte-parity with upstream
  (`TestNativePluginCheckSparsityRegion`); the BED/TSV fix is validated against
  the real upstream binary by feeding upstream the colon-equivalent region-list
  it can parse (`TestNativePluginCheckSparsityRegionBEDFixOnPort`).

#### bcftools-som-write-map

- **bcftools `som --train` always fails / `--classify` is unusable.**
  `vcfsom.c:170` (`som_write_map`) writes the map header with
  `fwrite("SOMv1",5,1,fp)` and checks the result `!=5`. `fwrite` returns
  the number of **elements** written (1 here, since `nmemb==1`), not the
  byte count, so the comparison is always true and the code calls
  `error("Failed to write 5 bytes\n")`, which `exit()`s with status 255
  after the file has been truncated to those 5 bytes. As a result every
  `bcftools som --train -p PREFIX ...` invocation aborts before writing a
  usable `PREFIX.som`, and `bcftools som --classify -p PREFIX` then fails
  with "Could not parse PREFIX.som" (also exit 255). `som --train` on a
  missing input file segfaults (exit 139). The `som` subcommand is
  therefore effectively dead upstream.

  Disposition: **fix-on-port (done).** We register `som` and fix the
  upstream write-map bug. The Go port's write path
  (`writeMaps` in `tools/bcftools/pkg/bcftools/som.go`) serialises the
  full map and validates the byte count of every write, so `--train`
  produces a usable map and the `--train`→`--classify` pipeline works
  end to end. Because upstream's on-disk `SOMv1` format is unusable, the
  port defines its own clean, versioned binary format (magic `SOMGO1`;
  documented in `som.go`) rather than reproducing the broken one. Two
  further deliberate divergences from a byte-exact port: (a) the SOM
  reads INFO annotations straight out of a VCF/BCF rather than a
  pre-extracted `annots.tab.gz`, and (b) weight initialisation uses Go's
  `math/rand` (deterministic per seed) instead of glibc's `random()`
  (TYPE_3 additive-feedback PRNG) — matching glibc's PRNG byte-for-byte
  only mattered for reproducing a tool that crashes. Validated by a
  train→classify round-trip, a map-file round-trip, and hand-checkable
  BMU / SOM-update unit tests in
  `tools/bcftools/pkg/bcftools/som_test.go` (no live oracle exists,
  since upstream crashes before writing a map). See
  `docs/PARITY_ROADMAP.md` (the `som` status note in the bcftools
  section).

#### vcftools-site-pi

- **vcftools `--site-pi` formula — RESOLVED: not a bug; byte-for-byte
  parity achieved.** Earlier notes (and issue #24) suspected upstream
  computed some per-genotype pairwise quantity that diverged from the
  textbook nucleotide diversity `(n² − Σ cₐ²) / (n(n-1))`. Reading the
  actual source settles it.

  Upstream `output_per_site_nucleotide_diversity`
  (`reference_code/vcftools/src/cpp/variant_file_output.cpp:3870`) does:

  ```cpp
  e->get_allele_counts(allele_counts, N_non_missing_chr);
  unsigned int total_alleles = accumulate(allele_counts...);   // == N_non_missing_chr
  unsigned int N_alleles = e->get_N_alleles();                 // ALT.size()+1
  int mismatches = 0;
  for (allele = 0; allele < N_alleles; allele++)
      mismatches += allele_counts[allele] * (total_alleles - allele_counts[allele]);
  int    pairs = total_alleles * (total_alleles - 1);
  double pi    = mismatches / (double) pairs;
  ```

  `get_allele_counts` (`entry_getters.cpp:395`) only ever reads the two
  diploid slots per sample, so `total_alleles == N_non_missing_chr == n`
  and `Σₐ cₐ·(n − cₐ) = n² − Σₐ cₐ²`. Therefore upstream's `pi` is
  **exactly** the textbook `(n² − Σ cₐ²) / (n(n-1))` — there is no
  per-genotype variant and **no bug**. The only things a naive port can
  get wrong are formatting and site selection:

  1. **Formatting.** Upstream writes `pi` straight to a default C++
     `ostream` (`out << pi`), i.e. `std::defaultfloat` precision 6 — so
     `0.6` not `0.600000`, and `0` not `0.000000`. Our port reproduces
     this with `formatFreq` (`%g`, 6 significant digits).
  2. **Site selection.** Upstream emits exactly the sites for which
     `entry::is_diploid()` is true (`entry_getters.cpp:94`), warning once
     "sitePi: Only using fully diploid sites." It therefore KEEPS a
     fully-missing diploid site (e.g. `20:1235237`, pi=0) and DROPS any
     site that has even one haploid included sample (e.g. chrX male
     calls). Our port applies the identical `siteIsDiploid` gate, and
     `n < 2` sites are dropped (division by zero / no pairs).

  Disposition: **fix-on-port complete, parity verified.** The default
  output matches the upstream 0.1.18 binary byte-for-byte. No
  `--site-pi-vcftools-compat` flag is needed because the default already
  equals both the upstream value and the published textbook definition.

  Parity test status: `TestParity_SitePi` (and the synthetic
  `TestParity_SitePi_EdgeCases`) build the upstream binary from the
  `reference_code/vcftools` submodule and assert byte-for-byte equality
  against it — a hard `t.Fatalf`, never `t.Skip`. A separate
  `TestParity_SitePi_Formula` spot-checks four hand-computed values
  offline.

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

#### vcftools spill-file outputs: 1-byte stack buffer overflow aborts the binary

Every output path that spills to a temporary file — the LD / chi-square
family `output_haplotype_r2` (`--hap-r2`), `output_genotype_r2`
(`--geno-r2`), `output_genotype_chisq` (`--geno-chisq`),
`output_interchromosomal_genotype_r2` (`--interchrom-geno-r2`),
`output_interchromosomal_haplotype_r2` (`--interchrom-hap-r2`), the
SNP-list-vs-all variants, AND the `--012` matrix writer
`output_as_012_matrix` (`variant_file_format_convert.cpp:377-405`) —
builds its temp-file name like this (`variant_file_output.cpp:1441-1443`
and the seven-plus sibling sites):

```cpp
string new_tmp = params.temp_dir+"/vcftools.XXXXXX";
char tmpname[new_tmp.size()];          // size N, NO room for the NUL
strcpy(tmpname, new_tmp.c_str());      // copies N+1 bytes → overflow
```

`char tmpname[new_tmp.size()]` is a variable-length array sized to the
string length, but `strcpy` writes `new_tmp.size() + 1` bytes (the
terminating NUL goes one past the end). On any toolchain that compiles
with `_FORTIFY_SOURCE` (the default for modern GCC/glibc, including the
build produced by the repo's own `./configure && make`), `__strcpy_chk`
detects the 1-byte overflow and aborts:

```
*** buffer overflow detected ***: terminated
```

**Severity:** the upstream 0.1.18 binary, built with its own autotools
config, **cannot run any `--hap-r2` / `--geno-r2` / `--geno-chisq` /
`--interchrom-*` / `--012` analysis at all** — it crashes before writing
a single data row (only the header / `.012.indv` prefix is flushed).
Reproduce: `vcftools --vcf any.vcf --hap-r2 --ld-window-bp 1000000` or
`vcftools --vcf any.vcf --012` → exit 134.

**Fixed in port** (this PR). The Go port computes these outputs in-memory
(no temp file, so no VLA/`strcpy` hazard) and emits the same column
layout upstream's writer would have produced, with the same C++
`defaultfloat` precision-6 formatting for the R²/D/Dprime/chi-square
columns. Because the upstream binary aborts, these modes cannot be
byte-validated against it; the port's LD and `--012` outputs are pinned
by the existing in-package unit tests (`ld_test.go`,
`ld_interchrom_test.go`, and the format-conversion tests) instead of a
live oracle.

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
  range. The BCF 2.2 spec allows it; in practice it just enlarges
  records and forces every consumer to add a 64-bit path. Not a bug per
  se but worth raising with the htslib maintainers about whether the
  encoder should clamp to int32 when it can. Tracked here so the next
  port can decide whether to round-trip 64-bit verbatim or downcast on
  read (we currently downcast — see `pkg/htsgo/bcf/typed.go`).

  Note: an earlier `bcftools index` CSI parity test was skipped citing
  this int64 descriptor as a blocker. That was inaccurate — the downcast
  means our reader handles upstream BCF (including int64 FORMAT counters)
  fine, so `TestParityIndex_CSIForBCF` now runs and asserts **functional**
  CSI parity (region reads return identical records on an upstream-produced
  BCF). Byte-identical `.csi` is a separate, deliberate non-target: htslib's
  `bcf_index` adapts the CSI depth to the longest contig and stores no aux
  block for BCF, whereas our CSI carries a small tabix-style aux so the
  reader can self-resolve contig names to ref IDs.

<a id="bcf-fmt-keys-missing"></a>

- **Our BCF reader drops per-record FORMAT keys on htslib-produced
  input** — RESOLVED (htsgo-gzi-bcf PR). The historical symptom was
  per-record `FmtKeys` decoding as `[<resolved>, -1, -1, ...]` (only the
  first key resolving, the rest as `MissingInt32`). This is no longer
  reproducible: `decodeIndiv` + `DecodeTypedInt` correctly read each
  FORMAT key (whether the dictionary index is int8- or int16-encoded)
  and `DecodeFormatTyped` advances the offset correctly across every
  per-sample value type, so subsequent keys resolve. Verified by live
  `bcftools view -O b` → our reader → VCF-text round-trips in
  `pkg/htsgo/bcf/fmtkey_parity_test.go` (`TestBCF_FormatKeyParity`),
  which exercises many numeric keys, string FORMAT fields, mixed ploidy,
  ragged int/float vectors with vector-end padding, phased GT,
  all-missing samples, ragged per-sample strings, and int16-encoded keys
  (dictionary padded past 127 entries). The reader no longer requires the
  VCF/VCF.gz workaround. (The `splitPerSample` TypeChar path was also
  hardened to slice per-sample slots by the descriptor's per-sample
  width `tv.Length` rather than `len(s)/nSample`.)

<a id="bcf-info-order"></a>

- **Our BCF writer does not preserve `InfoOrder` on encode** —
  RESOLVED. `bcf.Writer.encodeRecord` (typed_write.go) now emits INFO
  fields in the variant's recorded `InfoOrder` (then any stragglers in a
  stable name order), so a VCF→BCF→VCF round-trip preserves INFO key
  order. Pinned by `TestWriterInfoOrderDeterministic`
  (pkg/htsgo/bcf/writer_test.go).

- **Our BCF writer mis-encoded three value cases** — RESOLVED. (1) A
  missing INFO/FORMAT *integer* (`.`) was bit-truncated to `0` when the
  column narrowed to int8/int16 (`byte(int8(MissingInt32))` == 0) instead
  of the width's missing sentinel (`0x80`/`0x8000`); `narrowInt8`/
  `narrowInt16` now map the sentinels. (2) A missing GT allele was stored
  as the integer missing sentinel and only round-tripped because of that
  same truncation bug; `parseGT` now emits `bcf_gt_missing == 0` directly.
  (3) A Flag was encoded as a count-1 int8 of value 1, so htslib rendered
  it `TAG=1`; it is now the count-0 typed-int8 descriptor (`0x01`) htslib
  uses, rendering bare. A full VCF→our-BCF→VCF cycle now preserves
  records, and **upstream `bcftools` reads our BCF byte-equivalently**
  (developer cross-check). Pinned by `TestParityView_RoundTrip_OurBCF`
  and the updated `parseGT`/`encodeInfoValue` unit tests.

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

## samtools consensus <a id="samtools-consensus"></a>

### `samtools consensus --het-only` is a dead option (parsed, never read) <a id="samtools-consensus-het-only"></a>

**Severity:** documented contract violated — `--help` / the man page
advertise `--het-only` as "Only output heterozygous sites", but the flag
has no effect whatsoever on the output.

**Upstream behaviour (the bug).** In the vendored samtools
(`reference_code/samtools` at submodule commit
`e406d9eb0f7051b04c98a00f9aeb8f3de11b85dc`,
`bam_consensus.c`):

- the option is declared on the `consensus_opts` struct —
  `bam_consensus.c:232` (`int het_only;`),
- default-initialised to 0 — `bam_consensus.c:2997`
  (`.het_only = 0,`),
- registered in the long-option table — `bam_consensus.c:3049`
  (`{"het-only", no_argument, NULL, 6}`),
- and SET on parse — `bam_consensus.c:3097`
  (`case 6: opts.het_only = 1; break;`).

But `opts.het_only` is **never read** anywhere else: a repo-wide grep for
the symbol returns exactly those four sites (declaration, default, parse —
the getopt entry matches the literal string `"het-only"`, not the
variable). The consensus calling and output paths never branch on it, so
`samtools consensus --het-only` produces byte-for-byte identical output to
`samtools consensus` without the flag. This is a classic dead option: the
flag is accepted (no error) but silently inert, violating the documented
contract.

**Why it's a bug, not a design choice.** The flag's name and the man-page
text promise the consensus will be restricted to heterozygous-called
positions; a user passing `--het-only` reasonably expects homozygous and
no-call positions to be filtered out, and gets a full consensus instead —
silently wrong output for the user's intent.

**Fixed in port** (PR #221). Our Go `samtools consensus` implements the
intended behaviour: `--het-only` restricts output to HETEROZYGOUS-called
positions, suppressing homozygous and no-call positions.

- FASTA/FASTQ: suppressed positions render as `N` (coordinates
  preserved); leading/trailing non-het runs trim away like uncovered
  positions, unless `-a`/`-aa` forces full-length emission.
- pileup: suppressed rows are omitted entirely.
- Het-ness is determined INDEPENDENTLY of `--ambig`: in simple mode a
  position is het when `score2 >= het_fract*score1` on a confidently-
  called position; in bayesian mode when the het log-odds is positive on
  a confident call (depth/cutoff gated).

The fix lives in
`tools/samtools/pkg/samtools/consensus.go` (the `HetOnly` option, the
`consensusCall.isHet` field computed in `callConsensus` and
`callConsensusBayesian`, and the emit-loop suppression) and is wired on
the CLI in `tools/samtools/cmd/samtools/cmds_tail.go`.

**Tests.** Unit tests
(`TestConsensus_HetOnly_SuppressesHomozygous`,
`TestConsensus_HetOnly_AllHomozygous`,
`TestConsensus_HetOnly_OffIsUnaffected` in
`tools/samtools/pkg/samtools/consensus_test.go`) cover both calling modes
and both the `--ambig` and non-ambig paths. A live upstream test
(`TestConsensus_HetOnlyUpstreamBug` in
`tools/samtools/pkg/samtools/consensus_upstream_test.go`) builds the
vendored samtools and DEMONSTRATES the bug: it asserts upstream
`consensus --het-only` output is identical to upstream WITHOUT the flag
(proving upstream ignores it), and that our Go `--het-only` output DIFFERS
from our own no-flag output (the intentional, correct divergence). It also
confirms our baseline (no-flag) consensus still matches upstream's, so the
divergence is confined to the flag.

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

- **Symptom (historical).** Upstream `seqtk sample` uses a seeded
  reservoir sampler (krand / MT19937-64, default seed 11) and
  produces a deterministic, fraction-correct subset for any input.
  Our Go port originally implemented a "deterministic every-Nth-record"
  sampler that did not match upstream's selection regardless of seed.

- **Root cause.** Upstream's algorithm is the streaming reservoir
  pass over `n` records (`y = n-1 < num ? n-1 : r*n`, replacing slot
  `y` when `y < num`) with `keep = (kr_drand() < frac)` per record in
  fraction mode, and a two-pass random subset in `-2` mode. The old
  port short-circuited to `written/count < fraction` per record.

- **Disposition — RESOLVED (byte parity).** `SampleN` /
  `SampleFraction` (`tools/seqtk/pkg/seqtk/sample.go`) now port
  upstream's krand RNG and `stk_sample` exactly: streaming-reservoir
  fixed-number mode, two-pass (`-2`) mode, and per-record fraction
  mode, with `stk_printseq`'s no-wrap output. The CLI is wired to
  these. `TestParity_Seqtk_Sample_UpstreamByteParity` now asserts
  byte-for-byte against upstream-generated fixtures (fraction 0.3/0.5,
  number 5 at seeds 11 and 42, two-pass, and a FASTA case) and passes;
  it is no longer skipped. `TestParity_Seqtk_Sample_StructuralInvariants`
  remains as a complementary check of the legacy `Sample` helper.

### seqtk-randbase-rng <a id="seqtk-randbase-rng"></a>

- **Symptom (historical).** Upstream `seqtk randbase` uses `drand48()`
  (with the glibc default seed-0 state) and is therefore deterministic
  across runs but not seed-controllable. Our Go port used `math/rand`
  with a caller-supplied seed, so the output sequences differed.

- **Root cause.** Different RNGs.

- **Disposition — RESOLVED (byte parity).** `RandbaseUpstream`
  (`tools/seqtk/pkg/seqtk/mutations.go`) reimplements glibc's
  `drand48` (the 48-bit LCG `X = (0x5DEECE66D*X + 0xB) mod 2^48`,
  default `X = 0`) and reproduces `stk_randbase`'s output layout
  exactly: only 2-base IUPAC codes are drawn (`m = drand48() < 0.5`),
  3/4-base codes and N pass through, the comment is dropped, and the
  sequence is wrapped at 60 columns. The CLI uses this path by default
  (no `-s`); `-s INT` selects a seeded `math/rand` extension.
  `TestParity_Seqtk_Randbase_UpstreamByteParity` now asserts
  byte-for-byte against upstream fixtures and passes (no longer
  skipped); the structural-invariants test is retained for the seeded
  helper.

### seqtk-trimfq-algorithm <a id="seqtk-trimfq-algorithm"></a>

- **Symptom (historical).** Upstream `seqtk trimfq` runs a modified
  Mott algorithm with an error-rate threshold (default `-q 0.05`) and
  a `-l 30` minimum-length floor. The legacy `TrimQuality` helper did
  a simple Phred-quality threshold trim, producing different cuts.

- **Disposition — RESOLVED (byte parity).** `TrimFQ`
  (`tools/seqtk/pkg/seqtk/trimfq.go`) ports `stk_trimfq`'s modified
  Mott algorithm (including the `q_int2real` table, the
  `[36,127]`-clamped per-base sum, and the window-based fallback when
  the Mott window is shorter than `-l`) and the fixed-offset path
  (`-b`/`-e`/`-L`). The CLI default path is now Mott.
  `TestParity_Seqtk_Trimfq_UpstreamByteParity` asserts byte-for-byte
  against upstream fixtures (short-read pass-through, Mott on 40 bp
  reads, `-b/-e`, and `-L`) and passes (no longer skipped). The
  legacy `TrimQuality` Phred-threshold helper remains available but is
  off the default CLI path.

## prinseq

The prinseq parity audit (PR for `prinseq-seqtk-parity-validation`)
fixed three discrepancies on **our side** (not upstream); they are
listed under
[tools/PARITY_VALIDATION.md > prinseq parity validation](../tools/PARITY_VALIDATION.md).
No upstream bugs surfaced — PRINSEQ-lite's documented behaviour
agreed with the corpus we tested for every option we exercised
(see the prinseq table for the 18 cases).

### `--seq_id` trailing-comment divergence (resolved on our side)

**Severity:** behavioural (our port only); **Status:** fixed.

Our earlier `--seq_id` implementation rewrote the header to just
`<prefix><N>`, dropping any trailing FASTA/FASTQ comment. Upstream
`prinseq-lite.pl:3685-3704` emits `$sid.($header ? ' '.$header : '')`,
i.e. it replaces only the identifier token and re-appends the original
comment. The PR
`claude/festive-planck-n9o2lm-prinseq-transforms-misc` resolves this:
`renameDescription` now rewrites only the id and `joinHeader` re-attaches
the comment, so `>read1 sample=A` with `--seq_id S_` becomes
`>S_1 sample=A`, matching upstream byte-for-byte (verified by
`TestParityTransforms/seq_id_fasta` and `/seq_id_fastq`). The mapping TSV
written by `--seq_id_mappings` still keys on the bare identifier (`$sid`),
matching upstream line 3646.

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

## bcftools gtcheck: `-e EXPR` silently forces integer scoring

Upstream `vcfgtcheck.c`'s option parser, in `case 'e'`, runs
`args->gt_err = strtol(optarg, &tmp, 10)` before deciding whether the
argument is a filter expression or the deprecated `-e` error-probability
value. For a non-numeric filter such as `-e 'INFO/AC<4'`, `strtol`
returns 0 and the code then (correctly) treats the argument as an
exclude filter — but the `gt_err` side effect has already happened, so
the error probability is left at 0. A `gt_err` of 0 selects the
integer-mismatch scoring path (equivalent to `-E 0`), so the discordance
column switches from the abstract floating-point score to a raw integer
mismatch count purely as a side effect of using `-e` with a filter. The
`qry:`/`gt:`-prefixed form and `-i` (include) are unaffected — only the
bare `-e EXPR` branch runs the `strtol`.

**Our behaviour:** `-e EXPR` applies the exclude filter without touching
the error-probability / scoring mode; the scoring stays whatever `-E` /
`-u` (or the default) selected, which is what a user reasonably expects.
A user who genuinely wants integer scoring passes `-E 0` explicitly. The
live parity tests use the `qry:`-prefixed `-e` form so both sides agree;
the bare-`-e` scoring-mode flip is the documented, intentional deviation.

## bedtools merge/groupby `distinct_only`: spurious leading delimiter

`KeyListOpsMethods::getDistinctOnly()` (reference_code/bedtools, KeyListOps)
walks the value-string-sorted `freqMap` and appends a delimiter before every
element except the FIRST ENTRY OF THE MAP — not the first element actually
emitted:

```cpp
for (; _freqIter != _freqMap.end(); _freqIter++) {
    if (_freqIter->second > 1) continue;            // skip repeated values
    if (_freqIter != _freqMap.begin()) _retStr += _delimStr;
    _retStr.append(_freqIter->first);
}
```

When the first map key has frequency > 1 it is skipped, but the very next
emitted (frequency-1) value still sees `_freqIter != _freqMap.begin()` as true
and so is prefixed with a delimiter. The result has a leading comma, e.g. for
values `3,1,10,3,1` (counts 1→2, 3→2, 10→1) upstream prints `,10` instead of
`10`.

**Our behaviour:** `distinct_only` emits only the genuine frequency-1 values
(value-string sorted) with no spurious leading delimiter — `10` for the example
above. This is the documented, intentional deviation; all other KeyListOps
operations (`distinct`, `concat`, `distinct_sort_num[_desc]`, `freqasc`,
`freqdesc`, …) match the live upstream binary byte-for-byte.

## bedtools bamtobed `-tag` + `-split`: malformed extra column

`bamToBed.cpp` `PrintBed()` has two split-mode branches. When no `-tag` is
given it correctly prints `chrom start end name mapq strand` (6 columns, one
per block). But the `-tag` branch prints an extra `bam.Position` column and
streams the block start/end after it:

```cpp
cout << chrom << "\t"
     << bam.Position << "\t"   // <- spurious extra column
     << curr.start << "\t"
     << curr.end << "\t"
     << name << "\t"
     << PrintTag(bam, bamTag) << "\t"
     << strand << endl;
```

So `bamtobed -tag NM -split` emits a 7-column line `chrom pos start end name
tag strand` instead of the intended BED6 `chrom start end name tag strand`.

**Our behaviour:** reproduced exactly for byte-for-byte parity — pipelines
that already consume this 7-column shape keep working. `bedbamtobed` only
emits the spurious column in the `-tag`+`-split` combination, matching the
live upstream binary; the plain `-split` and non-split `-tag` paths are clean
BED6.

## bedtools bedtobam: out-of-genome chromosome silently maps to ref id 0

`bedToBam.cpp` `ConvertBedToBam()` resolves the reference id with
`bam.RefID = chromToId[bed.chrom];`. Because `chromToId` is a `std::map`,
`operator[]` on a missing key inserts a default-constructed value of `0` — so a
BED record whose chromosome is absent from the `-g` genome file is written
against the FIRST `@SQ` reference rather than being rejected or skipped. For
example, with a genome listing only `1`, an input on `chrX` produces a record
on `1` (ref id 0), not an error.

**Our behaviour:** reproduced for byte-for-byte parity — an unknown chromosome
is emitted against the first reference in the genome file (requires a non-empty
genome). `bedtobam` matches the live upstream binary here; an empty genome file
is an explicit error rather than a crash.

## bcftools +parental-origin: NULL deref (segfault) on a site-only `-e`/`-i` whose expression matches no site

`plugins/parental-origin.c` `process_record()` requests the per-sample mask
from `filter_test(args->filter, rec, &smpl_pass)`. For a *site-only*
expression (e.g. `QUAL<10`, no FORMAT term) htslib leaves `smpl_pass == NULL`.
In the `FLT_EXCLUDE` branch, when the site does **not** pass the expression
(`pass_site == 0`) the code falls into

```c
else
    for (i=0; i<3; i++) smpl_pass[args->trio.idx[i]] = 1;
```

which dereferences the NULL `smpl_pass` and crashes. So
`bcftools +parental-origin -p P,F,M -t dup -e 'QUAL<10' file.vcf` (where every
site has `QUAL>=10`, so nothing matches the exclude expression) **segfaults**
in upstream 1.23.1. The symmetric `-i` site-only case where everything matches
is fine; only the exclude-with-no-match (and include-with-no-match, which hits
the same NULL write through a different path) crashes.

**Our behaviour:** fixed-on-port. A site-only expression that selects no
per-sample mask is treated as a whole-site verdict: under `-e`, a site that
does not match the expression is *kept* with all three trio members included
(the intended semantics), and under `-i` a non-matching site is dropped. No
NULL write occurs, so the run completes and prints its summary instead of
crashing. The byte-parity tests therefore exercise the *non-crashing* filter
forms (`-i QUAL>10`, per-sample `FMT/GQ` expressions) against the live
upstream binary; the crashing form is asserted only to *not* crash in our
port, since upstream produces no output to compare against.
