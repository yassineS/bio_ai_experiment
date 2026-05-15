# Parity roadmap

**Goal:** **1:1 feature parity** with the upstream tool for every Go port in
this repo. This file is the authoritative gap list per tool.

The project's stated goal is to make these tools faster, better tested, and
better documented than their originals — which requires that we actually
implement the same features. "Working subset" is a milestone, not the
destination.

> The companion file [`UPSTREAM_BUGS.md`](UPSTREAM_BUGS.md) tracks bugs
> we identify in the originals as we go. We don't carry those over; we
> either fix on port or note and skip with a clear pointer.

The companion file [`../analysis/tool_ranking_2026.md`](../analysis/tool_ranking_2026.md)
ranks **next** tools to port. It is **not** a "deprioritise existing tools"
filter — existing tools all get carried to 1:1.

## Definition of "1:1"

A tool is **1:1** when:

1. Every subcommand the upstream binary exposes is present in the Go port.
2. Every documented flag/option is recognised (either implemented or
   gracefully rejected with a clear error pointing at this roadmap).
3. For every input shape upstream accepts, our port produces the same
   **logical result**, modulo:
   - documented intentional deviations (recorded in the tool's README +
     this roadmap),
   - upstream bugs we chose to fix (recorded in [`UPSTREAM_BUGS.md`](UPSTREAM_BUGS.md)),
   - **RNG / stochastic divergence** — see the policy below.
4. The validated-parity test suite (runs the upstream test corpus through
   our port) passes for every supported case, with explicit `t.Skip()` for
   each documented exception.

### RNG / stochastic-output policy (2026-05-15)

For subcommands that use randomness (`bedshuffle`, `bedsample`, `seqtk sample`,
`seqtk randbase`, `samtools view -s subsample`, etc.) the parity bar is:

- **Reproducibility within our tool**: same seed + same input → same output
  byte-for-byte, every time, across Go versions and platforms.
- **Logical equivalence with upstream**: our output must satisfy the same
  invariants upstream's output does (correct sampling fraction, no
  duplicates without replacement, strand filters honoured, etc.).
- **NOT** byte-identical with upstream's output. Upstream uses its own
  C/C++ RNG (typically a Mersenne Twister from libstdc++); Go uses
  `math/rand`. Porting upstream's RNG would be ~200 lines of
  bit-twiddling per tool with no functional benefit — the user
  explicitly opted out of that work in favour of focusing on real
  feature parity.

The parity-test infrastructure handles this by either:

- structural-invariant assertions (e.g. "every shuffled interval has the
  same length as the input; every shuffled interval is on a chrom in the
  genome"), or
- a documented `t.Skip("RNG byte-parity, see PARITY_ROADMAP.md#rng-policy")`
  with a pointer to this section.

We're not there yet for any tool. The bedtools subset (PR #55) is the
closest — 127 parity tests, 85 passing, 42 documented `t.Skip` — and even
there we have ~30 subcommands not yet started.

---

## Per-tool gap list

Numbers reflect state at 2026-05-14 (post-#71). Update when each gap is
closed.

### `seqtk`

**Status:** 11 of ~20 subcommands. ~55%.

Missing subcommands:

- `gap` — find gap regions in FASTA.
- `listhet` — extract heterozygous sites from VCF/BCF.
- `gc` — find GC-rich regions.
- `fqchk` — FASTQ quality check report.
- `seqshuf` — shuffle FASTA/Q records.
- `pair` — pair up R1/R2 from interleaved input.
- `dropse` — drop unpaired reads.
- `hpc-bg` — homopolymer-compress with mismatch tolerance.
- `kfreq` — k-mer frequency analysis.
- `cnregion` — find regions of constant base.
- `gcdc` — GC depth count.

Option-tail gaps (per existing subcommand):

- `comp` — missing `-r REGION` to restrict to a BED region.
- `seq` — missing `-A` (force ASCII output), `-C` (mask sequence with N), `-M FILE` (mask regions), the `-T int` trim option.
- `sample` — missing `-2` (output two paired files).
- `trimfq` — missing `-L int` (max length cap), `-B int` (min base quality).
- `subseq` — missing the regex-name mode.
- `mutfa` — missing the inverse `--inverse` mode.

**Validation:** no upstream-test-suite run yet.

### `prinseq-lite`

**Status:** 2 subcommands (`stats`, `filter`) covering most common knobs.

Missing flags (per the upstream `-help` listing — incomplete inventory; see
`reference_code/prinseq-lite/` for the full Perl source):

- `--graph_data` and the corresponding HTML/PNG report generation.
- `--out_format` for FASTA/FASTQ/QUAL/SOL conversion combinations.
- `--seq_id_mappings` for ID-mapping output.
- `--ns_max_p` (percentage form of `--ns_max_n`).
- `--noniupac` strict-IUPAC filtering.
- `--phreds` quality-score output format.

**Validation:** no upstream-test-suite run yet. Upstream PRINSEQ-lite has no
formal test suite — we'd need to construct one from the documented examples.

### `sickle`

**Status:** 2 subcommands (`se`, `pe`). Claimed "Complete" but never
validated byte-for-byte against the C original.

To do:

- Run upstream sickle on a corpus of test FASTQ and diff our output.
- Verify the quality-encoding auto-detect (PR #33) tracks upstream's
  heuristic exactly.
- Check `-q`/`-l` boundary semantics for off-by-one with upstream.
- Confirm gzip output level matches upstream's `gzip -6` default.

### `skewer`

**Status:** 2 subcommands (`se`, `pe`). Claimed "Complete" but never
validated byte-for-byte against the C++ original.

To do:

- Run upstream skewer on the adapter-trimming corpus from the original
  paper and diff our output.
- Confirm log-format compatibility (`-x`, `--quiet`, the `.log` file
  format).
- Check that `--mode` (`any`, `head`, `tail`) matches upstream
  semantics for each.

### `fastp`

**Status:** Single `fastp` command with sliding-window cut, auto adapter
detection, HTML+JSON reports, duplication evaluation, UMI processing.

Missing:

- **Overrepresented sequence analysis** (`--overrepresentation_analysis`,
  `--overrepresentation_sampling`).
- **Base-correction in PE overlap** (`--correction`).
- **Quality-trimming overlap mode** (`--overlap_len_require`,
  `--overlap_diff_limit`, `--overlap_diff_percent_limit`).
- **PolyG/polyX more knobs**: `--poly_g_min_len`, `--poly_x_min_len`.
- **Splitting output**: `--split`, `--split_by_lines`, `--split_prefix_digits`.
- **Adapter list output to FASTA**: `--adapter_fasta`.
- **Disable adapter trimming**: `--disable_adapter_trimming`.
- **JSON schema completeness**: a few sub-fields under
  `before_filtering`/`after_filtering` and the per-cycle base content
  arrays are present but a handful of additional keys upstream emits
  are still missing. Run upstream `fastp` on a sample input and diff.

#### Validated-parity audit (this PR)

15-case test corpus at `tools/fastp/pkg/fastp/parity_test.go` against
upstream fastp 1.0.1. 11 PASS, 4 SKIP. See
[tools/PARITY_VALIDATION.md#fastp-parity-validation](../tools/PARITY_VALIDATION.md#fastp-parity-validation)
for the case list.

Bugs in the Go port surfaced + fixed inline by this audit:

- **UMI tag format** was unconditionally `":UMI_<umi>"`. Upstream uses
  `":<umi>"` (no prefix) or `":<prefix>_<umi>"` (with prefix). Fixed.
- **Low-complexity definition** was "unique 2-mers / total 2-mers".
  Upstream uses "fraction of adjacent positions where seq[i] !=
  seq[i+1]". Fixed.
- **`low_complexity_reads` JSON counter** was missing. Added.

Bugs in the Go port we **identified but did NOT fix in this PR** (skipped
parity cases pointing back here):

- **PolyG mismatch tolerance**: upstream's `trimPolyG` tolerates 1
  mismatch per 8 bases scanned (capped at 5 total) and anchors on the
  last-G position (`reference_code/fastp/src/polyx.cpp::trimPolyG`).
  Our Go port does a strict consecutive-G count. Follow-up needed:
  port the upstream algorithm verbatim — it's ~20 lines.
- **Sliding-window boundary** (`cut_front` / `cut_tail` / `cut_right`):
  three off-by-1..2 issues in `slidingWindowCut` vs upstream's
  `filter.cpp::trimAndCut`. cut_right needs to keep the high-Q prefix
  of the offending window; cut_front needs to skip past trailing N's
  at the cut; window-iteration bounds differ by one. Follow-up: port
  the upstream algorithm; ~50 lines total across the three modes.
- **SE adapter auto-detect**: upstream builds a kmer overlap-tree from
  the first 10000 reads (`evaluator.cpp`). We do a simple substring
  search against a small built-in adapter table. Different algorithm,
  different results. Bigger fix; tracked here.

**Validation:** **16-case parity test suite, 12 passing, 4
documented `t.Skip`** (post this PR).

### bedtools (35 subcommands ported)

**Status:** 35 of ~40 subcommands (~88%). 141 passing parity tests
against the upstream test suite (across PR #55 + Phase-3 wave 1 + wave
2 simple + wave 2 algo) + 17 new cases from wave 3 (PR #87) + 6 cases
from the reldist/fisher full-parity wave (PR #90) + 6 cases from the
nuc/annotate wave (PR #91) + 7 cases from the multicov/multiinter
wave (PR #92) + 6 cases from the igv/links wave (PR #93) + 10 cases
from the BEDPE pair-ops wave (PR #94, five each for `bedpairtobed`
and `bedpairtopair`, derived from the upstream source as neither
subcommand ships a `test/<name>/` subdir) + **9 newly passing**
parity cases from the column-ops + discrepancies wave (this PR —
`jaccard.t02/t03/t05/t06/t10/t11`, `merge.t15`, `map.t11`,
`map.t13`).
Phase-3 wave 1 (PR #78) added `bedgroupby`/`bed12tobed6`/`bedmakewindows`;
wave 2 simple (PR #80) added `bedexpand`/`bedgetfasta`/`bedsample`/`bedspacing`;
wave 2 algo added `bedcoverage`/`bedmap`/`bedshuffle`; wave 3 tail
(PR #87) added `bedcluster`/`bedsplit`/`bedsummary`/`bedtag`/`bedwindow`;
the reldist/fisher wave (PR #90) added `bedreldist`/`bedfisher`; the
nuc/annotate wave (PR #91) added `bednuc`/`bedannotate`; the
multicov/multiinter wave (PR #92) closed the last two of the six
originally-planned algorithmic subcommands; the igv/links wave
(PR #93) landed the two pure-format converters; **the BEDPE pair-ops
wave (this PR) closes the last two missing subcommands — bedtools now
has no missing subcommands**.

Missing subcommands: none. All algorithmic, BEDPE, and converter
subcommands are ported. Remaining bedtools work is option-tail polish.

Resolved in the column-ops + discrepancies wave:

- `bedjaccard` now pre-merges A and B before computing intersection /
  union, matching upstream's `setUseMergedIntervals(true)`. Previously
  skipped parity cases `jaccard.t02 / t03 / t05 / t06 / t10 / t11` now
  pass byte-for-byte.
- `bedmerge -s` now matches upstream's per-strand merge semantics: `.`
  (unknown) strand records are dropped under `-s` and `+` / `-` groups
  merge independently before the (chrom, start) re-merge. Previously
  skipped `merge.t15` now passes.
- The full column-op vocabulary from upstream (`stdev`, `sstdev`,
  `absmin`, `absmax`, `cat`, `cat_uniq`) is wired into the shared
  `bedmerge.ApplyOp`, unblocking `bedgroupby` (`TestGroup_StdevSstdev`,
  `TestGroup_AdditionalOps`) and `bedmap` (`map.t11`, `map.t13`).

Option-tail gaps on the wave-2 additions:

- `bedgetfasta` — `-fullHeader` is now implemented: contigs are indexed
  by the full FASTA header line (whitespace included) via
  `pkg/bioformats/fasta.BuildIndexFullHeader` /
  `OpenRandomAccessFullHeader`, and `bedgetfasta -fullHeader` flows the
  flag through to the index build. Upstream `getfasta.t06` (the
  `-fullHeader` two-line case) and `t07` (the no-`-fullHeader` warning
  case) now both pass byte-for-byte. BGZF FASTA input is also wired
  through (this PR): `pkg/bioformats/fasta` now sniffs the BGZF magic
  in `OpenRandomAccess` / `OpenRandomAccessFullHeader` and routes to a
  new `OpenRandomAccessBGZF` that fully decompresses the payload
  in-memory and reuses the existing FAI index path. The `.gzi` sidecar
  (when present) is parsed for early validation via a stdlib-only
  little-endian reader in `pkg/bioformats/fasta/bgzf.go`; a samtools-
  compatible `.fa.gz.fai` is honoured when present, otherwise the
  index is rebuilt from the decompressed payload. Upstream
  `getfasta.t18` (BGZF FASTA + `-split` BED12) now passes byte-for-byte
  using the upstream `t.fa.gz` fixture. Partial-decompression seek via
  `.gzi` is a future optimization; the in-memory path is sufficient
  for parity and for the reference genomes bedtools is typically used
  against.
- `bedsort` — `-header` is now implemented (this PR): leading
  `#`-prefixed comment, `track`, and `browser` directive lines are
  buffered and emitted verbatim ahead of the sorted body. Upstream
  `sort.t09` now passes byte-for-byte.
- `bedsample` — output PRNG is Go `math/rand` and is not byte-compatible
  with upstream's C++ sampler. Seeded runs are deterministic within
  `bedsample` (same seed → same output) but cross-tool record-for-record
  parity with upstream is not feasible without porting upstream's
  `random_shuffle`.
- `bedmulticov` — <a id="bedmulticov-bam"></a>BAM input is wired through
  via `pkg/bioformats/sam.NewBAMReader`; primary alignments contribute
  one interval each over their reference span, and `-q` MAPQ filter +
  `-D` per-A-interval depth cap are honoured. Upstream `multicov.t1`
  through `t4` and `t10` pass byte-for-byte.
  **`-split`** block-aware coverage on BAM CIGAR `N` ops is now
  implemented (this PR): the BAM index pass walks each alignment's
  CIGAR and emits one block per contiguous reference-consuming op-run
  (M/=/X, with D extending the current block — matching upstream's
  `breakOnDeletionOps=false`), skipping any `N`-op gap. Each alignment
  is counted at most once per A interval. When combined with
  `-f`, the threshold is applied to `total_block_overlap /
  sum_of_BAM_block_lengths` using strict `>` — a quirk of bedtools 2.x
  preserved here for byte-for-byte parity (mirrors
  `multiBamCov.cpp::FindBlockedOverlaps`). Upstream `multicov.t5`
  through `t9` now pass.
  **CRAM** input remains deferred — see `docs/CRAM_DESIGN.md`; the CLI
  surfaces a clear error for `.cram`.
- `bedmultiinter` — VCF/GFF input not implemented (upstream autodetects
  these via `BedFile`). Input is assumed sorted; out-of-order records
  within a single file are tolerated only because each file is
  re-sorted and merged before the sweep.

Column-op closure: the shared `bedmerge.ApplyOp` (used by `bedmerge`,
`bedgroupby`, `bedmap`, and `bedcoverage`) now supports the full
upstream KeyListOps vocabulary: `stdev`, `sstdev`, `absmin`, `absmax`,
`cat`, `cat_uniq` (in addition to the previously-implemented sum,
min/max, mean, median, mode/antimode, count, count_distinct, distinct,
collapse, first, last). Done; no remaining gaps.

### `vcftools`

**Status:** ~70 of ~147 options (~48%) after long-tail wave 1.

Closed in wave 1 (this PR):

- **Inter-chromosomal LD**: `--interchrom-geno-r2`, `--interchrom-hap-r2` ✅
- **Chi-square LD**: `--geno-chisq` ✅
- **Relatedness**: `--relatedness` (Yang 2010), `--relatedness2` (KING-robust) ✅
- **Runs of homozygosity**: `--LROH` (+ `--LROH-min-variants`) ✅
- **Phased blocks**: `--phased-blocks` ✅
- **FILTER tag include/exclude**: `--remove-filtered`, `--keep-filtered` ✅
- **INFO selection in recode**: `--keep-INFO TAG`, `--remove-INFO TAG` ✅
- **INFO extraction**: `--get-INFO TAG[,TAG]` → `.INFO` ✅

Remaining gaps:

- **Mendelian inheritance checks**: `--mendel`.
- **Diff family extensions**: `--diff-indv-map`, `--diff-discordance-matrix`,
  `--diff-switch-error`, `--gzdiff` (already implicit via iohelper).
- **Output formats**: missing `--ldhat`, `--ldhat-geno`, `--ldhelmet`,
  `--IMPUTE`, `--phase` output paths.
- **Haplotype analyses**: `--haploid` (`--phased-blocks` done).
- **Per-individual output**: `--missing-per-ind` (we have `--missing-indv`
  but the per-individual `.imiss` row layout has fields we don't emit).
- **Other**: `--FILTER-PASS-summary`, `--remove-INFO-all` (use
  `--keep-INFO`/`--remove-INFO`), `--non-ref-af*` family, `--pca` family.

**Validation:** wave 1 adds header byte-for-byte parity tests for the new
output files; full upstream-test-suite run still pending.

### `bgzip`

**Status:** 1 / 1 command, most flags.

Missing:

- **Multi-threaded compression** (`-t / --threads N` is accepted but
  single-threaded; BGZF is trivially parallel per block).
- **Output-rename to follow upstream conventions on stdin**: minor.

**Validation:** round-trips through `tabix` work; no full upstream-test
suite run yet.

### `tabix`

**Status:** 1 / 1 command, most flags.

Missing:

- **`--reheader FILE`** — replace bgzipped file's header lines in place.
- **`--targets` strictness** — currently behaves as `-R`; needs to be a
  true post-filter that only emits records strictly inside the targets.

**Validation:** no full upstream-test-suite run yet.

### `samtools`

**Status:** 21 of ~25 subcommands (~84%). `view`, `sort`, `index`, `depth`,
`fastq`, `flagstat`, **`mpileup`** (wave-1 + tail wiring), PR #88's
wave-1 tail (`merge`, `coverage`, `idxstats`, `cat`, `reheader`,
`addreplacerg`, `fixmate`, `dict`, `split`, `quickcheck`), the
heavy-hitter pair `markdup` + `stats`, and the pair landed in
the calmd/import PR: **`calmd`** + **`import`**.

Missing subcommands (in rough priority order):

- **`consensus`** — base-level consensus.
- **`phase`** — phase reads with their mates.
- **`targetcut`** — cut targeted regions.
- **`tview`** — terminal viewer (likely a deliberate skip; ~no-one uses it).
- **`view` flag-tail**: `-X` (custom-index input). `-L bed` landed as a
  linear-scan BED-region filter; `-M`/`--use-multi-region-iterator` is
  accepted but treated as a no-op since we always run the full
  intersection. `-d/-D` (tag-value filter) and `-N` (qname file) landed
  in the view-d-D-N PR.
- **`mpileup` tail** beyond PR #88 wiring: BCF output, `-g/-u` genotype-
  likelihood mode. `-aa` zero-fill of empty contigs is implemented (see
  `TestMpileup_AA_ZeroFillTableDriven`).

Plus:

- **CRAM** read/write throughout. Multi-month effort on its own; the
  rANS codec layer is the gating piece. Owner has OK'd third-party
  dependencies for CRAM codecs (see `CLAUDE.md#documented-exception-cram`).
  Design doc to follow.
- **`.csi`** for samtools (BAI is fine for chromosomes ≤512Mb).
- **Multi-threading** in `sort`, `index`, `view` (`-@`).

**`markdup` deferred features** (deliberately skipped in v1, all flag
slots are accepted on the CLI for compat):

- Optical-duplicate detection (`-d/--max-dist` + (x,y) parsing of Illumina
  qnames). v1 marks PCR duplicates only; nonzero `-d` triggers a stderr
  warning.
- Per-read-group keying (upstream's `-S` flag). v1 folds all read groups
  into a single namespace, so fixture
  `reference_code/samtools/test/markdup/17_read_group.sam` is a
  documented partial-parity skip.
- Barcode regex / barcode-tag keying.
- The `dt:Z:` "duplicate-type" aux tag (SQ / LB / OQ). The 0x400 flag
  bit is set correctly; only the typed aux is missing.

**`markdup -l/--max-len` is a no-op-by-design.** Upstream uses `-l` solely
as the streaming buffer flush window in `bam_markdup.c:1949`; it does NOT
affect key construction or scoring. Our two-pass implementation buffers
per-bucket state in memory, so output is identical for any `-l` value.
The flag is accepted on the CLI and the option is preserved on
`MarkdupOptions.MaxLen` for forward compatibility if we ever move to a
single-pass streaming model.

**`stats` deferred sections** (also documented in
`PARITY_VALIDATION.md`):

- COV/COV2 coverage histograms (require reference + BAI).
- GCD/GCT/GCC/GCL GC distributions (require reference bases).
- FFQ/LFQ per-cycle quality matrices and OXC oxidation-context counts.
- `--target-regions BED` restriction.
- The leading CHK checksum block (CRC32 reduction of read names /
  sequences / qualities).
- BWA-style quality trimming (`-q/--trim-quality`). The SN field
  `bases trimmed` is reported as 0; upstream also reports 0 when the
  flag is not passed, so byte parity holds for the default invocation.
- Mate-tracking memory cap: upstream's `cleanup_overlaps` periodically
  evicts stale `mates` entries. Our `mates` map currently grows
  unbounded — fine for the typical workload but worth fixing before
  running `stats` on multi-billion-record BAMs.

v1 emits byte-faithful **SN** (Summary Numbers) and useful **RL / MAPQ /
IS** rollups; the unsupported sections are quietly omitted (or, under
`--sparse`, all histogram blocks are suppressed entirely).

**Validation:** upstream fixtures from `reference_code/samtools/test/markdup/`
and `.../test/stat/` are vendored under
`tools/samtools/testdata/parity/{markdup,stat}/`. The byte-exact /
flag-exact / SN-byte cases are exercised in
`tools/samtools/pkg/samtools/markdup_test.go` and `stats_test.go`.

**`calmd` deferred features** (accepted as CLI flags, behaviour partial):

- **BAQ recalculation** (`-r`, `-E`, `-A`). Upstream's `bam_md.c` calls
  `sam_prob_realn` to recompute the BQ (base-alignment quality) aux
  tag and to drop MAPQ for low-quality reads. v1 fills in MD + NM
  correctly but does not touch BQ or MAPQ. ~200 lines of HMM-style
  alignment math; deferred per owner steer.
- **`-h` HASH_QNM** (hash-based query-name binarisation) — niche
  upstream-only optimisation; not implemented.
- **`-d` DROP_TAG** (drop all aux but RG) and **`-q` BIN_QUAL** (round
  qualities to 0/7) — flag-recognised in the CLI driver but not
  threaded through; safe to add later as small wrappers around the
  current calmd pipeline.
- **`-n` max-NM cap** — would mask high-mismatch reads with bin-quality;
  trivial follow-up once BIN_QUAL lands.
- **`-N` clear-MD/NM-bits**, **`-C` capQ**, **`--no-PG`** —
  CLI-accepted-and-ignored stubs.

**`import` deferred features**:

- **`--i1` / `--i2`** index-read inputs (the index-as-aux BC/QT shape).
  v1 wires `-0/-1/-2/-s` and the positional shapes; index files would
  attach a BC:Z and QT:Z tag computed from a separate index FASTQ. The
  parser scaffolding is in place; just needs a third walker.
- **`-i` CASAVA** parsing (extract barcode from CASAVA-style headers).
  We do parse the description tail for SAM aux fields directly, which
  covers the common case where the FASTQ was produced by `samtools
  fastq`; CASAVA-format input is a follow-up.
- **`--barcode-tag` / `--quality-tag`** renaming of the BC/QT tag pair.
  Not exposed in the v1 CLI (defaults to BC/QT).
- **`-O` / `--output-fmt`** and **`-@` / `--threads`** —
  CLI-accepted-and-ignored stubs. v1 picks output format from the
  output-path extension (`.sam` vs `.bam`) and is single-threaded.

**Validation:** small hand-built fixtures live under
`tools/samtools/testdata/parity/{calmd,import}/` covering the four
calmd code paths (match, mismatch, deletion, insertion+softclip) and
the six import shapes (-0, -1/-2, -s, single positional, two
positionals, -T aux extraction, --order, -R/-r RG). The upstream
`bam_md.c` / `bam_import.c` regression cases are marked as
`t.Skip(...)` parity stubs because upstream's BGZF output isn't
byte-identical with ours (different libdeflate). Logical correctness
is covered by hand-computed expected values in the table tests.

### `bcftools`

**Status:** 13 of ~30 subcommands (~43%). `view`, `index`, `stats`, `query`,
`concat`, `norm`, `call` (consensus + biallelic multi-allelic), and PR #86
wave-1 tail: `annotate`, `head`, `isec`, `merge`, `reheader`, `sort`.

Missing subcommands (priority order):

- **`mpileup`** — base-level pileup; required upstream input to
  `bcftools call`. Large port.
- **`csq`** — predict variant consequences against a GFF.
- **`filter`** — soft-filter records (different from view's hard-filter).
- **`consensus`** — apply variants to a FASTA.
- **`convert`** — between formats (GVCF / TSV / hapmap / PLINK).
- **`mendelian`** / **`mendelian2`** — Mendelian-inheritance checks.
- **`roh`** — runs of homozygosity.
- **`polysomy`** — copy-number estimation.
- **`cnv`** — CNV calling.
- **`gtcheck`** — genotype concordance.
- **`+plugins`** — full plugin system (substantial).

Option-tail gaps on the wave-1 additions (PR #86):

- `annotate --set-id '+%CHROM_%POS'` macro expansion is not implemented;
  `-x ID` / `-x INFO/TAG` / `-x FORMAT/TAG` removal works.
- `isec`: `--collapse some` (REF match + any-ALT-in-common) is approximated
  via strict tuple match; deeper semantics deferred.
- `merge`: pre-sort assumption is enforced; no automatic CHROM/POS sort.
- `reheader`: in-place rewrite (`-i`) currently emits to stdout — caller
  is responsible for the swap.
- `sort`: `-m/--max-mem` and `-T/--tmpdir` are accepted but always
  in-memory.

Plus:

- **`bcftools view`** more flags: `--regions-overlap`, `--targets-overlap`,
  `--no-version`, `--write-index`, `--phased`.
- **CSI seek** for region queries: today we validate via the index then
  linear-scan. Real chunk-seek is the natural follow-up.

Subcommand-tail gaps on `bcftools call`:

- **Full multi-allelic caller (`-m` on >2 ALT sites).** The v1 port
  falls back to the consensus model for sites with more than one ALT;
  upstream iterates over every allele combination.
- **BCF input.** Today `call` rejects BCF input with a roadmap-pointer
  error. The BCF reader's FORMAT-key reconstruction
  (`docs/UPSTREAM_BUGS.md`, `bcf-fmt-keys-missing`) is the prerequisite.
- **`--ploidy GRCh37 / GRCh38`.** Accepted by the CLI parser but
  rejected at runtime — the per-contig sex-chromosome overrides need a
  ploidy registry that's not yet wired in.
- **Index-backed region queries** (`-r` reuses the post-filter path).
- **`--gvcf` block-emit mode** (banded reference blocks).
- **`-C alleles --constrain`** family.

**Validation:** no upstream-test-suite run yet.

### `mosdepth`

**Status:** 1 / 1 command, most flags.

Missing:

- **`.csi`** output (currently emits `.tbi`).
- **D4 output** (`-d/--d4`).
- **Multi-threading** (`-t/--threads N`).
- **`--mapq` 0-only fast-path** — upstream has a special fast loop.

**Validation:** no upstream-test-suite run yet.

---

## Phase plan

### Phase 1 (this PR)

- Truth pass on `PORTING_STATUS.md`, `tools/README.md`, and `analysis/tool_ranking_2026.md`.
- Initial scaffold for `UPSTREAM_BUGS.md` and this file.

### Phase 2 — validated parity audits

For each tool that has an upstream test suite (samtools, bcftools, mosdepth,
vcftools, plus sickle and skewer which have small test corpora), set up
the same kind of validated-parity rig that bedtools got in PR #55:

- Submodule-init the upstream repo.
- Pull representative test cases (5-20 per subcommand) plus their
  expected outputs.
- Diff our output against expected; pass or `t.Skip()` with a documented
  reason.
- Fix small bugs we find inline; record upstream bugs in `UPSTREAM_BUGS.md`.

### Phase 3 — systematic gap closure

Tool-by-tool, dispatched as parallel agent waves of 3-4 PRs each:

1. **bedtools long tail** (~30 small subcommands, parallelisable).
2. **vcftools long tail** (~87 options, parallelisable).
3. **seqtk + fastp + sickle + skewer + prinseq closure** (small individual
   tools, all parallelisable).
4. **samtools mpileup + bcftools call + bcftools annotate** (the big
   variant-calling loop).
5. **`samtools markdup`/`idxstats`/`stats`/`merge`** (workhorse
   utilities).
6. **bcftools merge / isec / sort / annotate** (everyday set-ops).
7. **CRAM** (last; multi-week effort on its own).

Phase 1 lands first. Phases 2 and 3 are several sessions of parallel-agent
waves each.
