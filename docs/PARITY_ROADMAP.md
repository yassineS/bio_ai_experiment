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

**Status:** 13 of ~20 subcommands. ~65%.

Missing subcommands:

- `listhet` — extract heterozygous sites from VCF/BCF.
- `fqchk` — FASTQ quality check report.
- `seqshuf` — shuffle FASTA/Q records.
- `pair` — pair up R1/R2 from interleaved input.
- `dropse` — drop unpaired reads.
- `hpc-bg` — homopolymer-compress with mismatch tolerance.
- `kfreq` — k-mer frequency analysis.
- `gcdc` — GC depth count.

Note: `cnregion` was listed here before but is **not** an upstream seqtk
subcommand (verified against `reference_code/seqtk/seqtk.c` v1.5: the only
`stk_*` functions registered in `main()` are `seq`, `comp`, `sample`,
`subseq`, `mergefa`, `mutfa`, `mergepe`, `randbase`, `hety`, `gc`, `fqchk`,
`hrun`/`hpc`, `listhet`, `famask`, `trimfq`, `hpc-bg`/`hpc`, `seq`, `cutN`,
`gap`, and `kfreq` — no `cnregion`). Dropped from the gap list.

Option-tail gaps (per existing subcommand):

- `comp` — missing `-r REGION` to restrict to a BED region.
- `seq` — missing `-A` (force ASCII output), `-C` (mask sequence with N), `-M FILE` (mask regions), the `-T int` trim option.
- `sample` — missing `-2` (output two paired files).
- `trimfq` — missing `-L int` (max length cap), `-B int` (min base quality).
- `subseq` — missing the regex-name mode.
- `mutfa` — missing the inverse `--inverse` mode.
- `gap` — full upstream surface implemented (`-l` only). Note: upstream's
  "gap" is any non-ACGT byte (via `seq_nt6_table`), not just literal N,
  so IUPAC ambiguity codes (R, Y, S, W, K, M, B, D, H, V, N) all count.
  We match that byte-for-byte against `reference_code/seqtk` v1.5.
- `gc` — full upstream surface implemented (`-w`, `-f`, `-l`, `-x`). The
  algorithm is the upstream X-dropoff scan, not a sliding window. Output
  is BED4 (`chrom\tstart\tend\thits`). Byte-parity verified against
  `reference_code/seqtk` v1.5 on `gc_small.fa` fixtures.

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

**Status:** 24 of ~25 subcommands (~96%). `view`, `sort`, `index`, `depth`,
`fastq`, `flagstat`, **`mpileup`** (wave-1 + tail wiring), PR #88's
wave-1 tail (`merge`, `coverage`, `idxstats`, `cat`, `reheader`,
`addreplacerg`, `fixmate`, `dict`, `split`, `quickcheck`), the
heavy-hitter pair `markdup` + `stats`, the calmd/import pair
(**`calmd`** + **`import`**), the niche pair landed in the
phase/targetcut PR (**`phase`** + **`targetcut`**), and now
**`consensus`** (simple-mode FASTA/FASTQ/pileup; bayesian falls back
with a stderr warning).

Missing subcommands (in rough priority order):

- **`tview`** — terminal viewer. **Deliberate skip** (interactive
  curses UI; near-zero pipeline usage and would require an ncurses
  dependency). Not on the roadmap.
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

**`phase` deferred features** (accepted on the CLI, behaviour partial):

- **MCMC chimera repair**. Upstream's `phase.c` runs a
  Markov-chain-Monte-Carlo loop (`phase_core`) that flips read-cluster
  assignments to maximise haplotype consistency and resolve chimeric
  reads at junctions. The v1 Go port replaces this with a greedy
  same-vs-opposite vote between adjacent het sites. Tied junctions
  emit label `0` (ambiguous) rather than being repaired by MCMC.
  Tracked here; the upstream `FLAG_FIX_CHIMERA` flag is implicitly
  disabled in v1.
- **`-b STR` per-haplotype BAM split.** v1 emits the phased TSV
  stream to `-o`/stdout but does not yet split the input BAM into
  per-haplotype output BAMs (`<prefix>.0.bam` / `<prefix>.1.bam`
  / `<prefix>.chimera.bam` in upstream). The flag is accepted on the
  CLI and stored in `PhaseOptions.OutputPrefix` for a follow-up
  wiring pass.
- **`-F` use-full-read** is accepted on the CLI but is a no-op in v1
  (we always walk the aligned slice as decoded from the CIGAR).
- **`-A` mark-drop-in-chimera-output** is also a no-op pending the
  `-b` split landing.
- **`-e`/`-l` site-list mode** (only-phase-listed-sites). The
  upstream `loadpos` path is not implemented; the Go port always
  discovers hets from the pileup. Upstream itself comments `-e` and
  `-l` out of the usage block, so the omission is a small loss.

**`targetcut` scope reduction.** The user-facing spec for the Go port
is "cut the aligned slice from each read and emit FASTA". Upstream's
`cut_target.c` actually does something quite different —
HMM-based consensus calling over fosmid pools, emitting one consensus
SAM record per identified region. The HMM consensus mode is **not**
implemented; the upstream tool is rarely used outside fosmid
workflows. The simple aligned-slice FASTA mode landed here covers the
"cut a read down to its aligned bases" use case that users typically
mean when they reach for the name. The `-Q` flag is wired through to
the per-base quality filter as documented.

**Validation:** hand-built SAM fixtures in
`tools/samtools/pkg/samtools/phase_test.go` and
`tools/samtools/pkg/samtools/targetcut_test.go`. Phase tests cover
single-block chaining (consistent & label-flipping orderings),
ambiguous-label fall-back when reads don't bridge two hets, and the
MinMAPQ filter. Targetcut tests cover soft-clip flank stripping,
insertion retention, deletion handling, unmapped/secondary skipping,
SEQ='*' skipping, and `-Q` per-base filtering. There is no upstream
regression-test fixture for either tool in
`reference_code/samtools/test/` so byte-parity against upstream is
not pursued.

**`consensus` deferred features** (accepted as CLI flags, behaviour
partial). Upstream `bam_consensus.c` ships five modes — `simple` and
four bayesian flavours (`bayesian_r` aka "bayesian", `bayesian_m`,
`bayesian_p`, `bayesian_116`). v1 only implements `simple`. Because
upstream defaults to `MODE_RECALL` (a bayesian mode) at
`bam_consensus.c:2983`, the v1 binary's default invocation lands on
the bayesian branch, emits a single-line stderr warning, and falls
back to `simple`. The deferred surface:

- **Bayesian (Gap5-derived) mode.** All variants (`bayesian`,
  `bayesian_r`, `bayesian_m`, `bayesian_p`, `bayesian_116`) and their
  knobs (`-C/--cutoff`, `--P-het`, `--P-indel`, `--het-scale`,
  `--adj-qual`, `--use-MQ`, `--adj-MQ`, `--NM-halo`, `--SC-cost`,
  `--scale-MQ`, `--low-MQ`, `--high-MQ`, `-p/--homopoly-fix`,
  `--homopoly-score`, `--homopoly-redux`, `-t/--qual-calibration`,
  `-X/--config`) are accepted on the CLI but not yet implemented.
- **Pileup-mode insertion rows.** Upstream's default `--show-ins yes`
  emits extra rows with `nth>0` for each column of an inserted
  sequence. v1 emits only `nth=0` rows (one per reference position);
  insertion columns are folded into the FASTA/FASTQ stream when
  `--show-ins yes` (the default) but not into the pileup stream.
- **Mate-overlap dedup.** `--ignore-overlaps` is accepted but is a
  no-op; v1 counts each mate independently in the pileup walker.
- **Reference-aware modes.** `-T/--reference`, `--ref-qual`, and
  `--default-qual` are accepted but unused; the simple scoring path
  doesn't need a reference, and the bayesian path that does is
  deferred.
- **Threading.** `-@/--threads`, `-Z/--block-size`, and
  `--input-fmt-option` are accepted but ignored; v1 is single-pass
  and single-threaded.
- **Read-flag filtering.** `--rf/--incl-flags` and `--ff/--excl-flags`
  are accepted as text/int but ignored. v1's filter set is fixed
  (drop UNMAP|SECONDARY|QCFAIL|DUP, matching upstream's default
  `excl_flags`).
- **`--het-only`** suppression of homozygous calls is accepted but
  not implemented.
- **`--verbosity`** is accepted and ignored.

**`consensus` correctness model.** v1 mirrors upstream's
`calculate_consensus_simple` (`bam_consensus.c:1900-2006`) bit-for-bit
where it matters:

- One fraction gate only: `used_score < call_fract * tscore`
  (`bam_consensus.c:1988-1994`). There is **no** separate
  "min-fraction on the dominant base alone" gate; an earlier PR
  fabricated one and is corrected here.
- Heterozygous condition is `score2 >= het_fract * score1 && ambig`
  with no `score1 > 0` guard (`bam_consensus.c:1982`).
- `use_qual=0` by default (`bam_consensus.c:2984`) — bases score by
  frequency, not quality, until `-q/--use-qual` is set.
- `--show-del` is honoured in pileup mode too — rows whose call is
  `'*'` are suppressed when `--show-del no`
  (`bam_consensus.c:2244`).
- Insertion gating uses `MinCallFraction`, the same knob as the
  per-position gate.

**Validation:** table-driven hand-built SAM fixtures in
`tools/samtools/pkg/samtools/consensus_test.go` covering:
all-match FASTA/FASTQ/pileup, mixed at the 0.75 boundary,
`--show-del no` in pileup (no `*` rows),
`--show-del yes` in pileup (keeps `*` rows),
the canonical **50/30/20 + `-A`** fixture (must land `M`, not `N`),
frequency-only counting (low-Q vs high-Q gives identical output when
`UseQual=false`),
the `UseQual=true` flip (high-Q minority beats low-Q majority),
multi-contig, `-a` zero-fill, `--min-depth`, line-len wrapping,
insertion include/suppress + `--mark-ins`,
and the default-bayesian fallback emitting a stderr warning.
Coverage of the `pkg/samtools` package after this PR is ~80%.

### `bcftools`

**Status:** 24 of ~30 subcommands (~80%). `view`, `index`, `stats`, `query`,
`concat`, `norm`, `call` (consensus + biallelic multi-allelic), the PR #86
wave-1 tail (`annotate`, `head`, `isec`, `merge`, `reheader`, `sort`), the
convert/mendelian PR (`convert`, `mendelian`), the gtcheck/roh PR (`gtcheck`,
`roh`), the filter/consensus PR (`filter`, `consensus`), the
mendelian2/polysomy PR (`mendelian2`, `polysomy`), the cnv/csq PR
(`cnv` + `csq`), and the mpileup PR (**`mpileup`**).

Missing subcommands (priority order):

- **`+plugins`** — full plugin system (substantial).

Option-tail gaps on `gtcheck` (PR #107, simple-mode):

- `--cluster N,N` (HMM-style sample clustering), `--distinctive-sites`,
  `--n-matches` — accepted-and-rejected with PARITY_ROADMAP pointer;
  bayesian-mode follow-up.
- `-u PL` — PL/GL-based scoring; v1 only does hard-GT Hamming.
- `-O z` — bgzip output; v1 only emits tab-text (`-O t`).
- `[5]Average -log P(HWE)` column is zeroed until a real per-site HWE
  estimator from panel AF lands.
- Index-backed `-r/-R` seek (post-filter only in v1).
- Multi-allelic input is rejected (matches upstream's
  `bcftools norm -m -` requirement).

Option-tail gaps on `roh` (PR #107, simple-mode):

- `-b/--buffer-size`, `-e/--estimate-AF`, `-m/--genetic-map`,
  `-M/--rec-rate`, `-V/--viterbi-training` — accepted-and-rejected
  with PARITY_ROADMAP pointer.
- `-O z` — bgzip output; v1 only emits tab-text.
- Transition defaults are upstream's literal per-bp magnitudes
  (`6.7e-8` / `5e-9`) but NOT scaled by physical inter-marker
  distance, so RG quality scores are NOT comparable to upstream's
  until distance scaling lands.

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

Option-tail gaps on the convert/mendelian PR:

- `convert`: v1 covers only the pass-through round-trip
  (VCF↔BCF↔VCF.gz) with sample/region filtering and -i/-e expressions.
  The full upstream `vcfconvert.c` covers many extra shapes
  (`--gvcf2vcf`, `--haplegendsample2vcf`, `--hapsample2vcf`,
  `--tsv2vcf`, `--gensample2vcf`, `--gvcf`, PLINK / GEN / HAP).
  These are explicit follow-ups; the CLI emits a usage block that
  lists them under "Deferred output paths".
- `mendelian`: the v1 port detects Mendelian inconsistencies for
  PED-style trios (one or more `-t CHILD,FATHER,MOTHER` flags, or
  `-T trio-file`), emits `INFO/MERR` per record, and supports the
  five upstream modes `a|c|x|d|+`. The `+` mode is treated as a
  synonym for annotate (the upstream-only `##bcftools_PG=` provenance
  header is deferred). The `--rules FILE` ploidy specification is
  accepted but currently only the chrX heuristic is honoured; full
  per-contig ploidy override (PAR boundaries, mitochondrial
  haploidy) is a follow-up.

Option-tail gaps on `filter` (this PR, simple-mode):

- `--mask [^]REGION` and `-M/--mask-file [^]FILE` — accepted but
  hard-rejected at runtime with a roadmap pointer. The CLI parses the
  flag (so a downstream automation that always passes `-M ""` doesn't
  break); the underlying BED-driven soft-filter logic is deferred.
- `--mask-overlap 0|1|2` — accepted; v1 ignores (always treats POS-in-region).
- `-W/--write-index[=FMT]` — accepted; v1 never auto-indexes outputs.
- `-v/--verbosity INT` — accepted; v1 ignores.
- `--regions-overlap` / `--targets-overlap` — accepted; v1 always
  uses POS-in-region semantics.
- `-g/--SnpGap INT:TYPE` — the `:TYPE` qualifier (indel|mnp|bnd|other|overlap)
  is parsed but always treated as "indel" in v1.
- BCF output (`-O b|u`) round-trips through the shared `pkg/bioformats/bcf`
  writer; CSI auto-indexing is the `-W` follow-up above.

Option-tail gaps on `mendelian2` (PR #109, simple-mode):

- `--rules ASSEMBLY` — predefined inheritance rules (GRCh37 / GRCh38
  / `list?`). Accepted but rejected at runtime with a roadmap pointer;
  v1 uses the chrX heuristic from the legacy `mendelian` port instead.
- `--rules-file FILE` — custom inheritance rules file. Accepted; v1
  rejects at runtime.
- `-W/--write-index[=FMT]` — auto-index output. Accepted; v1 never
  auto-indexes.
- `--regions-overlap 0|1|2`, `--targets-overlap 0|1|2` — accepted;
  v1 always uses POS-in-region semantics.
- `-v/--verbosity INT`, `--no-version` — accepted; v1 ignores both.
- `-r/-R/-t/-T` — region / target post-filter is wired through the
  CLI but the BCF synced-reader region-jump path is not used in v1;
  filtering happens after the records are read.
- `sites_not_diploid` counter never goes up in v1 because our
  `vcf.Variant` decoder coerces non-diploid GTs to missing rather
  than tracking ploidy. Tracked alongside the broader BCF FORMAT
  reconstruction in `docs/UPSTREAM_BUGS.md`.
- The PED-row sort order is by child name (deterministic);
  upstream sorts by min(sample-index) for sequential VCF reads. The
  sort is a performance optimisation only — the set of reported
  trios and per-trio counters are identical.
- `-i/-e EXPR` are accepted at the library boundary but only
  applied as record-level filters in v1 (no per-sample mask). The
  `sites_fail` counter is therefore always 0 in v1.

Option-tail gaps on `polysomy` (this PR, simple-mode):

- The full Gaussian-mixture peak-fit (`peakfit.c` + GSL) is
  **deferred**. v1 emits CN calls from a median-deviation heuristic:
  CN1 when n_het == 0, CN2 when |median(BAF) - 0.5| ≤ MinBafDev (with
  CnPenalty scaling), CN3 otherwise. The upstream algorithm fits CN2,
  CN3, and CN4 Gaussian mixtures and picks the lowest CN whose fit
  beats `(1 - cn_penalty) * previous_fit`.
- `-b/--peak-size FLOAT`, `-f/--fit-th FLOAT`, `-i/--include-aa`,
  `-p/--peak-symmetry FLOAT`, `--nbins INT`, `--smooth INT`,
  `--ra-rr-scaling` — accepted at the CLI for parity but inert in v1
  (the heuristic doesn't run a peak fit).
- `-o/--output-dir PATH` — accepted but ignored; v1 writes the
  per-chromosome TSV to stdout (no per-chromosome PNG plots).
- `--regions-overlap`, `--targets-overlap` — accepted; v1 always uses
  POS-in-region.
- `-v/--verbosity INT` — accepted; v1 ignores.
- `-n/--include-noise` — accepted; v1 always emits every chromosome
  (we don't classify any as noise / `?`).
- `--force-cn INT` (hidden upstream option) — implemented as
  per-chromosome override.
- BAF source: upstream requires FORMAT/BAF; v1 also accepts
  FORMAT/AD = REF,ALT as a fallback (synthesises BAF as
  `ALT / (REF + ALT)` at het sites).
- Per-record `-i/-e` are NOT in upstream `polysomy.c:main_polysomy`
  and we follow upstream's surface exactly (no invented flags).

Option-tail gaps on `cnv` (this PR, v1 heuristic):

- **The v1 algorithm is NOT the upstream HMM.** Upstream's vcfcnv.c
  runs a 5-state HMM (CN0/CN1/CN2/CN3/CN4) over each contig with
  joint BAF + LRR Gaussian emissions and a configurable transition
  matrix. The v1 port replaces this with a per-sample × per-chrom
  median-BAF + mean-LRR heuristic that classifies each chromosome
  into one of the same 5 CN states. The full Viterbi sweep is the
  natural follow-up; the CLI surface is already parity-clean for
  it. EVERY HMM tuning knob (`-a/--aberrant`, `-b/--BAF-weight`,
  `-e/--err-prob`, `-l/--LRR-weight`, `-L/--LRR-smooth-win`,
  `-O/--optimize`, `-P/--same-prob`, `-W/--baum-welch`,
  `-x/--xy-prob`, `--AF-file`) is parsed and stored in `CNVOptions`
  but the heuristic does NOT consume them. Only `-d/--BAF-dev` and
  `-k/--LRR-dev` (the per-sample expected std-dev floors) actually
  drive the v1 thresholds.
- `-o/--output-dir` — upstream writes per-sample / per-region plot
  data into this directory; v1 always streams a single summary TSV
  to stdout regardless of the path (the flag is still required for
  CLI parity).
- `-p/--plot-threshold` — accepted; v1 emits no plots.
- `--regions-overlap` / `--targets-overlap` — accepted; v1 always
  uses POS-in-region semantics.
- `-v/--verbosity` — accepted; v1 ignores.
- BCF / VCF.gz output — v1 always emits the summary TSV; the
  upstream `-O b|u|z|v` selector does not apply (the upstream tool
  produces several per-region files; v1 produces one summary).
- Indel / non-SNP records — the BAF/LRR signals are typically
  per-marker SNP data; v1 honours upstream's behaviour (treat each
  record as one marker regardless of REF/ALT).

Option-tail gaps on `csq` (this PR, v1 SNP-only):

- **The v1 classifier is NOT haplotype-aware.** Upstream's csq.c
  phases variants per haplotype, walks the GFF transcripts, and
  reports the per-haplotype consequence chain (including
  compound-het effects). v1 instead classifies one SNP at a time
  against the GFF CDS exons and emits `INFO/BCSQ` per-transcript
  for the matching position. Indels, splice-site disruption,
  start-gain, stop-retained, and compound-het bookkeeping are all
  deferred.
- `-p/--phase a|m|r|R|s` — parsed and stored; the per-record SNP
  classifier ignores phasing because consequences are computed
  position-by-position. Will become load-bearing when haplotype-
  aware phasing lands.
- `-i/--include` / `-e/--exclude` — accepted; v1 ignores (every
  input record runs through the classifier). The expression
  evaluator already exists in `pkg/bcftools`; the wire-up is a
  trivial follow-up.
- `-s/--samples` / `-S/--samples-file` — accepted; v1 does not
  subset (consequences are position-driven). Once haplotype-aware
  phasing lands the sample list will gate which haplotypes are
  walked.
- `-n/--ncsq INT` — accepted; v1 emits every matching transcript
  without a cap (one BCSQ entry per transcript).
- `-B/--trim-protein-seq INT` — accepted; v1 does not truncate
  predictions.
- `-b/--brief-predictions` — hard-rejected with a roadmap pointer
  (upstream deprecates this flag itself).
- `-C/--genetic-code INT|l` — only `0` (standard) is accepted in
  v1; other tables are hard-rejected with a roadmap pointer.
- `-l/--local-csq` — accepted; v1 always operates per-record.
- `--unify-chr-names LIST` — only `0` (no rewriting) is accepted in
  v1; non-zero specs are hard-rejected.
- `--dump-gff` — hard-rejected; v1 has no GFF dump path.
- `-O b|u|z|t` — only `-O v` (VCF text) is supported in v1; the
  others are hard-rejected with a roadmap pointer.
- `--threads`, `-v/--verbosity`, `-W/--write-index`, `--force`,
  `--no-version`, `-q/--quiet` — accepted; v1 ignores.
- The minimal GFF3 parser (`pkg/bioformats/gff`) understands `gene`,
  `mRNA` / `transcript`, `CDS`, and `exon` rows. Other feature
  types are silently skipped — fine for the v1 SNP classifier
  but the parser will need extension for splice-site / UTR work.

Option-tail gaps on `mpileup` (this PR, v1 SNP + uniform-error):

- **The v1 likelihood model is NOT the upstream MAQ model.** Upstream's
  `bam2bcf.c::glfgen` reads per-base error probabilities from the
  Heng Li MAQ recalibrator with BAQ adjustments; v1 instead uses the
  simpler samtools-0.1.19-style uniform-error binomial: e = 10^(-Q/10)
  per base, summed in log10 across reads, then phred-scaled and
  rebased to min=0 for the [0/0, 0/1, 1/1] triple. The MAQ port is
  the natural follow-up; the CLI surface and FORMAT/PL layout are
  parity-clean for it.
- **No BAQ recalibration.** `-B/--no-BAQ` is the v1 default (the flag
  is accepted as a no-op); `-D/--full-BAQ` is accepted but inert;
  `-E/--redo-BAQ` is hard-rejected with a roadmap pointer because
  silently skipping a recalibration step a downstream caller asked
  for would yield misleading PLs.
- **No indel calling.** The full upstream indel realigner
  (`bam2bcf_indel.c`) and the consensus indel mode
  (`bam2bcf_edlib.c`) are deferred. Every knob that drives the indel
  model — `-e/--ext-prob`, `-F/--gap-frac`, `-h/--tandem-qual`,
  `--indel-bias`, `--indel-size`, `-I/--skip-indels`,
  `-L/--max-idepth`, `-m/--min-ireads`, `-M/--max-read-len`,
  `--open-prob`, `--indels-cns`, `--indels-2.0`, `--no-indels-cns`,
  `--ar-prob`, `--ambig-reads / --ar`, `--del-bias`, `--poly-mqual`,
  `--no-poly-mqual`, `--score-vs-ref`, `--seqq-offset` — is accepted
  at the CLI for parity but inert in v1. The v1 emit path is
  equivalent to running upstream with `-I/--skip-indels` set.
- **No multi-allelic FORMAT/PL grid.** The PL we emit is always
  biallelic [PL(0/0), PL(0/1), PL(1/1)] against ALT[0]. Sites with
  multiple ALTs still parse upstream-style in `bcftools call`, but
  the PL grid for the 2nd / 3rd ALT genotypes is treated as 0
  (uninformative). The full j(j+1)/2 + i grid lands with the MAQ
  port.
- **No BAI seek.** `-r/--regions` and `-R/--regions-file` are
  post-filters applied after a linear scan of every input BAM; the
  BAI-seek fast path lives in `pkg/bioformats/sam` but is not wired
  through `mpileup` in v1. Tracked as a follow-up — perf only, no
  output difference.
- **No per-read group filtering.** `-G/--read-groups` is parsed and
  stored; v1 includes every record whose @RG passes the standard
  filters. `-Z/--ignore-RG` is accepted but inert.
- **No gVCF blocking.** `-g/--gvcf` is accepted; v1 always emits one
  VCF record per variant site (REF-only sites are skipped, matching
  upstream when `--gvcf` is unset).
- **`-a/--annotate LIST` is accepted but inert.** v1 always emits the
  default `INFO/DP`, `INFO/I16`, `FORMAT/PL` set. The
  `INFO/AD,ADF,ADR,SP,SCR,IDV,IMF`, `FORMAT/AD,ADF,ADR,DP,DV,DPR,SP,SCR,QS`
  tags will land alongside the per-tag stream when called from
  `bcftools call`.
- **`-O u|b` (BCF output) is hard-rejected.** The BCF writer in
  `pkg/bioformats/bcf` can handle generic records, but mpileup carries
  custom INFO/I16 typing rules; the wire-up is a follow-up. `-O v`
  (text VCF) is the default; `-O z` (gzipped VCF) is accepted at the
  CLI but currently streams text — gzip-wrap-stdout from a follow-up
  CLI shim will close that gap.
- **`--threads`, `-v/--verbosity`, `-W/--write-index`, `--no-version`,
  `-A/--count-orphans`, `-x/--ignore-overlaps`, `-d/--max-depth`,
  `-q/--min-MQ`, `-Q/--min-BQ`, `--max-bq`** — fully implemented in v1.
- `-X/--config STR` (presets like `1.12`, `2.1`, `ultima`,
  `pacbio-ccs-1.20`) — accepted; v1 ignores. Most presets toggle the
  indel-model knobs which are inert above.
- `-6/--illumina1.3+` — accepted; v1 ignores (input BAMs are
  Phred+33 across the board).
- `-C/--adjust-MQ INT` (MAPQ tail adjustment) — accepted; v1 ignores.
- The 5-flag mask flags (`--skip-any-unset`, `--skip-all-unset`,
  `--skip-any-set`, `--skip-all-set`, `--ls`) — accepted and stored;
  v1 honours only the standard `--ff` defaults (UNMAP, SECONDARY,
  QCFAIL, DUP, SUPPLEMENTARY) baked into `mpileupKeepRecord`.
- `--seed`, `--delta-BQ` — accepted; v1 ignores.

Option-tail gaps on `consensus` (this PR, simple-mode):

- `-c/--chain FILE` — liftover chain file. Accepted; v1 rejects with a
  roadmap pointer at runtime.
- `-H NpIu` — phased-index / unphased-IUPAC encoding. Accepted; v1
  rejects with a roadmap pointer.
- `--regions-overlap 0|1|2` — accepted; v1 ignores (no synced-reader
  region jump path).
- `-v/--verbosity INT` — accepted; v1 ignores.
- `-r/-R/-t/-T` — upstream `consensus` does NOT advertise these flags
  (only `--regions-overlap`); the v1 port follows upstream exactly.
- Multi-sample apply: `-s LIST` accepts a comma list but v1 honours only
  the first name (the other entries are silently ignored, mirroring
  upstream's single-sample focus when -H is unset).
- Complex MNP and SV (BND, breakend, `<DEL>` etc.) records are not yet
  applied to the reference. The current v1 covers SNPs and simple
  REF/ALT length-difference indels.
- Overlapping variants: first wins (left-to-right). Upstream emits a
  warning and folds them into a longer ALT; v1 doesn't.

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
