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
   bytes, modulo:
   - documented intentional deviations (recorded in the tool's README +
     this roadmap),
   - upstream bugs we chose to fix (recorded in [`UPSTREAM_BUGS.md`](UPSTREAM_BUGS.md)).
4. The validated-parity test suite (runs the upstream test corpus through
   our port) passes for every supported case, with explicit `t.Skip()` for
   each documented exception.

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

### bedtools (20 subcommands ported)

**Status:** 20 of ~40 subcommands (~50%). 135 passing parity tests against
upstream test suite (across PR #55 + Phase-3 wave 1 + wave 2 simple + wave 2
algo), 53 documented `t.Skip`, 2 known discrepancies. Phase-3 wave 1
(PR #78) added `bedgroupby`/`bed12tobed6`/`bedmakewindows`; wave 2 simple
(PR #80) added `bedexpand`/`bedgetfasta`/`bedsample`/`bedspacing`; wave 2
algo (this PR) adds `bedcoverage`/`bedmap`/`bedshuffle`.

Missing subcommands:

- `multicov` — multi-sample coverage.
- `multiinter` — multi-way intersection.
- `reldist` — relative distance distribution.
- `fisher` — Fisher's exact for overlap.
- `nuc` — nucleotide content per interval.
- `window` — overlap A with a window around B.
- `split` — split into approximately equal-sized files.
- `summary` — per-chrom summary.
- `tag` — annotate A with B's name.
- `igv` — generate IGV launch URLs.
- `links` — generate UCSC links.
- `cluster` — cluster overlapping intervals.
- `pairtopair`, `pairtobed` — paired BEDPE operations.
- `annotate` — annotate one BED with multiple BEDs.

Skipped parity cases from PR #55 to revisit:

- `bedjaccard` auto-merge behaviour.
- `bedmerge` `.`-strand fan-out.

Option-tail gaps on the wave-2 additions:

- `bedgetfasta` — `-fullHeader` (whitespace-aware contig name parsing) and
  BGZF FASTA input via `.gzi` are not yet implemented. `pkg/bioformats/fasta`
  needs a `.gzi` reader before BGZF random-access fetch is feasible.
- `bedsample` — output PRNG is Go `math/rand` and is not byte-compatible
  with upstream's C++ sampler. Seeded runs are deterministic within
  `bedsample` (same seed → same output) but cross-tool record-for-record
  parity with upstream is not feasible without porting upstream's
  `random_shuffle`.

Column-op gaps (shared between `bedmerge`, `bedgroupby`, and the future
`bedmap`/`bedcoverage`):

- `stdev` / `sstdev` — standard deviation / sample stdev.
- `absmin` / `absmax` — min/max by absolute value.
- `cat` / `cat_uniq` — concatenate / concatenate-unique.

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

**Status:** 6 of ~25 subcommands (~25%). `view`, `sort`, `index`, `depth`,
`fastq`, `flagstat`.

Missing subcommands (in rough priority order for variant-calling
workflows):

- **`mpileup`** — base-level pileup; required for any caller. Large port.
- **`markdup`** — mark/remove PCR duplicates.
- **`idxstats`** — per-reference read counts.
- **`stats`** — exhaustive per-file statistics (different from `flagstat`).
- **`merge`** — combine sorted BAMs.
- **`coverage`** — region-level coverage summary (different from `depth`).
- **`addreplacerg`** — add/replace read group.
- **`fixmate`** — fill in mate-read coordinates.
- **`calmd`** — compute MD/NM tags.
- **`reheader`** — replace header in place.
- **`cat`** — concatenate sorted BAMs without re-sorting.
- **`quickcheck`** — fast format sanity check.
- **`dict`** — emit a sequence dictionary.
- **`split`** — split by read group.
- **`consensus`** — base-level consensus.
- **`import`** — convert FASTQ/SAM to BAM.
- **`phase`** — phase reads with their mates.
- **`targetcut`** — cut targeted regions.
- **`tview`** — terminal viewer (likely a deliberate skip; ~no-one uses it).
- **`view` flag-tail**: `-L bed`, `-M` (read-id list), `-d/-D` (tag-value
  filter), `-N` (qname file), `-X` (custom-index input).

Plus:

- **CRAM** read/write throughout. Big.
- **`.csi`** for samtools (BAI is fine for chromosomes ≤512Mb).
- **Multi-threading** in `sort`, `index`, `view` (`-@`).

**Validation:** no upstream-test-suite run yet.

### `bcftools`

**Status:** 7 of ~30 subcommands (~23%). `view`, `index`, `stats`, `query`,
`concat`, `norm`, `call` (consensus + biallelic multi-allelic).

Missing subcommands (priority order):

- **`mpileup`** — base-level pileup; required upstream input to
  `bcftools call`. Large port.
- **`annotate`** — annotate records from a tab-indexed table.
- **`csq`** — predict variant consequences against a GFF.
- **`isec`** — set operations on VCF/BCF files.
- **`merge`** — combine VCFs from different samples.
- **`filter`** — soft-filter records (different from view's hard-filter).
- **`sort`** — sort VCF/BCF.
- **`reheader`** — replace header.
- **`consensus`** — apply variants to a FASTA.
- **`head`** — emit just the header (different from view's `-h`).
- **`convert`** — between formats (GVCF / TSV / hapmap / PLINK).
- **`mendelian`** / **`mendelian2`** — Mendelian-inheritance checks.
- **`roh`** — runs of homozygosity.
- **`polysomy`** — copy-number estimation.
- **`cnv`** — CNV calling.
- **`gtcheck`** — genotype concordance.
- **`+plugins`** — full plugin system (substantial).

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
