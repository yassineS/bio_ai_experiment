# Project status

**Last refreshed:** 2026-06-11 (post the ~70-PR wave: `call` modes,
`convert` GEN/HAP/TSV modes, `annotate`, `consensus` chain, `mendelian2`
rules, `csq` slices 1–4, mpileup legacy `bam2bcf_indel` + `--indels-cns`,
`phase`, `gtcheck`, view/sort/markdup threading, samtools `coverage -A` /
`markdup -d/-s/-S` / `calmd -C/-e/-u` / consensus indel calling, bcftools
`query %INFO/%SAMPLE` / `roh -Oz` / `filter -M` / `cnv --AF-file`,
mosdepth `--fragment-mode/--quantize/-t/--use-median`, bedtools BAM/VCF/GFF
input + `bedclosest` direction flags + `intersect -c`, CRAM bzip2 encode)

This is the honest "distance to done" snapshot for the started tools. It is
**the summary source of truth**; the authoritative gap list per tool is
[`docs/PARITY_ROADMAP.md`](docs/PARITY_ROADMAP.md), and the docs map in
[`docs/README.md`](docs/README.md) says which document owns what. This file
is the skimmable summary on top of the roadmap — keep status numbers here and
in the roadmap, and link to them from everywhere else rather than re-stating.

## What "done" means

A tool is **done** (1:1) when every upstream subcommand and documented flag
is present, every supported input produces the same logical result (modulo
documented intentional deviations, fixed-on-port upstream bugs, and the
RNG-byte-parity carve-out), and the validated-parity suite passes with an
explicit `t.Skip()` for each documented exception. "Working subset" is a
milestone, not the finish line. The percentages below are an honest,
evidence-based estimate of remaining surface — not a rosy reading.

## Per-tool completion

| Tool | Implemented | Remaining surface | % | Effort to finish |
|------|-------------|-------------------|---:|------------------|
| **seqtk** | All 24 subcommands; byte-parity vs v1.5 | none | **100%** | done |
| **prinseq** | `stats` / `filter` / `graph_data`; all in-scope flags + `--graph_data` | ~16 niche knobs (`--range_len`, `--trim_to_len`, `--seq_case`, `--line_width`, …); PNG report flow (out of scope) | **~95%** | small |
| **sickle** | `se` / `pe`; 15/15 parity cases | none (gzip-output level untested) | **100%** | done |
| **skewer** | `se` / `pe`; 14/14 parity cases | none | **100%** | done |
| **fastp** | single cmd; sliding-window, auto-adapter, HTML+JSON, dup-eval, UMI, PE base `--correction`, overrepresentation (`-p/-P`), `--split*`, merge writer (`-m`), `--adapter_fasta`, `--poly_x_min_len`, `--disable_adapter_trimming`; 16/16 + 9 tail parity | multi-thread `--split` file-boundary distribution; `merged_and_filtered` JSON block (merged FASTQ bytes are byte-identical) | **~97%** | small |
| **bedtools** | 37 bed* tools; no missing subcommands; 141+ parity tests; BAM input (`bedintersect`/`bedmulticov`), VCF/GFF input (`bedintersect`/`bedmultiinter`), `bedclosest` direction flags, `intersect -c` | scattered option-tail polish; CRAM input deferred | **~95%** | small (long tail) |
| **vcftools** | single cmd; **146/146 upstream long flags (100%)**, incl. BCF I/O, PCA, LD, RoH, relatedness, `--freq2`/`--counts2` (schema complete) | per-output column-set polish only; `--max-indv` uses deterministic truncation not upstream RNG shuffle (RNG-policy non-goal) | **~98%** | small |
| **bcftools** | 24 subcommands (all present); mpileup MAQ SNP model (slices 1–4) + legacy `bam2bcf_indel` + `--indels-cns` (edlib realigner) + BAQ + bias tags; full multi-allelic `call` (`-m`/`-c`/`--gvcf`/`-C alleles`/`-G`/`--ploidy GRCh37/38`/`--ploidy-file`); `convert` GEN/HAP/TSV/gVCF modes; `gtcheck`/`mendelian2`/`consensus` (chain+iupac)/`annotate`; `filter -M`/`cnv --AF-file`/`roh -Oz`/`query %INFO/%SAMPLE`; csq slices 1–4 (FORMAT/TBCSQ, --unify-chr-names, -O b\|u\|z, --dump-gff); full HMM `roh`/`cnv`/`polysomy`; subprocess plugin system; `csq -l/--local-csq` (test_cds_local) | `gtcheck -c/--cluster` + filter exprs; `convert` PLINK exporters (phantom — commented out upstream); `query %N_ALT` (non-goal); `som`/`tview` (non-goals) | **~97%** | small |
| **samtools** | 24 functional subcommands; CRAM r/w + bzip2 encode done; `.csi` done; consensus `--het-only` + indel calling; `coverage -A`; `markdup -d/-s/-S`; `calmd -C/-e/-u`; `phase` (full upstream-schema emit); mpileup MAQ **BCF/VCF emit (slices 1–4) + legacy indel + --indels-cns** done; `-@` threading **done for view/sort/markdup** | consensus pileup `-a` placeholder rows (deletion/zero-cov); `-@` parallel *input* BGZF/CRAM decode + non-BAM-writing subcommands; `tview` (deliberate skip) | **~97%** | small |
| **mosdepth** | single cmd; ALL flags wired incl. `-d/--d4`, `--fragment-mode`, `--quantize`, `-t/--threads`, `--use-median`, `--mapq` fast-path; emits `.csi` | CRAM input (`-f/--fasta` accepted but BAM-only decode) | **~95%** | small |
| **bgzip** | 1/1 cmd, most flags | multi-threaded compression (`-t`) | **~92%** | small |
| **tabix** | 1/1 cmd, most flags | `--reheader`; true `--targets` post-filter strictness | **~92%** | small |
| **htsfile** | 1/1 cmd | none (intentional `-c` omission) | **~98%** | done |
| **CRAM / htsgo formats** | CRAM read+write (rANS 4x8/4x16 in-tree, LZMA via ulikunitz/xz); BGZF, BAI/CSI, SAM/BAM, BCF v2.2, tabix index | partial-decompress seek via `.gzi` (perf); BCF FORMAT-key reconstruction edge cases | **~90%** | medium |

Percentages weight *remaining feature surface and effort*, not flag counts
alone — e.g. vcftools is 100% on flags but ~98% because a few outputs still
need column polish. The former bcftools/samtools "boulders" (full
multi-allelic `call`, mpileup SNP slices 1–4, **both** the legacy
`bam2bcf_indel` and the `--indels-cns` edlib indel callers, the csq
haplotype engine through slice 4, `convert` GEN/HAP/TSV modes, the
roh/cnv/polysomy HMMs) have all landed and are live-oracle validated.

## Genuinely-remaining real gaps

The list below is the *definitive* remaining-gap set (each is small and
individually scoped). Everything not on this list is either done or a
documented **non-goal** (see "Non-goals" below).

1. **htsgo `hfile` cloud I/O** — remote (S3 / http(s)) alignment/variant
   file access is not implemented (deferred P3 in `docs/HTSGO_ROADMAP.md`).
   Local files, BGZF, and `.gzi`/`.crai`/`.csi` seek all work. *Medium.*
2. **bcftools `gtcheck -c/--cluster`** (dendrogram, which upstream itself
   errors "to be implemented") and `gtcheck` filter expressions. *Small.*
3. **CRAM long-tail correctness + perf** — some BCF FORMAT-key
   reconstruction edge cases and the network REF_PATH/EBI reference fetch
   (an unresolvable reference is surfaced as a clear MD5 error); CRAM v4.0
   awaits a final spec. *Medium.*
4. **Scattered option-tail polish** — vcftools per-output column sets;
   prinseq niche knobs; a handful of bedtools per-subcommand flag tails.
   Individually *small*.

Cross-cutting: **multi-threading (`-@`/`-t`)** has landed for the dominant
BAM-output path — `samtools view`/`sort`/`markdup` and mosdepth `-t` drive a
parallel BGZF compressor / decompressor (decode-equal across thread counts)
— and remains a deferred no-op for parallel *input* BGZF/CRAM decode in
several non-BAM-writing subcommands and bgzip. It's a single deferred
parallel-read pass, not per-tool feature gaps.

## Non-goals (not gaps)

These are deliberately not ported and should not be counted against parity:

- **`bcftools som`** — the upstream source has an `fwrite`-return bug that
  reproduces as a crash when ported; documented in
  [`docs/UPSTREAM_BUGS.md`](docs/UPSTREAM_BUGS.md).
- **`bcftools`/`samtools tview`** — interactive ncurses viewer, no pipeline
  use, would require an ncurses dependency.
- **`query %N_ALT`** and **`import --skipBamQ`** — **not** upstream flags.
- **RNG byte-parity** for `seqtk sample`/`randbase`, `vcftools --max-indv`,
  `bedshuffle`/`bedsample` — policy is structural invariants + within-tool
  reproducibility, **not** byte-identity with upstream's C RNG (see the RNG
  policy section in `docs/PARITY_ROADMAP.md`). (`seqtk sample` and
  `bedsample` additionally now port the upstream RNG and *are* byte-exact.)
- **prinseq PNG report** (`prinseq-graphs.pl` graphics flow) — out of scope;
  the `graph_data`/`report` subcommands cover the equivalent surface.
- **The ~30 bundled upstream bcftools `.so` plugins** — the plugin *system*
  (a VCF-on-stdin/stdout subprocess protocol) is implemented; re-porting
  upstream's plugin catalogue is explicit non-goal scope.

## Where to look next

| Question | Doc |
|---|---|
| Which doc owns which info? | [`docs/README.md`](docs/README.md) (the docs map) |
| Authoritative per-tool gap list? | [`docs/PARITY_ROADMAP.md`](docs/PARITY_ROADMAP.md) |
| Per-subcommand feature lists? | [`tools/PORTING_STATUS.md`](tools/PORTING_STATUS.md) |
| Upstream bugs we fixed on port? | [`docs/UPSTREAM_BUGS.md`](docs/UPSTREAM_BUGS.md) |
| Which tools to port next? | [`analysis/tool_ranking_2026.md`](analysis/tool_ranking_2026.md) |
| How is the repo organised? | [`CLAUDE.md`](CLAUDE.md) |
| How do I use tool X? | `tools/<tool>/README.md` |
