# Project status

**Last refreshed:** 2026-06-10 (post PRs #220–#252: `call` modes, `convert`,
`annotate`, `consensus`, `mendelian2`, `csq`, `mpileup --indels-cns`,
`phase`, `gtcheck`, and view/sort/markdup threading)

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
| **fastp** | single cmd; sliding-window, auto-adapter, HTML+JSON, dup-eval, UMI; 16/16 parity | overrepresentation analysis, PE base `--correction`, overlap-trim knobs, `--split*`, `--adapter_fasta`, JSON sub-fields | **~85%** | medium |
| **bedtools** | 35 of ~40 subcommands; no missing subcommands; 141+ parity tests | option-tail polish only (`bedmultiinter` VCF/GFF input, `bedsample` RNG byte-parity, scattered flag tails); CRAM input deferred | **~90%** | medium (long tail) |
| **vcftools** | single cmd; **146/146 upstream long flags (100%)**, incl. BCF I/O, PCA, LD, RoH, relatedness | per-output column-set polish only; `--max-indv` uses deterministic truncation not upstream shuffle | **~97%** | small |
| **bcftools** | 24 subcommands (all present); mpileup MAQ SNP model + `--indels-cns` (edlib realigner) + BAQ + bias tags; `call` modes (`-m`/`-c`/`--gvcf`/`-C alleles`/`-G`/`--ploidy`); `convert`/`gtcheck`/`mendelian2`/`consensus`/`annotate`; csq slices 1–3; full HMM `roh`/`cnv`/`polysomy`; plugin system (subprocess design) | **classic `bam2bcf_indel` indel caller** (the `--indels-cns` CNS model landed; the legacy non-CNS indel path is deferred); **convert's ~18 extra modes**; **csq slice 4** (FORMAT/TBCSQ, --unify-chr-names, -O non-text); `gtcheck` PL/cluster tail; assorted `-O z`/`-W`/region-seek tails | **~70%** | large |
| **samtools** | 24 functional subcommands; CRAM r/w done; `.csi` done; consensus `--het-only` done; `phase` (full upstream-schema emit); mpileup MAQ genotype-likelihood **BCF/VCF emit (slices 1–4) done**; multi-threading (`-@`) **done for view/sort/markdup** BAM output | classic indel calling (`bam2bcf_indel`, shared deferral with bcftools); `-@` parallel *input* BGZF/CRAM decode + non-BAM-writing subcommands; `tview` (deliberate skip) | **~88%** | medium–large |
| **mosdepth** | single cmd; most flags; `-d/--d4` done (byte-identical) | `.csi` output (emits `.tbi`); multi-threading (`-t`); `--mapq 0` fast-path | **~85%** | small–medium |
| **bgzip** | 1/1 cmd, most flags | multi-threaded compression (`-t`) | **~92%** | small |
| **tabix** | 1/1 cmd, most flags | `--reheader`; true `--targets` post-filter strictness | **~92%** | small |
| **htsfile** | 1/1 cmd | none (intentional `-c` omission) | **~98%** | done |
| **CRAM / htsgo formats** | CRAM read+write (rANS 4x8/4x16 in-tree, LZMA via ulikunitz/xz); BGZF, BAI/CSI, SAM/BAM, BCF v2.2, tabix index | partial-decompress seek via `.gzi` (perf); BCF FORMAT-key reconstruction edge cases | **~90%** | medium |

Percentages weight *remaining feature surface and effort*, not flag counts
alone — e.g. vcftools is 100% on flags but ~97% because a few outputs still
need column polish, while bcftools is at ~70% because the remaining items
(indel calling, convert modes, csq slice 4) are individually large.

## Biggest remaining boulders

1. **bcftools `convert`'s ~18 modes** — only the pass-through round-trip
   (VCF↔BCF↔VCF.gz) is implemented. `--gvcf2vcf`, `--hapsample2vcf`,
   `--tsv2vcf`, `--gensample2vcf`, PLINK/GEN/HAP, etc. are all deferred.
   *Large.*
2. **Classic mpileup indel caller (`bam2bcf_indel.c`)** — the SNP MAQ model,
   BAQ, bias tags, the samtools-side genotype-likelihood/BCF emit path
   (slices 1–4), and bcftools `mpileup --indels-cns` (the CNS realigner,
   pure-Go `pkg/htsgo/edlib`) have all landed. What remains is the *legacy*
   (non-CNS) indel-calling path; the classic indel-model knobs that don't
   route through `--indels-cns` are accepted-but-inert today.
   *Large; shared between samtools and bcftools.*
3. **bcftools `csq` slice 4** — the haplotype engine (slices 1–3) is done
   and INFO/BCSQ is byte-parity, but `FORMAT/TBCSQ` query expansion,
   `--unify-chr-names`, `--dump-gff`, `-l/--local-csq`, and `-O b|u|z`
   non-text output remain. *Medium–large.* (csq `-b`/`-C` standard-table
   landed in #225's wave.)
4. **CRAM long-tail correctness + perf** — CRAM r/w landed, but
   partial-decompress seek via `.gzi` and some BCF FORMAT-key
   reconstruction edge cases remain. *Medium.*
5. **bedtools & vcftools option long tails** — both have no missing
   subcommands/flags but a scattered tail of option polish (bedtools:
   `bedmultiinter` VCF/GFF input, RNG byte-parity; vcftools: per-output
   column sets, `--max-indv` shuffle). Individually *small*, collectively
   *medium*.

Cross-cutting: **multi-threading (`-@`/`-t`)** has landed for the dominant
BAM-output path — `samtools view`/`sort`/`markdup` drive a parallel BGZF
compressor (decode-equal across thread counts) — and remains a deferred
no-op for parallel *input* BGZF/CRAM decode, mosdepth `-t`, bgzip's
remaining paths, and several non-BAM-writing subcommands. It's a single
deferred parallel-read pass, not per-tool feature gaps.

## Prioritized path to completion

1. **Close the small tails first** (fast wins): vcftools column polish,
   bgzip/tabix flag tails, mosdepth `.csi` + `--mapq` fast-path, prinseq
   niche knobs. Brings 5 tools to ~100%.
2. **fastp medium wave**: overrepresentation analysis, `--correction`,
   overlap-trim knobs, `--split*`, JSON schema completeness.
3. **bedtools option-tail wave**: `bedmultiinter` VCF/GFF input and the
   scattered per-subcommand flag tails.
4. **The variant-calling boulder** (highest effort, highest value):
   mpileup indel model + samtools `-g/-u` BCF output, then bcftools
   `call` full multi-allelic `-m` + BCF input.
5. **bcftools `convert` modes** and **csq slice 4** as parallel large items.
6. **Cross-cutting multi-threading pass** once feature parity is in place.

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
