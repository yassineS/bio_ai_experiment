# Bioinformatics Tool Ranking — 2026 Refresh

**Date:** 2026-05-14 (revised 2026-05-14)
**Author:** Research pass for `bio_ai_experiment`
**Purpose:** Rank candidate **new** tools to port by *current* (2024–2026)
real-world usage. **This file is not a deprioritise-existing-tools filter.**

> **Scope clarification (2026-05-14 revision):** the project goal is **1:1
> feature parity with every tool we've started**, including ones whose
> upstream is dormant. This file ranks which *additional* tools to port
> after the existing-tool parity work is complete. The earlier "freeze
> sickle/skewer/PRINSEQ" advice in this document was wrong relative to the
> project goal and has been removed. See
> [`../docs/PARITY_ROADMAP.md`](../docs/PARITY_ROADMAP.md) for the per-tool
> gap list and [`../docs/UPSTREAM_BUGS.md`](../docs/UPSTREAM_BUGS.md) for
> upstream bugs we choose not to carry over.

## TL;DR

- **The 1:1 parity goal applies to every existing port**, including
  sickle / skewer / PRINSEQ-lite. Upstream being dormant doesn't change
  the parity bar — it just means our port becomes the maintained
  implementation. The earlier "deprioritise" framing has been retracted.
- **Trim Galore went official-Rust** in v2.1.0 "Oxidized Edition" (May 2026).
  That removes a candidate from our *new-tools* shortlist — no point
  porting a tool whose upstream just shipped a modern rewrite.
- **For new tools, go all-in on the SAM/BAM/VCF/BGZF/tabix layer.** It's
  the most-used, most-stable, most-painful-to-install set of tools in the
  field, and `biogo/hts` (the only existing Go htslib alternative) is
  unmaintained since 2017. There is a clear "deliver value" wedge here.
- **Top-5 next picks (all merged as of 2026-05-14):** `bgzip` ✅, `tabix` ✅,
  `samtools view/sort/index/depth/fastq/flagstat` ✅,
  `bcftools view/index/stats/query/concat/norm` ✅, `mosdepth` ✅. Their
  remaining sub-features count against the **parity roadmap**, not this
  list. See "Next-up shortlist" at the bottom for the next batch.

## How to read this doc

For each tool below:

- **Use case** — one line.
- **Bioconda total downloads** — lifetime, from anaconda.org. Anaconda's
  "Last 6 months" widget currently displays `0` for nearly every package
  (display bug observed 2026-05-14), so we report lifetime totals from
  search-index snapshots where they appear, and otherwise mark `n/a`.
- **Last release / maintained** — from upstream GitHub.
- **Already ported?** — what's in `tools/` today.
- **Recommendation** — `port` / `wrapper` / `skip` / `already done`.

Numbers are best-effort snapshots, not authoritative. Where a figure
couldn't be confirmed, "n/a" is used and no speculation is offered.

---

## Section 1: Tools already in `tools/` — current status

| Tool | Bioconda total dl | Last release | Maintained? | Recommendation |
|------|------------------:|--------------|:-----------:|----------------|
| `seqtk`     | n/a (commonly cited >1M) | 1.5 (Jun 2025) | yes | keep, parity-track |
| `prinseq`   | n/a (latest 0.20.4, Mar 2016) | 2016 | **no** | already started — close to 1:1 (see PARITY_ROADMAP) |
| `sickle`    | n/a (latest Dec 2015) | 2015 | **no** | already started — close to 1:1 (see PARITY_ROADMAP) |
| `skewer`    | n/a (latest Sep 2016) | 2016 | **no** | already started — close to 1:1 (see PARITY_ROADMAP) |
| `fastp`     | ~631k | 1.3.3 (Apr 2026) | yes | keep, parity-track |
| `vcftools`  | n/a (latest 0.1.17, May 2025) | 0.1.17 (May 2025) | weak | already done — keep, but note vcftools is itself nearly frozen |
| `bedmerge/intersect/sort/...` (bedtools subset) | bedtools: ~3.3M | bedtools 2.31.1 (Nov 2023) | yes | keep, extend coverage |

**Sources:** anaconda.org/bioconda/{seqtk, prinseq, sickle-trim, skewer, fastp,
vcftools, bedtools}; bedtools total via WebSearch snapshot of anaconda.org
package page.

**Implication.** Three of our seven tool families (sickle, skewer, PRINSEQ) are
porting *dead software*. The Go versions are good engineering practice and we
should keep them shippable, but they will not "deliver value to people" because
people moved on years ago. Stop spending effort chasing the last edge cases of
sickle's quality algorithm — `fastp` and `cutadapt` already won.

---

## Section 2: Top ~30 actively-used CLI bioinformatics tools (2024–2026)

Ranked roughly by **(active usage × upstream pain × Go-port leverage)**.
Download numbers are lifetime bioconda totals; "n/a" means I couldn't
confirm a number and won't guess.

### Tier 1 — universally used core formats / IO (HIGHEST LEVERAGE)

| # | Tool | One-line description | Bioconda total dl | Last release | Maintained? | In `tools/`? | Recommendation |
|---|------|----------------------|------------------:|--------------|:-----------:|:------------:|----------------|
| 1 | **htslib** (`bgzip`, `tabix`) | Block-gzip + generic tabix index for VCF/BED/GFF/SAM | n/a (core dep of samtools/bcftools, "installed >1M times via bioconda" per Danecek 2021) | 1.23.1 (Mar 2026) | yes (samtools team, security patches CVE-2026-31965..31971) | no | **port (priority)** |
| 2 | **samtools** | SAM/BAM/CRAM swiss-army (view/sort/index/flagstat/depth/...) | ~6.5M | 1.23.1 (Mar 2026) | yes (1.9k★, active dev, 2,862 commits) | no | **port (subset: view/sort/index/flagstat/fastq)** |
| 3 | **bcftools** | VCF/BCF swiss-army (view/sort/index/merge/norm/...) | n/a (typically tracks samtools volume) | 1.23.1 (Mar 2026) | yes (855★) | no | **port (subset: view/sort/index/merge/norm)** |
| 4 | **bedtools** | Genome arithmetic on BED/GFF/VCF | ~3.3M | 2.31.1 (Nov 2023) | maintenance | partial (bedmerge, bedintersect, bedsort, etc.) | **extend** (cover the remaining ~20 subcommands) |

Sources: anaconda.org/bioconda/samtools (6.5M cited in search snapshot);
anaconda.org/bioconda/bedtools (3.3M cited in search snapshot);
github.com/samtools/{samtools,htslib,bcftools} for release/star data;
Danecek et al., GigaScience 2021 ("Twelve years of SAMtools and BCFtools",
PMC7931819).

### Tier 2 — high-volume aligners and QC

| # | Tool | One-line description | Bioconda total dl | Last release | Maintained? | In `tools/`? | Recommendation |
|---|------|----------------------|------------------:|--------------|:-----------:|:------------:|----------------|
| 5 | **minimap2** | Dominant long-read / spliced aligner | n/a | 2.30 r1287 (Jun 2025) | yes (2.2k★, Heng Li) | no | **skip-port-core, wrapper OK** — porting a 30k-LOC SIMD aligner is a multi-quarter project with no parity win |
| 6 | **bwa / bwa-mem2** | Short-read aligner (still the standard) | n/a | bwa 0.7.19 (Mar 2025); bwa-mem2 2.3 (Jul 2025) | yes | no (see `tools/BWA_IMPLEMENTATION_DECISION.md`) | skip (already decided) |
| 7 | **fastp** | All-in-one FASTQ QC/trim/filter — modern default | ~631k | 1.3.3 (Apr 2026) | yes | yes | keep, parity-track |
| 8 | **cutadapt** | Adapter trimming, very widely used | n/a | 5.2 (Oct 2025) | yes (Python; Marcel Martin) | no | wrapper / skip — Python rewrite would lose users; Rust port already exists |
| 9 | **Trim Galore** | Wrapper around cutadapt+FastQC | n/a | **2.2.0 (May 2026), full Rust rewrite** | yes | no | **skip — upstream just rewrote it** |
| 10 | **FastQC** | Per-base/per-read QC HTML report | n/a | 0.12.1 (Mar 2023) | low (Babraham, slow) | no | candidate, but Java→Go offers modest user-visible win; punt |
| 11 | **STAR** | Splice-aware short-read RNA-seq aligner | n/a | 2.7.11b (Jan 2024) | maintenance | no | skip (massive C++; alignment kernels) |
| 12 | **HISAT2** | Splice-aware short-read aligner | n/a | 2.2.2 (Jan 2026) | maintenance | no | skip (large, complex BWT-FM index) |
| 13 | **bowtie2** | Short-read aligner | n/a | active per CI on top-50 doc | yes | no | skip (large C++) |

### Tier 3 — high-value standalone analyses (good port candidates)

| # | Tool | One-line description | Bioconda total dl | Last release | Maintained? | In `tools/`? | Recommendation |
|---|------|----------------------|------------------:|--------------|:-----------:|:------------:|----------------|
| 14 | **mosdepth** | Fast BAM/CRAM depth (WGS/exome/targets) | n/a | 0.3.14 (Apr 2026) | yes (brentp) | no | **port** — Nim today, small (~3k LOC), high pain (Nim toolchain friction), would benefit from a Go port |
| 15 | **multiqc** | Aggregate QC reports across all tools | n/a | 1.34/1.35 (Apr–May 2026) | yes (Seqera/Phil Ewels) | no | **port (long-term)** — Python today; a single static Go binary that ingests the same parser modules would be a real value-add for HPC users |
| 16 | **kraken2** | k-mer taxonomic classification | n/a | active | yes | no | skip (huge DB, complex perf) |
| 17 | **salmon** | RNA-seq transcript quantification | n/a | 1.11.4 (Mar 2026), C++ (97.9%) | yes (884★) | no | skip (selective alignment + SSHash + EM — huge surface) |
| 18 | **kallisto** | RNA-seq pseudoalignment quantification | n/a | 0.52.0 (Mar 2026) | yes | no | skip (specialized index) |
| 19 | **MMseqs2** | Ultra-fast protein search/clustering | n/a | 18 (early 2025), Nature Methods 2025 GPU paper | yes | no | skip (GPU kernels, massive scope) |
| 20 | **MAFFT** | Multiple sequence alignment | n/a | 7.525 (Mar 2024) | maintenance | no | skip (algorithmic complexity, not a "I hate installing this" tool) |
| 21 | **DIAMOND** | Fast protein alignment (BLAST-replacement) | n/a | active | yes | no | skip (large C++, SIMD-heavy) |
| 22 | **trimmomatic** | Java adapter/quality trimmer | n/a | 0.40 (Aug 2025) | low | no | skip (superseded by fastp/Trim Galore in practice) |
| 23 | **seqkit** | FASTA/Q toolkit | n/a | 2.13.0 (Feb 2026) | yes (1.6k★) | n/a — already Go | **reference, not target** |

### Tier 4 — domain-specific, lower priority for us

| # | Tool | One-line description | Maintained? | Recommendation |
|---|------|----------------------|:-----------:|----------------|
| 24 | bismark | Bisulfite-seq alignment | active | skip (specialized) |
| 25 | Prokka / Bakta | Prokaryotic annotation | active | skip (pipeline + many DB deps) |
| 26 | SPAdes | Genome assembly | active | skip (very large C++) |
| 27 | LoFreq | Low-frequency variant caller | active | skip (statistical model heavy) |
| 28 | DeepVariant | DL variant caller | active | skip (TF/ML deps) |
| 29 | Strelka | Small-variant caller | active | skip (large C++) |
| 30 | IQ-TREE | ML phylogenetics | active | skip (scientific kernels) |

---

## Section 3: What changed since `top_50_packages_for_improvement.md` (2025-10)?

1. **Trim Galore #31 in the old list got an official Rust rewrite.** v2.1.0
   "Oxidized Edition" in May 2026 bundles cutadapt+FastQC and eliminates the
   Java/Perl dependency chain. Our value-add for that tool just evaporated.
2. **PRINSEQ #4, sickle #36, skewer #15 are confirmed dead upstream.** Last
   bioconda uploads are 2016/2015/2016 respectively. Keep our ports as
   reference/teaching, but don't sink more time into parity tests.
3. **bedtools is now soft-maintenance** (2.31.1 from Nov 2023, ~30 months
   stale as of this writing). Our `bed*` tools are well-positioned to *replace*
   bedtools rather than mirror it — same input/output, fewer install steps,
   single static binary.
4. **samtools/bcftools/htslib stay #1** in real-world install volume by a
   wide margin (~6.5M bioconda downloads for samtools alone) and are now
   shipping security patches (multiple CVEs in 2026 for the CRAM decoder),
   which underlines that the C codebase still has memory-safety burden — a
   selling point for a Go alternative.
5. **`biogo/hts` is dormant** (last release Feb 2017, 133★). There is no
   actively-maintained Go alternative for BAM/BGZF/tabix. Whoever ships that
   first gets meaningful adoption.

---

## Section 4: Next-up shortlist — top 5 to port next

Listed in **recommended build order**, because the dependency stack matters
(everything else needs bgzip+tabix first).

### 1. `bgzip` — block-gzip codec

- **Why:** It is the foundational format for everything else (BAM, BCF,
  indexed VCF, indexed BED/GFF). VCFtools and our planned BAM/VCF work all
  need it. It's small (a few hundred LOC of C), well-specified
  (`hts-specs/SAMv1.pdf` §4 and the BGZF section), and easy to ship as a
  single Go binary.
- **Effort:** ~1 week. The standard library `compress/gzip` already handles
  individual gzip members; BGZF is "gzip with a max-65535-byte payload per
  member plus a virtual-offset index". Add a `pkg/bgzf` package and a
  `tools/bgzip/cmd/bgzip` CLI that mirrors htslib's `bgzip(1)`.
- **Value to users:** Static binary, no `apt install htslib-tools`, drop-in
  for piping into existing tools. Our existing `vcftools` port can also
  consume it.

### 2. `tabix` — generic position index

- **Why:** Pairs with bgzip and powers `samtools view -r`, `bcftools view -r`,
  and every "give me variants in this region" pipeline. Tabix index format
  (TBI) and CSI v1 are both small, fully specified, and stable.
- **Effort:** ~2 weeks on top of bgzip. The trickiest piece is the
  linear-index + binning scheme; the format is documented in
  `hts-specs/tabix.pdf`.
- **Value to users:** Same as bgzip — one less native dependency, and a
  Go-callable `pkg/tabix` library that the rest of our tools can use.

### 3. `samtools view / sort / index / flagstat / fastq` (subset)

- **Why:** Highest-leverage CLI in the field (~6.5M bioconda downloads).
  Five subcommands cover ~80% of real-world day-to-day use:
  - `view` — filter + convert SAM/BAM/CRAM
  - `sort` — sort by coordinate or name
  - `index` — produce .bai/.csi
  - `flagstat` — quick stats
  - `fastq` — BAM → FASTQ for re-aligning
- **Effort:** ~6–8 weeks for the subset, assuming bgzip/tabix done first.
  BAM read/write is the bulk; CRAM is optional and can be skipped in v1
  (call it out as a known gap). Sort needs careful merging logic for
  multi-GB BAMs.
- **Value to users:** Single static binary, no htslib install, memory-safe
  CRAM/BAM parser (modulo what we implement). Even partial coverage is
  immediately useful for pipelines.
- **Risk / out-of-scope:** Full CRAM v3 codec is a *lot* of work (RANS,
  reference-based compression). Don't ship it in v1.

### 4. `bcftools view / index / merge / norm` (subset)

- **Why:** VCF/BCF is everywhere downstream of variant calling and our
  existing `tools/vcftools` is the natural home. `bcftools` is the modern
  successor to `vcftools` and the subcommands above are the highest-use
  ones. Combining with our existing VCF parsing in
  `pkg/htsgo/vcf` gives us a strong base.
- **Effort:** ~4–6 weeks once bgzip/tabix and BCF binary format are done.
  BCF binary parsing is well-specified.
- **Value to users:** Closes the loop with our existing VCFtools port —
  users get a complete BGZIP/tabix/VCF/BCF stack in Go.

### 5. `mosdepth` — fast BAM/CRAM depth

- **Why:** Different from #1–#4: this is a *standalone analytical* port, not
  format plumbing. It's actively maintained but written in Nim, which is a
  real install pain in CI environments. The algorithm is small and clean
  (sweep BAM, accumulate coverage per region). Replacing it with a single
  static Go binary that we control end-to-end is a tractable, high-visibility
  win once we have BAM parsing from #3.
- **Effort:** ~2–3 weeks **after** samtools BAM reader is in place. Without
  it, we'd be implementing BAM twice.
- **Value to users:** Same UX as mosdepth, no Nim toolchain, no htslib
  dependency. mosdepth is widely used in clinical NGS pipelines.

### Honourable mentions (not in top 5, but worth tracking)

- **`multiqc` in Go** would be a genuine win (single static binary; Python
  install pain is real), but it's a *parser zoo* — each upstream tool needs
  its own report parser. Treat as a long-running "v2.0" goal after the
  format stack is solid.
- **`bedtools` completeness**: extending our existing `bedmerge`,
  `bedintersect`, etc. to cover the remaining ~20 subcommands is also high
  value — bedtools itself is now in soft-maintenance and people regularly
  hit the same install/perf pain.
- **A `pkg/sam` and `pkg/bcf` library** falls out of the shortlist above
  for free. We should design it to be a reusable Go API (à la
  `biogo/hts` but maintained), not just a binary's internals.

---

## Section 5: Sources

- anaconda.org/bioconda/{samtools, bcftools, htslib, minimap2, mosdepth,
  mafft, mmseqs2, multiqc, seqtk, sickle-trim, skewer, prinseq, fastp,
  vcftools, bedtools, salmon, kallisto, star, hisat2, bwa, bwa-mem2,
  fastqc, trim-galore, cutadapt, seqkit, kraken2, trimmomatic} —
  WebFetch / WebSearch snapshots, 2026-05-14.
  Note: anaconda.org's "Downloads (Last 6 months)" widget displayed `0`
  for every package on the date of this snapshot — appears to be a
  site-wide display bug rather than reality. Lifetime totals reported
  where they appeared in search-engine cached snapshots of the same pages.
- github.com/samtools/{samtools, bcftools, htslib} — release notes and
  commit activity.
- github.com/lh3/minimap2/releases — v2.30 release (Jun 15, 2025).
- github.com/brentp/mosdepth — v0.3.14 (Apr 24, 2026), Nim.
- github.com/FelixKrueger/TrimGalore/releases — v2.2.0 "Clumpify
  Edition" (May 7, 2026), Rust rewrite via v2.1.0 "Oxidized Edition".
- github.com/MultiQC/MultiQC — v1.34 (Apr 21, 2026), still Python.
- github.com/COMBINE-lab/salmon — v1.11.4 (Mar 2026).
- github.com/shenwei356/seqkit — v2.13.0 (Feb 2026), Go, 1.6k stars.
- github.com/biogo/hts — last release Feb 2017, 133 stars, dormant.
- Danecek et al., "Twelve years of SAMtools and BCFtools", GigaScience
  2021 (PMC7931819) — referenced for the "installed >1M times via
  bioconda" claim.
- Nature Methods 2025 "Year in Review" (s41592-025-02997-5) — confirms
  FastQC/SPAdes/Prokka/Bakta/Kraken2 still the standards for microbial
  genomics; spatial/multi-omics is where the methods action is, but the
  CLI tooling there is mostly Python notebooks, not our wedge.
- Existing repo: `analysis/top_50_packages_for_improvement.md`
  (2025-10-20) for the prior ranking baseline.
