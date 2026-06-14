# Project status

**Last refreshed:** 2026-06-14 (post the parity-completion wave: a shared
`pkg/htsgo/alnbed` BAM/SAM→BED12 layer brought **BAM/SAM input** to
`bedgenomecov` (`-ibam`, plus `-pc`/`-fs`), `bedjaccard`, `bedcoverage`,
`bedspacing`, and `bedgroupby`; **`--split`** (BED12 block-aware) to
`bedjaccard`/`bedgenomecov`/`bedcoverage`; **GFF/VCF input** to `bedmerge` and
**GFF** to `bedmap`; `bedclosest` multi-database (`-b … -names/-filenames/-mdb`)
and `-s`/`-S`/`-N`; corrected flag bugs in `bedsubtract -N`, `bedmerge -S`,
`bedjaccard -S`, `bedcomplement -L`, `bedmerge --delim`; `bcftools view -s`
AC/AN recompute (`-I`) and `-v/-V` type selectors, `query` position tokens
(`%POS0/%END/%END0/%FIRST_ALT/%IS_TS`), and three **BCF writer** encoding
fixes (missing-value sentinels, GT-missing, Flag) giving full upstream BCF
interop; `samtools` mpileup BAQ, `fastq -T '*'` + QNAME-based pairing, and
depth `-a`/`-b` parity. Every change is parity-validated against the vendored
upstream binaries.)

Parity-skip elimination wave (2026-06-14, PRs #286–#296): every remaining
`t.Skip` that masked a **feature-parity** gap was driven to a real
byte-for-byte (or value-exact) assertion against the vendored upstream
binaries, surfacing and fixing **7 genuine port bugs** along the way —
`bedgenomecov` per-base `-dz` 0-based offset and CIGAR-`D` block splitting,
`samtools import` `FMUNMAP` on `/1`,`/2` FASTQ suffixes, `samtools calmd`
MD/NM aux re-append ordering, `bcftools norm -c` letter semantics
(`s`=set/fix, `x`=exclude), and `bcftools concat -a` contig ordering. Closed
in this wave: `bcftools stats` (all 7 sections), `concat -a/-D`, `norm -f/-c`,
`query %N_ALT`, CSI; `mosdepth` overlap-pair detection + frag-len/`--chrom`
validation; `seqtk` `sample`/`randbase` (glibc drand48) + Mott `trimfq`;
`samtools` `import`/`calmd`/`mpileup -aa`/consensus/CRAM-region; `bedfisher`
GFF input; `bedgenomecov` BAM/SAM/CRAM input; `bedreldist`/`bedsplit`/`bednuc
-fullHeader`/`bedslop` float-precision/`bedcluster`/`bedjaccard` fixtures;
`bedpairtopair -ss`/`bedpairtobed -slop`/`bedsample` CLI rejection; and the
whole CRAM/htscodecs compliance + `prinseq` live-parity corpus (env guards
converted to hard `t.Fatalf`-with-init-hint per the parity-rig policy). The
**only** runtime skip left in the repo is the `hfile` live third-party-network
backend test (opt-in via `HFILE_NETTEST=1`; verified passing against live
https/s3/gcs public objects) — test hygiene, not a parity gap.

Previous wave (2026-06-11, ~70 PRs): `call` modes, `convert` GEN/HAP/TSV,
`annotate`, `consensus` chain, `mendelian2`, `csq` slices 1–4, mpileup
`bam2bcf_indel` + `--indels-cns`, `phase`, `gtcheck`, threading, samtools
`coverage -A`/`markdup`/`calmd`/consensus, bcftools `query %INFO/%SAMPLE`/
`roh -Oz`/`filter -M`/`cnv --AF-file`, mosdepth flags, bedtools input formats,
CRAM bzip2 encode.

This is the honest "distance to done" snapshot for the started tools. It is
**the summary source of truth**; the authoritative gap list per tool is
[`docs/PARITY_ROADMAP.md`](docs/PARITY_ROADMAP.md), and the docs map in
[`docs/README.md`](docs/README.md) says which document owns what. This file
is the skimmable summary on top of the roadmap — keep status numbers here and
in the roadmap, and link to them from everywhere else rather than re-stating.

## Project focus (current)

**No new tools; 100% parity on the ones we have.** Nothing is out of scope —
there are no "non-goals." Every upstream feature, flag, input format, and edge
case is in scope. The only scope is 100% feature parity + bug fixes + better
docs + unit/parity tests + drop-in POSIX CLIs. Remaining work is the parity
gaps in the existing tools (`docs/PARITY_ROADMAP.md`), not new ports.

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
| **prinseq** | `stats` / `filter` / `graph_data`; all in-scope flags incl. `--graph_data`, `--range_len`, `--trim_to_len`, `--seq_case`, `--line_width`, `--range_gc` (all implemented) | PNG report flow only (out of scope — needs an image renderer) | **~98%** | done (bar PNG) |
| **sickle** | `se` / `pe`; 15/15 parity cases | none (gzip-output level untested) | **100%** | done |
| **skewer** | `se` / `pe`; 14/14 parity cases | none | **100%** | done |
| **fastp** | single cmd; sliding-window, auto-adapter, HTML+JSON, dup-eval, UMI, PE base `--correction`, overrepresentation (`-p/-P`), `--split*`, merge writer (`-m`), `--adapter_fasta`, `--poly_x_min_len`, `--disable_adapter_trimming`, `--merge` (`merged_and_filtered` JSON block); 16/16 + 9 tail parity | multi-thread `--split` file-boundary distribution (perf) | **~98%** | small |
| **bedtools** | 37 bed* tools; no missing subcommands; 141+ parity tests; BAM **and CRAM** input (`bedintersect`/`bedmulticov`), VCF/GFF input (`bedintersect`/`bedmultiinter`), `bedclosest` direction flags, `intersect -c` | scattered option-tail polish | **~96%** | small (long tail) |
| **vcftools** | single cmd; **146/146 upstream long flags (100%)**, incl. BCF I/O, PCA, LD, RoH, relatedness, `--freq2`/`--counts2` (schema complete) | per-output column-set polish only; `--max-indv` uses deterministic truncation not upstream RNG shuffle (RNG-policy non-goal) | **~98%** | small |
| **bcftools** | 24 subcommands (all present); mpileup MAQ SNP model (slices 1–4) + legacy `bam2bcf_indel` + `--indels-cns` (edlib realigner) + BAQ + bias tags; full multi-allelic `call` (`-m`/`-c`/`--gvcf`/`-C alleles`/`-G`/`--ploidy GRCh37/38`/`--ploidy-file`); `convert` GEN/HAP/TSV/gVCF modes + PLINK exporters (`--plink`/`--tped`/`--plink-bed`, PLINK1 spec); `gtcheck`/`mendelian2`/`consensus` (chain+iupac)/`annotate`; `filter -M`/`cnv --AF-file`/`roh -Oz`/`query %INFO/%SAMPLE`; csq slices 1–4 (FORMAT/TBCSQ, --unify-chr-names, -O b\|u\|z, --dump-gff); full HMM `roh`/`cnv`/`polysomy`; subprocess plugin system; `csq -l/--local-csq` (test_cds_local); `gtcheck` `-i`/`-e` filter expressions (qry:/gt: scope); `concat --ligate`; remote URL inputs (`view -r`/`query -r` via hfile); `som` train/classify (upstream write bug fixed); `gtcheck -c/--cluster` clustering (own design — upstream is an error stub) | `query %N_ALT` (non-goal) | **~98%** | small |
| **samtools** | 25 functional subcommands; CRAM r/w (v2/v3/**v4.0**) + bzip2 encode; `.csi`; consensus `--het-only` + indel calling + **`-a` placeholder rows** (ref-skip/del/zero-cov, live parity); `coverage -A`; `markdup -d/-s/-S`; `calmd -C/-e/-u`; `phase`; mpileup MAQ **BCF/VCF emit (slices 1–4) + legacy indel + --indels-cns**; **`tview` text/HTML/interactive `-d C`** (byte-for-byte text/HTML, pure-Go termios for `-d C`); remote URL inputs; `-@` threading for view/sort/markdup | `-@` parallel *input* BGZF/CRAM decode for non-BAM-writing subcommands (perf only) | **~99%** | small (perf) |
| **mosdepth** | single cmd; ALL flags wired incl. `-d/--d4`, `--fragment-mode`, `--quantize`, `-t/--threads`, `--use-median`, `--mapq` fast-path, **CRAM input** (`-f/--fasta`, via alnio); emits `.csi` | none | **~99%** | done |
| **bgzip** | 1/1 cmd, all flags incl. parallel block compression (`-@`/`-t`) + `--test` integrity check | none | **~99%** | done |
| **tabix** | 1/1 cmd, all flags incl. `--reheader`, strict `--targets` post-filter, remote URL region queries | none | **~99%** | done |
| **htsfile** | 1/1 cmd | none (intentional `-c` omission) | **~98%** | done |
| **CRAM / htsgo formats** | CRAM read+write **v2/v3 + v4.0** (uint7-varint + 64-bit pos; v4 transform codecs XPACK/XRLE/XDELTA decode; rANS 4x8/4x16 in-tree, LZMA via ulikunitz/xz); writer output (v3 **and** v4) decodes through upstream samtools (live cross-checks); BGZF, BAI/CSI, SAM/BAM, BCF v2.2, tabix index | partial-decompress seek via `.gzi` (perf only) | **~95%** | small |

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

1. **htsgo `hfile` cloud I/O — DONE.** The `pkg/htsgo/hfile` backend
   (HTTP(S)/S3/GCS, stdlib-only, hand-rolled AWS SigV4 + GCS bearer token,
   read-ahead-buffered `OpenSeekable`) is complete and **validated against
   live infrastructure**: HTTPS, anonymous S3 (`s3://1000genomes`), and public
   GCS (`gs://gcp-public-data--broad-references`) backends; a `tabix` region
   query over the 214 MB 1000 Genomes chr22 VCF on S3; and a `samtools view`
   region count over a 188 MB exome BAM on S3 (both fetch only the index + the
   chunk). Remote URLs flow through the streaming opens (`iohelper`/`alnio` —
   whole-file ops, incl. live CRAM decode) and **every** indexed region-query
   path: `samtools view`/`idxstats`/`mpileup`, `tabix`, `bcftools view -r` /
   `query -r`, and `.crai`-indexed CRAM `samtools view` (seek-based container
   reader, live-upstream-parity validated).
2. **CRAM long-tail correctness + perf.** **CRAM v4.0 is implemented for read AND write**
   (uint7-varint + 64-bit positions, version-threaded; the v4 transform
   codecs XPACK/XRLE/XDELTA decode; v4.0 decode is 15/15-record byte-for-byte
   vs an upstream-written v4.0 CRAM, and our v4.0 **and** v3.0 writer output
   decodes field-for-field through upstream `samtools view`). (`samtools
   view -T` over CRAM regenerates `MD:Z`/`NM:i` via `pkg/htsgo/mdnm`; the
   network REF_PATH/EBI reference fetch is implemented — both live-parity /
   httptest validated.) *Small.*
3. **Scattered option-tail polish — closed.** bcftools `concat --ligate`
   (live-parity), mendelian2 `sites_not_diploid`, and vcftools's BCF-binary
   I/O family (`--bcf`/`--recode-bcf`/`--diff-bcf`/`--contigs`,
   roundtrip-tested) are all implemented; vcftools has no long-flag gaps.
   `convert` PLINK exporters (PLINK1 spec), `som` train/classify (upstream
   write bug fixed), `gtcheck -c/--cluster` clustering, and `samtools tview`
   text/HTML (byte-for-byte) are now all implemented. No actionable
   option-tail gaps remain.

Cross-cutting: **multi-threading (`-@`/`-t`)** has landed for the dominant
BAM-output path — `samtools view`/`sort`/`markdup` and mosdepth `-t` drive a
parallel BGZF compressor / decompressor (decode-equal across thread counts)
— and remains a deferred no-op for parallel *input* BGZF/CRAM decode in
several non-BAM-writing subcommands and bgzip. It's a single deferred
parallel-read pass, not per-tool feature gaps.

## Non-goals (not gaps)

These are deliberately not ported and should not be counted against parity:

- **`samtools tview`** — the non-interactive text (`-d T`) and HTML (`-d H`)
  viewer modes are **implemented** (byte-for-byte parity with upstream, see
  `tools/samtools/README.md`); only the interactive ncurses **`-d C`** mode is
  a deliberate non-goal (no pipeline use, would require a TTY UI dependency).
- **`query %N_ALT`** and **`import --skipBamQ`** — **not** upstream flags.
- **`bcftools gtcheck -c/--cluster`** — the cluster/dendrogram option is
  **commented out** in upstream's usage (`vcfgtcheck.c`); it is unadvertised
  dead surface, so there is nothing to port. (`gtcheck`'s `-i`/`-e` filter
  expressions, the real feature, are implemented.)
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
