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
documented intentional deviations and fixed-on-port upstream bugs), and the
validated-parity suite passes with an explicit `t.Skip()` only for the few
genuinely untestable cases. "Working subset" is a
milestone, not the finish line. The percentages below are an honest,
evidence-based estimate of remaining surface — not a rosy reading.

## Per-tool completion

| Tool | Implemented | Remaining surface | % | Effort to finish |
|------|-------------|-------------------|---:|------------------|
| **seqtk** | All 24 subcommands; byte-parity vs v1.5 | none | **100%** | done |
| **prinseq** | `stats` / `filter` / `graph_data` / `graph_png`; all in-scope flags incl. `--graph_data`, `--range_len`, `--trim_to_len`, `--seq_case`, `--line_width`, `--range_gc`; the `prinseq-graphs.pl` PNG report flow is implemented (`graph_png`) | none (PNG byte-identity N/A — pure-Go renderer, not Perl Cairo/GD) | **100%** | done |
| **sickle** | `se` / `pe`; 15/15 parity cases | none (gzip-output level untested) | **100%** | done |
| **skewer** | `se` / `pe`; 14/14 parity cases | none | **100%** | done |
| **fastp** | single cmd; sliding-window, auto-adapter, HTML+JSON, dup-eval, UMI, PE base `--correction`, overrepresentation (`-p/-P`), `--split*`, merge writer (`-m`), `--adapter_fasta`, `--poly_x_min_len`, `--disable_adapter_trimming`, `--merge` (`merged_and_filtered` JSON block); 16/16 + 9 tail parity | multi-thread `--split` file-boundary distribution (perf) | **~98%** | small |
| **bedtools** | 37 bed* tools; no missing subcommands; 141+ parity tests; BAM **and CRAM** input (`bedintersect`/`bedmulticov`), VCF/GFF input (`bedintersect`/`bedmultiinter`), `bedclosest` direction flags, `intersect -c`; **`bedcoverage --split` over a blocked query (`-a`)** — a BED12 record or a spliced (`N`-CIGAR) BAM alignment — across every mode (default/`--counts`/`--depth`/`--hist`/`--mean`, `-s`/`-S`), byte-validated vs upstream 2.31.1 (intron bases excluded from overlap but kept in the full-span length/depth vector; B straddling an intron counted per overlapped block) | scattered option-tail polish | **~96%** | small (long tail) |
| **vcftools** | single cmd; **146/146 upstream long flags (100%)**, incl. BCF I/O, PCA, LD, RoH, relatedness, `--freq2`/`--counts2` (schema complete) | per-output column-set polish only; `--max-indv` ports upstream's glibc `rand()` + `std::random_shuffle` and is byte-exact for `--max-indv-seed` (upstream itself is time-seeded/non-reproducible) | **~98%** | small |
| **bcftools** | 24 subcommands (all present); mpileup MAQ SNP model (slices 1–4) + legacy `bam2bcf_indel` + `--indels-cns` (edlib realigner) + BAQ + bias tags; full multi-allelic `call` (`-m`/`-c`/`--gvcf`/`-C alleles`/`-G`/`--ploidy GRCh37/38`/`--ploidy-file`); `convert` GEN/HAP/TSV/gVCF modes + PLINK exporters (`--plink`/`--tped`/`--plink-bed`, PLINK1 spec); `gtcheck`/`mendelian2`/`consensus` (chain+iupac)/`annotate`; `filter -M`/`cnv --AF-file`/`roh -Oz`/`query %INFO/%SAMPLE`; csq slices 1–4 (FORMAT/TBCSQ, --unify-chr-names, -O b\|u\|z\|t, --dump-gff); full HMM `roh`/`cnv`/`polysomy`; subprocess plugin system; `csq -l/--local-csq` (test_cds_local); `gtcheck` `-i`/`-e` filter expressions (qry:/gt: scope); FORMAT/GT/sample-level filter engine (`view`/`filter` `-i`/`-e`, e.g. `GT="het"`, `FMT/DP>10`, with per-sample masks); `concat --ligate`; remote URL inputs (`view -r`/`query -r` via hfile); `som` train/classify (upstream write bug fixed); `gtcheck -c/--cluster` clustering (own design — upstream is an error stub); **complete native plugin catalogue — all 41 upstream `+<name>` plugins reimplemented in pure Go (CLI-to-CLI byte-parity vs 1.23.1), exec subprocess fallback retained**, including the native-plugin `-i`/`-e` site/sample pre-filter modes (guess-ploidy, smpl-stats, indel-stats, contrast, trio-stats; split per-output; scatter accepts-and-ignores, matching upstream's no-op), plus the **curly-brace multi-threshold `-i/-e` expansion** for smpl-stats/indel-stats/trio-stats (`-i 'GQ>{10,20,30}'` runs the stats once per expanded threshold, cartesian over multiple `{...}` groups, each in its own `FLT*` report section — byte-parity vs 1.23.1); **native-plugin `-r/-R/-t/-T` region/target selection** — a shared host-side filter (`region_target.go`) applies `-r`/`-R` (span-overlap) and `-t`/`-T` (start-position, with `^` negation) before records reach each plugin, byte-validated vs upstream 1.23.1 across check-sparsity, remove-overlaps, prune, smpl-stats, indel-stats, contrast, guess-ploidy (`-r/-R` only — its `-t` is `--tag`), mendelian2, trio-stats, isecGT (both streams), split, scatter; check-sparsity keeps its own region-labelled grouping + fix-on-port accepts BED/TSV `-R` files that upstream silently drops (see `docs/UPSTREAM_BUGS.md#bcftools-check-sparsity-regions-file`); **native-plugin `-W/--write-index[=csi\|tbi]` output auto-indexing** — writes a CSI (default) or TBI index next to each indexable `.vcf.gz`/`.bcf` output (non-indexable plain-VCF/stdout errors exactly as upstream), byte-validated vs `bcftools index` across contrast, isecGT, mendelian2, split and scatter (multi-output plugins index every emitted file); the CSI writer now emits the trailing `n_no_coor` field even when zero to match htslib byte-for-byte | `query %N_ALT`; `parental-origin` (dup/del, byte-parity via in-tree kfunc `kf_betai`) and `color-chrs` (HMM Viterbi `.dat`, byte-parity) are ported; `trio-dnm3` is now **complete** — `--use-NAIVE` (GT-only Mendelian DNM/VA, byte-parity) plus the three float models `--use-DMM`/`--use-ALM`/`--use-DNG` (Dirichlet-multinomial / allele-likelihood / DeNovoGear de-novo scores over FORMAT/AD/PL/QS/QM/SP, with `--dng-priors`, chrX/PAR, `-n`, the `--pn`/`--pns`/`--phi`/`--max-QM`/`--min-vaf`/`--mrate`/`--noise-prior`/`--strand-bias` knobs, `>4`-allele trimming, and `DNM:log`/`phred`/`prob` + `VA` + `VAF` outputs). The float scores go through the bit-stable in-tree `kfBetai`/`kfLgamma` and Go's `math`, validated with a tolerance-aware ("proximity") comparison (string fields exact, scores within ~6 sig-figs of upstream — the libm last-ULP boundary); on linux/amd64 they land byte-for-byte; split-vep `-i`/`-e` (filters over split-vep's derived per-transcript CSQ columns, not plain VCF fields); **split-vep now closes its last five modes — the `primary`/`pick`/`mane`/EXPRESSION transcript selectors, the `:worst` PRN qualifier, `-g/--gene-list [+]FILE` restrict/prioritise (+`--gene-list-fields`), `-S/--severity -\|FILE` custom scale, and `--columns-types -\|FILE` type table — all byte-validated vs 1.23.1 (see `docs/PARITY_ROADMAP.md`)**; **`+setGT` now closes its last three modes — `-t b:TAG CMP VAL` binomial (two-tailed `binom.test` over FORMAT/AD via the in-tree `kfBetai`, bit-exact), `-t r:FLOAT` random with `-s/--seed` (in-tree deterministic drand48 LCG → byte-parity for any fixed seed, forced serial), and `-n X` largest-FORMAT/AD allele — all byte-validated vs 1.23.1**; **`+prune` now closes all seven previously-unsupported modes — `rand` selection + `--random-seed`/`--randomize-missing` (reusing the in-tree drand48 → byte-parity), `-m R2=/LD=/RD=` LD thresholding and `-a r2/LD/RD/count` annotation (the LD math a byte-exact `+-*/`/`sqrt` port of `calc_ld`), `-f LABEL` soft-filter, `--keep-sites`, and default `maxAF` without `--AF-tag` — all byte-validated vs 1.23.1**; **`+fill-tags` now closes its last four modes — `-l/--list-tags` (exact stderr table), `-S/--samples-file` population grouping (per-pop `_GROUP`-suffixed AN/AC/AF/MAF/NS/HWE/ExcHet + the summary `ALL` pop), the full built-in tag set (incl. HWE/ExcHet via the in-tree Wigginton exact test and sites-only AF-from-AN/AC), and the custom expression `TAG[:Number]=[int\|float](EXPR)` via an in-tree evaluator (sum/avg/max/min/median/stdev + smpl_* variants, arithmetic, ABS/PHRED, F_MISSING/N_MISSING/F_PASS/N_PASS) — all byte-validated vs 1.23.1**; **the stats/contrast plugins' remaining flags are closed — `-o/--output FILE` (smpl-stats/indel-stats/trio-stats; report bytes identical to stdout), indel-stats `-p/--ped` de-novo indel mode (+`--alt2ref-DNM`; fix-on-port: reports AD-less PED indel VCFs that upstream aborts on, see `docs/UPSTREAM_BUGS.md`), trio-stats `-a/--alt-trios` deferred singleton/doubleton accounting, and contrast `-f/--max-allele-freq` rare-allele enrichment (the extra `max_AC/PASSOC/FASSOC/NASSOC:` stderr summary) — all byte-validated vs 1.23.1 with binary-free `TestUnit*` coverage of the pure helpers**; **the format/output-mode tails are closed — `+ad-bias --clean-vcf` (allele-subset VCF via an in-tree `bcf_remove_allele_set` port) + `-f` query-format report column, `+remove-overlaps --missing 0\|DP` (the min(QUAL) DP coverage heuristic) + `-Ot/-Otz` text-list output, `+tag2tag --LXX-to-XX`/`--LPL-to-PL`/`--LAD-to-AD` localized-allele expansion (LAA-mapped Number=G PL / Number=R AD, with `-d`/`-s`), `+guess-ploidy -g/--genome` b37/b38/hg19/hg38 non-PAR-chrX presets (→`-r REGION`), and `+af-dist -p/-d` bins-from-file — all byte-validated vs 1.23.1 (one remove-overlaps ring-wrap stale-mark corner fixed-on-port, see `docs/UPSTREAM_BUGS.md`) with binary-free `TestUnit*` coverage**; **`+gvcfz` and `+frameshifts` are now genuinely ported (previously registration stubs that errored from Init): `+gvcfz` resizes gVCF blocks by the `-g FILTER:EXPR[;…]` groups (each a full FORMAT/GT filter expression via the in-tree filter engine — single `&`/`|` now accepted alongside `&&`/`||`), merging INFO/END + FORMAT DP/GQ\|RGQ/PL across the collapsed block, with `-i/-e` pre-filter and `-a`; `+frameshifts` annotates INFO/OOF over an exon BED/region-list through an in-tree port of htslib's `bcf_sr_regions_overlap` cursor — both byte-validated vs 1.23.1 (plain-BED and tabixed-BED cursor paths), with binary-free `TestUnit*` coverage. The frameshifts default reproduces upstream's dead-code `OOF=-1` exactly (see `docs/UPSTREAM_BUGS.md#bcftools-frameshifts-oof-dead-code`); the corrected in-frame/out-of-frame computation is available via the opt-in `--fix-oof` flag**; **`+fixref` now closes its last mode — `id`/`--use-id FILE` (MODE_USE_ID): the REF allele is determined from a separate dbSNP VCF keyed by the ID column (per-input-chromosome rsID→{pos,ref} map; REF match → `none`, ALT match → `swap`+GT flip, otherwise → unresolved/`skip`; position-correction + the one-shot unsorted-VCF warning reproduced), byte-validated vs 1.23.1 on BOTH the corrected VCF and the stderr stats summary, with binary-free `TestUnitFixref*` coverage. Fix-on-port: our reader also accepts a plain (un-bgzipped/un-indexed) dbSNP VCF that upstream's synced reader refuses (see `docs/UPSTREAM_BUGS.md#bcftools-fixref-id-plain-vcf`)**; **`+vrfs` (variant read frequency score) is now genuinely ported (previously a registration stub that errored from Init): it piles up the BAM/CRAM alignment list against the FASTA at each indexed site, counts per-sample ref/alt supporting reads, bins the VAF and emits the `SITE`/`MEAN`/`VAR2` profile **byte-for-byte** vs 1.23.1 (streaming + `--use-index`, `--min-depth`, `--nbins` rescaling, `--recalc hc\|data\|file`, `--batch I/N`/`k=N`, `--merge-batches`/`--merge-files`, `-O t\|z`). Parity hinges on vrfs running mpileup2 in `LEGACY_MODE`, whose realignment is an upstream stub — no BAQ, no base-quality filter — so the count reduces to read-flag filtering + a direct CIGAR walk (incl. htslib `is_del`/`indel` column semantics); byte-validated in `native_plugin_vrfs_oracle_test.go` with binary-free `TestUnitVrfs*` helper coverage** | **~99%** | small |
| **samtools** | 25 functional subcommands; CRAM r/w (v2/v3/**v4.0**) + bzip2 encode; `.csi`; consensus `--het-only` + indel calling + **`-a` placeholder rows** (ref-skip/del/zero-cov, live parity) + **`-T/--reference` no-coverage ref-base fill** (pileup/FASTA/FASTQ + `--ref-qual`, live parity); **`stats` command-line positional region arguments** (`in.bam chr1:100-200 chr2`, same overlap filter + target SN lines as `-t`, live parity); `coverage -A`; `markdup -d/-s/-S`; `calmd -C/-e/-u`; `phase` (**deterministic Viterbi DP + `fragphase` chimera repair — no MCMC upstream; full `PS`/`FL`/`M`/`EV` stream byte-identical to upstream across 300 randomized complex/chimera fixtures via the in-place khash kick-out rehash**, incl. **`-b` per-haplotype BAM split with `-A`/`-F`**, byte-exact `dump_aln` routing via in-tree glibc-`drand48` port, live-oracle validated; remaining gap is het *calling* LOD precision at very low `-q`, not phasing); mpileup MAQ **BCF/VCF emit (slices 1–4) + legacy indel + --indels-cns**; **`tview` text/HTML/interactive `-d C`** (byte-for-byte text/HTML, pure-Go termios for `-d C`); remote URL inputs; `-@` threading for view/sort/markdup; **`-@` parallel *input* BGZF decode for view/flagstat/idxstats/stats/depth**; **`-@` parallel *input* CRAM container/slice decode for view/flagstat/idxstats/stats/depth** (records emitted in strict file order — byte-identical for any thread count) | — | **~99%** | — |
| **mosdepth** | single cmd; ALL flags wired incl. `-d/--d4`, `--fragment-mode`, `--quantize`, `-t/--threads`, `--use-median`, `--mapq` fast-path, **CRAM input** (`-f/--fasta`, via alnio); emits `.csi` | none | **~99%** | done |
| **bgzip** | 1/1 cmd, all flags incl. parallel block compression (`-@`/`-t`) + `--test` integrity check | none | **~99%** | done |
| **tabix** | 1/1 cmd, all flags incl. `--reheader`, strict `--targets` post-filter, remote URL region queries | none | **~99%** | done |
| **htsfile** | 1/1 cmd | none (intentional `-c` omission) | **~98%** | done |
| **CRAM / htsgo formats** | CRAM read+write **v2/v3 + v4.0** (uint7-varint + 64-bit pos; v4 transform codecs XPACK/XRLE/XDELTA decode; rANS 4x8/4x16 in-tree, LZMA via ulikunitz/xz); writer output (v3 **and** v4) decodes through upstream samtools (live cross-checks); **lossy read names (`lossy_names=1`) decode — detached mate name read from the mate block + dropped duplicate names reconstructed as `<prefix>:<n>` per htslib `cram_to_bam`**; **writer `=`/`X` CIGAR (`--eqx`) folded to per-base features like M (live `samtools` parity); `B` back-step rejected to match htslib**; BGZF, BAI/CSI, SAM/BAM, BCF v2.2, tabix index | partial-decompress seek via `.gzi` (perf only) | **~95%** | small |

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
parallel BGZF compressor — and for the **input** path: `samtools`
`view`/`flagstat`/`idxstats`/`stats`/`depth` and mosdepth now route BGZF-wrapped
BAM input through `bgzf.MultiReader`, inflating blocks across `-@` worker
goroutines with output that is byte-identical across thread counts. The only
remaining deferral is parallel *input* **CRAM** slice decode (CRAM carries its
own container framing, not BGZF blocks) and bgzip's read path — perf-only, not
per-tool feature gaps.

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
- **RNG byte-parity** — no longer a carve-out. `seqtk sample`/`randbase`
  (glibc drand48 / krand MT), `bedsample` and `bedshuffle`
  (`std::mt19937_64`), and `vcftools --max-indv` (glibc `rand()` +
  `std::random_shuffle`, via the new `--max-indv-seed`) all port the exact
  upstream RNG and are **byte-for-byte identical** to the upstream binary for
  a given seed. The sole genuine exception is that `vcftools --max-indv`
  upstream seeds from `srand(time(NULL))` with no seed flag, so a *plain*
  upstream run is non-reproducible by construction; our seeded path matches
  the algorithm exactly. See the RNG policy section in
  `docs/PARITY_ROADMAP.md`.
- **prinseq PNG report** (`prinseq-graphs.pl` graphics flow) — implemented as
  the `prinseq graph_png` subcommand. PNG byte-identity is N/A (pure-Go
  stdlib renderer, not Perl Cairo/GD); the graph set + plotted data series
  are the asserted parity surface.
- **The bundled upstream bcftools `.so` plugins** — **done.** All 41 in-tree
  plugins are reimplemented in pure Go (`+<name>` dispatch ahead of the
  subprocess fallback), each driven to CLI-to-CLI byte-parity against the real
  upstream `bcftools` 1.23.1 binary, with every previously-"unsupported" mode
  closed (region/target + overlap modes, write-index, multi-threshold filters,
  split-vep, setGT/prune RNG & LD, fill-tags, the stats/format tails, gvcfz,
  frameshifts, trio-dnm3 float models, gtisec polyploid, fixref `--use-id`, and
  the `vrfs` pileup plugin). The VCF-on-stdin/stdout subprocess protocol is
  retained as the fallback for user-supplied executables in any language.

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
