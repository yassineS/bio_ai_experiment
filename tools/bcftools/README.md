# bcftools (pure-Go)

A pure-Go reimplementation of selected [bcftools](https://samtools.github.io/bcftools/)
subcommands. The current implementation ships:

- `bcftools view` — filter / project / convert VCF and BCF. Allele-count /
  frequency selectors (`-c/-C/-q/-Q`) and the private filters (`-x/-X`)
  recompute and append `INFO/AC` and `INFO/AN` like upstream's `calc_ac`
  (suppressed by `-I/--no-update`), even when no sample subset is applied.
- `bcftools query` — format-string output for VCF/BCF records.
- `bcftools concat` — concatenate VCF/BCF files.
- `bcftools norm` — left-align indels, split/join multi-allelics, atomize, dedup.
  The `-m+` join groups biallelics at one position by variant type (or all for
  `-m+any`), builds a common padded REF, and merges FORMAT/GT with upstream's
  rule (a conflicting donor allele goes to the first free strand, so
  `0/1`+`0/1` → `2/1`).
- `bcftools call` — variant calling from per-position genotype likelihoods
  (consensus `-c` + full multi-allelic `-m`).
- `bcftools index` — build a `.csi` (or `.tbi`) index for a BCF / VCF.gz.
- `bcftools stats` — sectioned summary numbers compatible with
  `plot-vcfstats`.
- `bcftools convert` — re-emit VCF/BCF in a different format (`-O v|z|b|u`)
  with optional sample / region filtering, plus the GEN/HAP/TSV conversion
  modes: `--gvcf2vcf` (and its upstream prefix abbreviation `--gvcf`),
  `--tsv2vcf`, `--gensample`/`--gensample2vcf`,
  `--hapsample`/`--hapsample2vcf`, `--haplegendsample`/`--haplegendsample2vcf`,
  and the PLINK exporters `-p`/`--plink` (`.ped`+`.map`), `--tped`
  (`.tped`+`.tfam`) and `--plink-bed` (binary `.bed`+`.bim`+`.fam`).
  See the PLINK export notes below.
- `bcftools mendelian` — detect Mendelian-inconsistent genotypes given
  one or more `CHILD,FATHER,MOTHER` trios. Emits `INFO/MERR` per
  record (annotate mode), a TSV trio rollup (`-c`), or a filtered
  VCF (`-d`). The X-chromosome mode (`-m x`) treats the father as
  haploid on `X`/`chrX`.
- `bcftools mendelian2` — the rewritten upstream plugin. Adds PED-file
  ingestion (`-P/--ped`), `-p/--pfm [1X:|2X:]P,F,M` shortcut, and the
  richer mode bitmask `-m c|[adeEgmMS]` (a=annotate, d=set offending
  trio GTs to ./., e=list err sites, E=drop err sites, g=list good
  sites, m=list missing sites, M=drop missing sites, S=drop skipped
  sites). Multiple letters can be combined. The annotate mode (`-m a`)
  adds the full per-site INFO quartet upstream emits — `MERR` (trios
  with a Mendelian error), `MGOOD` (evaluable, consistent trios),
  `MMISS` (trios with missing/unusable genotypes), and `MNORULE`
  (trios with no applicable inheritance rule) — with the verbatim
  `##INFO` definitions and values, byte-validated against upstream
  1.23.1. The `+mendelian2` plugin form shares this engine.
  *Deliberate divergence (faithful upstream match):* the inheritance-rule
  grammar (`--rules`/`-r`) rejects a per-parent ploidy greater than 2
  ("ploidy > 2 is not supported"), exactly as upstream does — Mendelian
  inheritance is defined only for haploid/diploid genotypes, so there is
  no rule to express for triploid+; this is not a gap.
- `bcftools gtcheck` — sample-identity check by hard-GT Hamming
  concordance (with `-g panel`). Emits the upstream `tsv`-format
  `DC` / `INFO` tables as plain text (`-O t`, default) or BGZF-compressed
  (`-O z`, with an optional `-O z<N>` compression level). Supports
  `--n-matches` (top-N per sample) and `--distinctive-sites`.
- `bcftools roh` — runs-of-autozygosity detection via a 2-state HMM
  (`HW` / `AZ`) over hard-GT input (`-G/--GTs-only`). The ST/RG tables are
  emitted chromosome-major then sample-major (header order), matching
  upstream's per-chromosome synced-reader flush.
- `bcftools merge` — combine multiple per-sample VCF/BCF files. Duplicate
  sample names across inputs error unless `--force-samples` is given, which
  de-dupes by prefixing the clashing name from input *i* with `<i+1>:`
  (e.g. `A + A → A, 2:A`), matching upstream `vcfmerge.c`.
- `bcftools concat` (with `-a/--allow-overlaps`) emits records at a shared
  position in upstream's synced-reader order: grouped by REF>ALT, the groups
  ordered by descending pre-dedup record count, ties by first-appearance.
- `bcftools annotate -x` strips the matching `##FILTER`/`##INFO`/`##FORMAT`
  header lines as well as the column values (bare `FILTER`/`INFO` drop all;
  `FORMAT`/`FMT` keep GT; `FILTER/NAME` rewrites an emptied record to PASS).
- `bcftools filter` — soft-filter records by include / exclude
  expression. Failing records keep their place in the output but
  have FILTER set to the `-s/--soft-filter NAME` (or appended via
  `-m +`); optional `-S/--set-GTs .|0` rewrites failing samples'
  GTs. Supports `-g/--SnpGap` and `-G/--IndelGap` clustering.
- `bcftools consensus` — apply VCF variants (SNPs + simple indels)
  to a reference FASTA. Supports `-s/--samples`, `-H/--haplotype`
  (R/A/I/LR/LA/SR/SA + numeric), `-I/--iupac-codes`, `-m/--mask`
  with `--mask-with`, `--mark-ins / --mark-snv / --mark-del`,
  `-p/--prefix`, `-a/--absent`, `-M/--missing`. When the VCF carries
  samples and neither `-H` nor an allele pick is given, het sites are
  emitted as IUPAC ambiguity codes (upstream's `iupac_GTs` default),
  folded across the `-s` sample or all samples; sites-only VCFs apply
  the first ALT.
- `bcftools polysomy` — estimate chromosomal copy number from BAF.
  Emits a per-sample × per-chromosome TSV
  (`sample / chrom / n_het / mean_baf / median_baf / cn_call`). The
  v1 algorithm is a median-deviation heuristic (CN1 when no hets,
  CN2 when |median - 0.5| ≤ threshold, CN3 otherwise); the full
  Gaussian-mixture peak fit is tracked in
  `docs/PARITY_ROADMAP.md#bcftools`. Reads BAF from FORMAT/BAF
  when present, falls back to FORMAT/AD = REF,ALT.
- `bcftools cnv` — copy-number variation caller. v1 ships a
  heuristic per-sample × per-chromosome median-BAF + mean-LRR
  classifier (CN0..CN4); the full upstream HMM Viterbi is tracked
  in `docs/PARITY_ROADMAP.md`. Output is a TSV with columns
  `sample, chrom, n_sites, median_baf, mean_lrr, cn_call`.
- `bcftools csq` — predict variant consequences against a GFF3
  annotation, via the full haplotype-aware engine. Output is a VCF
  with an `INFO/BCSQ` tag of the form
  `consequence|gene|transcript|biotype|strand|aa_change|dna_change`
  (and a per-sample `FORMAT/BCSQ` bitmask). All upstream output
  containers are supported: `-O v|z|b|u` (VCF / BGZF-VCF / BCF /
  uncompressed-BCF) plus `-O t`, the streaming tab-delimited text
  form (upstream `FT_TAB_TEXT`) that emits one
  `CSQ<TAB>sample<TAB>haplotype<TAB>chrom<TAB>pos<TAB>consequence`
  row per (sample, haplotype) consequence, byte-for-byte with
  upstream `text_print_vcsq` (the leading `#`-comment version/command
  provenance lines aside). Remaining gaps (e.g. `-s -` sample
  dropping) are tracked in `docs/PARITY_ROADMAP.md`.
- `bcftools mpileup` — per-position genotype likelihoods from BAM
  input, the upstream input to `bcftools call`. v1 ships the
  SNP-only, uniform-error binomial likelihood model: each base
  contributes log10(P(base | g)) for g ∈ {0/0, 0/1, 1/1} based on
  e = 10^(-Q/10), then sums across reads and emits Phred-scaled
  `FORMAT/PL`. Reads `pkg/htsgo/sam` for BAM ingestion and
  `pkg/htsgo/fasta` for the reference. Output is a streaming
  VCF with `INFO/DP`, `INFO/I16`, and biallelic `FORMAT/PL`. BAQ
  recalibration, indel calling, and the full MAQ likelihood model
  are tracked in `docs/PARITY_ROADMAP.md`.
- `bcftools som` — Self-Organizing Map (Kohonen map) variant
  classifier. `--train` reads a VCF/BCF, extracts the per-site INFO
  annotation vector (`-t/--training-annots`, default
  `QUAL,MQ,MQ0F,BQB,MQB,RPB,SGB`, min/max-normalised onto [0,1]),
  trains a 2-D map of the given `-s/--size` (default 20) with the usual
  learning-rate / radius decay, and writes a usable `<prefix>.som`.
  `--classify` loads that map and prints one SOM score per site
  (`1 - dist/sqrt(kdim)`, higher = more training-like). The on-disk map
  format is **our own** clean, versioned binary format (magic `SOMGO1`)
  because upstream's `SOMv1` layout is unusable: upstream's
  `som_write_map` (`vcfsom.c:170`) has an `fwrite`-return bug that aborts
  `--train` after truncating the map to 5 bytes (see
  `docs/UPSTREAM_BUGS.md#bcftools-som-write-map`); this port fixes it so
  the train→classify pipeline works. Upstream's `-f/--nfold` /
  `-m/--merge` cross-validation knobs are accepted as surface only (v1
  trains a single map).
- `bcftools plugin` — run a user-supplied plugin over a VCF/BCF, also
  reachable as the `bcftools +<name>` shorthand. Unlike upstream — which
  loads plugins as native shared objects via `dlopen` — this port runs a
  plugin as an ordinary **child process**: it streams the input VCF as
  uncompressed text to the plugin's stdin and reads the plugin's stdout
  back as VCF, so a plugin is "a filter from VCF on stdin to VCF on
  stdout" and can be written in any language. Plugins are discovered via
  the `BCFTOOLS_PLUGINS` colon-separated directory list, exactly as
  upstream. `bcftools plugin -l` (`-lv` for verbose) lists them. The
  host applies `-o`/`-O` formatting around the plugin's output. The
  ~30 upstream bundled plugins are intentionally **not** ported — only
  the mechanism is, so users can write their own. The full contract is
  in `docs/PLUGIN_PROTOCOL.md`, with a reference example plugin under
  `tools/bcftools/plugins/example/`.
- `pkg/htsgo/bcf` — reader and writer for the BCF v2.2 binary format.

All pieces share the existing `pkg/htsgo/vcf` types so downstream
consumers see records as familiar `vcf.Variant` values.

## `bcftools plugin` / `bcftools +<name>`

A plugin is any executable found by name in one of the colon-separated
directories in the `BCFTOOLS_PLUGINS` environment variable. The host
pipes uncompressed VCF text to the plugin's stdin and reads VCF text from
its stdout; plugin arguments are passed as the plugin's `argv`, a literal
`--` separates them from the host input file. A non-zero plugin exit is
surfaced as an error with the plugin's stderr. See
`docs/PLUGIN_PROTOCOL.md` for the complete specification.

```bash
# Discover plugins.
export BCFTOOLS_PLUGINS=/path/to/plugins
bcftools plugin -l            # names, one per line
bcftools plugin -lv           # path + --about description

# Run a plugin (these two forms are equivalent).
bcftools +example   -- input.vcf.gz
bcftools plugin example -O z -o out.vcf.gz -- input.vcf.gz chr1:1-1000
```

### Region / target selection for native plugins (`-r/-R/-t/-T`)

The native (pure-Go) plugins honour the same `-r/--regions`,
`-R/--regions-file`, `-t/--targets` and `-T/--targets-file` options upstream
does, applied by a shared host-side filter before any record reaches the
plugin:

- `-r`/`-R` is **span-overlap** based — a record is kept if `[POS,
  POS+len(REF)-1]` overlaps the region.
- `-t`/`-T` is **record-start** based — a record is kept if `POS` falls in the
  target window; a leading `^` negates (exclude matches). This is upstream's
  exact `-r` vs `-t` difference: an indel at `POS=100` spanning `100..104` is
  kept by `-r chr:102-102` but dropped by `-t chr:102-102`.
- `-R`/`-T` files follow the htslib synced-reader format: a `.bed` path is
  0-based half-open; any other file is 1-based (`chr<TAB>pos` or
  `chr<TAB>beg<TAB>end`). Inline `-r`/`-t` strings keep the `chr:beg-end` colon
  syntax.

Honoured by `+check-sparsity`, `+remove-overlaps`, `+prune`, `+smpl-stats`,
`+indel-stats`, `+contrast`, `+guess-ploidy` (only `-r/-R`; its `-t` is
`--tag`), `+mendelian2`, `+trio-stats`, `+isecGT` (applied to both inputs),
`+split` and `+scatter`. `+check-sparsity` additionally groups and labels its
report per region and — unlike upstream, which silently drops a BED/TSV `-R`
line (see `docs/UPSTREAM_BUGS.md`) — accepts BED/TSV `-R` files as a fix-on-port
while keeping colon region-list lines byte-identical to upstream.

### Curly-brace multi-threshold `-i/-e` expansion (stats plugins)

`+smpl-stats`, `+indel-stats` and `+trio-stats` support upstream's curly-brace
multi-threshold filter syntax: an `-i`/`-e` expression may contain a
`{a,b,c}` list, which is expanded into one concrete filter per element (the
braces replaced by that element), each tallied into its own `FLT*`/`SITE*`
(and, for indel-stats, `SN*`/`DVAF*`/`DLEN*`/`DFRAC*`/`NFRAC*`) report section.
For example `bcftools +smpl-stats -i 'FMT/GQ>{10,20,30}' file.bcf` reports the
per-sample stats three times, once at each GQ threshold. Multiple `{...}`
groups combine as a cartesian product, in upstream's exact ordering; an empty
list (`{}`) collapses to the single default "all" filter, and an unmatched `{`
is a hard parse error — all byte-validated against upstream 1.23.1.

### Stats / contrast plugin remaining flags (`-o`, `-p`, `-a`, `-f`)

The four stats/contrast plugins now match upstream on every flag:

- **`-o`/`--output FILE`** (`+smpl-stats`, `+indel-stats`, `+trio-stats`) writes
  the report to FILE instead of stdout; the bytes are identical to the stdout
  form (the `CMD` line echoes the verbatim argv in both).
- **`+indel-stats -p`/`--ped FILE`** restricts the stats to de-novo indels in
  each PED trio's child (with `--alt2ref-DNM` widening the DNM definition); the
  SN* "number of samples" column then reports the trio count. (Fix-on-port:
  upstream aborts on a PED indel VCF without FORMAT/AD; our port reports the
  AD-independent counts anyway — see `docs/UPSTREAM_BUGS.md`.)
- **`+trio-stats -a`/`--alt-trios INT`** applies the deferred singleton/doubleton
  transmission-rate accounting: an allele counts only when present in at most
  `INT` alternate trios at the site (`0` = unlimited, the default).
- **`+contrast -f`/`--max-allele-freq NUM`** adds rare-allele enrichment: the
  per-site VCF output is unchanged, and a second stderr summary line
  `max_AC/PASSOC/FASSOC/NASSOC:` reports the region-wide Fisher's exact
  probability and control/case non-REF fractions over the pooled minor alleles
  (an integer `NUM` is an allele-count threshold; a float in `[0,1]` is an
  allele-frequency threshold scaled by the sample count). The
  `--regions-overlap`/`--targets-overlap` region-matching modes remain
  unsupported.

All byte-validated against upstream 1.23.1.

### Format / output-mode tails (`+ad-bias`, `+remove-overlaps`, `+tag2tag`, `+guess-ploidy`, `+af-dist`)

The remaining per-plugin format/output modes now match upstream 1.23.1:

- **`+ad-bias --clean-vcf`/`-c`** emits the VCF subset to only the ALT alleles
  whose Fisher p-value passes `-t` (dropping sites where nothing passes), with
  `AC`/`AD`/`PL`/`GT` and other Number=A/R/G fields remapped via a faithful port
  of htslib's `bcf_remove_allele_set`. **`+ad-bias -f`/`--format`** appends a
  `bcftools query`-style format column (evaluated once per record) to every `FT`
  report line; `-f` and `-c` are mutually exclusive, as upstream. (The `-c`
  short form takes no argument; the `--clean-vcf` long form consumes and ignores
  one, reproducing an upstream `getopt_long` quirk.)
- **`+remove-overlaps -m 'min(QUAL)' --missing DP`** resolves overlaps using a
  coverage heuristic for missing-QUAL records (scale the window's maximum QUAL by
  `INFO/DP`: `maxQUAL*DP/maxQUAL_DP`); `--missing 0` is the explicit scalar
  default. **`-Ot`/`-Otz`** emits a plain (or bgzip-framed) `chr<TAB>pos` list
  instead of the VCF. (Fix-on-port: a deletion-window-boundary + ring-wrap corner
  where upstream leaks a stale overlap mark across windows is corrected — see
  `docs/UPSTREAM_BUGS.md#bcftools-remove-overlaps-minqual-stale-mark`.)
- **`+tag2tag --LXX-to-XX`** (and the partial `--LPL-to-PL` / `--LAD-to-AD`)
  expand the localized FORMAT tags (`LAA` + `LPL`/`LAD`) back into the standard
  Number=G `PL` and Number=R `AD`, mapping each sample's localized indices via
  `FORMAT/LAA`. `-d`/`--defaults` supplies the value for untouched cells and
  `-s`/`--skip-nalt` skips sites above an allele threshold. The reverse
  direction (`--XX-to-LXX`) is a `todo` upstream too and is rejected with the
  same restriction.
- **`+guess-ploidy -g`/`--genome`** is the `b37`/`b38`/`hg19`/`hg38` shortcut for
  the non-PAR chrX region, expanded to the equivalent `-r CHR:BEG-END`
  (`X:2699521-154931043` etc.) before the shared region filter runs.
- **`+af-dist -p`/`-d`** bin lists may be read from a file (one boundary per
  line) as well as inline (a comma-separated list), exactly as upstream's
  `bin_init`/`hts_readlist` decides.

All byte-validated against upstream 1.23.1.

`-W`/`--write-index` (bare, or `=csi`/`=tbi`) writes a CSI (default) or TBI
index next to each indexable (`.vcf.gz`/`.bcf`) output; plain-VCF and stdout
outputs are non-indexable and error exactly as upstream does. Honoured by
`+contrast`, `+isecGT`, `+mendelian2`, `+split` and `+scatter` (multi-output
plugins index every file). The produced index is byte-validated against
`bcftools index` over the same data.

### `+gvcfz` — resize gVCF blocks (native)

`gvcfz` is implemented natively (pure Go). It groups consecutive gVCF reference
blocks (records whose only ALT is `<NON_REF>`/`<*>`) by the
`-g/--group-by FILTER:EXPR[; FILTER:EXPR …]` clauses, where each `EXPR` is a full
bcftools filter expression evaluated per record (FORMAT/GT predicates included,
e.g. `GQ>60 & DP<20`, `GT!="alt"`; a `-` expression is the catch-all). The first
matching group selects the block a record belongs to; consecutive same-group
records merge into the first (representative) record, with INFO/END extended to
the block end, FORMAT/DP set to the minimum MIN_DP (or DP), FORMAT/GQ (or RGQ) to
the minimum, and FORMAT/PL to the element-wise minimum. A real variant flushes the
current block and passes through verbatim. A non-PASS group label adds a
`##FILTER` line whose Description is the verbatim `-g` string (with `"` rewritten
to `'`). `-i/-e` apply a record-level pre-filter; `-a/--trim-alt-alleles` is
accepted; `-o/-O/-W` are handled by the host pipeline.

```bash
# Resize blocks by GQ and DP; non-PASS labels mark the block FILTER.
bcftools +gvcfz input.bcf -g'PASS:GQ>60 & DP<20; PASS:GQ>40 & DP<15; Flt1:GQ>20; Flt2:-'
# Collapse all non-variant sites into one block, removing unused ALTs.
bcftools +gvcfz input.bcf -a -g'PASS:GT!="alt"'
```

### `+fixref` — fix REF strand orientation (native)

`+fixref` determines and fixes REF/ALT strand orientation against a FASTA
reference. All upstream modes are ported: `stats` (default; collect+print the
strand-convention stats, no VCF output), `ref-alt`, `swap`, `flip`, `flip-all`,
`top` (Illumina TOP → fwd with ambiguous-pair sequence walking), and **`id` /
`--use-id`** — which determines the correct REF allele from a separate dbSNP
VCF keyed by the **ID (rsID) column** instead of from strand convention. Each
converted record is annotated with `INFO/FIXREF` (`-t` renames the tag)
recording the change (`none`/`swap`/`flip`/`GT`/`skip`/`err`).

```bash
# Match REF/ALT to a dbSNP VCF by rsID; discard sites with no dbSNP match.
bcftools +fixref input.vcf -- -f ref.fa -i dbsnp.vcf.gz -d
# Equivalent: -m id with an explicit dbSNP file.
bcftools +fixref input.vcf -- -f ref.fa -m id -i dbsnp.vcf.gz
```

In `id` mode each input record's ID is looked up in a per-chromosome
rsID→{position, ref-base} map built from the dbSNP file (skipping non-SNPs,
non-`[ACGT]` REF and missing `.` IDs; the first record wins on a duplicate ID).
If the input REF already equals the dbSNP REF the site is left unchanged
(`none`); if the input ALT equals the dbSNP REF, REF/ALT are swapped and every
sample genotype is flipped (`swap`); a missing/unknown ID or neither-allele
match is left unresolved (`skip`, or dropped with `-d/--discard`). When the
dbSNP record sits at a different position the input position is corrected (and a
`fixed pos` is counted). Both the corrected VCF and the end-of-run stats summary
match upstream 1.23.1 byte-for-byte. Upstream requires the dbSNP file to be
bgzip-compressed and tabix/CSI-indexed; this port additionally accepts a plain
(un-indexed) `.vcf`/`.vcf.gz` as a fix-on-port robustness superset (see
`docs/UPSTREAM_BUGS.md#bcftools-fixref-id-plain-vcf`).

### `+frameshifts` — annotate frameshift indels (native)

`frameshifts` is implemented natively (pure Go). It reads exons from
`-e/--exons FILE` (a BED or region-list, optionally bgzipped+tabixed) into an
in-tree port of htslib's `bcf_sr_regions_overlap` cursor and adds INFO/OOF (one
Integer per ALT allele) to every indel record that overlaps an exon.

By default it reproduces the **shipped upstream behaviour byte-for-byte**: the
plugin's per-allele in-frame/out-of-frame computation is dead code in the real
binary (the `var[i].type != VCF_INDEL` guard is always true under modern htslib),
so every exon-overlapping indel allele is annotated `OOF=-1` ("not applicable")
and the intended length-mod-3 result is never produced — see
`docs/UPSTREAM_BUGS.md#bcftools-frameshifts-oof-dead-code`. The **corrected**
computation (trim the inserted/deleted length against the exon, then take it
mod 3: out-of-frame ⇒ 1, in-frame ⇒ 0) is available via the opt-in `--fix-oof`
flag, which deviates from drop-in parity on purpose.

```bash
bcftools +frameshifts in.vcf -- -e exons.bed.gz          # upstream-exact OOF=-1
bcftools +frameshifts in.vcf -- -e exons.bed.gz --fix-oof # corrected in/out-of-frame
```

Both are byte-validated against upstream 1.23.1 (plain-BED and tabixed-BED cursor
paths), with binary-free `TestUnit*` coverage of the pure helpers.

### `+prune` — prune/annotate by LD or window density (native)

`prune` is implemented natively (pure Go, no subprocess) with **all** upstream
modes, byte-validated against bcftools 1.23.1:

- `-n/--nsites-per-win N` keeps at most `N` sites per `-w` window, selecting by
  `-N maxAF` (default — biggest allele frequency from `--AF-tag` or, without it,
  from INFO/AC+AN or the genotypes), `-N 1st` (first encountered), or `-N rand`
  (random; `--random-seed INT` makes it reproducible).
- `-m count=N` removes clusters of more than `N` sites within the window;
  `-m R2=/LD=/RD=FLOAT` (or a bare number == r2) discards sites whose linkage
  disequilibrium with a kept upstream site exceeds the threshold. The three LD
  measures (correlation r², Lewontin's D', Ragsdale's RD) are computed exactly
  as upstream's `calc_ld` (`+ - * /` and `sqrt` only, so byte-identical after
  htslib's float32 narrowing).
- `-f LABEL` soft-filters instead of discarding (sets the FILTER column,
  requires `-m`); `-a count|r2|LD|RD` annotates each site with the cluster size
  or the maximum LD value and the partner site's position (`R2`/`POS_R2`, …).
- `-w INT[bp|kb|Mb]` sets a bp window (suffix) or a site-count window (bare
  integer); `-k/--keep-sites` leaves `-i`/`-e`-filtered sites in place;
  `--randomize-missing` fills missing genotypes from the site allele frequency
  via the same deterministic drand48 stream.

Two upstream quirks are reproduced for parity: `maxAF` ranks by alt/ref (so a
monomorphic-ALT site sorts *lowest*) and the soft-filter header renders a bp
window as "within 0kb"; see
`docs/UPSTREAM_BUGS.md#bcftools-prune-maxaf-ranks-by-altref-not-allele-frequency`.
The oracle lives in `native_plugin_prune_oracle_test.go`.

### `+trio-dnm3` — de-novo mutation screening (native)

`trio-dnm3` is implemented natively (pure Go, no subprocess). It screens
trios for de-novo mutations and writes `FORMAT/DNM` (the score),
`FORMAT/VA` (the de-novo allele) and `FORMAT/VAF` (percent ALT reads).
All four upstream models are supported:

- `--use-NAIVE` — GT-only Mendelian-incompatibility flag. Pure integer
  table lookup; **byte-exact** vs upstream (`DNM`/`VA` are integers).
- `--use-DNG` — the original DeNovoGear likelihood over `FORMAT/PL`
  (implies `--dng-priors`).
- `--use-ALM` — the allele-likelihood model over `FORMAT/QS` (or `PL`
  with `--with-pPL`, or fake-QS-from-AD with `--with-pAD`).
- `--use-DMM` — the default Dirichlet-multinomial model over
  `FORMAT/AD`+`QM` (and `PL`, unless `--with-cAD`).

The supported knobs match upstream: `--dnm-tag TAG[:log|phred|prob|flag]`,
`--va`, `--vaf`, `-n/--strictly-novel`, `--mrate`, `--pn`/`--pns`,
`--phi`, `--max-QM`, `--min-vaf`, `--noise-prior`/`--np`, `--strand-bias`/`--sb`,
`--allelic-dropout`/`--ad`, `-X/--chrX` (GRCh37/GRCh38 PAR presets), `-m/--min-score`,
the `-i`/`-e` per-trio filters, and `>4`-allele trimming.

**libm-tolerance boundary.** The NAIVE verdict is integer-exact. The
DMM/ALM/DNG **scores** are long `log`/`exp`/`pow`/`lgamma` reductions; the
incomplete-beta and log-gamma kernels go through the bit-stable in-tree
`kfBetai`/`kfLgamma` port (the same AS245 code upstream's `kfunc.c` uses),
while the remaining transcendentals use Go's `math`. Because libm
transcendentals are only guaranteed to the last ULP, the de-novo score may
differ from the C build in the last printed digit (e.g. `-46.0521` vs
`-46.0522`) after htslib narrows it to a 32-bit float and prints it with
`%g`. That is the floating-point reproducibility boundary, not a bug:
byte parity is **not** the contract for the float scores. Parity is
asserted with a field-aware, tolerance-aware comparison (string fields
exact; numeric `DNM`/`VA`/`VAF` fields equal within ~6 significant figures
or a small relative/absolute epsilon — see `numeric_parity_test.go`). On
linux/amd64 the scores in fact land byte-for-byte.

### `+split-vep` — query VEP/CSQ subfields (native)

`split-vep` is implemented natively (pure Go, no subprocess). It splits a
structured `INFO/CSQ` (Ensembl VEP) or `INFO/BCSQ` (`bcftools csq`)
annotation into its pipe-delimited subfields and supports the **full**
upstream surface:

- **`-c/--columns`** extracts named/indexed subfields (ranges, `:TYPE`
  suffixes, `-` for all) into new `INFO` tags; **`-f/--format`** prints a
  `bcftools query`-style line; `-A/--all-fields` expands `%CSQ`; `-d`
  duplicates per transcript; `-p/--annot-prefix`, `-x/-X`, `-u`.
- **`-s/--select TR:CSQ:PRN`** transcript and consequence selection:
  - `TR` — `all`, `worst` (most severe, see `-S`), `primary`
    (`CANONICAL=YES`), `pick` (`PICK=1`), `mane` (`MANE_SELECT!=""`), or an
    arbitrary `<FIELD><OP><VALUE>` EXPRESSION with the `=`, `!=`, `~`, `!~`
    operators (the value may be double-quoted).
  - `CSQ` — `any` or a severity term with the `+` (this-or-more-severe) / `-`
    (this-or-less-severe) modifiers. Terms are matched case-sensitively
    against the lowercased scale (so `:MISSENSE` is rejected, exactly as
    upstream).
  - `PRN` — `all` (print every consequence term) or `worst` (rewrite the
    printed Consequence to its single worst `&`-joined term).
- **`-g/--gene-list [+]FILE`** restricts to transcripts whose gene appears in
  `FILE` (one gene per line); a leading `+` *prioritises* instead — all
  transcripts are kept but the listed-gene ones are moved to the front.
  `--gene-list-fields LIST` chooses which subfields are matched (default
  `SYMBOL,Gene,gene`).
- **`-S/--severity -|FILE`** overrides the built-in consequence severity
  scale from a file (one tier per line, whitespace-separated synonyms);
  `-S -`/`-S ?` print the default scale.
- **`--columns-types -|FILE`** overrides the auto-detected column types via a
  regex→type table (each pattern is anchored `^…$`); `-` prints the default
  table. It drives both the emitted `##INFO` `Type` and the numeric
  re-parsing of values.
- **`-i/--include` / `-e/--exclude`** evaluate over the derived
  per-transcript CSQ subfields (auto-registered as `INFO` tags), not plain
  VCF fields.

Every mode is byte-validated against the upstream binary 1.23.1 via the
CLI-to-CLI oracle (`TestNativePluginSplitVep*` in
`native_plugin_batch7_oracle_test.go`). One upstream quirk is preserved
deliberately: `PRN :worst` ranks `&`-joined terms by an *exact* scale lookup
(not the substring matcher), so a compound term whose parts aren't exact
scale tokens keeps its first term — see
`docs/UPSTREAM_BUGS.md#bcftools-split-vep-prn-worst-exact-match`.

### `+setGT` — set genotypes (native)

The native `+setGT` port now covers the **full** upstream target/new-gt grammar,
including the modes that used to be rejected:

- **`-t b:TAG CMP VAL`** — set diploid heterozygous genotypes whose two-tailed
  binomial test over a FORMAT integer tag (typically `AD`) satisfies the
  comparison. `CMP` is one of `<`, `<=`, `>`, `>=`, `==`/`=`; the p-value is
  `binom.test(nAlt, nRef+nAlt, 0.5)` computed via the same regularized
  incomplete-beta function htslib uses, so it is bit-exact. Example:
  `bcftools +setGT in.vcf -- -t 'b:AD<1e-3' -n 0`.
- **`-t r:FLOAT` with `-s/--seed INT`** — act on a random proportion `FLOAT`
  (0<FLOAT<1) of the targeted genotypes. Upstream seeds htslib's *deterministic*
  drand48 PRNG from `-s` (default 0) and nothing else, so the result is fully
  reproducible; the port reimplements the same 48-bit LCG and is byte-identical
  to upstream for any fixed seed (and across thread counts). Used alone, `-t r`
  implicitly targets all genotypes.
- **`-n X`** — set every allele of the genotype to the allele with the largest
  FORMAT/AD value for that sample (also usable inside a `c:` template, e.g.
  `-n c:0/X`). Requires a FORMAT/AD header.

All three modes — plus every operator and several new-gt targets — are
byte-validated against the upstream binary 1.23.1 via the CLI-to-CLI oracle
(`TestNativePluginSetGTBinom`, `TestNativePluginSetGTRandom`,
`TestNativePluginSetGTReadDepth`), and the drand48 port is pinned to the
canonical POSIX sequence by `TestDrand48KnownVectors`.

### `+fill-tags` — recompute INFO/FORMAT annotations (native)

The native `+fill-tags` port covers the **full** upstream surface:

- **Every built-in tag.** `AN`, `AC`, `AC_Hom`, `AC_Het`, `AC_Hemi`, `AF`,
  `MAF`, `NS`, `HWE`, `ExcHet`, `END`, `TYPE`, `FORMAT/VAF`, `FORMAT/VAF1`, and
  the `F_MISSING` expression. Counts follow upstream's exact het/hom/hemi/half
  classification (including `-d/--drop-missing`); `HWE`/`ExcHet` use the in-tree
  Wigginton (PMID:15789306) exact test, and `AF` falls back to `INFO/AN,AC` for
  sites-only records.
- **`-t LIST`** selection with `INFO/`/`FORMAT/` qualifiers and the `all`
  keyword (which excludes `END`/`TYPE`, as upstream).
- **`-S/--samples-file FILE`** population grouping. The file lists
  `SAMPLE  GROUP[,GROUP2,...]` per line; each distinct group gets its own
  `_GROUP`-suffixed tags (and `## ... in GROUP` headers), alongside the summary
  `ALL` population. Example:
  `bcftools +fill-tags in.vcf -- -S groups.txt -t AN,AC,AF,HWE`.
- **Custom expression `TAG[:Number]=[int|integer|float](EXPR)`.** An in-tree
  evaluator supports INFO/FORMAT tag references, the aggregations
  `sum`/`avg`(`mean`)/`max`/`min`/`median`/`stdev` and their per-sample
  `smpl_*` (`sSUM`, `sMEAN`, ...) variants, arithmetic `+ - * /`, unary minus,
  `abs`, `phred`, and the genotype reductions `F_MISSING`/`N_MISSING`/
  `F_PASS(COND)`/`N_PASS(COND)`. `int()`/`integer()` produces an Integer field
  (C `round()` half-away), `float()`/bare produces Float; `:Number` fixes the
  count, otherwise `Number=.`. Examples:
  `-t 'DP:1=int(sum(FORMAT/DP))'`, `-t 'FORMAT/VD:1=int(smpl_sum(FORMAT/AD))'`,
  `-t 'good=N_PASS(GT="het")'`.
- **`-l/--list-tags`** prints the available-tag table to stderr and exits.

Byte-validated against the upstream binary 1.23.1 via the CLI-to-CLI oracle
(`TestNativePluginFillTagsPops`, `TestNativePluginFillTagsListTags`,
`TestNativePluginFillTags`), with binary-free `TestUnitFillTags*` unit tests for
the pure helpers (formula calculators, tag-list/samples-file parsers, expression
evaluators).

### `+vrfs` — variant read frequency score (native)

`vrfs` is implemented natively (pure Go). It assesses site noisiness from a
large set of unaffected alignments: given an alignment list (`-a/--alns`, one
BAM/CRAM path per line), a FASTA reference (`-f/--fasta-ref`) and a
tab-delimited sites file (`-s/--sites`, `chr pos ref alt`), it piles up every
alignment at each indexed site, counts the per-sample ref/alt supporting reads,
bins the variant-allele fraction into a per-site histogram, and emits the
`SITE`/`MEAN`/`VAR2` profile.

```bash
# Typical run (streaming the BAMs)
bcftools +vrfs -f ref.fa -a bams.txt -s sites.txt -o scores.txt

# Use the BAM index to jump to the sites (faster with few sites)
bcftools +vrfs -f ref.fa -a bams.txt -s sites.txt -i -o scores.txt

# Batch a large list, then merge
bcftools +vrfs -a bams.txt --batch k=3                       # prints the batch count
bcftools +vrfs -f ref.fa -a bams.txt -s sites.txt --batch 1/2 -o s1.txt
bcftools +vrfs -f ref.fa -a bams.txt -s sites.txt --batch 2/2 -o s2.txt
bcftools +vrfs --merge-batches list.txt                       # list.txt holds s1.txt, s2.txt
```

Supported options: `-d/--min-depth`, `-n/--nbins` (with hard-coded-profile
rescaling), `-r/--recalc hc|data|file:PATH`, `-i/--use-index`, `-b/--batch
I/N` and `k=N`, `-m/--merge-batches`/`-M/--merge-files`, and
`-o`/`-O t|z[0-9]` output.

**Parity boundary (byte-for-byte).** vrfs runs the htslib mpileup2 engine in
`LEGACY_MODE`, whose realignment step is **stubbed out** upstream
(`mpileup2/mpileup.c`): there is **no BAQ and no base-quality adjustment**, and
the vrfs count loop never reads base quality. The pileup therefore reduces to
read-level flag filtering (drop unmapped/secondary/qcfail/dup; no MAPQ floor, no
supplementary drop, no orphan filtering) plus a direct CIGAR walk. Because no
ambiguous indel realignment is ever performed, the output is **byte-exact** vs
upstream 1.23.1 across SNV *and* indel sites — there is no proximity tolerance
in play. The CIGAR walk reproduces htslib `bam_pileup1_t` `is_del`/`indel`
column semantics (a deletion column reads the post-deletion base; the last
aligned base before an I/D op is classified as a generic indel). The
empty-profile `MEAN` line is emitted as `-nan` to match glibc `printf` of
`0.0/0`. See `docs/UPSTREAM_BUGS.md#bcftools-vrfs-legacy-mode-no-realign` for
the boundary detail.

Byte-validated against the upstream binary 1.23.1 via the CLI-to-CLI oracle
(`native_plugin_vrfs_oracle_test.go`: 10 profiling cases + 2 merge cases over
upstream-samtools-built BAM fixtures), with binary-free `TestUnitVrfs*` unit
tests covering the VAF binning, the sites/aln-list parsers, the per-(sample,
site) accumulator, the profile mean/var aggregation, the variance rescaling, the
C-style float formatting and the CIGAR→column walk.

## Quick start

```bash
# Build everything from the repo root.
go build ./tools/bcftools/cmd/bcftools

# Pass-through (VCF → VCF) — defaults to stdout, format -O v.
./bcftools view input.vcf

# Keep only PASS records, drop the per-sample columns, write gzipped VCF.
./bcftools view -f PASS -G -O z -o filtered.vcf.gz input.vcf

# Restrict samples and apply an expression.
./bcftools view -s sample1,sample3 -i 'INFO/DP>30 && FILTER="PASS"' input.vcf

# Region query, .tbi-backed when the sibling index exists.
./bcftools view -r chr1:100-200 input.vcf.gz

# Decode BCF and emit VCF text.
./bcftools view input.bcf

# Summary stats — emits the same tab-prefixed multi-section format that
# `plot-vcfstats` expects.
./bcftools stats input.vcf
./bcftools stats -s sample1,sample3 -i 'FILTER="PASS"' input.vcf
./bcftools stats -r chr1:100-200 -d 0,100,10 input.vcf
```

## POSIX short-flag handling (all subcommands)

Every subcommand parses its arguments through `pkg/cliflag`'s
`cliflag.Parse`, so it accepts the same getopt-style short-flag forms as
upstream bcftools in addition to GNU long flags:

- **Bundling**: `view -hG` is equivalent to `view -h -G`; boolean short
  flags can be clustered (`call -mv` == `call -m -v`).
- **Value concatenation**: `-Ob` == `-O b`, `norm -m-` == `norm -m -`,
  `stats -s-` == `stats -s -`.
- `--` ends option parsing; a bare `-` means stdin/stdout.

A handful of upstream legacy/compat short flags are also accepted so old
command lines (including bundled ones) keep working: `norm -D` (alias of
`-d exact`), `call -f` (alias of `-a`/`--annotate`), `call -Y` (deprecated
alias of `--ploidy Y`), `call -N` (omit-REF-N, the default), and
`plugin -lv`/`-l -v` (verbose plugin listing). Like upstream getopt, the
short `-O=v` *equals* form is not accepted; use `-O v`, `-Ov`, or the long
`--output-type=v` form.

## Supported flags

| Short | Long                  | Meaning |
| ----- | --------------------- | ------- |
| `-O`  | `--output-type`       | `v` (default), `z` (BGZF VCF), `u` (uncompressed BCF), `b` (BGZF BCF). All four output types are implemented (BCF via the in-tree `pkg/htsgo/bcf` writer). |
| `-o`  | `--output`            | Output path (defaults to stdout). |
| `-h`  | `--header-only`       | Emit only the header. Help is on `-?` / `--help`. |
| `-H`  | `--no-header`         | Drop the header. |
| `-G`  | `--drop-genotypes`    | Strip the FORMAT and per-sample columns. |
| `-c`  | `--min-ac`            | Minimum non-reference allele count. |
| `-C`  | `--max-ac`            | Maximum non-reference allele count. |
| `-q`  | `--min-af`            | Minimum allele frequency. |
| `-Q`  | `--max-af`            | Maximum allele frequency. |
| `-i`  | `--include`           | Keep records matching expression. |
| `-e`  | `--exclude`           | Drop records matching expression. |
| `-f`  | `--apply-filters`     | Comma list of FILTER names to keep. |
| `-r`  | `--regions`           | Region list (`chr:beg-end[,...]`). Uses `.tbi` if present. |
| `-R`  | `--regions-file`      | BED-like regions file (`CHROM \t BEG \t END`). |
| `-t`  | `--targets`           | Like `-r` but always a post-filter (no index needed). |
| `-T`  | `--targets-file`      | BED-like targets file (post-filter). |
| `-s`  | `--samples`           | Restrict samples to this comma list. |
| `-S`  | `--samples-file`      | File of sample IDs (one per line). |
| `-x`  | `--private`           | Keep only sites whose non-reference alleles are exclusive (private) to the subset samples. Requires a subset (`-s`/`-S`). |
| `-X`  | `--exclude-private`   | Inverse of `-x`: drop sites private to the subset samples. |
| `-l`  | `--compression-level` | gzip level for `-O z` output. |
|       | `--threads`           | Accepted; v1 is single-threaded. |
| `-?`  | `--help`              | Show help. |
|       | `--version`           | Show version. |

### Note on `-h` semantics

Upstream `bcftools` uses `-h` for "header-only" instead of "help" (a divergence
from the rest of the GNU world). This implementation follows the upstream
convention so existing scripts keep working. Help is exposed on `-?` and
`--help`.

## Expression language (`-i` / `-e`)

The expression evaluator is a small recursive-descent parser. The supported
grammar is:

```text
expr     := or_expr
or_expr  := and_expr ("||" and_expr)*
and_expr := unary ("&&" unary)*
unary    := "!" unary | primary
primary  := "(" expr ")" | comparison | value
comparison := value op value
value    := "INFO/" IDENT | "FILTER" | NUMBER | STRING | IDENT
op       := "==" | "=" | "!=" | "<" | "<=" | ">" | ">="
```

Examples:

```bash
bcftools view -i 'INFO/DP > 30'
bcftools view -i 'INFO/DP > 30 && FILTER="PASS"'
bcftools view -i 'INFO/AF >= 0.05 || FILTER="LowQual"'
bcftools view -e 'INFO/MQ < 40'
bcftools view -e '!(INFO/H2)'
```

INFO field values are coerced to numbers when both sides of a comparison parse
as such; otherwise the comparison is lexical. Multi-value INFO entries take
the first comma-separated element (matching the upstream default).

## `bcftools stats`

`bcftools stats` emits the same nine tab-prefixed sections that the upstream
binary produces. Each section starts with a `# <ID>, <description>` comment,
the column-header comment, and then data rows whose first column is the
section short name:

| Section | Meaning |
| ------- | ------- |
| `SN`    | Summary numbers — record / SNP / MNP / indel / multi-allele totals. |
| `AF`    | Counts binned by non-reference allele frequency. |
| `QUAL`  | Counts binned by `QUAL`. |
| `IDD`   | Indel-length distribution. |
| `ST`    | Substitution-type counts (A>C, A>G, …). |
| `DP`    | Depth distribution (sites and per-sample GTs). |
| `PSC`   | Per-sample counts (RefHom / NonRefHom / Hets / Ts / Tv / Indels / avgDP). |
| `PSI`   | Per-sample indel counts. |
| `HWE`   | Hardy-Weinberg-equilibrium statistic per AF bucket. |

Supported flags (all accept POSIX short + GNU long forms):

| Short | Long                  | Meaning |
| ----- | --------------------- | ------- |
| `-s`  | `--samples`           | Restrict to a comma list of samples. |
| `-S`  | `--samples-file`      | Sample IDs, one per line. |
| `-r`  | `--regions`           | Region list (`chr:beg-end[,…]`) — post-filter (no index). |
| `-R`  | `--regions-file`      | BED-like regions file. |
| `-t`  | `--targets`           | Like `-r`, always a post-filter. |
| `-T`  | `--targets-file`      | BED-like targets file. |
| `-i`  | `--include`           | Keep records matching expression (same syntax as `view`). |
| `-e`  | `--exclude`           | Drop records matching expression. |
| `-f`  | `--apply-filters`     | Keep only PASS or named filters. |
| `-d`  | `--depth`             | `MIN,MAX,STEP` depth bins (default `0,500,1`). |
| `-a`  | `--af-bins`           | Bin edges (default `0,0.1,…,0.9,0.99,1.0`). |
| `-c`  | `--collapse`          | Accepted; v1 always treats each ALT separately. |
| `-1`  | `--1st-allele-only`   | Count only the first ALT allele. |
|       | `--af-tag TAG`        | Read AF from `INFO/TAG` instead of computing it from GT. |
| `-o`  | `--output`            | Output file (default stdout). |
|       | `--threads`           | Accepted; v1 is single-threaded. |

### Intentional deviations from upstream output

- The IDD section has columns `[length] [count] [nGenotypes] [meanVAF]`; the
  per-genotype and mean-VAF columns are placeholders (always `0` / `0.00`)
  because they require the upstream call cache.
- The AF section emits ten data rows (lower bin edges 0.0 … 0.9, plus the
  `0.99` bucket). The trailing `[8]repeat-consistent`, `[9]repeat-inconsistent`,
  and `[10]not applicable` columns are reported as zero — those come from
  upstream's local-realignment step, which we do not perform.
- Mixed sites (a SNP plus an indel at the same record) are counted in
  *both* the SNP and the indel SN counters once; this matches the upstream
  default for multi-allelic non-collapsed input.

## `bcftools query`

`bcftools query` writes records through a format string. Tokens:

| Token            | Meaning                                              |
| ---------------- | ---------------------------------------------------- |
| `%CHROM`         | Contig name.                                         |
| `%POS`           | 1-based position.                                    |
| `%REF`           | Reference allele.                                    |
| `%ALT`           | Comma-joined ALT alleles.                            |
| `%QUAL`          | Quality score (or `.` when missing).                 |
| `%ID`            | Variant ID (or `.`).                                 |
| `%FILTER`        | `;`-joined FILTER column.                            |
| `%TYPE`          | `SNP` / `MNP` / `INDEL` / `OTHER`.                   |
| `%INFO/<TAG>`    | INFO field by name. Missing keys render as `.`.      |
| `%INFO`          | The entire INFO column, in record order.             |
| `%GT`            | Raw GT (sample context only).                        |
| `%TGT`           | Translated genotype (`0/1` -> `A/T`).                |
| `%SAMPLE`        | Sample name (sample context only).                   |
| `%FMT/<TAG>`     | Sample FORMAT field by name (incl. `Character` tags). |
| `[%TOKEN ...]`   | Sample-repeated: emit inner once per sample, tab-joined. |
| `\n`, `\t`       | Literal newline / tab.                               |

### Flags

| Short | Long                | Meaning |
| ----- | ------------------- | ------- |
| `-f`  | `--format`          | Format string (required unless `-l`). |
| `-H`  | `--print-header`    | Emit a header row derived from the format string. |
| `-l`  | `--list-samples`    | Print one sample per line and exit. |
| `-s`  | `--samples`         | Comma list narrowing per-sample expansion. |
| `-S`  | `--samples-file`    | File of sample IDs (one per line). |
| `-r`  | `--regions`         | Region list (`chr:beg-end[,...]`). Uses `.tbi` / `.csi` when present. |
| `-R`  | `--regions-file`    | BED-like regions file. |
| `-t`  | `--targets`         | Like `-r` but always a post-filter (no index needed). |
| `-T`  | `--targets-file`    | BED-like targets file (post-filter). |
| `-i`  | `--include`         | Keep records matching expression. |
| `-e`  | `--exclude`         | Drop records matching expression. |
| `-F`  | `--apply-filters`   | Comma list of FILTER names to keep. (Upstream's `-f`; `-f` is the format string here.) |
| `-o`  | `--output`          | Output path (default stdout). |
|       | `--threads`         | Accepted; v1 is single-threaded. |
| `-?`  | `--help`            | Show help. |
|       | `--version`         | Show version. |

### Examples

```bash
# TSV of variant coordinates.
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' input.vcf

# Per-sample GT and DP, tab-separated within each row.
bcftools query -f '%CHROM\t%POS\t[%GT\t%DP]\n' input.vcf

# Filtered by depth, with a -H header row.
bcftools query -H -f '%CHROM\t%POS\t%INFO/DP\n' -i 'INFO/DP>30' input.vcf

# Sample names from a CSI-indexed BCF.
bcftools query -l input.bcf
```

## `bcftools concat`

`bcftools concat` concatenates two or more VCF/BCF inputs. By default it
assumes the inputs are sorted and non-overlapping; pass `-a` for a sort-merge
and `-D` to collapse duplicate (chrom, pos, ref, alt) records.

### Flags

| Short | Long                 | Meaning |
| ----- | -------------------- | ------- |
| `-a`  | `--allow-overlaps`   | Sort-merge across inputs. |
| `-D`  | `--remove-duplicates`| Drop adjacent duplicate records. |
| `-f`  | `--file-list`        | Read inputs from a file (one path per line). |
| `-O`  | `--output-type`      | `v` (default), `z`, `u`, `b`. |
| `-o`  | `--output`           | Output path (default stdout). |
| `-l`  | `--ligate`           | Ligate overlapping phased chunks (chunked imputation output): the overlap is emitted once, phase is reconciled across chunks, and `FORMAT/PS` and `FORMAT/PQ` are added. |
| `-q`  | `--min-PQ`           | Break the phase set when a sample's phasing quality at a boundary is below this value (default 30). |
|       | `--ligate-force`     | Ligate even non-overlapping chunks, keeping all sites. |
|       | `--ligate-warn`      | Drop sites in imperfect overlaps instead of erroring. |
|       | `--compression-level`| gzip level for `-O z` output. |
|       | `--threads`          | Accepted; v1 is single-threaded. |
| `-?`  | `--help`             | Show help. |
|       | `--version`          | Show version. |

### Header merging rules

- `##fileformat` is taken from the first input.
- `##contig` lines are union-merged in first-seen order.
- `##INFO`, `##FORMAT`, `##FILTER` lines are union-merged by ID. If the same
  ID appears with a different definition the command errors out with a
  message naming the conflicting line.
- Other meta lines are de-duplicated by exact-string equality.
- The sample sets of all inputs must match (same names in the same order);
  mismatches abort with an error.

### Examples

```bash
# Plain concat of two pre-sorted, non-overlapping VCFs.
bcftools concat a.vcf b.vcf -o joined.vcf

# Sort-merge across overlapping chunks (e.g. per-chromosome shards).
bcftools concat -a chr*.vcf.gz -O z -o all.vcf.gz

# Read the input list from a file.
bcftools concat -f files.txt -O b -o all.bcf

# Drop adjacent duplicates after sort-merge.
bcftools concat -a -D split.*.vcf.gz -o uniq.vcf
```

## `bcftools norm`

`bcftools norm` performs the standard pre-analysis fix-ups: left-align indels
against a reference FASTA, split multiallelics into biallelics (or vice
versa), drop duplicate records, and atomize complex variants.

```bash
# Left-align indels.
bcftools norm -f ref.fa in.vcf > out.vcf

# Split multi-allelic SNPs and indels into biallelic records, then left-align.
bcftools norm -f ref.fa -m -both in.vcf > out.vcf

# Skip left-alignment but split.
bcftools norm -N -m -any in.vcf > out.vcf

# Drop byte-for-byte duplicate records, write gzipped VCF.
bcftools norm -d exact -O z -o clean.vcf.gz in.vcf
```

### Supported `norm` flags

| Short | Long                | Meaning |
| ----- | ------------------- | ------- |
| `-f`  | `--fasta-ref`       | Reference FASTA. Builds an in-memory `.fai` on the fly if no sidecar exists. |
|       | `--check-ref`       | `e` (error, default), `w` (warn), `s` (skip) on REF / FASTA mismatch. |
| `-m`  | `--multiallelics`   | `-snps`, `-indels`, `-both`, `-any` to split; `+snps`, ... to join. |
| `-d`  | `--rm-dup`          | Drop duplicates: `snps`, `indels`, `both`, `all`, `exact`, `none` (default). |
| `-a`  | `--atomize`         | Decompose complex variants (same-length REF/ALT > 1bp) into single-base events. |
| `-N`  | `--do-not-normalize`| Skip left-alignment (useful when combined with `-m`). |
| `-s`  | `--strict-filter`   | Apply `--apply-filters` before splitting (default: after). |
| `-r`  | `--regions`         | Region(s) `chr:beg-end[,...]` (post-filter; no index required). |
| `-R`  | `--regions-file`    | BED-like regions file. |
| `-t`  | `--targets`         | Like `-r` but always a post-filter. |
| `-T`  | `--targets-file`    | BED-like targets file (post-filter). |
|       | `--apply-filters`   | Comma list of FILTER values to keep. |
| `-O`  | `--output-type`     | `v` / `z` / `u` / `b` (same semantics as `view`). |
| `-o`  | `--output`          | Output file (default stdout). |
| `-l`  | `--compression-level` | gzip level for `-O z`. |
|       | `--threads`         | Accepted; v1 is single-threaded. |
| `-h`  | `--help`            | Show help. |
|       | `--version`         | Show version. |

### Left-alignment algorithm

We implement the classical Tan-Abecasis-Durbin algorithm used by `bcftools` and
`vt normalize`:

```text
repeat:
  if all alleles end with the same base AND not all are length 1:
    trim the trailing base from every allele
  if any allele is now empty:
    prepend the upstream reference base from the FASTA
  stop when no change is made.
```

A final pass trims any shared leading bases beyond the single VCF anchor.
Left-alignment is bounded by the chromosome start; an indel that would have
to walk past position 1 is reported as an error.

### `.fai` index

`pkg/htsgo/fasta` ships a small `.fai` reader. The format is the same
five-column tab-separated layout `samtools faidx` produces:

```text
NAME    LENGTH    OFFSET    LINEBASES    LINEWIDTH
```

If no sidecar `<fasta>.fai` exists alongside the FASTA, `bcftools norm`
calls `fasta.BuildIndex` to scan and build the equivalent on the fly. This
needs the FASTA to have uniform line widths per contig (the same constraint
`samtools faidx` enforces).

### Multi-allelic split / join

Splitting (`-m -*`) explodes a record with N ALTs into N biallelic records.
For each child record:

- `INFO/AC` and `INFO/AF` are narrowed to the chosen allele.
- `FORMAT/GT` is remapped: the chosen ALT becomes "1"; other ALTs collapse
  to "0".

Joining (`-m +*`) walks adjacent biallelics sharing CHROM/POS/REF and
collapses them into one multiallelic. `AC`/`AF` become comma lists; `GT`
allele indices are renumbered to match each donor record's new position.

## `bcftools call`

`bcftools call` decides whether each input site is a variant given the
per-position genotype likelihoods produced upstream by `samtools mpileup
-g/--BCF`. The two operate as a pair: `mpileup` emits `FORMAT/PL` (Phred-
scaled likelihoods) per site per sample, and `call` aggregates those
likelihoods, applies a Hardy-Weinberg + mutation-rate prior, picks the
most-likely genotype per sample, and emits the standard variant record.

```bash
# Consensus caller (Li 2011), drop all-reference sites.
bcftools call -c -v input.vcf > out.vcf

# Multi-allelic caller writing BGZF-compressed VCF.
bcftools call -m -v -O z -o out.vcf.gz input.vcf

# Treat samples as haploid (e.g. for chrY).
bcftools call -c --ploidy 1 -v input.vcf

# Keep every declared ALT, including those with zero supporting reads.
bcftools call -c -A input.vcf

# Tighten the call rate.
bcftools call -c -v -p 0.01 -P 1e-4 input.vcf
```

### Supported `call` flags

| Short | Long                    | Meaning |
| ----- | ----------------------- | ------- |
| `-c`  | `--consensus-caller`    | Old Li-2011 consensus caller. |
| `-m`  | `--multiallelic-caller` | Full multi-allelic caller (the `-m` PL grid over >2 ALTs pairs with the mpileup MAQ port). |
| `-A`  | `--keep-alts`           | Emit every declared ALT, even those with zero supporting reads. |
| `-v`  | `--variants-only`       | Drop all-reference sites. |
| `-P`  | `--prior`               | Mutation rate prior (default `1.1e-3`). |
| `-p`  | `--pval-threshold`      | Variant-posterior threshold (default `0.5`). |
|       | `--ploidy`              | `2` (default), `1`, `GRCh37`, `GRCh38`, or `--ploidy-file`. All implemented. |
| `-X`  | `--chromosome-X`        | Legacy alias for `--ploidy 1`. |
| `-O`  | `--output-type`         | `v`, `z`, `u`, `b` (same semantics as `view`). |
| `-o`  | `--output`              | Output file (default stdout). |
| `-r`  | `--regions`             | Region(s) `chr:beg-end[,...]` (v1: post-filter only). |
| `-R`  | `--regions-file`        | BED-like regions file. |
| `-t`  | `--targets`             | Like `-r` but always a post-filter. |
| `-T`  | `--targets-file`        | BED-like targets file (post-filter). |
| `-s`  | `--samples`             | Restrict to these samples. |
| `-S`  | `--samples-file`        | File of sample IDs. |
|       | `--threads`             | Accepted; v1 is single-threaded. |
| `-?`  | `--help`                | Show help. |
|       | `--version`             | Show version. |

### Calling algorithm (consensus, `-c`)

For each input site:

1. For each sample, find the most-likely genotype from `FORMAT/PL` (the
   index with PL=0, or the smallest PL when there is no zero).
2. Aggregate ALT-allele support across samples (`INFO/AC` / `INFO/AN`).
3. Compute a per-site variant-posterior by summing the per-sample
   log-likelihood ratio "best-non-ref vs ref" (capped at +/-200 to keep
   the maths well-behaved on extreme PLs), then combining with the
   `-P` mutation-rate prior in a one-vs-rest logistic.
4. The site is "variant" iff `posterior > 1 - pvalThreshold` or any
   sample's best genotype is non-reference. Without `-v` every site is
   still emitted (with the called GTs); with `-v` non-variant sites are
   dropped (unless `-A` overrides).
5. Compute `QUAL` as `-10 * log10(1 - posterior)`, clamped to `[0, 999]`.
6. Rewrite `FORMAT/GT` per sample from the chosen genotype index.

### Calling algorithm (multi-allelic, `-m`)

The full multi-allelic caller iterates over every possible allele
combination and picks the most likely combined genotype across all
samples. This is now implemented: the `-m` PL grid over >2 ALTs pairs
with the mpileup MAQ port (the consensus `-c` path also remains
available). See `docs/PARITY_ROADMAP.md` under the `bcftools call`
entry for the validation detail.

### `convert` PLINK export notes

Upstream's `vcfconvert.c` leaves its `--plink`/`--tped`/`--bin` option
block commented out (no `case 'p'`, no implementation), so there is no
upstream binary to diff against. These exporters are implemented directly
to the [PLINK1 file-format
spec](https://www.cog-genomics.org/plink/1.9/formats):

- `-p`/`--plink <prefix>` writes a PLINK1 text fileset: `<prefix>.ped`
  (one line per sample: 6 mandatory columns `FID IID PAT MAT SEX PHENO`
  with `FID=IID=sample`, `PAT=MAT=0`, `SEX=0`, `PHENO=-9`, then the two
  REF/ALT allele letters per variant, `0 0` for a missing genotype) and
  `<prefix>.map` (`CHROM SNP_ID 0 BP`).
- `--tped <prefix>` writes the transposed text fileset: `<prefix>.tped`
  (one line per variant: `CHROM SNP_ID 0 BP` then the two alleles of every
  sample) and `<prefix>.tfam` (the 6 `.ped` columns, one line per sample).
- `--plink-bed <prefix>` writes the PLINK1 **binary** fileset:
  `<prefix>.bed` (3-byte magic `0x6c 0x1b 0x01`, SNP-major; per variant
  `ceil(nsamples/4)` bytes, 2 bits per sample little-endian within the
  byte), `<prefix>.bim` (`CHROM SNP_ID 0 BP A1 A2`) and `<prefix>.fam`.

Conventions (chosen because upstream has none to match):

- **`A1=ALT`, `A2=REF`** in the `.bim`, matching `plink --make-bed
  --keep-allele-order` from a VCF. The on-disk `.bed` 2-bit code is the
  ALT-allele count: `00`=hom-REF (two A2), `10`=het, `11`=hom-ALT (two A1),
  `01`=missing.
- **Chromosome codes**: numeric `1..22` pass through; `X→23`, `Y→24`,
  `XY→25`, `MT`/`M→26`; an optional `chr`/`CHR` prefix is stripped first.
  Any other contig name is written verbatim (PLINK1.9 `--allow-extra-chr`).
- **Biallelic only**: PLINK is a biallelic format, so records with more
  than one ALT (or no ALT) are skipped and counted, with a one-line stderr
  warning on the first multi-allelic site. Split with `bcftools norm -m-`
  first to retain them. The `SNP_ID` is the VCF `ID`, or `CHROM:POS` when
  the ID is `.`.

The exporters accept a bare `<prefix>` (canonical suffixes) or an explicit
comma-separated file list (2 names for `--plink`/`--tped`, 3 for
`--plink-bed`).

### Deviations from upstream

- **BCF input edge cases.** Some htslib-produced BCF `FORMAT`-key
  reconstruction edge cases remain (see `docs/UPSTREAM_BUGS.md`
  `bcf-fmt-keys-missing`); for those, convert with
  `bcftools view in.bcf > in.vcf` first. Tracked in the roadmap.
- **No index-backed region queries** for `-r`; the flag is honoured as
  a post-filter. The view-side CSI seek will move here in a follow-up.

## Scope and deferred work

All 24 bcftools subcommands now ship: `view`, `index`, `stats`, `query`,
`concat`, `norm`, `call`, `merge`, `annotate`, `isec`, `sort`, `reheader`,
`convert`, `mendelian`, `mendelian2`, `gtcheck`, `roh`, `filter`,
`consensus`, `mpileup`, `cnv`, `polysomy`, `csq`, and `plugin`. BCF
read + write (`-O u`/`-O b`), VCF and BGZF VCF I/O (`-O v`/`-O z`),
CSI/TBI indexing, and multi-threaded (`-@`) output compression are all
in. The full multi-allelic `-m` caller, `--ploidy GRCh37`/`GRCh38`,
mpileup SNP genotype likelihoods (slices 1–4) plus both the legacy
`bam2bcf_indel` and the `--indels-cns` (edlib) indel callers, the csq
haplotype engine (slices 1–4), and the `roh`/`cnv`/`polysomy` HMMs are
implemented and live-oracle validated.

The `convert` PLINK exporters (`--plink`/`--tped`/`--plink-bed`) are now
implemented (see the PLINK export notes above); the GEN/HAP/TSV/gVCF modes
were already done.

`gtcheck` is now feature-complete against upstream: GT and PL scoring,
the probability/integer discordance paths, the HWE column, `--n-matches`,
`--distinctive-sites`, `-i/-e` filter expressions (with the `qry:`/`gt:`
scope prefix), the `-c/--cluster` clustering (this port's own design;
upstream leaves it an error stub), and both `-O` output containers
(`t` text and `z` BGZF, with an optional `z<N>` compression level).

What remains **genuinely open** (see `docs/PARITY_ROADMAP.md` for the
authoritative gap list):

- `query %N_ALT` (**not** an upstream `query` token — non-goal).
- `bcftools tview` is a deliberate skip (interactive ncurses; no pipeline
  use).

(`bcftools som` is now **implemented** — the upstream `fwrite`-return write
bug is fixed so train→classify works; see `docs/UPSTREAM_BUGS.md` and the
`som` entry above.)
- The ~30 bundled upstream `.so` plugins are non-goal scope (the plugin
  *system* — a VCF-on-stdin/stdout subprocess protocol — is implemented).

## How records flow

```text
on disk (.vcf, .vcf.gz, .bcf)
        │
        ▼
pkg/htsgo/iohelper.OpenReader    (auto-detects BGZF / gzip / plain)
        │
        ▼
streaming dispatcher in pkg bcftools
   ├── BCF? → pkg/htsgo/bcf decoder → ToVariant()
   └── VCF? → pkg/htsgo/vcf decoder
        │
        ▼
filter pipeline (region / -t / -f / -c / -C / -q / -Q / -i / -e / -s / -G)
        │
        ▼
vcf.Writer  →  stdout / -o file / BGZF wrap when -O z / BCF when -O b|u
```

## Tests

Run from the repo root:

```bash
go test ./pkg/htsgo/bcf/... ./tools/bcftools/...
go test -race -cover ./pkg/htsgo/bcf/... ./tools/bcftools/...
```

Coverage targets:

- `pkg/htsgo/bcf` ≥ 80% (BCF parsing is fiddly; we hit ~82% in this
  slice).
- `tools/bcftools/pkg/bcftools` ≥ 85% (we hit ~86% with the `stats` slice
  in place; `stats.go`-only coverage is ~93%).
