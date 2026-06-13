# bcftools (pure-Go)

A pure-Go reimplementation of selected [bcftools](https://samtools.github.io/bcftools/)
subcommands. The current implementation ships:

- `bcftools view` — filter / project / convert VCF and BCF.
- `bcftools query` — format-string output for VCF/BCF records.
- `bcftools concat` — concatenate VCF/BCF files.
- `bcftools norm` — left-align indels, split/join multi-allelics, atomize, dedup.
- `bcftools call` — variant calling from per-position genotype likelihoods
  (consensus `-c` + full multi-allelic `-m`).
- `bcftools index` — build a `.csi` (or `.tbi`) index for a BCF / VCF.gz.
- `bcftools stats` — sectioned summary numbers compatible with
  `plot-vcfstats`.
- `bcftools convert` — re-emit VCF/BCF in a different format (`-O v|z|b|u`)
  with optional sample / region filtering, plus the GEN/HAP/TSV conversion
  modes: `--gvcf2vcf`, `--tsv2vcf`, `--gensample`/`--gensample2vcf`,
  `--hapsample`/`--hapsample2vcf`, `--haplegendsample`/`--haplegendsample2vcf`.
  Only the PLINK exporters remain unported (tracked in
  `docs/PARITY_ROADMAP.md`).
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
  sites). Multiple letters can be combined.
- `bcftools gtcheck` — sample-identity check by hard-GT Hamming
  concordance (with `-g panel`). Emits the upstream `tsv`-format
  `DC` / `INFO` tables.
- `bcftools roh` — runs-of-autozygosity detection via a 2-state HMM
  (`HW` / `AZ`) over hard-GT input (`-G/--GTs-only`).
- `bcftools filter` — soft-filter records by include / exclude
  expression. Failing records keep their place in the output but
  have FILTER set to the `-s/--soft-filter NAME` (or appended via
  `-m +`); optional `-S/--set-GTs .|0` rewrites failing samples'
  GTs. Supports `-g/--SnpGap` and `-G/--IndelGap` clustering.
- `bcftools consensus` — apply VCF variants (SNPs + simple indels)
  to a reference FASTA. Supports `-s/--samples`, `-H/--haplotype`
  (R/A/I/LR/LA/SR/SA + numeric), `-I/--iupac-codes`, `-m/--mask`
  with `--mask-with`, `--mark-ins / --mark-snv / --mark-del`,
  `-p/--prefix`, `-a/--absent`, `-M/--missing`.
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
  annotation. v1 ships only the protein-coding SNP classifier
  (missense / synonymous / stop\_gained / stop\_lost / start\_lost);
  indels, splice-site, compound-het, and haplotype-aware phasing
  are tracked in `docs/PARITY_ROADMAP.md`. Output is a VCF with an
  `INFO/BCSQ` tag of the form
  `consequence|gene|transcript|biotype|strand|aa_change|dna_change`.
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

What remains **genuinely open** (see `docs/PARITY_ROADMAP.md` for the
authoritative gap list):

- `convert` PLINK exporters (the GEN/HAP/TSV/gVCF modes are done).
- `gtcheck -c/--cluster` dendrogram (upstream itself errors "to be
  implemented") and `gtcheck` filter expressions.
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
