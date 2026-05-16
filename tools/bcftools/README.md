# bcftools (pure-Go)

A pure-Go reimplementation of selected [bcftools](https://samtools.github.io/bcftools/)
subcommands. The current implementation ships:

- `bcftools view` — filter / project / convert VCF and BCF.
- `bcftools query` — format-string output for VCF/BCF records.
- `bcftools concat` — concatenate VCF/BCF files.
- `bcftools norm` — left-align indels, split/join multi-allelics, atomize, dedup.
- `bcftools call` — variant calling from per-position genotype likelihoods
  (consensus + biallelic multi-allelic).
- `bcftools index` — build a `.csi` (or `.tbi`) index for a BCF / VCF.gz.
- `bcftools stats` — sectioned summary numbers compatible with
  `plot-vcfstats`.
- `bcftools convert` — re-emit VCF/BCF in a different format (`-O v|z|b|u`)
  with optional sample / region filtering. Upstream's exotic shapes
  (gVCF↔VCF, HAPLEGEND/HAPSAMPLE, PLINK) are tracked in
  `docs/PARITY_ROADMAP.md` and deferred.
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
- `pkg/bioformats/bcf` — reader and writer for the BCF v2.2 binary format.

All pieces share the existing `pkg/bioformats/vcf` types so downstream
consumers see records as familiar `vcf.Variant` values.

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

## Supported flags

| Short | Long                  | Meaning |
| ----- | --------------------- | ------- |
| `-O`  | `--output-type`       | `v` (default), `z`, `u`, `b`. `u`/`b` (BCF output) is **NOT YET implemented**; see scope below. |
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
| `%GT`            | Raw GT (sample context only).                        |
| `%TGT`           | Translated genotype (`0/1` -> `A/T`).                |
| `%FMT/<TAG>`     | Sample FORMAT field by name.                         |
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
| `-q`  | `--min-PQ`           | Accepted; no-op in v1. |
| `-l`  | `--ligate`           | Accepted; no-op in v1 (chunked imputation). |
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

`pkg/bioformats/fasta` ships a small `.fai` reader. The format is the same
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

# Multi-allelic caller (v1: biallelic-only) writing gzipped VCF.
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
| `-m`  | `--multiallelic-caller` | Multi-allelic caller (v1: biallelic-only; falls back to consensus on multi-allelic sites). |
| `-A`  | `--keep-alts`           | Emit every declared ALT, even those with zero supporting reads. |
| `-v`  | `--variants-only`       | Drop all-reference sites. |
| `-P`  | `--prior`               | Mutation rate prior (default `1.1e-3`). |
| `-p`  | `--pval-threshold`      | Variant-posterior threshold (default `0.5`). |
|       | `--ploidy`              | `2` (default), `1`, or `GRCh37` / `GRCh38` (deferred). |
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

Upstream's full multi-allelic caller iterates over every possible allele
combination and picks the most likely combined genotype across all
samples. **Our v1 implementation handles biallelic sites with the same
maths as `-c` and falls back to `-c` for multi-allelic sites.** Sites
with >2 alleles still produce a valid call (the most likely per-sample
genotype is selected from the full PL vector) but the posterior is
computed as if the site were biallelic with the most-supported ALT.
Tracked in `docs/PARITY_ROADMAP.md` under the `bcftools call` entry.

### Deviations from upstream

- **No BCF input.** The current decoder doesn't reconstruct
  `FORMAT/PL` from htslib-produced BCF (see `docs/UPSTREAM_BUGS.md`
  `bcf-fmt-keys-missing`). Convert with `bcftools view in.bcf > in.vcf`
  first. Tracked in the roadmap.
- **`--ploidy GRCh37` / `GRCh38` are rejected** at runtime; v1 only
  supports fixed ploidies. The CLI parses the spec so existing scripts
  fail early with a clear error.
- **No index-backed region queries** for `-r`; the flag is honoured as
  a post-filter. The view-side CSI seek will move here in a follow-up.
- **Multi-allelic caller** is biallelic-only, as described above.

## Scope and deferred work

What ships in this slice:

- `view`, `index`, `stats`, `query`, `concat`, `norm`, `call`,
  `merge`, `annotate`, `isec`, `sort`, `reheader`, `convert`,
  `mendelian`, `gtcheck`, `roh` subcommands.
- BCF reader + writer (CHROM/POS/REF/ALT/QUAL/FILTER/INFO + per-sample FORMAT).
- VCF and BGZF-wrapped VCF input.
- VCF and gzip-VCF (`-O v`, `-O z`) output for most subcommands.
- CSI / TBI index reading and writing.

What is **deferred** to a follow-up PR:

- Other subcommands (`csq`, `filter`, `mpileup`, `consensus`,
  `polysomy`, `cnv`, `mendelian2`, `+plugins`, ...).
- The plugin system (`bcftools plugin`).
- `bcftools stats -E exons.tab.gz` (upstream's exon-overlap section).
- BCF input for `bcftools call`, the full multi-allelic caller, and
  GRCh37 / GRCh38 ploidy specs.

## How records flow

```text
on disk (.vcf, .vcf.gz, .bcf)
        │
        ▼
pkg/bioformats/iohelper.OpenReader    (auto-detects BGZF / gzip / plain)
        │
        ▼
streaming dispatcher in pkg bcftools
   ├── BCF? → pkg/bioformats/bcf decoder → ToVariant()
   └── VCF? → pkg/bioformats/vcf decoder
        │
        ▼
filter pipeline (region / -t / -f / -c / -C / -q / -Q / -i / -e / -s / -G)
        │
        ▼
vcf.Writer  →  stdout / -o file / gzip wrap when -O z
```

## Tests

Run from the repo root:

```bash
go test ./pkg/bioformats/bcf/... ./tools/bcftools/...
go test -race -cover ./pkg/bioformats/bcf/... ./tools/bcftools/...
```

Coverage targets:

- `pkg/bioformats/bcf` ≥ 80% (BCF parsing is fiddly; we hit ~82% in this
  slice).
- `tools/bcftools/pkg/bcftools` ≥ 85% (we hit ~86% with the `stats` slice
  in place; `stats.go`-only coverage is ~93%).
