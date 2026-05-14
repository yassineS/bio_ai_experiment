# bcftools (pure-Go)

A pure-Go reimplementation of selected [bcftools](https://samtools.github.io/bcftools/)
subcommands. Currently shipping:

- `bcftools view` — filter / project / convert VCF and BCF.
- `bcftools query` — format-string output for VCF/BCF records.
- `bcftools concat` — concatenate VCF/BCF files.
- `bcftools index` — build a CSI / .tbi index for a BCF / VCF.gz file.
- `pkg/bioformats/bcf` — BCF v2.2 reader and writer.

Both pieces share the existing `pkg/bioformats/vcf` types so downstream
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

## Scope and deferred work

What ships today:

- BCF v2.2 reader and writer.
- VCF, BGZF-wrapped VCF and BCF input.
- VCF, gzip-VCF, BCF (`-O v`, `-O z`, `-O u`, `-O b`) output.
- `view`, `query`, `concat`, `index` subcommands with the flags listed in
  this README.
- CSI (BCF) and .tbi (VCF.gz) index reading and writing for region queries.

What is **deferred** to follow-up PRs:

- Other subcommands (`stats`, `norm`, `merge`, ...).
- The plugin system (`bcftools plugin`).
- `concat --ligate` (specialised imputation merging) and `--min-PQ` are
  accepted on the CLI but currently no-ops.

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
- `tools/bcftools/pkg/bcftools` ≥ 85% (we hit ~85%).
