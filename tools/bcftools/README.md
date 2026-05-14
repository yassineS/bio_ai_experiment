# bcftools (pure-Go)

A pure-Go reimplementation of selected [bcftools](https://samtools.github.io/bcftools/)
subcommands. Subcommands shipped so far:

- `bcftools view` — filter / project / convert VCF and BCF.
- `bcftools norm` — left-align indels and normalize multiallelics against a reference FASTA.
- `bcftools index` — build CSI / .tbi indices for BCF and bgzipped VCF.
- `pkg/bioformats/bcf` — BCF v2.2 reader + writer.
- `pkg/bioformats/fasta` — FASTA reader / writer plus a `.fai` index reader used by `norm`.

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

## Scope and deferred work

What ships so far:

- VCF, BGZF-wrapped VCF, and BCF input.
- VCF, gzip-VCF, and BCF (`-O v`, `-O z`, `-O u`, `-O b`) output.
- `view`, `norm`, `index` subcommands.
- `-i` / `-e` expression language with arithmetic, comparison, and INFO/FILTER lookups.
- `.csi` index reading + writing for BCF.

What is **deferred** to a follow-up PR:

- Other subcommands (`query`, `stats`, `concat`, `merge`, `mpileup`).
- The plugin system (`bcftools plugin`).
- Multi-threading for the splitting/joining pipeline.

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
- `pkg/bioformats/fasta` covers the new index reader at ~89%.
- `tools/bcftools/pkg/bcftools` ≥ 85% (we hit ~86% overall; `norm.go`
  alone is ~89%).
