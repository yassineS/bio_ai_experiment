# bedgetfasta

Pure-Go reimplementation of `bedtools getfasta`. For each BED interval, pulls
the corresponding subsequence from a FAI-indexed FASTA reference. Output is
FASTA by default; TSV with `-tab`.

## Usage

```bash
# Default: emit FASTA, header = chrom:start-end
bedgetfasta -fi genome.fa -bed peaks.bed > peaks.fa

# Use BED name column as header
bedgetfasta -fi genome.fa -bed peaks.bed -name > peaks.fa

# Reverse-complement '-' strand intervals (case + IUPAC preserved)
bedgetfasta -fi genome.fa -bed peaks.bed -s -name > peaks.fa

# Split BED12 records and concatenate the per-block sequence
bedgetfasta -fi genome.fa -bed transcripts.bed12 -split -s > transcripts.fa
```

## Flags

| Short | Long | Notes |
|---|---|---|
| `-fi` | `--fasta` | FASTA reference. Required. `.fai` built on demand. |
| `-bed` | `--bed` | BED file. Required. `-` = stdin. Transparent gzip. |
| `-fo` | `--output` | Output file (default stdout). `-` = stdout. |
| `-name` | `--name` | Header is `<name>::<chrom>:<start>-<end>`. |
| `-name+` | | Deprecated alias of `-name` (identical header). |
| `-nameOnly` | `--nameOnly` | Header is just `<name>`. |
| `-tab` | `--tab` | TSV output (`<header>\t<seq>`). |
| `-bedOut` | `--bedOut` | Re-emit the BED record with a trailing sequence column (tab-delimited) instead of FASTA. |
| `-s` | `--strand` | Reverse-complement `-` strand intervals (IUPAC + case preserved). |
| `-split` | `--split` | Split BED12 records into their constituent blocks. |
| `-rna` | `--rna` | Emit `U/u` in place of `T/t` (applied after `-s`). |
| `-fullHeader` | `--full-header` | Index FASTA contigs by the full header line (whitespace included). |
| `-h` | `--help` | |
| `-v` | `--version` | |

## Behaviour

- **Case-preserving fetch.** Where the shared `pkg/htsgo/fasta` random
  access uppercases output for case-insensitive downstream comparison
  (`bcftools norm` etc.), `bedgetfasta` ships its own case-preserving
  `FetchPreserveCase` so IUPAC codes round-trip exactly.
- **Strand suffix.** `-s` toggles the `(+)`/`(-)` suffix in the header (the
  default behaviour without `-s` omits it).
- **Name headers with no name column.** Matching upstream, `-name` on a BED
  row without a name column emits `>::chrom:start-end` (empty name, *not* a
  fall-back to `chrom:start-end`); `-nameOnly` emits an empty header (`>`).
- **`-bedOut`.** Re-emits the original BED columns followed by a trailing
  sequence column. The sequence still honours `-s`, `-split` and `-rna`.
  Columns beyond 6 are preserved verbatim (matches upstream's
  `reportBedTab`).
- **Missing chromosome.** A BED interval whose chromosome is not in the FAI
  produces a warning on stderr (`WARNING. chromosome (X) was not found in
  the FASTA file. Skipping.`) and the record is dropped, matching upstream.
- **Out-of-range coordinates.** A feature extending past the contig length is
  skipped with `Feature (chrom:start-end) beyond the length of chrom size
  (N bp).  Skipping.` (note the two spaces before `Skipping`, matching
  upstream) rather than aborting.
- **Zero-length features.** A record with `start == end` is skipped with
  `Feature (chrom:start-end) has length = 0, Skipping.` and produces no
  output, matching upstream.
- **Stale index warning.** When a sibling `.fai` exists but is older than the
  FASTA file, `Warning: the index file is older than the FASTA file.` is
  emitted on stderr (matches upstream/htslib `getfasta.t10`).
- **`-split` + `-s`.** Blocks are extracted in genomic order, concatenated,
  and the whole sequence is then reverse-complemented — exactly as upstream
  `ReportSeq` does. Matches upstream `getfasta.t05`.
- **`-rna`.** Applied after any reverse-complement: T→U and t→u (case
  preserved). Other bases pass through.

## BGZF input

`-fi *.fa.gz` (bgzipped FASTA) is supported transparently: the file is
sniffed for the BGZF magic and routed through `pkg/htsgo/fasta`'s
`OpenRandomAccessBGZF`, which fully decompresses the payload in-memory
and reuses the standard FAI index path. A samtools-style sibling
`<path>.fa.gz.fai` is honoured when present; otherwise the index is
rebuilt from the decompressed payload. A `.gzi` sidecar (when present)
is parsed for validation only — partial-decompression seek via `.gzi`
is a future optimisation, not required for parity. Upstream
`getfasta.t18` passes byte-for-byte.

## Not yet implemented

(No outstanding gaps in the upstream `getfasta` test corpus.)

See [`../PARITY_VALIDATION.md`](../PARITY_VALIDATION.md) for the validated
parity matrix.
