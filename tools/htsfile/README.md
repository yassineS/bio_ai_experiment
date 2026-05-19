# htsfile — bioinformatics file format identifier

A pure-Go re-implementation of htslib's `htsfile` binary. Sniffs each
input file and prints a one-line summary of its detected format,
container, and (where known) version.

## What it identifies

| Container         | Payloads                                  |
|-------------------|-------------------------------------------|
| Plain text        | SAM, VCF, FASTA, FASTQ, BED, GFF, generic |
| Plain gzip        | All of the above                          |
| BGZF              | All of the above + BAM, BCF               |
| Binary (no gzip)  | BAM, BCF, CRAM                            |

For SAM / VCF / GFF / CRAM the version is extracted from the file's
own version field (`@HD VN:`, `##fileformat=VCFv`, `##gff-version`,
the CRAM `major.minor` magic-byte tail). For BCF the version comes
from the magic bytes (`BCF\2\2` → 2.2, `BCF\2\1` → 2.1). BAM has no
in-band version field.

## Usage

```
htsfile FILE [FILE...]
```

Each FILE is identified independently. Pass `-` to sniff stdin.

```
$ htsfile alignments.bam variants.vcf.gz reads.fastq.gz
alignments.bam: BAM BGZF-compressed sequence data
variants.vcf.gz: VCF version 4.2 BGZF-compressed variant calling data
reads.fastq.gz: FASTQ gzip-compressed sequence data
```

`-h` / `--help` and `-v` / `--version` are also accepted.

## Differences from upstream

- We don't link against libhts — the sniff is a pure-Go peek that
  never decompresses more than the first BGZF block.
- The `-c` "concatenate to stdout" mode of upstream `htsfile` is
  **not implemented** in v1. The scope is identification only;
  callers wanting the raw bytes can use shell redirection.
- The detection result for ambiguous inputs is conservative:
  rather than guess between two formats with overlapping prefixes,
  we report the most-specific match (e.g. `@read1\nACGT\n+\n!!!!\n`
  is classified as FASTQ, not SAM, even though both start with `@`).
- The "mostly text" fallback uses an ASCII-only heuristic: bytes
  in the printable range (0x20–0x7E) plus `\t`, `\n`, `\r` count
  as text; NUL forces a "binary" classification. A UTF-8 file
  with bytes ≥ 0x80 (e.g. a comment containing `é`) will fall
  into PayloadBinary — same limitation as upstream htsfile.

## Implementation

The sniff is split into:

- `tools/htsfile/pkg/htsfile/sniff.go` — `Identify(path)` and
  `IdentifyReader(io.Reader)` produce a `*Format` describing
  compression + payload + version. The `Format.Describe()`
  method emits the one-line summary string.
- `tools/htsfile/cmd/htsfile/main.go` — the CLI that iterates
  arguments and prints `path: <description>` per file.

`pkg/htsfile` is intentionally kept inside the tool tree (not
under `pkg/htsgo/`) because the heuristics are CLI-shaped — they
look at a fixed-size prefix and produce a text summary, not a
streaming reader. The `pkg/htsgo` packages stay format-primitive.
