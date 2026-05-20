# tabix — pure-Go generic random-access index

`tabix` is a pure-Go reimplementation of htslib's `tabix` command and the
`.tbi` index format. It builds a position-index (`<file>.tbi`) for a
bgzipped tab-delimited file (VCF, BED, GFF, SAM, or any custom layout) and
then supports fast region queries against the indexed file:

```bash
tabix -p vcf calls.vcf.gz                  # build calls.vcf.gz.tbi
tabix calls.vcf.gz chr1:1,000,000-2,000,000   # print matching records
```

The implementation builds on `pkg/htsgo/bgzf` (the in-tree BGZF
codec) and has **no third-party dependencies** — pure Go standard library.

It is pick #2 of the 2026 next-up list in `analysis/tool_ranking_2026.md`
because tabix's bin/linear-index machinery is what the rest of the
htslib-style ecosystem (`samtools view`, `bcftools view`, BAI indices)
re-uses.

## What is the `.tbi` format?

A `.tbi` file is a bgzipped binary index with this layout:

1. 4-byte magic `TBI\1`
2. `n_ref`, `format`, `col_seq`, `col_beg`, `col_end`, `meta`, `skip` as
   little-endian `int32`s
3. The length-prefixed concatenation of NUL-terminated chromosome names
4. For each reference: an array of `(bin, chunks…)` records followed by
   a linear-index array of `uint64` virtual offsets, one per 16-kbp tile
5. Optional `n_no_coor` `uint64` trailer for records without coordinates

Each chunk is a half-open `[Beg, End)` pair of 64-bit *virtual offsets*:
the top 48 bits index into the compressed file (a BGZF block start), the
bottom 16 bits give the byte offset inside that block's decompressed
payload. With BGZF blocks capped at 64 KiB the low 16 bits are enough.

Records are placed in the UCSC binning scheme — five levels of fixed-size
bins from `2^14` bp at the leaves up to `2^29` bp at the root. A region
query consults every bin that overlaps the requested range, then trims
chunks against the linear-index lower bound `Linear[beg >> 14]` and merges
adjacent chunks before reading the underlying BGZF.

## Build

```bash
go build ./tools/tabix/cmd/tabix
```

## Usage

```text
tabix [options] file.gz                      # build index → file.gz.tbi
tabix [options] file.gz REGION [REGION...]   # query records by region
```

A `REGION` is `CHROM`, `CHROM:START-END`, or `CHROM:START` (all 1-based,
inclusive). The CLI translates 1-based positions to the 0-based half-open
form used internally.

### Flags

| Short | Long              | Description                                      |
|-------|-------------------|--------------------------------------------------|
| `-p`  | `--preset`        | Standard preset: `gff`, `bed`, `sam`, `vcf`.     |
| `-s`  | `--seq-col N`     | 1-based column of the sequence (chrom) name.     |
| `-b`  | `--begin-col N`   | 1-based column of the begin position.            |
| `-e`  | `--end-col N`     | 1-based column of the end (0 = use begin only).  |
| `-S`  | `--skip-lines N`  | Number of header lines to skip past (default 0). |
| `-c`  | `--meta-char C`   | Comment-line prefix character (default `#`).     |
| `-0`  | `--zero-based`    | 0-based half-open coordinates (BED-style).       |
| `-f`  | `--force`         | Overwrite an existing `.tbi` index.              |
| `-R`  | `--regions FILE`  | Read regions from a BED-like file.               |
| `-T`  | `--targets FILE`  | Restrict output to records overlapping FILE.     |
| `-l`  | `--list-chroms`   | Print chromosome names recorded in the index.    |
| `-h`  | `--print-header`  | Also emit header lines from the queried file.    |
|       | `--only-header`   | Emit only the header lines.                      |
| `-D`  |                   | Do not save the index (for `build` mode only).   |
|       | `--help`          | Show help and exit.                              |
| `-v`  | `--version`       | Show version and exit.                           |

### Presets

| Preset | Format          | seq | beg | end | meta | skip |
|--------|-----------------|-----|-----|-----|------|------|
| `gff`  | generic         | 1   | 4   | 5   | `#`  | 0    |
| `bed`  | generic, 0-based| 1   | 2   | 3   | `#`  | 0    |
| `sam`  | SAM             | 3   | 4   | 0   | `@`  | 0    |
| `vcf`  | VCF             | 1   | 2   | 0   | `#`  | 0    |

### Examples

Index and query a VCF:

```bash
bgzip my.vcf
tabix -p vcf my.vcf.gz                    # writes my.vcf.gz.tbi
tabix my.vcf.gz chr1:1000000-2000000      # records in that region
tabix -h my.vcf.gz chr1:1000000-2000000   # include header
tabix --only-header my.vcf.gz             # header only
tabix -l my.vcf.gz                        # list chromosomes
```

Index a BED file:

```bash
bgzip features.bed
tabix -p bed features.bed.gz
tabix features.bed.gz chr1:1-1000000
```

Custom column layout (chrom in col 5, single position in col 6):

```bash
tabix -s 5 -b 6 -c '#' weird.tsv.gz
```

## Library

The CLI is a thin wrapper around `pkg/htsgo/tabix`:

```go
import "github.com/yassineS/bio_ai_experiment/pkg/htsgo/tabix"

cfg, _ := tabix.PresetConfig(tabix.PresetVCF)
idx, _ := tabix.Build("calls.vcf.gz", cfg)
idx.WriteFile("calls.vcf.gz.tbi")

reloaded, _ := tabix.ReadFile("calls.vcf.gz.tbi")
records, _ := reloaded.QueryBytes("calls.vcf.gz", "chr1", 999_999, 2_000_000)
for _, rec := range records {
    fmt.Printf("%s\n", rec)
}
```

Exported types:

- `Index`, `Config` — the parsed `.tbi` and its column/format spec.
- `Bin`, `Chunk`, `RefIndex` — the bin / linear-index data structures.
- `VOffset`, `MakeVOffset` — 48-bit-coff + 16-bit-uoff packed offsets.
- `Reg2bin(beg, end) int`, `Reg2bins(beg, end) []int`,
  `LinearTile(pos) int` — the UCSC binning math, exposed because BAI
  re-uses the same scheme.

## Deviations from upstream

This v1 implementation matches upstream `tabix` for the documented build /
query / list / header flags with the following intentional differences:

1. **`-r` / `--reheader` is not implemented.** htslib's `tabix --reheader`
   replaces the header lines of an existing bgzipped file in place; that
   is a write-path operation orthogonal to indexing and queries. It will
   land in a follow-up PR.
2. **`-T` / `--targets` is parsed but currently behaves as `-R`.** A real
   targets filter (post-filter records to only those overlapping a set of
   target regions) is a thin layer on top of the existing query path and
   will be tightened up alongside the same PR as `--reheader`.
3. **Linear-index "no record yet" sentinel.** Internally the builder uses
   `^uint64(0)` while accumulating per-tile minimum virtual offsets,
   replacing the sentinel with `0` (or the carry-forward of the last
   recorded offset, per the htslib convention) in `finalize()` before
   serialisation. This is invisible on disk — the byte layout matches
   htslib byte-for-byte.

## Layout

```text
tools/tabix/
├── cmd/tabix/main.go              # CLI entry point
├── pkg/tabix/
│   ├── doc.go                     # package overview
│   ├── binning.go                 # Reg2bin, Reg2bins, LinearTile
│   ├── voffset.go                 # VOffset packing
│   ├── tabix.go                   # Index, Build / Read / Write / QueryBytes
│   ├── binning_test.go
│   ├── voffset_test.go
│   ├── tabix_test.go
│   └── tabix_extra_test.go
└── README.md
```

## What this unblocks

With tabix and its bin / linear-index machinery in hand, the natural
follow-ups (NOT part of this PR) are:

- `samtools view / sort / index` — BAI files re-use exactly the same
  binning scheme; the only structural difference is that BAI computes
  end-of-record from the CIGAR string rather than from a column index.
- `bcftools view` random-access queries — same `.tbi` layer.
- Wiring tabix-backed lookup into `pkg/htsgo/iohelper` so the rest
  of the repo gets transparent VCF region access.

## References

- Heng Li, *Tabix: fast retrieval of sequence features from generic
  TAB-delimited files*, Bioinformatics 27 (5), 2011.
- The SAM/BAM/Tabix specification, section "The tabix index file format".
- htslib sources [`htslib/hts.c`](https://github.com/samtools/htslib),
  `tabix.c`.
