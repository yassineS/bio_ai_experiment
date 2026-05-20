# htsgo

A pure-Go bioinformatics format / index / compression library, scoped to
the surface htslib covers in C: BGZF, BAM, BAI, CRAM, CRAI, SAM, VCF,
BCF, Tabix/CSI, FASTA, FASTA index, FASTQ, BED, GFF, region parsing.

**Status:** **in-flight migration.** The legacy home was
`pkg/bioformats/`; htsgo is the consolidation target. See
`docs/HTSGO_ROADMAP.md` for the full plan and PR sequence (A → I).

## Current contents (after PR-B)

| Sub-package         | Source of truth                | Notes                                  |
|---------------------|--------------------------------|----------------------------------------|
| `iohelper`          | `pkg/htsgo/iohelper/`          | Transparent gzip + stdin/stdout + BGZF sniffing |
| `fasta`             | `pkg/htsgo/fasta/`             | Reader/Writer + faidx index + random access     |
| `fastq`             | `pkg/htsgo/fastq/`             | Reader/Writer with Phred33/64 enums    |
| `bed`               | `pkg/htsgo/bed/`               | BED + BEDPE + interval-tree                     |
| `gff`               | `pkg/htsgo/gff/`               | GFF/GTF reader                                  |
| `sam`               | `pkg/htsgo/sam/`               | SAM + BAM reader/writer, flag + cigar consts    |
| `vcf`               | `pkg/htsgo/vcf/`               | Text VCF reader/writer                          |
| `bcf`               | `pkg/htsgo/bcf/`               | Binary BCF v2.2 reader/writer + dict helpers    |
| `bgzf`              | `pkg/htsgo/bgzf/`              | BGZF reader/writer + `.gzi` index               |
| `bam`               | `pkg/htsgo/bam/`               | BAI (.bai) index format + `BuildBAI(*sam.BAMReader)` |
| `tabix`             | `pkg/htsgo/tabix/`             | Tabix (.tbi) + CSI (.csi) + binning helpers (Reg2bin, LinearTile) |
| `region`            | `pkg/htsgo/region/`            | `chr:start-end` parser + ResolveRegions (BAI-agnostic) |

All format packages live under `pkg/htsgo/`. The legacy
`pkg/bioformats/` directory, the `tools/bgzip/pkg/bgzip/` /
`tools/tabix/pkg/tabix/` shim paths, and the two in-samtools
sub-shims (`bai_shim.go`, `region_shim.go`) have all been
removed; every importer uses the `pkg/htsgo/` path directly.

Note: the BAM reader/writer themselves currently live in
`pkg/htsgo/sam/` (`sam.BAMReader`, `sam.BAMWriter`) — they moved
into htsgo via PR-A/B as part of the `sam` migration. The
`pkg/htsgo/bam` package added here only owns the **BAI index**
format primitives and the BAI builder. A future PR can split BAM
out of `sam/` into `bam/` if a cleaner SAM/BAM separation
becomes worthwhile (htslib itself keeps them together).

## Coming in subsequent PRs

| PR    | Brings in                                                       |
|-------|-----------------------------------------------------------------|
| ~~PR-B~~ | ~~`sam`, `vcf`, `bcf` move~~ **landed**                       |
| ~~PR-C~~ | ~~`bgzf` extracted from `tools/bgzip/pkg/bgzip/`~~ **landed** |
| ~~PR-D~~ | ~~`bam` (BAI primitives) extracted from `tools/samtools/`~~ **landed** |
| ~~PR-E~~ | ~~`tabix` (+ `.csi`) extracted from `tools/tabix/pkg/tabix/`~~ **landed** |
| ~~PR-F~~ | ~~`region/` — single source for `chr1:100-200` parsing~~ **landed** |
| PR-G  | already landed inline as the wave-21/22/23 vcftools BCF wiring  |
| PR-H (partial) | `tools/htsfile` CLI landed; `faidx` polish + region-iterator API still TBD |
| ~~PR-I~~ | ~~drop `pkg/bioformats/` and tool-package re-export shims~~ **landed** |
| PR-J+ | CRAM read/write (see `docs/CRAM_DESIGN.md`)                     |
| PR-K+ | `hfile` virtual filesystem (HTTP / S3 / GCS ranged reads)       |

## Conventions

- **Pure Go** — no cgo, no calls to libhts. CRAM codecs are the only
  third-party-dep zone (`pkg/htsgo/cram/codec/`), and only for rANS /
  LZMA primitives that aren't in the standard library. See
  `docs/CRAM_DESIGN.md` and the sanctioned-dep section of `CLAUDE.md`.
- **Tests move with code.** The shim packages have no tests; all
  coverage lives under `pkg/htsgo/<pkg>/*_test.go`.
- **Idiomatic Go API.** htslib's C API is the behavioural spec, not
  a literal one.
- Documented exported identifiers follow the project-wide
  complete-sentence docstring rule from `CLAUDE.md`.

## Treating htsgo as a tool

Per the convention from `tools/PORTING_STATUS.md`, htsgo gets:

- This `README.md` (the canonical inventory; supersedes
  `pkg/bioformats/README.md`).
- A row in `tools/PORTING_STATUS.md` once the migration completes
  (PR-I).
- CLI wrappers ported into `tools/`:
  - `tools/bgzip/` (exists) becomes a thin wrapper after PR-C.
  - `tools/tabix/` (exists) the same after PR-E.
  - `tools/htsfile/` (new in PR-H) is htslib's format-sniffer CLI.

## When to import from here

Always import from `pkg/htsgo/<pkg>` — the legacy
`pkg/bioformats/<pkg>` and tool-package shim paths
(`tools/bgzip/pkg/bgzip`, `tools/tabix/pkg/tabix`) were deleted
in PR-I.

bgzip-specific note: callers of the relocated BGZF code use
`bgzip "github.com/.../pkg/htsgo/bgzf"` as the canonical import
form — the htsgo target package is named `bgzf` but call sites
across the tree use the `bgzip.X` qualifier, so the alias
keeps every existing reference at zero diff cost. New code
should follow the same pattern.

The samtools subcommands now import `pkg/htsgo/bam` and
`pkg/htsgo/region` directly and qualify every reference
(`bam.BAIIndex`, `region.ParseRegion`, …) — the bare-name
sub-shims that PR-I had deferred were removed in a follow-up.
