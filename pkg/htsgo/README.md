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

Old `pkg/bioformats/{iohelper,fasta,fastq,bed,gff,sam,vcf,bcf}/`
paths still resolve via tiny re-export shims; new code should import
the htsgo path directly. The shims will be deleted in PR-I.

## Coming in subsequent PRs

| PR    | Brings in                                                       |
|-------|-----------------------------------------------------------------|
| ~~PR-B~~ | ~~`sam`, `vcf`, `bcf` move~~ **landed**                       |
| PR-C  | `bgzf` extracted from `tools/bgzip/pkg/bgzip/`                  |
| PR-D  | `bam` (+ `.bai`) extracted from `tools/samtools/pkg/samtools/` |
| PR-E  | `tabix` (+ `.csi`) extracted from `tools/tabix/pkg/tabix/`     |
| PR-F  | `region/` — single source for `chr1:100-200` parsing            |
| PR-G  | already landed inline as the wave-21/22/23 vcftools BCF wiring  |
| PR-H  | `faidx` polish + region-iterator API + `tools/htsfile` CLI      |
| PR-I  | drop `pkg/bioformats/` and tool-package re-export shims         |
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

## When to import from here vs. tool-package paths

All format packages (`iohelper`, `fasta`, `fastq`, `bed`, `gff`,
`sam`, `vcf`, `bcf`) are now under `pkg/htsgo/<pkg>` — prefer those
import paths; the matching `pkg/bioformats/<pkg>` paths are
deprecated re-export shims that PR-I will delete.

The in-tool packages still at their original locations
(`tools/bgzip/pkg/bgzip`, `tools/samtools/pkg/samtools`,
`tools/tabix/pkg/tabix`) migrate in PRs C–E per the table above.
