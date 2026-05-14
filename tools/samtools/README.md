# samtools — pure-Go reimplementation (first slice)

`samtools` is a pure-Go reimplementation of the canonical
[samtools](https://github.com/samtools/samtools) command. This first slice
ships the foundational pieces needed by everything downstream:

- a shared SAM/BAM I/O library at
  `pkg/bioformats/sam` (text SAM reader/writer plus a binary BAM
  reader/writer on top of the in-tree
  `tools/bgzip/pkg/bgzip`);
- `samtools view` — print, filter, and convert SAM↔BAM records; and
- `samtools flagstat` — the classic 16-line alignment summary.

Subsequent PRs will add `sort`, `index`, `fastq`, `depth`, and `mpileup`.

`samtools` is pick #3 of the 2026 next-up list in
`analysis/tool_ranking_2026.md` — the most widely-used CLI in genomics, with
~6.5M bioconda downloads over its lifetime.

The implementation has **no third-party dependencies** — pure Go standard
library plus our existing in-tree libraries (`pkg/bioformats/{sam,iohelper}`,
`tools/bgzip/pkg/bgzip`, `tools/tabix/pkg/tabix`, `pkg/cliflag`).

## Build

```bash
go build ./tools/samtools/cmd/samtools
```

## Usage

```text
samtools view     [options] <in.bam|in.sam> [region ...]
samtools flagstat <in.bam|in.sam>
samtools help
samtools version
```

`in.bam`, `in.sam` and `in.sam.gz`/`in.bam` (BGZF-wrapped) are all auto-
detected; pass `-` to read from stdin.

### `samtools view`

| Short | Long                      | Description                                   |
|-------|---------------------------|-----------------------------------------------|
| `-b`  | `--bam`                   | Output BAM (default text SAM).                |
| `-h`  | `--with-header`           | Include the header alongside records.         |
| `-H`  | `--header-only`           | Print only the header.                        |
| `-c`  | `--count`                 | Print just the count of matching records.     |
| `-f`  | `--include-flags N`       | Keep records where ALL bits in N are set.     |
| `-F`  | `--exclude-flags N`       | Drop records where ANY bit in N is set.       |
| `-G`  | `--exclude-flags-all N`   | Drop only when ALL bits in N are set.         |
| `-q`  | `--min-mapq N`            | Minimum MAPQ.                                 |
| `-r`  | `--read-group ID`         | Keep records matching this RG.                |
| `-R`  | `--read-groups-file F`    | File of RG IDs (one per line).                |
| `-L`  | `--regions-file F`        | BED of regions (deferred — see Deviations).   |
| `-s`  | `--subsample F`           | Keep fraction `F`, or `<seed>.<frac>`.        |
| `-o`  | `--output PATH`           | Output file (default stdout).                 |
| `-T`  | `--reference FASTA`       | Accepted; CRAM is not supported in v1.        |
| `-@`  | `--threads N`             | Accepted; single-threaded in v1.              |
|       | `--no-PG`                 | Suppress `@PG` line emission.                 |
|       | `--help`                  | Show help.                                    |
|       | `--version`               | Show version.                                 |

Examples:

```bash
# Count primary mapped reads.
samtools view -c -F 0x900 -f 0x000 in.bam

# Pull a read group out and re-emit as BAM.
samtools view -b -r rg1 -o out.bam in.bam

# Print just the header.
samtools view -H in.bam
```

### `samtools flagstat`

Emits the classic 16-line summary:

```text
N + 0 in total (QC-passed reads + QC-failed reads)
N + 0 primary
N + 0 secondary
N + 0 supplementary
N + 0 duplicates
N + 0 primary duplicates
N + 0 mapped (XX.XX% : N/A)
N + 0 primary mapped (XX.XX% : N/A)
N + 0 paired in sequencing
N + 0 read1
N + 0 read2
N + 0 properly paired (XX.XX% : N/A)
N + 0 with itself and mate mapped
N + 0 singletons (XX.XX% : N/A)
N + 0 with mate mapped to a different chr
N + 0 with mate mapped to a different chr (mapQ>=5)
```

Each counter is split QC-passed + QC-failed via the 0x200 flag bit.

## Deviations from upstream samtools

This is the **first slice**; the following are intentionally out of scope and
will land in follow-up PRs.

- **Region queries** (`chr:start-end` on the CLI, or `-L regions.bed`) require
  a `.bai` index. BAI indexing is not yet implemented; passing a region or
  `-L` causes a clear error ("region-query support requires .bai indexing —
  not yet implemented"). Linear streaming filters (flags, MAPQ, RG) work
  fine without an index.
- **CRAM** is not supported. The `-T/--reference` flag is accepted (so
  pipelines passing it through do not break) but has no effect.
- **Multi-threading.** `-@/--threads` is accepted but the v1 pipeline is
  single-threaded.
- **`--no-PG`** is accepted but has no observable effect — the v1 view
  pipeline never injects a `@PG` line into the header on its own.
- `flagstat` always reports `N + 0` (i.e. all records counted on the
  QC-passed side); QC-failed totals are derived from the 0x200 flag bit so
  this is correct, but the `0` after `+` reflects that no extra QC-fail
  filtering is performed in v1.

## What this unblocks

`pkg/bioformats/sam` is the foundation for everything htslib-shaped that
follows:

- `samtools sort` (in-memory and external-merge variants).
- `samtools index` (BAI) — and by extension, region queries in
  `samtools view`.
- `samtools depth` and `mosdepth` (the latter is pick #9 of the 2026
  ranking).
- `samtools fastq` and `bam2fq` conversions.
- Downstream callers (`bcftools mpileup`, GATK-style metrics, etc.).

## References

- [SAM/BAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf).
- `reference_code/samtools/sam.c`, `bam.c` (vendored upstream — read-only).
