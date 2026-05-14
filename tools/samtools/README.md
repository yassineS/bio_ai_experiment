# samtools — pure-Go reimplementation

`samtools` is a pure-Go reimplementation of the canonical
[samtools](https://github.com/samtools/samtools) command. Slices landed so far:

- a shared SAM/BAM I/O library at
  `pkg/bioformats/sam` (text SAM reader/writer plus a binary BAM
  reader/writer on top of the in-tree
  `tools/bgzip/pkg/bgzip`);
- `samtools view` — print, filter, convert SAM↔BAM records, **with region
  queries (`chr:start-end`) backed by a sibling `.bai` index when one
  exists**;
- `samtools sort` — external-merge sort by coordinate, name, natural name,
  or aux tag;
- `samtools index` — build a BAI index for a coordinate-sorted BAM; and
- `samtools flagstat` — the classic 16-line alignment summary.

Subsequent PRs will add `fastq`, `depth`, `mpileup`, and CSI indexing.

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
samtools sort     [options] <in.bam|in.sam>
samtools index    [options] <in.sorted.bam>
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

### `samtools sort`

External-merge sort that bounds memory use by spilling chunked-and-sorted
runs to temporary BAM files, then doing a k-way heap merge for the final
output. The sort is stable.

| Short | Long                       | Description                              |
|-------|----------------------------|------------------------------------------|
| `-o`  | `--output PATH`            | Output path (default stdout).            |
| `-O`  | `--output-fmt {bam,sam}`   | Force output format.                     |
| `-n`  | `--by-name`                | Sort by QName (lexicographic).           |
| `-N`  | `--by-natural-name`        | Sort by QName with numeric-run ordering. |
| `-t`  | `--by-tag TAG`             | Sort by an Aux tag value.                |
| `-m`  | `--max-mem N[K|M|G]`       | Per-shard memory budget (default 768M).  |
| `-T`  | `--tmpdir PREFIX`          | Tmpfile prefix for spill files.          |
| `-l`  | `--compress-level N`       | BGZF deflate level 0..9.                 |
| `-@`  | `--threads N`              | Accepted; single-threaded in v1.         |
|       | `--no-PG`                  | No `@PG` injection (v1 never injects).   |

```bash
samtools sort -o out.bam in.bam            # coordinate sort
samtools sort -n -o by-name.bam in.bam     # query-name sort
samtools sort -t NM -o by-nm.bam in.bam    # sort by NM tag value
```

The output's `@HD SO:` field is rewritten to match the chosen order
(`coordinate`, `queryname`, or `unknown` for tag-sorted output).

### `samtools index`

Builds a `.bai` index for a coordinate-sorted BAM. The on-disk format is
exactly the BAI binary layout used by htslib: per-reference bin lists with
two-uint64 chunks, a metadata pseudo-bin (37450) carrying the per-ref
(firstVOff, lastVOff) and (mapped, unmapped) tallies, a linear index of
16-Kbp tile entries, and an optional `n_no_coor` trailer for unplaced
records.

| Short | Long                  | Description                                  |
|-------|-----------------------|----------------------------------------------|
| `-b`  | `--bai`               | Emit `.bai` (default).                       |
| `-c`  | `--csi`               | Emit `.csi` — **NOT YET IMPLEMENTED**.       |
|       | `--csi-min-shift N`   | CSI shift (accepted with `-c`).              |
| `-o`  | `--output PATH`       | Index output path (default `<input>.bai`).   |
| `-@`  | `--threads N`         | Accepted; single-threaded.                   |

```bash
samtools sort -o sorted.bam in.bam
samtools index sorted.bam       # writes sorted.bam.bai
samtools view sorted.bam chr1:1000-2000
```

When `samtools view` sees a region argument and finds a sibling `.bai` it
uses the BAI bins and linear index to seek directly to the relevant BGZF
chunks; otherwise it falls back to a full linear scan with a stderr
warning.

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

The following are intentionally out of scope or deferred and will land in
follow-up PRs.

- **CSI indexing.** `samtools index -c/--csi` is accepted on the CLI but
  surfaces a clear error: "CSI output (-c/--csi) is not yet implemented;
  v1 emits BAI only". CSI is only needed for chromosomes longer than
  ~512 Mb, which excludes every common reference genome.
- **`-L regions.bed`.** Region querying via a BED file is still deferred;
  `samtools view` will print a warning to stderr and fall back to a
  whole-file scan when `-L` is provided. Single-region specifiers
  (`chr1:1000-2000`) are fully supported.
- **CRAM** is not supported. The `-T/--reference` flag is accepted (so
  pipelines passing it through do not break) but has no effect.
- **Multi-threading.** `-@/--threads` is accepted by every subcommand but
  the v1 pipelines are single-threaded.
- **`--no-PG`** is accepted but has no observable effect — none of our
  subcommands inject a `@PG` line into the header.
- `samtools sort` produces deterministic output (stable sort + QName
  tie-break) which differs subtly from upstream samtools, which uses an
  unstable parallel sort and may permute records that share their primary
  key.
- `flagstat` always reports `N + 0` (i.e. all records counted on the
  QC-passed side); QC-failed totals are derived from the 0x200 flag bit so
  this is correct, but the `0` after `+` reflects that no extra QC-fail
  filtering is performed in v1.

## What this unblocks

With `sort` + `index` + region queries landed, the natural next slices are:

- `samtools mpileup` and `bcftools mpileup` (need indexed seek to walk a
  region per ref).
- `samtools fastq` / `bam2fq` (works on coordinate- or name-sorted input).
- `samtools depth` and `mosdepth` (the latter is pick #9 of the 2026
  ranking — uses sorted+indexed BAM throughout).
- CSI indexing for chromosomes longer than 512 Mb.
- Multi-threaded sort and index, once a profiling pass shows where the
  single-threaded baseline bottlenecks.

## References

- [SAM/BAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf).
- `reference_code/samtools/sam.c`, `bam.c` (vendored upstream — read-only).
