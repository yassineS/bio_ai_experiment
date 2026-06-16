# samtools — pure-Go reimplementation

`samtools` is a pure-Go reimplementation of the canonical
[samtools](https://github.com/samtools/samtools) command. Slices landed so far:

- a shared SAM/BAM I/O library at
  `pkg/htsgo/sam` (text SAM reader/writer plus a binary BAM
  reader/writer on top of the in-tree
  `pkg/htsgo/bgzf`);
- `samtools view` — print, filter, convert SAM↔BAM records, **with region
  queries (`chr:start-end`) backed by a sibling `.bai` index when one
  exists**;
- `samtools sort` — external-merge sort by coordinate, name, natural name,
  or aux tag;
- `samtools index` — build a `.bai` (default) or `.csi` (`-c`) index for
  a coordinate-sorted BAM;
- `samtools flagstat` — the classic 16-line alignment summary;
- `samtools depth` — per-position depth across one or more BAMs; and
- `samtools fastq` (and the `bam2fq` alias) — convert SAM/BAM to FASTQ.

Subsequent PRs will continue to fill in the `mpileup` subcommand tail.

`samtools` is pick #3 of the 2026 next-up list in
`analysis/tool_ranking_2026.md` — the most widely-used CLI in genomics, with
~6.5M bioconda downloads over its lifetime.

The implementation has **no third-party dependencies** — pure Go standard
library plus our existing in-tree libraries (`pkg/htsgo/{sam,iohelper}`,
`pkg/htsgo/bgzf`, `pkg/htsgo/tabix`, `pkg/cliflag`).

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
samtools depth    [options] <in1.bam> [<in2.bam> ...]
samtools fastq    [options] <in.bam|in.sam>
samtools bam2fq   [options] <in.bam|in.sam>   # alias for fastq
samtools markdup  [options] <in.bam> <out.bam>
samtools stats    [options] <in.bam|in.sam>
samtools tview    [options] <in.bam|in.cram> [ref.fasta]
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
| `-L`  | `--regions-file F`        | Keep records overlapping any BED interval.    |
| `-M`  | `--use-multi-region-iterator` | Accepted (we always run the full intersection). |
| `-s`  | `--subsample F`           | Keep fraction `F`, or `<seed>.<frac>`.        |
| `-o`  | `--output PATH`           | Output file (default stdout).                 |
| `-T`  | `--reference FASTA`       | Reference FASTA (used for CRAM resolution).   |
| `-@`  | `--threads N`             | Parallel BGZF compression of the BAM output.  |
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
| `-c`  | `--csi`               | Emit `.csi` (for references > 512 Mb).       |
| `-m`  | `--min-shift N`       | CSI bin-hierarchy `min_shift` (default 14;   |
|       |                       | `--csi-min-shift` is accepted as an alias).  |
| `-o`  | `--output PATH`       | Index output path (default `<input>.bai`,    |
|       |                       | or `<input>.csi` with `-c`).                 |
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

### `samtools depth`

Prints per-position depth (1-based coordinates) across one or more BAM/SAM
inputs. The output is `chrom\tpos\tdepth1[\tdepth2 ...]\n`, with the
depths of each input file in positional order. All inputs must share the
same `@SQ` ordering.

| Short | Long                 | Description                                  |
|-------|----------------------|----------------------------------------------|
| `-a`  | `--all`              | Emit zero-depth positions inside covered regions. |
| `-A`  | `--all-trans`        | Emit every position of every reference.      |
| `-r`  | `--region chr[:S-E]` | Limit to region (repeatable).                |
| `-b`  | `--bed FILE`         | Limit to BED regions.                        |
| `-q`  | `--min-mapq N`       | Skip reads with MAPQ below N.                |
| `-Q`  | `--min-baseq N`      | Skip bases with quality below N.             |
| `-l`  | `--min-readlen N`    | Skip reads shorter than N query bases.       |
| `-f`  | `--include-flags N`  | Require ALL flag bits in N to be set.        |
| `-F`  | `--exclude-flags N`  | Drop reads with ANY of these flag bits (default `0x4`). |
| `-d`  | `--max-depth N`      | Cap reported depth (`0` = no cap).           |
| `-o`  | `--output PATH`      | Output path (default stdout).                |
| `-@`  | `--threads N`        | Accepted; single-threaded.                   |

Depth is incremented for every reference base covered by a `M`, `=`, or
`X` CIGAR operation; `I`/`S`/`H`/`N`/`P` ops do not contribute. The
default `-F 0x4` filters out unmapped reads, matching upstream.

```bash
samtools depth -a -r chr1:1000-2000 sorted.bam
samtools depth -b regions.bed sample1.bam sample2.bam
```

### `samtools fastq` (alias `bam2fq`)

Converts a SAM/BAM file to FASTQ. For paired output, name-sorted input is
required — coordinate-sorted input falls back to writing every record to
`-o` (or singletons) with a stderr warning.

| Short | Long                   | Description                                |
|-------|------------------------|--------------------------------------------|
| `-1`  | `--read1 FILE`         | Output for first-in-pair reads.            |
| `-2`  | `--read2 FILE`         | Output for second-in-pair reads.           |
| `-0`  | `--read-orphan FILE`   | Paired reads where 0x40/0x80 are both set or both unset. |
| `-s`  | `--singleton FILE`     | Output for unpaired reads.                 |
| `-o`  | `--output FILE`        | Default sink (interleaved if `-1/-2` unset). |
| `-N`  | `--output-name`        | Always append `/1` or `/2` to read names.  |
| `-n`  | `--no-suffix`          | Never append `/1` `/2`.                    |
| `-f`  | `--include-flags N`    | Required flag bits.                        |
| `-F`  | `--exclude-flags N`    | Excluded flag bits (default `0x900`).      |
| `-G`  | `--exclude-flags-all N`| Drop only when ALL bits match.             |
| `-T`  | `--add-tags TAGS`      | Comma-separated aux tags to append to the read description. |
| `-t`  | `--no-CO`              | Accepted; we never emit `@CO` lines.       |
| `-c`  | `--compress-level N`   | Gzip level for `.gz` outputs.              |
| `-O`  | `--use-qq`             | Use `OQ` aux tag for quality when present. |
|       | `--threads N`          | Accepted; single-threaded.                 |

Reverse-strand records (`FLAG & 0x10`) have their SEQ reverse-complemented
and QUAL reversed back to original-read orientation. Paired suffixes
(`/1`, `/2`) are appended unless `-n` is set or the QNAME already ends
with them.

```bash
# Split paired data into two files (name-sorted input):
samtools fastq -1 r1.fq -2 r2.fq -s singleton.fq -0 orphan.fq paired.bam

# Interleaved FASTQ on stdout:
samtools fastq paired.bam | bgzip > paired.fq.gz

# Compressed split outputs:
samtools fastq -1 r1.fq.gz -2 r2.fq.gz -c 6 paired.bam

# Append the NM tag to each read description:
samtools fastq -T NM paired.bam > tagged.fq
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

### `samtools markdup`

Mark or remove PCR duplicates using a two-pass `pair_key`/`single_key`
streaming algorithm that mirrors upstream's `bam_markdup.c`. The
algorithm walks the input twice: pass 1 builds the duplicate buckets,
pass 2 re-emits records with the 0x400 flag set on the non-chosen
members of each bucket. The "best" record is the one whose
`calc_score(read) + ms_tag(read)` is highest (`calc_score` = sum of base
qualities ≥ 15).

| Short | Long             | Description                                |
|-------|------------------|--------------------------------------------|
| `-r`  | `--remove-dups`  | Drop duplicates from output (vs flag).     |
| `-d`  | `--max-dist N`   | Optical-dup distance (accepted, not impl). |
| `-s`  | `--mode {t/s/tp}`| Key mode: template (default) / sequence.   |
| `-T`  | `--tmpdir PATH`  | Accepted; v1 streams in memory.            |
| `-l`  | `--max-len N`    | Max read length considered (default 300).  |
|       | `--include-flags N` | Require ALL bits set.                   |
|       | `--exclude-flags N` | Drop records with ANY bit set.          |
| `-c`  | `--clear-tags`   | Strip pre-existing `do`/`dt`/`mc` tags.    |
| `-t`  | `--add-tag`      | Write `do:Z:<winner-qname>` on duplicates. |
| `-d`  | `--max-dist N`   | Optical-duplicate max pixel distance.       |
| `-s`  | `--stats`        | Emit the Picard-style duplicate-stats report. |
| `-S`  |                  | Mark supplementary/secondary of a duplicate too. |
| `-@`  | `--threads N`    | Parallel BGZF compression of the BAM output. |
| `-o`  | `--output PATH`  | Output BAM (default stdout).               |
|       | `--no-PG`        | Suppress `@PG` line emission.              |

Byte-parity validated against upstream's `test/markdup/5_markdup.sam` and
`6_remove_dups.sam`; flag-parity on `18_primary_duplicate_count.sam`.
Optical-duplicate detection (`-d`), the stats report (`-s`/`-S`), and
supplementary marking (`-S`) are implemented; see `PARITY_VALIDATION.md`.

### `samtools stats`

Emit the upstream-compatible Summary Numbers (SN) block plus useful
read-length, MAPQ, and insert-size histograms. The SN block is
byte-faithful with upstream for the 6 fixtures we exercise in
`stats_test.go` (fixtures 1, 2, 5, 7, 8, 10 from
`reference_code/samtools/test/stat/`).

| Short | Long                  | Description                                 |
|-------|-----------------------|---------------------------------------------|
| `-r`  | `--ref-seq FASTA`     | Reference FASTA; GCD GC + MPC from the ref. |
| `-c`  | `--coverage SPEC`     | Coverage bin spec (default `1,1000,1`).     |
|       | `--GC-depth N`        | GC-depth bin width (default 20000).         |
| `-l`  | `--required-flag N`   | Require ALL bits set.                       |
| `-F`  | `--filtering-flag N`  | Drop records with ANY bit set.              |
| `-d`  | `--max-depth N`       | Cap depth (placeholder).                    |
| `-q`  | `--trim-quality N`    | BWA-style 3'-end quality-trim threshold.    |
|       | `--min-mapq N`        | Skip records with MAPQ < N.                 |
|       | `--remove-dups`       | Drop duplicate-flagged records.             |
|       | `--remove-overlaps`   | Accepted; no-op in v1.                      |
| `-i`  | `--insert-size N`     | Max insert size for IS section (8000).      |
| `-x`  | `--sparse`            | Suppress IS rows that have no insertions.   |
| `-t`  | `--target-regions F`  | Restrict stats to a target-regions file.    |
| `-g`  | `--cov-threshold N`   | Coverage threshold for the target SN line.  |
|       | `--ref-stats`         | Emit the RFS reference-statistics section.  |
|       | `--ref-stats-chunk N` | RFS reference-fetch chunk width (MB).       |
| `-@`  | `--threads N`         | Accepted; single-threaded.                  |
| `-o`  | `--output PATH`       | Output path (default stdout).               |

Sections emitted: the **CHK** checksum block, **SN**, **FFQ/LFQ**,
**MPC** (with `--ref-seq`), **GCF/GCL**, **GCC/GCT**, **FBC/FTC/LBC/LTC**,
the per-barcode **`<tag>C`/`<tag>Q`** tables (emitted automatically when
the BC/CR/OX/RX aux tags are present), **RL/FRL/LRL**, **MAPQ**, **IS**,
**IC/ID**, **COV**, **GCD** and **RFS** (with `--ref-stats`). Every
output section is implemented; `samtools stats` is at full 1:1 section
parity. `-x/--sparse` thins all-zero IS rows only.

### `samtools tview`

An alignment viewer with three display modes. `-d T` (plain text) and `-d H`
(HTML) are pipeline-friendly; `-d C` (and the bare default on a terminal) is the
**interactive** viewer. All three render the same character grid upstream's
`bam_tview.c` builds for a region — a ruler line, the reference bases, a
per-position consensus, then one row per read with non-overlapping reads packed
greedily onto shared rows. The text mode emits the grid as plain characters (no
ANSI escapes — upstream only colourises a TTY); the HTML mode wraps each cell in
a coloured `<span>`, a direct port of `bam_tview_html.c`. Both non-interactive
modes are verified **byte-for-byte** against the vendored upstream binary
(`TestTviewLiveParity`).

```bash
samtools tview -d T -p chr1:10000 -w 160 aln.bam ref.fa   # text grid
samtools tview -d H -p chr1:10000 aln.bam ref.fa > view.html
samtools tview aln.bam ref.fa                             # interactive (TTY)
```

| Short | Long              | Description                                                      |
|-------|-------------------|-----------------------------------------------------------------|
| `-d`  | `--display T\|H\|C` | Text / HTML / interactive (curses) output. Default: interactive on a TTY. |
| `-p`  | `--position REG`  | Start at this region (`chr` or `chr:pos`).                       |
| `-w`  | `--width INT`     | Display width in columns (default: terminal width, or 80 for T/H).|
| `-s`  | `--sample STR`    | Show only reads from this `@RG` sample (`SM:`) or read-group ID. |
| `-T`  | `--reference FA`  | Reference FASTA (also accepted positionally as `ref.fasta`).     |
| `-i`  | `--hide-inserts`  | Hide insertion columns (default expands them).                  |

Rendering matches upstream exactly: a base equal to the reference shows `.`
(forward) / `,` (reverse); a mismatch shows the read base (UPPER forward / lower
reverse); deletions show `*`; reference skips `>` (forward) / `<` (reverse).
The consensus line uses the same bam2bcf `bcf_call_glfgen` + MAQ error-model
call upstream uses, and the read-row packing reproduces the `bam_lpileup.c`
level-pool algorithm (a freed row is reusable after the upstream two-column gap,
lowest free row first). BAM and CRAM inputs are accepted; CRAM uses `-T` as the
decode reference.

**Interactive mode (`-d C`)** is implemented in **pure Go (Linux, no
ncurses)**. It puts the terminal into raw (cbreak/no-echo) mode via the
`TCGETS`/`TCSETS` termios ioctls, queries the window size with `TIOCGWINSZ`
(falling back to 80×40), and reuses the exact same frame renderer as `-d T`,
redrawing with ANSI escapes after each keystroke. The original termios is
restored on exit (and on panic, via `defer`). Key bindings (ported from
`bam_tview_curses.c`):

| Key | Action |
|-----|--------|
| `←`/`h`, `→`/`l` | scroll one column |
| `↑`/`k`, `↓`/`j` | scroll read rows |
| `H` / `L` | page 20 columns left / right |
| space / backspace | page one screen right / left |
| `Ctrl-H` / `Ctrl-L` | jump 1000 columns left / right |
| `0` / Home | jump to the start of the contig |
| `g` (or `/`) | prompt for a `chr:pos` region and jump |
| `m` / `b` / `n` / `N` | colour mode: map-q / base-q / nucleotide / none |
| `.` | toggle base-vs-dot |
| `i` | toggle insertion columns |
| `r` | toggle by-read-name colour |
| `?` | help screen (any key returns) |
| `q` / `Esc` / `Ctrl-C` | quit |

The OS-specific termios/ioctl code is isolated in `tview_tty_linux.go`; the
key→action map and action→state machine are pure functions (`tview_interactive.go`)
unit-tested without a TTY. **Piped `-d C` exits with a clear message** (use
`-d T` / `-d H` in pipelines) rather than garbling the stream. On **non-Linux
platforms** `-d C` reports that interactive mode requires Linux (the rest of the
tool still builds and runs).

## Deviations from upstream samtools

The following are intentionally out of scope or deferred and will land in
follow-up PRs.

- **`-L regions.bed`.** `samtools view -L` is supported: a per-chromosome
  interval tree is built from the BED and records are kept only when their
  `[Pos, Pos+refLen)` half-open range overlaps at least one BED interval
  on the record's reference. The walk is always a linear scan (no
  `.bai` shortcut yet); `-M`/`--use-multi-region-iterator` is accepted but
  produces an identical record set, so we always run the full intersection.
- **CRAM** read and write are supported throughout (rANS 4x8/4x16 in-tree
  pure Go; `ulikunitz/xz` confined to the LZMA block codec; X_EXT bzip2
  encode in-tree). The `-T/--reference` FASTA is consumed for CRAM
  reference resolution. v2.1 decode and v3.0/v3.1 read+write are validated
  byte-for-byte against live `samtools`; the only residuals are the network
  REF_PATH/EBI fetch and CRAM v4.0 (spec not final).
- **Multi-threading.** `-@/--threads` drives a genuine parallel BGZF
  compressor for the BAM-writing subcommands (`view`, `sort`, `markdup`)
  **and a parallel BGZF decompressor on the input path** for `view`,
  `flagstat`, `idxstats` (no-index scan), `stats`, and `depth`: a
  BGZF-wrapped BAM input is inflated block-by-block across `-@` worker
  goroutines via `bgzf.MultiReader`, which reorders the decoded blocks so the
  byte stream — and every record and output line — is byte-identical for any
  thread count; only decode throughput changes. The remaining deferral is
  parallel *input* **CRAM** slice decode (CRAM uses its own container framing,
  not BGZF blocks, so it is decoded single-threaded) — perf-only, not a
  feature gap.
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
- `samtools depth` requires identical `@SQ` ordering across all inputs and
  surfaces an explicit error when they differ — matching upstream's
  positional-output contract. Depth contributions follow the conservative
  rule "only `M`/`=`/`X` CIGAR ops count for reference coverage"; `D`/`N`
  advance the reference but do not add depth (consistent with upstream).
- `samtools fastq` requires **name-sorted** input when `-1`/`-2` are used
  for paired-split output. Coordinate-sorted input is detected from the
  `@HD SO:` field and falls back to writing every record through `-o` (or
  to `-s` if only that is configured) with a warning to stderr. A second
  pass that pairs records on disk for coordinate-sorted input may follow
  in a later slice.

## What this unblocks

With `view`, `sort`, `index`, `flagstat`, `depth`, and `fastq` landed,
the natural next slices are:

- `samtools mpileup` and `bcftools mpileup` (depth's per-position scan is
  most of the iteration scaffold).
- `mosdepth` (pick #9 of the 2026 ranking — coverage-summary tool that
  builds on top of indexed BAM iteration).
- A coordinate-sorted code-path for `fastq -1/-2` that pairs records on
  disk (currently this configuration falls back to interleaved output).
- Multi-threaded sort, index, depth, and fastq once a profiling pass
  shows where the single-threaded baselines bottleneck.

## References

- [SAM/BAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf).
- `reference_code/samtools/sam.c`, `bam.c` (vendored upstream — read-only).
