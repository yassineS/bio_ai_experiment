# bedintersect - BED Interval Intersection Finder

A fast, drop-in re-implementation of `bedtools intersect`, written in Go.

## Features

- **Full `bedtools intersect` flag parity**: `-wa/-wb/-wo/-wao/-loj/-u/-c/-C/-v`,
  `-f/-F/-r/-e`, `-s/-S`, `-split`, multiple `-b` databases
  (`-names`/`-filenames`/`-sortout`), `-sorted`/`-g`, `-header`, `-nonamecheck` —
  validated byte-for-byte against the live upstream binary.
- **Multiple input formats**: BED, GFF, VCF (incl. structural-variant spans),
  BAM and CRAM, with gzip/BGZF auto-decompression (including on stdin).
- **Fast performance**: Chromosome indexing with an optional interval tree for
  very large B files.
- **bedintersect extensions**: distance to nearest feature (`-d`), closest
  feature (`-k`), and a `--stats` summary.
- **Built-in gzip support**: Automatic handling of `.gz` files.
- **Streaming**: File A is processed in streaming fashion (memory-efficient).

## Installation

### From Source

```bash
cd tools/bedintersect
go build ./cmd/bedintersect
```

### Using Go Install

```bash
go install github.com/yassineS/bio_ai_experiment/tools/bedintersect/cmd/bedintersect@latest
```

## Usage

### Basic Usage

```bash
bedintersect -a genes.bed -b peaks.bed > overlaps.bed
```

### Options

- `-a FILE` - Input file A (BED/GFF/VCF/BAM; `-` for stdin). Required.
- `-b FILE [FILE ...]` - Input file(s) B. `-b` may be followed by multiple
  database files (each B record is then prefixed with a DB-id column in
  `-wb`/`-loj`/`-wao` output). Required.
- `-abam FILE` / `-ibam FILE` - Aliases for `-a` with a BAM file.
- `-o, --output FILE` - Output file (default: stdout)
- `-f, --fraction-a NUM` - Minimum overlap as a fraction of A, in `(0.0, 1.0]`
  (default ~1bp).
- `-F, --fraction-b NUM` - Minimum overlap as a fraction of B, in `(0.0, 1.0]`
  (default ~1bp).
- `-r, --reciprocal` - Require the fraction overlap be reciprocal for A AND B.
- `-e` - Require the minimum fraction be satisfied for A OR B (rather than both).
- `-s, --strand` - Only report hits on the same strand.
- `-S` - Only report hits on the opposite strand. (Matches `bedtools intersect`;
  note this differs from older versions of this tool, where `-S` meant
  "statistics" — that is now the long flag `--stats` only.)
- `-v, --invert` - Report A entries with NO overlap with B.
- `-wa, --write-a` - Write the original A entry for each overlap.
- `-wb, --write-b` - Write the original B entry for each overlap (with `-wa`:
  write A and B side-by-side per overlap).
- `-u` - Write each original A entry once if any overlap is found.
- `-c, --count` - Report the total count of B overlaps for each A.
- `-C` - Report the overlap count with each B file separately (one line per B
  file; with multiple B files a DB-id column is included).
- `-loj` - Left outer join: report every A; B is null when there is no overlap.
- `-wo` - Write A, B and the number of overlapping bases per overlap.
- `-wao` - Like `-wo`, but also report A (with null B and overlap `0`) when there
  is no overlap.
- `-split` - Treat split BAM/BED12 entries as distinct intervals (block-aware
  overlap and overlap-base counting).
- `-names NAME [NAME ...]` - Aliases for each B file (printed instead of a
  numeric file id in the DB-id column).
- `-filenames` - Print each B file's name instead of a numeric file id.
- `-sortout` - Sort each A record's DB hits by position across all B files.
- `-header` - Print the A file's header lines verbatim before the results.
- `-sorted` - Validate that the inputs are coordinate-sorted (errors otherwise),
  mirroring upstream. The output order is unchanged (it already matches
  upstream's bin order whether or not `-sorted` is given).
- `-g FILE` - Genome file fixing the chromosome order for `-sorted` validation.
- `-nonamecheck` - Suppress the chromosome naming-convention warning.
- `-bed` - With BAM/CRAM input, write output as BED instead of the default
  binary alignments. By default a BAM/CRAM query (`-a`/`-abam`/`-ibam`) writes
  the intersecting alignments back out as BAM, or as CRAM when a CRAM reference
  is available (see BAM/CRAM output below); `-bed` forces BED12 text output
  instead.
- `--cram-ref FILE` - CRAM reference FASTA. A CRAM query writes **CRAM** output
  (rather than BAM) only when this flag, or the `CRAM_REFERENCE` environment
  variable, names a reference, matching upstream.
- `-ubam` - Accepted for compatibility (requests uncompressed BAM). Upstream's
  format choice still follows the CRAM reference, so a CRAM query with a
  reference stays CRAM under `-ubam`; our BAM output is always BGZF-compressed.
- `-m, --min-overlap INT` - Minimum overlap in bp (bedintersect extension,
  default 1).
- `-d, --distance` - Report distance to nearest B feature (bedintersect
  extension).
- `-k, --closest` - Output closest B feature for each A (bedintersect extension).
- `-t, --tree` - Use an interval tree for large B files (bedintersect extension).
- `--stats` - Print summary statistics to stderr (bedintersect extension; single
  B file only).
- `-h, --help` - Show help message.
- `--version` - Show version information.

### Examples

#### Find overlapping regions (default: intersection)

```bash
bedintersect -a genes.bed -b peaks.bed > overlaps.bed
```

Output: The overlapping portion of each pair:

```
chr1 150 200
chr1 350 400
```

#### Report original A entries

```bash
bedintersect -a genes.bed -b peaks.bed -wa > genes_with_peaks.bed
```

Output: Original gene coordinates that have peaks:

```
chr1 100 200
chr1 300 400
```

#### Report B entries that overlap A

```bash
bedintersect -a genes.bed -b peaks.bed -wb > peaks_in_genes.bed
```

#### Count overlaps per A interval

```bash
bedintersect -a genes.bed -b peaks.bed -c > gene_peak_counts.bed
```

Output: Each gene with count in name field:

```
chr1 100 200 3
chr1 300 400 1
```

#### Find A intervals with no B overlap

```bash
bedintersect -a genes.bed -b peaks.bed -v > genes_without_peaks.bed
```

#### Require minimum overlap

```bash
bedintersect -a genes.bed -b peaks.bed -m 50 > overlaps.bed
```

Only reports overlaps of at least 50bp.

#### Require fractional overlap of A

```bash
bedintersect -a genes.bed -b peaks.bed -f 0.8 > overlaps.bed
```

Requires that at least 80% of each gene overlaps a peak.

#### Require fractional overlap of B

```bash
bedintersect -a genes.bed -b peaks.bed -F 0.5 > overlaps.bed
```

Requires that at least 50% of each peak overlaps a gene.

#### Strand-specific intersection

```bash
bedintersect -a genes.bed -b peaks.bed -s > overlaps.bed
```

Only reports overlaps on the same strand (requires strand column in both files).

#### Require reciprocal overlap

```bash
bedintersect -a genes.bed -b peaks.bed -f 0.5 -F 0.5 -r > overlaps.bed
```

Requires that at least 50% of both the gene and the peak overlap each other.

#### Report distance to nearest feature

```bash
bedintersect -a genes.bed -b peaks.bed -d > distances.bed
```

Outputs each gene with the distance to its nearest peak in the name field (0 for overlapping).

#### Report closest feature

```bash
bedintersect -a genes.bed -b peaks.bed -k > closest_peaks.bed
```

Outputs the closest peak for each gene.

#### Use interval tree for large files

```bash
bedintersect -a genes.bed -b large_database.bed -t > overlaps.bed
```

Uses an interval tree data structure for improved performance with very large B files.

#### Show statistics

```bash
bedintersect -a genes.bed -b peaks.bed --stats > overlaps.bed
```

(`--stats` is a bedintersect extension; in older versions this was `-S`, but
`-S` now matches `bedtools intersect`'s opposite-strand meaning.)

Outputs to stderr:

```
Intervals in A: 1000
Intervals in B: 500
A intervals with hits: 450
A intervals with no hits: 550
Total overlaps: 680
```

#### Gzip support

```bash
bedintersect -a genes.bed.gz -b peaks.bed.gz > overlaps.bed
bedintersect -a genes.bed -b peaks.bed | gzip > overlaps.bed.gz
```

## Input Format

Standard BED format with at least 3 columns:

```
chr1 100 200
chr1 300 400
```

- Tab-delimited
- Minimum 3 fields: chromosome, start, end
- Optional: name, score, strand, etc.
- 0-based, half-open coordinates [start, end)
- Does not need to be sorted

## Output Format

Depends on options:

### Default (intersection coordinates)

```
chr1 150 200
chr1 350 400
```

### With -wa (original A)

```
chr1 100 200
chr1 300 400
```

### With -wb (B entries)

```
chr1 150 250
chr1 350 450
```

### With -c (counts)

```
chr1 100 200 2
chr1 300 400 1
```

## Algorithm

1. **Read B file** completely into memory
2. **Index B intervals** by chromosome
3. **Build interval trees** (optional, with -t flag) for O(log n) query time
4. **For each A interval** (streaming):
   - Find candidate B intervals on same chromosome
   - Use interval tree or linear search to find overlaps
   - Check for overlap considering options
   - Output according to mode (-wa, -wb, -c, -d, -k, default)

With interval tree enabled (-t), query complexity is O(log n + k) where k is the number of results.
Without interval tree, query complexity is O(n) where n is the number of B intervals per chromosome.

## Performance

Benchmarked on 100,000 intervals in each file:

- **Time**: ~0.3 seconds
- **Memory**: ~50 MB
- **Speedup**: Comparable to bedtools intersect

Performance is similar to bedtools intersect for most use cases.

## Use Cases

### Find genes overlapping peaks

```bash
bedintersect -a genes.bed -b peaks.bed -wa > genes_with_peaks.bed
```

### Count binding sites per gene

```bash
bedintersect -a genes.bed -b binding_sites.bed -c > gene_binding_counts.bed
```

### Find genes without enhancers

```bash
bedintersect -a genes.bed -b enhancers.bed -v > genes_no_enhancers.bed
```

### Find significant overlaps

```bash
bedintersect -a regions.bed -b features.bed -m 100 -f 0.5 > significant.bed
```

### Strand-specific analysis

```bash
bedintersect -a genes.bed -b reads.bed -s -c > sense_coverage.bed
```

### Find reciprocal best hits

```bash
bedintersect -a genes.bed -b orthologs.bed -f 0.8 -F 0.8 -r > reciprocal_hits.bed
```

### Find nearest regulatory elements

```bash
bedintersect -a genes.bed -b enhancers.bed -d > gene_enhancer_distances.bed
```

### Get closest transcription factor binding site

```bash
bedintersect -a genes.bed -b tfbs.bed -k > closest_tfbs.bed
```

### Process very large feature databases

```bash
bedintersect -a queries.bed -b large_database.bed -t > results.bed
```

## Comparison with bedtools

### Similarities

- Same basic intersection algorithm
- Compatible input/output formats
- Supports most common options

### Differences

| Feature | bedintersect | bedtools intersect |
|---------|--------------|-------------------|
| Language | Go | C++ |
| Installation | Single binary | External dependency |
| Memory usage | Lower | Higher |
| Speed | Comparable | Comparable |
| Built-in gzip | Yes | No |
| Interval tree | Yes (optional) | No |
| Output modes (`-wa/-wb/-wo/-wao/-loj/-u/-c/-C/-v`) | Yes | Yes |
| Fractions (`-f/-F/-r/-e`) | Yes | Yes |
| Strand (`-s/-S`) | Yes | Yes |
| `-split` block-aware overlap | Yes | Yes |
| Multiple `-b` (`-names/-filenames/-sortout`) | Yes | Yes |
| `-sorted`/`-g` order validation | Yes | Yes |
| BED/GFF/VCF/BAM/CRAM input | Yes | Yes |
| BAM output (BAM/CRAM query, default) | Yes | Yes |
| CRAM binary output (CRAM query + reference) | Yes | Yes |
| Distance mode (`-d`) | Yes (extension) | No |
| Closest feature (`-k`) | Yes (extension) | Via separate tool |

The output of every shared option is validated byte-for-byte against the live
upstream `bedtools intersect` binary over the upstream test fixtures (see the
`*_parity_test.go` files), so it is a drop-in replacement for text (BED/GFF/VCF)
workflows.

**Use bedtools intersect when:**

- You're already using the bedtools suite end-to-end. BAM output (from a BAM or
  CRAM query) and CRAM output (from a CRAM query with a reference) are both
  supported here and match upstream (see BAM/CRAM output below).

## Testing

Run unit tests:

```bash
cd tools/bedintersect
go test ./pkg/bedintersect
```

Run with coverage:

```bash
go test -cover ./pkg/bedintersect
```

## Implementation Details

- Written in pure Go using standard library
- Uses existing `pkg/htsgo/bed` parser
- Chromosome-based indexing for efficiency
- Linear search within chromosome (fast for typical datasets)

## BAM/CRAM output

When the query (`-a`/`-abam`/`-ibam`) is a BAM or CRAM file and `-bed` is **not**
given, the surviving alignments are written back out as binary (the original
header plus the original alignment records, in input order), matching upstream
`bedtools intersect`'s default behaviour. The output framing follows upstream's
gating exactly:

- A **BAM** query writes **BAM**.
- A **CRAM** query writes **CRAM** only when a CRAM reference is available — via
  the `--cram-ref FILE` flag or the `CRAM_REFERENCE` environment variable.
  Without a reference a CRAM query writes **BAM**, exactly like upstream (whose
  writer opens htslib mode `wc` when a reference is set and `wb` otherwise). The
  `-ubam` flag does not change this — upstream selects the format from the
  reference alone — so a CRAM query with a reference stays CRAM under `-ubam`.

The BAM output is decoded and validated byte-for-byte (as SAM) against the live
upstream binary over the upstream BAM fixtures
(`cmd/bedintersect/bam_output_parity_test.go`). CRAM bytes are not identical
across encoders (block layout and codec choices differ), so the CRAM output is
validated by decoding **both** upstream's and our CRAM back to SAM (with the
fixture reference) and diffing the alignment records
(RNAME/POS/CIGAR/FLAG/MAPQ/SEQ/QUAL/aux) byte-for-byte; the header is compared as
a line set — see `cmd/bedintersect/cram_output_parity_test.go`.

The alignment-level flags behave exactly as upstream gates them for a BAM query
without `-bed`:

- **Produce BAM:** default, `-wa`, `-u` (each A alignment with ≥1 overlap, once)
  and `-v` (each A alignment with no overlap). `-C` is *not* an error: it falls
  through to the default BAM-output selection (upstream keeps it "printable").
  Unmapped reads never overlap, so they are absent under the default mode and
  reported under `-v`, exactly as upstream's `printUnmapped` path.
- **Error (require `-bed`):** `-c` (`writeCount`) and `-wo`/`-wao`
  (`writeOverlap`/`writeAllOverlap`) print the same `***** ERROR: … is not valid
  with BAM query input, unless bed output is specified with -bed option. *****`
  banner and exit 1.
- **Warn and ignore:** `-wb`/`-loj` and `-header` print the same stderr warning
  upstream does and then proceed as default BAM output (the flags have no BED
  columns to add in BAM mode).

`-bed` forces BED12 text output for a BAM/CRAM query instead (matching upstream's
`-bed` BED output byte-for-byte). The output BAM is BGZF-compressed; `-ubam`
(uncompressed BAM) is accepted but always emits compressed BAM (the `sam` writer
has no uncompressed-BAM mode). CRAM output is reference-free CRAM v3.0 (a
self-contained file that decodes without an external reference); its SEQ/QUAL is
carried verbatim from the records read out of the query, so a decode of our CRAM
matches a decode of upstream's record-for-record.

## Limitations

- **Uncompressed BAM (`-ubam`).** Accepted but always emits BGZF-compressed BAM;
  the `pkg/htsgo/sam` writer has no uncompressed-BAM mode. Upstream's `-ubam`
  compression hook is itself a no-op, and the decoded records are identical
  either way, so this does not affect record-level parity.
- Loads the B file(s) completely into memory (necessary for random access).
- The chromosome naming-convention warning is emitted before the data rather
  than interleaved into it; when stdout and stderr are captured separately the
  content matches upstream exactly, but a combined `2>&1` stream orders the
  warning differently.

## Recent Enhancements

- ✅ Interval tree for very large B files (use -t flag)
- ✅ Reciprocal overlap mode (use -r flag with -f and -F)
- ✅ Distance to nearest feature (use -d flag)
- ✅ Output closest feature (use -k flag)
- ✅ Streaming mode for file A (always enabled)
- ✅ Left outer join (`-loj`), write-overlap (`-wo`/`-wao`), `-wa -wb`
  side-by-side, and `-split` block-aware overlap — all validated
  byte-for-byte against the upstream `bedtools intersect` binary (BED3–BED12
  null shapes, zero-length intervals, B-file order, and `-s` UNKNOWN-strand
  handling).
- ✅ Full upstream flag parity: `-u` (unique), `-C` (per-B-file counts),
  `-S` (opposite strand), `-e` (either-fraction), multiple `-b` databases with
  `-names`/`-filenames`/`-sortout`, `-header`, `-sorted`/`-g` order validation,
  and `-nonamecheck` — each validated byte-for-byte against the live upstream
  binary over the upstream fixtures.
- ✅ Bin-order output: default (non-sorted) hit ordering reproduces upstream's
  UCSC-bin traversal order, so nested/overlapping B records print in the same
  order as `bedtools intersect`.
- ✅ Input-format parity: BED/GFF/VCF (including structural-variant END/SVLEN
  spans), BAM/CRAM (with `/1`,`/2` mate suffixes and unmapped-read placeholders),
  and gzip/BGZF-compressed text on stdin.
- ✅ BAM binary output: a BAM/CRAM query without `-bed` writes the intersecting
  alignments back out as BAM (original header + records, in input order),
  matching upstream's default and its `-u`/`-v`/`-wa` selection, its
  `-c`/`-wo`/`-wao` "requires `-bed`" errors, and its `-wb`/`-loj`/`-header`
  warn-and-ignore behaviour — validated byte-for-byte (decoded to SAM) against
  the live upstream binary.
- ✅ CRAM binary output: a CRAM query without `-bed` and **with a CRAM reference**
  (`--cram-ref` / `CRAM_REFERENCE`) writes the intersecting alignments back out
  as CRAM, matching upstream's reference-gated format selection (no reference →
  BAM) and its `-u`/`-v`/`-wa`/`-C` selection — validated against the live
  upstream binary by decoding both CRAM outputs to SAM and diffing the records.

These join/overlap modes echo the original A and B input columns verbatim,
in the original B-file order, matching upstream. BAM/VCF/GFF inputs and the
`bedclosest` directional flags (`-id`/`-iu`/`-fu`/`-fd`) remain out of scope —
see `docs/PARITY_ROADMAP.md` for the documented remainder.

## Future Enhancements

- [ ] Sorted file optimization to reduce memory usage for B file
- [ ] Parallel processing for multi-core systems

## Contributing

Contributions welcome! Please:

- Add tests for new features
- Follow Go coding standards
- Update documentation

## License

Apache License 2.0 - See LICENSE file for details.

## See Also

- [bedtools](https://bedtools.readthedocs.io/) - Comprehensive genomic interval toolkit
- [BED format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) - Format specification
- Other tools: bedmerge, seqtk, prinseq, sickle, skewer, fastp
