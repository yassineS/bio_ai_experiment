# bedintersect - BED Interval Intersection Finder

A fast tool to find intersecting intervals between two BED files, implemented in Go.

## Features

- **Fast performance**: Efficient interval searching with chromosome indexing and optional interval tree
- **Flexible output**: Report intersections, original entries, counts, distance, or closest feature
- **Multiple filter options**: Minimum overlap, fraction overlap, reciprocal overlap, strand-specific
- **Advanced modes**: Distance to nearest feature, closest feature reporting
- **Interval tree**: Optional interval tree data structure for very large B files
- **Standard BED format**: Compatible with bedtools intersect
- **Built-in gzip support**: Automatic handling of .gz files
- **Statistics**: Optional intersection statistics output
- **Streaming**: File A is processed in streaming fashion (memory-efficient)

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

- `-a, --input-a FILE` - Input BED file A (required)
- `-b, --input-b FILE` - Input BED file B (required)
- `-o, --output FILE` - Output file (default: stdout)
- `-m, --min-overlap INT` - Minimum overlap required (default: 1)
- `-f, --fraction-a NUM` - Minimum fraction of A that must overlap (0.0-1.0)
- `-F, --fraction-b NUM` - Minimum fraction of B that must overlap (0.0-1.0)
- `-r, --reciprocal` - Require reciprocal overlap (both -f and -F must be satisfied)
- `-s, --strand` - Only report hits on same strand
- `-v, --invert` - Report A entries with NO overlap with B
- `-wa, --write-a` - Write original A entry (default: write intersection)
- `-wb, --write-b` - Write B entry instead of A (with `-wa`: write A and B
  side-by-side per overlap)
- `-loj` - Left outer join: report every A; B is null when there is no overlap
- `-wo` - Write A, B and the number of overlapping bases per overlap
- `-wao` - Like `-wo`, but also report A (with null B and overlap `0`) when there
  is no overlap
- `-split` - Treat split BED12 entries as distinct intervals (block-aware
  overlap and overlap-base counting)
- `-c, --count` - Report count of B overlaps for each A
- `-d, --distance` - Report distance to nearest B feature
- `-k, --closest` - Output closest B feature for each A
- `-t, --tree` - Use interval tree for large B files (better performance)
- `-S, --stats` - Print statistics to stderr
- `-h, --help` - Show help message

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
bedintersect -a genes.bed -b peaks.bed -S > overlaps.bed
```

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
| Reciprocal mode | Yes | Yes |
| Distance mode | Yes | No |
| Closest feature | Yes | Via separate tool |
| Advanced options | Growing | Extensive |
| Sorted requirement | No | No |

**Use bedintersect when:**

- You need a simple, standalone tool
- You want built-in gzip support
- You prefer Go tools
- You use common intersection operations
- You need distance or closest feature functionality
- You're working with very large B files (use -t flag)
- You need reciprocal overlap filtering

**Use bedtools intersect when:**

- You need advanced options (wao, sorted optimization, etc.)
- You're already using bedtools suite
- You need perfect output compatibility with existing pipelines

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

## Limitations

- Loads B file completely into memory (necessary for random access)
- No sorted file optimization for memory-constrained environments

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
