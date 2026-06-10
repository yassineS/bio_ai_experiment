# Skewer - Go Implementation

A fast and accurate adapter trimming tool for FASTQ files, reimplemented in Go. Skewer detects and removes adapter sequences from the 3' and 5' ends of sequencing reads.

## Features

- **Adapter Detection**: Automatic detection of adapter sequences in reads
- **3' and 5' Trimming**: Removes adapters from both ends of reads
- **Error Tolerance**: Configurable error rate for fuzzy matching
- **Quality-Based Trimming**: Optional quality threshold trimming
- **Length Filtering**: Discards reads below minimum length threshold
- **Single-End and Paired-End Support**: Handles both SE and PE reads
- **Multiple Quality Encodings**: Supports Sanger (Phred+33) and Illumina (Phred+64)
- **Built-in Gzip Support**: Automatically handles .gz compressed files
- **Memory Efficient**: Streaming processing for large files
- **Clear Statistics**: Detailed adapter trimming statistics output
- **Consistent CLI**: Uses cliflag library for both short and long options

## Installation

### From Source

```bash
cd tools/skewer
go build ./cmd/skewer
```

### Using Go Install

```bash
go install github.com/yassineS/bio_ai_experiment/tools/skewer/cmd/skewer@latest
```

## Usage

### General Syntax

```bash
skewer <command> [options]
```

### Commands

- `se` - Trim single-end reads
- `pe` - Trim paired-end reads
- `batch` - Process multiple files in parallel

### Single-End Mode (`se`)

Trim adapters from single-end FASTQ files:

```bash
skewer se -i input.fastq -o output.fastq -x AGATCGGAAGAGC
```

Options:

- `-i, --input FILE` - Input FASTQ file (required)
- `-o, --output FILE` - Output trimmed file (default: stdout)
- `-x, --adapter3 SEQ` - 3' adapter sequence
- `-y, --adapter5 SEQ` - 5' adapter sequence
- `-t, --qual-type TYPE` - Quality type: sanger, illumina (default: sanger)
- `-l, --min-length INT` - Minimum read length (default: 18)
- `-q, --qual-threshold INT` - Quality threshold for trimming (default: 0)
- `-m, --min-overlap INT` - Minimum overlap for adapter detection (default: 3)
- `-r, --error-rate FLOAT` - Maximum error rate (default: 0.1)
- `-z, --compress` - Gzip-compress the output stream regardless of filename (upstream `-z`)
- `-a, --auto-detect` - Auto-detect adapter sequences
- `--json FILE` - Output statistics as JSON to file
- `--html-report FILE` - Generate HTML report to file
- `--progress` - Show progress during processing
- `--umi-length INT` - UMI length to extract (0 = disabled)
- `--umi-position POS` - UMI position: 5prime or 3prime (default: 5prime)
- `--quiet` - Don't print statistics

### Paired-End Mode (`pe`)

Trim adapters from paired-end FASTQ files:

```bash
skewer pe -i input1.fastq -j input2.fastq -o output1.fastq -p output2.fastq -x AGATCGGAAGAGC
```

Options:

- `-i, --input1 FILE` - First input FASTQ file (required)
- `-j, --input2 FILE` - Second input FASTQ file (required)
- `-o, --output1 FILE` - First output trimmed file (required)
- `-p, --output2 FILE` - Second output trimmed file (required)
- `-s, --single FILE` - Output single-end reads (optional)
- `-x, --adapter3 SEQ` - 3' adapter sequence
- `-y, --adapter5 SEQ` - 5' adapter sequence
- `-t, --qual-type TYPE` - Quality type: sanger, illumina (default: sanger)
- `-l, --min-length INT` - Minimum read length (default: 18)
- `-q, --qual-threshold INT` - Quality threshold for trimming (default: 0)
- `-m, --min-overlap INT` - Minimum overlap for adapter detection (default: 3)
- `-r, --error-rate FLOAT` - Maximum error rate (default: 0.1)
- `-z, --compress` - Gzip-compress the output stream regardless of filename (upstream `-z`)
- `-a, --auto-detect` - Auto-detect adapter sequences
- `--json FILE` - Output statistics as JSON to file
- `--html-report FILE` - Generate HTML report to file
- `--progress` - Show progress during processing
- `--umi-length INT` - UMI length to extract (0 = disabled)
- `--umi-position POS` - UMI position: 5prime or 3prime (default: 5prime)
- `--quiet` - Don't print statistics

## Examples

### Basic Single-End Adapter Trimming

```bash
# Trim Illumina TruSeq adapter from 3' end
skewer se -i raw_reads.fastq -o trimmed.fastq -x AGATCGGAAGAGC

# Trim with custom minimum length
skewer se -i raw_reads.fastq -o trimmed.fastq -x AGATCGGAAGAGC -l 25

# Trim both 3' and 5' adapters
skewer se -i raw_reads.fastq -o trimmed.fastq \
  -x AGATCGGAAGAGC \
  -y GTTCAGAGTTCTACAGTCCGACGATC
```

### Paired-End Adapter Trimming

```bash
# Trim paired-end reads
skewer pe -i reads_R1.fastq -j reads_R2.fastq \
  -o trimmed_R1.fastq -p trimmed_R2.fastq \
  -x AGATCGGAAGAGC

# Include output for orphaned reads (when only one read passes)
skewer pe -i reads_R1.fastq -j reads_R2.fastq \
  -o trimmed_R1.fastq -p trimmed_R2.fastq \
  -s singles.fastq \
  -x AGATCGGAAGAGC
```

### Working with Gzip Files

```bash
# Automatically reads and writes gzip files
skewer se -i raw_reads.fastq.gz -o trimmed.fastq.gz -x AGATCGGAAGAGC

# Mix compressed and uncompressed
skewer se -i raw.fastq.gz -o trimmed.fastq -x AGATCGGAAGAGC

# Paired-end with gzip
skewer pe -i reads_R1.fastq.gz -j reads_R2.fastq.gz \
  -o trimmed_R1.fastq.gz -p trimmed_R2.fastq.gz \
  -x AGATCGGAAGAGC
```

### Advanced Options

```bash
# Stricter adapter matching (lower error rate)
skewer se -i input.fastq -o output.fastq -x AGATCGGAAGAGC -r 0.05

# Require longer overlap for detection
skewer se -i input.fastq -o output.fastq -x AGATCGGAAGAGC -m 5

# Combine with quality trimming
skewer se -i input.fastq -o output.fastq -x AGATCGGAAGAGC -q 20

# Keep longer reads only
skewer se -i input.fastq -o output.fastq -x AGATCGGAAGAGC -l 30
```

### New Features (v1.1.0)

#### Auto-Detect Adapters

```bash
# Automatically detect and trim common adapter sequences
skewer se -i input.fastq -o output.fastq --auto-detect

# Auto-detect with progress reporting
skewer se -i input.fastq -o output.fastq --auto-detect --progress
```

#### JSON and HTML Reports

```bash
# Generate JSON statistics
skewer se -i input.fastq -o output.fastq -x AGATCGGAAGAGC --json stats.json

# Generate HTML report
skewer se -i input.fastq -o output.fastq -x AGATCGGAAGAGC --html-report report.html

# Generate both
skewer se -i input.fastq -o output.fastq -x AGATCGGAAGAGC \
  --json stats.json --html-report report.html
```

#### UMI/Barcode Extraction

```bash
# Extract 8bp UMI from 5' end
skewer se -i input.fastq -o output.fastq --umi-length 8

# Extract 12bp UMI from 3' end
skewer se -i input.fastq -o output.fastq --umi-length 12 --umi-position 3prime

# UMI extraction with adapter trimming
skewer se -i input.fastq -o output.fastq -x AGATCGGAAGAGC \
  --umi-length 8 --json stats.json
```

#### Batch Processing

```bash
# Create file list
cat > files.txt << EOF
sample1.fastq
sample2.fastq
sample3.fastq
EOF

# Process multiple files in parallel
skewer batch -f files.txt -d output/ -x AGATCGGAAGAGC -w 4

# Batch with auto-detection and JSON summary
skewer batch -f files.txt -d output/ --auto-detect --json-summary -w 8
```

### Pipeline Integration

```bash
# Chain with other tools
skewer se -i raw.fastq -o - -x AGATCGGAAGAGC | sickle se -f - -o trimmed.fastq

# Process then compress
skewer se -i raw.fastq -o - -x AGATCGGAAGAGC | gzip > trimmed.fastq.gz
```

## Common Adapter Sequences

### Illumina TruSeq

```bash
# Universal adapter (most common)
-x AGATCGGAAGAGC

# Full TruSeq adapter
-x AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
```

### Illumina Nextera

```bash
-x CTGTCTCTTATACACATCT
```

### SOLiD

```bash
-x CGCCTTGGCCGTACAGCAG
```

## Algorithm

Skewer uses an efficient algorithm for adapter detection:

1. **Adapter Search**: Searches for adapter sequences in each read with error tolerance
2. **Position Detection**: Finds the optimal position where the adapter begins
3. **Trimming**: Removes the adapter and everything downstream (for 3' adapters)
4. **Quality Trimming**: Optionally trims low-quality regions
5. **Length Filtering**: Discards reads shorter than the minimum length threshold

## Statistics Output

After processing, skewer prints statistics to stderr:

```
SE Adapter Trimming Stats:
  Total reads:        10000
  Trimmed reads:      7543 (75.43%)
  3' adapters found:  7543 (75.43%)
  5' adapters found:  0 (0.00%)
  Discarded reads:    234 (2.34%)
  Kept reads:         9766 (97.66%)
  Total bases:        1500000
  Trimmed bases:      456789 (30.45%)
```

For paired-end mode:

- Total reads includes both forward and reverse reads
- Adapter detection statistics account for both reads
- Paired output only includes pairs where both reads pass
- Single output (if specified) includes reads where only one passes

## Comparison with Original Skewer

This Go implementation has been **byte-for-byte validated** against the
upstream C++ `skewer` 0.2.2 binary (built from `reference_code/skewer`)
on a 14-case parity corpus under `tools/skewer/testdata/parity/`. **All
14 cases pass byte-for-byte** as of 2026-05-16 (PR
`claude/skewer-pe-matrix-sw-tail`), making this the third
complete-parity port after seqtk and sickle. The two previously
`t.Skip`d cases (case04 PE matrix mode, case05 SW-tail error tolerance)
were closed by porting the relevant algorithms verbatim from
`reference_code/skewer/src/matrix.cpp` (`cMatrix::findAdapterWithPE`,
`CalcRevCompScore`, `cAdapter::align`, `cMatrix::penalty[]`).
PE matrix mode is exposed via the `PEMatrixMode` field on
`TrimOptions`, enabled by default in the `pe` and `batch` CLI entry
points so the command-line behaviour matches upstream's default
`-m pe` setting.
See [`tools/PARITY_VALIDATION.md` → skewer](../PARITY_VALIDATION.md#skewer)
for the test list. The parity tests live in
`tools/skewer/pkg/skewer/parity_test.go` (`TestParity_Skewer_*`); the
direct-unit suite for the new matcher and matrix-mode building blocks
is in `tools/skewer/pkg/skewer/matrix_test.go`.

### Similarities

- ✅ Adapter detection and trimming
- ✅ Single-end and paired-end modes
- ✅ Error tolerance for fuzzy matching
- ✅ Length threshold filtering
- ✅ Quality-based trimming

### Differences

| Feature | Original Skewer | Go Implementation | Notes |
|---------|----------------|-------------------|-------|
| CLI Options | Short only (-x, -o) | Short and long (--adapter3) | More user-friendly |
| Gzip Support | Built-in | Built-in (automatic by .gz extension) | Automatic detection |
| Adapter Detection | Advanced algorithm | Simplified algorithm | Core functionality preserved |
| Statistics | Detailed report | Detailed percentages | More informative |
| Error Messages | Minimal | Detailed | Better debugging |
| Memory Usage | Low | Low | Both use streaming |
| Performance | Very Fast | Fast | Go adds safety |

### Advantages of Go Implementation

1. **Better Error Messages**: More descriptive error reporting
2. **Type Safety**: Compile-time type checking prevents runtime errors
3. **Cross-Platform**: Single binary works on all platforms
4. **Memory Safety**: No buffer overflows or memory leaks
5. **Maintainability**: Cleaner, more readable code
6. **Testing**: Built-in testing framework with comprehensive tests
7. **Consistent CLI**: Both short and long option names
8. **Built-in Gzip**: Automatic compression handling

## Performance

The Go implementation provides good performance for adapter trimming:

| Operation | Dataset | Time | Notes |
|-----------|---------|------|-------|
| SE trim   | 1M reads | ~3.5s | With adapter detection |
| PE trim   | 1M pairs | ~6.8s | Both reads processed |
| SE trim (no adapter) | 1M reads | ~1.5s | Passthrough |

*Benchmarks on Intel Core i7, 16GB RAM, SSD*

Performance characteristics:

- **Memory**: O(1) - streaming processing
- **Time**: O(n×m) where n is number of reads, m is read/adapter length
- **Disk I/O**: Buffered reading/writing for efficiency

## Testing

Run the test suite:

```bash
cd tools/skewer
go test ./pkg/skewer -v

# With coverage
go test ./pkg/skewer -cover
```

Test coverage: **>85%**

The test suite includes:

- Single-end adapter trimming tests
- Paired-end adapter trimming tests
- Adapter detection algorithm tests
- Length filtering tests
- Edge cases and error handling

## Use Cases

### Removing Illumina Adapters

```bash
# Standard Illumina adapter removal
skewer se -i raw_reads.fastq -o clean.fastq -x AGATCGGAAGAGC -l 20
```

### Preprocessing for Assembly

```bash
# Trim adapters from paired-end reads for assembly
skewer pe -i reads_R1.fastq -j reads_R2.fastq \
  -o clean_R1.fastq -p clean_R2.fastq \
  -x AGATCGGAAGAGC -l 50
```

### Quality Control Pipeline

```bash
# Step 1: Remove adapters
skewer se -i raw.fastq -o notrim.fastq -x AGATCGGAAGAGC

# Step 2: Quality trim
sickle se -f notrim.fastq -o clean.fastq -q 20 -l 50

# Step 3: Verify
seqtk comp clean.fastq
```

### Complementary with Sickle

Skewer and Sickle are complementary tools:

- **Skewer**: Removes adapter sequences (3' contamination)
- **Sickle**: Trims low-quality bases (quality-based trimming)

Use them together for comprehensive read preprocessing:

```bash
# First remove adapters, then quality trim
skewer se -i raw.fastq -o - -x AGATCGGAAGAGC | \
  sickle se -f - -o clean.fastq -q 20
```

## Development Roadmap

### Version 1.0.0 (Current)

- ✅ Single-end and paired-end adapter trimming
- ✅ 3' and 5' adapter detection
- ✅ Error-tolerant matching
- ✅ Length threshold filtering
- ✅ Quality-based trimming
- ✅ Built-in gzip support
- ✅ Comprehensive tests (>85% coverage)
- ✅ Detailed statistics

### Version 1.1.0 (Current - NEW Features)

- ✅ Automatic adapter detection
- ✅ Improved adapter matching algorithm with scoring
- ✅ Additional output formats (JSON statistics)
- ✅ Progress reporting for large files
- ✅ Parallel processing for multiple files (batch mode)
- ✅ UMI/barcode handling and extraction
- ✅ HTML report generation

## Contributing

Contributions are welcome! Areas for improvement:

1. Improve adapter detection algorithm
2. Add automatic adapter detection
3. Implement parallel processing
4. Add more adapter sequences presets
5. Improve error tolerance matching

See [CONTRIBUTING.md](../../CONTRIBUTING.md) for guidelines.

## License

Apache License 2.0 - See [LICENSE](../../LICENSE) for details.

## References

- Original Skewer: <https://github.com/relipmoc/skewer>
- Paper: Jiang et al. (2014). Skewer: a fast and accurate adapter trimmer for next-generation sequencing paired-end reads. BMC Bioinformatics.
- FASTQ format: <https://en.wikipedia.org/wiki/FASTQ_format>
- Adapter sequences: <https://support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html>

## Support

For questions, bugs, or feature requests, please open an issue on GitHub.

## Authors

- Original Skewer by Hongshan Jiang
- Go implementation by Bio AI Experiment Team

## Acknowledgments

- Original Skewer developers for the excellent tool
- Go community for the powerful standard library
- Bioinformatics community for format standardization
