# CLI Design Conventions for Bio AI Experiment Tools

This document outlines the command-line interface (CLI) design conventions used across all tools in the Bio AI Experiment repository. These conventions ensure consistency and usability across all reimplemented bioinformatics tools.

## General Principles

1. **POSIX-style flags**: Use both short (single dash, single letter) and long (double dash, full word) options
2. **Consistency**: Use the same option names for similar functionality across all tools
3. **Discoverability**: Provide helpful error messages and comprehensive help text
4. **Backward compatibility**: Where possible, maintain compatibility with original tool options

## Option Naming Conventions

### Short Options (Single Letter)

Short options use a single dash followed by a single letter:

- `-i` for input files
- `-o` for output files
- `-q` for quality-related options
- `-l` for length-related options
- `-g` for GC content options
- `-t` for trimming options
- `-d` for duplicate removal options
- `-h` for help
- `-v` for version

### Long Options (Full Word)

Long options use a double dash followed by a descriptive word or hyphenated phrase:

- `--input` for input files
- `--output` for output files
- `--quality` for quality thresholds
- `--min-length` for minimum length
- `--max-length` for maximum length
- `--help` for help
- `--version` for version

### Paired Options

When tools support paired-end reads, use numbered suffixes:

- `-i1, --input1` for first read file
- `-i2, --input2` for second read file
- `-o1, --output1` for first output file
- `-o2, --output2` for second output file

## Standard Options

All tools should support these standard options:

- `-h, --help`: Display help information
- `-v, --version`: Display version information
- `-i, --input`: Primary input file (use `-` for stdin)
- `-o, --output`: Primary output file (use `-` for stdout if not specified)

## File Format Options

For tools that support multiple formats:

- `--fasta`: Input is FASTA format
- `--fastq`: Input is FASTQ format
- `--format`: Specify format explicitly (fasta, fastq, etc.)

## Common Filtering Options

Standard options across filtering tools:

- `-l, --min-length INT`: Minimum sequence length
- `-L, --max-length INT`: Maximum sequence length
- `-g, --min-gc FLOAT`: Minimum GC content percentage
- `-G, --max-gc FLOAT`: Maximum GC content percentage
- `-q, --min-quality FLOAT`: Minimum quality score
- `-Q, --max-quality FLOAT`: Maximum quality score
- `-n, --max-ns INT`: Maximum number of N bases
- `-N, --max-ns-percent FLOAT`: Maximum percentage of N bases

## Trimming Options

Standard trimming options:

- `--trim-left INT`: Trim N bases from 5' end
- `--trim-right INT`: Trim N bases from 3' end
- `--trim-qual-left INT`: Quality-based trimming from 5' end
- `--trim-qual-right INT`: Quality-based trimming from 3' end
- `--trim-n-left INT`: Trim poly-N from 5' end
- `--trim-n-right INT`: Trim poly-N from 3' end

## Duplicate Removal Options

- `-d, --derep MODE`: Duplicate removal mode (1=exact, 4=revcomp, 5=both)
- `--derep-min INT`: Minimum occurrences to keep

## Output Options

- `-o, --output FILE`: Main output file
- `--out-good FILE`: Passing sequences (alternative naming)
- `--out-bad FILE`: Rejected sequences
- `--stats FILE`: Statistics output file

## Help Text Format

Help text should follow this structure:

```
tool-name - Brief description

Usage:
  tool-name <command> [options]

Commands:
  command1    Description
  command2    Description

Options:
  -s, --short-name TYPE    Description
                           Additional description lines indented

Examples:
  tool-name command1 -i input.fastq -o output.fastq
  tool-name command2 --min-length 100 input.fasta
```

## Examples

### Good Example (PRINSEQ filter)

```bash
# Using short options
prinseq filter -i reads.fastq -o filtered.fastq -l 100 -q 20

# Using long options
prinseq filter --input reads.fastq --output filtered.fastq --min-length 100 --min-quality 20

# Mixed short and long options
prinseq filter -i reads.fastq --min-length 100 --trim-qual-left 20
```

### Good Example (paired-end)

```bash
# Short options
prinseq filter -i1 R1.fastq -i2 R2.fastq -o1 out_R1.fastq -o2 out_R2.fastq

# Long options
prinseq filter --input1 R1.fastq --input2 R2.fastq --output1 out_R1.fastq --output2 out_R2.fastq
```

## Implementation Guidelines

### Go Implementation

For Go tools using the standard `flag` package, we implement a custom flag type that accepts both short and long options:

```go
// Define both short and long versions
fs.String("i", "", "")
fs.String("input", "", "Input file (use '-' for stdin)")

// Or use a helper function
addFlag(fs, "i", "input", "", "Input file (use '-' for stdin)")
```

### Validation

- Validate that required options are provided
- Provide clear error messages when options conflict
- Show relevant help text when errors occur

## Version History

- 2025-10-20: Initial version established for PRINSEQ tool
- Future: May be extended as more tools are ported

## References

- POSIX Utility Conventions: <https://pubs.opengroup.org/onlinepubs/9699919799/basedefs/V1_chap12.html>
- GNU Coding Standards: <https://www.gnu.org/prep/standards/html_node/Command_002dLine-Interfaces.html>
