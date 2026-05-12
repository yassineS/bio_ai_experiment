# Bio AI Experiment

Deploying agents to recode and modernise bioinformatics tools.

## Overview

This repository is an experimental project that leverages AI agents and automated tools to improve the bioinformatics software ecosystem. The primary focus is to identify, analyze, recode, test, and document popular bioinformatics tools, making them more robust, performant, and easier to use.

## Project Goals

The main objectives of this repository are to use AI agents to:

1. **Build a Comprehensive Tool Database**
   - Compile a list of the top 100 most used bioinformatics tools
   - Document each tool's name, link, citation, and download statistics
   - Create a structured database for easy reference

2. **Identify Areas for Improvement**
   - Analyze tools for poor code quality, performance issues, and documentation gaps
   - Identify edge cases and use cases that are not well-handled
   - Document findings systematically

3. **Recode Tools in GoLang**
   - Reimplement selected tools using Go for improved performance and maintainability
   - Maintain compatibility with original tool functionality
   - Leverage Go's concurrency and performance features

4. **Comprehensive Testing**
   - Provide test data for each tool
   - Write unit tests for main functionalities
   - Write unit tests for edge cases
   - Ensure robust test coverage

5. **Documentation**
   - Document all code comprehensively
   - Create user guides and API documentation
   - Document known issues and workarounds
   - Maintain up-to-date documentation throughout iterations

6. **Iterative Improvement**
   - Repeat the analysis-recode-test-document cycle until tools are robust
   - Continuously refine based on findings and feedback

7. **MCP Integration**
   - Write Model Context Protocol (MCP) servers for each tool
   - Simplify integration of bioinformatics tools with Large Language Models
   - Enable easier access to bioinformatics capabilities through AI interfaces

## Ultimate Goals

- **Improve Usability**: Make bioinformatics tools more accessible and easier to use
- **Enhance Documentation**: Provide clear, comprehensive documentation for all tools
- **Boost Performance**: Leverage modern programming practices and languages for better performance
- **Document AI Agent Utility**: Track and document the effectiveness (or lack thereof) of coding agents in this process

## Repository Structure

```
bio_ai_experiment/
├── .github/
│   └── agents/          # Agent configuration files
├── tools/               # Directory for recoded tools
│   └── [tool-name]/
│       ├── src/         # Go source code
│       ├── tests/       # Unit tests
│       ├── testdata/    # Test data
│       └── docs/        # Tool-specific documentation
├── analysis/            # Tool analysis and findings
├── mcp-servers/         # MCP server implementations
└── docs/                # General documentation

```

## Getting Started

### Prerequisites

- Go 1.21 or later
- Git
- Basic understanding of bioinformatics tools

### Installation

```bash
git clone https://github.com/yassineS/bio_ai_experiment.git
cd bio_ai_experiment
```

### Building Tools

Individual tools can be built from their respective directories:

```bash
cd tools/[tool-name]
go build ./cmd/[tool-name]
```

Example:

```bash
cd tools/seqtk
go build ./cmd/seqtk
./seqtk help
```

### Running Tests

```bash
# Run all tests
go test ./...

# Run tests with coverage
go test -cover ./...

# Run tests for a specific tool
cd tools/seqtk
go test ./pkg/seqtk
```

### Using Shared Libraries

The bioformats library provides reusable parsers for common file formats:

```bash
# View library documentation
go doc github.com/yassineS/bio_ai_experiment/pkg/bioformats/fasta
```

## Current Status

### ✅ Completed

**Shared Libraries:**

- FASTA format parser/writer with validation and utilities
- FASTQ format parser/writer with Phred33/64 support
- VCF format parser/writer with genotype methods
- BED format parser/writer with interval operations
- **cliflag** library for consistent CLI flag handling (short and long options)

**Tools:**

- **seqtk v1.0.0** - Complete reimplementation
  - comp, fq2fa, seq -r, sample, trimfq commands
  - 85.7% test coverage, all tests passing
  - Performance comparable to original C version
  - Zero external dependencies
  - **New**: Consistent CLI with both short and long option flags

- **prinseq v1.0.0** - Core functionality implemented
  - Statistics calculation (reads, bases, lengths, GC, N content, quality)
  - Multi-criteria filtering (length, GC, N content, quality)
  - 90.2% test coverage, all tests passing
  - 20-26% faster than original Perl version
  - Zero external dependencies
  - Consistent CLI with both short and long option flags

**Documentation:**

- Comprehensive Go implementation guide
- Bioformats library documentation
- Tool analysis and comparisons (seqtk, prinseq)
- Best practices and patterns

### 📋 In Progress

- Additional seqtk commands
- PRINSEQ trimming and duplicate removal
- Compressed file support (gzip)
- Performance benchmarking suite

See [docs/GO_IMPLEMENTATION_SUMMARY.md](docs/GO_IMPLEMENTATION_SUMMARY.md) for detailed status.

## Available Tools

### seqtk - FASTA/Q Sequence Processor

Fast and efficient sequence processing tool with consistent CLI flags.

```bash
cd tools/seqtk
go build ./cmd/seqtk

# Get statistics
./seqtk comp sequences.fasta

# Convert FASTQ to FASTA (short options)
./seqtk fq2fa reads.fastq > reads.fasta

# Convert FASTQ to FASTA (long options)
./seqtk fq2fa --output reads.fasta reads.fastq

# Reverse complement
./seqtk seq -r sequences.fasta > rev_comp.fasta
./seqtk seq --reverse sequences.fasta > rev_comp.fasta

# Sample 10% of reads
./seqtk sample reads.fastq 0.1 > sample.fastq

# Trim low-quality bases
./seqtk trimfq -q 20 reads.fastq > trimmed.fastq
./seqtk trimfq --quality 20 reads.fastq > trimmed.fastq
```

See [tools/seqtk/README.md](tools/seqtk/README.md) for complete documentation.

### prinseq - Sequence Quality Control

Sequence quality control and preprocessing tool.

```bash
cd tools/prinseq
go build ./cmd/prinseq

# Get statistics
./prinseq stats -fastq reads.fastq

# Filter by length
./prinseq filter -fastq reads.fastq -min_len 100 > filtered.fastq

# Filter by GC content
./prinseq filter -fastq reads.fastq -min_gc 40 -max_gc 60 > filtered.fastq

# Filter by quality
./prinseq filter -fastq reads.fastq -min_qual_mean 20 > filtered.fastq

# Combined filters
./prinseq filter -fastq reads.fastq \
  -min_len 100 \
  -min_gc 40 \
  -max_gc 60 \
  -min_qual_mean 20 \
  -ns_max_p 5 \
  -out_good filtered.fastq
```

See [tools/prinseq/README.md](tools/prinseq/README.md) for complete documentation.

## Documentation

- [Go Implementation Guide](docs/GOLANG_GUIDE.md) - Best practices and patterns
- [Bioformats Library](pkg/bioformats/README.md) - Format parser documentation
- [Implementation Summary](docs/GO_IMPLEMENTATION_SUMMARY.md) - Project overview and metrics
- [Tool Analyses](docs/tools/) - Detailed tool comparisons

## Project Structure

```
bio_ai_experiment/
├── pkg/bioformats/      # Shared format libraries
│   ├── fasta/          # FASTA parser/writer
│   ├── fastq/          # FASTQ parser/writer
│   ├── vcf/            # VCF parser/writer
│   └── bed/            # BED parser/writer
├── tools/              # Go tool implementations
│   └── seqtk/         # First reference implementation
├── docs/               # Documentation
│   ├── GOLANG_GUIDE.md
│   ├── GO_IMPLEMENTATION_SUMMARY.md
│   └── tools/         # Tool analyses
├── analysis/           # Analysis scripts and data
├── reference_code/     # Original tool repositories
└── mcp-servers/        # MCP server implementations
```

## Contributing

Contributions are welcome! This project is designed to be collaborative and benefits from diverse perspectives.

### How to Contribute

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add some amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

### Contribution Guidelines

- Follow Go best practices and conventions
- Include tests for new functionality
- Update documentation as needed
- Ensure all tests pass before submitting PR

See [CONTRIBUTING.md](CONTRIBUTING.md) for detailed contribution guidelines.

### Getting Help

- 💬 [Discussions](https://github.com/yassineS/bio_ai_experiment/discussions) - Ask questions, share ideas
- 🐛 [Issues](https://github.com/yassineS/bio_ai_experiment/issues) - Report bugs, request features
- 📖 [Documentation](docs/) - Guides and references
- 📋 [Analysis](analysis/) - Tool analyses and evaluations

## License

This project is licensed under the Apache License 2.0 - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- The bioinformatics community for creating the original tools
- AI/LLM technologies that make this experiment possible
- All contributors to this project

## Project Status

This is an active experimental project. Progress and findings will be documented regularly.

## Contact

For questions, suggestions, or collaboration opportunities, please open an issue in this repository.
