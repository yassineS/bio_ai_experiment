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
   - **POSIX-compliant CLI**: once a tool reaches full feature parity with its
     original, its command-line interface must also be POSIX-compliant — POSIX
     short options (`-i`, `-q`, ...), GNU-style long options (`--input`,
     `--quality`, ...), `--` to end option parsing, `-` to mean stdin/stdout,
     and clean exit codes. See [docs/CLI_CONVENTIONS.md](docs/CLI_CONVENTIONS.md).

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

7. **Robust, Drop-in POSIX CLIs** (the project's interface direction)
   - Expose every tool as a standalone, POSIX-compliant command-line program
   - Keep retro-compatibility with the upstream tool's flags so each port is a
     drop-in replacement (getopt-style bundling, short + long options) — see
     [docs/CLI_CONVENTIONS.md](docs/CLI_CONVENTIONS.md)
   - The earlier plan to ship per-tool Model Context Protocol (MCP) servers has
     been **descoped**; standalone CLIs are the supported interface (see
     [mcp-servers/README.md](mcp-servers/README.md))

## Ultimate Goals

- **Improve Usability**: Make bioinformatics tools more accessible and easier to use
- **Enhance Documentation**: Provide clear, comprehensive documentation for all tools
- **Boost Performance**: Leverage modern programming practices and languages for better performance
- **Feature Parity + POSIX-Compliant CLIs**: each ported tool's terminal target
  is 100% feature parity with its original *and* a POSIX-compliant CLI (see
  [docs/CLI_CONVENTIONS.md](docs/CLI_CONVENTIONS.md)). Until a tool gets there
  it is documented in [tools/PORTING_STATUS.md](tools/PORTING_STATUS.md) as
  "Partial".
- **Document AI Agent Utility**: Track and document the effectiveness (or lack thereof) of coding agents in this process

## Repository structure

Single Go module at the root (`github.com/yassineS/bio_ai_experiment`); no
per-tool `go.mod`, no third-party Go dependencies. Reference upstream sources
are vendored as git submodules under `reference_code/` for parity work.

```text
bio_ai_experiment/
├── go.mod                 # single root module
├── pkg/                   # shared libraries
│   ├── bioformats/        # fasta, fastq, vcf, bed, sam, bcf, iohelper
│   └── cliflag/           # POSIX short + GNU long flag helpers
├── tools/                 # tool ports, one subdir per tool
│   ├── PORTING_STATUS.md  # per-tool feature status
│   ├── PARITY_VALIDATION.md  # byte-for-byte audit results
│   └── <tool>/
│       ├── cmd/<tool>/main.go     # CLI entry
│       ├── pkg/<tool>/            # logic + tests
│       └── README.md
├── docs/
│   ├── PARITY_ROADMAP.md      # authoritative gap list
│   ├── UPSTREAM_BUGS.md       # bugs in originals we do not carry
│   ├── CLI_CONVENTIONS.md     # CLI rules
│   └── archive/               # historical Phase 0/1 docs
├── analysis/                  # tool ranking + research
├── reference_code/            # upstream sources as submodules (parity work)
└── .github/agents/            # AI-agent role descriptions
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
go doc github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta
```

## Current status

**No port in this repo is at 1:1 feature parity yet.** Every tool below has
a *working subset* of upstream functionality; the authoritative gap list is
[`docs/PARITY_ROADMAP.md`](docs/PARITY_ROADMAP.md), and upstream bugs we
identify but choose not to carry over are tracked in
[`docs/UPSTREAM_BUGS.md`](docs/UPSTREAM_BUGS.md).

### Shared libraries (`pkg/`)

- `pkg/htsgo/{fasta,fastq,vcf,bed,gff,sam,bcf,bgzf,tabix,bam,region,iohelper}` — parsers, writers,
  and a transparent gzip/BGZF I/O helper (BGZF auto-detected via the
  `BC` extra-subfield magic).
- `pkg/cliflag` — POSIX short + GNU long flag wiring on a standard
  `flag.FlagSet`.

### Tools ported

~50 tool ports as of 2026-06. The canonical completion table lives in
[`PROJECT_STATUS.md`](PROJECT_STATUS.md); [`tools/README.md`](tools/README.md)
indexes every tool, [`tools/PORTING_STATUS.md`](tools/PORTING_STATUS.md) is the
per-subcommand feature inventory, and
[`docs/PARITY_ROADMAP.md`](docs/PARITY_ROADMAP.md) is the authoritative gap
list. See [`docs/README.md`](docs/README.md) for the docs map (which doc owns
what).

Highlights:

- **htslib core**: `bgzip`, `tabix`, `htsfile`, `mosdepth`, `samtools`
  (24 functional subcommands incl. CRAM r/w, `mpileup`, `phase`,
  `consensus`), `bcftools` (24 subcommands incl. `call`, `mpileup
  --indels-cns`, `convert`, `gtcheck`, `csq`, `roh`).
- **bedtools surface**: ~37 `bed*` tools (`bedmerge`, `bedintersect`,
  `bedmap`, `bedcoverage`, `bedfisher`, …) covering the bedtools subcommand
  set.
- **Sequence preprocessing**: `seqtk` (byte-parity vs v1.5), `fastp`,
  `prinseq`, `sickle`, `skewer`.
- **VCF**: `vcftools` (146/146 upstream long flags).

### Validated parity against upstream test suites

Where we've run the upstream regression tests against our binaries and
diffed the output, the pass rates are:

| Audit | PR | Tests | Pass | Skip | Bugs fixed in our code |
|-------|----|------:|-----:|-----:|-----------------------:|
| bedtools (10 subcmds) | #55 | 127 | 85 | 42 | 7 |
| sickle + skewer | #73 | 27 | 22 | 3 | 3 |
| samtools (6 subcmds) | #75 | 43 | 34 | 9 | 3 |
| bcftools (6 subcmds) | #74 | 52 | 32 | 20 | 9 |
| mosdepth + vcftools | #76 | 65 | 50 | 15 | 9 |
| **Total** |  | **314** | **223** | **89** | **31** |

The 89 documented `t.Skip()`s are all unimplemented features, each
cross-referenced to `docs/PARITY_ROADMAP.md`. The 31 bug fixes were real
divergences between our ports and upstream surfaced by the audits and
corrected on the way (see each PR for details).

`prinseq`, `seqtk`, and `fastp` have **not yet** been parity-validated; the
common-path tests pass but 1:1 byte-equivalence with upstream is untested.

### CI

The CI workflow (`.github/workflows/ci.yml`) is currently a no-op
(`workflow_dispatch`-only) while the project is in heavy iteration; the full
job config (`gofmt -l`, `go vet ./...`, `go test -race -cover ./...`,
`go build ./...`, markdown lint over `**/*.md`) is kept commented in the file
ready to re-enable. Run those checks locally before pushing and record the
output in your PR.

## Documentation

- [`tools/README.md`](tools/README.md) — per-tool table + quick start
- [`tools/PORTING_STATUS.md`](tools/PORTING_STATUS.md) — feature status + coverage per tool
- [`tools/PARITY_VALIDATION.md`](tools/PARITY_VALIDATION.md) — byte-for-byte audit results vs upstream
- [`docs/PARITY_ROADMAP.md`](docs/PARITY_ROADMAP.md) — authoritative gap list
- [`docs/UPSTREAM_BUGS.md`](docs/UPSTREAM_BUGS.md) — bugs in upstream we deliberately do not carry
- [`docs/CLI_CONVENTIONS.md`](docs/CLI_CONVENTIONS.md) — CLI flag rules (POSIX short + GNU long)
- [`docs/GOLANG_GUIDE.md`](docs/GOLANG_GUIDE.md) — Go best practices used here
- [`analysis/tool_ranking_2026.md`](analysis/tool_ranking_2026.md) — ranking for *next* tools to port

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

## Contact

For questions, suggestions, or collaboration opportunities, please open an issue in this repository.
