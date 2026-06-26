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
     been **descoped**; standalone CLIs are the supported interface

## Ultimate Goals

- **Improve Usability**: Make bioinformatics tools more accessible and easier to use
- **Enhance Documentation**: Provide clear, comprehensive documentation for all tools
- **Boost Performance**: Leverage modern programming practices and languages for better performance
- **Feature Parity + POSIX-Compliant CLIs**: each ported tool's terminal target
  is 100% feature parity with its original *and* a POSIX-compliant CLI (see
  [docs/CLI_CONVENTIONS.md](docs/CLI_CONVENTIONS.md)). Until a tool gets there
  it is tracked in [docs/PARITY_ROADMAP.md](docs/PARITY_ROADMAP.md) as a
  remaining gap.
- **Document AI Agent Utility**: Track and document the effectiveness (or lack thereof) of coding agents in this process

## Repository structure

Single Go module at the root (`github.com/yassineS/bio_ai_experiment`); no
per-tool `go.mod`. Third-party Go dependencies are kept to an absolute minimum
— the standard library and the in-tree `pkg/` libraries are preferred, and the
only sanctioned externals are `gonum` (linear algebra for the vcftools PCA
family), `klauspost/compress` (the DEFLATE backend for BGZF/gzip I/O), and
`ulikunitz/xz` (LZMA decode confined to the CRAM codec layer). Reference
upstream sources are vendored as git submodules under `reference_code/` for
parity work.

```text
bio_ai_experiment/
├── go.mod                 # single root module (Go 1.24.x)
├── pkg/                   # shared libraries
│   ├── htsgo/             # fasta, fastq, vcf, bed, gff, sam, bcf, bam,
│   │                      #   bgzf, tabix, cram, region, iohelper, …
│   ├── cliflag/           # POSIX short + GNU long flag helpers
│   └── cppsort/           # C++-compatible sort helpers for byte-exact parity
├── tools/                 # tool ports, one subdir per tool
│   └── <tool>/
│       ├── cmd/<tool>/main.go     # CLI entry
│       ├── pkg/<tool>/            # logic + tests
│       └── README.md
├── pipeline/              # validation harnesses: GIAB concordance,
│                          #   differential fuzzing, conformance + edge cases
├── docs/
│   ├── PARITY_ROADMAP.md      # authoritative gap list
│   ├── UPSTREAM_BUGS.md       # bugs in originals we do not carry
│   ├── CLI_CONVENTIONS.md     # CLI rules
│   └── manuscript/            # manuscript plan + labeled bug corpus
├── analysis/                  # tool ranking + research
├── scripts/                   # repo scripts (e.g. recompute-metrics.sh)
├── reference_code/            # upstream sources as submodules (parity work)
└── .github/agents/            # AI-agent role descriptions
```

> The legacy `pkg/bioformats/` package no longer exists — all format parsing
> moved to `pkg/htsgo/` (the migration completed across PRs A–I).

## Getting Started

### Prerequisites

- Go 1.24+ (the exact toolchain is pinned in `go.mod`)
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

The `pkg/htsgo` library provides reusable parsers/writers for common file
formats (FASTA/FASTQ/VCF/BCF/BED/GFF/SAM/BAM/CRAM, plus BGZF/tabix and a
transparent gzip/BGZF + stdin/stdout I/O helper):

```bash
# View library documentation
go doc github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta
```

## Current status

The project is **no longer taking on new tools**; the focus is driving the
tools already started to **100% upstream feature parity** (every flag, input
format, and edge case), with bug fixes, better docs, and parity tests for
everything. Many tools are byte-exact against upstream today (the htslib core
and `seqtk` among them) and the live `upstream-parity` CI job re-checks this on
every run; the per-tool/per-feature state is tracked in
[`PROJECT_STATUS.md`](PROJECT_STATUS.md) (summary) and
[`docs/PARITY_ROADMAP.md`](docs/PARITY_ROADMAP.md) (the authoritative gap
list). Upstream bugs we identify but choose not to carry over are tracked in
[`docs/UPSTREAM_BUGS.md`](docs/UPSTREAM_BUGS.md).

### Shared libraries (`pkg/`)

- `pkg/htsgo/{fasta,fastq,vcf,bed,gff,sam,bcf,bam,bgzf,tabix,cram,region,iohelper}`
  — parsers, writers, and a transparent gzip/BGZF + stdin/stdout I/O helper
  (BGZF auto-detected via the `BC` extra-subfield magic).
- `pkg/cliflag` — POSIX short + GNU long flag wiring on a standard
  `flag.FlagSet`.

### Tools ported

53 drop-in CLIs as of 2026-06. The canonical completion table lives in
[`PROJECT_STATUS.md`](PROJECT_STATUS.md), and
[`docs/PARITY_ROADMAP.md`](docs/PARITY_ROADMAP.md) is the authoritative
per-tool/per-subcommand gap list. See [`docs/README.md`](docs/README.md) for the
docs map (which doc owns what); each tool also has its own
`tools/<tool>/README.md`.

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

> **Note — these are the original audit snapshots (PRs #55–#76).** The tool
> surface and parity have since expanded well beyond them: `bcftools` now ships
> 24 subcommands plus native (pure-Go) reimplementations of all 41 in-tree
> plugins, and the remaining parity tail — host-level region/target selection
> (`-r/-R/-t/-T`), `-W/--write-index`, multi-threshold `{…}` filter expansion,
> the full `split-vep` machinery, the `setGT`/`prune` RNG & LD modes (byte-exact
> via an in-tree `drand48` port), the `trio-dnm3` float de-novo models, and the
> htsgo uncompressed-BAM / canonical-CRAM-header writer paths — was closed in a
> later wave. For the current per-tool/per-feature state always trust
> `PROJECT_STATUS.md` and `docs/PARITY_ROADMAP.md`, not this historical table.

Since those original audits, live byte-for-byte parity harnesses have been
added for far more of the surface — including `seqtk` (glibc `drand48`
`sample`/`randbase` + Mott `trimfq`) and the `prinseq` corpus — and the
`pipeline/` validation suites (GIAB concordance, differential fuzzing,
htslib/htscodecs conformance, and a silent-corruption edge-case battery) extend
the gate beyond the upstream regression tests. `fastp` remains validated mainly
on the common path. Per-tool parity status is tracked in
[`PROJECT_STATUS.md`](PROJECT_STATUS.md) and
[`docs/PARITY_ROADMAP.md`](docs/PARITY_ROADMAP.md).

### CI

The CI workflow (`.github/workflows/ci.yml`) runs on every push to `main` and
every PR targeting `main`. Jobs:

- **gofmt + vet** — `gofmt -l` and `go vet ./...`.
- **test + cover** — `go test -coverprofile … ./...` over the whole module.
  The tool suites are not hermetic (their live-parity helpers compare against
  the real upstream binary), so this job pre-builds htslib/bcftools/samtools/
  bedtools once, serially, before the parallel test run.
- **race (pkg)** — `go test -race ./pkg/...` for the concurrent code (the
  parallel BGZF reader/writer and threaded readers).
- **build** — `go build ./...`.
- **macOS (build + short tests)** — `go build`/`go vet`/`go test -short ./...`
  on `macos-latest`: the tools are pure-Go (cgo-free), so this guards Darwin
  portability while `-short` skips the live-upstream parity tests (whose
  byte-exact goldens are Linux-built).
- **markdown lint** — `markdownlint-cli2` over `**/*.md` (excluding
  `reference_code`).
- **upstream parity (live)** — the independent, non-self-reported re-execution
  of the byte-exact gate: builds htslib/bcftools/samtools from the submodules
  and runs the `*Upstream*` parity tests against them.

Run `gofmt -w .`, `go vet ./...`, and `go test ./...` locally before pushing.

## Documentation

- [`PROJECT_STATUS.md`](PROJECT_STATUS.md) — top-level per-tool completion table
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
