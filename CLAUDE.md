# CLAUDE.md

Guidance for AI assistants (Claude Code and others) working in this repository.

## What this project is

`bio_ai_experiment` is an experiment in using AI agents to re-implement popular
bioinformatics command-line tools in Go — making them faster, better tested, and
better documented than their originals (typically C, C++, or Perl). The original
tools are vendored as git submodules under `reference_code/` for reference.

## Repository layout

```
.
├── go.mod                    # single Go module: github.com/yassineS/bio_ai_experiment
├── pkg/                      # shared libraries used by all tools
│   ├── bioformats/           # parsers/writers: fasta, fastq, vcf, bed, iohelper
│   └── cliflag/              # helper for flags that accept both short (-i) and long (--input) forms
├── tools/                    # Go re-implementations, one subdir per tool
│   └── <tool>/
│       ├── cmd/<tool>/main.go        # CLI entry point
│       └── pkg/<tool>/*.go           # tool logic + *_test.go
├── reference_code/           # git submodules: upstream sources of the original tools
├── analysis/                 # Python scripts + CSV/JSON/Markdown ranking the top ~200 tools
├── mcp-servers/              # MCP server implementations (planned; currently just a README)
├── docs/                     # project docs (see below)
└── .github/agents/           # role descriptions for the AI agents that build this repo
```

Implemented tools so far: `seqtk`, `prinseq`, `sickle`, `skewer`, `fastp`,
`bedmerge`, `bedintersect`, `vcftools`. Each has its own `README.md` under
`tools/<tool>/`.

### Important: it is ONE Go module

Despite what some older docs (`tools/README.md`, `CONTRIBUTING.md`,
`.github/agents/*`) say, there is **no per-tool `go.mod`**. Everything lives in
the single root module. Third-party dependencies are kept to an absolute
minimum — prefer the standard library and the shared `pkg/` libraries.
Adding a new external dependency needs its own conversation; the currently
sanctioned set is enumerated below.

#### Sanctioned third-party deps

Adding a new external dependency outside these zones requires explicit
owner approval.

1. **`gonum.org/v1/gonum`** (BSD-3) — used by the vcftools `--pca`
   family for symmetric eigendecomposition of the N×N Genomic
   Relatedness Matrix. Upstream calls LAPACK `dgeev_`; gonum's pure-Go
   `mat.SymEigen` is the equivalent without dragging in cgo / libblas.
   Confined to `tools/vcftools/pkg/vcftools/pca.go` today; future
   stats-heavy tools (relatedness PCs, Fst, ADMIXTURE-style models)
   are pre-approved to reuse the same dep.
2. **CRAM codec layer** (when we get there). CRAM uses several custom
   compression codecs (rANS 4x8, rANS 4x16) that have no Go-stdlib
   equivalent and would otherwise require ~1,500 lines of careful
   bit-twiddling to port from htslib's C. The owner has explicitly
   OK'd accepting third-party deps for these primitives rather than
   a from-scratch rANS port. The dep must be confined to a single
   sub-package under `pkg/bioformats/cram/codec/` so the rest of the
   repo can still claim "stdlib + gonum only" for non-CRAM workflows.
   See `docs/CRAM_DESIGN.md`.

Preference order:

1. Standard library.
2. In-tree implementation (the bgzip / tabix / sam / bcf packages all
   went this route).
3. Sanctioned third-party dep (gonum for linalg, future CRAM codec
   libraries).

## Common commands (run from repo root)

```bash
go build ./...                                   # build everything
go test ./...                                    # run all tests
go test -race -cover ./...                        # race detector + coverage (matches CI)
go test ./tools/seqtk/...                          # test one tool
go test -bench=. -benchmem ./tools/<tool>/...      # benchmarks
gofmt -l .                                        # list unformatted files (CI fails if any)
gofmt -w .                                        # format
go vet ./...                                      # vet
go build ./tools/seqtk/cmd/seqtk                  # build a single tool binary
go run ./tools/seqtk/cmd/seqtk comp file.fasta    # run a tool without installing
```

CI (`.github/workflows/ci.yml`) runs: `gofmt -l`, `go vet ./...`,
`go test -race -coverprofile=... ./...`, `go build ./...`, and markdown lint
(`**/*.md`). Keep markdown well-formed. Note: CI pins Go 1.21 while `go.mod`
declares 1.24.9 — stick to language/stdlib features available in 1.21.

## Conventions

### Code

- Idiomatic Go, `gofmt`-clean, passes `go vet`. Small focused functions, meaningful names.
- Document all exported identifiers with complete-sentence doc comments.
- Tool logic goes in `tools/<tool>/pkg/<tool>/`; `cmd/<tool>/main.go` only does
  argument parsing, wiring, and exit codes.
- Reuse `pkg/bioformats/*` for file format parsing instead of re-implementing it.
  Use `pkg/bioformats/iohelper` for transparent gzip and stdin/stdout (`-`) handling.

### CLI design

- Tools should accept both POSIX short flags (`-i`, `-o`, `-q`, `-l`, ...) and
  GNU long flags (`--input`, `--output`, `--min-length`, ...). Use the
  `pkg/cliflag` helpers (`cliflag.StringVar`, `cliflag.IntVar`, `cliflag.Float64Var`,
  `cliflag.BoolVar`) which register both names on a `flag.FlagSet`.
- `-` means stdin/stdout. Always provide `-h/--help` and `-v/--version`.
- Where practical, keep compatibility with the original tool's flags.
- Full rules: `docs/CLI_CONVENTIONS.md`.

### Tests

- Every new behavior gets a test. Prefer table-driven tests. Aim for >80% coverage.
- Tests live next to the code as `*_test.go` in `tools/<tool>/pkg/<tool>/`.
- Add benchmarks for performance-sensitive paths; this project cares about being
  faster than the originals (see `PERFORMANCE_BENCHMARKS.md`).

### Commits & PRs

- Imperative, concise commit subjects (e.g. "Add reverse complement to seqtk seq").
- Branch off `main`. PRs target `main`; there's a PR template in `.github/`.
- Run `go test ./...`, `go vet ./...`, and `gofmt -w .` before committing.

## Submodules

`reference_code/` contains the upstream tool sources as submodules. They're
read-only references for behavior/parity checking — never modify them. Init with
`git submodule update --init reference_code/<tool>` only if you actually need a
specific one; you usually don't.

## Where to look

- `README.md` — project overview and per-tool quick start.
- `docs/GOLANG_GUIDE.md` — Go patterns/best practices used here.
- `docs/CLI_CONVENTIONS.md` — the canonical CLI flag spec.
- `pkg/bioformats/README.md` — format library docs.
- `tools/PORTING_STATUS.md`, `tools/IMPLEMENTATION_SUMMARY.md` — tool-by-tool status.
- `tools/<tool>/README.md` — per-tool usage and parity notes.
- `.github/agents/*.md` — the agent roles (tool-analysis, golang-recoding, testing,
  documentation, mcp-server) and how the work is divided.

## Caveats / known stale docs

Some documents predate the current code and describe an aspirational structure
that wasn't followed:

- `tools/README.md` claims the `tools/` dir is empty and that each tool has its
  own `go.mod`/`go.sum` and a `tests/` + `testdata/` + `docs/` subtree. In
  reality tools are populated, share the root module, and keep tests inline.
- `PROJECT_STATUS.md` says "0 tools implemented" — outdated.

When in doubt, trust the actual code over these older Markdown files, and feel
free to update the stale docs as part of related work.
