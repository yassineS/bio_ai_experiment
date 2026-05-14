# Bioinformatics tools — Go ports

Go re-implementations of popular bioinformatics CLI tools. Each is in its own
subdirectory; they all live in the same Go module (`github.com/yassineS/bio_ai_experiment`),
no per-tool `go.mod`, no third-party dependencies.

## Tools in this directory

| Tool | Maps to | Purpose | Coverage |
|------|---------|---------|---------:|
| [`seqtk`](seqtk/) | seqtk | FASTA/Q processing (comp/fq2fa/seq/sample/trimfq/subseq/mergepe/cutN/mutfa/randbase/hpc) | ~70% |
| [`prinseq`](prinseq/) | PRINSEQ-lite | FASTA/Q stats + filtering | 99.9% |
| [`sickle`](sickle/) | sickle | Quality trimming | ~82% |
| [`skewer`](skewer/) | skewer | Adapter trimming | 100% |
| [`fastp`](fastp/) | fastp | All-in-one preprocessor (cut/sliding-window/HTML+JSON reports/auto-adapter) | ~76% |
| [`bedmerge`](bedmerge/) | bedtools merge | Merge overlapping intervals | 100% |
| [`bedintersect`](bedintersect/) | bedtools intersect | Interval intersection | 100% |
| [`bedsort`](bedsort/) | bedtools sort | Sort BED files | 92% |
| [`bedslop`](bedslop/) | bedtools slop | Extend intervals | 95% |
| [`bedcomplement`](bedcomplement/) | bedtools complement | Gaps over a genome | 95% |
| [`bedsubtract`](bedsubtract/) | bedtools subtract | A − B intervals | 94% |
| [`bedflank`](bedflank/) | bedtools flank | Flanking intervals | 92% |
| [`bedclosest`](bedclosest/) | bedtools closest | Nearest-neighbour intervals | 93% |
| [`bedgenomecov`](bedgenomecov/) | bedtools genomecov | Genome-wide coverage | 94% |
| [`bedjaccard`](bedjaccard/) | bedtools jaccard | Jaccard statistic on intervals | 96% |
| [`vcftools`](vcftools/) | vcftools | VCF filtering/stats/conversion + LD analysis | ~68% |
| [`bgzip`](bgzip/) | htslib `bgzip` | Block-gzip codec used by `.vcf.gz`, BAM, BCF, tabix | 90% |
| [`tabix`](tabix/) | htslib `tabix` | Region index for bgzipped VCF/BED/GFF/SAM | 86% |
| [`samtools`](samtools/) | htslib `samtools` | SAM/BAM `view`/`sort`/`index`/`depth`/`fastq`/`flagstat` | 87% |
| [`bcftools`](bcftools/) | htslib `bcftools` | VCF/BCF `view` (filtering, conversion, expression evaluator) | 85% |

For up-to-date per-tool feature lists and migration notes, see each tool's
own `README.md`, [`PORTING_STATUS.md`](PORTING_STATUS.md), and the
authoritative gap list in [`../docs/PARITY_ROADMAP.md`](../docs/PARITY_ROADMAP.md).

## Project goal: 1:1 feature parity

Every tool in this directory targets **1:1 feature parity** with its
upstream, validated byte-for-byte against the upstream's own test suite
where one exists. No port here is yet "done" by that bar; the gap list
per tool is in [`../docs/PARITY_ROADMAP.md`](../docs/PARITY_ROADMAP.md), and
bugs we identify in the upstream (which we do not carry over) are tracked
in [`../docs/UPSTREAM_BUGS.md`](../docs/UPSTREAM_BUGS.md).

The companion file [`../analysis/tool_ranking_2026.md`](../analysis/tool_ranking_2026.md)
ranks which **new** tools to port next; it is not a deprioritise-existing
filter.

## Building, testing, running

Run everything from the repo root.

```bash
# Build all tools
go build ./tools/...

# Build a single tool's binary
go build ./tools/seqtk/cmd/seqtk

# Run a tool without installing
go run ./tools/seqtk/cmd/seqtk comp file.fasta

# Test one tool
go test ./tools/seqtk/...

# Test one tool with race detector + coverage
go test -race -cover ./tools/seqtk/...

# Benchmark one tool
go test -bench=. -benchmem ./tools/seqtk/...
```

## Where the shared code lives

- [`pkg/bioformats/`](../pkg/bioformats/) — parsers + writers for FASTA, FASTQ,
  VCF, BED, **SAM/BAM**, **BCF**, plus `iohelper` for transparent gzip / BGZF
  (auto-detected by magic-byte sniff) and `-` (stdin/stdout).
- [`pkg/cliflag/`](../pkg/cliflag/) — POSIX short + GNU long flag wiring on a
  standard `flag.FlagSet`.

## Standard tool layout

```text
tools/<tool>/
├── cmd/<tool>/main.go         # CLI entry point
├── pkg/<tool>/                # tool logic + tests
│   ├── <tool>.go
│   └── <tool>_test.go
└── README.md
```

Older docs sometimes describe a deeper structure (`tests/unit/`, `testdata/`,
`docs/`, per-tool `go.mod`/`go.sum`). That structure was never used; trust
the layout above.

## CLI conventions

All tools must support **POSIX short flags** (`-i`, `-o`, ...) **and** GNU
long flags (`--input`, `--output`, ...) for the same option. See
[`../docs/CLI_CONVENTIONS.md`](../docs/CLI_CONVENTIONS.md) for the full
rules.

## CI

The CI workflow is currently disabled (manual-only via
`workflow_dispatch`). Contributors run `gofmt -l`, `go vet`, `go test -race
-cover`, `go build`, and `markdownlint` locally and document the output in
each PR description.
