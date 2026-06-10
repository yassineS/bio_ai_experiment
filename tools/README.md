# Bioinformatics tools — Go ports

Go re-implementations of popular bioinformatics CLI tools. Each is in its own
subdirectory; they all live in the same Go module (`github.com/yassineS/bio_ai_experiment`),
no per-tool `go.mod`. Third-party deps are kept to the sanctioned minimum
(`gonum` for linalg; `ulikunitz/xz` confined to the CRAM codec layer) — see
[`../CLAUDE.md`](../CLAUDE.md).

> **Source of truth for status.** The parity percentages below are a quick
> index only. The canonical, skimmable completion table lives in
> [`../PROJECT_STATUS.md`](../PROJECT_STATUS.md), and the authoritative
> per-tool gap list lives in [`../docs/PARITY_ROADMAP.md`](../docs/PARITY_ROADMAP.md).
> When the numbers here disagree, those two files win.

## Tools in this directory

| Tool | Maps to | Purpose | Parity |
|------|---------|---------|---------:|
| [`seqtk`](seqtk/) | seqtk | FASTA/Q processing (all 24 subcommands; byte-parity vs v1.5) | 100% |
| [`prinseq`](prinseq/) | PRINSEQ-lite | FASTA/Q stats + filtering | ~95% |
| [`sickle`](sickle/) | sickle | Quality trimming | 100% |
| [`skewer`](skewer/) | skewer | Adapter trimming | 100% |
| [`fastp`](fastp/) | fastp | All-in-one preprocessor (cut/sliding-window/HTML+JSON reports/auto-adapter) | ~85% |
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
| [`bedgroupby`](bedgroupby/) | bedtools groupby | Group records + apply column ops | 94% |
| [`bed12tobed6`](bed12tobed6/) | bedtools bed12tobed6 | Split BED12 records into BED6 blocks | 91% |
| [`bedmakewindows`](bedmakewindows/) | bedtools makewindows | Partition intervals into windows | 91% |
| [`bedexpand`](bedexpand/) | bedtools expand | Expand comma-separated list columns into rows | 92% |
| [`bedgetfasta`](bedgetfasta/) | bedtools getfasta | Pull FASTA subsequences for BED intervals | 86% |
| [`bedsample`](bedsample/) | bedtools sample | Random subsample (reservoir, without replacement) | 91% |
| [`bedspacing`](bedspacing/) | bedtools spacing | Distance to previous interval per chromosome | 85% |
| [`bedcoverage`](bedcoverage/) | bedtools coverage | Per-A coverage (count / bp / fraction / histogram / per-base) | 92% |
| [`bedmap`](bedmap/) | bedtools map | Apply column ops to A using B's values | 92% |
| [`bedshuffle`](bedshuffle/) | bedtools shuffle | Randomly relocate intervals across a genome | 92% |
| [`bedcluster`](bedcluster/) | bedtools cluster | Cluster overlapping intervals + tag cluster ID | ~90% |
| [`bedsplit`](bedsplit/) | bedtools split | Split a BED into N approximately equal shards | ~90% |
| [`bedsummary`](bedsummary/) | bedtools summary | Per-chromosome interval-length summary stats | ~90% |
| [`bedtag`](bedtag/) | bedtools tag | Annotate A records with tags from overlapping B | ~90% |
| [`bedwindow`](bedwindow/) | bedtools window | Overlap A with an expanded window around B | ~90% |
| [`bedreldist`](bedreldist/) | bedtools reldist | Relative-distance distribution between two BED files | ~91% |
| [`bedfisher`](bedfisher/) | bedtools fisher | Fisher's exact test of overlap enrichment | ~89% |
| [`bednuc`](bednuc/) | bedtools nuc | Nucleotide content profile per BED interval | ~95% |
| [`bedannotate`](bedannotate/) | bedtools annotate | Annotate A intervals with overlap stats from N BED files | ~89% |
| [`bedmulticov`](bedmulticov/) | bedtools multicov | Per-interval overlap counts against N BED files | ~88% |
| [`bedmultiinter`](bedmultiinter/) | bedtools multiinter | Multi-way intersection across N BED files | ~87% |
| [`bedigv`](bedigv/) | bedtools igv | Emit an IGV batch-mode script (one snapshot per interval) | 90% |
| [`bedlinks`](bedlinks/) | bedtools links | Emit a UCSC Genome Browser HTML link table | 87% |
| [`bedpairtobed`](bedpairtobed/) | bedtools pairtobed | Overlap BEDPE pairs against a BED (either/both/neither/notboth/xor) | 90% |
| [`bedpairtopair`](bedpairtopair/) | bedtools pairtopair | Overlap two BEDPE files (both/either/neither/notboth + slop) | 91% |
| [`vcftools`](vcftools/) | vcftools | VCF filtering/stats/conversion + LD/RoH/relatedness/PCA (146/146 long flags) | ~97% |
| [`bgzip`](bgzip/) | htslib `bgzip` | Block-gzip codec used by `.vcf.gz`, BAM, BCF, tabix | ~92% |
| [`tabix`](tabix/) | htslib `tabix` | Region index for bgzipped VCF/BED/GFF/SAM | ~92% |
| [`htsfile`](htsfile/) | htslib `htsfile` | Identify file format / compression | ~98% |
| [`mosdepth`](mosdepth/) | mosdepth | Fast BAM/CRAM depth over windows/BED (`-d/--d4` byte-identical) | ~85% |
| [`samtools`](samtools/) | htslib `samtools` | SAM/BAM/CRAM tools — 24 functional subcommands: `view`/`sort`/`index`/`depth`/`fastq`/`flagstat`/`merge`/`coverage`/`idxstats`/`cat`/`reheader`/`addreplacerg`/`fixmate`/`dict`/`split`/`quickcheck`/`mpileup`/`markdup`/`stats`/`calmd`/`import`/`phase`/`targetcut`/`consensus` | ~88% |
| [`bcftools`](bcftools/) | htslib `bcftools` | VCF/BCF tools — 24 subcommands: `view`/`index`/`stats`/`query`/`concat`/`norm`/`call`/`annotate`/`head`/`isec`/`merge`/`reheader`/`sort`/`filter`/`consensus`/`convert`/`mendelian`/`mendelian2`/`gtcheck`/`roh`/`polysomy`/`cnv`/`csq`/`mpileup` (+ `plugin`) | ~70% |

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

- [`pkg/htsgo/`](../pkg/htsgo/) — parsers + writers for FASTA, FASTQ,
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

The CI workflow (`.github/workflows/ci.yml`) is currently disabled
(manual-only via `workflow_dispatch`) while the project iterates; the full
check set is kept commented in the file. Run it locally before pushing and
document the output in the PR:

```bash
gofmt -l .
go vet ./...
go test -race -cover ./...
go build ./...
```

Stick to Go 1.21 language/stdlib features (the commented CI config pins 1.21).
