# Statement coverage over the parity-exercised code

**Headline: 64.25 % of statements** (63 319 / 98 545) are covered across the
tool packages (`tools/*/pkg/*` + `cmd`), the shared format libraries
(`pkg/htsgo/*`), and the parity runner (`pipeline/runner`) — i.e. the code the
parity sweep drives, including every `*Upstream*` / `*Parity` test that compares
our output against the vendored real binaries.

Go's coverage tooling measures **statement** coverage, which is the standard
proxy for branch coverage (see the interpretation note below).

## Reproducing command

Run from the repo root:

```bash
go test \
  -coverpkg=./tools/...,./pkg/htsgo/...,./pipeline/runner \
  -coverprofile=cover.out \
  ./tools/... ./pkg/htsgo/... ./pipeline/runner
go tool cover -func=cover.out | tail -1   # -> total: (statements) 64.2%
```

Notes:

- `-coverpkg` instruments the union of the tool packages, the htsgo format
  libraries, and the parity runner, so a line counts as covered whenever **any**
  test in the run exercises it — the parity/`*Upstream*` tests are a subset of
  that run, and they are the tests that drive the real-data code paths (BAM/CRAM
  decode, VCF/BCF, BED, SAM records) against the live upstream oracles.
- Because `-coverpkg` names many packages, every per-package test binary
  re-emits the full instrumented block set, so the raw `cover.out` contains one
  block record **per test binary**. `go tool cover -func` deduplicates by
  merging identical block keys (OR-ing the covered flag); that merged view is
  the authoritative 64.25 %. A naïve `awk` sum over the raw profile double-counts
  and must merge by block key first (the per-package table below does exactly
  that, and its TOTAL reconciles to `-func`'s 64.25 %).
- One environment-dependent CRAM test (`pkg/htsgo/cram` bzip2 parity) can FAIL
  when the vendored upstream `samtools` does not select a BZIP2 external block;
  it does not affect the emitted profile or the parity comparator and is
  orthogonal to this measurement.

## Per-tool / per-package coverage (statement-weighted)

Rolled up by tool family (the 37 `bed*` tools are grouped):

| Group | Statements covered | Coverage |
|---|---|---|
| `pkg/htsgo/*` (format libs) | 14 124 / 16 511 | **85.5 %** |
| `tools/mosdepth` | 1 103 / 1 290 | **85.5 %** |
| `tools/vcftools` | 4 052 / 4 929 | **82.2 %** |
| `tools/prinseq` | 3 318 / 4 575 | 72.5 % |
| `tools/tabix` | 148 / 224 | 66.1 % |
| `tools/fastp` | 1 752 / 2 627 | 66.7 % |
| `tools/bed*` (37 bedtools) | 7 795 / 12 172 | 64.0 % |
| `tools/samtools` | 8 169 / 12 899 | 63.3 % |
| `tools/htsfile` | 143 / 243 | 58.8 % |
| `tools/skewer` | 816 / 1 432 | 57.0 % |
| `tools/seqtk` | 1 808 / 3 185 | 56.8 % |
| `tools/bcftools` | 19 490 / 37 053 | 52.6 % |
| `tools/bgzip` | 171 / 252 | 67.9 % |
| `tools/sickle` | 222 / 666 | 33.3 % |
| `pipeline/runner` | 208 / 487 | 42.7 % |
| **TOTAL** | **63 319 / 98 545** | **64.25 %** |

Selected `pkg/htsgo/*` breakdown (the format libraries the real-data paths lean
on hardest):

| Package | Coverage |
|---|---|
| `pkg/htsgo/errmod` | 96.4 % |
| `pkg/htsgo/cram` | 88.5 % |
| `pkg/htsgo/baq` | 88.7 % |
| `pkg/htsgo/hfile` | 89.0 % |
| `pkg/htsgo/tabix` | 85.4 % |
| `pkg/htsgo/bam` | 84.2 % |
| `pkg/htsgo/bgzf` | 84.3 % |
| `pkg/htsgo/bcf` | 83.6 % |
| `pkg/htsgo/sam` | 82.0 % |
| `pkg/htsgo/vcf` | 71.4 % |
| `pkg/htsgo/fasta` | 76.6 % |
| `pkg/htsgo/fastq` | 62.2 % |

The full per-package numbers (all 75 packages) are reproducible from
`go tool cover -func=cover.out`.

## Interpretation and caveats

The headline 64.25 % is the fraction of executable **statements** in the
parity-exercised code that at least one test runs. It is the sharpest available
answer to "which input regions are untested": the uncovered ~36 % concentrates
in (a) rarely-hit error/EOF branches and defensive `return err` paths, (b)
seldom-used flag combinations on the largest surfaces (`bcftools` at 52.6 % and
`samtools` at 63.3 % carry the two biggest bodies of code and the longest tail
of sub-command flags), and (c) `sickle`/`bedshift`/`bedslop`-style tools whose
crafted-input tests cover the mainline but not every edge branch.

Two caveats matter for reading this as a *branch*-coverage figure:

1. **Statement coverage is a proxy for branch coverage, not a substitute.** Go
   instruments basic blocks, so a fully-covered `if`/`switch` statement is
   counted covered even if only one arm was taken; true branch coverage is
   therefore **≤** this number. The figure is an upper bound on branch coverage,
   useful for locating untested *regions* rather than certifying every decision
   edge.
2. **The parity oracle drives the real-data paths, which is where coverage
   matters most.** The `*Upstream*`/`*Parity` tests execute the actual BAM/CRAM,
   VCF/BCF, BED and SAM parsing/formatting code against live upstream binaries,
   so the covered statements are disproportionately the correctness-critical
   ones. Coverage here is a floor of confidence, complemented by the byte-exact
   parity gate (which asserts *behaviour*, not merely *execution*): a covered
   line that produced a wrong byte would still fail parity, so the two metrics
   are orthogonal and together stronger than either alone.
