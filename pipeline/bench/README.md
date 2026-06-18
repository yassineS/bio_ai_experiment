# Performance & scalability suite (`pipeline/bench`)

This package is the **performance** half of the parity pipeline. Where the
parity runner (`pipeline/cmd/parity-pipeline`) proves the Go ports produce the
*same output* as the vendored upstream binaries, this suite measures what they
*cost* — and how that cost scales with input size.

## What it measures

For every benchmark cell it times **our** binary against the **upstream** binary
and records three resource axes for each side, from `wait4(2)`'s `struct rusage`
(no `/usr/bin/time` dependency):

| axis | source | meaning |
|---|---|---|
| **wall-clock** | monotonic clock around the run | latency a user sees |
| **CPU time** | `ru_utime + ru_stime` | total compute (user + kernel) |
| **max RSS** | `ru_maxrss` | peak resident memory |

Wall and CPU are reduced to the **minimum** over repetitions (least
scheduler/IO noise — the standard "how fast can it go"); RSS to the **maximum**
(true peak). The report prints `ratio = ours / upstream` per axis (**< 1.0 means
we use less**).

## Scalability axis

The same matrix is swept across the fixture **scale tiers** (see
`pipeline/fixtures/scale.go`). Each tier scales read count, coverage, variant
count and interval count together, so running multiple tiers turns the point
measurements into **size-vs-resource curves**:

| tier | reference | reads | variants | intervals | ~footprint |
|---|---|---|---|---|---|
| smoke | 0.04 Mb | 2 k | 400 | 500 | <1 MB |
| small | 1 Mb | 40 k | 8 k | 6 k | ~5 MB |
| medium | 16 Mb | 300 k | 60 k | 40 k | ~50 MB |
| large | 192 Mb | 2.5 M | 400 k | 250 k | ~500 MB |

> **Reading the numbers.** At `smoke`/`small`, Go process startup, runtime init
> and GC dominate, so the ratios overstate our cost — the ports look *slower*
> simply because the work is over before the fixed overhead is amortised. The
> `medium` and `large` tiers are where steady-state throughput and memory
> behaviour show through; those are the tiers to quote.

## Coverage

Cells span every full-scale input format and the load-bearing operations,
including the same-file vs different-file overlap cases and the two-input set
operations:

- **BAM** — `samtools` view (BAM→BAM), sort, flagstat, stats, depth, mpileup
- **CRAM** — `samtools` view BAM→CRAM (reference-compressed encode) and CRAM→BAM (decode)
- **VCF** — `bcftools` view, norm, stats, query, call, **isec** (two indexed VCFs)
- **BED** — `bedtools` **intersect (same file)**, **intersect (different files)**, merge, coverage, genomecov, sort
- **FASTQ** — `seqtk seq`, `sickle se`

## Running it

```bash
# quick local sanity pass
go run ./pipeline/bench/cmd/parity-bench -scales small -reps 2

# the manuscript sweep (heavy — run on a fat node / HPC; large needs ~1 GB scratch)
go run ./pipeline/bench/cmd/parity-bench -scales medium,large -reps 5

# narrow to one format/cell while iterating
go run ./pipeline/bench/cmd/parity-bench -scales medium -group BED
go run ./pipeline/bench/cmd/parity-bench -scales medium -cell mpileup
```

Reports are written to `pipeline/.fixtures/<lastScale>/bench/bench.{json,md}`:
`bench.json` for downstream plotting/aggregation, `bench.md` for a
human-readable per-format table plus a scalability table of wall-time across the
swept tiers.

Fixtures are generated (and cached) on demand through the shared `fixtures`
package, so a bench run reuses whatever the parity runner produced for the same
tier and seed.

## Relationship to `bench_test.go`

`bench_test.go` holds `go test -bench` micro-benchmarks of the same operations
(handy under `go test`/CI and for `benchstat`). The `cmd/parity-bench` runner
adds the CPU + RSS axes, the cross-tier scalability sweep, and the JSON/Markdown
report that the micro-benchmark harness does not produce.
