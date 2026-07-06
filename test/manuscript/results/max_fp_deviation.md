# Max absolute / relative floating-point deviation per tool

The honest "worst-case how far off" — beyond the pass/fail within-ε verdict.
This aggregates, **per tool**, the largest absolute (`|ours − upstream|`) and
largest relative (`|ours − upstream| / max(|ours|,|upstream|)`) numeric-field
deviation observed across the parity matrix and the real-data parity battery.

## How the deltas are computed

The similarity comparator lives in `pipeline/runner/compare.go`:

- `CompareSimilarity(ours, upstream, eps)` tokenises both provenance-stripped
  streams field-by-field, parses numeric fields, and records **both** the max
  relative deviation (`relDev`, `MaxDeviation`) and the max absolute deviation
  (`math.Abs(a-b)`, `MaxAbsDeviation`). Non-numeric fields must match exactly.
- The default relative tolerance is `similarityEpsilon = 1e-6`; a matrix entry
  may widen it via `Entry.Tolerance` (`resolveEpsilon`).
- Byte-exact tools go through `CompareByteExact` / `CompareDigests` (md5 of the
  provenance-stripped stream). Any float divergence there is a **DIVERGE**, so
  their accepted max deviation is **exactly 0.0** by construction.
- `pipeline/cmd/realparity` (the real-data samtools/bcftools battery) is
  **strictly byte-exact** — it compares md5 digests of the provenance-stripped
  streams (`runner.CompareByteExact` semantics). It therefore contributes 0.0 to
  every tool it exercises; a last-ULP float difference would fail as a DIVERGE
  rather than be absorbed.

`MaxAbsDeviation` was added to `CompareResult` / `Result` (JSON key
`max_abs_deviation`) alongside the pre-existing relative `MaxDeviation` so both
are recorded per cell; the per-tool maxima below are the max over each tool's
cells.

## Per-tool max deviation

Only three matrix cells across the whole sweep run in `Similarity` mode (every
other cell is byte-exact); all remaining tools are exactly 0.0.

| Tool | Max \|abs\| | Max \|rel\| | Exact (0.0)? | Notes |
|---|---|---|---|---|
| `bcftools` (`call -m`) | ~1e-4 phred | ~7.2e-6 | no | libm last-ULP QUAL residual (tol 2e-5) |
| `bedgenomecov` (histogram fraction) | 1.0e-6 | 9.7e-6 | no | `%g` round-half last-digit flip (tol 1e-5) |
| `samtools` (`consensus`, gap5 bayesian `cq`) | ≤4 (cq units); ±1 typical | discrete phred | no† | libm last-ULP; base/seq/qual byte-exact |
| `vcftools` (all stats/PCA/Fst cells) | 0.0 | 0.0 | **yes** | byte-exact |
| `mosdepth` | 0.0 | 0.0 | **yes** | byte-exact |
| `samtools` (all other subcommands) | 0.0 | 0.0 | **yes** | byte-exact |
| `bcftools` (all other subcommands) | 0.0 | 0.0 | **yes** | byte-exact |
| all other `bed*` tools | 0.0 | 0.0 | **yes** | byte-exact |
| `seqtk`, `fastp`, `prinseq`, `sickle`, `skewer` | 0.0 | 0.0 | **yes** | byte-exact |
| `bgzip`, `tabix`, `htsfile` | 0.0 | 0.0 | **yes** | byte-exact |

† `samtools consensus` in the default gap5 Bayesian mode is byte-exact for the
base / sequence / quality bytes; only the derived `cq` score column carries the
residual (it is not a `Similarity` matrix cell — it is validated by the
`consensus`/roundtrip parity tests, and `--mode simple` is fully byte-exact).

## Which tools carry a residual, and why

Three residuals exist, all of the **same class**: Go's `math.Log` / `math.Exp` /
`math.Pow` (and the C++ `ostream %g` rounding) differ from C glibc `libm` in the
last ULP, and glibc itself is not bit-stable across versions. These are accepted
proximity-parity residuals, not bugs:

1. **`bcftools call -m` — QUAL.** `QUAL = -4.343·(ref_lk − logsumexp2(...))` is
   an accumulation of `log`/`exp` and a `pl2p` pow table that are not
   bit-identical to glibc. Values are correct to ~6 significant figures; the
   observed worst relative deviation is **7.2e-6** (matrix tolerance 2e-5, ≈3×
   headroom, still orders of magnitude below any real-bug field move). Max
   absolute wobble on a QUAL near 15.7 is ~**1e-4** phred. Source:
   `pipeline/matrix/bcftools.go` (the `bcftools_call` entry) and
   `pipeline/edgecases/qual_pl_ulp_test.go` (which additionally asserts the
   GT/FILTER *decision* is byte-identical and |ΔQUAL| < 0.01 phred).

2. **`bedGenomeCov` histogram — depth-fraction column.** The `chrom / depth /
   count / genomeSize` columns are byte-identical; only the `%g`-printed
   `fraction` column flips its last significant digit on exact round-half values
   (Go rounds half-to-even, C++ `ostream %g` rounds half-up — e.g. 0.103500 vs
   0.103501). Worst relative deviation **9.7e-6**, worst absolute **1.0e-6**
   (matrix tolerance 1e-5). Source: `pipeline/matrix/bedtools.go` (the
   `base` / `flagmax_5` / `flagstrand_+` genomecov entries).

3. **`samtools consensus` — gap5 Bayesian `cq` score.** Whole-contig-20, ours vs
   upstream differs at **68 / 10 416 641 positions, in the `cq` column only** —
   every base, sequence and quality byte is byte-exact. The differences are
   **bidirectional ±1** (ruling out a fixable rounding-mode bias; up to ±2–4 at a
   handful of ultra-deep columns), the signature of `log`/`exp` last-ULP noise
   landing the final phred on the other side of an integer boundary. Same class
   the project accepts for the `bcftools trio-dnm3` float scores. Closing it
   would require a bit-for-bit `libm` reimplementation (out of scope); `--mode
   simple` is byte-exact. Source: `tools/samtools/pkg/samtools/consensus_bayesian.go`
   (the ACCEPTED RESIDUAL note) and `docs/PARITY_ROADMAP.md`.

**Everything else is byte-exact (max abs = max rel = 0.0).** The float-output
tools that one might expect to carry residuals — `vcftools` stats/Fst/PCA,
`mosdepth`, `samtools stats`, `bcftools stats` — are all byte-exact under the
provenance-stripped comparison; their numeric output matches upstream to the
last printed digit.

## Reproducing command

The per-cell abs/rel deltas are emitted by the runner (JSON keys
`max_deviation`, `max_abs_deviation`). For the two `Similarity` matrix cells the
documented worst-case deviations reproduce directly through the comparator:

```go
// bedGenomeCov histogram worst case: rel 9.66e-6, abs 1e-6
CompareSimilarity([]byte("chr1\t1\t2\t0.103500\n"),
                  []byte("chr1\t1\t2\t0.103501\n"), 1e-5)
// bcftools call QUAL last-ULP: rel ~6.4e-6, abs ~1e-4 phred
CompareSimilarity([]byte("x\t20\t.\tA\tG\t15.6999\n"),
                  []byte("x\t20\t.\tA\tG\t15.6998\n"), 2e-5)
```

To drive the full matrix / real-data battery (requires the upstream oracles from
`test/nextflow/Dockerfile`):

```bash
# real-data samtools/bcftools parity battery (byte-exact -> 0.0 deviation)
go run ./pipeline/cmd/realparity -samtools-ours=bin/ours/samtools \
    -samtools-up=bin/upstream/samtools -bam=<bam> -ref=<ref> -vcf=<vcf>

# the matrix Similarity cells (bcftools_call, bedgenomecov base/flag*) run via
# the runner's *Upstream* parity tests; their per-cell max_deviation /
# max_abs_deviation appear in the emitted JSON Result.
go test ./tools/bcftools/... ./tools/bedgenomecov/... -run Parity
```
