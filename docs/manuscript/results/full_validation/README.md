# Full validation + performance round

Run after the real-data bug-fix rounds, on the `linux/arm64` container
([`../hardware.md`](../hardware.md)), **sequentially / uncontended** so the
performance numbers are not skewed by competing load. Two stages:

1. **Full `go test ./...`** — whole-repo regression.
2. **`full-validation -scales=smoke,small,medium -reps=10`** — the parity matrix
   (all 53 CLIs) + bidirectional round-trip interop + the wall/CPU/RSS bench.
   Reports here: `report.{json,md}`, `roundtrip.md`, `bench.{json,md}`.

## 1. `go test ./...` — 102 packages OK, 7 FAIL (all pre-existing arm64 artifacts)

The 7 failures are **not** regressions from any fix — they are the known
`arm64`-platform-sensitive parity tests (byte-exact on `amd64`/CI), in files no
fix touched:

- `vcftools` — `TestParity_SiteMeanDepth`, `TestParity_OutputModes_Upstream`,
  `TestParity_WeirFst_Upstream` (arm64 `-nan` / last-ULP FP).
- `samtools` — `TestConsensus_AllPositionsUpstreamParity`,
  `TestLiveCoverageHistogram`.
- `bcftools` plugins — `TestVrfsOracleParity`, `TestNativePluginColorChrs`
  (arm64 `-nan`).

## 2. full-validation

| scale | parity | round-trip interop | perf bench |
|---|---|---|---|
| smoke | 382 pass / 4 similar / 3 diverge / 0 error | **14/14 PASS** | 22 cells |
| small | 381 pass / 4 similar / 4 diverge / 0 error | **14/14 PASS** | 22 cells |
| medium | 382 pass / 4 similar / 3 diverge / 0 error | **14/14 PASS** | 22 cells (peak heap 10.7 GB) |

**The diverges are all explained, none are real port defects:**

- `fastp_detect_adapter_pe_heavy` — a **stale test-matrix false positive**: the
  cell hard-coded the *old* non-upstream fastp CLI for our side, and the new
  drop-in CLI made `-I` mean read-2, scrambling the I/O. **Fixed** in the matrix
  (commit `c34070e`); verified the cell is byte-exact (r1+r2 md5-identical to
  upstream). This run pre-dates that fix, so it still appears here.
- `vcftools_site_mean_depth`, `vcftools_hardy`, `vcftools_hap_r2` — the documented
  **arm64 `-nan` / FP** cells (byte-exact on `amd64`).
- The 4 `SIMILAR` cells are the documented floating-point cells (`bcftools call`;
  three `bedtools` similarity cells), accepted by the gate.

So the genuine parity result is **all cells pass except the arm64-FP vcftools
cells**, with **round-trip interop 14/14** at every scale (BGZF/BAM/CRAM/VCF.gz/
BCF/FASTQ + BAI/CSI/TBI region queries, both directions).

## 3. Performance (`bench.md`, medium, reps=10) — `ratio = ours/upstream`

| cell | wall× | note |
|---|---|---|
| `bed_intersect_pair` | 0.61 | faster |
| `bed_intersect_self` | 0.63 | faster |
| `sam_sort` | 0.68 | faster |
| `sam_view_bam2bam` | 0.70 | faster |
| `sam_view_bam2cram` | 0.73 | faster |
| `bed_coverage` | 0.73 | faster |
| `sam_depth` | 1.01 | par |
| `bed_genomecov` | 1.03 | par |
| `bed_sort` | 1.05 | par |
| `sam_stats` | 1.07 | par (CPU× 1.64 — more threads) |
| `bcf_view` | 1.12 | slower |
| `sam_flagstat` | 1.31 | slower |
| `sam_view_cram2bam` | 1.32 | slower |
| `bcf_norm` | 1.42 | slower |
| `bcf_query` | 1.56 | slower |
| `bcf_call` | 1.75 | **slow** |
| `bed_merge` | 2.05 | **slow** |
| `sam_mpileup` | 2.40 | **slow** (CPU× 3.35) |
| `bcf_isec` | 2.99 | **slowest** (CPU× 4.88) |

I/O-bound conversions and `bedtools` intersect/coverage are **faster** than
upstream; the compute-heavy variant cells (`mpileup`, `call`, `isec`) are
**slower** and reported plainly.

**Caveats (read before citing):**

- Medium ran under `GOMEMLIMIT=8GiB` to fit the 12 GB VM, which adds GC pressure
  that **penalises the heavy cells** (`mpileup`/`call`/`isec` would look better
  unconstrained). The light/fast cells (well under the limit) are unaffected.
- The `RSS×` column reads the **orchestrator's** RSS, not per-subprocess, so it is
  a uniform 1.00 and is **not** a meaningful memory comparison.
- Timings reduce to min-over-reps; the median/IQR/ratio-CI upgrade (H1a) is still
  a separate task.
- `arm64` container, not a fat node — `mpileup` etc. are absolute-time slow here.
