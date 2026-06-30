# Full validation + performance round

Run after the real-data bug-fix rounds, on the `linux/arm64` container
([`../hardware.md`](../hardware.md)), **sequentially / uncontended** so the
performance numbers are not skewed by competing load. Two stages:

1. **Full `go test ./...`** — whole-repo regression.
2. **`full-validation -scales=smoke,small,medium -reps=10`** — the parity matrix
   (all 53 CLIs) + bidirectional round-trip interop + the wall/CPU/RSS bench.
   Reports here: `report.{json,md}`, `roundtrip.md`, `bench.{json,md}`, plus the
   raw `run.log`.

> This is the **clean rerun** under the fixed `fastp` test-matrix (commit
> `c34070e`) and the new robust perf statistics (H1a, commit `e34e62d`). The
> earlier round's `fastp_detect_adapter_pe_heavy` divergence was a stale
> test-matrix false positive and is **gone** — that cell now passes at every
> scale. The only remaining diverges are the documented `arm64`-platform
> floating-point `vcftools` cells (byte-exact on `amd64`/CI).

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

(The H1a change added `pipeline/stats` quantile/IQR/bootstrap-CI tests, all
green, and left every tool package untouched.)

## 2. full-validation (clean rerun)

| scale | parity | round-trip interop | perf bench | peak heap |
|---|---|---|---|---|
| smoke | 383 pass / 4 similar / 2 diverge / 0 error | **14/14 PASS** | 22 cells | 168 MB |
| small | 382 pass / 4 similar / 3 diverge / 0 error | **14/14 PASS** | 22 cells | 1.8 GB |
| medium | 383 pass / 4 similar / 2 diverge / 0 error | **14/14 PASS** | 22 cells | 11.5 GB |

**Every diverge is a documented `arm64`-platform FP cell — none are port defects:**

- `vcftools_site_mean_depth`, `vcftools_hardy` — arm64 `-nan` vs `nan` sign on a
  division-by-zero field (smoke + medium).
- small adds `vcftools_hap_r2` — a last-ULP `1.4e-17` vs `0` FP difference.

All three are byte-exact on `amd64`/CI. The 4 `SIMILAR` cells are the documented
floating-point cells (`bcftools call`; three `bedtools` similarity cells),
accepted by the gate.

So the genuine parity result is **all cells pass except the arm64-FP vcftools
cells**, with **round-trip interop 14/14** at every scale (BGZF/BAM/CRAM/VCF.gz/
BCF/FASTQ + BAI/CSI/TBI region queries, both directions).

## 3. Performance (`bench.md`, medium, reps=10) — robust stats

`ratio = ours / upstream` (`< 1.0` = we are faster). The headline is now the
**median ratio** with its **95% percentile-bootstrap CI** over the 10 reps (H1a);
the legacy min-over-reps point estimate is retained in `bench.md` below each
table. `±` is the inter-quartile range of the per-side wall time.

| cell | wall× (median) | 95% CI | note |
|---|---|---|---|
| `sickle_se` | 0.51 | [0.48, 0.56] | faster |
| `bed_intersect_pair` | 0.66 | [0.62, 0.70] | faster |
| `sam_sort` | 0.68 | [0.67, 0.69] | faster |
| `bed_intersect_self` | 0.69 | [0.67, 0.74] | faster |
| `sam_view_bam2bam` | 0.71 | [0.71, 0.72] | faster |
| `bed_coverage` | 0.75 | [0.74, 0.76] | faster |
| `sam_view_bam2cram` | 0.76 | [0.74, 0.80] | faster |
| `bed_genomecov` | 1.04 | [1.03, 1.08] | par |
| `bcf_stats` | 1.09 | [0.99, 1.12] | par |
| `bed_sort` | 1.09 | [0.98, 1.17] | par |
| `bcf_view` | 1.11 | [1.09, 1.13] | par |
| `sam_depth` | 1.12 | [1.06, 1.13] | par |
| `sam_stats` | 1.14 | [1.13, 1.16] | slower (CPU× 1.75 — more threads) |
| `sam_view_cram2bam` | 1.31 | [1.30, 1.33] | slower |
| `sam_flagstat` | 1.33 | [1.29, 1.34] | slower |
| `seqtk_seq` | 1.33 | [1.28, 1.37] | slower |
| `bcf_norm` | 1.46 | [1.42, 1.51] | slower |
| `bcf_query` | 1.57 | [1.52, 1.65] | slower |
| `bcf_call` | 1.72 | [1.70, 1.77] | **slow** |
| `bed_merge` | 2.15 | [1.85, 2.41] | **slow** |
| `sam_mpileup` | 2.40 | [2.37, 2.41] | **slow** (CPU× 3.39) |
| `bcf_isec` | 3.19 | [2.81, 3.40] | **slowest** (CPU× 5.34) |

I/O-bound conversions, `bedtools` intersect/coverage, and the `sickle` trimmer
are **faster** than upstream; the compute-heavy variant cells (`mpileup`,
`call`, `isec`) and `bed_merge` are **slower** and reported plainly.

**Caveats (read before citing):**

- Medium ran under `GOMEMLIMIT=8GiB` to fit the 12 GB VM, which adds GC pressure
  that **penalises the heavy cells** (`mpileup`/`call`/`isec` would look better
  unconstrained). The light/fast cells (well under the limit) are unaffected.
- The `RSS×` column reads the **orchestrator's** RSS, not per-subprocess, so it is
  a uniform 1.00 and is **not** a meaningful memory comparison.
- `arm64` container, not a fat node — `mpileup` etc. are absolute-time slow here.
