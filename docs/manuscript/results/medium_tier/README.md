# Medium tier — partial (8 GB VM memory wall)

The medium tier was run on the local laptop container described in
[`../hardware.md`](../hardware.md). It is **partial**: the 8 GB Docker VM cannot
hold the heaviest medium cells, so two of the four tool groups were OOM-killed.
This is a resource limit of the box, **not** a parity failure — see
[`../LOCAL_RUN_NOTES.md`](../LOCAL_RUN_NOTES.md) §3 for the mechanism.

## What ran

The matrix was split into tool groups (`-tools=…`) run sequentially so each
would fit; the harness runs cells sequentially and buffers each cell's full
ours+upstream stdout in RAM, so peak memory ≈ the single heaviest cell.

| Group | Result | Notes |
|---|---|---|
| `samtools` ([samtools/](samtools/)) | **72 / 72 PASS**, 0 DIVERGE, 0 ERROR | byte-exact at medium scale |
| QC + htslib ([qc/](qc/)) | **90 PASS, 2 DIVERGE, 10 SKIP, 0 ERROR** | DIVERGE = arm64-FP vcftools (below); the 10 SKIP are `mosdepth` (run separately, below) |
| `mosdepth` ([mosdepth/](mosdepth/)) | **10 / 10 PASS**, 0 DIVERGE, 0 ERROR | run in a `linux/amd64` container (below) |
| `bcftools` | **OOM-killed** | a single heavy cell (`mpileup_heavy`/`call`) exceeds the 8 GB VM |
| `bedtools` family | **OOM-killed** | a single heavy cell (e.g. `genomecov -d` per-base over the medium genome) exceeds the 8 GB VM |

The `samtools` group was run with `-reps=1 -skip-bench` (parity only); the QC
group likewise. No medium-tier performance bench is reported here — the
performance numbers are taken from the small tier (see `../small_tier/`), which
completed in full.

### mosdepth (unskipped via a `linux/amd64` container)

`mosdepth` ships only a `linux/amd64` release binary, and the harness SKIPs it on
`arm64`. To actually run it, the v0.3.14 binary was pulled and the `mosdepth`
cells were re-run in an **emulated `linux/amd64` container** (`MOSDEPTH_BIN` set;
our pure-Go `mosdepth` rebuilt for `amd64`). Result: **10/10 PASS, byte-exact**,
including the `--by` BED/window/threshold/region-dist/summary paths. The cells'
timing ratios in [mosdepth/report.md](mosdepth/report.md) are **emulation
artifacts** (both sides run under `amd64` emulation on the M2) and are **not**
performance results. The round-trip stage of that run reports failures only
because the upstream non-mosdepth binaries on disk are `arm64` and cannot exec in
the `amd64` container — irrelevant to mosdepth parity; the canonical interop
result is the `arm64` 14/14 PASS in `../small_tier/roundtrip.md`.

## The 2 QC-group DIVERGEs are arm64 floating-point artifacts

Both are `vcftools` statistics that print a NaN, and both differ only in the
**sign of the NaN**, which is `arm64`-specific (they are byte-exact on `amd64`,
which is why CI is green):

- `vcftools_site_mean_depth` — ours `-nan`, upstream `nan`
- `vcftools_hardy` — ours `-nan`, upstream `nan`

See the run notes §2 for the mechanism (Go's NaN carries a set sign bit on
`arm64` for this code path; the same expression yields an unsigned NaN on
`amd64`). These are **not** port defects and were **not** masked.

## To complete the medium (and large) tier

Run on the runbook's intended **`linux/amd64` fat node** (32–64 GB RAM,
≥30 GB scratch), where (a) the heavy cells fit, (b) `amd64` floating-point makes
the vcftools cells byte-exact, and (c) `mosdepth` has a native upstream binary:

```bash
go run ./pipeline/cmd/full-validation -scales=medium,large -reps=10 \
  -out="$PWD/docs/manuscript/results/large_tier"
```
