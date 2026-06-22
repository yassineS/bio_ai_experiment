# Local heavy-tier run — notes, findings, and honest limits

This records an attempt to run the heavy-tier validation
([`HEAVY_TIER_RUNBOOK.md`](HEAVY_TIER_RUNBOOK.md) / [`LOCAL_AGENT_PROMPT.md`](LOCAL_AGENT_PROMPT.md))
on a **developer laptop** (Apple M2, 16 GB, macOS) rather than the Linux fat
node the runbook targets. The honest headline:

- **What this box could validate cleanly:** the **small tier** in full (parity +
  bidirectional interop + performance), plus the **samtools** and **QC/htslib**
  groups at the **medium** tier.
- **What it could not:** the heavy **medium** groups (`bcftools`, `bedtools`) and
  the entire **large** tier — they exceed the 8 GB Docker-VM memory wall. The
  GIAB real-data step (H2a) and the K-run study (H3) were out of scope for this
  session.

Three platform findings below are worth carrying into the manuscript's
threats-to-validity, because they show the byte-exact parity bar is
**platform-sensitive** in ways the `linux/amd64` CI does not exercise.

Environment: [`hardware.md`](hardware.md). Reports: [`small_tier/`](small_tier/),
[`medium_tier/`](medium_tier/).

## 1. macOS is invalid for byte-exact parity (`libc++` vs `libstdc++`)

The first attempt ran the harness natively on macOS. It produced 6 `DIVERGE`
cells, **all** in `bedsort`/`bedcluster`/`bedsplit`, and **all** were the same
shape: records that tie on the sort key (identical `start`, or identical size)
emitted in a different order.

Root cause, confirmed two ways:

- The macOS upstream `bedtools` links **`libc++`** (`otool -L` shows
  `libc++.1.dylib`).
- The repo's `cppsort` package is, by its own source comment, *"a byte-faithful
  port of [**`libstdc++`**'s] introsort"* — it reproduces GNU `std::sort`'s
  (unstable) tie-ordering so it matches `bedtools` exactly.

`libc++` and `libstdc++` break `std::sort` ties differently, so on macOS our
`cppsort` (faithful to `libstdc++`) diverges from the `libc++`-built upstream.
On Linux (`libstdc++`, where CI runs) these cells are byte-exact. **Switching to
a `linux/arm64` container made all 6 vanish.**

**Implication:** drop-in byte-exact parity against `bedtools` (and any C++ tool
whose output order derives from `std::sort`) is only meaningful where upstream is
built with the *same* C++ standard library the port targets — i.e. `libstdc++`.
This should be stated as a portability caveat in the runbook and manuscript.

## 2. `arm64` vs `amd64` floating-point (NaN sign + FMA)

Under the (correct) `linux/arm64` container, all bed-tool divergences disappeared
and `ERROR` dropped to 0, but **3 `vcftools` cells** diverge — and only on
floating-point formatting:

| cell | ours | upstream |
|---|---|---|
| `vcftools_site_mean_depth` | `-nan` | `nan` |
| `vcftools_hardy` | `-nan` | `nan` |
| `vcftools_hap_r2` | `1.38778e-17` | `0` |

Mechanism (verified in-container):

- Plain C `0.0/0.0` on this `arm64` box yields `nan` (sign bit clear), so upstream
  prints `nan`. Go's NaN on this code path carries a **set** sign bit on `arm64`
  → our port prints `-nan`. The *same* Go expression yields an unsigned NaN on
  `amd64`, which is why these cells are byte-exact on CI.
- `1.38778e-17` vs `0` is sub-`2^-53` round-off, consistent with `arm64` FMA
  contraction of `a*b+c` differing from `amd64`.

These are **not** port defects and were deliberately **not** "fixed" (forcing an
unsigned NaN to match `arm64` upstream would risk breaking `amd64` parity and
violates the preserve-logic rule). They are byte-exact on the manuscript's
validated platform (`amd64`). Worth a sentence in threats-to-validity:
*statistical text output can be `arm64`/`amd64`-sensitive at the NaN-sign / last-ULP level.*

## 3. The medium/large memory wall (8 GB VM)

`pipeline/runner.RunEntry` runs matrix cells **sequentially** but buffers each
cell's **entire** ours-and-upstream stdout in memory (`runner.go:167-168`,
`bytes.Buffer` in `timedRun`) to byte-compare them; `Result` does not retain the
bytes, so peak RAM ≈ the **single heaviest cell**, not a cumulative sum.

On the medium fixtures that single-cell peak exceeds the 8 GB VM for the heavy
cells:

- `samtools` group (72 cells): peaks ~3.8 GB → **fits**, 72/72 PASS.
- `bcftools` group: OOM-killed on a heavy cell (`mpileup_heavy`/`call`).
- `bedtools` family: OOM-killed on a heavy cell (e.g. `genomecov -d` per-base
  over the medium genome → multi-GB ours+upstream buffers).
- QC/htslib group: fits (90 PASS / 2 arm64-FP DIVERGE / 10 mosdepth SKIP).

The large tier (≈19–30 GB fixtures + bigger per-cell buffers) is out of reach
here entirely. The Docker VM was **not** enlarged because the host was running
other containers and restarting Docker Desktop would have killed them; the
harness was **not** modified to stream-compare (a legitimate future improvement,
but out of scope for "run the validation"). **Medium/large belong on the
runbook's fat node (32–64 GB).**

## What passed — small tier (the canonical complete result)

[`small_tier/`](small_tier/) — full matrix, all tools, `-reps=10`:

- **Parity:** 400 cells — **382 PASS, 4 SIMILAR, 3 DIVERGE, 0 ERROR**, 11 SKIP.
  - The 3 DIVERGE are exactly the arm64-FP `vcftools` cells from §2.
  - The 4 SIMILAR are the documented floating-point cells (`bcftools call`; three
    `bedtools` similarity cells) — accepted by the runbook gate.
  - 0 ERROR (every upstream binary built and resolved).
- **Round-trip interop:** **14/14 PASS**, including bidirectional
  ours↔upstream interop for BGZF, BAM, CRAM, VCF.gz, BCF, FASTQ and BAI/CSI/TBI
  region-query index interop (`roundtrip.md`).
- **Performance:** `bench.md` / `bench.json`, `-reps=10`. Honest mixed picture
  (`ratio = ours/upstream`, wall-clock, min-over-reps):

  | cell | wall× | note |
  |---|---|---|
  | `sam_view_bam2cram` | 0.60 | faster |
  | `bed_intersect_pair` | 0.63 | faster |
  | `sam_sort` | 0.74 | faster |
  | `sam_view_bam2bam` | 0.75 | faster |
  | `bed_coverage` | 0.80 | faster |
  | `sam_depth` | 0.99 | par (but CPU× 1.96 — more threads) |
  | `sam_stats` | 1.18 | slower |
  | `bed_genomecov` | 1.26 | slower |
  | `sam_flagstat` | 1.32 | slower |
  | `sam_view_cram2bam` | 1.45 | slower |
  | `bed_merge` | 1.51 | slower |
  | `sam_mpileup` | 2.56 | **slow** (CPU× 3.64) — reported plainly |

  Two caveats on these numbers: (a) the RSS column reads the **orchestrator's**
  RSS, not the per-subprocess RSS, so `RSS×` is uniformly 1.00 and is **not**
  meaningful; (b) timings are the **min over reps**, which is precisely the H1a
  limitation — the manuscript wants **median + IQR + a ratio CI**. That stats
  upgrade (record raw per-rep samples; add `MedianIQR` + bootstrap/Hodges-Lehmann
  ratio-CI to `pipeline/stats`) is a separate code task, not done here.

## Not done (and why)

- **`bcftools`/`bedtools` medium, all of large** — 8 GB VM memory wall (§3). Fat node.
- **H2a GIAB real-data parity** — deferred by request ("large-tier only for now");
  needs multi-GB downloads.
- **H1a perf-stats upgrade** (median/IQR/ratio-CI) — code change, not attempted
  in this session.
- **H3 K-run reproducibility** — needs LLM-API budget; out of scope.

## Recommendations

1. Run medium/large (and H2a) on a **`linux/amd64`** fat node per the runbook;
   `amd64` also makes the §2 vcftools cells byte-exact and gives `mosdepth` a
   native binary.
2. Add the §1 (`libstdc++`-only) and §2 (`arm64`/`amd64` FP) caveats to the
   runbook's prerequisites and the manuscript's threats-to-validity.
3. Land the H1a perf-stats upgrade before citing performance ratios.
