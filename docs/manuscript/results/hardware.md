# Hardware & toolchain (local heavy-tier run)

These results were produced on a developer laptop, **not** the fat node the
heavy-tier runbook targets. The byte-exact parity and performance numbers must
be read with the platform caveats in
[`LOCAL_RUN_NOTES.md`](LOCAL_RUN_NOTES.md) (macOS `libc++` vs `libstdc++`
`std::sort` tie-order; `arm64` vs `amd64` floating-point; 8 GB VM memory wall).

## Figure data provenance — which machine each tier came from (2026-06-29)

> **The performance figures mix two machines, by tier — read them accordingly.**
> Each *tier* is internally one machine (no cell within a tier mixes hosts), but
> the small/medium tiers and the large tier were measured on different hosts:

| Figure tier | Machine | When | Notes |
|---|---|---|---|
| **smoke / small / medium** | **darwin/arm64 host** (Apple M2, 16 GiB; macOS, native) | **2026-06-29** | Re-measured to reflect the CRAM `view -b -T` decode-RSS work; all 26 bench cells, `reps=5`, upstream binaries built **natively** on this host (Mach-O arm64). |
| **large** | linux/arm64 container (the pinned bench node below) | prior run | **Not** re-measured (too memory/disk-heavy for the laptop — the large fixtures are ~500 MB and the `mpileup`/`call` cells write ~17 GB temp). Left as the prior pinned-Linux data. |

Why mixed: the `darwin` `ru_maxrss` cross-platform **unit bug** (bytes vs KiB)
was fixed in `pipeline/bench/measure.go` (build-tagged `measure_linux.go` /
`measure_notlinux.go`), so peak-RSS is now reported in the same unit on both
hosts and the darwin re-measurement is valid. The large tier was left on the
prior linux/arm64 numbers rather than mixing hosts *within* a tier. The
small/medium → large discontinuity visible in some `fig_memory` bars (e.g. the
streaming-histogram and FASTQ cells, leaner at large) is partly the larger
inputs amortising the Go runtime's RSS floor and partly the host change — it is
**not** a code regression. The wall/CPU figures are robust to the host change
(both are arm64); the RSS figures carry the host caveat above.

> **CRAM `view` cells (`sam_view_cram2bam`, `sam_view_bam2cram`):** the
> small/medium bars now reflect the post-improvement decode path measured on
> this darwin host. NB the bench fixtures are synthetic (not the
> `up_chr20.cram` / `up_small.cram` corpus quoted in `METRICS.md`), so the
> absolute ratios differ from the realdata 1.85×/2.07× figures — the synthetic
> fixtures are smaller, so the Go runtime's fixed RSS floor weighs more heavily.

## Figure re-measurement host — smoke/small/medium tiers (2026-06-29)

The smoke/small/medium figure cells were re-measured directly on the **macOS
host** below (native, no container) after the `ru_maxrss` unit fix made
darwin RSS comparable to Linux. This is a *performance* re-measurement, not a
byte-exact parity run (parity stays on the Linux jobs, see the run notes).

| Property | Value |
|---|---|
| Model | Apple Mac14,2 (MacBook Air, Apple M2) |
| CPU | Apple M2, 8 cores |
| RAM | 16 GiB |
| OS | macOS 26.5.1 (Darwin kernel 25.5.0, `arm64`) |
| Go toolchain | `go1.26.4 darwin/arm64` (toolchain auto-upgrade over the `go.mod` 1.24.9 directive) |
| C/C++ compiler | Apple clang 21.0.0 |
| Reps | 5 per side, all 26 bench cells |
| CRAM mem-limit | `BIOAI_CRAM_MEMLIMIT_MIB` unset → default 36 MiB |

Upstream parity binaries used for these tiers were built **natively** on this
host (verified `Mach-O 64-bit executable arm64`):

| Tool | Version (native arm64 build) |
|---|---|
| samtools | 1.22.1-19-ge406d9e |
| bcftools | 1.23.1-73-ge0ec6ab0 |
| bedtools | v2.31.1 |
| seqtk | 1.5-r133 |
| sickle | 1.33 |

## Execution environment (where the linux/arm64 harness ran — large tier + parity)

The validation ran inside a **native `linux/arm64` Docker container** (Debian
bookworm) so that the upstream binaries link `libstdc++` — the C++ standard
library the ports were tuned against. Running the harness on the macOS host
directly is invalid for byte-exact parity (see the run notes). The **large**
figure tier and all byte-exact parity numbers come from this container.

| Property | Value |
|---|---|
| Container OS / kernel | Debian 12 (bookworm), Linux `6.12.76-linuxkit` |
| Architecture | `aarch64` (`linux/arm64`, native — not emulated) |
| vCPUs visible to container | 8 |
| Container/VM RAM | 7.7 GiB (Docker Desktop LinuxKit VM) |
| Go toolchain | `go1.24.9` (matches `go.mod`) |
| C/C++ compiler | gcc/g++ 12.2.0 |
| C library | GNU libc (glibc) 2.36 |
| C++ standard library | `libstdc++` (gcc 12) |

## Host (the physical machine under the VM)

| Property | Value |
|---|---|
| Model | Apple Mac14,2 (MacBook Air, Apple M2) |
| CPU | Apple M2, 8 cores |
| RAM | 16 GiB |
| OS | macOS 26.5.1 (Darwin kernel 25.5.0, `arm64` / T8112) |

> The Docker VM is capped at ~8 GB on this 16 GB host; the host was also running
> other containers, so the VM was **not** enlarged. This memory cap is the
> binding constraint on the medium/large tiers (see run notes §3).

## Upstream parity oracle (built from the vendored submodules)

All upstream binaries were built from `reference_code/` in the container.
**htslib was configured `--without-libdeflate`** (verified `HAVE_LIBDEFLATE`
undefined in `config.h`) so the BGZF byte stream matches the project's goldens.

| Tool | Version |
|---|---|
| samtools | 1.22.1-19-ge406d9eb |
| bcftools | 1.23.1-73-ge0ec6ab0 |
| htslib (bgzip/tabix/htsfile) | 1.23.1-32-gcdf22929 |
| bedtools | v2.31.1 |
| vcftools | 0.1.18 |
| fastp | 1.0.1 |
| seqtk / sickle / skewer / prinseq | from the pinned submodule commits |

> `mosdepth` is **not** built: upstream publishes a `linux/amd64`-only release
> binary (it is a Nim project, not built from source here), so every `mosdepth`
> cell SKIPs on `arm64` — mirroring the per-tool parity test's own guard.
