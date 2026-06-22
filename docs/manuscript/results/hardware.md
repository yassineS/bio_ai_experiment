# Hardware & toolchain (local heavy-tier run)

These results were produced on a developer laptop, **not** the fat node the
heavy-tier runbook targets. The byte-exact parity and performance numbers must
be read with the platform caveats in
[`LOCAL_RUN_NOTES.md`](LOCAL_RUN_NOTES.md) (macOS `libc++` vs `libstdc++`
`std::sort` tie-order; `arm64` vs `amd64` floating-point; 8 GB VM memory wall).

## Execution environment (where the harness actually ran)

The validation ran inside a **native `linux/arm64` Docker container** (Debian
bookworm) so that the upstream binaries link `libstdc++` — the C++ standard
library the ports were tuned against. Running the harness on the macOS host
directly is invalid for byte-exact parity (see the run notes).

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
