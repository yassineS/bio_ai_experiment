# Reference Code Directory

This directory vendors, as git submodules, the **upstream source of the tools this
project has ported** — kept as read-only references for behaviour/parity checking.
It is *not* a general tool collection: the earlier phase-1 survey submodules (BWA,
STAR, salmon, kraken2, …) that were never ported have been removed to keep the
checkout small and the `reference_code/` tree honest about scope.

## Vendored upstream sources (the ports' oracles)

| Submodule | Provides (our port) |
|---|---|
| `htslib` | the `bgzip`/`tabix`/`htsfile` utilities + the SAM/BAM/CRAM/VCF/BCF library reference |
| `htscodecs` | CRAM custom codec (rANS/fqzcomp/tok3) test vectors — byte-exactness oracle for `pkg/htsgo/cram/codec` |
| `samtools` | `samtools` (24 subcommands) |
| `bcftools` | `bcftools` (24 subcommands) |
| `bedtools` | the `bed*` tool family |
| `seqtk` | `seqtk` |
| `prinseq` | `prinseq` |
| `sickle` | `sickle` |
| `skewer` | `skewer` |
| `fastp` | `fastp` |
| `vcftools` | `vcftools` |

Exact pinned versions/SHAs are in `SUBMODULES.csv` (and `SUBMODULES.md`).

## Usage

These are read-only references — **never modify them**. Initialise one only when you
actually need it (e.g. to build the upstream binary for a live parity test):

```bash
git submodule update --init reference_code/<tool>
# htslib pulls its own htscodecs sub-submodule:
git submodule update --init --recursive reference_code/htslib
```

The parity harness (`pipeline/`) builds the upstream binaries from these submodules
on demand; CI's `upstream-parity` job builds `htslib`/`bcftools`/`samtools` and runs
the `*Upstream*` tests against them.

`patches/` holds the small local patches applied to an upstream source where needed
(e.g. a build fix); see `patches/README.md`.

> Historical note: `reference_code/` previously held ~50 submodules from the phase-1
> "top-200 tools" survey (`analysis/`). Per the project's current scope — finishing
> the already-ported tools, not broadening the set — the ~40 unported survey tools
> were removed. The survey data itself still lives under `analysis/`.
