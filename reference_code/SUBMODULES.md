# Git Submodules — upstream sources of the ported tools

**Total submodules**: 11 (the upstream source of every tool this project ports;
the unported phase-1 survey tools were removed — see `README.md`).
Pinned versions/SHAs: `SUBMODULES.csv`.

| Tool | Language | Provides (our port) | Path |
|------|----------|---------------------|------|
| htslib | C | `bgzip`/`tabix`/`htsfile` + SAM/BAM/CRAM/VCF/BCF library reference | `reference_code/htslib` |
| htscodecs | C | CRAM codec (rANS/fqzcomp/tok3) byte-exactness test vectors | `reference_code/htscodecs` |
| samtools | C | `samtools` | `reference_code/samtools` |
| bcftools | C | `bcftools` | `reference_code/bcftools` |
| BEDTools | C++ | the `bed*` tool family | `reference_code/bedtools` |
| seqtk | C | `seqtk` | `reference_code/seqtk` |
| PRINSEQ | Perl | `prinseq` | `reference_code/prinseq` |
| Sickle | C | `sickle` | `reference_code/sickle` |
| Skewer | C++ | `skewer` | `reference_code/skewer` |
| fastp | C++ | `fastp` | `reference_code/fastp` |
| VCFtools | C++/Perl | `vcftools` | `reference_code/vcftools` |

These are read-only parity references — never modify them. Init on demand:
`git submodule update --init reference_code/<tool>` (use `--recursive` for `htslib`,
which carries its own `htscodecs` sub-submodule).
