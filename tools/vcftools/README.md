# vcftools - VCF File Utilities

A partial Go reimplementation of vcftools for working with VCF (Variant Call
Format) files.

> **Scope:** this is **not** a drop-in replacement for upstream vcftools. It
> implements roughly 70 of vcftools' ~147 options — the commonly used filtering,
> per-site statistics, basic and inter-chromosomal LD analysis, chi-square LD
> tests, Weir & Cockerham Fst, KING / Yang relatedness, runs of homozygosity,
> phased-block reporting, INFO-tag extraction, and a few format conversions
> (listed below). Many other options (`--mendel`, `--ldhat`, `--IMPUTE`, full
> `--diff` extensions …) are **rejected with an error** rather than silently
> ignored. See [ROADMAP.md](ROADMAP.md) and
> [FEATURE_COMPARISON.md](FEATURE_COMPARISON.md).

## Overview

This port can:

- Filter VCF files (position, SNP ID, quality, allele frequency/count, variant type, genotype-level filters)
- Calculate per-site and per-individual statistics
- Convert to a few other formats (012 matrix, PLINK PED/MAP, PLINK TPED/TFAM)
- Handle both plain and gzipped VCF files

## Features

### Filtering Options

- **Position filtering**: by chromosome (`--chr`/`--not-chr`), position range (`--from-bp`/`--to-bp`), position list (`--positions`/`--exclude-positions`), or BED intervals (`--bed`/`--exclude-bed`)
- **SNP-ID filtering**: `--snp`, `--snps`, `--exclude`, `--exclude-snps`, `--thin`
- **Variant type filtering**: `--keep-only-indels`, `--remove-indels`, `--min-alleles`/`--max-alleles`
- **Quality filtering**: `--minQ`, `--remove-filtered-all`
- **Allele frequency / count filtering**: `--maf`/`--max-maf`, `--mac`/`--max-mac`
- **Genotype-level filtering**: `--max-missing`, `--max-missing-count`, `--hwe`, `--min-meanDP`/`--max-meanDP`, `--minDP`/`--maxDP`, `--minGQ`
- **Sample filtering**: `--indv`, `--remove-indv`, `--keep`, `--remove`

### Statistics Output

- Allele frequency and counts (`--freq`, `--counts`, `--freq2`, `--counts2`)
- Depth: per individual (`--depth` → `.idepth`), per genotype (`--geno-depth` → `.gdepth`), per site summed (`--site-depth`) and mean (`--site-mean-depth`)
- Site quality (`--site-quality`)
- Missingness: per individual (`--missing-indv`) and per site (`--missing-site`)
- Hardy-Weinberg equilibrium (`--hardy`)
- Heterozygosity / F per individual (`--het`); singletons (`--singletons`)
- Transition/transversion ratios: `--TsTv-summary`, `--TsTv N`, `--TsTv-by-count`, `--TsTv-by-qual` (→ `.TsTv.qual`: Ts/Tv counts and ratios cumulative below and at-or-above each distinct QUAL threshold)
- Nucleotide diversity: per site (`--site-pi`) and windowed (`--window-pi`, `--window-pi-step`)
- Tajima's D in non-overlapping windows (`--TajimaD N` → `.Tajima.D`)
- Weir & Cockerham 1984 Fst per biallelic SNP (`--weir-fst-pop` → `.weir.fst`) plus optional windowed output (`--fst-window-size`, `--fst-window-step` → `.windowed.weir.fst`); the per-site mean and weighted Fst summary is printed to stderr.
- Indel length histogram (`--hist-indel-len` → `.indel.hist`)
- FILTER summary (`--FILTER-summary`); SNP density (`--SNPdensity N`)
- Site-pass / site-fail trace (`--kept-sites` → `.kept.sites`,
  `--removed-sites` → `.removed.sites`): two-column `CHROM\tPOS` TSVs
  listing which sites pass or fail every filter, in input order.

`--site-pi` uses the standard per-site formula `(n² − Σ cₐ²) / (n(n−1))` over
non-missing chromosomes; `--window-pi` reports the sum of per-site π over each
window. `--TajimaD` uses `D = (π − θ_W) / sqrt(e₁S + e₂S(S−1))` per window,
with the chromosome count taken as the modal value among the window's SNPs
(exact for complete data). (Earlier builds of this port reported a different,
incorrect quantity for `--site-pi` and silently ignored `--TajimaD`.)

### Format Conversion

- `--012` (0/1/2 genotype matrix), `--plink` (PED/MAP), `--plink-tped` (TPED/TFAM), with `--chrom-map`
- `--BEAGLE-GL` → `<prefix>.BEAGLE.GL` — log10 genotype likelihoods derived
  from FORMAT/PL. Biallelic SNPs only; sites without a PL field are skipped
  with a one-time stderr warning.
- `--BEAGLE-PL` → `<prefix>.BEAGLE.PL` — same selection as `--BEAGLE-GL` but
  with the raw Phred PL triplets.

### VCF Comparison (--diff family)

Compare two VCFs site-by-site and per-individual. The second VCF is loaded
fully into memory; the first file is streamed and respects any normal
filtering flags (`--chr`, `--bed`, `--minQ`, ...). The second file is
*not* filtered, matching upstream vcftools.

Flags:

- `--diff FILE` — second VCF to compare against (`.gz` auto-detected).
- `--diff-site` → `<prefix>.diff.sites_in_files` with columns
  `CHROM POS1 POS2 IN_FILE REF1 REF2 ALT1 ALT2` where `IN_FILE ∈ {1, 2, B}`.
- `--diff-indv` → `<prefix>.diff.indv_in_files` listing every sample tagged
  with `1`, `2`, or `B`.
- `--diff-site-discordance` → `<prefix>.diff.sites` with columns
  `CHROM POS FILES MATCHING_ALLELES N_COMMON_CALLED N_DISCORD DISCORDANCE`,
  matching upstream's 7-column layout. Sites unique to one file are
  still listed (FILES = `1` or `2`) with zero counts and `-nan`
  discordance, mirroring upstream's libstdc++ output.
- `--diff-indv-discordance` → `<prefix>.diff.indv` with columns
  `INDV N_COMMON_CALLED N_DISCORD DISCORDANCE` over the *union* of
  file-1 and (post-map) file-2 samples in alphabetical order.
  Samples that appear in only one file get `0 0 -nan`.
- `--diff-indv-map FILE` — two-column whitespace-separated table that
  renames file-2 sample IDs before matching against file-1 (mirrors
  upstream `variant_file_diff.cpp:11-34`). Lines starting with `#` and
  blank lines are skipped. Composes with all other `--diff-*` outputs.
- `--diff-discordance-matrix` → `<prefix>.diff.discordance_matrix` —
  the 4x4 genotype-by-genotype contingency table (`0/0`, `0/1`, `1/1`,
  `./.`) for biallelic loci with matching ALT alleles in both files
  (mirrors upstream `variant_file_diff.cpp:944`). Byte-for-byte parity
  vs upstream is covered in `diff_parity_test.go`.

Discordance compares unphased, sorted allele indices restricted to REF/first
ALT; samples with multi-allelic calls at a given site are treated as missing
for that site (mirroring upstream's default behaviour). `--gzdiff` is wired
as a plain alias for `--diff` (the `iohelper` reader auto-sniffs gzip from
the magic bytes, so there's no distinction at the parser layer).

### Haplotype counts (`--hapcount`)

`--hapcount BED` writes `<prefix>.hapcount` — per-BED-bin tallies of unique
kept-individual haplotypes (columns:
`#CHROM BIN_START BIN_END N_SNP N_UNIQ_HAPS N_GROUPS {MULTIPLICITY:FREQ}...`).
Implies `--phased`. Diploid-only per-site. Bins must be non-overlapping.

This port fixes three upstream defects in `output_haplotype_count`
(`variant_file_output.cpp:1169-1401`):

1. `prev_bin_idx` shift on within-chromosome bin transitions silently
   overwrote one bin's counts with the next bin's.
2. End-of-stream read-after-free silently dropped (or zeroed) the last
   chromosome's rows.
3. The BED first line was unconditionally skipped; header-less BEDs
   silently lost their first bin.

See `docs/UPSTREAM_BUGS.md#fix-on-port-resolved` for the full writeup.

### Misc

`--temp DIR` is accepted for CLI parity with upstream (parameters.cpp:341)
but has no effect — this port doesn't spill to disk. A stderr note is
printed when the flag is supplied.

### Linkage Disequilibrium

Genotype-based (`--geno-r2`) and haplotype-based (`--hap-r2`) pairwise LD
between sites on the same chromosome. Window size, minimum distance, and an
r² threshold can all be configured.

Flags:

- `--geno-r2` → `<prefix>.geno.ld` (columns: `CHR POS1 POS2 N_INDV R^2`).
  Uses per-sample diploid allele counts `g_i ∈ {0,1,2}` over individuals
  non-missing at both sites; multi-allelic sites use only the first ALT.
- `--hap-r2` → `<prefix>.hap.ld` (columns: `CHR POS1 POS2 N_CHR R^2 D Dprime`).
  Requires phased GTs (`a|b`); unphased samples are skipped.
- `--geno-r2-positions FILE`, `--hap-r2-positions FILE` — only emit pairs
  where at least one endpoint is listed in `FILE` (tab-separated `chrom pos`).
- `--ld-window INT`, `--ld-window-bp INT` — maximum SNP and bp distance
  between paired sites (default: unbounded).
- `--ld-window-min INT`, `--ld-window-bp-min INT` — minimum SNP and bp
  distance (default 0).
- `--min-r2 FLOAT` — drop pairs below this r².

Examples:

```bash
# Genotype-based LD within 1 kb windows.
vcftools --vcf input.vcf --geno-r2 --ld-window-bp 1000 --out ld

# Haplotype-based LD (phased data) keeping only strong LD.
vcftools --vcf phased.vcf --hap-r2 --min-r2 0.5 --out hapld
```

Upstream byte-for-byte parity hasn't been validated yet; see the follow-up
issue in `ROADMAP.md`.

#### Inter-chromosomal LD and chi-square

- `--interchrom-geno-r2` → `<prefix>.interchrom.geno.ld`
  (columns: `CHR1 POS1 CHR2 POS2 N_INDV R^2`). Buffers all kept variants in
  memory and emits every cross-chromosome pair.
- `--interchrom-hap-r2` → `<prefix>.interchrom.hap.ld`
  (columns: `CHR1 POS1 CHR2 POS2 N_CHR R^2 D Dprime`).
- `--geno-chisq` → `<prefix>.geno.chisq`
  (columns: `CHR1 POS1 CHR2 POS2 N_INDV CHI^2 DF P-VALUE`). Pearson
  chi-square test of association on the 3×3 contingency table of diploid
  ALT counts, with degrees of freedom restricted to non-empty rows/columns.
  P-values via the regularised upper incomplete gamma function.

### Relatedness and runs of homozygosity

- `--relatedness` → `<prefix>.relatedness`
  (columns: `INDV1 INDV2 RELATEDNESS_AJK`). Yang et al. 2010 unadjusted
  per-pair A_jk estimator averaged over biallelic SNPs.
- `--relatedness2` → `<prefix>.relatedness2`
  (columns: `INDV1 INDV2 N_AaAa N_AAaa N1_Aa N2_Aa RELATEDNESS_PHI`).
  KING-robust kinship coefficient (Manichaikul et al. 2010).
- `--LROH` → `<prefix>.LROH`
  (columns: `CHROM AUTO_START AUTO_END N_VARIANTS INDV`). Contiguous runs
  of homozygous (`0/0` or `1/1`) diploid genotypes per individual. Default
  minimum run length is 10 variants; override with `--LROH-min-variants N`.
- `--phased-blocks` → `<prefix>.blocks`
  (columns: `CHROM BLOCK_START BLOCK_END N_VARIANTS INDV`). Per-individual
  contiguous runs of phased (`a|b`) diploid genotypes; runs shorter than
  two variants are not emitted.

### FILTER / INFO selection

- `--remove-filtered NAME[,NAME...]` — drop sites listing any of the named
  FILTERs. Composes with `--remove-filtered-all`.
- `--keep-filtered NAME[,NAME...]` — keep only sites listing at least one
  named FILTER (sites with `PASS` or `.` are dropped).
- `--keep-INFO TAG` — SITE FILTER. Keep only sites where the named
  Flag-type INFO tag is present. Errors at runtime if the tag is not
  declared `Type=Flag` in the header. Repeatable; multiple tags
  compose via OR (any present keeps the site). Mirrors upstream
  `entry_filters.cpp:1033`.
- `--recode-INFO TAG` — recode-column selector. Restrict the INFO
  column in `.recode.vcf` output to the listed tag(s). Repeatable.
  Mirrors upstream `parameters.cpp:319`.
- `--remove-INFO TAG` — strip a tag from the INFO column of recoded
  output. Repeatable. NOTE: upstream defines this flag as a SITE
  FILTER (drop sites where the named Flag IS present); the port
  currently implements the column-stripping semantic and tracks the
  divergence in `docs/PARITY_ROADMAP.md`.
- `--get-INFO TAG[,TAG...]` → `<prefix>.INFO`
  (columns: `CHROM POS REF ALT TAG1 TAG2 ...`). Missing values emit `.`.

### Principal Component Analysis

- `--pca` → `<prefix>.pca` (rows: `INDV` header, an `EIGENVALUE` row,
  then one row per kept sample with that sample's loadings on each
  principal component). Implements upstream's
  `variant_file_output.cpp:4871-5042` recipe: build the per-individual
  centred / variance-normalised genotype matrix `M`, form the
  `(1/n)·M·Mᵀ` GRM, eigendecompose, sort by `|eigenvalue|` descending.
  Eigendecomposition uses `gonum.org/v1/gonum/mat`'s symmetric solver
  (the project's second sanctioned third-party dep after the CRAM
  codec carveout — see `CLAUDE.md`).
- `--pca-no-norm` — implies `--pca`, but skips the per-SNP `1/√(p(1-p))`
  normalisation (still mean-centres). Mirrors upstream
  `parameters.cpp:298`.
- `--pca-snp-loadings INT` → `<prefix>.pca.loadings` — per-site
  projections of raw genotypes onto the first `INT` principal
  components (`GAMMA_0..GAMMA_{INT-1}`). Mirrors upstream
  `variant_file_output.cpp:5044-5246`. Does NOT imply `--pca` (a
  separate output, independent of `.pca` per upstream
  `vcftools.cpp:110-111`).

Eigenvector signs are arbitrary in both LAPACK and gonum; the port
applies a "first nonzero component positive" canonicalisation so the
output is deterministic across LAPACK builds. Parity tests use a
per-column sign-tolerant comparison.

Wave 19 also fixes a latent upstream bug — `output_PCA` reads past
the end of the per-individual `M[i]` vectors when any kept individual
has a missing genotype; the port drops such sites entirely. This is
a deliberate divergence from both upstream's buggy "skip the push,
advance the index" behaviour AND Patterson, Price & Reich 2006's
mean-imputation recipe — the port chooses the most conservative
correct behaviour (drop the site) so the GRM is built only from
fully-observed sites. See `docs/UPSTREAM_BUGS.md` for the upstream
bug write-up; if mean-imputation is preferred for a downstream
analysis, regenerate the input with imputation applied first.

### Not implemented

These options are recognised but **rejected with an error** (older builds
accepted them and produced nothing): the BCF binary I/O family
(`--bcf`, `--diff-bcf`, `--recode-bcf`, `--contigs`; blocked on htsgo
PR-G — see `docs/HTSGO_ROADMAP.md`), and a long tail of less-used
upstream options. `--mendel` and `--diff-switch-error` landed in
wave 5; `--hwe` and `--max-missing-count` landed in wave 8; the PCA
family landed in wave 19.
See `ROADMAP.md`.

### Format Support

- VCF v4.0, v4.1, and v4.2
- Automatic gzip handling
- Stdin/stdout support

## Installation

```bash
cd tools/vcftools
go build ./cmd/vcftools
```

## Usage

### Basic Usage

```bash
# Show help
./vcftools --help

# Read from VCF file
./vcftools --vcf input.vcf [options]

# Read from gzipped VCF
./vcftools --gzvcf input.vcf.gz [options]

# Read from stdin
cat input.vcf | ./vcftools --stdin [options]
```

### Filtering Examples

```bash
# Filter by chromosome
./vcftools --vcf input.vcf --chr chr1 --recode --out chr1_only

# Filter by position range
./vcftools --vcf input.vcf --chr chr1 --from-bp 1000000 --to-bp 2000000 --recode --out region

# Filter by quality
./vcftools --vcf input.vcf --minQ 30 --recode --out high_quality

# Filter by minor allele frequency
./vcftools --vcf input.vcf --maf 0.05 --recode --out common_variants

# Remove indels
./vcftools --vcf input.vcf --remove-indels --recode --recode-INFO-all --out snps_only

# Keep only indels
./vcftools --vcf input.vcf --keep-only-indels --recode --out indels_only

# Filter by missing data (keep sites with <10% missing)
./vcftools --vcf input.vcf --max-missing 0.9 --recode --out low_missing

# Combine multiple filters
./vcftools --vcf input.vcf \
  --chr chr1 \
  --minQ 30 \
  --maf 0.05 \
  --max-missing 0.9 \
  --remove-indels \
  --recode --recode-INFO-all \
  --out filtered
```

### Statistics Examples

```bash
# Calculate allele frequencies
./vcftools --vcf input.vcf --freq --out output

# Calculate allele counts
./vcftools --vcf input.vcf --counts --out output

# Calculate mean depth per site
./vcftools --vcf input.vcf --site-mean-depth --out output

# Calculate site depth
./vcftools --vcf input.vcf --site-depth --out output

# Calculate site quality
./vcftools --vcf input.vcf --site-quality --out output

# Calculate site missingness
./vcftools --vcf input.vcf --missing-site --out output

# Calculate individual missingness
./vcftools --vcf input.vcf --missing-indv --out output

# Test for Hardy-Weinberg equilibrium
./vcftools --vcf input.vcf --hardy --out output

# Calculate Ts/Tv ratio
./vcftools --vcf input.vcf --TsTv-summary --out output

# Calculate Ts/Tv in 10kb bins
./vcftools --vcf input.vcf --TsTv 10000 --out output

# Calculate nucleotide diversity
./vcftools --vcf input.vcf --site-pi --out output

# Calculate multiple statistics
./vcftools --vcf input.vcf \
  --freq \
  --site-mean-depth \
  --missing-site \
  --hardy \
  --TsTv-summary \
  --out stats
```

### Sample Filtering Examples

```bash
# Keep specific individuals
./vcftools --vcf input.vcf --indv sample1 --indv sample2 --recode --out subset

# Remove specific individuals
./vcftools --vcf input.vcf --remove-indv sample1 --recode --out filtered

# Keep individuals from file (one per line)
./vcftools --vcf input.vcf --keep samples_to_keep.txt --recode --out subset

# Remove individuals from file
./vcftools --vcf input.vcf --remove samples_to_remove.txt --recode --out filtered
```

### Position Filtering Examples

```bash
# Filter by positions in a file
# File format: tab-delimited with chromosome and position
# Example file content:
#   chr1  1000
#   chr1  2000
#   chr2  3000

./vcftools --vcf input.vcf --positions positions.txt --recode --out subset

# Exclude positions from file
./vcftools --vcf input.vcf --exclude-positions exclude.txt --recode --out filtered

# Keep only sites inside intervals in a BED file (0-based half-open).
./vcftools --vcf input.vcf --bed regions.bed --recode --out in_regions

# Drop sites inside the regions file instead.
./vcftools --vcf input.vcf --exclude-bed exclude.bed --recode --out outside_regions
```

### VCF Comparison Example

```bash
# Compare callsets — emits .diff.sites_in_files, .diff.sites,
# .diff.indv_in_files, and .diff.indv next to <prefix>.
./vcftools --vcf callerA.vcf \
  --diff callerB.vcf.gz \
  --diff-site --diff-site-discordance \
  --diff-indv --diff-indv-discordance \
  --out compare
```

### BEAGLE Output Example

```bash
# Write BEAGLE-format genotype likelihoods (skips non-SNP sites and sites
# without a FORMAT/PL field, with a one-time stderr warning).
./vcftools --vcf input.vcf --BEAGLE-GL --out for_beagle
./vcftools --vcf input.vcf --BEAGLE-PL --out for_beagle
```

## Output Files

### Recode Output

- `<prefix>.recode.vcf` - Filtered VCF file

### Statistics Outputs

- `<prefix>.frq` - Allele frequencies
- `<prefix>.frq.count` - Allele counts
- `<prefix>.ldepth.mean` - Mean depth per site
- `<prefix>.ldepth` - Sum depth per site
- `<prefix>.lqual` - Quality per site
- `<prefix>.lmiss` - Site missingness
- `<prefix>.imiss` - Individual missingness
- `<prefix>.hwe` - Hardy-Weinberg equilibrium test results
- `<prefix>.TsTv.summary` - Ts/Tv ratio summary
- `<prefix>.TsTv` - Ts/Tv ratios by bin
- `<prefix>.sites.pi` - Nucleotide diversity per site

## Performance

vcftools (Go) provides comparable performance to the original C++/Perl implementation with the following advantages:

- Single binary with no external dependencies
- Consistent behavior across platforms
- Automatic gzip handling
- Memory-efficient streaming processing

## Differences from Original vcftools

### Feature Coverage

This implementation provides approximately **80% of commonly used vcftools features**, covering:

- ✅ Core filtering options (position, quality, allele frequency, variant type)
- ✅ Essential statistics (frequency, depth, missingness, HWE, Ts/Tv, pi)
- ✅ VCF recoding with filtering
- ✅ Sample filtering and management

See [FEATURE_COMPARISON.md](FEATURE_COMPARISON.md) for a detailed feature-by-feature comparison.

### Major Missing Features

The following commonly-used features are not yet implemented:

**High Priority:**

- ❌ Linkage disequilibrium (LD) calculations (`--geno-r2`, `--hap-r2`)
- ❌ PLINK format conversion (`--plink`, `--plink-tped`)

**Medium Priority:**

- ❌ VCF comparison/diff operations (`--diff`, `--diff-site-discordance`)
- ❌ Windowed statistics (`--window-pi`, `--TajimaD`)
- ❌ Additional format conversions (BEAGLE, IMPUTE, LDhat, 012 matrix)

**Low Priority:**

- ❌ Advanced filtering (BED files, SNP IDs, masks)
- ❌ Relatedness and PCA analysis
- ❌ BCF format support

These features may be added in future versions. Contributions are welcome!

### Improvements

- Simpler command-line interface
- Automatic handling of gzipped files
- No need for separate Perl and C++ tools
- Single binary distribution

## Feature Coverage Summary

This Go implementation includes **~35 of the original ~147 vcftools options** (24% total coverage), but covers approximately **80% of commonly used features** in real-world analyses. The implementation prioritizes:

1. **Core filtering operations** - Essential for data quality control
2. **Basic statistics** - Most frequently requested outputs
3. **VCF manipulation** - Recoding and sample management

For a complete feature-by-feature comparison with the original vcftools, see [FEATURE_COMPARISON.md](FEATURE_COMPARISON.md).

## Implementation Details

### VCF Format Support

- VCF v4.0, v4.1, v4.2
- Handles standard VCF fields: CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, samples
- Supports both phased (|) and unphased (/) genotypes

### Filtering Logic

Filters are applied in the following order:

1. Position filters (chromosome, position range, position list)
2. Variant type filters (indels, allele count)
3. Quality filters (minQ, filter status)
4. Allele frequency filters (MAF, MAC)
5. Genotype filters (missing data, depth)

All filters are combined with AND logic (all must pass).

### Statistics Calculations

**Allele Frequency**: Calculated as count of allele / total non-missing alleles

**Minor Allele Frequency (MAF)**: Frequency of the less common allele

**Hardy-Weinberg Equilibrium**: Chi-square test comparing observed vs expected genotype frequencies

**Ts/Tv Ratio**: Ratio of transitions (A↔G, C↔T) to transversions (all other changes)

**Nucleotide Diversity (π)**: Average pairwise difference between sequences

## Testing

```bash
# Run tests
go test ./pkg/vcftools

# Run tests with coverage
go test -cover ./pkg/vcftools
```

## Examples

See the `examples/` directory for sample VCF files and usage examples.

## Contributing

Contributions are welcome! Areas for improvement:

- Additional statistics calculations
- LD analysis
- Format conversion
- Performance optimizations
- Additional filtering options

## License

This project is licensed under the Apache License 2.0 - see the LICENSE file for details.

## Citation

If you use this tool in your research, please cite the original vcftools paper:

> Danecek P, Auton A, Abecasis G, Albers CA, Banks E, DePristo MA, Handsaker RE, Lunter G, Marth GT, Sherry ST, McVean G, Durbin R; 1000 Genomes Project Analysis Group. The variant call format and VCFtools. Bioinformatics. 2011 Aug 1;27(15):2156-8. doi: 10.1093/bioinformatics/btr330

## Acknowledgments

Based on the original vcftools by:

- Adam Auton
- Petr Danecek
- Anthony Marcketta

## Contact

For questions or issues, please open an issue on the GitHub repository.
