# vcftools - VCF File Utilities

A partial Go reimplementation of vcftools for working with VCF (Variant Call
Format) files.

> **Scope:** this is **not** a drop-in replacement for upstream vcftools. It
> implements roughly 40 of vcftools' ~147 options — the commonly used filtering,
> per-site statistics, and a few format conversions (listed below). Options that
> are not implemented (LD analysis, Fst, ...) are **rejected with an error** rather
> than silently ignored. See [ROADMAP.md](ROADMAP.md) and
> [FEATURE_COMPARISON.md](FEATURE_COMPARISON.md).

## Overview

This port can:

- Filter VCF files (position, SNP ID, quality, allele frequency/count, variant type, genotype-level filters)
- Calculate per-site and per-individual statistics
- Convert to a few other formats (012 matrix, PLINK PED/MAP, PLINK TPED/TFAM)
- Handle both plain and gzipped VCF files

## Features

### Filtering Options

- **Position filtering**: by chromosome (`--chr`/`--not-chr`), position range (`--from-bp`/`--to-bp`), or position list (`--positions`/`--exclude-positions`)
- **SNP-ID filtering**: `--snp`, `--snps`, `--exclude`, `--exclude-snps`, `--thin`
- **Variant type filtering**: `--keep-only-indels`, `--remove-indels`, `--min-alleles`/`--max-alleles`
- **Quality filtering**: `--minQ`, `--remove-filtered-all`
- **Allele frequency / count filtering**: `--maf`/`--max-maf`, `--mac`/`--max-mac`
- **Genotype-level filtering**: `--max-missing`, `--min-meanDP`/`--max-meanDP`, `--minDP`/`--maxDP`, `--minGQ`
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

`--site-pi` uses the standard per-site formula `(n² − Σ cₐ²) / (n(n−1))` over
non-missing chromosomes; `--window-pi` reports the sum of per-site π over each
window. `--TajimaD` uses `D = (π − θ_W) / sqrt(e₁S + e₂S(S−1))` per window,
with the chromosome count taken as the modal value among the window's SNPs
(exact for complete data). (Earlier builds of this port reported a different,
incorrect quantity for `--site-pi` and silently ignored `--TajimaD`.)

### Format Conversion

- `--012` (0/1/2 genotype matrix), `--plink` (PED/MAP), `--plink-tped` (TPED/TFAM), with `--chrom-map`

### Not implemented

These options are recognised but **rejected with an error** (older builds
accepted them and produced nothing): all LD analysis
(`--geno-r2`, `--hap-r2`, ...).

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
