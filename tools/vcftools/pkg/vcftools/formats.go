package vcftools

import (
	"fmt"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/vcf"
)

// Format conversion functions for vcftools

// output012Matrix outputs genotypes in 012 matrix format
// 0 = homozygous reference
// 1 = heterozygous
// 2 = homozygous alternate
// -1 = missing
func output012Matrix(variants []*vcf.Variant, header *vcf.Header, prefix string) error {
	// Output .012 file (genotype matrix)
	f012, err := iohelper.OpenWriter(prefix + ".012")
	if err != nil {
		return err
	}
	defer f012.Close()

	// Output .012.indv file (sample names)
	fIndv, err := iohelper.OpenWriter(prefix + ".012.indv")
	if err != nil {
		return err
	}
	defer fIndv.Close()

	// Output .012.pos file (variant positions)
	fPos, err := iohelper.OpenWriter(prefix + ".012.pos")
	if err != nil {
		return err
	}
	defer fPos.Close()

	// Write sample names
	for _, sample := range header.Samples {
		fmt.Fprintln(fIndv, sample)
	}

	// Restrict to biallelic sites — upstream emits a one-off warning then
	// drops multi-allelic loci because the 0/1/2 encoding has no room for
	// alt2/alt3.
	biallelic := make([]*vcf.Variant, 0, len(variants))
	for _, v := range variants {
		if len(v.Alt) == 1 {
			biallelic = append(biallelic, v)
		}
	}

	// Write variant positions (.012.pos)
	for _, v := range biallelic {
		fmt.Fprintf(fPos, "%s\t%d\n", v.Chrom, v.Pos)
	}

	// Write genotype matrix (one row per sample). Upstream prefixes each
	// row with the 0-based sample index, NOT the sample name. See
	// reference_code/vcftools/src/cpp/variant_file_format_convert.cpp
	// output_as_012_matrix().
	for sampleIdx := range header.Samples {
		fmt.Fprintf(f012, "%d", sampleIdx)

		for _, v := range biallelic {
			genotype := -1 // missing by default

			if sampleIdx < len(v.Samples) {
				sampleData := v.Samples[sampleIdx]
				gt, ok := sampleData.Data["GT"]

				if ok && !strings.Contains(gt, ".") {
					alleles := strings.FieldsFunc(gt, func(r rune) bool {
						return r == '/' || r == '|'
					})

					if len(alleles) == 2 {
						a1, a2 := alleles[0], alleles[1]
						if a1 == "0" && a2 == "0" {
							genotype = 0 // homozygous reference
						} else if a1 != a2 {
							genotype = 1 // heterozygous
						} else {
							genotype = 2 // homozygous alternate
						}
					}
				}
			}

			fmt.Fprintf(f012, "\t%d", genotype)
		}
		fmt.Fprintln(f012)
	}

	return nil
}

// outputPlink outputs genotypes in PLINK PED/MAP format
func outputPlink(variants []*vcf.Variant, header *vcf.Header, prefix string, chromMap map[string]int) error {
	// Output .ped file (genotypes)
	fPed, err := iohelper.OpenWriter(prefix + ".ped")
	if err != nil {
		return err
	}
	defer fPed.Close()

	// Output .map file (variant map)
	fMap, err := iohelper.OpenWriter(prefix + ".map")
	if err != nil {
		return err
	}
	defer fMap.Close()

	// Write MAP file (chr, variant ID, genetic distance, position)
	for _, v := range variants {
		chromNum := getChromNumber(v.Chrom, chromMap)
		varID := v.ID
		if varID == "" || varID == "." {
			varID = fmt.Sprintf("%s:%d", v.Chrom, v.Pos)
		}
		// Genetic distance is 0 (unknown)
		fmt.Fprintf(fMap, "%d\t%s\t0\t%d\n", chromNum, varID, v.Pos)
	}

	// Write PED file (family, individual, father, mother, sex, phenotype, genotypes)
	for _, sample := range header.Samples {
		// PLINK PED format: FID IID PAT MAT SEX PHENOTYPE genotypes...
		// We use sample name for both FID and IID, unknowns for others
		fmt.Fprintf(fPed, "%s\t%s\t0\t0\t0\t-9", sample, sample)

		for _, v := range variants {
			allele1, allele2 := "0", "0" // missing by default

			// Find this sample's genotype
			for _, sampleData := range v.Samples {
				if sampleData.Name == sample {
					gt, ok := sampleData.Data["GT"]

					if ok && !strings.Contains(gt, ".") {
						alleles := strings.FieldsFunc(gt, func(r rune) bool {
							return r == '/' || r == '|'
						})

						if len(alleles) == 2 {
							// Convert allele indices to actual alleles
							idx1, idx2 := alleles[0], alleles[1]
							allele1 = getAllele(v, idx1)
							allele2 = getAllele(v, idx2)
						}
					}
					break
				}
			}

			fmt.Fprintf(fPed, "\t%s\t%s", allele1, allele2)
		}
		fmt.Fprintln(fPed)
	}

	return nil
}

// outputPlinkTped outputs genotypes in PLINK transposed format
func outputPlinkTped(variants []*vcf.Variant, header *vcf.Header, prefix string, chromMap map[string]int) error {
	// Output .tped file (transposed genotypes)
	fTped, err := iohelper.OpenWriter(prefix + ".tped")
	if err != nil {
		return err
	}
	defer fTped.Close()

	// Output .tfam file (family information)
	fTfam, err := iohelper.OpenWriter(prefix + ".tfam")
	if err != nil {
		return err
	}
	defer fTfam.Close()

	// Write TFAM file (FID IID PAT MAT SEX PHENOTYPE)
	for _, sample := range header.Samples {
		fmt.Fprintf(fTfam, "%s\t%s\t0\t0\t0\t-9\n", sample, sample)
	}

	// Write TPED file (chr, variant ID, genetic distance, position, genotypes)
	for _, v := range variants {
		chromNum := getChromNumber(v.Chrom, chromMap)
		varID := v.ID
		if varID == "" || varID == "." {
			varID = fmt.Sprintf("%s:%d", v.Chrom, v.Pos)
		}

		fmt.Fprintf(fTped, "%d\t%s\t0\t%d", chromNum, varID, v.Pos)

		// Write genotypes for all samples
		for _, sample := range header.Samples {
			allele1, allele2 := "0", "0" // missing by default

			// Find this sample's genotype
			for _, sampleData := range v.Samples {
				if sampleData.Name == sample {
					gt, ok := sampleData.Data["GT"]

					if ok && !strings.Contains(gt, ".") {
						alleles := strings.FieldsFunc(gt, func(r rune) bool {
							return r == '/' || r == '|'
						})

						if len(alleles) == 2 {
							idx1, idx2 := alleles[0], alleles[1]
							allele1 = getAllele(v, idx1)
							allele2 = getAllele(v, idx2)
						}
					}
					break
				}
			}

			fmt.Fprintf(fTped, "\t%s\t%s", allele1, allele2)
		}
		fmt.Fprintln(fTped)
	}

	return nil
}

// getChromNumber converts chromosome name to number for PLINK
func getChromNumber(chrom string, chromMap map[string]int) int {
	if chromMap != nil {
		if num, ok := chromMap[chrom]; ok {
			return num
		}
	}

	// Try to extract number from chrom name
	chrom = strings.TrimPrefix(chrom, "chr")
	chrom = strings.TrimPrefix(chrom, "Chr")
	chrom = strings.TrimPrefix(chrom, "CHR")

	// Handle special cases
	switch chrom {
	case "X", "x":
		return 23
	case "Y", "y":
		return 24
	case "XY", "xy":
		return 25
	case "MT", "mt", "M", "m":
		return 26
	}

	// Try to parse as number
	var num int
	if _, err := fmt.Sscanf(chrom, "%d", &num); err == nil {
		return num
	}

	return 0 // unknown
}

// getAllele gets the actual allele string from variant
func getAllele(v *vcf.Variant, idx string) string {
	if idx == "." {
		return "0" // missing
	}

	var alleleIdx int
	if _, err := fmt.Sscanf(idx, "%d", &alleleIdx); err != nil {
		return "0"
	}

	if alleleIdx == 0 {
		return v.Ref
	}

	if alleleIdx > 0 && alleleIdx <= len(v.Alt) {
		return v.Alt[alleleIdx-1]
	}

	return "0" // unknown
}

// loadChromMap loads chromosome name to number mapping from file
func loadChromMap(filename string) (map[string]int, error) {
	if filename == "" {
		return nil, nil
	}

	f, err := iohelper.OpenReader(filename)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	chromMap := make(map[string]int)
	var chrom string
	var num int

	for {
		_, err := fmt.Fscanf(f, "%s\t%d\n", &chrom, &num)
		if err != nil {
			break
		}
		chromMap[chrom] = num
	}

	return chromMap, nil
}
