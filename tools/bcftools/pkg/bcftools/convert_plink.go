package bcftools

// This file implements the PLINK exporters of the `bcftools convert`
// subcommand. Upstream's vcfconvert.c advertises a --plink/--tped/--bin
// family but leaves the option block commented out (lines ~1697-1699)
// with no implementation, so there is nothing to diff against upstream.
// The exporters below follow the well-documented PLINK1 file-format spec
// (https://www.cog-genomics.org/plink/1.9/formats):
//
//   - --plink     <prefix>  -> <prefix>.ped + <prefix>.map  (PLINK1 text)
//   - --tped      <prefix>  -> <prefix>.tped + <prefix>.tfam (transposed text)
//   - --plink-bed <prefix>  -> <prefix>.bed + <prefix>.bim + <prefix>.fam
//                              (PLINK1 BINARY, SNP-major)
//
// Conventions chosen (documented in tools/bcftools/README.md):
//
//   - Allele order. For the binary .bed/.bim we set A1 = ALT and A2 = REF,
//     which matches `plink --make-bed --keep-allele-order` applied to a VCF
//     (PLINK1's default A1 is the minor allele; --keep-allele-order pins it
//     to the VCF ALT). The 2-bit codes are therefore:
//       00 = hom-A2 (REF/REF)   <- two A2 (major) alleles
//       01 = missing
//       10 = het   (REF/ALT)
//       11 = hom-A1 (ALT/ALT)   <- two A1 (minor) alleles
//     i.e. the on-disk count is the number of A1 (ALT) alleles, exactly as
//     PLINK1 stores it. .bim columns are "CHROM SNP_ID 0 BP A1 A2" with
//     A1=ALT, A2=REF.
//
//   - Chromosome codes. Numeric contigs 1..22 pass through; X->23, Y->24,
//     XY->25, MT/M->26 (with an optional "chr"/"CHR" prefix accepted and
//     stripped). Anything else is written verbatim, which PLINK1.9 accepts
//     under --allow-extra-chr.
//
//   - Multi-allelic sites. PLINK1 is strictly biallelic. Records with more
//     than one ALT allele are skipped and counted; the first such record
//     emits a one-line warning to stderr (mirroring the hap exporters).
//     Records with no ALT allele are likewise skipped.

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// PlinkConvertOptions controls the PLINK exporters of `bcftools convert`
// (--plink, --tped and --plink-bed). The sample/region/expression filters
// mirror the other convert modes.
type PlinkConvertOptions struct {
	// Prefix is the value passed to --plink / --tped / --plink-bed. A bare
	// prefix expands to the canonical per-format suffixes; an explicit
	// comma-separated file list is also accepted (2 names for --plink and
	// --tped, 3 for --plink-bed).
	Prefix string

	// Samples / SamplesFile / ForceSamples restrict the per-sample columns
	// of the input VCF before export.
	Samples      []string
	SamplesFile  string
	ForceSamples bool

	// Regions / RegionsFile / Targets / TargetsFile apply CHROM[:beg-end]
	// post-filters to the input VCF.
	Regions     []string
	RegionsFile string
	Targets     []string
	TargetsFile string

	// IncludeExpr / ExcludeExpr are the standard -i / -e filter
	// expressions applied to the input VCF.
	IncludeExpr string
	ExcludeExpr string
}

// asHapOptions adapts a PlinkConvertOptions into the HapConvertOptions
// shape consumed by loadFilteredVCF (the shared input-loading path used by
// the IMPUTE2 hap exporters).
func (o *PlinkConvertOptions) asHapOptions() *HapConvertOptions {
	return &HapConvertOptions{
		Samples:      o.Samples,
		SamplesFile:  o.SamplesFile,
		ForceSamples: o.ForceSamples,
		Regions:      o.Regions,
		RegionsFile:  o.RegionsFile,
		Targets:      o.Targets,
		TargetsFile:  o.TargetsFile,
		IncludeExpr:  o.IncludeExpr,
		ExcludeExpr:  o.ExcludeExpr,
	}
}

// plinkOutputNames resolves a prefix (or explicit comma-separated name
// list) into the per-format output file names. suffixes are the canonical
// extensions applied to a bare prefix; their count also sets how many
// comma-separated names an explicit list must supply.
func plinkOutputNames(prefix string, suffixes []string) ([]string, error) {
	if strings.IndexByte(prefix, ',') < 0 {
		names := make([]string, len(suffixes))
		for i, s := range suffixes {
			names[i] = prefix + s
		}
		return names, nil
	}
	parts := strings.Split(prefix, ",")
	if len(parts) != len(suffixes) {
		return nil, fmt.Errorf("error parsing PLINK filenames: %s (expected %d comma-separated names)", prefix, len(suffixes))
	}
	return parts, nil
}

// plinkChrom maps a VCF contig name to the PLINK1 chromosome code. Numeric
// 1..22 pass through; X/Y/XY/MT/M map to 23/24/25/26; an optional
// "chr"/"CHR" prefix is stripped first. Any other name is returned
// verbatim (accepted by PLINK1.9 under --allow-extra-chr).
func plinkChrom(chrom string) string {
	c := chrom
	if len(c) >= 3 && strings.EqualFold(c[:3], "chr") {
		c = c[3:]
	}
	switch strings.ToUpper(c) {
	case "X":
		return "23"
	case "Y":
		return "24"
	case "XY":
		return "25"
	case "MT", "M":
		return "26"
	}
	// Numeric 1..22 (and any other plain positive integer) pass through.
	if n, err := strconv.Atoi(c); err == nil && n > 0 {
		return c
	}
	// Non-numeric, non-special: keep the original name verbatim.
	return chrom
}

// plinkVariantID returns the SNP id for a variant: its VCF ID, or
// "CHROM:POS" when the ID is missing (".").
func plinkVariantID(v *vcf.Variant) string {
	if v.ID != "" && v.ID != "." {
		return v.ID
	}
	return fmt.Sprintf("%s:%d", v.Chrom, v.Pos)
}

// plinkSkips tracks the two reasons a record is dropped by the PLINK
// exporters: no ALT allele, or more than one ALT (non-biallelic).
type plinkSkips struct {
	noAlt        int
	nonBiallelic int
}

// plinkExportable filters variants to the biallelic, ALT-bearing subset the
// PLINK formats accept, warning once on the first multi-allelic record. The
// returned slice preserves input order.
func plinkExportable(variants []*vcf.Variant, stderr io.Writer) ([]*vcf.Variant, plinkSkips) {
	var skips plinkSkips
	warnedMulti := false
	kept := make([]*vcf.Variant, 0, len(variants))
	for _, v := range variants {
		if len(v.Alt) == 0 || (len(v.Alt) == 1 && (v.Alt[0] == "." || v.Alt[0] == "")) {
			skips.noAlt++
			continue
		}
		if len(v.Alt) > 1 {
			if !warnedMulti && stderr != nil {
				fmt.Fprintln(stderr, "Warning: non-biallelic records are skipped (PLINK is biallelic). Consider 'bcftools norm -m-' first.")
				warnedMulti = true
			}
			skips.nonBiallelic++
			continue
		}
		kept = append(kept, v)
	}
	return kept, skips
}

// reportPlinkSkips prints the per-conversion summary line to stderr.
func reportPlinkSkips(stderr io.Writer, nok int, s plinkSkips) {
	if stderr == nil {
		return
	}
	fmt.Fprintf(stderr, "%d records written, %d skipped: %d/%d no-ALT/non-biallelic\n",
		nok, s.noAlt+s.nonBiallelic, s.noAlt, s.nonBiallelic)
}

// pedAlleles renders one sample's two .ped/.tped allele letters for a
// biallelic variant: the actual REF/ALT letters keyed by the GT allele
// indices, or "0 0" for a missing genotype. Diploid calls use both
// alleles; a haploid call is rendered as a homozygote.
func pedAlleles(gt, ref, alt string) string {
	alleles := splitGT(gt)
	letter := func(idx int) string {
		if idx == 0 {
			return ref
		}
		return alt
	}
	switch len(alleles) {
	case 0:
		return "0 0"
	case 1:
		a0, m0 := parseAllele(alleles[0])
		if m0 {
			return "0 0"
		}
		l := letter(a0)
		return l + " " + l
	default:
		a0, m0 := parseAllele(alleles[0])
		a1, m1 := parseAllele(alleles[1])
		// A partially-missing genotype (e.g. "0/.") is treated as fully
		// missing, matching PLINK's all-or-nothing missing convention.
		if m0 || m1 {
			return "0 0"
		}
		return letter(a0) + " " + letter(a1)
	}
}

// bedCode returns the 2-bit PLINK1 .bed code for one sample's genotype
// under the A1=ALT, A2=REF convention. The byte stores the number of A1
// (ALT) alleles: 00 = two A2 (hom-REF), 11 = two A1 (hom-ALT), 10 = het,
// 01 = missing.
func bedCode(gt string) byte {
	alleles := splitGT(gt)
	switch len(alleles) {
	case 0:
		return 0b01
	case 1:
		a0, m0 := parseAllele(alleles[0])
		if m0 {
			return 0b01
		}
		if a0 == 0 {
			return 0b00 // hom A2 (REF)
		}
		return 0b11 // hom A1 (ALT)
	default:
		a0, m0 := parseAllele(alleles[0])
		a1, m1 := parseAllele(alleles[1])
		if m0 || m1 {
			return 0b01
		}
		n1 := 0
		if a0 != 0 {
			n1++
		}
		if a1 != 0 {
			n1++
		}
		switch n1 {
		case 0:
			return 0b00 // hom A2 (REF/REF)
		case 1:
			return 0b10 // het
		default:
			return 0b11 // hom A1 (ALT/ALT)
		}
	}
}

// gtFor returns sample i's GT string, or "" when the variant carries no GT
// FORMAT field for that sample.
func gtFor(v *vcf.Variant, i int) string {
	if i >= len(v.Samples) {
		return ""
	}
	return v.Samples[i].Data["GT"]
}

// requireGT confirms the variant declares a FORMAT/GT field.
func requireGT(v *vcf.Variant) error {
	if formatIndex(v.Format, "GT") < 0 {
		return fmt.Errorf("FORMAT/GT tag not present at %s:%d", v.Chrom, v.Pos)
	}
	return nil
}

// requireGTAll confirms every kept variant declares a FORMAT/GT field.
func requireGTAll(variants []*vcf.Variant) error {
	for _, v := range variants {
		if err := requireGT(v); err != nil {
			return err
		}
	}
	return nil
}

// writePlinkMap writes the .map file: one "CHROM SNP_ID 0 BP" line per
// variant.
func writePlinkMap(name string, variants []*vcf.Variant) error {
	f, err := os.Create(name)
	if err != nil {
		return err
	}
	defer f.Close()
	bw := bufio.NewWriter(f)
	for _, v := range variants {
		if _, err := fmt.Fprintf(bw, "%s\t%s\t0\t%d\n", plinkChrom(v.Chrom), plinkVariantID(v), v.Pos); err != nil {
			return err
		}
	}
	if err := bw.Flush(); err != nil {
		return err
	}
	return f.Close()
}

// writePlinkFam writes a PLINK .fam / .tfam file: six columns per sample
// "FID IID PAT MAT SEX PHENO" with FID=IID=sample name, PAT=MAT=0,
// SEX=0 (unknown), PHENO=-9 (missing).
func writePlinkFam(name string, samples []string) error {
	f, err := os.Create(name)
	if err != nil {
		return err
	}
	defer f.Close()
	bw := bufio.NewWriter(f)
	for _, s := range samples {
		if _, err := fmt.Fprintf(bw, "%s %s 0 0 0 -9\n", s, s); err != nil {
			return err
		}
	}
	if err := bw.Flush(); err != nil {
		return err
	}
	return f.Close()
}

// VCFToPlink implements `bcftools convert --plink`: it reads the VCF/BCF at
// path, applies the filters in opts, and writes a PLINK1 text fileset
// (<prefix>.ped + <prefix>.map). It returns the number of variants written.
func VCFToPlink(path string, opts PlinkConvertOptions, stderr io.Writer) (int, error) {
	hdr, variants, err := loadFilteredVCF(path, opts.asHapOptions())
	if err != nil {
		return 0, err
	}
	names, err := plinkOutputNames(opts.Prefix, []string{".ped", ".map"})
	if err != nil {
		return 0, err
	}
	pedName, mapName := names[0], names[1]

	kept, skips := plinkExportable(variants, stderr)
	if err := requireGTAll(kept); err != nil {
		return 0, err
	}

	if mapName != "" {
		if err := writePlinkMap(mapName, kept); err != nil {
			return 0, err
		}
	}
	if pedName != "" {
		if err := writePlinkPed(pedName, hdr.Samples, kept); err != nil {
			return 0, err
		}
	}
	reportPlinkSkips(stderr, len(kept), skips)
	return len(kept), nil
}

// writePlinkPed writes the .ped file: one line per sample, six mandatory
// columns then the two allele letters of every variant.
func writePlinkPed(name string, samples []string, variants []*vcf.Variant) error {
	f, err := os.Create(name)
	if err != nil {
		return err
	}
	defer f.Close()
	bw := bufio.NewWriter(f)
	for si, s := range samples {
		// FID IID PAT MAT SEX PHENO with FID=IID=name, SEX=0, PHENO=-9.
		if _, err := fmt.Fprintf(bw, "%s %s 0 0 0 -9", s, s); err != nil {
			return err
		}
		for _, v := range variants {
			if err := bw.WriteByte(' '); err != nil {
				return err
			}
			if _, err := bw.WriteString(pedAlleles(gtFor(v, si), v.Ref, v.Alt[0])); err != nil {
				return err
			}
		}
		if err := bw.WriteByte('\n'); err != nil {
			return err
		}
	}
	if err := bw.Flush(); err != nil {
		return err
	}
	return f.Close()
}

// VCFToPlinkTransposed implements `bcftools convert --tped`: it reads the
// VCF/BCF at path, applies the filters in opts, and writes a transposed
// PLINK text fileset (<prefix>.tped + <prefix>.tfam). It returns the number
// of variants written.
func VCFToPlinkTransposed(path string, opts PlinkConvertOptions, stderr io.Writer) (int, error) {
	hdr, variants, err := loadFilteredVCF(path, opts.asHapOptions())
	if err != nil {
		return 0, err
	}
	names, err := plinkOutputNames(opts.Prefix, []string{".tped", ".tfam"})
	if err != nil {
		return 0, err
	}
	tpedName, tfamName := names[0], names[1]

	kept, skips := plinkExportable(variants, stderr)
	if err := requireGTAll(kept); err != nil {
		return 0, err
	}

	if tfamName != "" {
		if err := writePlinkFam(tfamName, hdr.Samples); err != nil {
			return 0, err
		}
	}
	if tpedName != "" {
		if err := writePlinkTped(tpedName, hdr.Samples, kept); err != nil {
			return 0, err
		}
	}
	reportPlinkSkips(stderr, len(kept), skips)
	return len(kept), nil
}

// writePlinkTped writes the .tped file: one line per variant,
// "CHROM SNP_ID 0 BP" then the two allele letters of every sample.
func writePlinkTped(name string, samples []string, variants []*vcf.Variant) error {
	f, err := os.Create(name)
	if err != nil {
		return err
	}
	defer f.Close()
	bw := bufio.NewWriter(f)
	for _, v := range variants {
		if _, err := fmt.Fprintf(bw, "%s %s 0 %d", plinkChrom(v.Chrom), plinkVariantID(v), v.Pos); err != nil {
			return err
		}
		for si := range samples {
			if err := bw.WriteByte(' '); err != nil {
				return err
			}
			if _, err := bw.WriteString(pedAlleles(gtFor(v, si), v.Ref, v.Alt[0])); err != nil {
				return err
			}
		}
		if err := bw.WriteByte('\n'); err != nil {
			return err
		}
	}
	if err := bw.Flush(); err != nil {
		return err
	}
	return f.Close()
}

// plinkBedMagic is the 3-byte PLINK1 .bed header: 0x6c 0x1b then 0x01 for
// SNP-major (variant-major) ordering.
var plinkBedMagic = []byte{0x6c, 0x1b, 0x01}

// VCFToPlinkBinary implements `bcftools convert --plink-bed`: it reads the
// VCF/BCF at path, applies the filters in opts, and writes a PLINK1 BINARY
// fileset (<prefix>.bed + <prefix>.bim + <prefix>.fam). It returns the
// number of variants written.
func VCFToPlinkBinary(path string, opts PlinkConvertOptions, stderr io.Writer) (int, error) {
	hdr, variants, err := loadFilteredVCF(path, opts.asHapOptions())
	if err != nil {
		return 0, err
	}
	names, err := plinkOutputNames(opts.Prefix, []string{".bed", ".bim", ".fam"})
	if err != nil {
		return 0, err
	}
	bedName, bimName, famName := names[0], names[1], names[2]

	kept, skips := plinkExportable(variants, stderr)
	if err := requireGTAll(kept); err != nil {
		return 0, err
	}

	if famName != "" {
		if err := writePlinkFam(famName, hdr.Samples); err != nil {
			return 0, err
		}
	}
	if bimName != "" {
		if err := writePlinkBim(bimName, kept); err != nil {
			return 0, err
		}
	}
	if bedName != "" {
		if err := writePlinkBed(bedName, len(hdr.Samples), kept); err != nil {
			return 0, err
		}
	}
	reportPlinkSkips(stderr, len(kept), skips)
	return len(kept), nil
}

// writePlinkBim writes the .bim file: one "CHROM SNP_ID 0 BP A1 A2" line
// per variant with A1=ALT (minor) and A2=REF (major).
func writePlinkBim(name string, variants []*vcf.Variant) error {
	f, err := os.Create(name)
	if err != nil {
		return err
	}
	defer f.Close()
	bw := bufio.NewWriter(f)
	for _, v := range variants {
		// A1=ALT, A2=REF.
		if _, err := fmt.Fprintf(bw, "%s\t%s\t0\t%d\t%s\t%s\n",
			plinkChrom(v.Chrom), plinkVariantID(v), v.Pos, v.Alt[0], v.Ref); err != nil {
			return err
		}
	}
	if err := bw.Flush(); err != nil {
		return err
	}
	return f.Close()
}

// writePlinkBed writes the SNP-major PLINK1 .bed body: the 3-byte magic
// header then, for each variant, ceil(nsamples/4) bytes packing 2 bits per
// sample little-endian within the byte (sample 0 -> bits 0-1).
func writePlinkBed(name string, nsamples int, variants []*vcf.Variant) error {
	f, err := os.Create(name)
	if err != nil {
		return err
	}
	defer f.Close()
	bw := bufio.NewWriter(f)
	if _, err := bw.Write(plinkBedMagic); err != nil {
		return err
	}
	bytesPerVariant := (nsamples + 3) / 4
	block := make([]byte, bytesPerVariant)
	for _, v := range variants {
		for i := range block {
			block[i] = 0
		}
		for si := 0; si < nsamples; si++ {
			code := bedCode(gtFor(v, si))
			block[si/4] |= code << uint((si%4)*2)
		}
		if _, err := bw.Write(block); err != nil {
			return err
		}
	}
	if err := bw.Flush(); err != nil {
		return err
	}
	return f.Close()
}
