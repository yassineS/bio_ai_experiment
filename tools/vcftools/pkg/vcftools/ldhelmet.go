package vcftools

import (
	"bufio"
	"fmt"
	"strconv"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/vcf"
)

// LDhelmet output writer, ported from upstream vcftools
// src/cpp/variant_file_format_convert.cpp:982-1157 (output_as_LDhelmet).
//
// CLI flag is `--ldhelmet`; upstream registers it in parameters.cpp:275 as:
//
//	output_as_ldhelmet = true; phased_only = true; remove_indels = true;
//	num_outputs++;
//
// So `--ldhelmet` implies both --phased and --remove-indels. Like --ldhat /
// --ldhat-geno, it also requires exactly one chromosome (--chr) — enforced at
// parameters.cpp:717.
//
// Output layout:
//
//	<prefix>.ldhelmet.snps  — FASTA-style haplotype rows (2 per kept sample),
//	                          using the actual nucleotide (REF or ALT[i]) for
//	                          each kept site or 'N' for missing.
//	<prefix>.ldhelmet.pos   — one integer POS per kept site, in record order.
//
// Header of .snps is implicit (FASTA-style, one ">indv-k" per haplotype).
//
// Differences from --ldhat / --ldhat-geno:
//
//   - LDhelmet emits the nucleotide character(s) rather than '0'/'1'/'2'/'?',
//     so multi-allelic loci can contribute meaningful output (alleles[geno]
//     indexes into REF + ALT).
//   - LDhelmet writes positions as bare integers, not "%.4f" kb. Upstream
//     does set ios::fixed/precision(4) on the stream but writes an int via
//     `<<`, which produces plain integer formatting (see line 1119-1127).
//   - No biallelic-only filter; upstream just looks up alleles[geno] for
//     whatever geno is in the GT slot. (--remove-indels is on by default so
//     tri-allelic SNPs would still be in scope; the existing
//     params.MaxAlleles guard in passFilters defaults to 0 = unbounded so
//     such sites pass.)
//
// The upstream implementation streams via temporary files (one per
// haplotype) to keep memory bounded. We buffer in memory instead — same
// per-haplotype string concatenation, just without the temp-file dance.
// Output bytes are identical.

// ldhelmetRunner buffers per-site haplotype contributions for --ldhelmet
// and emits the paired (.ldhelmet.snps, .ldhelmet.pos) bundle at flush.
type ldhelmetRunner struct {
	prefix  string
	samples []string
	// haplotypes[2*i+k] is the kept-site allele string for sample i,
	// haplotype k. Each appended entry is a single allele substring (one
	// nucleotide for SNPs, longer for non-SNV alleles — but --ldhelmet
	// implies --remove-indels so indels are dropped upstream of us).
	// Missing/out-of-range alleles append "N".
	haplotypes [][]string
	positions  []int
}

// newLDhelmetRunner constructs the runner. samples is the *filtered* sample
// list (so --indv / --keep / --remove are already honoured upstream of us).
func newLDhelmetRunner(prefix string, samples []string) *ldhelmetRunner {
	r := &ldhelmetRunner{
		prefix:     prefix,
		samples:    append([]string(nil), samples...),
		haplotypes: make([][]string, 2*len(samples)),
	}
	return r
}

// addVariant buffers one variant's contribution to the LDhelmet output.
func (r *ldhelmetRunner) addVariant(v *vcf.Variant) {
	if r == nil {
		return
	}
	// alleles[0] = REF, alleles[1..] = ALT[0..]. Matches upstream's
	// get_alleles_vector (entry_getters.cpp:147-153).
	alleles := make([]string, 0, len(v.Alt)+1)
	alleles = append(alleles, v.Ref)
	alleles = append(alleles, v.Alt...)

	r.positions = append(r.positions, v.Pos)
	for i := range r.samples {
		a1, a2, _ := parseGTForLDhat(getGT(v, i))
		r.haplotypes[2*i] = append(r.haplotypes[2*i], ldhelmetAllele(a1, alleles))
		r.haplotypes[2*i+1] = append(r.haplotypes[2*i+1], ldhelmetAllele(a2, alleles))
	}
}

// ldhelmetAllele looks up the actual allele string by index. Mirrors
// upstream's switch at variant_file_format_convert.cpp:1100-1114:
//
//	if (geno >= 0 && include_genotype[ui]==true) -> alleles[geno]
//	else -> "N"
//
// We map our parseGTForLDhat negative codes (-1 = ".", -2 = haploid second
// allele absent) to "N" too. Numeric indices out of range likewise yield
// "N" (parseGTForLDhat preserves the int; if the GT references an allele
// that doesn't exist we cannot emit a nucleotide).
func ldhelmetAllele(allele int, alleles []string) string {
	if allele < 0 || allele >= len(alleles) {
		return "N"
	}
	return alleles[allele]
}

// close writes the .ldhelmet.snps and .ldhelmet.pos files. Safe to call on a
// nil receiver.
func (r *ldhelmetRunner) close() error {
	if r == nil {
		return nil
	}

	posPath := r.prefix + ".ldhelmet.pos"
	snpsPath := r.prefix + ".ldhelmet.snps"

	posW, err := iohelper.OpenWriter(posPath)
	if err != nil {
		return fmt.Errorf("opening %s: %w", posPath, err)
	}
	defer posW.Close()
	posBuf := bufio.NewWriter(posW)
	for _, p := range r.positions {
		posBuf.WriteString(strconv.Itoa(p))
		posBuf.WriteByte('\n')
	}
	if err := posBuf.Flush(); err != nil {
		return fmt.Errorf("writing %s: %w", posPath, err)
	}

	snpsW, err := iohelper.OpenWriter(snpsPath)
	if err != nil {
		return fmt.Errorf("opening %s: %w", snpsPath, err)
	}
	defer snpsW.Close()
	snpsBuf := bufio.NewWriter(snpsW)
	for i, name := range r.samples {
		for k := 0; k < 2; k++ {
			fmt.Fprintf(snpsBuf, ">%s-%d\n", name, k)
			for _, a := range r.haplotypes[2*i+k] {
				snpsBuf.WriteString(a)
			}
			snpsBuf.WriteByte('\n')
		}
	}
	if err := snpsBuf.Flush(); err != nil {
		return fmt.Errorf("writing %s: %w", snpsPath, err)
	}
	return nil
}
