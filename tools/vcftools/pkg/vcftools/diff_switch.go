// Phase-switch error analysis for the vcftools --diff family.
//
// CLI flag: --diff-switch-error.
//
// Upstream source: reference_code/vcftools/src/cpp/variant_file_diff.cpp:1207-
// 1507 (variant_file::output_switch_error). Registered in
// parameters.cpp:208 as:
//
//	else if (in_str == "--diff-switch-error") {
//	    diff_switch_error = true; num_outputs++;
//	}
//
// Output files (mirroring upstream byte layout):
//
//	<prefix>.diff.switch        per-event log; header
//	                            "CHROM\tPOS_START\tPOS_END\tINDV". One row per
//	                            detected switch event (file2_hap1 disagrees
//	                            with both file1 haplotypes when comparing the
//	                            transition from the previous phased het to the
//	                            current site).
//
//	<prefix>.diff.indv.switch   per-individual summary; header
//	                            "INDV\tN_COMMON_PHASED_HET\tN_SWITCH\tSWITCH".
//	                            SWITCH = N_SWITCH / N_COMMON_PHASED_HET when
//	                            the denominator is > 0; otherwise 0.
//
// Both files participate in the regular --diff sample-intersection logic.
// When --diff-indv-map is supplied the mapping is honoured before forming the
// intersection (upstream return_indv_union does the same).
//
// Algorithm (per individual; mirrors lines 1396-1469 of upstream):
//
//   1. Iterate file-1 variants in their VCF-file order; for every common
//      (file-1, file-2) site that is heterozygous *and* phased in both files,
//      treat it as one "common phased het site".
//   2. The first such site initialises prev_geno/prev_pos for the individual
//      but emits no event.
//   3. For each subsequent common phased het site, build
//
//        file1_hap1 = (prev_geno_file1.first,  curr_geno_file1.first)
//        file1_hap2 = (prev_geno_file1.second, curr_geno_file1.second)
//        file2_hap1 = (prev_geno_file2.first,  curr_geno_file2.first)
//
//      If file2_hap1 matches neither file1 haplotype the phases disagree and
//      we emit one row in <prefix>.diff.switch.
//
//   4. Increment N_SWITCH for the individual on every emitted event;
//      increment N_COMMON_PHASED_HET on every common phased het site.
//
// On stream EOF, write the per-individual summary in commonPairs order, so
// the output is stable (file-1 sample order, mirroring upstream's
// combined_individuals iteration which is map-ordered).

package vcftools

import (
	"bufio"
	"fmt"
	"io"
	"strconv"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// switchEventOut bundles the <prefix>.diff.switch writer with its file handle.
type switchEventOut struct {
	f io.WriteCloser
	w *bufio.Writer
}

// switchRunner accumulates phase-switch state per individual.
type switchRunner struct {
	prefix string
	pairs  []commonPair

	// Per-individual previous phased het state. Indices line up with
	// `pairs`.
	prevAllele1F1 []string // first-haplotype allele in file1 at prev het
	prevAllele2F1 []string // second-haplotype allele in file1 at prev het
	prevAllele1F2 []string // first-haplotype allele in file2 at prev het
	prevAllele2F2 []string // second-haplotype allele in file2 at prev het
	prevChromF1   []string
	prevPosF1     []int
	prevChromF2   []string
	prevPosF2     []int
	hasPrev       []bool

	// Per-individual aggregates.
	nCommonPhasedHet []int
	nSwitch          []int

	out *switchEventOut
}

// newSwitchRunner opens <prefix>.diff.switch and prepares the per-individual
// counters. pairs is the file-1-ordered intersection of (file-1 sample,
// file-2 sample) names — same source as the rest of the diff runners.
func newSwitchRunner(prefix string, pairs []commonPair) (*switchRunner, error) {
	if len(pairs) == 0 {
		return nil, fmt.Errorf("--diff-switch-error: no overlapping individuals can be found")
	}
	path := prefix + ".diff.switch"
	f, err := iohelper.OpenWriter(path)
	if err != nil {
		return nil, fmt.Errorf("opening %s: %w", path, err)
	}
	w := bufio.NewWriter(f)
	if _, err := w.WriteString("CHROM\tPOS_START\tPOS_END\tINDV\n"); err != nil {
		f.Close()
		return nil, fmt.Errorf("writing header to %s: %w", path, err)
	}
	r := &switchRunner{
		prefix:           prefix,
		pairs:            pairs,
		prevAllele1F1:    make([]string, len(pairs)),
		prevAllele2F1:    make([]string, len(pairs)),
		prevAllele1F2:    make([]string, len(pairs)),
		prevAllele2F2:    make([]string, len(pairs)),
		prevChromF1:      make([]string, len(pairs)),
		prevPosF1:        make([]int, len(pairs)),
		prevChromF2:      make([]string, len(pairs)),
		prevPosF2:        make([]int, len(pairs)),
		hasPrev:          make([]bool, len(pairs)),
		nCommonPhasedHet: make([]int, len(pairs)),
		nSwitch:          make([]int, len(pairs)),
		out:              &switchEventOut{f: f, w: w},
	}
	return r, nil
}

// addVariant processes a file-1 variant against its matching file-2 record (if
// any). v is the post-filter file-1 variant; rec is the file-2 record at the
// same (chrom,pos), or nil if file-2 has no entry there.
func (r *switchRunner) addVariant(v *vcf.Variant, rec *diffRecord) {
	if r == nil || rec == nil {
		return
	}
	for i, pair := range r.pairs {
		gt1, ok1 := findSampleGT(v, pair.f1Name)
		gt2, ok2 := rec.genotypes[pair.f2Name]
		if !ok1 || !ok2 {
			continue
		}
		a11, a12, phased1 := splitPhasedAlleles(gt1)
		a21, a22, phased2 := splitPhasedAlleles(gt2)
		// Both must be non-missing diploid calls.
		if a11 == "" || a12 == "" || a21 == "" || a22 == "" {
			continue
		}
		if a11 == "." || a12 == "." || a21 == "." || a22 == "." {
			continue
		}
		// Matching genotype check (set equality with both orderings).
		matched :=
			(a11 == a21 && a12 == a22) ||
				(a11 == a22 && a12 == a21)
		if !matched {
			continue
		}
		// Heterozygote check: the two alleles in file1 must differ.
		if a11 == a12 {
			continue
		}
		// Both must be phased ('|') — upstream only enters the inner
		// branch when phase1 == '|' && phase2 == '|'.
		if !phased1 || !phased2 {
			continue
		}

		// We are in a "common phased het" site.
		r.nCommonPhasedHet[i]++

		if r.hasPrev[i] {
			// file1_hap1 = (prevA1F1, curA1F1)
			// file1_hap2 = (prevA2F1, curA2F1)
			// file2_hap1 = (prevA1F2, curA1F2)
			f2h1a := r.prevAllele1F2[i]
			f2h1b := a21
			f1h1a := r.prevAllele1F1[i]
			f1h1b := a11
			f1h2a := r.prevAllele2F1[i]
			f1h2b := a12
			eq1 := f2h1a == f1h1a && f2h1b == f1h1b
			eq2 := f2h1a == f1h2a && f2h1b == f1h2b
			if !eq1 && !eq2 {
				// Switch error event.
				r.nSwitch[i]++
				// POS_START / POS_END: same-chromosome adjacency
				// only. Upstream lines 1453-1459: if the two prev
				// chromosomes match, pick the smaller of the two
				// prev positions for POS_START. We require both
				// prev sites and the current site to share a
				// chromosome to emit (file2 is keyed at the same
				// (chrom,pos) as file1, so prev_chrom_f1 ==
				// prev_chrom_f2 by construction).
				if r.prevChromF1[i] == r.prevChromF2[i] {
					startPos := r.prevPosF1[i]
					if r.prevPosF2[i] < startPos {
						startPos = r.prevPosF2[i]
					}
					chrom := r.prevChromF1[i]
					// POS_END is the *current* file-1 position
					// (upstream uses POS1).
					fmt.Fprintf(r.out.w, "%s\t%d\t%d\t%s\n",
						chrom, startPos, v.Pos, pair.f1Name)
				}
			}
		}

		// Update prev to current.
		r.prevAllele1F1[i] = a11
		r.prevAllele2F1[i] = a12
		r.prevAllele1F2[i] = a21
		r.prevAllele2F2[i] = a22
		r.prevChromF1[i] = v.Chrom
		r.prevPosF1[i] = v.Pos
		r.prevChromF2[i] = v.Chrom // file2 is keyed by same (chrom,pos)
		r.prevPosF2[i] = v.Pos
		r.hasPrev[i] = true
	}
}

// close flushes <prefix>.diff.switch and writes <prefix>.diff.indv.switch.
func (r *switchRunner) close() error {
	if r == nil {
		return nil
	}
	// Flush the per-event log first.
	var firstErr error
	if err := r.out.w.Flush(); err != nil {
		firstErr = err
	}
	if err := r.out.f.Close(); err != nil && firstErr == nil {
		firstErr = err
	}
	if firstErr != nil {
		return firstErr
	}

	// Per-individual summary.
	idiscPath := r.prefix + ".diff.indv.switch"
	f, err := iohelper.OpenWriter(idiscPath)
	if err != nil {
		return fmt.Errorf("opening %s: %w", idiscPath, err)
	}
	defer f.Close()
	w := bufio.NewWriter(f)
	if _, err := w.WriteString("INDV\tN_COMMON_PHASED_HET\tN_SWITCH\tSWITCH\n"); err != nil {
		return fmt.Errorf("writing header to %s: %w", idiscPath, err)
	}
	for i, pair := range r.pairs {
		var rate float64
		if r.nCommonPhasedHet[i] > 0 {
			rate = float64(r.nSwitch[i]) / float64(r.nCommonPhasedHet[i])
		}
		if _, err := fmt.Fprintf(w, "%s\t%d\t%d\t%s\n",
			pair.f1Name, r.nCommonPhasedHet[i], r.nSwitch[i],
			formatSwitchRate(rate)); err != nil {
			return fmt.Errorf("writing %s: %w", idiscPath, err)
		}
	}
	if err := w.Flush(); err != nil {
		return fmt.Errorf("flushing %s: %w", idiscPath, err)
	}
	return nil
}

// formatSwitchRate formats a switch-error rate the way upstream's
// `ostream << double` does — C++ iostream's default `operator<<(double)`
// uses precision 6 with %g-style rounding and trims trailing zeros from
// the significand. Go's `strconv.FormatFloat(v, 'g', 6, 64)` matches that
// exactly for non-special values; the rate is in [0, 1] here so we never
// hit the exponent-formatting threshold.
//
// Special case: exact zero is printed as "0" (not "0.000000") to match
// upstream where N_phased_het_sites == 0 sets switch_error = 0 and the
// ostream prints it as "0".
func formatSwitchRate(v float64) string {
	if v == 0 {
		return "0"
	}
	return strconv.FormatFloat(v, 'g', 6, 64)
}

// splitPhasedAlleles parses a diploid GT string into its two allele strings
// and a phased flag. We keep the alleles as strings so we can compare them by
// value without committing to integer parsing (allele indices > 9 work too).
//
// Returns ("", "", false) for empty / haploid / malformed GT. A haploid call
// (no '/' or '|' separator) is not considered a phased het site by upstream
// here — the inner branch requires `phase1 == '|' && phase2 == '|'` *and*
// the genotype-pair set-match plus diploid heterozygosity, which a haploid
// site can't satisfy.
//
// "." or empty in either slot is returned as that literal string so the
// caller can drop the site via the missing-allele guard.
func splitPhasedAlleles(gt string) (a1, a2 string, phased bool) {
	if gt == "" {
		return "", "", false
	}
	sep := -1
	for i := 0; i < len(gt); i++ {
		if gt[i] == '/' || gt[i] == '|' {
			sep = i
			break
		}
	}
	if sep < 0 {
		// Haploid: not relevant for switch-error.
		return "", "", false
	}
	phased = gt[sep] == '|'
	return gt[:sep], gt[sep+1:], phased
}
