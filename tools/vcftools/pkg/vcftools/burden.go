package vcftools

import (
	"bufio"
	"fmt"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// indvBurdenRunner accumulates per-individual diploid genotype counts and
// writes `<prefix>.iburden`. Mirrors upstream
// variant_file_output.cpp:378-498 (variant_file::output_indv_burden):
//
//   - Only fully diploid sites contribute (upstream's
//     `if (e->is_diploid() == false) continue;` at line 429-433). A
//     site is "fully diploid" when every kept-individual GT is either
//     missing-haploid or `a/b` / `a|b` with two visible alleles.
//   - With `--derived`, the site's INFO/AA tag picks the ancestral
//     allele index; sites where AA is missing (".", "?", absent) or
//     does not match any REF/ALT allele are skipped (upstream's
//     `continue` branches at lines 440-462). Counts then split into
//     hom-anc / het / hom-der instead of hom-ref / het / hom-alt.
//   - For each kept individual at a kept site, exactly one of
//     hom-ref / het / hom-alt / missing is incremented, matching
//     upstream's chained `if/else if/else` at lines 471-483.
//   - Output: TSV with header
//     `INDV\tN_HOM_REF\tN_HET\tN_HOM_ALT\tN_MISS` (or `N_HOM_ANC` /
//     `N_HOM_DER` under --derived; upstream lines 402-405). One data
//     row per kept individual, in the input VCF sample order.
type indvBurdenRunner struct {
	samples []string
	derived bool
	homRef  []int
	het     []int
	homAlt  []int
	missing []int
}

// newIndvBurdenRunner constructs a runner sized to len(samples). Pass the
// post-filter sample list (the same one fed to the per-variant addVariant
// call). When derived is true the header columns are renamed
// HOM_ANC / HOM_DER and per-site allele counts are interpreted with the
// ancestral allele in the "reference" slot.
func newIndvBurdenRunner(samples []string, derived bool) *indvBurdenRunner {
	n := len(samples)
	return &indvBurdenRunner{
		samples: append([]string(nil), samples...),
		derived: derived,
		homRef:  make([]int, n),
		het:     make([]int, n),
		homAlt:  make([]int, n),
		missing: make([]int, n),
	}
}

// addVariant updates the per-individual counters for v. The variant's
// Samples slice is assumed to be aligned with the runner's samples slice
// (this is the contract our other runners use). Upstream warns once when
// it sees a non-diploid site and skips it; we do the same shape (skip),
// without re-emitting the warning byte-for-byte.
func (r *indvBurdenRunner) addVariant(v *vcf.Variant) {
	if r == nil || len(r.samples) == 0 {
		return
	}
	// Diploid check across the kept individuals — upstream's
	// is_diploid() returns false if any kept individual has a non-diploid
	// GT (haploid or missing-haploid). We mirror the same predicate.
	if !isFullyDiploid(v) {
		return
	}

	// Derive the ancestral-allele index. Upstream's loop at lines 446-456
	// matches AA (upper-cased) against e->get_allele(ui) for each allele
	// — REF and every ALT. For a non-derived run aaIdx stays 0
	// (REF), matching upstream's initialisation at line 410.
	aaIdx := 0
	if r.derived {
		idx, ok := ancestralAlleleIndex(v)
		if !ok {
			return
		}
		aaIdx = idx
	}

	for i := range r.samples {
		if i >= len(v.Samples) {
			r.missing[i]++
			continue
		}
		a1, a2, ok := diploidAlleles(getGT(v, i))
		if !ok {
			// Missing-haploid genotype contributes to the missing
			// counter for that individual (upstream's `else` branch
			// at line 481).
			r.missing[i]++
			continue
		}
		switch {
		case a1 == aaIdx && a2 == aaIdx:
			r.homRef[i]++
		case a1 >= 0 && a2 >= 0 && a1 != a2:
			r.het[i]++
		case a1 >= 0 && a2 >= 0 && a1 == a2:
			r.homAlt[i]++
		default:
			r.missing[i]++
		}
	}
}

// writeOutput flushes the accumulated counts to `<prefix>.iburden`. Safe
// to call on a nil receiver. The format is upstream-byte-identical
// (variant_file_output.cpp:402-497): header row, one row per kept
// individual, fields tab-separated and rows newline-terminated.
func (r *indvBurdenRunner) writeOutput(prefix string) error {
	if r == nil {
		return nil
	}
	f, err := iohelper.OpenWriter(prefix + ".iburden")
	if err != nil {
		return err
	}
	defer f.Close()
	w := bufio.NewWriter(f)
	defer w.Flush()
	if r.derived {
		if _, err := w.WriteString("INDV\tN_HOM_ANC\tN_HET\tN_HOM_DER\tN_MISS\n"); err != nil {
			return err
		}
	} else {
		if _, err := w.WriteString("INDV\tN_HOM_REF\tN_HET\tN_HOM_ALT\tN_MISS\n"); err != nil {
			return err
		}
	}
	for i, s := range r.samples {
		if _, err := fmt.Fprintf(w, "%s\t%d\t%d\t%d\t%d\n", s, r.homRef[i], r.het[i], r.homAlt[i], r.missing[i]); err != nil {
			return err
		}
	}
	return nil
}

// indvFreqBurdenRunner accumulates a per-individual × per-allele-count
// matrix and writes `<prefix>.ifreqburden`. Mirrors upstream
// variant_file_output.cpp:501-627 (variant_file::output_indv_freq_burden).
//
// The matrix dimensions are N_kept × (2*N_kept + 1). For each kept site
// (diploid only, and with a resolvable ancestral allele when --derived
// is supplied), upstream computes the per-allele count vector across the
// kept individuals (upstream's `get_allele_counts`, lines 590-591) and
// for each kept individual's diploid genotype (a1, a2):
//
//	if (a1 != aaIdx && a1 >= 0) burden_matrix[indv][allele_counts[a1]]++
//	if (doubleCountHomAlt == 0 || a1 != a2)
//	  if (a2 != aaIdx && a2 >= 0) burden_matrix[indv][allele_counts[a2]]++
//
// (lines 603-609). `--indv-freq-burden` sets doubleCountHomAlt=0;
// `--indv-freq-burden2` sets it to 1 (so hom-alt contributes 1 instead
// of 2). When --derived is not supplied aaIdx is 0 (the REF slot).
//
// Output header is `INDV\t0\t1\t...\t(2*N)` then one row per kept
// individual.
//
// Upstream BUG NOTE (FIXED in this port; tracked in docs/UPSTREAM_BUGS.md):
// upstream variant_file_output.cpp:621 emits `meta_data.indv[indv_count]`
// for the per-row INDV label, where `indv_count` is the kept-position
// index, NOT the original-index `ui`. That decouples the INDV label
// from the burden values when --remove-indv drops a non-trailing
// sample (e.g. removing S2 from [S1,S2,S3,S4] yields labels
// `[S1,S2,S3]` for the rows that actually contain S1/S3/S4 data).
// The companion `output_indv_burden` (lines 488-497) uses
// `meta_data.indv[ui]` and is unaffected. Per CLAUDE.md ("don't
// replicate upstream bugs"), we emit the CORRECT kept-sample label.
type indvFreqBurdenRunner struct {
	samples        []string // kept (post-filter) sample names, in kept order
	n              int      // == len(samples)
	maxChrCount    int      // 2 * n
	burden         [][]int  // [n][maxChrCount+1]
	doubleCountHom bool
	derived        bool
}

// newIndvFreqBurdenRunner sets up the matrix. `keptSamples` is the
// post-sample-filter list — rows of the matrix, the genotypes we read
// per variant, AND the leading INDV column label (correcting the
// upstream `meta_data.indv[indv_count]` bug described above).
func newIndvFreqBurdenRunner(keptSamples []string, doubleCountHomAlt, derived bool) *indvFreqBurdenRunner {
	n := len(keptSamples)
	maxChr := 2 * n
	burden := make([][]int, n)
	for i := range burden {
		burden[i] = make([]int, maxChr+1)
	}
	return &indvFreqBurdenRunner{
		samples:        append([]string(nil), keptSamples...),
		n:              n,
		maxChrCount:    maxChr,
		burden:         burden,
		doubleCountHom: doubleCountHomAlt,
		derived:        derived,
	}
}

// addVariant updates the per-individual frequency-burden matrix for v.
// The variant's Samples slice is assumed to align with r.samples.
func (r *indvFreqBurdenRunner) addVariant(v *vcf.Variant) {
	if r == nil || r.n == 0 {
		return
	}
	if !isFullyDiploid(v) {
		return
	}
	aaIdx := 0
	if r.derived {
		idx, ok := ancestralAlleleIndex(v)
		if !ok {
			return
		}
		aaIdx = idx
	}
	// Per-allele non-missing chromosome counts across the kept
	// individuals — matches upstream `e->get_allele_counts(allele_counts,
	// N_non_missing_chr)` at line 590 with `include_indv` mask applied.
	nAlleles := 1 + len(v.Alt)
	alleleCounts := make([]int, nAlleles)
	for i := 0; i < r.n; i++ {
		if i >= len(v.Samples) {
			continue
		}
		a1, a2, ok := diploidAlleles(getGT(v, i))
		if !ok {
			continue
		}
		if a1 >= 0 && a1 < nAlleles {
			alleleCounts[a1]++
		}
		if a2 >= 0 && a2 < nAlleles {
			alleleCounts[a2]++
		}
	}

	for i := 0; i < r.n; i++ {
		if i >= len(v.Samples) {
			continue
		}
		a1, a2, ok := diploidAlleles(getGT(v, i))
		if !ok {
			continue
		}
		// First allele: count when it is not the ancestral allele and
		// not missing (upstream line 603).
		if a1 != aaIdx && a1 >= 0 && a1 < nAlleles {
			c := alleleCounts[a1]
			if c >= 0 && c <= r.maxChrCount {
				r.burden[i][c]++
			}
		}
		// Second allele: count when not hom-alt OR doubleCountHomAlt is
		// off (upstream line 605: `(double_count_hom_alt == 0) ||
		// (geno.first != geno.second)`). Note that doubleCountHomAlt=0
		// means: always count the second allele (no skip); the
		// `--indv-freq-burden2` flag sets it to 1, which means: skip
		// the second-allele count when both alleles match. Read
		// upstream carefully — the variable name is inverted from how
		// the if-condition reads.
		if !r.doubleCountHom || a1 != a2 {
			if a2 != aaIdx && a2 >= 0 && a2 < nAlleles {
				c := alleleCounts[a2]
				if c >= 0 && c <= r.maxChrCount {
					r.burden[i][c]++
				}
			}
		}
	}
}

// writeOutput flushes the matrix to `<prefix>.ifreqburden`. Safe to call
// on a nil receiver.
func (r *indvFreqBurdenRunner) writeOutput(prefix string) error {
	if r == nil {
		return nil
	}
	f, err := iohelper.OpenWriter(prefix + ".ifreqburden")
	if err != nil {
		return err
	}
	defer f.Close()
	w := bufio.NewWriter(f)
	defer w.Flush()
	// Header: INDV\t0\t1\t...\t(2*N)
	var sb strings.Builder
	sb.WriteString("INDV")
	for i := 0; i <= r.maxChrCount; i++ {
		fmt.Fprintf(&sb, "\t%d", i)
	}
	sb.WriteByte('\n')
	if _, err := w.WriteString(sb.String()); err != nil {
		return err
	}
	for i := 0; i < r.n; i++ {
		if _, err := w.WriteString(r.samples[i]); err != nil {
			return err
		}
		for j := 0; j <= r.maxChrCount; j++ {
			if _, err := fmt.Fprintf(w, "\t%d", r.burden[i][j]); err != nil {
				return err
			}
		}
		if _, err := w.WriteString("\n"); err != nil {
			return err
		}
	}
	return nil
}

// isFullyDiploid mirrors upstream's entry::is_diploid for the kept
// individuals: returns true when every kept sample has a diploid GT
// (`a/b` or `a|b` with two visible allele slots). Haploid calls (no
// separator) and empty/missing-dot GTs disqualify the site. Upstream
// emits a one-off warning and skips such sites
// (variant_file_output.cpp:429-433).
func isFullyDiploid(v *vcf.Variant) bool {
	for i := range v.Samples {
		gt := getGT(v, i)
		if gt == "" {
			return false
		}
		sep := -1
		for j := 0; j < len(gt); j++ {
			if gt[j] == '/' || gt[j] == '|' {
				sep = j
				break
			}
		}
		if sep < 0 {
			return false
		}
	}
	return true
}

// diploidAlleles parses a diploid GT into (a1, a2, ok). Missing slots
// (`.`) parse to -1. ok=false when the GT is not diploid (haploid or
// junk) so the caller can treat the site as missing. Mirrors
// upstream's parser for the two-allele branch
// (vcf_entry_setters.cpp:67-101) with the convention that "." in either
// slot yields a -1 (handled here via parseLDhatAllele).
func diploidAlleles(gt string) (a1, a2 int, ok bool) {
	if gt == "" {
		return -1, -1, false
	}
	sep := -1
	for i := 0; i < len(gt); i++ {
		if gt[i] == '/' || gt[i] == '|' {
			sep = i
			break
		}
	}
	if sep < 0 {
		return -1, -1, false
	}
	left, lok := parseLDhatAllele(gt[:sep])
	right, rok := parseLDhatAllele(gt[sep+1:])
	if !lok {
		left = -1
	}
	if !rok {
		right = -1
	}
	return left, right, true
}

// ancestralAlleleIndex implements upstream's `aa_idx` resolution from
// variant_file_output.cpp:437-462. The INFO/AA tag is upper-cased then
// compared (case-insensitive) against REF and each ALT. Returns the
// matching allele index (0 for REF) and ok=true on success. Returns
// ok=false when AA is missing, `.`, `?`, empty, or does not match any
// allele.
func ancestralAlleleIndex(v *vcf.Variant) (int, bool) {
	aa, present := v.Info["AA"]
	if !present || aa == "" || aa == "." || aa == "?" {
		return 0, false
	}
	aaUp := strings.ToUpper(aa)
	if strings.ToUpper(v.Ref) == aaUp {
		return 0, true
	}
	for i, alt := range v.Alt {
		if strings.ToUpper(alt) == aaUp {
			return i + 1, true
		}
	}
	return 0, false
}
