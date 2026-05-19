package vcftools

import (
	"bufio"
	"fmt"
	"os"
	"strconv"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// IMPUTE output writer, ported from upstream vcftools
// src/cpp/variant_file_format_convert.cpp:504-614 (output_as_IMPUTE).
//
// CLI flag is `--IMPUTE` (case-sensitive). Upstream registers it in
// parameters.cpp:255 as:
//
//	output_as_IMPUTE = true; phased_only = true;
//	min_site_call_rate = 1.0; min_alleles = 2; max_alleles = 2;
//	num_outputs++;
//
// Effects:
//   - phased_only ON (filter_sites_by_phase drops any site with an unphased
//     kept-individual GT; haploid GTs are considered phased upstream).
//   - min_alleles == max_alleles == 2 (biallelic only).
//   - min_site_call_rate == 1.0 (no missing genotypes allowed at the site,
//     enforced inside output_as_IMPUTE itself rather than relying on the
//     generic filter; see lines 556-584).
//
// Output is a three-file bundle:
//
//	<prefix>.impute.legend    header "ID pos allele0 allele1\n" + one line
//	                          per site. ID is the VCF ID column if not ".",
//	                          else "CHROM-POS".
//	<prefix>.impute.hap       one space-separated row per kept site: two
//	                          integers per kept sample (allele indices, 0 or
//	                          1). Tokens are space-joined; no leading space.
//	<prefix>.impute.hap.indv  one kept-sample name per line (no header).
//
// Upstream errors out if --IMPUTE is paired with --stdout
// (parameters.cpp:734-735). We don't error here because --stdout in our
// port routes through OpenWriter same as --out; each file is a separate
// path.

// imputeRunner buffers per-site rows for --IMPUTE. We write the three
// files at flush time so the per-site missing/unphased check (which
// upstream applies inline in output_as_IMPUTE rather than via the generic
// filter pipeline) can short-circuit before we emit anything.
type imputeRunner struct {
	prefix  string
	samples []string

	// legendEntries collects one "id pos allele0 allele1" line per kept
	// site. We write them lazily at flush so we don't have to keep three
	// open writers across the variant-stream loop. (We could keep them
	// open; buffering is just as cheap and matches the other runners in
	// this package — see ldhatRunner / ldhelmetRunner.)
	legendEntries []string

	// hapRows collects one row per kept site. Each row is the space-joined
	// list of "<allele1> <allele2>" pairs across kept samples.
	hapRows []string

	// warnedMulti tracks whether we've emitted the upstream one-off
	// warning for non-biallelic sites. Upstream emits it via
	// LOG.one_off_warning so it fires at most once per run.
	warnedMulti bool
}

// newImputeRunner constructs the runner. samples is the *filtered* sample
// list (so --indv / --keep / --remove are already honoured upstream of us).
func newImputeRunner(prefix string, samples []string) *imputeRunner {
	return &imputeRunner{
		prefix:  prefix,
		samples: append([]string(nil), samples...),
	}
}

// addVariant buffers one variant's contribution to the IMPUTE output. The
// upstream per-site checks (biallelic; every kept sample has a phased,
// non-missing diploid GT) are enforced here, mirroring lines 550-584 of
// variant_file_format_convert.cpp.
func (r *imputeRunner) addVariant(v *vcf.Variant) {
	if r == nil {
		return
	}

	// Biallelic-only. Upstream emits a one-off warning via
	// LOG.one_off_warning("\tIMPUTE: Only outputting biallelic loci.")
	// We surface the same warning at most once per run.
	if len(v.Alt) != 1 {
		if !r.warnedMulti {
			fmt.Fprintln(os.Stderr, "\tIMPUTE: Only outputting biallelic loci.")
			r.warnedMulti = true
		}
		return
	}

	// Per-site missingness / phase check: any kept sample with a missing
	// or unphased GT disqualifies the site (upstream lines 556-584).
	// Note: this is upstream's filter_sites_by_phase + the inline
	// "missing == true" branch combined. We already passed
	// filter_sites_by_phase via params.Phased (set by parameters.cpp:255
	// via phased_only); this loop covers the inline allele.first < 0 /
	// alleles.second < 0 / PHASE != '|' guard.
	alleles := make([]int, 0, 2*len(r.samples))
	for i := range r.samples {
		a1, a2, phased := parseGTForLDhat(getGT(v, i))
		if !phased {
			return
		}
		if a1 < 0 || a2 < 0 {
			return
		}
		// Upstream only writes alleles 0 or 1 (biallelic enforced
		// above), but defensively guard so we never emit a garbage
		// index. Multi-allele indices on a biallelic site indicate a
		// malformed VCF; we treat that as missing-data and drop the
		// site, matching upstream's get_indv_GENOTYPE_ids returning a
		// negative for out-of-range.
		if a1 > 1 || a2 > 1 {
			return
		}
		alleles = append(alleles, a1, a2)
	}

	// Build the legend line. ID defaults to "CHROM-POS" if the VCF ID
	// column is "." (upstream lines 586-591).
	id := v.ID
	if id == "" || id == "." {
		id = fmt.Sprintf("%s-%d", v.Chrom, v.Pos)
	}
	r.legendEntries = append(r.legendEntries,
		fmt.Sprintf("%s %d %s %s", id, v.Pos, v.Ref, v.Alt[0]))

	// Build the hap row: "a1 a2 a1 a2 ..." with single spaces.
	var row []byte
	for j, a := range alleles {
		if j > 0 {
			row = append(row, ' ')
		}
		row = strconv.AppendInt(row, int64(a), 10)
	}
	r.hapRows = append(r.hapRows, string(row))
}

// close writes the three IMPUTE output files. Safe to call on a nil
// receiver.
func (r *imputeRunner) close() error {
	if r == nil {
		return nil
	}

	legendPath := r.prefix + ".impute.legend"
	hapPath := r.prefix + ".impute.hap"
	indvPath := r.prefix + ".impute.hap.indv"

	// Legend file: header + one line per kept site.
	legendW, err := iohelper.OpenWriter(legendPath)
	if err != nil {
		return fmt.Errorf("opening %s: %w", legendPath, err)
	}
	defer legendW.Close()
	legendBuf := bufio.NewWriter(legendW)
	legendBuf.WriteString("ID pos allele0 allele1\n")
	for _, entry := range r.legendEntries {
		legendBuf.WriteString(entry)
		legendBuf.WriteByte('\n')
	}
	if err := legendBuf.Flush(); err != nil {
		return fmt.Errorf("writing %s: %w", legendPath, err)
	}

	// Hap file: one row per kept site (no header).
	hapW, err := iohelper.OpenWriter(hapPath)
	if err != nil {
		return fmt.Errorf("opening %s: %w", hapPath, err)
	}
	defer hapW.Close()
	hapBuf := bufio.NewWriter(hapW)
	for _, row := range r.hapRows {
		hapBuf.WriteString(row)
		hapBuf.WriteByte('\n')
	}
	if err := hapBuf.Flush(); err != nil {
		return fmt.Errorf("writing %s: %w", hapPath, err)
	}

	// Indv file: one kept sample name per line.
	indvW, err := iohelper.OpenWriter(indvPath)
	if err != nil {
		return fmt.Errorf("opening %s: %w", indvPath, err)
	}
	defer indvW.Close()
	indvBuf := bufio.NewWriter(indvW)
	for _, name := range r.samples {
		indvBuf.WriteString(name)
		indvBuf.WriteByte('\n')
	}
	if err := indvBuf.Flush(); err != nil {
		return fmt.Errorf("writing %s: %w", indvPath, err)
	}

	return nil
}
