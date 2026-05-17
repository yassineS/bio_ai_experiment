package vcftools

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/vcf"
)

// LDhat / LDhat-geno output writers, ported from upstream vcftools
// src/cpp/variant_file_format_convert.cpp:
//   - output_as_LDhat_phased   (--ldhat, lines 616-797)
//   - output_as_LDhat_unphased (--ldhat-geno, lines 799-980)
//
// Both flags write a paired (<prefix>.ldhat.sites, <prefix>.ldhat.locs) bundle.
// Upstream requires exactly one chromosome (--chr) to be selected; this is
// enforced in src/cpp/parameters.cpp:717. `--ldhat` additionally implies
// --phased (phased_only=true) which we model in passFilters.
//
// Phased layout (sites file):
//
//	<2*N_indv>\t<n_sites>\t1
//	>S1-0
//	<n_sites alleles>
//	>S1-1
//	<n_sites alleles>
//	...
//
// Each haplotype emits one character per kept site: '0', '1', or '?' for
// missing or out-of-range alleles (and for unphased GTs that survive the
// site-level phase filter via --ldhat's phased_only).
//
// Unphased layout (sites file):
//
//	<N_indv>\t<n_sites>\t2
//	>S1
//	<n_sites alleles>
//	...
//
// Per-sample mapping of (a1,a2) -> char follows upstream's switch in
// variant_file_format_convert.cpp:903-936:
//
//	(0,0) -> '0'              homozygous ref
//	(0,1) or (1,0) -> '2'     heterozygote
//	(1,1) -> '1'              homozygous alt
//	(0,.) phased haploid -> '0'
//	(1,.) phased haploid -> '1'
//	anything else -> '?'      missing / out-of-range
//
// The locs file emits in both cases:
//
//	<n_sites>\t<max_pos/1000.0>\tL\n
//	<pos_1/1000.0>
//	<pos_2/1000.0>
//	...
//
// All decimal numbers use 4 fractional digits ("%.4f" via std::fixed +
// precision(4) in upstream).

// ldhatMode selects between the phased and unphased layouts.
type ldhatMode int

const (
	ldhatPhased ldhatMode = iota
	ldhatUnphased
)

// ldhatSite is a per-site row buffered while streaming. Each entry stores
// one or two characters per kept sample (depending on mode).
type ldhatSite struct {
	pos int
	// For ldhatPhased: 2 chars per sample (samples flattened).
	// For ldhatUnphased: 1 char per sample.
	chars []byte
}

// ldhatRunner buffers per-site genotype rows for --ldhat / --ldhat-geno and
// emits the paired (.ldhat.sites, .ldhat.locs) files at flush time.
type ldhatRunner struct {
	mode    ldhatMode
	prefix  string
	samples []string
	sites   []ldhatSite
	// warned tracks whether we've emitted the upstream one-off warning for
	// non-biallelic sites. Upstream emits it via LOG.one_off_warning so it
	// fires at most once per run.
	warned bool
}

// newLDhatRunner constructs the runner. samples is the *filtered* sample
// list (so --indv / --keep / --remove are already honoured upstream of us).
func newLDhatRunner(mode ldhatMode, prefix string, samples []string) *ldhatRunner {
	return &ldhatRunner{
		mode:    mode,
		prefix:  prefix,
		samples: append([]string(nil), samples...),
	}
}

// addVariant buffers one variant's contribution to the LDhat output. Only
// biallelic loci contribute; multi-allelic sites are silently dropped with
// a single warning to stderr, matching upstream.
func (r *ldhatRunner) addVariant(v *vcf.Variant) {
	if r == nil {
		return
	}
	if len(v.Alt) != 1 {
		if !r.warned {
			fmt.Fprintln(os.Stderr, "\tLDhat: Only outputting biallelic loci.")
			r.warned = true
		}
		return
	}
	switch r.mode {
	case ldhatPhased:
		row := make([]byte, 0, 2*len(r.samples))
		for i := range r.samples {
			a1, a2, phased := parseGTForLDhat(getGT(v, i))
			row = append(row, ldhatPhasedChar(a1, phased))
			row = append(row, ldhatPhasedChar(a2, phased))
		}
		r.sites = append(r.sites, ldhatSite{pos: v.Pos, chars: row})
	case ldhatUnphased:
		row := make([]byte, 0, len(r.samples))
		for i := range r.samples {
			a1, a2, phased := parseGTForLDhat(getGT(v, i))
			row = append(row, ldhatUnphasedChar(a1, a2, phased))
		}
		r.sites = append(r.sites, ldhatSite{pos: v.Pos, chars: row})
	}
}

// close writes the .ldhat.sites and .ldhat.locs files. Safe to call on a
// nil receiver.
func (r *ldhatRunner) close() error {
	if r == nil {
		return nil
	}

	nSites := len(r.sites)
	maxPos := -1
	for _, s := range r.sites {
		if s.pos > maxPos {
			maxPos = s.pos
		}
	}

	locsPath := r.prefix + ".ldhat.locs"
	sitesPath := r.prefix + ".ldhat.sites"

	// Match upstream's empty-output behaviour: max_pos starts at -1, so the
	// "max/1000.0" prefix on the empty case becomes -0.0010 (see
	// variant_file_format_convert.cpp:758 and :944 for the precision-4 format).
	maxPosKb := float64(maxPos) / 1000.0

	locsW, err := iohelper.OpenWriter(locsPath)
	if err != nil {
		return fmt.Errorf("opening %s: %w", locsPath, err)
	}
	defer locsW.Close()
	locsBuf := bufio.NewWriter(locsW)
	// Header: "<n_sites>\t<max/1000>\tL\n"
	fmt.Fprintf(locsBuf, "%d\t%.4f\tL\n", nSites, maxPosKb)
	for _, s := range r.sites {
		fmt.Fprintf(locsBuf, "%.4f\n", float64(s.pos)/1000.0)
	}
	if err := locsBuf.Flush(); err != nil {
		return fmt.Errorf("writing %s: %w", locsPath, err)
	}

	sitesW, err := iohelper.OpenWriter(sitesPath)
	if err != nil {
		return fmt.Errorf("opening %s: %w", sitesPath, err)
	}
	defer sitesW.Close()
	sitesBuf := bufio.NewWriter(sitesW)

	// Sites header: phased uses "2*N\tn_sites\t1" (one row per haplotype),
	// unphased uses "N\tn_sites\t2" (one row per individual, three states).
	var nHapsOrIndv int
	var stateCount int
	if r.mode == ldhatPhased {
		nHapsOrIndv = 2 * len(r.samples)
		stateCount = 1
	} else {
		nHapsOrIndv = len(r.samples)
		stateCount = 2
	}
	fmt.Fprintf(sitesBuf, "%d\t%d\t%d\n", nHapsOrIndv, nSites, stateCount)

	switch r.mode {
	case ldhatPhased:
		// Two haplotype rows per sample. For each sample we walk the buffered
		// sites and pluck out the matching position in `chars`.
		for i, name := range r.samples {
			for k := 0; k < 2; k++ {
				fmt.Fprintf(sitesBuf, ">%s-%d\n", name, k)
				for _, s := range r.sites {
					sitesBuf.WriteByte(s.chars[2*i+k])
				}
				sitesBuf.WriteByte('\n')
			}
		}
	case ldhatUnphased:
		for i, name := range r.samples {
			fmt.Fprintf(sitesBuf, ">%s\n", name)
			for _, s := range r.sites {
				sitesBuf.WriteByte(s.chars[i])
			}
			sitesBuf.WriteByte('\n')
		}
	}

	if err := sitesBuf.Flush(); err != nil {
		return fmt.Errorf("writing %s: %w", sitesPath, err)
	}
	return nil
}

// ldhatPhasedChar mirrors the upstream phased branch
// (variant_file_format_convert.cpp:736-750): emit the digit when the allele
// is callable, '?' otherwise. We also enforce upstream's phased-only invariant
// here: an unphased genotype that somehow survives passFilters (e.g. the
// `--ldhat` site-level filter found exactly zero unphased samples but a
// particular sample is a singleton mismatch — not possible in upstream's
// model but defensive here) yields '?'.
func ldhatPhasedChar(allele int, phased bool) byte {
	if !phased {
		return '?'
	}
	switch allele {
	case 0:
		return '0'
	case 1:
		return '1'
	default:
		return '?'
	}
}

// ldhatUnphasedChar mirrors the upstream switch on alleles.first /
// alleles.second in variant_file_format_convert.cpp:903-936.
//
//	(0,0) -> '0'              homozygous ref
//	(0,1) -> '2', (1,0) -> '2' heterozygote
//	(1,1) -> '1'              homozygous alt
//	(0,-1) phased -> '0'      haploid ref
//	(1,-1) phased -> '1'      haploid alt
//	(0,-2) -> '0'             haploid ref (truly absent second allele)
//	(1,-2) -> '1'             haploid alt (truly absent second allele)
//	anything else -> '?'
//
// We encode "second allele absent" with secondAbsent: in upstream that's
// alleles.second == -2 (haploid call: GT has no separator).
func ldhatUnphasedChar(a1, a2 int, phased bool) byte {
	// First allele branch must be valid (0 or 1) for anything other than '?'.
	switch a1 {
	case 0:
		switch a2 {
		case 0:
			return '0'
		case 1:
			return '2'
		case -1:
			// alleles.second == -1 is "second allele present but unknown"
			// (e.g. trailing dot in "0|."). Upstream maps this to ref iff
			// phased.
			if phased {
				return '0'
			}
			return '?'
		case -2:
			// Haploid: no second allele in the GT string at all.
			return '0'
		default:
			return '?'
		}
	case 1:
		switch a2 {
		case 0:
			return '2'
		case 1:
			return '1'
		case -1:
			if phased {
				return '1'
			}
			return '?'
		case -2:
			return '1'
		default:
			return '?'
		}
	default:
		return '?'
	}
}

// parseGTForLDhat parses a VCF GT string into (allele1, allele2, phased).
// The semantics follow upstream vcftools' parser
// (vcf_entry_setters.cpp:67-101):
//
//   - "a/b" or "a|b" with two visible alleles: a1=a, a2=b, phased iff '|'.
//   - "a" with no separator: haploid, a1=a, a2=-2, phased=true (upstream
//     sets PHASE='|' for the haploid branch at line 90).
//   - A single character "." (missing dot) is handled here as haploid
//     missing: a1=-1, a2=-2, phased=true (upstream's haploid branch with the
//     missing allele).
//   - "." in a slot is interpreted as -1 (missing); a numeric value is
//     parsed as the integer index. Allele indices >1 are kept as their
//     integer values so the upstream switch can map them to '?'.
//
// On parse error (junk in an allele slot) returns (-1,-1,false). Empty
// GT string is the same as missing.
func parseGTForLDhat(gt string) (a1, a2 int, phased bool) {
	if gt == "" || gt == "." {
		// Upstream: pos==string::npos branch -> ploidy=1, PHASE='|',
		// allele1=int(".")=-1, allele2=-2.
		return -1, -2, true
	}
	sep := -1
	for i := 0; i < len(gt); i++ {
		if gt[i] == '/' || gt[i] == '|' {
			sep = i
			break
		}
	}
	if sep < 0 {
		// Haploid call: "0", "1", "2", or just digits. Upstream sets
		// PHASE='|' here unconditionally.
		v, ok := parseLDhatAllele(gt)
		if !ok {
			return -1, -2, true
		}
		return v, -2, true
	}
	phased = gt[sep] == '|'
	left, lOK := parseLDhatAllele(gt[:sep])
	right, rOK := parseLDhatAllele(gt[sep+1:])
	if !lOK {
		left = -1
	}
	if !rOK {
		right = -1
	}
	return left, right, phased
}

// parseLDhatAllele parses a single allele slot from a GT field. Returns
// (-1, true) for "." and (n, true) for numeric strings. Trailing whitespace
// or any leftover after the digit yields (-1, false).
func parseLDhatAllele(s string) (int, bool) {
	if s == "." || s == "" {
		return -1, true
	}
	n, err := strconv.Atoi(strings.TrimSpace(s))
	if err != nil {
		return -1, false
	}
	return n, true
}

// getGT pulls the GT value for the i'th sample, returning "" if the sample
// or its GT entry is missing.
func getGT(v *vcf.Variant, i int) string {
	if i < 0 || i >= len(v.Samples) {
		return ""
	}
	gt, ok := v.Samples[i].Data["GT"]
	if !ok {
		return ""
	}
	return gt
}

// isPhasedSite mirrors upstream's filter_sites_by_phase
// (entry_filters.cpp:989-1010): a site is "phased" iff every kept-individual
// genotype's PHASE byte is '|' (i.e. the separator is '|', or the genotype
// is haploid which upstream codes as PHASE='|' too).
//
// Returns true to keep the site, false to drop it.
func isPhasedSite(v *vcf.Variant) bool {
	for i := range v.Samples {
		gt := getGT(v, i)
		// Use the same parser as parseGTForLDhat so we are consistent.
		// parseGTForLDhat sets phased=true for haploid calls (incl. missing
		// "."), matching upstream's PHASE='|' assignment in
		// vcf_entry_setters.cpp:90. Diploid "a/b" yields phased=false.
		_, _, phased := parseGTForLDhat(gt)
		if !phased {
			return false
		}
	}
	return true
}
