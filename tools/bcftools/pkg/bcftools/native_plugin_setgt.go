// Native port of the upstream `setGT` plugin (plugins/setGT.c). It sets
// genotypes matching a target rule (-t) to a new value (-n).
//
// Supported target rules (-t): "." (partially or completely missing),
// "./." (completely missing), "./x" (partially missing), and "a" (all).
// Supported new-gt rules (-n): "." (missing), "0" (reference), "M" (major
// allele), "m" (minor allele), "p" (phase), "u" (unphase + sort), "i" (invert
// phase), "c:GT" (custom genotype), and the "0p"/"Mp" phased combinations.
//
// The filter-expression (-t q with -i/-e), binomial (-t b:...), random
// (-t r:FLOAT), and read-depth ("X" / FMT/AD) modes require the bcftools filter
// engine and are reported as unsupported in batch 1.
package bcftools

import (
	"fmt"
	"io"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() { registerNativePlugin("setGT", func() NativePlugin { return &setGTPlugin{} }) }

// setGT target-genotype bits (subset of setGT.c's GT_* masks).
const (
	sgtMissing = 1 << iota // ./.
	sgtPartial             // ./x
	sgtAll                 // a
)

// setGT new-genotype bits (subset of setGT.c's GT_* masks).
const (
	sngMissing  = 1 << iota // .
	sngRef                  // 0
	sngMajor                // M
	sngMinor                // m
	sngPhased               // p (or trailing p on 0/M/m)
	sngUnphased             // u
	sngInvPhase             // i
	sngCustom               // c:GT
)

// setGTPlugin implements setGT. It is per-record and parallel; the major/minor
// allele is computed from each record's own FORMAT/GT.
type setGTPlugin struct {
	tgtMask int
	newMask int
	custom  customGT
	stderr  io.Writer
}

// customGT holds a parsed `c:GT` new-genotype template.
type customGT struct {
	alleles []int  // per-position allele code (>=0 literal, or one of the c* sentinels)
	phased  []bool // separator phase preceding each position
	ploidy  int
}

// Custom-allele sentinels mirror setGT.c's MINOR/MAJOR/MISSING_ALLELE codes.
const (
	cMinor   = -1
	cMajor   = -2
	cMissing = -4
)

// Name returns the plugin name.
func (p *setGTPlugin) Name() string { return "setGT" }

// About returns the one-line description, matching setGT.c about().
func (p *setGTPlugin) About() string {
	return "Set genotypes: partially missing to missing, missing to ref/major allele, etc."
}

// Parallel reports true.
func (p *setGTPlugin) Parallel() bool { return true }

// SetStderr wires the host stderr for the end-of-run summary.
func (p *setGTPlugin) SetStderr(w io.Writer) { p.stderr = w }

// Init parses -t and -n, rejecting the modes that need the filter engine.
func (p *setGTPlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	var tgt, ngt string
	for i := 0; i < len(args); i++ {
		a := args[i]
		switch a {
		case "-t", "--target-gt":
			if i+1 >= len(args) {
				return nil, fmt.Errorf("setGT: -t requires an argument")
			}
			i++
			tgt = args[i]
		case "-n", "--new-gt":
			if i+1 >= len(args) {
				return nil, fmt.Errorf("setGT: -n requires an argument")
			}
			i++
			ngt = args[i]
		case "-i", "--include", "-e", "--exclude":
			return nil, fmt.Errorf("setGT: -i/-e filter expressions are not supported in the native plugin")
		case "-s", "--seed":
			return nil, fmt.Errorf("setGT: -s (random seed / -t r) is not supported in the native plugin")
		default:
			if strings.HasPrefix(a, "-t") && len(a) > 2 {
				tgt = a[2:]
				continue
			}
			if strings.HasPrefix(a, "-n") && len(a) > 2 {
				ngt = a[2:]
				continue
			}
			return nil, fmt.Errorf("setGT: unsupported option %q", a)
		}
	}
	if ngt == "" {
		return nil, fmt.Errorf("setGT: expected -n option")
	}
	if tgt == "" {
		return nil, fmt.Errorf("setGT: expected -t option")
	}
	if err := p.parseTarget(tgt); err != nil {
		return nil, err
	}
	if err := p.parseNew(ngt); err != nil {
		return nil, err
	}

	// Add FORMAT/GT to the header if missing (matches setGT.c init()).
	out := &vcf.Header{Samples: hdr.Samples}
	out.MetaInfo = append(out.MetaInfo, hdr.MetaInfo...)
	if !hasFormatHeader(out.MetaInfo, "GT") {
		out.MetaInfo = appendInfoHeader(out.MetaInfo,
			`##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">`)
	}
	return out, nil
}

// parseTarget parses the -t value into the target bitmask.
func (p *setGTPlugin) parseTarget(s string) error {
	switch s {
	case ".":
		p.tgtMask = sgtMissing | sgtPartial
	case "./x":
		p.tgtMask = sgtPartial
	case "./.":
		p.tgtMask = sgtMissing
	case "a":
		p.tgtMask = sgtAll
	case "q", "?":
		return fmt.Errorf("setGT: -t q (filter query) is not supported in the native plugin")
	default:
		if strings.HasPrefix(s, "r:") {
			return fmt.Errorf("setGT: -t r:FLOAT (random) is not supported in the native plugin")
		}
		if strings.Contains(s, "b") {
			return fmt.Errorf("setGT: -t b:... (binomial) is not supported in the native plugin")
		}
		return fmt.Errorf("setGT: unknown -t value %q", s)
	}
	return nil
}

// parseNew parses the -n value into the new-gt bitmask and any custom template.
func (p *setGTPlugin) parseNew(s string) error {
	if strings.HasPrefix(s, "c:") {
		p.newMask |= sngCustom
		return p.parseCustom(s[2:])
	}
	for _, c := range s {
		switch c {
		case '.':
			p.newMask |= sngMissing
		case '0':
			p.newMask |= sngRef
		case 'M':
			p.newMask |= sngMajor
		case 'm':
			p.newMask |= sngMinor
		case 'p':
			p.newMask |= sngPhased
		case 'u':
			p.newMask |= sngUnphased
		case 'i':
			p.newMask |= sngInvPhase
		case 'X':
			return fmt.Errorf("setGT: -n X (FMT/AD read depth) is not supported in the native plugin")
		default:
			return fmt.Errorf("setGT: unknown -n value %q", s)
		}
	}
	if p.newMask == 0 {
		return fmt.Errorf("setGT: unknown -n value %q", s)
	}
	return nil
}

// parseCustom parses a `c:GT` template such as "0/0", "m/M", or "./.".
func (p *setGTPlugin) parseCustom(s string) error {
	i := 0
	for i < len(s) {
		c := s[i]
		switch {
		case c == 'm':
			p.custom.alleles = append(p.custom.alleles, cMinor)
			p.custom.phased = append(p.custom.phased, false)
			p.newMask |= sngMinor
			i++
		case c == 'M':
			p.custom.alleles = append(p.custom.alleles, cMajor)
			p.custom.phased = append(p.custom.phased, false)
			p.newMask |= sngMajor
			i++
		case c == '.':
			p.custom.alleles = append(p.custom.alleles, cMissing)
			p.custom.phased = append(p.custom.phased, false)
			i++
		case c == '/' || c == '|':
			if len(p.custom.phased) == 0 {
				return fmt.Errorf("setGT: could not parse custom genotype %q", s)
			}
			p.custom.phased[len(p.custom.phased)-1] = c == '|'
			i++
		case c >= '0' && c <= '9':
			j := i
			for j < len(s) && s[j] >= '0' && s[j] <= '9' {
				j++
			}
			n, err := strconv.Atoi(s[i:j])
			if err != nil {
				return fmt.Errorf("setGT: could not parse custom genotype %q", s)
			}
			p.custom.alleles = append(p.custom.alleles, n)
			p.custom.phased = append(p.custom.phased, false)
			i = j
		default:
			return fmt.Errorf("setGT: could not parse custom genotype %q", s)
		}
	}
	p.custom.ploidy = len(p.custom.alleles)
	// The phasing sign comes before the allele; shift it to the right one,
	// matching setGT.c's post-parse fixup.
	for k := p.custom.ploidy - 1; k > 0; k-- {
		p.custom.phased[k] = p.custom.phased[k-1]
	}
	if p.custom.ploidy > 0 {
		p.custom.phased[0] = false
	}
	return nil
}

// Process applies the configured target/new-gt rule to each sample genotype.
func (p *setGTPlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	if len(v.Samples) == 0 {
		return []*vcf.Variant{v}, nil
	}
	nals := 1 + len(v.Alt)

	// Resolve major/minor allele if needed (per record, from FORMAT/GT).
	majorAllele, minorAllele := 0, 0
	if p.newMask&(sngMajor|sngMinor) != 0 {
		ac := majorAlleleCounts(v, nals)
		majorAllele = argmax(ac)
		minorAllele = secondArgmax(ac)
	}

	for i := range v.Samples {
		gt, ok := sampleGT(v, i)
		if !ok {
			continue
		}
		ploidy := gt.ploidy()
		nmiss := gt.nMissing()

		doSet := false
		switch {
		case p.tgtMask&sgtAll != 0:
			doSet = true
		case p.tgtMask&sgtPartial != 0 && nmiss > 0:
			doSet = true
		case p.tgtMask&sgtMissing != 0 && ploidy == nmiss:
			doSet = true
		}
		if !doSet {
			continue
		}

		newGT, changed := p.applyNewGT(gt, majorAllele, minorAllele)
		if changed {
			v.Samples[i].Data["GT"] = newGT.String()
		}
	}
	return []*vcf.Variant{v}, nil
}

// Destroy releases resources (none held).
func (p *setGTPlugin) Destroy() error { return nil }

// applyNewGT transforms one genotype per the new-gt rule, returning the new
// genotype and whether it changed.
func (p *setGTPlugin) applyNewGT(gt genotype, major, minor int) (genotype, bool) {
	switch {
	case p.newMask&sngUnphased != 0:
		return unphaseGT(gt)
	case p.newMask == sngPhased:
		return phaseGT(gt)
	case p.newMask&sngCustom != 0:
		return p.applyCustom(gt, major, minor)
	case p.newMask&sngInvPhase != 0:
		return invertPhaseGT(gt)
	default:
		// Set every present (non-vector-end) allele to the target allele.
		allele := 0
		isMissing := false
		switch {
		case p.newMask&sngMissing != 0:
			isMissing = true
		case p.newMask&sngRef != 0:
			allele = 0
		case p.newMask&sngMajor != 0:
			allele = major
		case p.newMask&sngMinor != 0:
			allele = minor
		}
		phased := p.newMask&sngPhased != 0
		return setAllAlleles(gt, allele, isMissing, phased)
	}
}

// setAllAlleles sets every allele in gt to allele (or missing), applying the
// phase flag to alleles after the first. It mirrors setGT.c set_gt().
func setAllAlleles(gt genotype, allele int, missing, phased bool) (genotype, bool) {
	changed := false
	out := genotype{
		alleles: append([]int(nil), gt.alleles...),
		phased:  append([]bool(nil), gt.phased...),
	}
	// htslib replaces each allele with new_gt, an encoded value that carries
	// its own phase bit: bcf_gt_unphased(a) unless -n p was given. Setting the
	// allele therefore also overwrites the separator phase (and missing alleles
	// are always unphased).
	newPhase := phased && !missing
	for j := range out.alleles {
		newA := allele
		if missing {
			newA = missingAllele
		}
		if out.alleles[j] != newA || (j > 0 && out.phased[j] != newPhase) {
			changed = true
		}
		out.alleles[j] = newA
		if j > 0 {
			out.phased[j] = newPhase
		}
	}
	return out, changed
}

// phaseGT adds phasing to every allele after the first, mirroring phase_gt().
func phaseGT(gt genotype) (genotype, bool) {
	out := genotype{
		alleles: append([]int(nil), gt.alleles...),
		phased:  append([]bool(nil), gt.phased...),
	}
	changed := false
	for j := 1; j < len(out.phased); j++ {
		if !out.phased[j] {
			out.phased[j] = true
			changed = true
		}
	}
	return out, changed
}

// unphaseGT removes phasing and sorts the alleles ascending, mirroring
// unphase_gt() (which insertion-sorts the allele values after unphasing).
func unphaseGT(gt genotype) (genotype, bool) {
	out := genotype{
		alleles: append([]int(nil), gt.alleles...),
		phased:  append([]bool(nil), gt.phased...),
	}
	changed := false
	for j := 1; j < len(out.phased); j++ {
		if out.phased[j] {
			out.phased[j] = false
			changed = true
		}
	}
	// Sort the allele values ascending; missing (-1) sorts first, which matches
	// htslib's encoding where bcf_gt_missing is the smallest GT value.
	sort.Ints(out.alleles)
	return out, changed
}

// invertPhaseGT inverts the phase of a diploid genotype, mirroring
// invert_phase_gt() (only applied for ploidy 2).
func invertPhaseGT(gt genotype) (genotype, bool) {
	if gt.ploidy() != 2 {
		return gt, false
	}
	// Upstream invert_phase_gt proceeds even when an allele is missing; it only
	// bails on a vector-end, which in the textual model means ploidy != 2. It
	// swaps the two allele values and clears the phase on the first.
	out := genotype{
		alleles: []int{gt.alleles[1], gt.alleles[0]},
		phased:  []bool{false, gt.phased[1]},
	}
	return out, true
}

// applyCustom applies a parsed c:GT template, mirroring set_gt_custom().
func (p *setGTPlugin) applyCustom(gt genotype, major, minor int) (genotype, bool) {
	out := genotype{}
	changed := false
	for i := 0; i < p.custom.ploidy; i++ {
		var a int
		switch p.custom.alleles[i] {
		case cMinor:
			a = minor
		case cMajor:
			a = major
		case cMissing:
			a = missingAllele
		default:
			a = p.custom.alleles[i]
		}
		out.alleles = append(out.alleles, a)
		out.phased = append(out.phased, p.custom.phased[i])
		changed = true
	}
	// Compare with the original to decide if anything actually changed.
	if sameGT(gt, out) {
		return gt, false
	}
	return out, changed
}

// sameGT reports whether two genotypes render identically.
func sameGT(a, b genotype) bool { return a.String() == b.String() }

// argmax returns the index of the maximum value (first on ties).
func argmax(v []int) int {
	imax := 0
	for i := 1; i < len(v); i++ {
		if v[i] > v[imax] {
			imax = i
		}
	}
	return imax
}

// secondArgmax returns the index of the second-largest value, mirroring the
// minor-allele selection in setGT.c (imax2 search).
func secondArgmax(v []int) int {
	if len(v) <= 1 {
		return 0
	}
	imax := argmax(v)
	imax2 := 0
	if imax == 0 {
		if len(v) > 1 {
			imax2 = 1
		}
	}
	for i := 0; i < len(v); i++ {
		if i != imax && v[imax2] < v[i] {
			imax2 = i
		}
	}
	return imax2
}

// hasFormatHeader reports whether a ##FORMAT line for id is present.
func hasFormatHeader(meta []string, id string) bool {
	for _, m := range meta {
		if headerKind(m) == "##FORMAT" && headerID(m) == id {
			return true
		}
	}
	return false
}
