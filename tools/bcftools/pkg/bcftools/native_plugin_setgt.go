// Native port of the upstream `setGT` plugin (plugins/setGT.c). It sets
// genotypes matching a target rule (-t) to a new value (-n).
//
// Supported target rules (-t): "." (partially or completely missing),
// "./." (completely missing), "./x" (partially missing), and "a" (all).
// Supported new-gt rules (-n): "." (missing), "0" (reference), "M" (major
// allele), "m" (minor allele), "X" (allele with the largest FORMAT/AD), "p"
// (phase), "u" (unphase + sort), "i" (invert phase), "c:GT" (custom genotype),
// and the "0p"/"Mp" phased combinations.
//
// The filter-expression mode (-t q with -i/-e) is supported via the native
// filter engine: -i sets matching genotypes, -e sets the complement (matching
// setGT.c's FLT_INCLUDE / FLT_EXCLUDE per-sample handling).
//
// The binomial (-t b:TAG CMP VAL), random (-t r:FLOAT with -s SEED) and
// read-depth (-n X / FMT/AD) modes are all supported natively:
//
//   - binomial: a two-tailed binomial test over a FORMAT integer tag (usually
//     AD) for each diploid heterozygous genotype, computed via the same
//     regularized incomplete-beta function (kfBetai) htslib uses, so it is
//     bit-exact.
//   - random: htslib's deterministic drand48 PRNG (see native_drand48.go),
//     seeded from -s (default 0), gives byte-for-byte parity with upstream.
//   - read-depth: -n X picks, per sample, the allele with the largest
//     FORMAT/AD value and sets every allele of the genotype to it.
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
	sgtQuery               // q: select genotypes via -i/-e filter expression
	sgtBinom               // b: heterozygous GTs failing a two-tailed binomial test
	sgtRand                // r: a random proportion of genotypes
)

// setGT filter-logic values, mirroring setGT.c's FLT_INCLUDE / FLT_EXCLUDE.
const (
	sgtFilterInclude = 1 // -i
	sgtFilterExclude = 2 // -e
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
	sngXVAF                 // X: allele with the largest FORMAT/AD
)

// setGTPlugin implements setGT. It is per-record and parallel; the major/minor
// allele is computed from each record's own FORMAT/GT.
type setGTPlugin struct {
	tgtMask     int
	newMask     int
	custom      customGT
	stderr      io.Writer
	filter      *Filter // compiled -i/-e expression (nil unless -t q)
	filterLogic int     // sgtFilterInclude / sgtFilterExclude
	filterExpr  string  // raw expression text (for error messages)

	// Binomial mode (-t b:TAG CMP VAL).
	binomTag string  // FORMAT integer tag (e.g. "AD")
	binomVal float64 // threshold VAL
	binomCmp func(a, b float64) bool

	// Random mode (-t r:FLOAT with -s SEED).
	randFrac float64 // proportion of genotypes to act on (0<frac<1)
	randSeed int64   // -s value (default 0)
	rng      *drand48
}

// customGT holds a parsed `c:GT` new-genotype template.
type customGT struct {
	alleles []int  // per-position allele code (>=0 literal, or one of the c* sentinels)
	phased  []bool // separator phase preceding each position
	ploidy  int
}

// Custom-allele sentinels mirror setGT.c's MINOR/MAJOR/X_VAF/MISSING_ALLELE codes.
const (
	cMinor   = -1
	cMajor   = -2
	cXVAF    = -3
	cMissing = -4
)

// Name returns the plugin name.
func (p *setGTPlugin) Name() string { return "setGT" }

// About returns the one-line description, matching setGT.c about().
func (p *setGTPlugin) About() string {
	return "Set genotypes: partially missing to missing, missing to ref/major allele, etc."
}

// Parallel reports whether the plugin is safe to run through the per-record
// worker pool. The random mode (-t r) consumes a single shared drand48 stream
// that must advance in strict input order for byte-parity with upstream, so it
// forces serial execution; every other mode is per-record and concurrency-safe.
func (p *setGTPlugin) Parallel() bool { return p.tgtMask&sgtRand == 0 }

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
			if i+1 >= len(args) {
				return nil, fmt.Errorf("setGT: %s requires an argument", a)
			}
			if p.filterExpr != "" {
				return nil, fmt.Errorf("setGT: only one -i or -e expression can be given, and they cannot be combined")
			}
			i++
			p.filterExpr = args[i]
			if a == "-e" || a == "--exclude" {
				p.filterLogic = sgtFilterExclude
			} else {
				p.filterLogic = sgtFilterInclude
			}
		case "-s", "--seed":
			if i+1 >= len(args) {
				return nil, fmt.Errorf("setGT: -s requires an argument")
			}
			i++
			n, err := strconv.ParseInt(args[i], 10, 64)
			if err != nil {
				return nil, fmt.Errorf("setGT: could not parse: -s %s", args[i])
			}
			p.randSeed = n
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

	// Random mode (-t r): if used on its own it implicitly targets all
	// genotypes, and the deterministic drand48 PRNG is seeded once from -s
	// (default 0), exactly as setGT.c init() (lines 269-273).
	if p.tgtMask&sgtRand != 0 {
		if p.tgtMask == sgtRand {
			p.tgtMask |= sgtAll
		}
		p.rng = newDrand48(p.randSeed)
	}

	// The binomial FORMAT tag must be declared in the header (setGT.c:177).
	if p.tgtMask&sgtBinom != 0 && !hasFormatHeader(hdr.MetaInfo, p.binomTag) {
		return nil, fmt.Errorf("setGT: the FORMAT tag \"%s\" is not present in the VCF", p.binomTag)
	}
	// -n X requires FORMAT/AD (setGT.c:322-323).
	if p.newMask&sngXVAF != 0 && !hasFormatHeader(hdr.MetaInfo, "AD") {
		return nil, fmt.Errorf("setGT: the FORMAT/AD annotation does exist, cannot run with --new-gt %s", ngt)
	}

	// -t q must pair with exactly one -i/-e and vice versa (setGT.c:314-316).
	if p.filterExpr != "" && p.tgtMask&sgtQuery == 0 {
		return nil, fmt.Errorf("setGT: expected -t q with -i/-e")
	}
	if p.filterExpr == "" && p.tgtMask&sgtQuery != 0 {
		return nil, fmt.Errorf("setGT: expected -i/-e with -t q")
	}
	if p.filterExpr != "" {
		f, err := CompileFilterWithHeader(p.filterExpr, hdr)
		if err != nil {
			return nil, fmt.Errorf("setGT: %w", err)
		}
		p.filter = f
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

// parseTarget parses the -t value into the target bitmask. It mirrors
// setGT.c init()'s case 't', which OR's together every pattern it recognises
// in the string (so the random "r:FLOAT" selector can be combined with another
// mode, and any 'b' triggers the binomial parser) and only errors when nothing
// matched.
func (p *setGTPlugin) parseTarget(s string) error {
	switch s {
	case ".":
		p.tgtMask |= sgtMissing | sgtPartial
	case "./x":
		p.tgtMask |= sgtPartial
	case "./.":
		p.tgtMask |= sgtMissing
	case "a":
		p.tgtMask |= sgtAll
	case "q", "?":
		p.tgtMask |= sgtQuery
	}
	if strings.HasPrefix(s, "r:") {
		f, err := strconv.ParseFloat(s[2:], 64)
		if err != nil {
			return fmt.Errorf("setGT: could not parse: -t %s", s)
		}
		if f <= 0 || f >= 1 {
			return fmt.Errorf("setGT: expected value between 0 and 1 with -t")
		}
		p.randFrac = f
		p.tgtMask |= sgtRand
	}
	if i := strings.IndexByte(s, 'b'); i >= 0 {
		if err := p.parseBinomExpr(s[i:]); err != nil {
			return err
		}
	}
	if p.tgtMask == 0 {
		return fmt.Errorf("setGT: unknown -t value %q", s)
	}
	return nil
}

// parseBinomExpr parses a "b:TAG CMP VAL" binomial selector, mirroring
// setGT.c's parse_binom_expr(). str begins at the 'b'. TAG is a FORMAT tag,
// CMP one of <, <=, >, >=, ==, =, and VAL a floating-point threshold.
func (p *setGTPlugin) parseBinomExpr(str string) error {
	errExpr := func() error {
		return fmt.Errorf("setGT: error parsing the expression: %s\n"+
			"expected TAG CMP VAL, where TAG is a FORMAT tag, "+
			"CMP one of <, <=, >, >=, and VAL a value", str)
	}
	if len(str) < 2 || str[1] != ':' {
		return errExpr()
	}
	rest := str[2:]
	// Skip leading whitespace before the tag.
	i := 0
	for i < len(rest) && isSpaceByte(rest[i]) {
		i++
	}
	beg := i
	for i < len(rest) {
		c := rest[i]
		if isSpaceByte(c) || c == '<' || c == '=' || c == '>' {
			break
		}
		i++
	}
	if i >= len(rest) || i == beg {
		return errExpr()
	}
	p.binomTag = rest[beg:i]
	// Skip whitespace between tag and operator.
	for i < len(rest) && isSpaceByte(rest[i]) {
		i++
	}
	if i >= len(rest) {
		return errExpr()
	}
	switch {
	case strings.HasPrefix(rest[i:], "<="):
		p.binomCmp = cmpLE
		i += 2
	case strings.HasPrefix(rest[i:], ">="):
		p.binomCmp = cmpGE
		i += 2
	case strings.HasPrefix(rest[i:], "=="):
		p.binomCmp = cmpEQ
		i += 2
	case rest[i] == '<':
		p.binomCmp = cmpLT
		i++
	case rest[i] == '>':
		p.binomCmp = cmpGT
		i++
	case rest[i] == '=':
		p.binomCmp = cmpEQ
		i++
	default:
		return errExpr()
	}
	// Skip whitespace before the value.
	for i < len(rest) && isSpaceByte(rest[i]) {
		i++
	}
	if i >= len(rest) {
		return errExpr()
	}
	valStr, consumed := scanFloatPrefix(rest[i:])
	if consumed == 0 {
		return errExpr()
	}
	v, err := strconv.ParseFloat(valStr, 64)
	if err != nil {
		return errExpr()
	}
	// Anything other than trailing whitespace after the value is an error,
	// matching strtod's end-pointer check in parse_binom_expr.
	tail := rest[i+consumed:]
	for j := 0; j < len(tail); j++ {
		if !isSpaceByte(tail[j]) {
			return errExpr()
		}
	}
	p.binomVal = v
	p.tgtMask |= sgtBinom
	return nil
}

// scanFloatPrefix returns the leading run of bytes that strtod would consume as
// a floating-point literal (sign, digits, decimal point, exponent), and how
// many bytes it spans.
func scanFloatPrefix(s string) (string, int) {
	i := 0
	if i < len(s) && (s[i] == '+' || s[i] == '-') {
		i++
	}
	digits := 0
	for i < len(s) && s[i] >= '0' && s[i] <= '9' {
		i++
		digits++
	}
	if i < len(s) && s[i] == '.' {
		i++
		for i < len(s) && s[i] >= '0' && s[i] <= '9' {
			i++
			digits++
		}
	}
	if digits == 0 {
		return "", 0
	}
	if i < len(s) && (s[i] == 'e' || s[i] == 'E') {
		j := i + 1
		if j < len(s) && (s[j] == '+' || s[j] == '-') {
			j++
		}
		expDigits := 0
		for j < len(s) && s[j] >= '0' && s[j] <= '9' {
			j++
			expDigits++
		}
		if expDigits > 0 {
			i = j
		}
	}
	return s[:i], i
}

// Binomial comparison operators mirror setGT.c's cmp_* functions.
func cmpEQ(a, b float64) bool { return a == b }
func cmpLE(a, b float64) bool { return a <= b }
func cmpGE(a, b float64) bool { return a >= b }
func cmpLT(a, b float64) bool { return a < b }
func cmpGT(a, b float64) bool { return a > b }

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
			p.newMask |= sngXVAF
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
		case c == 'X':
			p.custom.alleles = append(p.custom.alleles, cXVAF)
			p.custom.phased = append(p.custom.phased, false)
			p.newMask |= sngXVAF
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

	// Resolve the -n X allele per sample (the allele with the largest FORMAT/AD).
	// If AD is absent or malformed for the record, setGT.c returns the record
	// unchanged (process(): "else return rec").
	var xAllele []int // per-sample: >=0 allele index, or missingAllele
	if p.newMask&sngXVAF != 0 {
		xv, ok := p.resolveXVAF(v, nals)
		if !ok {
			return []*vcf.Variant{v}, nil
		}
		xAllele = xv
	}

	// Binomial mode (-t b): only diploid heterozygous genotypes are considered.
	if p.tgtMask&sgtBinom != 0 {
		return p.processBinom(v, nals, majorAllele, minorAllele, xAllele)
	}

	// Filter-query mode (-t q): the samples to set are chosen by the per-sample
	// mask of the -i/-e expression rather than by missingness.
	smplPass, ok := p.querySampleMask(v)
	if p.tgtMask&sgtQuery != 0 && !ok {
		// The whole site was filtered out: emit it unchanged (setGT.c returns
		// the record without modifying any genotype).
		return []*vcf.Variant{v}, nil
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
		case p.tgtMask&sgtQuery != 0:
			doSet = smplPass == nil || smplPass[i]
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
		// The random draw is consumed only for samples already selected, and
		// exactly once, mirroring setGT.c's ordering for byte-reproducibility.
		if p.tgtMask&sgtRand != 0 && p.randomDraw() {
			continue
		}

		newGT, changed := p.applyNewGTX(gt, nals, majorAllele, minorAllele, xAllele, i)
		if changed {
			v.Samples[i].Data["GT"] = newGT.String()
		}
	}
	return []*vcf.Variant{v}, nil
}

// processBinom implements the -t b binomial branch of setGT.c process()
// (lines 530-568): for each diploid heterozygous genotype it computes a
// two-tailed binomial test over the FORMAT binom tag (typically AD) indexed by
// the two GT alleles, and sets the genotype when the comparison passes.
func (p *setGTPlugin) processBinom(v *vcf.Variant, nals, major, minor int, xAllele []int) ([]*vcf.Variant, error) {
	smplPass, _ := p.querySampleMask(v)
	for i := range v.Samples {
		// -i/-e filtering applies per sample inside the binomial loop too
		// (setGT.c:535-539). With include logic a non-passing sample is skipped;
		// with exclude logic a passing sample is skipped. querySampleMask already
		// inverts the exclude mask, so here a true mask entry means "set".
		if p.filter != nil && smplPass != nil && !smplPass[i] {
			continue
		}
		gt, ok := sampleGT(v, i)
		if !ok || gt.ploidy() < 2 {
			continue
		}
		a0, a1 := gt.alleles[0], gt.alleles[1]
		if a0 == missingAllele || a1 == missingAllele {
			continue
		}
		if a0 == a1 {
			continue // a hom
		}
		ad, parseOK := sampleADInts(v, i)
		if !parseOK {
			continue
		}
		if a0 >= len(ad) || a1 >= len(ad) {
			return nil, fmt.Errorf("setGT: the sample %s has incorrect number of %s fields at %s:%d",
				v.Samples[i].Name, p.binomTag, v.Chrom, v.Pos)
		}
		prob := calcBinomTwoSided(ad[a0], ad[a1], 0.5)
		if !p.binomCmp(prob, p.binomVal) {
			continue
		}
		if p.tgtMask&sgtRand != 0 && p.randomDraw() {
			continue
		}
		newGT, changed := p.applyNewGTX(gt, nals, major, minor, xAllele, i)
		if changed {
			v.Samples[i].Data["GT"] = newGT.String()
		}
	}
	return []*vcf.Variant{v}, nil
}

// resolveXVAF returns, per sample, the allele index with the largest FORMAT/AD
// value (or missingAllele when every AD entry is missing), mirroring setGT.c's
// GT_X_VAF block. The second return is false when the record lacks a usable
// FORMAT/AD (in which case setGT.c emits the record unchanged).
func (p *setGTPlugin) resolveXVAF(v *vcf.Variant, nals int) ([]int, bool) {
	out := make([]int, len(v.Samples))
	for i := range v.Samples {
		ad, ok := sampleADInts(v, i)
		if !ok || len(ad) != nals {
			// Upstream requires bcf_get_format_int32(...)==n_allele*n_sample for
			// the whole record; any sample short of n_allele AD values aborts the
			// X resolution and leaves the record unchanged.
			return nil, false
		}
		jmax := -1
		for j := 0; j < len(ad); j++ {
			if jmax == -1 || ad[jmax] < ad[j] {
				jmax = j
			}
		}
		out[i] = jmax
	}
	return out, true
}

// randomDraw mirrors setGT.c's random_draw(): it returns true (i.e. "skip this
// genotype") when the next drand48 value exceeds the requested fraction. The
// draw is reversed so a larger fraction keeps more genotypes.
func (p *setGTPlugin) randomDraw() bool {
	return p.rng.float64() > p.randFrac
}

// sampleADInts returns the integer FORMAT/AD values for sample i. The boolean
// is false when AD is absent or contains a missing entry, mirroring setGT.c's
// reliance on bcf_get_format_int32 (a missing value aborts use of the field).
func sampleADInts(v *vcf.Variant, i int) ([]int, bool) {
	ad, ok := v.Samples[i].Data["AD"]
	if !ok {
		return nil, false
	}
	return parseADList(ad)
}

// querySampleMask evaluates the -i/-e filter for the record and returns the
// per-sample mask of genotypes to set, together with whether the site is
// processed at all. It reproduces setGT.c's FLT_INCLUDE / FLT_EXCLUDE handling
// (setGT.c:570-590):
//
//   - With no filter (non-query mode) it returns (nil, true).
//   - With -i: the site is skipped entirely unless it passes; a passing site
//     with a per-sample mask sets only the matching samples (mask returned),
//     and a passing site-level expression with no per-sample component sets
//     every sample (nil mask).
//   - With -e: a site that does NOT pass sets every sample; a site that passes
//     has its per-sample mask inverted (samples that matched are left alone,
//     the rest are set), and if no sample remains after inversion the site is
//     skipped. A passing site-level expression with no per-sample component
//     skips the whole site.
func (p *setGTPlugin) querySampleMask(v *vcf.Variant) (mask []bool, process bool) {
	if p.filter == nil {
		return nil, true
	}
	passSite, smpl := p.filter.EvalSamples(v)
	if p.filterLogic == sgtFilterExclude {
		if !passSite {
			// Site excluded => set all samples.
			return nil, true
		}
		if smpl == nil {
			// Site-level expression passed: nothing to set.
			return nil, false
		}
		inv := make([]bool, len(smpl))
		anyLeft := false
		for i, b := range smpl {
			inv[i] = !b
			if inv[i] {
				anyLeft = true
			}
		}
		if !anyLeft {
			return nil, false
		}
		return inv, true
	}
	// Include logic.
	if !passSite {
		return nil, false
	}
	return smpl, true
}

// Destroy releases resources (none held).
func (p *setGTPlugin) Destroy() error { return nil }

// applyNewGTX transforms one genotype per the new-gt rule for sample i,
// resolving the -n X allele from xAllele (nil unless X is requested). It
// mirrors setGT.c process()'s per-sample apply chain (lines 553-567).
func (p *setGTPlugin) applyNewGTX(gt genotype, nals, major, minor int, xAllele []int, i int) (genotype, bool) {
	x := 0
	xMissing := false
	if xAllele != nil {
		if xAllele[i] == missingAllele {
			xMissing = true
		} else {
			x = xAllele[i]
		}
	}
	switch {
	case p.newMask&sngUnphased != 0:
		return unphaseGT(gt)
	case p.newMask == sngPhased:
		return phaseGT(gt)
	case p.newMask&sngCustom != 0:
		return p.applyCustom(gt, nals, major, minor, x, xMissing, xAllele != nil)
	case p.newMask&sngXVAF != 0:
		// -n X: set every allele to the per-sample largest-AD allele (unphased),
		// or missing when no AD entry was usable.
		return setAllAlleles(gt, x, xMissing, false)
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
// nals is REF+ALT count; x/xMissing carry the resolved -n X allele for this
// sample (used by an 'X' position inside the template) and xUsed reports
// whether X resolution ran at all. The custom template overrides ploidy.
func (p *setGTPlugin) applyCustom(gt genotype, nals, major, minor, x int, xMissing, xUsed bool) (genotype, bool) {
	out := genotype{}
	for i := 0; i < p.custom.ploidy; i++ {
		var newAllele int
		switch p.custom.alleles[i] {
		case cMinor:
			newAllele = minor
		case cMajor:
			newAllele = major
		case cXVAF:
			if !xUsed || xMissing {
				newAllele = nals // triggers the missing branch below
			} else {
				newAllele = x
			}
		case cMissing:
			newAllele = nals // triggers the missing branch below
		default:
			newAllele = p.custom.alleles[i]
		}
		if newAllele >= nals {
			// Requested index exceeds the alleles present: set missing (unphased).
			out.alleles = append(out.alleles, missingAllele)
			out.phased = append(out.phased, false)
		} else {
			out.alleles = append(out.alleles, newAllele)
			out.phased = append(out.phased, p.custom.phased[i])
		}
	}
	// Compare with the original to decide if anything actually changed.
	if sameGT(gt, out) {
		return gt, false
	}
	return out, true
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
