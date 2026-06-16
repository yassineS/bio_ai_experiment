// Native port of the upstream `tag2tag` plugin (plugins/tag2tag.c) for the
// common GL/PL/GP/GT conversions selected with --SRC-to-DST. These convert
// between the similar genotype-likelihood tags:
//
//	GL  : log10 likelihoods (Float, Number=G)
//	PL  : phred-scaled likelihoods (Integer, Number=G)
//	GP  : genotype probabilities (Float, Number=G)
//	GT  : hard-called genotype (String, Number=1)
//
// Options -r/--replace (drop the source tag) and -t/--threshold (GP->GT
// hard-call threshold) are supported, plus the FORMAT/QR,QA -> FORMAT/QS
// conversion (--QR-QA-to-QS), which concatenates the per-sample reference
// quality (QR, Number=1) and alternate quality sums (QA, Number=A) into the
// Number=R QS tag.
//
// The localized-allele expansion family (--LXX-to-XX, --LPL-to-PL, --LAD-to-AD)
// is supported: it expands the localized FORMAT tags (LAA + LPL/LAD) back into
// the standard Number=G PL and Number=R AD tags, using FORMAT/LAA to map each
// sample's localized indices to global allele indices, with -d/--defaults
// supplying the value placed in untouched cells and -s/--skip-nalt skipping
// sites above an allele threshold. This ports process_LXX() from tag2tag.c. The
// reverse direction (--XX-to-LXX, --PL-to-LPL, --AD-to-LAD) is unimplemented
// upstream too (process_XX is a stub that errors "todo"), so it is rejected with
// that same restriction.
package bcftools

import (
	"fmt"
	"math"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() { registerNativePlugin("tag2tag", func() NativePlugin { return &tag2tagPlugin{} }) }

// Localized-tag bitmask values, mirroring the enum tag in tag2tag.c (only the
// bits used by the supported LXX-to-XX conversions are needed).
const (
	t2tBitLAA = 1 << iota
	t2tBitLPL
	t2tBitLAD
	t2tBitPL
	t2tBitAD
)

// tag2tagPlugin implements the GL/PL/GP/GT conversions and the localized-allele
// expansion (LXX-to-XX). It is per-record and parallel.
type tag2tagPlugin struct {
	src, dst string
	dropSrc  bool
	gpThresh float64
	qrqa     bool // --QR-QA-to-QS mode

	// Localized-allele expansion state (--LXX-to-XX and friends). locExpand is
	// set when one of the supported localized-to-standard conversions was
	// selected; locSrc/locDst are bitmasks of the localized source tags consumed
	// and the standard destination tags produced. skipNalt and the dflt* values
	// mirror -s/--skip-nalt and -d/--defaults.
	locExpand bool
	locSrc    int
	locDst    int
	skipNalt  int
	dfltAD    string // "." (missing) by default, else the integer text
	dfltPL    string
}

// Name returns the plugin name.
func (p *tag2tagPlugin) Name() string { return "tag2tag" }

// About returns the one-line description, matching tag2tag.c about().
func (p *tag2tagPlugin) About() string {
	return "Convert between similar tags, such as GL,PL,GP or QR,QA,QS or localized alleles, eg PL,LPL."
}

// Parallel reports true.
func (p *tag2tagPlugin) Parallel() bool { return true }

// Init parses the --SRC-to-DST selector and -r/-t/-s/-d options.
func (p *tag2tagPlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	// Upstream tag2tag.c leaves gp_th calloc'd to 0 by default (unlike the
	// usage text's "[0.1]"), so the GP->GT hard call requires GP==1 unless -t
	// is given explicitly.
	p.gpThresh = 0
	// Defaults for --LXX-to-XX cells, mirroring args->dflt_AD/dflt_PL =
	// bcf_int32_missing (rendered "." in the VCF text model).
	p.dfltAD = "."
	p.dfltPL = "."
	for i := 0; i < len(args); i++ {
		a := args[i]
		switch {
		case a == "-r" || a == "--replace":
			p.dropSrc = true
		case a == "-t" || a == "--threshold":
			if i+1 >= len(args) {
				return nil, fmt.Errorf("tag2tag: -t requires an argument")
			}
			i++
			v, err := strconv.ParseFloat(args[i], 64)
			if err != nil || v < 0 || v > 1 {
				return nil, fmt.Errorf("tag2tag: expected value between 0 and 1 for -t, got %q", args[i])
			}
			p.gpThresh = v
		case a == "-s" || a == "--skip-nalt":
			if i+1 >= len(args) {
				return nil, fmt.Errorf("tag2tag: -s requires an argument")
			}
			i++
			n, err := strconv.Atoi(args[i])
			if err != nil || n < 0 {
				return nil, fmt.Errorf("tag2tag: could not parse: --skip-nalt %s", args[i])
			}
			p.skipNalt = n
		case a == "-d" || a == "--defaults":
			if i+1 >= len(args) {
				return nil, fmt.Errorf("tag2tag: -d requires an argument")
			}
			i++
			if err := p.parseDefaults(args[i]); err != nil {
				return nil, err
			}
		case strings.HasPrefix(a, "--") && strings.Contains(strings.ToUpper(a), "-TO-"):
			if err := p.parseSelector(a); err != nil {
				return nil, err
			}
		default:
			return nil, fmt.Errorf("tag2tag: unsupported option %q", a)
		}
	}
	if p.locExpand {
		return p.initLocalized(hdr)
	}
	if p.qrqa {
		for _, t := range []string{"QR", "QA"} {
			if !hasFormatHeader(hdr.MetaInfo, t) {
				return nil, fmt.Errorf("tag2tag: the source tag does not exist: %s", t)
			}
		}
		out := &vcf.Header{Samples: hdr.Samples}
		out.MetaInfo = append(out.MetaInfo, hdr.MetaInfo...)
		if p.dropSrc {
			out.MetaInfo = removeFormatHeader(out.MetaInfo, "QR")
			out.MetaInfo = removeFormatHeader(out.MetaInfo, "QA")
		}
		out.MetaInfo = appendInfoHeader(out.MetaInfo, `##FORMAT=<ID=QS,Number=R,Type=Integer,Description="Phred-score allele quality sum">`)
		return out, nil
	}

	if p.src == "" {
		return nil, fmt.Errorf("tag2tag: a conversion such as --PL-to-GT must be given")
	}
	if !hasFormatHeader(hdr.MetaInfo, p.src) {
		return nil, fmt.Errorf("tag2tag: the source tag does not exist: %s", p.src)
	}

	out := &vcf.Header{Samples: hdr.Samples}
	out.MetaInfo = append(out.MetaInfo, hdr.MetaInfo...)
	if p.dropSrc {
		out.MetaInfo = removeFormatHeader(out.MetaInfo, p.src)
	}
	out.MetaInfo = appendInfoHeader(out.MetaInfo, dstHeaderLine(p.dst))
	return out, nil
}

// parseSelector parses a --SRC-to-DST option into src/dst (or the localized /
// QR-QA modes), mirroring parse_ori2new_option() in tag2tag.c.
func (p *tag2tagPlugin) parseSelector(opt string) error {
	up := strings.ToUpper(opt[2:]) // strip leading --
	if up == "QR-QA-TO-QS" {
		p.qrqa = true
		return nil
	}
	// Localized-to-standard expansions (the supported direction): set the
	// loc_src/loc_dst masks exactly as parse_ori2new_option does.
	switch up {
	case "LXX-TO-XX":
		p.locExpand = true
		p.locSrc = t2tBitLPL | t2tBitLAD | t2tBitLAA
		p.locDst = t2tBitPL | t2tBitAD
		return nil
	case "LPL-TO-PL":
		p.locExpand = true
		p.locSrc = t2tBitLPL | t2tBitLAA
		p.locDst = t2tBitPL
		return nil
	case "LAD-TO-AD":
		p.locExpand = true
		p.locSrc = t2tBitLAD | t2tBitLAA
		p.locDst = t2tBitAD
		return nil
	case "XX-TO-LXX", "PL-TO-LPL", "AD-TO-LAD":
		// Upstream's process_XX is a stub: error("todo: --XX-to-LXX\n"). Reject
		// with the same restriction rather than silently diverging.
		return fmt.Errorf("tag2tag: todo: --XX-to-LXX")
	}
	parts := strings.Split(up, "-TO-")
	if len(parts) != 2 {
		return fmt.Errorf("tag2tag: could not parse conversion %q", opt)
	}
	src, dst := parts[0], parts[1]
	valid := map[string]bool{"GL": true, "PL": true, "GP": true, "GT": true}
	if !valid[src] || !valid[dst] || src == dst {
		return fmt.Errorf("tag2tag: the conversion is not supported: %s", opt)
	}
	if src == "GT" {
		return fmt.Errorf("tag2tag: GT cannot be a source tag")
	}
	p.src, p.dst = src, dst
	return nil
}

// Process converts the source tag to the destination tag for one record.
func (p *tag2tagPlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	if p.locExpand {
		return p.processLocalized(v)
	}
	if p.qrqa {
		return p.processQRQA(v)
	}
	if !formatHasTag(v, p.src) {
		return []*vcf.Variant{v}, nil
	}
	nals := 1 + len(v.Alt)

	for i := range v.Samples {
		raw, ok := v.Samples[i].Data[p.src]
		if !ok || raw == "." || raw == "" {
			continue
		}
		// Convert the per-sample source values to GL (log10 likelihoods).
		gl, ok := p.toGL(raw)
		if !ok {
			continue
		}
		switch p.dst {
		case "GL":
			v.Samples[i].Data["GL"] = formatGLList(gl)
		case "PL":
			v.Samples[i].Data["PL"] = glToPL(gl)
		case "GP":
			v.Samples[i].Data["GP"] = glToGP(gl)
		case "GT":
			v.Samples[i].Data["GT"] = glToGT(gl, nals, p.gpThresh)
		}
	}

	for _, tag := range []string{"GL", "PL", "GP", "GT"} {
		if tag == p.dst {
			ensureFormatTag(v, tag)
		}
	}
	if p.dropSrc {
		removeFormatTag(v, p.src)
	}
	return []*vcf.Variant{v}, nil
}

// parseDefaults parses the -d/--defaults LIST, e.g. "AD:0,PL:.", setting the
// per-tag placeholder used by --LXX-to-XX for cells with no localized source
// value. A "." selects the missing value. It ports parse_defaults() in
// tag2tag.c (which scans for "AD:"/"PL:" prefixes).
func (p *tag2tagPlugin) parseDefaults(opt string) error {
	ptr := opt
	for ptr != "" {
		var dst *string
		switch {
		case strings.HasPrefix(strings.ToUpper(ptr), "AD:"):
			dst = &p.dfltAD
		case strings.HasPrefix(strings.ToUpper(ptr), "PL:"):
			dst = &p.dfltPL
		default:
			// Upstream's while loop only advances on a recognised prefix; an
			// unrecognised character would loop forever. Treat it as a parse
			// error to fail fast instead.
			return fmt.Errorf("tag2tag: could not parse: --defaults %s", opt)
		}
		ptr = ptr[3:]
		// Read up to the next comma.
		end := strings.IndexByte(ptr, ',')
		var tok string
		if end < 0 {
			tok, ptr = ptr, ""
		} else {
			tok, ptr = ptr[:end], ptr[end+1:]
		}
		if tok == "." || tok == "" {
			*dst = "."
		} else {
			if _, err := strconv.Atoi(tok); err != nil {
				return fmt.Errorf("tag2tag: could not parse: --defaults %s", opt)
			}
			*dst = tok
		}
	}
	return nil
}

// initLocalized validates the localized source tags and builds the output
// header for the --LXX-to-XX family, mirroring the LXX branch of init().
func (p *tag2tagPlugin) initLocalized(hdr *vcf.Header) (*vcf.Header, error) {
	// All requested localized source tags must be declared FORMAT fields.
	for _, t := range []struct {
		bit  int
		name string
	}{{t2tBitLPL, "LPL"}, {t2tBitLAD, "LAD"}, {t2tBitLAA, "LAA"}} {
		if p.locSrc&t.bit != 0 && !hasFormatHeader(hdr.MetaInfo, t.name) {
			return nil, fmt.Errorf("tag2tag: the source tag does not exist: %s", t.name)
		}
	}

	out := &vcf.Header{Samples: hdr.Samples}
	out.MetaInfo = append(out.MetaInfo, hdr.MetaInfo...)

	// With -r/--replace remove the consumed source headers, but NOT when
	// -s/--skip-nalt is set (some records may retain the tags). LAA is removed
	// only when it is the sole remaining loc_src bit (drop_laa == 1<<LAA).
	if p.dropSrc && p.skipNalt == 0 {
		dropLAA := p.locSrc
		if p.locSrc&t2tBitLAD != 0 {
			out.MetaInfo = removeFormatHeader(out.MetaInfo, "LAD")
			dropLAA &^= t2tBitLAD
		}
		if p.locSrc&t2tBitLPL != 0 {
			out.MetaInfo = removeFormatHeader(out.MetaInfo, "LPL")
			dropLAA &^= t2tBitLPL
		}
		if dropLAA == t2tBitLAA {
			out.MetaInfo = removeFormatHeader(out.MetaInfo, "LAA")
		}
	}

	// Append the destination headers (AD before PL, matching tags_XX = {PL,AD}
	// header order is PL then AD, but the per-record update order is AD then PL;
	// htslib appends PL header first, then AD — reproduce that header order).
	if p.locDst&t2tBitPL != 0 {
		out.MetaInfo = appendInfoHeader(out.MetaInfo, `##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods">`)
	}
	if p.locDst&t2tBitAD != 0 {
		out.MetaInfo = appendInfoHeader(out.MetaInfo, `##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">`)
	}
	return out, nil
}

// processLocalized expands the localized tags into AD (Number=R) and/or PL
// (Number=G) for one record, porting process_LXX() in tag2tag.c. It is skipped
// for sites above the -s/--skip-nalt allele threshold, and for records without
// FORMAT/LAA.
func (p *tag2tagPlugin) processLocalized(v *vcf.Variant) ([]*vcf.Variant, error) {
	nals := 1 + len(v.Alt)
	if p.skipNalt != 0 && nals > p.skipNalt {
		return []*vcf.Variant{v}, nil
	}
	if !formatHasTag(v, "LAA") {
		return []*vcf.Variant{v}, nil
	}

	// Per-sample localized allele lists. A missing LAA ("." or empty) yields an
	// empty list (only REF is localized).
	laa := make([][]int, len(v.Samples))
	for i := range v.Samples {
		laa[i] = parseLAAList(v.Samples[i].Data["LAA"])
	}

	dropLAA := p.locSrc

	if p.locSrc&t2tBitLAD != 0 && formatHasTag(v, "LAD") {
		for i := range v.Samples {
			src := splitCSV(v.Samples[i].Data["LAD"])
			v.Samples[i].Data["AD"] = ladToAD(src, laa[i], nals, p.dfltAD)
		}
		ensureFormatTag(v, "AD")
		if p.dropSrc {
			removeFormatTag(v, "LAD")
			dropLAA &^= t2tBitLAD
		}
	}

	if p.locSrc&t2tBitLPL != 0 && formatHasTag(v, "LPL") {
		for i := range v.Samples {
			src := splitCSV(v.Samples[i].Data["LPL"])
			v.Samples[i].Data["PL"] = lplToPL(src, laa[i], nals, p.dfltPL)
		}
		ensureFormatTag(v, "PL")
		if p.dropSrc {
			removeFormatTag(v, "LPL")
			dropLAA &^= t2tBitLPL
		}
	}

	if p.dropSrc && dropLAA == t2tBitLAA {
		removeFormatTag(v, "LAA")
	}
	return []*vcf.Variant{v}, nil
}

// ladToAD expands a localized allelic-depth list into a Number=R AD vector,
// porting the LAD->AD block of process_LXX. dst[0] is the REF depth (src[0]);
// dst[1..nals-1] default to dflt; for each localized allele laa[j-1], dst at
// that global index takes src[j].
func ladToAD(src []string, laa []int, nals int, dflt string) string {
	dst := make([]string, nals)
	for j := range dst {
		dst[j] = dflt
	}
	if len(src) > 0 {
		dst[0] = src[0]
	}
	for j := 1; j < len(src); j++ {
		if j-1 >= len(laa) {
			break
		}
		a := laa[j-1]
		if a >= 0 && a < nals {
			dst[a] = src[j]
		}
	}
	return strings.Join(dst, ",")
}

// lplToPL expands a localized phred-likelihood list into a Number=G PL vector,
// porting the LPL->PL block of process_LXX. tmp_laa = [0, laa...]; dst has
// nals*(nals+1)/2 cells defaulting to dflt; for each pair (j>=k) over tmp_laa
// with tmp_laa[j] in [0,nals), dst[ tmp_laa[j]*(tmp_laa[j]+1)/2 + tmp_laa[k] ]
// consumes the next src value in order.
func lplToPL(src []string, laa []int, nals int, dflt string) string {
	ndst := nals * (nals + 1) / 2
	dst := make([]string, ndst)
	for j := range dst {
		dst[j] = dflt
	}
	tmpLAA := make([]int, 0, len(laa)+1)
	tmpLAA = append(tmpLAA, 0)
	tmpLAA = append(tmpLAA, laa...)
	si := 0
	for j := 0; j < len(tmpLAA); j++ {
		if !(tmpLAA[j] >= 0 && tmpLAA[j] < nals) {
			break
		}
		for k := 0; k <= j; k++ {
			idx := tmpLAA[j]*(tmpLAA[j]+1)/2 + tmpLAA[k]
			if idx >= 0 && idx < ndst && si < len(src) {
				dst[idx] = src[si]
			}
			si++
		}
	}
	return strings.Join(dst, ",")
}

// parseLAAList parses a FORMAT/LAA value into the per-sample list of localized
// allele indices, dropping a missing ("." / "") field to an empty list and
// stopping at the first missing element (htslib's vector-end semantics).
func parseLAAList(s string) []int {
	if s == "" || s == "." {
		return nil
	}
	out := make([]int, 0, 4)
	for _, tok := range strings.Split(s, ",") {
		if tok == "." || tok == "" {
			break
		}
		n, err := strconv.Atoi(tok)
		if err != nil {
			break
		}
		out = append(out, n)
	}
	return out
}

// splitCSV splits a per-sample value into its comma-separated tokens, returning
// an empty slice for a missing ("." / "") field.
func splitCSV(s string) []string {
	if s == "" || s == "." {
		return nil
	}
	return strings.Split(s, ",")
}

// processQRQA implements --QR-QA-to-QS: it concatenates FORMAT/QR (Number=1)
// and FORMAT/QA (Number=A) into FORMAT/QS (Number=R) per sample. The record is
// left unchanged if QR is absent, or (for multiallelic sites) if QA is absent.
func (p *tag2tagPlugin) processQRQA(v *vcf.Variant) ([]*vcf.Variant, error) {
	if !formatHasTag(v, "QR") {
		return []*vcf.Variant{v}, nil
	}
	nals := 1 + len(v.Alt)
	if nals > 1 && !formatHasTag(v, "QA") {
		return []*vcf.Variant{v}, nil
	}
	qs := make([]string, len(v.Samples))
	for i := range v.Samples {
		qr := v.Samples[i].Data["QR"]
		if nals == 1 {
			qs[i] = qr
			continue
		}
		qa := v.Samples[i].Data["QA"]
		qs[i] = qr + "," + qa
	}
	ensureFormatTag(v, "QS")
	for i := range v.Samples {
		v.Samples[i].Data["QS"] = qs[i]
	}
	if p.dropSrc {
		removeFormatTag(v, "QR")
		removeFormatTag(v, "QA")
	}
	return []*vcf.Variant{v}, nil
}

// Destroy releases resources (none held).
func (p *tag2tagPlugin) Destroy() error { return nil }

// toGL parses the source value and returns it as GL (log10 likelihoods),
// mirroring the "convert to GL" step in tag2tag.c process(). GP values are
// log10'd (0 -> -99); PL values are scaled by -0.1; GL is used as-is.
func (p *tag2tagPlugin) toGL(raw string) ([]float64, bool) {
	parts := strings.Split(raw, ",")
	out := make([]float64, len(parts))
	for i, s := range parts {
		if s == "." {
			out[i] = math.NaN()
			continue
		}
		f, err := strconv.ParseFloat(s, 64)
		if err != nil {
			return nil, false
		}
		switch p.src {
		case "GP":
			if f != 0 {
				out[i] = math.Log10(f)
			} else {
				out[i] = -99
			}
		case "PL":
			out[i] = -0.1 * f
		case "GL":
			out[i] = f
		}
	}
	return out, true
}

// formatGLList renders GL values, matching htslib's float %g (6 sig digits).
func formatGLList(gl []float64) string {
	parts := make([]string, len(gl))
	for i, x := range gl {
		if math.IsNaN(x) {
			parts[i] = "."
		} else {
			parts[i] = formatVCFFloat(x)
		}
	}
	return strings.Join(parts, ",")
}

// glToPL converts GL back to PL: PL = round(-10*GL), mirroring lroundf(-10*gl).
func glToPL(gl []float64) string {
	parts := make([]string, len(gl))
	for i, x := range gl {
		if math.IsNaN(x) {
			parts[i] = "."
		} else {
			parts[i] = strconv.Itoa(int(math.Round(-10 * x)))
		}
	}
	return strings.Join(parts, ",")
}

// glToGP converts GL to normalized GP: 10^GL then divide by the sum.
func glToGP(gl []float64) string {
	gp := make([]float64, len(gl))
	sum := 0.0
	for i, x := range gl {
		if math.IsNaN(x) {
			gp[i] = math.NaN()
			continue
		}
		gp[i] = math.Pow(10, x)
		sum += gp[i]
	}
	parts := make([]string, len(gl))
	for i := range gp {
		if math.IsNaN(gp[i]) {
			parts[i] = "."
		} else if sum > 0 {
			parts[i] = formatVCFFloat(gp[i] / sum)
		} else {
			parts[i] = formatVCFFloat(gp[i])
		}
	}
	return strings.Join(parts, ",")
}

// glToGT hard-calls a genotype from GL, mirroring the GP->GT branch of
// tag2tag.c: convert GL to normalized GP, pick the max, and apply the
// threshold (1 - gpThresh) below which the call becomes missing.
func glToGT(gl []float64, nals int, gpThresh float64) string {
	// Convert to GP (normalized 10^GL). htslib stores the likelihoods as C
	// floats, so the pow/sum/divide all happen in 32-bit precision; doing the
	// same here is what makes a near-certain call (e.g. PL 120,0,130) round to
	// GP==1 and get hard-called rather than falling just under the threshold.
	gp := make([]float64, len(gl))
	var sum float32
	for i, x := range gl {
		if math.IsNaN(x) {
			gp[i] = math.NaN()
			continue
		}
		f := float32(math.Pow(10, x))
		gp[i] = float64(f)
		sum += f
	}
	if sum > 0 {
		for i := range gp {
			if !math.IsNaN(gp[i]) {
				gp[i] = float64(float32(gp[i]) / sum)
			}
		}
	}
	if len(gp) == 0 || math.IsNaN(gp[0]) {
		return "./."
	}
	jmax := 0
	n1 := len(gp)
	for j := 1; j < n1; j++ {
		if math.IsNaN(gp[j]) {
			n1 = j
			break
		}
		if gp[j] > gp[jmax] {
			jmax = j
		}
	}

	// Haploid genotype: number of GP values == number of alleles.
	if n1 == nals {
		if gp[jmax] < 1-gpThresh {
			return "."
		}
		return strconv.Itoa(jmax)
	}
	if gp[jmax] < 1-gpThresh {
		return "./."
	}
	if jmax == 0 {
		return "0/0"
	}
	a, b := gt2alleles(jmax)
	return strconv.Itoa(a) + "/" + strconv.Itoa(b)
}

// gt2alleles maps the diploid GP/PL index back to its two alleles, mirroring
// htslib bcf_gt2alleles: index = b*(b+1)/2 + a with a <= b.
func gt2alleles(idx int) (int, int) {
	b := 0
	for (b+1)*(b+2)/2 <= idx {
		b++
	}
	a := idx - b*(b+1)/2
	return a, b
}

// dstHeaderLine returns the ##FORMAT header line for a destination tag.
func dstHeaderLine(tag string) string {
	switch tag {
	case "GP":
		return `##FORMAT=<ID=GP,Number=G,Type=Float,Description="Genotype probabilities">`
	case "GL":
		return `##FORMAT=<ID=GL,Number=G,Type=Float,Description="Genotype likelihoods">`
	case "PL":
		return `##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods">`
	case "GT":
		return `##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">`
	}
	return ""
}

// removeFormatHeader drops the ##FORMAT line for id from meta.
func removeFormatHeader(meta []string, id string) []string {
	out := meta[:0:0]
	for _, m := range meta {
		if headerKind(m) == "##FORMAT" && headerID(m) == id {
			continue
		}
		out = append(out, m)
	}
	return out
}

// removeFormatTag drops a FORMAT tag (and its per-sample values) from v.
func removeFormatTag(v *vcf.Variant, tag string) {
	idx := -1
	for i, f := range v.Format {
		if f == tag {
			idx = i
			break
		}
	}
	if idx < 0 {
		return
	}
	v.Format = append(v.Format[:idx], v.Format[idx+1:]...)
	for i := range v.Samples {
		delete(v.Samples[i].Data, tag)
	}
}
