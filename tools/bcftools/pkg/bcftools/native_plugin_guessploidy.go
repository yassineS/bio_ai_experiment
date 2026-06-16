// Native port of the upstream `guess-ploidy` plugin (plugins/guess-ploidy.c) for
// its whole-file (no region/genome jump) mode. It guesses each sample's sex by
// comparing, under Hardy-Weinberg equilibrium, the haploid and diploid
// likelihoods of the observed genotype likelihoods (PL/GL) or genotypes (GT)
// across the streamed sites — intended for the non-PAR region of chrX. The
// VCF/BCF output is suppressed; a per-sample table is printed.
//
// The -r/-R region selection is handled by the shared host region filter, and
// the -g/--genome shortcut expands to the equivalent -r CHR:BEG-END region for
// the four built-in genome presets (b37, b38, hg19, hg38) via RewriteArgs,
// mirroring guess-ploidy.c's `case 'g'` which simply sets `region` to the
// hardcoded non-PAR span of chrX for that genome. The filter expressions
// (--include/--exclude; note guess-ploidy maps -i to --include-indels and -e to
// --error-rate, so the filter is long-form only) are supported as a
// site/per-sample pre-filter via the shared filter engine, matching
// guess-ploidy.c's filter_init/filter_test usage. The default whole-file scan,
// the common case when piping `bcftools view -r chrX:... | bcftools
// +guess-ploidy`, is implemented.
package bcftools

import (
	"fmt"
	"io"
	"math"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() {
	registerNativePlugin("guess-ploidy", func() NativePlugin { return &guessPloidyPlugin{} })
}

// guess-ploidy tag-source selectors, matching GUESS_* in guess-ploidy.c.
const (
	guessGT = 1
	guessPL = 2
	guessGL = 4
)

// guessPloidyCount accumulates one sample's running log-likelihoods, mirroring
// count_t in guess-ploidy.c.
type guessPloidyCount struct {
	ncount     uint64
	phap, pdip float64
}

// guessPloidyPlugin implements the `guess-ploidy` plugin in whole-file mode.
type guessPloidyPlugin struct {
	hdr       *vcf.Header
	tag       int
	afTag     string
	afDflt    float64
	gtErrProb float64
	verbose   int
	indels    bool

	pl2p   []float64
	counts []guessPloidyCount

	filter *pluginFilter // compiled --include/--exclude pre-filter, nil if none

	out    io.Writer
	stderr io.Writer
	argv   []string
}

// SuppressVCF reports true: guess-ploidy emits only its textual report.
func (p *guessPloidyPlugin) SuppressVCF() bool { return true }

// SetStdout wires the host stdout writer the report is printed to.
func (p *guessPloidyPlugin) SetStdout(w io.Writer) { p.out = w }

// SetStderr wires the host stderr writer the PL/GL fallback warnings use.
func (p *guessPloidyPlugin) SetStderr(w io.Writer) { p.stderr = w }

// SetArgv records the upstream-equivalent argv for the verbose command-line
// header line.
func (p *guessPloidyPlugin) SetArgv(argv []string) { p.argv = argv }

// RunStyle reports that guess-ploidy is a run()-style plugin.
func (p *guessPloidyPlugin) RunStyle() bool { return true }

// FlagTakesValue reports whether one of guess-ploidy's value-taking flags
// consumes the following token.
func (p *guessPloidyPlugin) FlagTakesValue(flag string) bool {
	switch flag {
	case "--AF-tag", "--AF-dflt", "--exclude", "--include", "-e", "--error-rate",
		"-t", "--tag", "-g", "--genome", "-r", "--regions", "-R", "--regions-file":
		return true
	case "-v", "--verbosity", "--verbose":
		// Optional argument: do not unconditionally consume the next token.
		return false
	}
	return false
}

// guessPloidyGenomePresets maps the -g/--genome shortcut values to the
// hardcoded non-PAR chrX region upstream selects (guess-ploidy.c `case 'g'`).
// The coordinates are the GRCh37/GRCh38 chrX non-PAR boundaries; hg19/hg38 use
// the same coordinates with the "chr" contig prefix.
var guessPloidyGenomePresets = map[string]string{
	"b37":  "X:2699521-154931043",
	"b38":  "X:2781480-155701381",
	"hg19": "chrX:2699521-154931043",
	"hg38": "chrX:2781480-155701381",
}

// guessPloidyGenomeRegion returns the non-PAR chrX region for a -g/--genome
// preset value (case-insensitive, matching upstream's strcasecmp), or ok=false
// when the value is not one of b37/b38/hg19/hg38.
func guessPloidyGenomeRegion(value string) (region string, ok bool) {
	region, ok = guessPloidyGenomePresets[strings.ToLower(value)]
	return region, ok
}

// RewriteArgs expands -g/--genome PRESET into -r REGION before the host's
// region/target extraction runs, so the shared region filter restricts the
// stream to the genome's non-PAR chrX span exactly as upstream's `-g` shortcut
// does. An unrecognised preset is rejected with upstream's message.
func (p *guessPloidyPlugin) RewriteArgs(args []string) ([]string, error) {
	out := make([]string, 0, len(args))
	for i := 0; i < len(args); i++ {
		a := args[i]
		var value string
		switch {
		case a == "-g" || a == "--genome":
			if i+1 >= len(args) {
				return nil, fmt.Errorf("guess-ploidy: %s requires an argument", a)
			}
			i++
			value = args[i]
		case strings.HasPrefix(a, "-g") && len(a) > 2:
			value = a[2:]
		case strings.HasPrefix(a, "--genome="):
			value = a[len("--genome="):]
		default:
			out = append(out, a)
			continue
		}
		region, ok := guessPloidyGenomeRegion(value)
		if !ok {
			return nil, fmt.Errorf("guess-ploidy: the argument not recognised, expected --genome b37, b38, hg19 or hg38: %s", value)
		}
		out = append(out, "-r", region)
	}
	return out, nil
}

// RegionTargetCaps opts guess-ploidy into -r/-R region selection only. Its -t is
// --tag (the INFO field to read genotype likelihoods from), NOT targets, so the
// shared filter must leave -t/-T to guess-ploidy's own parser.
func (p *guessPloidyPlugin) RegionTargetCaps() regionTargetCaps {
	return regionsOnlyCaps
}

// Name returns the plugin name.
func (p *guessPloidyPlugin) Name() string { return "guess-ploidy" }

// About returns the one-line description, matching guess-ploidy.c about().
func (p *guessPloidyPlugin) About() string {
	return "Determine sample sex by checking genotype likelihoods in haploid regions.\n"
}

// Parallel reports false: per-sample likelihoods accumulate serially.
func (p *guessPloidyPlugin) Parallel() bool { return false }

// Init parses the options and rejects the region/filter modes, then validates
// the requested tag against the header (with the PL->GL->GT fallback warnings).
func (p *guessPloidyPlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	p.tag = guessPL
	p.gtErrProb = 1e-3
	p.afDflt = 0.5
	var filterExpr string
	var filterExclude, haveFilter bool
	for i := 0; i < len(args); i++ {
		a := args[i]
		next := func() (string, error) {
			if i+1 >= len(args) {
				return "", fmt.Errorf("guess-ploidy: %s requires an argument", a)
			}
			i++
			return args[i], nil
		}
		switch a {
		case "--AF-tag":
			v, err := next()
			if err != nil {
				return nil, err
			}
			p.afTag = v
		case "--AF-dflt":
			v, err := next()
			if err != nil {
				return nil, err
			}
			f, err := strconv.ParseFloat(v, 64)
			if err != nil {
				return nil, fmt.Errorf("guess-ploidy: could not parse: --AF-dflt %s", v)
			}
			p.afDflt = f
		case "--exclude", "--include":
			if haveFilter {
				return nil, fmt.Errorf("guess-ploidy: only one --include or --exclude expression can be given, and they cannot be combined")
			}
			v, err := next()
			if err != nil {
				return nil, err
			}
			filterExpr = v
			filterExclude = a == "--exclude"
			haveFilter = true
		case "-e", "--error-rate":
			v, err := next()
			if err != nil {
				return nil, err
			}
			f, err := strconv.ParseFloat(v, 64)
			if err != nil || f < 0 || f > 1 {
				return nil, fmt.Errorf("guess-ploidy: expected value from the interval [0,1]: -e %s", v)
			}
			p.gtErrProb = f
		case "-i", "--include-indels":
			p.indels = true
		case "-g", "--genome":
			// -g/--genome is rewritten to -r REGION by RewriteArgs before the
			// host's region extraction runs, so by the time Init sees the argv
			// the token is gone. Reaching it here means the rewrite consumed the
			// value but left the flag (it cannot happen); reject defensively.
			if _, err := next(); err != nil {
				return nil, err
			}
			return nil, fmt.Errorf("guess-ploidy: internal: -g/--genome was not rewritten to -r")
		case "-t", "--tag":
			v, err := next()
			if err != nil {
				return nil, err
			}
			switch strings.ToUpper(v) {
			case "GT":
				p.tag = guessGT
			case "PL":
				p.tag = guessPL
			case "GL":
				p.tag = guessGL
			default:
				return nil, fmt.Errorf("guess-ploidy: expected --tag GT, PL or GL: %s", v)
			}
		case "-v", "--verbosity", "--verbose":
			// Optional argument: a bare -v increments the level (getopt's
			// optional-argument handling). A separate following token is NOT
			// consumed by upstream (it stays a positional); the attached forms
			// (-v2, --verbosity=2) are handled in the default branch below.
			p.verbose++
		default:
			if lvl, ok := parseAttachedVerbosity(a); ok {
				if lvl < 0 {
					return nil, fmt.Errorf("guess-ploidy: could not parse argument: %s", a)
				}
				p.verbose = lvl
				continue
			}
			return nil, fmt.Errorf("guess-ploidy: unsupported option %q", a)
		}
	}
	if haveFilter {
		f, err := newPluginFilterWithHeader(filterExpr, filterExclude, hdr)
		if err != nil {
			return nil, fmt.Errorf("guess-ploidy: %w", err)
		}
		p.filter = f
	}
	p.hdr = hdr

	if p.afTag != "" && !hasInfoHeaderMeta(hdr.MetaInfo, p.afTag) {
		return nil, fmt.Errorf("guess-ploidy: no such INFO tag: %s", p.afTag)
	}
	if p.tag&guessPL != 0 && !hasFormatHeader(hdr.MetaInfo, "PL") {
		if p.stderr != nil {
			fmt.Fprint(p.stderr, "Warning: PL tag not found in header, switching to GL\n")
		}
		p.tag = guessGL
	}
	if p.tag&guessGL != 0 && !hasFormatHeader(hdr.MetaInfo, "GL") {
		if p.stderr != nil {
			fmt.Fprint(p.stderr, "Warning: GL tag not found in header, switching to GT\n")
		}
		p.tag = guessGT
	}
	if p.tag&guessGT != 0 && !hasFormatHeader(hdr.MetaInfo, "GT") {
		return nil, fmt.Errorf("guess-ploidy: GT tag not found in header")
	}

	if p.tag&guessPL != 0 {
		p.pl2p = make([]float64, 256)
		for i := 0; i < 256; i++ {
			p.pl2p[i] = math.Pow(10, -float64(i)/10)
		}
	}
	p.counts = make([]guessPloidyCount, len(hdr.Samples))

	if p.verbose > 0 && p.out != nil {
		fmt.Fprint(p.out, "# This file was produced by: bcftools +guess-ploidy(bio_ai_experiment+htslib-bio_ai_experiment)\n")
		fmt.Fprintf(p.out, "# The command line was:\tbcftools +%s\n", strings.Join(p.argv, " "))
		fmt.Fprint(p.out, "# [1]SEX\t[2]Sample\t[3]Predicted sex\t[4]log P(Haploid)/nSites\t[5]log P(Diploid)/nSites\t[6]nSites\t[7]Score: F < 0 < M ($4-$5)\n")
		if p.verbose > 1 {
			fmt.Fprint(p.out, "# [1]DBG\t[2]Chr\t[3]Pos\t[4]Sample\t[5]AF\t[6]pRR\t[7]pRA\t[8]pAA\t[9]P(Haploid)\t[10]P(Diploid)\n")
		}
	}
	return hdr, nil
}

// Process accumulates the per-sample haploid/diploid log-likelihoods for one
// record, mirroring process_region_guess() in guess-ploidy.c.
func (p *guessPloidyPlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	nAllele := len(v.Alt) + 1
	if nAllele == 1 {
		return nil, nil
	}
	if !p.indels && variantTypeMask(v)&vtSNP == 0 {
		return nil, nil
	}
	// Apply the -i/-e filter as a site/per-sample gate. A nil mask means the
	// expression is site-level (all samples accumulate when the site passes);
	// a non-nil mask drops the non-matching samples, mirroring guess-ploidy.c's
	// smpl_pass() guard in every per-sample loop.
	passSite, smplPass := p.filter.testSamples(v)
	if !passSite {
		return nil, nil
	}
	nsmp := len(v.Samples)
	freq := [2]float64{0, 0}
	tmp := make([][3]float64, nsmp)

	switch {
	case p.tag&guessGT != 0:
		if !formatHasTag(v, "GT") {
			return nil, nil
		}
		for s := 0; s < nsmp; s++ {
			if smplPass != nil && !smplPass[s] {
				tmp[s][0] = -1
				continue
			}
			gt, ok := sampleGT(v, s)
			if !ok || len(gt.alleles) == 0 || gt.alleles[0] == missingAllele {
				tmp[s][0] = -1
				continue
			}
			e := p.gtErrProb
			if len(gt.alleles) == 1 { // haploid
				if gt.alleles[0] == 0 {
					tmp[s][0] = 1 - 2*e
					tmp[s][1], tmp[s][2] = e, e
				} else {
					tmp[s][0], tmp[s][1] = e, e
					tmp[s][2] = 1 - 2*e
				}
				continue
			}
			a0, a1 := gt.alleles[0], gt.alleles[1]
			if a0 == 0 && a1 == 0 {
				tmp[s][0] = 1 - 2*e
				tmp[s][1], tmp[s][2] = e, e
			} else if a0 == a1 {
				tmp[s][0], tmp[s][1] = e, e
				tmp[s][2] = 1 - 2*e
			} else {
				tmp[s][1] = 1 - 2*e
				tmp[s][0], tmp[s][2] = e, e
			}
			freq[0] += 2*tmp[s][0] + tmp[s][1]
			freq[1] += tmp[s][1] + 2*tmp[s][2]
		}
	case p.tag&guessPL != 0:
		if !p.accumulatePL(v, nAllele, nsmp, tmp, &freq, smplPass) {
			return nil, nil
		}
	default: // GL
		if !p.accumulateGL(v, nAllele, nsmp, tmp, &freq, smplPass) {
			return nil, nil
		}
	}

	if p.afTag != "" {
		if af, ok := infoFloat0(v, p.afTag); ok {
			freq[0] = 1 - af
			freq[1] = af
		}
	}
	if freq[0] == 0 && freq[1] == 0 {
		freq[0] = 1 - p.afDflt
		freq[1] = p.afDflt
	}
	sum := freq[0] + freq[1]
	freq[0] /= sum
	freq[1] /= sum

	for s := 0; s < nsmp; s++ {
		if smplPass != nil && !smplPass[s] {
			continue
		}
		if tmp[s][0] < 0 {
			continue
		}
		phap := freq[0]*tmp[s][0] + freq[1]*tmp[s][2]
		pdip := freq[0]*freq[0]*tmp[s][0] + 2*freq[0]*freq[1]*tmp[s][1] + freq[1]*freq[1]*tmp[s][2]
		p.counts[s].phap += math.Log(phap)
		p.counts[s].pdip += math.Log(pdip)
		p.counts[s].ncount++
		if p.verbose > 1 && p.out != nil {
			fmt.Fprintf(p.out, "DBG\t%s\t%d\t%s\t%e\t%e\t%e\t%e\t%e\t%e\n",
				v.Chrom, v.Pos, p.hdr.Samples[s], freq[1], tmp[s][0], tmp[s][1], tmp[s][2], phap, pdip)
		}
	}
	return nil, nil
}

// accumulatePL fills tmp and freq from FORMAT/PL, mirroring the PL branch of
// process_region_guess(). It returns false when the record has no usable PL.
func (p *guessPloidyPlugin) accumulatePL(v *vcf.Variant, nAllele, nsmp int, tmp [][3]float64, freq *[2]float64, smplPass []bool) bool {
	pls, npl := parseFormatIntAll(v, "PL")
	if pls == nil || npl <= 0 {
		return false
	}
	ndipGT := nAllele * (nAllele + 1) / 2
	switch npl {
	case ndipGT: // diploid
		for s := 0; s < nsmp; s++ {
			if smplPass != nil && !smplPass[s] {
				tmp[s][0] = -1
				continue
			}
			ptr := pls[s]
			if len(ptr) < 3 || ptr[0] == intMissing || ptr[1] == intMissing || ptr[2] == intMissing {
				tmp[s][0] = -1
				continue
			}
			if ptr[0] == ptr[1] && ptr[0] == ptr[2] {
				tmp[s][0] = -1
				continue
			}
			for i := 0; i < 3; i++ {
				tmp[s][i] = guessPLToProb(p.pl2p, ptr[i])
			}
			sum := tmp[s][0] + tmp[s][1] + tmp[s][2]
			for i := 0; i < 3; i++ {
				tmp[s][i] /= sum
			}
			freq[0] += 2*tmp[s][0] + tmp[s][1]
			freq[1] += tmp[s][1] + 2*tmp[s][2]
		}
	case nAllele: // all haploid
		for s := 0; s < nsmp; s++ {
			if smplPass != nil && !smplPass[s] {
				tmp[s][0] = -1
				continue
			}
			ptr := pls[s]
			if len(ptr) < 2 || ptr[0] == intMissing || ptr[1] == intMissing {
				tmp[s][0] = -1
				continue
			}
			tmp[s][0] = guessPLToProb(p.pl2p, ptr[0])
			tmp[s][1] = p.pl2p[255]
			tmp[s][2] = guessPLToProb(p.pl2p, ptr[1])
			sum := tmp[s][0] + tmp[s][1] + tmp[s][2]
			for i := 0; i < 3; i++ {
				tmp[s][i] /= sum
			}
			freq[0] += tmp[s][0]
			freq[1] += tmp[s][2]
		}
	default:
		return false
	}
	return true
}

// accumulateGL fills tmp and freq from FORMAT/GL, mirroring the GL branch of
// process_region_guess(). It returns false when the record has no usable GL.
func (p *guessPloidyPlugin) accumulateGL(v *vcf.Variant, nAllele, nsmp int, tmp [][3]float64, freq *[2]float64, smplPass []bool) bool {
	gls, ngl := parseFormatFloatAll(v, "GL")
	if gls == nil || ngl <= 0 {
		return false
	}
	ndipGT := nAllele * (nAllele + 1) / 2
	switch ngl {
	case ndipGT: // diploid
		for s := 0; s < nsmp; s++ {
			if smplPass != nil && !smplPass[s] {
				tmp[s][0] = -1
				continue
			}
			ptr := gls[s]
			if len(ptr) < 3 || isFloatMissing(ptr[0]) || isFloatMissing(ptr[1]) || isFloatMissing(ptr[2]) {
				tmp[s][0] = -1
				continue
			}
			if ptr[0] == ptr[1] && ptr[0] == ptr[2] {
				tmp[s][0] = -1
				continue
			}
			for i := 0; i < 3; i++ {
				tmp[s][i] = math.Pow(10, ptr[i])
			}
			sum := tmp[s][0] + tmp[s][1] + tmp[s][2]
			for i := 0; i < 3; i++ {
				tmp[s][i] /= sum
			}
			freq[0] += 2*tmp[s][0] + tmp[s][1]
			freq[1] += tmp[s][1] + 2*tmp[s][2]
		}
	case nAllele: // all haploid
		for s := 0; s < nsmp; s++ {
			if smplPass != nil && !smplPass[s] {
				tmp[s][0] = -1
				continue
			}
			ptr := gls[s]
			if len(ptr) < 2 || isFloatMissing(ptr[0]) || isFloatMissing(ptr[1]) {
				tmp[s][0] = -1
				continue
			}
			tmp[s][0] = math.Pow(10, ptr[0])
			tmp[s][1] = 1e-26
			tmp[s][2] = math.Pow(10, ptr[1])
			sum := tmp[s][0] + tmp[s][1] + tmp[s][2]
			for i := 0; i < 3; i++ {
				tmp[s][i] /= sum
			}
			freq[0] += tmp[s][0]
			freq[1] += tmp[s][2]
		}
	default:
		return false
	}
	return true
}

// Destroy prints the per-sample sex predictions, mirroring guess-ploidy.c run().
func (p *guessPloidyPlugin) Destroy() error {
	if p.out == nil {
		return nil
	}
	for i := range p.counts {
		var phap, pdip float64 = 0.5, 0.5
		if p.counts[i].ncount != 0 {
			phap = p.counts[i].phap / float64(p.counts[i].ncount)
			pdip = p.counts[i].pdip / float64(p.counts[i].ncount)
		}
		sex := 'U'
		if phap > pdip {
			sex = 'M'
		} else if phap < pdip {
			sex = 'F'
		}
		if p.verbose > 0 {
			fmt.Fprintf(p.out, "SEX\t%s\t%c\t%f\t%f\t%d\t%f\n",
				p.hdr.Samples[i], sex, phap, pdip, p.counts[i].ncount, phap-pdip)
		} else {
			fmt.Fprintf(p.out, "%s\t%c\n", p.hdr.Samples[i], sex)
		}
	}
	return nil
}

// plToProb maps a PL value to a probability via the pl2p table, clamping out-of
// range values to pl2p[255], matching the (ptr<0||ptr>=256)?pl2p[255] guard.
func guessPLToProb(pl2p []float64, x int) float64 {
	if x < 0 || x >= 256 {
		return pl2p[255]
	}
	return pl2p[x]
}

// intMissing is the sentinel for a missing FORMAT integer value (".").
const intMissing = math.MinInt32

// parseFormatIntAll returns per-sample integer vectors for a FORMAT tag and the
// per-sample value count (the first sample's length), mirroring how
// bcf_get_format_int32 reports n/nsmpl. Missing entries become intMissing.
func parseFormatIntAll(v *vcf.Variant, tag string) ([][]int, int) {
	if !formatHasTag(v, tag) {
		return nil, 0
	}
	out := make([][]int, len(v.Samples))
	per := 0
	for i := range v.Samples {
		s := v.Samples[i].Data[tag]
		var row []int
		if s != "" && s != "." {
			for _, tok := range strings.Split(s, ",") {
				if tok == "." || tok == "" {
					row = append(row, intMissing)
					continue
				}
				n, err := strconv.Atoi(tok)
				if err != nil {
					row = append(row, intMissing)
					continue
				}
				row = append(row, n)
			}
		}
		out[i] = row
		if i == 0 {
			per = len(row)
		}
	}
	return out, per
}

// parseFormatFloatAll is the float counterpart of parseFormatIntAll for GL.
func parseFormatFloatAll(v *vcf.Variant, tag string) ([][]float64, int) {
	if !formatHasTag(v, tag) {
		return nil, 0
	}
	out := make([][]float64, len(v.Samples))
	per := 0
	for i := range v.Samples {
		s := v.Samples[i].Data[tag]
		var row []float64
		if s != "" && s != "." {
			for _, tok := range strings.Split(s, ",") {
				if tok == "." || tok == "" {
					row = append(row, math.NaN())
					continue
				}
				f, err := strconv.ParseFloat(tok, 64)
				if err != nil {
					row = append(row, math.NaN())
					continue
				}
				row = append(row, f)
			}
		}
		out[i] = row
		if i == 0 {
			per = len(row)
		}
	}
	return out, per
}

// isFloatMissing reports whether a parsed GL value is missing (NaN sentinel).
func isFloatMissing(f float64) bool { return math.IsNaN(f) }

// infoFloat0 returns the first float of an INFO tag, mirroring the
// bcf_get_info_float(...)>0 use of args->af[0].
func infoFloat0(v *vcf.Variant, tag string) (float64, bool) {
	s, ok := v.Info[tag]
	if !ok || s == "" || s == "." {
		return 0, false
	}
	first := s
	if idx := strings.IndexByte(s, ','); idx >= 0 {
		first = s[:idx]
	}
	f, err := strconv.ParseFloat(first, 64)
	if err != nil {
		return 0, false
	}
	return f, true
}

// parseAttachedVerbosity recognises getopt's attached optional-argument forms
// for the verbosity flag: -v<N>, --verbosity=<N> and --verbose=<N>. It returns
// the level and true when the token is one of those forms; a non-numeric value
// yields level -1 (a parse error for the caller to report).
func parseAttachedVerbosity(tok string) (int, bool) {
	var val string
	switch {
	case strings.HasPrefix(tok, "-v") && len(tok) > 2 && tok[2] != '-':
		val = tok[2:]
	case strings.HasPrefix(tok, "--verbosity="):
		val = tok[len("--verbosity="):]
	case strings.HasPrefix(tok, "--verbose="):
		val = tok[len("--verbose="):]
	default:
		return 0, false
	}
	n, err := strconv.Atoi(val)
	if err != nil || n < 0 {
		return -1, true
	}
	return n, true
}

// hasInfoHeaderMeta reports whether the header declares an INFO tag with the
// given ID.
func hasInfoHeaderMeta(meta []string, id string) bool {
	for _, m := range meta {
		if headerKind(m) == "##INFO" && headerID(m) == id {
			return true
		}
	}
	return false
}
