// Native port of the upstream `fill-tags` plugin (plugins/fill-tags.c) and its
// deprecated subset `fill-AN-AC` (plugins/fill-AN-AC.c).
//
// fill-tags recomputes INFO/FORMAT annotations from FORMAT/GT: AN, AC, AC_Hom,
// AC_Het, AC_Hemi, AF, MAF, NS, HWE, ExcHet, END, TYPE, and FORMAT/VAF, VAF1.
// The per-allele het/hom/hemi/half counting follows process_fmt's BRANCH_INT
// loop exactly, including the --drop-missing (-d) treatment of half-missing
// genotypes.
//
// This port supports the full upstream surface:
//   - the default tag set and -t LIST selection (incl. "all", the "-" drop
//     prefix, and INFO/FORMAT-qualified names);
//   - -S/--samples-file population grouping (per-population AN/AC/AF/... tags
//     suffixed with _GROUP plus the summary "ALL" population);
//   - the experimental custom expression TAG[:Number]=[int|float](EXPR), via
//     the in-tree fill-tags expression evaluator (native_plugin_filltags_expr.go);
//   - -l/--list-tags, which prints the available-tag table to stderr and exits.
package bcftools

import (
	"errors"
	"fmt"
	"io"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() {
	registerNativePlugin("fill-tags", func() NativePlugin { return &fillTagsPlugin{} })
	registerNativePlugin("fill-AN-AC", func() NativePlugin { return &fillTagsPlugin{anacOnly: true} })
}

// fill-tags tag selection bits, matching the SET_* flags in fill-tags.c.
const (
	setAN = 1 << iota
	setAC
	setACHom
	setACHet
	setACHemi
	setAF
	setNS
	setMAF
	setHWE
	setExcHet
	setEND
	setType
	setVAF
	setVAF1
)

// errFillTagsListed is the sentinel returned after -l/--list-tags has printed
// the available-tag table. It causes the host to exit non-zero with no record
// output, matching upstream's error()-based list_tags() behaviour.
var errFillTagsListed = errors.New("fill-tags: listed tags")

// IsListTagsError reports whether err is the fill-tags -l/--list-tags sentinel.
// The CLI uses it to exit non-zero (matching upstream's error() exit) without
// printing the generic "plugin failed" wrapper, since the available-tag table
// has already been written to stderr by Init.
func IsListTagsError(err error) bool {
	return errors.Is(err, errFillTagsListed)
}

// filltagsListText is the exact text upstream's list_tags() emits (printed to
// stderr via error()).
const filltagsListText = `INFO/AC        Number:A  Type:Integer  ..  Allele count in genotypes
INFO/AC_Hom    Number:A  Type:Integer  ..  Allele counts in homozygous genotypes
INFO/AC_Het    Number:A  Type:Integer  ..  Allele counts in heterozygous genotypes
INFO/AC_Hemi   Number:A  Type:Integer  ..  Allele counts in hemizygous genotypes
INFO/AF        Number:A  Type:Float    ..  Allele frequency from FMT/GT or AC,AN if FMT/GT is not present
INFO/AN        Number:1  Type:Integer  ..  Total number of alleles in called genotypes
INFO/ExcHet    Number:A  Type:Float    ..  Test excess heterozygosity; 1=good, 0=bad
INFO/END       Number:1  Type:Integer  ..  End position of the variant
INFO/F_MISSING Number:1  Type:Float    ..  Fraction of missing genotypes, synonymous with 'F_MISSING=F_PASS(GT="mis")'
INFO/HWE       Number:A  Type:Float    ..  HWE test (PMID:15789306); 1=good, 0=bad
INFO/MAF       Number:1  Type:Float    ..  Frequency of the second most common allele
INFO/NS        Number:1  Type:Integer  ..  Number of samples with data
INFO/TYPE      Number:.  Type:String   ..  The record type (REF,SNP,MNP,INDEL,etc)
FORMAT/VAF     Number:A  Type:Float    ..  The fraction of reads with the alternate allele, requires FORMAT/AD or ADF+ADR
FORMAT/VAF1    Number:1  Type:Float    ..  The same as FORMAT/VAF but for all alternate alleles cumulatively
TAG:Number=Type(EXPR)                  ..  Experimental support for user expressions such as DP:1=int(sum(DP))
               If Number and Type are not given (e.g. DP=sum(DP)), variable number (Number=.) of floating point
               values (Type=Float) will be used.
`

// ftfBinding is one compiled custom-expression tag for one population.
type ftfBinding struct {
	dstTag   string // tag name (without population suffix)
	suffix   string // population suffix ("" for ALL)
	isFormat bool   // FORMAT (per-sample) vs INFO (site) destination
	isInt    bool   // Integer (int()/integer()) vs Float destination
	fixedLen bool   // Number=N (fixed) vs Number=. (variable)
	count    int    // N when fixedLen
	expr     *fillExpr
	usmpl    []bool // population sample mask (nil for ALL)
}

// population is one sample group (a -S group, or the summary "ALL" group).
type population struct {
	name   string // group name ("" for ALL)
	suffix string // tag suffix ("" for ALL, "_NAME" otherwise)
	mask   []bool // per-sample membership mask (nil for ALL = every sample)
}

// alleleCounts accumulates, per allele, the het/hom/hemi/half-missing counts
// used to derive every fill-tags annotation. It mirrors fill-tags.c counts_t.
type alleleCounts struct {
	nhom, nhet, nhemi, nac int
}

// fillTagsPlugin implements fill-tags / fill-AN-AC.
type fillTagsPlugin struct {
	anacOnly    bool // true for the fill-AN-AC entry point
	tags        int
	dropMissing bool
	listTags    bool
	pops        []population
	bindings    []*ftfBinding // custom-expression tags, in parse order across pops
	stderr      io.Writer
}

// Name returns the plugin name.
func (p *fillTagsPlugin) Name() string {
	if p.anacOnly {
		return "fill-AN-AC"
	}
	return "fill-tags"
}

// About returns the one-line description, matching the upstream about().
func (p *fillTagsPlugin) About() string {
	if p.anacOnly {
		return "Fill INFO fields AN and AC. This plugin is DEPRECATED, use fill-tags instead."
	}
	return "Set INFO tags AF, AC, AC_Hemi, AC_Hom, AC_Het, AN, ExcHet, HWE, MAF, NS; FORMAT/VAF and more."
}

// Parallel reports true: each record is recomputed independently.
func (p *fillTagsPlugin) Parallel() bool { return true }

// SetStderr wires the host stderr so -l can print the available-tag table.
func (p *fillTagsPlugin) SetStderr(w io.Writer) { p.stderr = w }

// Init parses the plugin arguments and appends the relevant ##INFO/##FORMAT
// header lines in upstream's fixed order.
func (p *fillTagsPlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	tagsStr := "all"
	var samplesFile string
	if !p.anacOnly {
		for i := 0; i < len(args); i++ {
			a := args[i]
			switch a {
			case "-d", "--drop-missing":
				p.dropMissing = true
			case "-l", "--list-tags":
				p.listTags = true
			case "-t", "--tags":
				if i+1 >= len(args) {
					return nil, fmt.Errorf("fill-tags: -t requires an argument")
				}
				i++
				tagsStr = args[i]
			case "-S", "--samples-file":
				if i+1 >= len(args) {
					return nil, fmt.Errorf("fill-tags: -S requires an argument")
				}
				i++
				samplesFile = args[i]
			default:
				if strings.HasPrefix(a, "-t") && len(a) > 2 {
					tagsStr = a[2:]
					continue
				}
				if strings.HasPrefix(a, "-S") && len(a) > 2 {
					samplesFile = a[2:]
					continue
				}
				return nil, fmt.Errorf("fill-tags: unsupported option %q", a)
			}
		}
	}

	if p.listTags {
		if p.stderr != nil {
			fmt.Fprint(p.stderr, filltagsListText)
		}
		return nil, errFillTagsListed
	}

	// Build the populations: the -S groups first, then the summary "ALL"
	// population appended last (matching init_pops, which makes ALL the final
	// pop so its empty suffix sorts after the named groups in every tag loop).
	if samplesFile != "" {
		groups, err := parseFillTagsSamples(samplesFile, hdr, p.stderr)
		if err != nil {
			return nil, err
		}
		p.pops = append(p.pops, groups...)
	}
	p.pops = append(p.pops, population{name: "", suffix: "", mask: nil})

	out := &vcf.Header{Samples: hdr.Samples}
	out.MetaInfo = append(out.MetaInfo, hdr.MetaInfo...)

	flag, err := p.parseTags(tagsStr, hdr, out)
	if err != nil {
		return nil, err
	}
	p.tags = flag

	add := func(line string) { out.MetaInfo = appendInfoHeader(out.MetaInfo, line) }
	// Per-population header lines for the FORMAT-derived tags, in upstream's
	// fixed hdr_append order. The custom-expression tag headers were already
	// appended during parseTags (matching upstream's parse_func).
	hdrAppend := func(buildLine func(pop population) string) {
		for _, pop := range p.pops {
			add(buildLine(pop))
		}
	}
	desc := func(pop population, base string) string {
		if pop.name == "" {
			return base
		}
		return base + " in " + pop.name
	}
	if p.tags&setAN != 0 {
		hdrAppend(func(pop population) string {
			return fmt.Sprintf(`##INFO=<ID=AN%s,Number=1,Type=Integer,Description="%s">`, pop.suffix, desc(pop, "Total number of alleles in called genotypes"))
		})
	}
	if p.tags&setAC != 0 {
		hdrAppend(func(pop population) string {
			return fmt.Sprintf(`##INFO=<ID=AC%s,Number=A,Type=Integer,Description="%s">`, pop.suffix, desc(pop, "Allele count in genotypes"))
		})
	}
	if p.tags&setNS != 0 {
		hdrAppend(func(pop population) string {
			return fmt.Sprintf(`##INFO=<ID=NS%s,Number=1,Type=Integer,Description="%s">`, pop.suffix, desc(pop, "Number of samples with data"))
		})
	}
	if p.tags&setACHom != 0 {
		hdrAppend(func(pop population) string {
			return fmt.Sprintf(`##INFO=<ID=AC_Hom%s,Number=A,Type=Integer,Description="%s">`, pop.suffix, desc(pop, "Allele counts in homozygous genotypes"))
		})
	}
	if p.tags&setACHet != 0 {
		hdrAppend(func(pop population) string {
			return fmt.Sprintf(`##INFO=<ID=AC_Het%s,Number=A,Type=Integer,Description="%s">`, pop.suffix, desc(pop, "Allele counts in heterozygous genotypes"))
		})
	}
	if p.tags&setACHemi != 0 {
		hdrAppend(func(pop population) string {
			return fmt.Sprintf(`##INFO=<ID=AC_Hemi%s,Number=A,Type=Integer,Description="%s">`, pop.suffix, desc(pop, "Allele counts in hemizygous genotypes"))
		})
	}
	if p.tags&setAF != 0 {
		hdrAppend(func(pop population) string {
			return fmt.Sprintf(`##INFO=<ID=AF%s,Number=A,Type=Float,Description="%s">`, pop.suffix, desc(pop, "Allele frequency"))
		})
	}
	if p.tags&setMAF != 0 {
		hdrAppend(func(pop population) string {
			return fmt.Sprintf(`##INFO=<ID=MAF%s,Number=1,Type=Float,Description="%s">`, pop.suffix, desc(pop, "Frequency of the second most common allele"))
		})
	}
	if p.tags&setHWE != 0 {
		hdrAppend(func(pop population) string {
			base := "HWE test"
			if pop.name != "" {
				base += " in " + pop.name
			}
			return fmt.Sprintf(`##INFO=<ID=HWE%s,Number=A,Type=Float,Description="%s (PMID:15789306); 1=good, 0=bad">`, pop.suffix, base)
		})
	}
	if p.tags&setEND != 0 {
		add(`##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant">`)
	}
	if p.tags&setType != 0 {
		add(`##INFO=<ID=TYPE,Number=.,Type=String,Description="Variant type">`)
	}
	if p.tags&setExcHet != 0 {
		hdrAppend(func(pop population) string {
			base := "Test excess heterozygosity"
			if pop.name != "" {
				base += " in " + pop.name
			}
			return fmt.Sprintf(`##INFO=<ID=ExcHet%s,Number=A,Type=Float,Description="%s; 1=good, 0=bad">`, pop.suffix, base)
		})
	}
	if p.tags&setVAF != 0 {
		add(`##FORMAT=<ID=VAF,Number=A,Type=Float,Description="The fraction of reads with alternate allele (nALT/nSumAll)">`)
	}
	if p.tags&setVAF1 != 0 {
		add(`##FORMAT=<ID=VAF1,Number=1,Type=Float,Description="The fraction of reads with alternate alleles (nSumALT/nSumAll)">`)
	}
	return out, nil
}

// parseTags converts the comma-separated -t list into the SET_* bitmask and
// compiles any custom-expression tags (TAG=EXPR), appending their ##INFO /
// ##FORMAT header lines (per population) to out. Unknown tokens are an error,
// matching upstream parse_tags (which has no "-" drop syntax of its own).
func (p *fillTagsPlugin) parseTags(str string, hdr, out *vcf.Header) (int, error) {
	if p.anacOnly {
		// fill-AN-AC ignores -t and always fills AN,AC for the ALL population.
		return setAN | setAC, nil
	}
	flag := 0
	for _, raw := range strings.Split(str, ",") {
		t := strings.TrimSpace(raw)
		bit, isExpr, err := p.classifyTag(t, hdr, out)
		if err != nil {
			return 0, err
		}
		if isExpr {
			continue
		}
		flag |= bit
	}
	return flag, nil
}

// classifyTag maps a single tag token to its SET_* bit, or compiles it as a
// custom expression when it contains '='. isExpr is true for the expression
// case (no bit contributed).
func (p *fillTagsPlugin) classifyTag(t string, hdr, out *vcf.Header) (bit int, isExpr bool, err error) {
	norm := t
	norm = strings.TrimPrefix(norm, "INFO/")
	norm = strings.TrimPrefix(norm, "FORMAT/")
	switch strings.ToLower(norm) {
	case "all":
		return setAN | setAC | setACHom | setACHet | setACHemi | setAF | setNS | setMAF | setHWE | setExcHet | setVAF | setVAF1, false, p.addExpr("F_MISSING:1=F_MISSING", "F_MISSING", hdr, out)
	case "an":
		return setAN, false, nil
	case "ac":
		return setAC, false, nil
	case "ns":
		return setNS, false, nil
	case "ac_hom":
		return setACHom, false, nil
	case "ac_het":
		return setACHet, false, nil
	case "ac_hemi":
		return setACHemi, false, nil
	case "af":
		return setAF, false, nil
	case "maf":
		return setMAF, false, nil
	case "hwe":
		return setHWE, false, nil
	case "exchet":
		return setExcHet, false, nil
	case "end":
		return setEND, false, nil
	case "type":
		return setType, false, nil
	case "vaf":
		return setVAF, false, nil
	case "vaf1":
		return setVAF1, false, nil
	case "f_missing":
		return 0, true, p.addExpr("F_MISSING:1=F_MISSING", "F_MISSING", hdr, out)
	}
	if idx := strings.IndexByte(t, '='); idx >= 0 {
		return 0, true, p.addExpr(t, t[idx+1:], hdr, out)
	}
	return 0, false, fmt.Errorf("Error parsing \"--tags %s\": the tag \"%s\" is not supported", t, t)
}

// addExpr compiles a TAG[:Number]=[int|float](EXPR) custom expression for every
// population and appends its header line(s).
func (p *fillTagsPlugin) addExpr(tagExpr, expr string, hdr, out *vcf.Header) error {
	eq := strings.IndexByte(tagExpr, '=')
	if eq < 0 {
		return fmt.Errorf("fill-tags: could not parse the expression: %s", tagExpr)
	}
	dst := tagExpr[:eq]
	isFormat := false
	switch {
	case strings.HasPrefix(strings.ToLower(dst), "info/"):
		dst = dst[5:]
	case strings.HasPrefix(strings.ToLower(dst), "format/"):
		dst = dst[7:]
		isFormat = true
	case strings.HasPrefix(strings.ToLower(dst), "fmt/"):
		dst = dst[4:]
		isFormat = true
	}
	fixedLen := false
	count := 0
	if c := strings.IndexByte(dst, ':'); c >= 0 {
		numStr := dst[c+1:]
		dst = dst[:c]
		n, perr := strconv.Atoi(numStr)
		if perr != nil {
			return fmt.Errorf("fill-tags: could not parse the expression: %s", tagExpr)
		}
		count = n
		fixedLen = true
	}

	isInt := false
	inner := expr
	lower := strings.ToLower(expr)
	if strings.HasSuffix(expr, ")") {
		switch {
		case strings.HasPrefix(lower, "int("):
			inner = expr[4 : len(expr)-1]
			isInt = true
		case strings.HasPrefix(lower, "integer("):
			inner = expr[8 : len(expr)-1]
			isInt = true
		case strings.HasPrefix(lower, "float("):
			inner = expr[6 : len(expr)-1]
			isInt = false
		}
	}

	compiled, err := compileFillExpr(inner, hdr)
	if err != nil {
		return err
	}

	typeStr := "Float"
	if isInt {
		typeStr = "Integer"
	}
	numField := "."
	if fixedLen {
		numField = strconv.Itoa(count)
	}
	kind := "INFO"
	if isFormat {
		kind = "FORMAT"
	}

	for _, pop := range p.pops {
		b := &ftfBinding{
			dstTag:   dst,
			suffix:   pop.suffix,
			isFormat: isFormat,
			isInt:    isInt,
			fixedLen: fixedLen,
			count:    count,
			expr:     compiled,
			usmpl:    pop.mask,
		}
		p.bindings = append(p.bindings, b)

		descTail := ""
		if pop.name != "" {
			descTail = " in " + pop.name
		}
		line := fmt.Sprintf(`##%s=<ID=%s%s,Number=%s,Type=%s,Description="Added by +fill-tags expression %s%s">`,
			kind, dst, pop.suffix, numField, typeStr, escapeHeaderQuotes(tagExpr), descTail)
		out.MetaInfo = appendInfoHeader(out.MetaInfo, line)
	}
	return nil
}

// Process recomputes the requested annotations for a single record.
func (p *fillTagsPlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	nals := 1 + len(v.Alt)
	hasGT := formatHasTag(v, "GT")

	// Custom-expression tags run first (matching upstream, where the ftf
	// functions execute before process_fmt).
	for _, b := range p.bindings {
		p.fillExprTag(v, b)
	}

	if hasGT && p.tags&(setAN|setAC|setACHom|setACHet|setACHemi|setAF|setMAF|setNS|setHWE|setExcHet) != 0 {
		// Precompute per-population counts once.
		popCounts := make([][]alleleCounts, len(p.pops))
		popNS := make([]int, len(p.pops))
		for i, pop := range p.pops {
			popCounts[i], popNS[i] = p.popCounts(v, pop, nals)
		}
		if p.tags&setNS != 0 {
			for i, pop := range p.pops {
				setInfo(v, "NS"+pop.suffix, strconv.Itoa(popNS[i]))
			}
		}
		if p.tags&setAN != 0 {
			for i, pop := range p.pops {
				an := 0
				for j := 0; j < nals; j++ {
					an += alleleTotal(popCounts[i][j])
				}
				setInfo(v, "AN"+pop.suffix, strconv.Itoa(an))
			}
		}
		if p.tags&(setAF|setMAF) != 0 {
			for i, pop := range p.pops {
				an := 0
				for j := 0; j < nals; j++ {
					an += alleleTotal(popCounts[i][j])
				}
				p.fillAFMAF(v, pop, popCounts[i], nals, an)
			}
		}
		if p.tags&setAC != 0 {
			for i, pop := range p.pops {
				p.fillAC(v, pop, popCounts[i], nals)
			}
		}
		if p.tags&setACHet != 0 {
			for i, pop := range p.pops {
				p.fillPerAlt(v, "AC_Het"+pop.suffix, popCounts[i], nals, func(c alleleCounts) int { return c.nhet })
			}
		}
		if p.tags&setACHom != 0 {
			for i, pop := range p.pops {
				p.fillPerAlt(v, "AC_Hom"+pop.suffix, popCounts[i], nals, func(c alleleCounts) int { return c.nhom })
			}
		}
		if p.tags&setACHemi != 0 && nals > 1 {
			for i, pop := range p.pops {
				p.fillPerAlt(v, "AC_Hemi"+pop.suffix, popCounts[i], nals, func(c alleleCounts) int { return c.nhemi })
			}
		}
		if p.tags&(setHWE|setExcHet) != 0 {
			for i, pop := range p.pops {
				p.fillHWE(v, pop, popCounts[i], nals)
			}
		}
	}

	// Sites-only AF from AN/AC (process_info_af): only when there are no samples.
	if p.tags&setAF != 0 && len(v.Samples) == 0 {
		p.fillSitesOnlyAF(v, nals)
	}

	if p.tags&(setVAF|setVAF1) != 0 {
		p.fillVAF(v, nals)
	}
	if p.tags&setEND != 0 {
		setInfo(v, "END", strconv.Itoa(v.Pos+refLen(v)-1))
	}
	if p.tags&setType != 0 {
		setInfo(v, "TYPE", typeMacroValue(v))
	}
	return []*vcf.Variant{v}, nil
}

// Destroy releases resources (none held).
func (p *fillTagsPlugin) Destroy() error { return nil }

// popCounts returns the per-allele counts and the NS (samples-with-data) count
// for one population. Mirrors process_fmt's BRANCH_INT counting restricted to
// the population's sample set.
func (p *fillTagsPlugin) popCounts(v *vcf.Variant, pop population, nals int) ([]alleleCounts, int) {
	counts := make([]alleleCounts, nals)
	ns := 0
	for i := range v.Samples {
		if pop.mask != nil && (i >= len(pop.mask) || !pop.mask[i]) {
			continue
		}
		gt, ok := sampleGT(v, i)
		if !ok {
			continue
		}
		present := make([]bool, nals)
		nbits := 0
		ncalled := 0
		for _, a := range gt.alleles {
			if a == missingAllele {
				continue
			}
			if a < 0 || a >= nals {
				continue
			}
			ncalled++
			if !present[a] {
				nbits++
				present[a] = true
			}
		}
		if ncalled == 0 {
			continue
		}
		isHom := nbits == 1
		var isHemi, isHalf bool
		switch {
		case ncalled != gt.ploidy():
			if p.dropMissing {
				isHemi, isHalf = false, true
			} else {
				isHemi, isHalf = true, false
			}
		case ncalled == 1:
			isHemi, isHalf = true, false
		default:
			isHemi, isHalf = false, false
		}
		for a := 0; a < nals; a++ {
			if !present[a] {
				continue
			}
			switch {
			case isHalf:
				counts[a].nac++
			case !isHom:
				counts[a].nhet++
			case !isHemi:
				counts[a].nhom += 2
			default:
				counts[a].nhemi++
			}
		}
		ns++
	}
	return counts, ns
}

// alleleTotal is the total called count of allele j (het+hom+hemi+half).
func alleleTotal(c alleleCounts) int { return c.nhet + c.nhom + c.nhemi + c.nac }

func (p *fillTagsPlugin) fillAC(v *vcf.Variant, pop population, counts []alleleCounts, nals int) {
	parts := make([]string, 0, nals-1)
	for j := 1; j < nals; j++ {
		parts = append(parts, strconv.Itoa(alleleTotal(counts[j])))
	}
	setInfo(v, "AC"+pop.suffix, strings.Join(parts, ","))
}

func (p *fillTagsPlugin) fillPerAlt(v *vcf.Variant, key string, counts []alleleCounts, nals int, sel func(alleleCounts) int) {
	parts := make([]string, 0, nals-1)
	for j := 1; j < nals; j++ {
		parts = append(parts, strconv.Itoa(sel(counts[j])))
	}
	setInfo(v, key, strings.Join(parts, ","))
}

// fillAFMAF computes per-allele frequencies for one population.
func (p *fillTagsPlugin) fillAFMAF(v *vcf.Variant, pop population, counts []alleleCounts, nals int, an int) {
	freq := make([]float64, nals)
	if nals > 1 {
		for j := 0; j < nals; j++ {
			freq[j] = float64(alleleTotal(counts[j]))
		}
		if an != 0 {
			for j := 0; j < nals; j++ {
				freq[j] /= float64(an)
			}
		}
	}
	if p.tags&setAF != 0 {
		parts := make([]string, 0, nals-1)
		for j := 1; j < nals; j++ {
			if nals > 1 && an == 0 {
				parts = append(parts, ".")
			} else {
				parts = append(parts, formatVCFFloat(freq[j]))
			}
		}
		setInfo(v, "AF"+pop.suffix, strings.Join(parts, ","))
	}
	if nals > 1 && p.tags&setMAF != 0 {
		sorted := append([]float64(nil), freq...)
		if an != 0 {
			sort.Sort(sort.Reverse(sort.Float64Slice(sorted)))
		}
		maf := sorted[1]
		if an == 0 {
			setInfo(v, "MAF"+pop.suffix, ".")
		} else {
			setInfo(v, "MAF"+pop.suffix, formatVCFFloat(maf))
		}
	}
}

// fillHWE computes the HWE and ExcHet p-values per ALT allele for one
// population, porting calc_hwe (Wigginton 2005, PMID 15789306).
func (p *fillTagsPlugin) fillHWE(v *vcf.Variant, pop population, counts []alleleCounts, nals int) {
	hwe := make([]string, 0, nals-1)
	exc := make([]string, 0, nals-1)
	if nals > 1 {
		nrefTot := counts[0].nhom
		for j := 0; j < nals; j++ {
			nrefTot += counts[j].nhet
		}
		for j := 1; j < nals; j++ {
			nref := nrefTot - counts[j].nhet
			nalt := counts[j].nhet + counts[j].nhom
			nhet := counts[j].nhet
			var ph, pe float64 = 1, 1
			if nref > 0 && nalt > 0 {
				ph, pe = calcHWE(nref, nalt, nhet)
			}
			hwe = append(hwe, formatVCFFloat(ph))
			exc = append(exc, formatVCFFloat(pe))
		}
	}
	if p.tags&setHWE != 0 {
		setInfo(v, "HWE"+pop.suffix, strings.Join(hwe, ","))
	}
	if p.tags&setExcHet != 0 {
		setInfo(v, "ExcHet"+pop.suffix, strings.Join(exc, ","))
	}
}

// fillSitesOnlyAF computes INFO/AF from INFO/AN and INFO/AC when the record has
// no samples, mirroring process_info_af.
func (p *fillTagsPlugin) fillSitesOnlyAF(v *vcf.Variant, nals int) {
	anStr, ok := v.Info["AN"]
	if !ok {
		return
	}
	an, err := strconv.Atoi(strings.TrimSpace(anStr))
	if err != nil || an == 0 {
		return
	}
	acStr, ok := v.Info["AC"]
	if !ok {
		return
	}
	acParts := strings.Split(acStr, ",")
	if len(acParts) != nals-1 {
		return
	}
	parts := make([]string, 0, len(acParts))
	for _, s := range acParts {
		n, perr := strconv.Atoi(strings.TrimSpace(s))
		if perr != nil {
			return
		}
		parts = append(parts, formatVCFFloat(float64(n)/float64(an)))
	}
	setInfo(v, "AF", strings.Join(parts, ","))
}

// fillExprTag evaluates a custom-expression binding and writes its INFO or
// FORMAT value(s), porting ftf_filter_expr.
func (p *fillTagsPlugin) fillExprTag(v *vcf.Variant, b *ftfBinding) {
	res := b.expr.evaluate(v, b.usmpl)
	key := b.dstTag + b.suffix
	if !b.isFormat {
		var vals []float64
		if res.perSample {
			for i := range v.Samples {
				if b.usmpl != nil && (i >= len(b.usmpl) || !b.usmpl[i]) {
					continue
				}
				vals = res.values[i]
				break
			}
		} else {
			vals = res.site
		}
		if len(vals) == 0 && !b.fixedLen {
			return // nothing to set; INFO field stays absent
		}
		nfill := len(vals)
		if b.fixedLen {
			nfill = b.count
		}
		if nfill == 0 {
			return
		}
		parts := make([]string, nfill)
		for j := 0; j < nfill; j++ {
			if j < len(vals) {
				parts[j] = formatExprValue(vals[j], b.isInt)
			} else {
				parts[j] = "."
			}
		}
		setInfo(v, key, strings.Join(parts, ","))
		return
	}

	// FORMAT destination: one value (or fixed-count vector) per sample.
	nval1 := 1
	if b.fixedLen {
		nval1 = b.count
	} else {
		for i := range v.Samples {
			if res.perSample && i < len(res.values) && len(res.values[i]) > nval1 {
				nval1 = len(res.values[i])
			}
		}
	}
	ensureFormatTag(v, key)
	for i := range v.Samples {
		var vec []float64
		if res.perSample && i < len(res.values) {
			vec = res.values[i]
		} else {
			vec = res.site
		}
		parts := make([]string, nval1)
		for j := 0; j < nval1; j++ {
			if j < len(vec) {
				parts[j] = formatExprValue(vec[j], b.isInt)
			} else {
				parts[j] = "."
			}
		}
		v.Samples[i].Data[key] = strings.Join(parts, ",")
	}
}

// formatExprValue renders an expression value as Integer or Float, mapping the
// missing sentinel to ".". Integer rounding uses round-half-away-from-zero to
// match C's round() used by int32_from_double.
func formatExprValue(x float64, isInt bool) string {
	if isExprMissing(x) {
		return "."
	}
	if isInt {
		return strconv.FormatInt(int64(roundHalfAway(x)), 10)
	}
	return formatVCFFloat(x)
}

// roundHalfAway rounds to the nearest integer, ties away from zero, matching
// C's round().
func roundHalfAway(x float64) float64 {
	if x < 0 {
		return -floorAdd(-x)
	}
	return floorAdd(x)
}

func floorAdd(x float64) float64 {
	return float64(int64(x + 0.5))
}

// fillVAF computes FORMAT/VAF and VAF1 from FORMAT/AD, porting process_vaf.
func (p *fillTagsPlugin) fillVAF(v *vcf.Variant, nals int) {
	if nals <= 1 {
		return
	}
	if !formatHasTag(v, "AD") {
		return
	}
	doVAF := p.tags&setVAF != 0
	doVAF1 := p.tags&setVAF1 != 0
	var vafVals, vaf1Vals []string
	for i := range v.Samples {
		ad, ok := v.Samples[i].Data["AD"]
		adVals, parseOK := parseADList(ad)
		valid := ok && parseOK && len(adVals) == nals
		sum := 0
		if valid {
			for _, x := range adVals {
				sum += x
			}
		}
		if doVAF {
			if !valid {
				vafVals = append(vafVals, ".")
			} else {
				per := make([]string, nals-1)
				for j := 1; j < nals; j++ {
					if sum != 0 {
						per[j-1] = formatVCFFloat(float64(adVals[j]) / float64(sum))
					} else {
						per[j-1] = "0"
					}
				}
				vafVals = append(vafVals, strings.Join(per, ","))
			}
		}
		if doVAF1 {
			if !valid {
				vaf1Vals = append(vaf1Vals, ".")
			} else if sum != 0 {
				vaf1Vals = append(vaf1Vals, formatVCFFloat(float64(sum-adVals[0])/float64(sum)))
			} else {
				vaf1Vals = append(vaf1Vals, "0")
			}
		}
	}
	if doVAF {
		ensureFormatTag(v, "VAF")
		for i := range v.Samples {
			v.Samples[i].Data["VAF"] = vafVals[i]
		}
	}
	if doVAF1 {
		ensureFormatTag(v, "VAF1")
		for i := range v.Samples {
			v.Samples[i].Data["VAF1"] = vaf1Vals[i]
		}
	}
}

// calcHWE ports fill-tags.c calc_hwe: it returns (p_hwe, p_exc_het) for a
// biallelic configuration with nref reference alleles, nalt alt alleles, and
// nhet heterozygous genotypes (assuming diploid: total genotypes (nref+nalt)/2).
func calcHWE(nref, nalt, nhet int) (float64, float64) {
	ngt := (nref + nalt) / 2
	nrare := nref
	if nalt < nref {
		nrare = nalt
	}
	probs := make([]float64, nrare+1)

	mid := int(float64(nrare) * float64(nref+nalt-nrare) / float64(nref+nalt))
	if (nrare&1)^(mid&1) != 0 {
		mid++
	}
	het := mid
	homR := (nrare - mid) / 2
	homC := ngt - het - homR
	sum := 1.0
	probs[mid] = 1.0

	for het = mid; het > 1; het -= 2 {
		probs[het-2] = probs[het] * float64(het) * float64(het-1) / (4.0 * float64(homR+1) * float64(homC+1))
		sum += probs[het-2]
		homR++
		homC++
	}

	het = mid
	homR = (nrare - mid) / 2
	homC = ngt - het - homR
	for het = mid; het <= nrare-2; het += 2 {
		probs[het+2] = probs[het] * 4.0 * float64(homR) * float64(homC) / (float64(het+2) * float64(het+1))
		sum += probs[het+2]
		homR--
		homC--
	}

	for h := 0; h < nrare+1; h++ {
		probs[h] /= sum
	}

	pExcHet := probs[nhet]
	for h := nhet + 1; h <= nrare; h++ {
		pExcHet += probs[h]
	}

	pHWE := 0.0
	for h := 0; h <= nrare; h++ {
		if probs[h] > probs[nhet] {
			continue
		}
		pHWE += probs[h]
	}
	if pHWE > 1 {
		pHWE = 1
	}
	return pHWE, pExcHet
}

// formatHasTag reports whether tag appears in v.Format.
func formatHasTag(v *vcf.Variant, tag string) bool {
	for _, f := range v.Format {
		if f == tag {
			return true
		}
	}
	return false
}

// parseADList parses a comma-separated FORMAT/AD integer list. A "." entry (or
// an empty/missing field) marks the value as not fully present (false), to
// match upstream's "all values present or skip the sample" handling for AD.
func parseADList(s string) ([]int, bool) {
	if s == "" || s == "." {
		return nil, false
	}
	parts := strings.Split(s, ",")
	out := make([]int, len(parts))
	for i, pp := range parts {
		if pp == "." {
			return nil, false
		}
		n, err := strconv.Atoi(pp)
		if err != nil {
			return nil, false
		}
		out[i] = n
	}
	return out, true
}

// refLen returns the reference-span length of the record (len(REF) for simple
// records). This is used to compute INFO/END = POS + rlen - 1.
func refLen(v *vcf.Variant) int {
	if end, ok := v.Info["END"]; ok {
		if n, err := strconv.Atoi(end); err == nil {
			return n - v.Pos + 1
		}
	}
	return len(v.Ref)
}

// appendInfoHeader inserts an ##INFO/##FORMAT line into the meta lines if a
// definition for the same ID is not already present, preserving upstream's
// behaviour of not duplicating an existing header definition.
func appendInfoHeader(meta []string, line string) []string {
	id := headerID(line)
	if id != "" {
		for _, m := range meta {
			if headerID(m) == id && headerKind(m) == headerKind(line) {
				return meta
			}
		}
	}
	return append(meta, line)
}

// headerID extracts the ID=... value from a structured header line.
func headerID(line string) string {
	i := strings.Index(line, "ID=")
	if i < 0 {
		return ""
	}
	rest := line[i+3:]
	end := strings.IndexAny(rest, ",>")
	if end < 0 {
		return rest
	}
	return rest[:end]
}

// headerKind returns the header category prefix (e.g. "##INFO", "##FORMAT").
func headerKind(line string) string {
	if i := strings.Index(line, "=<"); i >= 0 {
		return line[:i]
	}
	return line
}
