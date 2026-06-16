// Native port of the upstream `indel-stats` plugin (plugins/indel-stats.c) for
// its default (no filtering expression, no PED) mode. It computes indel summary
// numbers, a per-length distribution (DLEN), a variant-allele-frequency
// distribution (DVAF) and the mean minor-allele fraction at het indels (DFRAC,
// with its support counts NFRAC) from FORMAT/AD, and prints a tab-separated
// report. The VCF/BCF output is suppressed.
//
// A single -i/--include or -e/--exclude filter expression is supported as a
// site/per-sample pre-filter via the shared filter engine, matching upstream's
// filter_init/filter_test usage in indel-stats.c (a FORMAT expression filters
// per-sample, a site expression filters whole records). The curly-brace
// multi-threshold expansion is supported too: such an expression is expanded
// into one filter per element (and a cartesian product across multiple groups),
// each tallied into its own SN* / DVAF* / DLEN* / DFRAC* / NFRAC* report
// section, matching upstream's parse_filters().
//
// The PED-restricted de-novo mode (-p/--ped FILE) is supported: the trios
// resolved from the PED file restrict the stats to de-novo indels in each
// trio's child (Mendelian-error genotypes that introduce an indel allele not
// inherited from either parent), exactly as indel-stats.c's process_record()
// does in its PED branch. With -p the SN* "number of samples" column reports
// the trio count, npass counts sites carrying a DNM, and npass_gt counts the
// de-novo indel genotypes. The --alt2ref-DNM toggle (accept 0/1 + 1/1 -> 0/0 as
// a valid DNM) is honoured. The report can be written to a file with
// -o/--output FILE; the bytes are identical to the stdout form.
package bcftools

import (
	"fmt"
	"io"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() { registerNativePlugin("indel-stats", func() NativePlugin { return &indelStatsPlugin{} }) }

// indelStatsCounters holds one filter's accumulators, mirroring stats_t/
// flt_stats_t in indel-stats.c.
type indelStatsCounters struct {
	nvafBins              []uint32
	nlen                  []uint32
	nfrac                 []uint32
	dfrac                 []float64
	npassGT               uint32
	npass                 uint32
	nsites                uint32
	nins, ndel            uint32
	nframeshift, ninframe uint32
}

// indelStatsFilter is one expanded -i/-e threshold: its compiled filter, the
// DEF-line label and the accumulators.
type indelStatsFilter struct {
	filter   *pluginFilter // compiled filter; nil for the default "all"
	label    string        // DEF-line label: "all" or the expanded expression
	counters indelStatsCounters
}

// indelStatsPlugin implements the `indel-stats` plugin in its default mode.
type indelStatsPlugin struct {
	hdr    *vcf.Header
	csqTag string
	maxLen int
	nvaf   int

	// PED de-novo mode (-p): when trios is non-nil the stats are restricted to
	// de-novo indels in the child of each resolved trio, and the SN* "number of
	// samples" column reports the trio count instead. allowAlt2RefDNM mirrors
	// upstream's --alt2ref-DNM (accept 0/1 + 1/1 -> 0/0 as a valid DNM).
	trios           []trioStatsTrio
	allowAlt2RefDNM bool

	// filters holds one entry per expanded threshold (one for the default "all"
	// or a single -i/-e expression, N for a curly-brace expansion).
	filters []*indelStatsFilter

	out        io.Writer
	stderr     io.Writer
	argv       []string
	outputFile string // -o/--output FILE; "" or "-" means stdout
}

// newCounters allocates the per-filter accumulators sized to the plugin's
// nvaf/maxLen bin counts.
func (p *indelStatsPlugin) newCounters() indelStatsCounters {
	return indelStatsCounters{
		nvafBins: make([]uint32, p.nvaf),
		nlen:     make([]uint32, p.maxLen*2+1),
		nfrac:    make([]uint32, p.maxLen*2+1),
		dfrac:    make([]float64, p.maxLen*2+1),
	}
}

// SetStderr wires the host stderr writer the "Collecting data ..." note uses.
func (p *indelStatsPlugin) SetStderr(w io.Writer) { p.stderr = w }

// SuppressVCF reports true: indel-stats emits only its textual report.
func (p *indelStatsPlugin) SuppressVCF() bool { return true }

// SetStdout wires the host stdout writer the report is printed to.
func (p *indelStatsPlugin) SetStdout(w io.Writer) { p.out = w }

// SetArgv records the upstream-equivalent argv for the CMD report line.
func (p *indelStatsPlugin) SetArgv(argv []string) { p.argv = argv }

// RunStyle reports that indel-stats is a run()-style plugin.
func (p *indelStatsPlugin) RunStyle() bool { return true }

// FlagTakesValue reports whether one of indel-stats' value-taking flags consumes
// the following token, used by the host to split the input-file positional.
func (p *indelStatsPlugin) FlagTakesValue(flag string) bool {
	switch flag {
	case "-c", "--csq-tag", "-i", "--include", "-e", "--exclude", "-o", "--output",
		"-p", "--ped", "-r", "--regions", "-R", "--regions-file", "-t", "--targets",
		"-T", "--targets-file", "-v", "--verbosity", "--max-len", "--nvaf":
		return true
	}
	return false
}

// Name returns the plugin name.
func (p *indelStatsPlugin) Name() string { return "indel-stats" }

// RegionTargetCaps opts indel-stats into the shared -r/-R/-t/-T region/target
// filter, applied to the records before the indel statistics are tallied.
func (p *indelStatsPlugin) RegionTargetCaps() regionTargetCaps { return allRegionTargetCaps }

// About returns the one-line description, matching indel-stats.c about().
func (p *indelStatsPlugin) About() string {
	return "Calculate indel stats scanning over a range of thresholds simultaneously.\n"
}

// Parallel reports false: the accumulators are updated serially.
func (p *indelStatsPlugin) Parallel() bool { return false }

// Init parses the supported options (-i/-e, -p, --alt2ref-DNM, -c, --max-len,
// --nvaf, -o, -v), resolving the PED trios for the de-novo mode when -p is given.
func (p *indelStatsPlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	p.csqTag = "CSQ"
	p.maxLen = 20
	p.nvaf = 20
	var filterExpr string
	var filterExclude, haveFilter bool
	var pedFile string
	for i := 0; i < len(args); i++ {
		a := args[i]
		next := func() (string, error) {
			if i+1 >= len(args) {
				return "", fmt.Errorf("indel-stats: %s requires an argument", a)
			}
			i++
			return args[i], nil
		}
		switch a {
		case "-i", "--include", "-e", "--exclude":
			if haveFilter {
				return nil, fmt.Errorf("indel-stats: only one of -i/--include or -e/--exclude can be given")
			}
			v, err := next()
			if err != nil {
				return nil, err
			}
			filterExpr = v
			filterExclude = a == "-e" || a == "--exclude"
			haveFilter = true
		case "-p", "--ped":
			v, err := next()
			if err != nil {
				return nil, err
			}
			pedFile = v
		case "-o", "--output":
			v, err := next()
			if err != nil {
				return nil, err
			}
			p.outputFile = v
		case "--alt2ref-DNM":
			p.allowAlt2RefDNM = true
		case "-c", "--csq-tag":
			v, err := next()
			if err != nil {
				return nil, err
			}
			p.csqTag = v
		case "--max-len":
			v, err := next()
			if err != nil {
				return nil, err
			}
			n, err := strconv.Atoi(v)
			if err != nil || n <= 0 {
				return nil, fmt.Errorf("indel-stats: could not parse: --max-len %s", v)
			}
			p.maxLen = n
		case "--nvaf":
			v, err := next()
			if err != nil {
				return nil, err
			}
			n, err := strconv.Atoi(v)
			if err != nil {
				return nil, fmt.Errorf("indel-stats: could not parse: --nvaf %s", v)
			}
			// Mirror upstream's (buggy) validation: it rejects any explicit
			// --nvaf outside the [0,1] interval even though the value is used as
			// a bin count, so only --nvaf 0 or 1 are accepted.
			if n < 0 || n > 1 {
				return nil, fmt.Errorf("indel-stats: expected value from the interval [0,1] with --nvaf")
			}
			p.nvaf = n
		case "-v", "--verbosity":
			if _, err := next(); err != nil {
				return nil, err
			}
		default:
			return nil, fmt.Errorf("indel-stats: unsupported option %q", a)
		}
	}
	p.hdr = hdr

	// Resolve the PED trios before the filters, matching indel-stats.c
	// init_data() which calls parse_ped() (and prints its stderr note) ahead of
	// parse_filters(). Unlike trio-stats.c, indel-stats.c's parse_ped does NOT
	// deduplicate trios or reject a child listed twice; it simply appends every
	// resolvable row (father, mother AND child all present in the header) and
	// sorts by minimum sample index.
	if pedFile != "" {
		trios, err := parseIndelStatsPED(pedFile, sampleIndex(hdr))
		if err != nil {
			return nil, err
		}
		if p.stderr != nil {
			fmt.Fprintf(p.stderr, "Identified %d complete trios in the VCF file\n", len(trios))
		}
		p.trios = trios
	}

	// Expand any curly-brace multi-threshold list into N concrete expressions,
	// matching upstream parse_filters(). When a filter is given, upstream always
	// prints the "Collecting data for N filtering expressions" note (N may be 0
	// when the braces collapse, e.g. "GQ>{}").
	var exprs []string
	if haveFilter {
		var err error
		exprs, err = expandPluginFilterExpr(filterExpr)
		if err != nil {
			return nil, fmt.Errorf("indel-stats: %s", err)
		}
		if p.stderr != nil {
			fmt.Fprintf(p.stderr, "Collecting data for %d filtering expressions\n", len(exprs))
		}
	}

	if len(exprs) == 0 {
		// No filter, or a brace list that collapsed to nothing: the single
		// default "all" filter.
		p.filters = []*indelStatsFilter{{label: "all", counters: p.newCounters()}}
		return hdr, nil
	}
	p.filters = make([]*indelStatsFilter, len(exprs))
	for i, expr := range exprs {
		f, err := newPluginFilterWithHeader(expr, filterExclude, hdr)
		if err != nil {
			return nil, fmt.Errorf("indel-stats: %w", err)
		}
		p.filters[i] = &indelStatsFilter{filter: f, label: pluginExprLabel(expr), counters: p.newCounters()}
	}
	return hdr, nil
}

// len2bin maps a net indel length to a DLEN bin index, clamping to the extremes.
func (p *indelStatsPlugin) len2bin(l int) int {
	if l < -p.maxLen {
		return 0
	}
	if l > p.maxLen {
		return 2 * p.maxLen
	}
	return p.maxLen + l
}

// vaf2bin maps a variant allele fraction in [0,1] to a DVAF bin index.
func (p *indelStatsPlugin) vaf2bin(vaf float64) int {
	b := int(vaf * float64(p.nvaf-1))
	if b < 0 {
		b = 0
	}
	if b >= p.nvaf {
		b = p.nvaf - 1
	}
	return b
}

// Process accumulates indel statistics for one record across every expanded
// filter, mirroring indel-stats.c run()'s loop over filters which calls
// process_record() once per filter.
func (p *indelStatsPlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	for _, flt := range p.filters {
		p.processOne(v, flt)
	}
	return nil, nil
}

// processOne tallies one record into one filter's accumulators, mirroring
// indel-stats.c process_record() in the default (no PED) path. Sites without any
// indel allele are skipped (matching the run-loop pre-filter), but nsites counts
// every indel-bearing record reaching the accumulator.
func (p *indelStatsPlugin) processOne(v *vcf.Variant, flt *indelStatsFilter) {
	c := &flt.counters
	if variantTypeMask(v)&vtINDEL == 0 {
		return
	}
	c.nsites++
	// Upstream increments nsites for every indel site (the run-loop pre-filter)
	// and only then applies the -i/-e filter; a site that fails the filter is
	// counted in nsites but contributes nothing else.
	passSite, smplPass := flt.filter.testSamples(v)
	if !passSite {
		return
	}
	nAllele := len(v.Alt) + 1

	haveGT := len(v.Samples) > 0 && formatHasTag(v, "GT")
	var adVals [][]int
	nad1 := 0
	if haveGT {
		if formatHasTag(v, "AD") {
			adVals = parseADAll(v, nAllele)
			nad1 = nAllele
		}
	}

	starAllele := starAlleleIndex(v)

	// PED de-novo mode: restrict the stats to de-novo indels in each trio's
	// child, mirroring indel-stats.c process_record()'s `args->ngt>0 && ntrio`
	// branch. A site with no DNM contributes only nsites (already counted) — the
	// CSQ / nins / ndel / npass tallies below are skipped via the early return.
	if haveGT && len(p.trios) > 0 {
		isDNM := false
		for ti := range p.trios {
			// In PED mode a trio is processed only when all three members pass the
			// per-sample filter mask. testSamples has already folded the
			// INCLUDE/EXCLUDE inversion into smplPass (true == member retained), so
			// a trio passes iff its three mask bits are all true — matching
			// indel-stats.c's per-trio pass bookkeeping for both logics.
			trio := p.trios[ti]
			if smplPass != nil {
				if !(smplPass[trio.child] && smplPass[trio.father] && smplPass[trio.mother]) {
					continue
				}
			}
			alsChild, kc := parseGenotypeAlleles(v, trio.child)
			if kc == gtMissing {
				continue
			}
			alsFather, kf := parseGenotypeAlleles(v, trio.father)
			if kf == gtMissing {
				continue
			}
			alsMother, km := parseGenotypeAlleles(v, trio.mother)
			if km == gtMissing {
				continue
			}

			// Is it a DNM? Same Mendelian logic as indel-stats.c.
			if !p.allowAlt2RefDNM && alsChild[0] == 0 && alsChild[1] == 0 {
				continue
			}
			if (alsChild[0] == alsFather[0] || alsChild[0] == alsFather[1]) &&
				(alsChild[1] == alsMother[0] || alsChild[1] == alsMother[1]) {
				continue
			}
			if (alsChild[1] == alsFather[0] || alsChild[1] == alsFather[1]) &&
				(alsChild[0] == alsMother[0] || alsChild[0] == alsMother[1]) {
				continue
			}
			if alsChild[0] == starAllele || alsChild[1] == starAllele {
				continue
			}
			if alsFather[0] == starAllele || alsFather[1] == starAllele {
				continue
			}
			if alsMother[0] == starAllele || alsMother[1] == starAllele {
				continue
			}

			childIsIndel := altVariantType(v, alsChild[0])&vtINDEL != 0 || altVariantType(v, alsChild[1])&vtINDEL != 0
			if !p.allowAlt2RefDNM {
				if !childIsIndel {
					continue
				}
			} else if !childIsIndel &&
				altVariantType(v, alsFather[0])&vtINDEL == 0 &&
				altVariantType(v, alsFather[1])&vtINDEL == 0 &&
				altVariantType(v, alsMother[0])&vtINDEL == 0 &&
				altVariantType(v, alsMother[1])&vtINDEL == 0 {
				continue // not an indel, in any sample
			}

			if childIsIndel && nad1 > 0 {
				p.updateIndelStats(c, v, trio.child, alsChild, adVals, starAllele)
			}
			c.npassGT++
			isDNM = true
		}
		if !isDNM {
			return
		}
	} else if haveGT && nad1 > 0 {
		for i := range v.Samples {
			if smplPass != nil && !smplPass[i] {
				continue
			}
			als, kind := parseGenotypeAlleles(v, i)
			if kind == gtMissing {
				continue
			}
			if altVariantType(v, als[0])&vtINDEL == 0 && altVariantType(v, als[1])&vtINDEL == 0 {
				continue
			}
			p.updateIndelStats(c, v, i, als, adVals, starAllele)
			c.npassGT++
		}
	}

	if csq, ok := v.Info[p.csqTag]; ok {
		if strings.Contains(csq, "inframe") {
			c.ninframe++
		}
		if strings.Contains(csq, "frameshift") {
			c.nframeshift++
		}
	}

	for i := 1; i < nAllele; i++ {
		if altVariantType(v, i)&vtINDEL == 0 {
			continue
		}
		n := indelAlleleLen(v, i)
		if n < 0 {
			c.ndel++
		} else if n > 0 {
			c.nins++
		}
		if !haveGT || nad1 == 0 {
			bin := p.len2bin(n)
			if bin >= 0 {
				c.nlen[bin]++
			}
		}
	}
	c.npass++
}

// updateIndelStats records the VAF bin, length bins and DFRAC contribution for a
// single sample's indel genotype into the given filter's counters, mirroring
// indel-stats.c update_indel_stats().
func (p *indelStatsPlugin) updateIndelStats(c *indelStatsCounters, v *vcf.Variant, ismpl int, als [2]int, adVals [][]int, starAllele int) {
	ad := adVals[ismpl]
	ntot := 0
	for _, x := range ad {
		if x < 0 {
			continue // missing
		}
		ntot += x
	}
	if ntot == 0 {
		return
	}
	al0, al1 := als[0], als[1]
	if altVariantType(v, al0)&vtINDEL == 0 {
		if altVariantType(v, al1)&vtINDEL == 0 {
			return
		}
		al0, al1 = als[1], als[0]
	} else if altVariantType(v, al1)&vtINDEL != 0 && al0 != al1 {
		// Both alleles are indels: select the more frequent one as al0.
		if ad[al0] < ad[al1] {
			al0, al1 = als[1], als[0]
		}
		// Record the length of the less-frequent indel allele too.
		bin := p.len2bin(indelAlleleLen(v, al1))
		if bin >= 0 {
			c.nlen[bin]++
		}
	}

	vaf := float64(ad[al0]) / float64(ntot)
	c.nvafBins[p.vaf2bin(vaf)]++

	lenBin := p.len2bin(indelAlleleLen(v, al0))
	if lenBin < 0 {
		return
	}
	c.nlen[lenBin]++

	if al0 != al1 {
		tot := ad[al0] + ad[al1]
		if tot != 0 {
			c.nfrac[lenBin]++
			c.dfrac[lenBin] += float64(ad[al0]) / float64(tot)
		}
	}
}

// Destroy prints the report, mirroring indel-stats.c report_stats(). With
// -o/--output the report is written to that file instead of stdout (the bytes
// are identical; the CMD line echoes the verbatim argv in both cases).
func (p *indelStatsPlugin) Destroy() error {
	if p.out == nil && (p.outputFile == "" || p.outputFile == "-") {
		return nil
	}
	fp, closeFn, err := statsReportWriter(p.outputFile, p.out)
	if err != nil {
		return fmt.Errorf("indel-stats: %w", err)
	}
	defer closeFn()
	fmt.Fprint(fp, indelStatsHeaderFor(p.nvaf, p.maxLen))
	fmt.Fprintf(fp, "CMD\t%s\n", strings.Join(p.argv, " "))
	for i, flt := range p.filters {
		fmt.Fprintf(fp, "DEF\tFLT%d\t%s\n", i, flt.label)
	}

	// Upstream's SN* "number of samples" column reports the trio count in PED
	// mode (`args->ntrio ? args->ntrio : args->nsmpl`).
	nsmp := len(p.hdr.Samples)
	if len(p.trios) > 0 {
		nsmp = len(p.trios)
	}
	var b strings.Builder
	for i, flt := range p.filters {
		c := &flt.counters
		b.WriteString(fmt.Sprintf("SN%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
			i, nsmp, c.nsites, c.npass, c.npassGT, c.nins, c.ndel, c.nframeshift, c.ninframe))

		b.WriteString(fmt.Sprintf("DVAF%d\t%d", i, p.nvaf))
		for _, x := range c.nvafBins {
			b.WriteString(fmt.Sprintf("\t%d", x))
		}
		b.WriteByte('\n')

		b.WriteString(fmt.Sprintf("DLEN%d\t%d", i, p.maxLen))
		for _, x := range c.nlen {
			b.WriteString(fmt.Sprintf("\t%d", x))
		}
		b.WriteByte('\n')

		b.WriteString(fmt.Sprintf("DFRAC%d\t%d", i, p.maxLen))
		for j := range c.dfrac {
			if c.nfrac[j] != 0 {
				b.WriteString(fmt.Sprintf("\t%.2f", c.dfrac[j]/float64(c.nfrac[j])))
			} else {
				b.WriteString("\t.")
			}
		}
		b.WriteByte('\n')

		b.WriteString(fmt.Sprintf("NFRAC%d\t%d", i, p.maxLen))
		for _, x := range c.nfrac {
			b.WriteString(fmt.Sprintf("\t%d", x))
		}
		b.WriteByte('\n')
	}
	fp.Write([]byte(b.String()))
	return nil
}

// indelStatsHeaderFor renders the report header, which embeds the NVAF/MAX_LEN
// column ranges, matching indel-stats.c report_stats().
func indelStatsHeaderFor(nvaf, maxLen int) string {
	var b strings.Builder
	b.WriteString("# CMD line shows the command line used to generate this output\n")
	b.WriteString("# DEF lines define expressions for all tested thresholds\n")
	b.WriteString("# SN* summary number for every threshold:\n")
	i := 0
	w := func(s string) { i++; b.WriteString(fmt.Sprintf("#   %d) %s\n", i, s)) }
	w("SN*, filter id")
	w("number of samples (or trios with -p)")
	w("number of indel sites total")
	w("number of indel sites that pass the filter (and, with -p, have a de novo indel)")
	w("number of indel genotypes that pass the filter (and, with -p, are de novo)")
	w("number of insertions (site-wise, not genotype-wise)")
	w("number of deletions (site-wise, not genotype-wise)")
	w("number of frameshifts (site-wise, not genotype-wise)")
	w("number of inframe indels (site-wise, not genotype-wise)")
	b.WriteString("#\n")
	i = 0
	b.WriteString("# DVAF* lines report indel variant allele frequency (VAF) distribution for every threshold,\n")
	b.WriteString("#   k-th bin corresponds to the frequency k/(nVAF-1):\n")
	w("DVAF*, filter id")
	w("nVAF, number of bins which split the [0,1] VAF interval.")
	b.WriteString(fmt.Sprintf("#   %d-%d) counts of indel genotypes in the VAF bin. For non-reference hets, the VAF of the less supported allele is recorded\n", i+1, i+nvaf))
	b.WriteString("#\n")
	i = 0
	b.WriteString("# DLEN* lines report indel length distribution for every threshold. When genotype fields are available,\n")
	b.WriteString("#   the counts correspond to the number of genotypes, otherwise the number of sites are given.\n")
	b.WriteString("#   The k-th bin corresponds to the indel size k-MAX_LEN, negative for deletions, positive for insertions.\n")
	b.WriteString("#   The first/last bin contains also all deletions/insertions larger than MAX_LEN:\n")
	w("DLEN*, filter id")
	w("maximum indel length")
	b.WriteString(fmt.Sprintf("#   %d-%d) counts of indel lengths (-max,..,0,..,max), all unique alleles in a genotype are recorded (alt hets increase the counters 2x, alt homs 1x)\n", i+1, i+maxLen*2+1))
	b.WriteString("#\n")
	i = 0
	b.WriteString("# DFRAC* lines report the mean minor allele fraction at HET indel genotypes as a function of indel size.\n")
	b.WriteString("#   The format is the same as for DLEN:\n")
	w("DFRAC*, filter id")
	w("maximum indel length")
	b.WriteString(fmt.Sprintf("#   %d-%d) mean fraction at indel lengths (-max,..,0,..,max)\n", i+1, i+maxLen*2+1))
	b.WriteString("#\n")
	i = 0
	b.WriteString("# NFRAC* lines report the number of indels informing the DFRAC distribution.\n")
	w("NFRAC*, filter id")
	w("maximum indel length")
	b.WriteString(fmt.Sprintf("#   %d-%d) counts at indel lengths (-max,..,0,..,max)\n", i+1, i+maxLen*2+1))
	b.WriteString("#\n")
	return b.String()
}

// parseADAll returns per-sample FORMAT/AD as integer slices of length nAllele.
// Missing values (".") become -1. A sample missing the AD field yields a slice
// of -1s.
func parseADAll(v *vcf.Variant, nAllele int) [][]int {
	out := make([][]int, len(v.Samples))
	for i := range v.Samples {
		row := make([]int, nAllele)
		for k := range row {
			row[k] = -1
		}
		s, ok := v.Samples[i].Data["AD"]
		if ok && s != "" && s != "." {
			for k, tok := range strings.Split(s, ",") {
				if k >= nAllele {
					break
				}
				if tok == "." || tok == "" {
					row[k] = -1
					continue
				}
				n, err := strconv.Atoi(tok)
				if err != nil {
					row[k] = -1
					continue
				}
				row[k] = n
			}
		}
		out[i] = row
	}
	return out
}
