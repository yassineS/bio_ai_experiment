// Native port of the upstream `trio-stats` plugin (plugins/trio-stats.c) for its
// default (no filtering expression) mode. Given a PED file (or a -P P,F,M trio)
// it reports, per trio, transmission and Mendelian-consistency statistics: the
// number of valid trio genotypes, non-reference trio GTs, DNMs / Mendelian
// errors, novel singleton alleles in the child, untransmitted and transmitted
// trio singletons, transitions / transversions and their ratio, plus homozygous
// and recurrent DNMs. The VCF/BCF output is suppressed; only the textual report
// is written.
//
// trio-stats is a run()-style plugin (its options precede the input file, no
// `--` separator) and is inherently stateful, so it runs through the serial
// pipeline and accumulates per-trio counters across the whole record stream.
//
// A single -i/--include or -e/--exclude filter expression is supported as a
// site/per-trio pre-filter via the shared filter engine, matching trio-stats.c's
// filter_init/filter_test usage: a site expression gates whole records, while a
// FORMAT expression is resolved per sample and then folded into a per-trio
// verdict (include keeps a trio only when all three members match; exclude keeps
// a trio only when none of the three match). The curly-brace multi-threshold
// expansion, the -a/--alt-trios accounting, region/target index-jumping and the
// -o file output remain unsupported. The default single filter — the common
// case — is implemented faithfully, including the -d/--debug MERR and
// TRANSMITTED line emission and the -v/--verbosity passthrough.
package bcftools

import (
	"fmt"
	"io"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() { registerNativePlugin("trio-stats", func() NativePlugin { return &trioStatsPlugin{} }) }

// trioStatsTrio is one (child, father, mother) tuple, resolved to sample
// indices, mirroring trio_t in trio-stats.c (idx[iCHILD/iFATHER/iMOTHER]).
type trioStatsTrio struct {
	child, father, mother int
}

// trioStatsCounters accumulates the per-trio statistics, mirroring
// trio_stats_t in trio-stats.c.
type trioStatsCounters struct {
	npass      uint32 // valid trio genotypes (all members non-missing)
	nnonRef    uint32 // non-reference trio GTs
	nmendelErr uint32 // DNMs / Mendelian errors
	nnovel     uint32 // novel singleton allele in the child
	nsingleton uint32 // untransmitted trio singletons
	ndoubleton uint32 // transmitted trio singletons
	nts, ntv   uint32 // transitions / transversions
	ndnmRecur  uint32 // recurrent DNMs / Mendelian errors
	ndnmHom    uint32 // homozygous DNMs / Mendelian errors
}

// trioStatsPlugin implements the `trio-stats` plugin in its default mode.
type trioStatsPlugin struct {
	hdr     *vcf.Header
	trios   []trioStatsTrio
	stats   []trioStatsCounters
	verbose int // VERBOSE_MENDEL (1) | VERBOSE_TRANSMITTED (2)
	out     io.Writer
	argv    []string

	filter    *pluginFilter // compiled -i/-e pre-filter, nil for the default "all"
	exprLabel string        // DEF-line label: "all" by default, else the expression
	stderr    io.Writer
}

// SetStderr wires the host stderr writer the "Collecting data ..." note uses.
func (p *trioStatsPlugin) SetStderr(w io.Writer) { p.stderr = w }

const (
	trioVerboseMendel      = 1
	trioVerboseTransmitted = 2
)

// SuppressVCF reports true: trio-stats emits only its textual report.
func (p *trioStatsPlugin) SuppressVCF() bool { return true }

// SetStdout wires the host stdout writer the report is printed to.
func (p *trioStatsPlugin) SetStdout(w io.Writer) { p.out = w }

// SetArgv records the upstream-equivalent argv for the CMD report line.
func (p *trioStatsPlugin) SetArgv(argv []string) { p.argv = argv }

// RunStyle reports that trio-stats is a run()-style plugin.
func (p *trioStatsPlugin) RunStyle() bool { return true }

// FlagTakesValue reports whether one of trio-stats' value-taking flags consumes
// the following token, used by the host to split the input-file positional.
func (p *trioStatsPlugin) FlagTakesValue(flag string) bool {
	switch flag {
	case "-a", "--alt-trios", "-d", "--debug", "-e", "--exclude",
		"-i", "--include", "-o", "--output", "-p", "--ped", "-P", "--pfm",
		"-r", "--regions", "-R", "--regions-file", "-t", "--targets",
		"-T", "--targets-file", "-v", "--verbosity":
		return true
	}
	return false
}

// Name returns the plugin name.
func (p *trioStatsPlugin) Name() string { return "trio-stats" }

// RegionTargetCaps opts trio-stats into the shared -r/-R/-t/-T region/target
// filter, applied to the records before the per-trio statistics are tallied.
func (p *trioStatsPlugin) RegionTargetCaps() regionTargetCaps { return allRegionTargetCaps }

// About returns the one-line description, matching trio-stats.c about().
func (p *trioStatsPlugin) About() string {
	return "Calculate transmission rate and other stats in trio children.\n"
}

// Parallel reports false: the per-trio accumulators are updated serially.
func (p *trioStatsPlugin) Parallel() bool { return false }

// Init parses the supported options (-p/-P, -d, -v), resolves the trios against
// the header and rejects the modes that need htslib machinery (filters,
// alt-trios accounting, region jumps, file output).
func (p *trioStatsPlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	p.exprLabel = "all"
	var pedFile, pfm string
	var filterExpr string
	var filterExclude, haveFilter bool
	for i := 0; i < len(args); i++ {
		a := args[i]
		switch a {
		case "-i", "--include", "-e", "--exclude":
			if haveFilter {
				return nil, fmt.Errorf("trio-stats: only one of -i/--include or -e/--exclude can be given")
			}
			if i+1 >= len(args) {
				return nil, fmt.Errorf("trio-stats: %s requires a value", a)
			}
			i++
			if strings.ContainsRune(args[i], '{') {
				return nil, fmt.Errorf("trio-stats: the curly-brace multi-threshold filter expansion is not supported by the native plugin")
			}
			filterExpr = args[i]
			filterExclude = a == "-e" || a == "--exclude"
			haveFilter = true
		case "-a", "--alt-trios":
			return nil, fmt.Errorf("trio-stats: the --alt-trios accounting (-a) is not supported by the native plugin")
		case "-o", "--output":
			return nil, fmt.Errorf("trio-stats: writing to a file (-o) is not supported by the native plugin; use stdout")
		case "-p", "--ped":
			if i+1 >= len(args) {
				return nil, fmt.Errorf("trio-stats: %s requires a value", a)
			}
			i++
			pedFile = args[i]
		case "-P", "--pfm":
			if i+1 >= len(args) {
				return nil, fmt.Errorf("trio-stats: %s requires a value", a)
			}
			i++
			pfm = args[i]
		case "-d", "--debug":
			if i+1 >= len(args) {
				return nil, fmt.Errorf("trio-stats: %s requires a value", a)
			}
			i++
			for _, feat := range strings.Split(args[i], ",") {
				switch strings.ToLower(strings.TrimSpace(feat)) {
				case "mendel-errors":
					p.verbose |= trioVerboseMendel
				case "transmitted":
					p.verbose |= trioVerboseTransmitted
				default:
					return nil, fmt.Errorf("trio-stats: the argument %q to option --debug is not recognised", feat)
				}
			}
		case "-v", "--verbosity":
			if i+1 < len(args) {
				i++
			}
		default:
			return nil, fmt.Errorf("trio-stats: unsupported option %q", a)
		}
	}
	if pedFile == "" && pfm == "" {
		return nil, fmt.Errorf("trio-stats: missing the -p or -P option")
	}
	if haveFilter {
		f, err := newPluginFilterWithHeader(filterExpr, filterExclude, hdr)
		if err != nil {
			return nil, fmt.Errorf("trio-stats: %w", err)
		}
		p.filter = f
		p.exprLabel = filterExpr
		if p.stderr != nil {
			fmt.Fprint(p.stderr, "Collecting data for 1 filtering expressions\n")
		}
	}
	p.hdr = hdr

	idx := sampleIndex(hdr)
	if pedFile != "" {
		trios, err := parseTrioStatsPED(pedFile, idx)
		if err != nil {
			return nil, err
		}
		p.trios = trios
	} else {
		t, err := parseTrioPFM(pfm, idx)
		if err != nil {
			return nil, fmt.Errorf("trio-stats: %w", err)
		}
		p.trios = []trioStatsTrio{t}
	}
	p.stats = make([]trioStatsCounters, len(p.trios))

	// Upstream prints the comment block and the CMD line in init_data(), before
	// any record is processed, so the streamed MERR / TRANSMITTED debug lines
	// (emitted during Process) appear after this header and before the final
	// DEF / FLT report written in Destroy.
	if p.out != nil {
		fmt.Fprint(p.out, trioStatsHeader)
		fmt.Fprintf(p.out, "CMD\t%s\n", strings.Join(p.argv, " "))
	}
	return hdr, nil
}

// trioFilterPass evaluates the -i/-e filter for one record and returns the
// per-trio pass mask (nil when no per-trio filtering applies, i.e. every trio is
// processed) and the site-level verdict (false ⇒ skip the whole record). It
// reproduces, verbatim, the FLT_INCLUDE / FLT_EXCLUDE bookkeeping in
// trio-stats.c process_record(): an INCLUDE trio passes only when all three
// members match; an EXCLUDE trio passes only when none of the three match; and
// the site is skipped when no trio passes.
func (p *trioStatsPlugin) trioFilterPass(v *vcf.Variant) (trioPass []bool, passSite bool) {
	siteMatch, mask, exclude := p.filter.rawSamples(v)
	if exclude {
		if !siteMatch {
			return nil, true // nothing matched: every trio passes
		}
		if mask == nil {
			return nil, false // site-only exclude matched: drop the record
		}
		pass := make([]bool, len(p.trios))
		any := false
		for ti, trio := range p.trios {
			ok := !(mask[trio.child] || mask[trio.father] || mask[trio.mother])
			pass[ti] = ok
			if ok {
				any = true
			}
		}
		if !any {
			return nil, false
		}
		return pass, true
	}
	// FLT_INCLUDE.
	if !siteMatch {
		return nil, false
	}
	if mask == nil {
		return nil, true // site-only include passed: every trio processed
	}
	pass := make([]bool, len(p.trios))
	any := false
	for ti, trio := range p.trios {
		ok := mask[trio.child] && mask[trio.father] && mask[trio.mother]
		pass[ti] = ok
		if ok {
			any = true
		}
	}
	if !any {
		return nil, false
	}
	return pass, true
}

// Process accumulates per-trio statistics for one record, mirroring
// process_record() in trio-stats.c with the default "all" filter and no
// --alt-trios accounting.
func (p *trioStatsPlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	trioPass, passSite := p.trioFilterPass(v)
	if !passSite {
		return nil, nil
	}

	nAllele := len(v.Alt) + 1
	ac := computeACWithRef(v, nAllele)

	// For ts/tv: numeric code of the reference allele, -1 for insertions
	// (multi-base REF). Matches `!rec->d.allele[0][1] ? bcf_acgt2int(...) : -1`.
	ref := -1
	if len(v.Ref) == 1 {
		ref = acgt2int(v.Ref[0])
	}
	starAllele := starAlleleIndex(v)

	acTrio := make([]int, nAllele)
	for ti := range p.trios {
		if trioPass != nil && !trioPass[ti] {
			continue
		}
		trio := p.trios[ti]
		stats := &p.stats[ti]

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

		stats.npass++

		als := [6]int{alsChild[0], alsChild[1], alsFather[0], alsFather[1], alsMother[0], alsMother[1]}

		hasStar, hasNonRef := false, false
		for j := range acTrio {
			acTrio[j] = 0
		}
		for _, a := range als {
			if a == starAllele {
				hasStar = true
				continue
			}
			if a != 0 {
				hasNonRef = true
			}
			if a >= 0 && a < nAllele {
				acTrio[a]++
			}
		}
		if !hasNonRef {
			continue
		}
		stats.nnonRef++

		// ts/tv: handles HetAA genotypes too.
		if ref != -1 {
			hasTs, hasTv := false, false
			for _, a := range als {
				if a == 0 || a == starAllele {
					continue
				}
				if a >= nAllele {
					return nil, fmt.Errorf("trio-stats: the GT index is out of range at %s:%d", v.Chrom, v.Pos)
				}
				// Only single-base ALT alleles contribute (rec->d.allele[a][1]==0).
				if a-1 >= len(v.Alt) || len(v.Alt[a-1]) != 1 {
					continue
				}
				alt := acgt2int(v.Alt[a-1][0])
				if abs(ref-alt) == 2 {
					hasTs = true
				} else {
					hasTv = true
				}
			}
			if hasTs {
				stats.nts++
			}
			if hasTv {
				stats.ntv++
			}
		}

		// Skip the remaining stats if the star allele is present (it was already
		// handled at the primary record).
		if hasStar {
			continue
		}

		// Detect Mendelian errors.
		a0F := alsChild[0] == alsFather[0] || alsChild[0] == alsFather[1]
		a1M := alsChild[1] == alsMother[0] || alsChild[1] == alsMother[1]
		if !a0F || !a1M {
			a0M := alsChild[0] == alsMother[0] || alsChild[0] == alsMother[1]
			a1F := alsChild[1] == alsFather[0] || alsChild[1] == alsFather[1]
			if !a0M || !a1F {
				stats.nmendelErr++

				dnmHom := false
				if alsChild[0] == alsChild[1] {
					stats.ndnmHom++
					dnmHom = true
				}

				var culprit int
				switch {
				case !a0F && !a0M:
					culprit = alsChild[0]
				case !a1F && !a1M:
					culprit = alsChild[1]
				case ac[alsChild[0]] < ac[alsChild[1]]:
					culprit = alsChild[0]
				default:
					culprit = alsChild[1]
				}

				dnmRecurrent := false
				if (!dnmHom && ac[culprit] > 1) || (dnmHom && ac[culprit] > 2) {
					stats.ndnmRecur++
					dnmRecurrent = true
				}

				if p.verbose&trioVerboseMendel != 0 {
					homStr, recStr := "-", "-"
					if dnmHom {
						homStr = "HOM"
					}
					if dnmRecurrent {
						recStr = "RECURRENT"
					}
					fmt.Fprintf(p.out, "MERR\t%s\t%d\t%s\t%s\t%s\t%s\t%s\n",
						v.Chrom, v.Pos,
						p.hdr.Samples[trio.child], p.hdr.Samples[trio.father], p.hdr.Samples[trio.mother],
						homStr, recStr)
				}
			}
		}

		// Singleton / doubleton classification.
		for j := 0; j < nAllele; j++ {
			if acTrio[j] == 0 {
				continue
			}
			if acTrio[j] == 1 { // singleton (in parent) or novel (in child)
				if alsChild[0] == j || alsChild[1] == j {
					stats.nnovel++
				} else {
					stats.nsingleton++
					if p.verbose&trioVerboseTransmitted != 0 {
						fmt.Fprintf(p.out, "TRANSMITTED\t%s\t%d\t%s\t%s\t%s\tNO\n",
							v.Chrom, v.Pos,
							p.hdr.Samples[trio.child], p.hdr.Samples[trio.father], p.hdr.Samples[trio.mother])
					}
				}
			} else if acTrio[j] == 2 { // possibly a doubleton
				if (alsChild[0] != j && alsChild[1] != j) || (alsChild[0] == j && alsChild[1] == j) {
					continue
				}
				if (alsFather[0] == j && alsFather[1] == j) || (alsMother[0] == j && alsMother[1] == j) {
					continue
				}
				stats.ndoubleton++
				if p.verbose&trioVerboseTransmitted != 0 {
					fmt.Fprintf(p.out, "TRANSMITTED\t%s\t%d\t%s\t%s\t%s\tYES\n",
						v.Chrom, v.Pos,
						p.hdr.Samples[trio.child], p.hdr.Samples[trio.father], p.hdr.Samples[trio.mother])
				}
			}
		}
	}
	return nil, nil
}

// Destroy prints the report, mirroring report_stats() in trio-stats.c.
func (p *trioStatsPlugin) Destroy() error {
	if p.out == nil {
		return nil
	}
	fp := p.out
	// The comment block and CMD line were already printed in Init (matching
	// upstream's init_data); Destroy emits only the DEF / FLT report, after any
	// streamed MERR / TRANSMITTED debug lines.
	fmt.Fprintf(fp, "DEF\tFLT0\t%s\n", p.exprLabel)
	for j := range p.trios {
		trio := p.trios[j]
		s := &p.stats[j]
		fmt.Fprintf(fp, "FLT0\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%d\t%d\n",
			p.hdr.Samples[trio.child], p.hdr.Samples[trio.father], p.hdr.Samples[trio.mother],
			s.npass, s.nnonRef, s.nmendelErr, s.nnovel, s.nsingleton, s.ndoubleton,
			s.nts, s.ntv, tstvStr(int(s.nts), int(s.ntv)), s.ndnmHom, s.ndnmRecur)
	}
	return nil
}

// trioStatsHeader is the verbatim comment block printed by trio-stats.c
// init_data() before the CMD line.
const trioStatsHeader = `# CMD line shows the command line used to generate this output
# DEF lines define expressions for all tested thresholds
# FLT* lines report numbers for every threshold and every trio:
#   1) filter id
#   2) child
#   3) father
#   4) mother
#   5) number of valid trio genotypes (all trio members pass filters, all non-missing)
#   6) number of non-reference trio GTs (at least one trio member carries an alternate allele)
#   7) number of DNMs/Mendelian errors
#   8) number of novel singleton alleles in the child (counted also as DNM / Mendelian error)
#   9) number of untransmitted trio singletons (one alternate allele present in one parent)
#   10) number of transmitted trio singletons (one alternate allele present in one parent and the child)
#   11) number of transitions, all distinct ALT alleles present in the trio are considered
#   12) number of transversions, all distinct ALT alleles present in the trio are considered
#   13) overall ts/tv, all distinct ALT alleles present in the trio are considered
#   14) number of homozygous DNMs/Mendelian errors (likely genotyping errors)
#   15) number of recurrent DNMs/Mendelian errors (non-inherited alleles present in other samples; counts GTs, not sites)
`
