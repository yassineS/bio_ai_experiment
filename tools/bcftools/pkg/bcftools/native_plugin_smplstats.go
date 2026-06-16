// Native port of the upstream `smpl-stats` plugin (plugins/smpl-stats.c) for its
// default (no filtering expression) mode. It computes basic per-sample genotype
// statistics — genotype-class counts, SNV/indel counts, singletons, ts/tv — plus
// a per-site summary, and prints a tab-separated report. The VCF/BCF output is
// suppressed.
//
// A single -i/--include or -e/--exclude filter expression is supported as a
// site/per-sample pre-filter via the shared filter engine, matching upstream's
// filter_init/filter_test usage in smpl-stats.c (a FORMAT expression filters
// per-sample, a site expression filters whole records). The curly-brace
// "{10,20,30}" multi-threshold expansion (which defines several FLT* filters at
// once) and the index/region jump options remain unsupported; the single
// default-or-explicit filter is the common case.
package bcftools

import (
	"fmt"
	"io"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() { registerNativePlugin("smpl-stats", func() NativePlugin { return &smplStatsPlugin{} }) }

// smplStats holds one accumulator (per sample, plus the per-site summary),
// matching stats_t in smpl-stats.c.
type smplStats struct {
	npass, nnonRef, nhomRR, nhomAA, nhemi, nhet  int
	nSNV, nIndel, nmissing, nsingleton, nts, ntv int
}

// smplStatsPlugin implements the `smpl-stats` plugin in its default mode.
type smplStatsPlugin struct {
	hdr       *vcf.Header
	stats     []smplStats // per-sample
	siteStats smplStats
	out       io.Writer
	stderr    io.Writer
	argv      []string

	filter    *pluginFilter // compiled -i/-e pre-filter, nil for the default "all"
	exprLabel string        // DEF-line label: "all" by default, else the expression
}

// SetStderr wires the host stderr writer the "Collecting data ..." note uses.
func (p *smplStatsPlugin) SetStderr(w io.Writer) { p.stderr = w }

// SuppressVCF reports true: smpl-stats emits only its textual report.
func (p *smplStatsPlugin) SuppressVCF() bool { return true }

// SetStdout wires the host stdout writer the report is printed to.
func (p *smplStatsPlugin) SetStdout(w io.Writer) { p.out = w }

// SetArgv records the upstream-equivalent argv for the CMD report line.
func (p *smplStatsPlugin) SetArgv(argv []string) { p.argv = argv }

// RunStyle reports that smpl-stats is a run()-style plugin (options precede the
// file, no `--` separator).
func (p *smplStatsPlugin) RunStyle() bool { return true }

// FlagTakesValue reports whether one of smpl-stats' value-taking flags consumes
// the following token, used by the host to split the input-file positional.
func (p *smplStatsPlugin) FlagTakesValue(flag string) bool {
	switch flag {
	case "-i", "--include", "-e", "--exclude", "-o", "--output",
		"-r", "--regions", "-R", "--regions-file", "-t", "--targets",
		"-T", "--targets-file", "-v", "--verbosity":
		return true
	}
	return false
}

// Name returns the plugin name.
func (p *smplStatsPlugin) Name() string { return "smpl-stats" }

// RegionTargetCaps opts smpl-stats into the shared -r/-R/-t/-T region/target
// filter, applied to the records before the per-sample statistics are tallied.
func (p *smplStatsPlugin) RegionTargetCaps() regionTargetCaps { return allRegionTargetCaps }

// About returns the one-line description, matching smpl-stats.c about().
func (p *smplStatsPlugin) About() string {
	return "Calculate basic per-sample stats scanning over a range of thresholds simultaneously.\n"
}

// Parallel reports false: the accumulators are updated serially.
func (p *smplStatsPlugin) Parallel() bool { return false }

// Init parses the supported options and rejects the filter/region modes that
// require htslib machinery the native pipeline does not provide.
func (p *smplStatsPlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	p.exprLabel = "all"
	var filterExpr string
	var filterExclude, haveFilter bool
	for i := 0; i < len(args); i++ {
		a := args[i]
		next := func() (string, error) {
			if i+1 >= len(args) {
				return "", fmt.Errorf("smpl-stats: %s requires an argument", a)
			}
			i++
			return args[i], nil
		}
		switch a {
		case "-i", "--include", "-e", "--exclude":
			if haveFilter {
				return nil, fmt.Errorf("smpl-stats: only one of -i/--include or -e/--exclude can be given")
			}
			v, err := next()
			if err != nil {
				return nil, err
			}
			if strings.ContainsRune(v, '{') {
				return nil, fmt.Errorf("smpl-stats: the curly-brace multi-threshold filter expansion is not supported by the native plugin")
			}
			filterExpr = v
			filterExclude = a == "-e" || a == "--exclude"
			haveFilter = true
		case "-o", "--output":
			return nil, fmt.Errorf("smpl-stats: writing to a file (-o) is not supported by the native plugin; use stdout")
		case "-v", "--verbosity":
			if i+1 < len(args) {
				i++
			}
		default:
			return nil, fmt.Errorf("smpl-stats: unsupported option %q", a)
		}
	}
	if haveFilter {
		f, err := newPluginFilterWithHeader(filterExpr, filterExclude, hdr)
		if err != nil {
			return nil, fmt.Errorf("smpl-stats: %w", err)
		}
		p.filter = f
		p.exprLabel = filterExpr
		if p.stderr != nil {
			fmt.Fprint(p.stderr, "Collecting data for 1 filtering expressions\n")
		}
	}
	p.hdr = hdr
	p.stats = make([]smplStats, len(hdr.Samples))
	return hdr, nil
}

// Process accumulates per-sample and per-site statistics for one record,
// mirroring smpl-stats.c process_record() with the default "all" filter.
func (p *smplStatsPlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	passSite, smplPass := p.filter.testSamples(v)
	if !passSite {
		return nil, nil
	}

	nAllele := len(v.Alt) + 1
	ac := computeACWithRef(v, nAllele)

	ref := -1
	if len(v.Ref) == 1 {
		ref = acgt2int(v.Ref[0])
	}
	starAllele := starAlleleIndex(v)

	var sitePass, siteSNV, siteIndel, siteHasTs, siteHasTv, siteSingleton int
	for i := range v.Samples {
		if smplPass != nil && !smplPass[i] {
			continue
		}
		stats := &p.stats[i]
		als, kind := parseGenotypeAlleles(v, i)
		if kind == gtMissing {
			stats.nmissing++
			continue
		}
		if kind == gtHemi {
			stats.nhemi++
		} else if als[0] != als[1] {
			stats.nhet++
		} else if als[0] == 0 {
			stats.nhomRR++
		} else {
			stats.nhomAA++
		}
		stats.npass++
		sitePass = 1

		hasNonRef := false
		for j := 0; j < 2; j++ {
			if als[j] == starAllele || als[j] == 0 {
				continue
			}
			hasNonRef = true
		}
		if !hasNonRef {
			continue
		}
		stats.nnonRef++

		hasTs, hasTv, hasSNV, hasIndel := false, false, false, false
		for j := 0; j < 2; j++ {
			a := als[j]
			if a == 0 || a == starAllele {
				continue
			}
			if a >= nAllele {
				return nil, fmt.Errorf("smpl-stats: GT index out of range at %s:%d", v.Chrom, v.Pos)
			}
			if a-1 < len(ac) && ac[a] == 1 {
				stats.nsingleton++
				siteSingleton = 1
			}
			vt := altVariantType(v, a)
			if vt&(vtSNP|vtMNP) != 0 {
				// ts/tv from the first differing base, mirroring upstream.
				alt := v.Alt[a-1]
				k := 0
				for k < len(v.Ref) && k < len(alt) {
					if v.Ref[k] == alt[k] {
						k++
						continue
					}
					altI := acgt2int(alt[k])
					if abs(ref-altI) == 2 {
						hasTs = true
					} else {
						hasTv = true
					}
					hasSNV = true
					k++
				}
			} else if vt == vtINDEL {
				hasIndel = true
			}
		}
		if hasTs {
			stats.nts++
			siteHasTs = 1
		}
		if hasTv {
			stats.ntv++
			siteHasTv = 1
		}
		if hasSNV {
			stats.nSNV++
			siteSNV = 1
		}
		if hasIndel {
			stats.nIndel++
			siteIndel = 1
		}
	}
	p.siteStats.npass += sitePass
	p.siteStats.nSNV += siteSNV
	p.siteStats.nIndel += siteIndel
	p.siteStats.nts += siteHasTs
	p.siteStats.ntv += siteHasTv
	p.siteStats.nsingleton += siteSingleton
	return nil, nil
}

// Destroy prints the report, mirroring smpl-stats.c report_stats().
func (p *smplStatsPlugin) Destroy() error {
	if p.out == nil {
		return nil
	}
	fp := p.out
	fmt.Fprint(fp, smplStatsHeader)
	fmt.Fprintf(fp, "CMD\t%s\n", strings.Join(p.argv, " "))
	fmt.Fprintf(fp, "DEF\tFLT0\t%s\n", p.exprLabel)
	for j := range p.stats {
		s := &p.stats[j]
		fmt.Fprintf(fp, "FLT0\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n",
			p.hdr.Samples[j], s.npass, s.nnonRef, s.nhomRR, s.nhomAA, s.nhet, s.nhemi,
			s.nSNV, s.nIndel, s.nsingleton, s.nmissing, s.nts, s.ntv, tstvStr(s.nts, s.ntv))
	}
	s := &p.siteStats
	fmt.Fprintf(fp, "SITE0\t%d\t%d\t%d\t%d\t%d\t%d\t%s\n",
		s.npass, s.nSNV, s.nIndel, s.nsingleton, s.nts, s.ntv, tstvStr(s.nts, s.ntv))
	return nil
}

// tstvStr formats the ts/tv ratio the way upstream does: "%.2f", or "inf" when
// there are transitions but no transversions (C's INFINITY printed by %.2f).
func tstvStr(nts, ntv int) string {
	if ntv == 0 {
		return "inf"
	}
	return fmt.Sprintf("%.2f", float64(nts)/float64(ntv))
}

const smplStatsHeader = `# CMD line shows the command line used to generate this output
# DEF lines define expressions for all tested thresholds
# FLT* lines report numbers for every threshold and every sample:
#   1) filter id
#   2) sample
#   3) number of genotypes which pass the filter
#   4) number of non-reference genotypes
#   5) number of homozygous ref genotypes (0/0 or 0)
#   6) number of homozygous alt genotypes (1/1, 2/2, etc)
#   7) number of heterozygous genotypes (0/1, 1/2, etc)
#   8) number of hemizygous genotypes (0, 1, etc)
#   9) number of SNVs
#   10) number of indels
#   11) number of singletons
#   12) number of missing genotypes (./., ., ./0, etc)
#   13) number of transitions (alt het genotypes such as "1/2" are counted twice)
#   14) number of transversions (alt het genotypes such as "1/2" are counted twice)
#   15) overall ts/tv
# SITE* lines report numbers for every threshold:
#   1) filter id
#   2) number of sites which pass the filter
#   3) number of SNVs
#   4) number of indels
#   5) number of singletons
#   6) number of transitions (counted at most once at multiallelic sites)
#   7) number of transversions (counted at most once at multiallelic sites)
#   8) overall ts/tv
`
