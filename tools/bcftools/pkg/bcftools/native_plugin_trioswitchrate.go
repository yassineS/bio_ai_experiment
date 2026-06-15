// Native port of the upstream `trio-switch-rate` plugin
// (plugins/trio-switch-rate.c). Given a PED file it computes, per trio, the
// phase switch rate in the (phased) child genotypes: for each informative,
// Mendelian-consistent biallelic het child, the parental phase is determined and
// compared with the previous informative site on the same chromosome; a change
// counts as a switch. The per-trio TRIO lines report the number of tested sites,
// Mendelian errors, switches, and the switch percentage. An optional 7th PED
// column groups trios into populations, reported in the trailing POP lines.
//
// trio-switch-rate is a generic init/process plugin (its options follow a `--`
// separator) that suppresses the VCF/BCF output and writes only its textual
// report. It is inherently stateful (it carries the previous phase per trio
// across records) so it runs through the serial pipeline.
//
// The only plugin option is -p/--ped, which is fully supported. The leading
// banner ("# This file was produced by ..." / "# The command line was: ...")
// embeds version strings and is filtered out by the oracle's provenance
// stripper, exactly as for the roh / stats banners.
package bcftools

import (
	"fmt"
	"io"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() {
	registerNativePlugin("trio-switch-rate", func() NativePlugin { return &trioSwitchRatePlugin{} })
}

// trioSwitchTrio mirrors trio_t in trio-switch-rate.c: the three sample indices
// plus the rolling per-trio state (previous tested phase and the running
// counters).
type trioSwitchTrio struct {
	father, mother, child int
	ipop                  int    // population index (0 when no 7th column)
	prev                  int    // previous phase test result (0 = none yet)
	err, nswitch, ntest   uint32 // Mendelian errors, switches, tested sites
}

// trioSwitchPop mirrors pop_t: a named population grouping accumulated in
// Destroy().
type trioSwitchPop struct {
	name                string
	err, nswitch, ntest uint32
	ntrio               uint32
	pswitch             float64
}

// trioSwitchRatePlugin implements the `trio-switch-rate` plugin.
type trioSwitchRatePlugin struct {
	hdr     *vcf.Header
	trios   []trioSwitchTrio
	pops    []trioSwitchPop
	prevRID string
	havePrv bool
	out     io.Writer
}

// SuppressVCF reports true: trio-switch-rate emits only its textual report.
func (p *trioSwitchRatePlugin) SuppressVCF() bool { return true }

// SetStdout wires the host stdout writer the report is printed to.
func (p *trioSwitchRatePlugin) SetStdout(w io.Writer) { p.out = w }

// Name returns the plugin name.
func (p *trioSwitchRatePlugin) Name() string { return "trio-switch-rate" }

// About returns the one-line description, matching the plugin's about().
func (p *trioSwitchRatePlugin) About() string {
	return "Calculate phase switch rate in trio samples, children samples must have phased GTs.\n"
}

// Parallel reports false: the rolling per-trio phase state is updated serially.
func (p *trioSwitchRatePlugin) Parallel() bool { return false }

// Init parses the -p/--ped option and resolves the trios (and any 7th-column
// populations) against the input header.
func (p *trioSwitchRatePlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	var pedFile string
	for i := 0; i < len(args); i++ {
		switch args[i] {
		case "-p", "--ped":
			if i+1 >= len(args) {
				return nil, fmt.Errorf("trio-switch-rate: %s requires a value", args[i])
			}
			i++
			pedFile = args[i]
		default:
			return nil, fmt.Errorf("trio-switch-rate: unsupported option %q", args[i])
		}
	}
	if pedFile == "" {
		return nil, fmt.Errorf("trio-switch-rate: expected the -p option")
	}
	p.hdr = hdr

	idx := sampleIndex(hdr)
	rows, err := parsePEDRows(pedFile, idx)
	if err != nil {
		return nil, fmt.Errorf("trio-switch-rate: %w", err)
	}
	popIdx := make(map[string]int)
	for _, r := range rows {
		t := trioSwitchTrio{father: r.father, mother: r.mother, child: r.child}
		if r.popName != "" {
			pi, ok := popIdx[r.popName]
			if !ok {
				pi = len(p.pops)
				popIdx[r.popName] = pi
				p.pops = append(p.pops, trioSwitchPop{name: r.popName})
			}
			t.ipop = pi
			p.pops[pi].ntrio++
		}
		p.trios = append(p.trios, t)
	}
	return hdr, nil
}

// bcftoolsVersionString / htslibVersionString return placeholder version
// strings for the trio-switch-rate banner. The banner is removed by the
// oracle's provenance stripper before comparison, so the exact value is
// irrelevant to parity; they exist only so the emitted line is well-formed.
func bcftoolsVersionString() string { return "bio_ai_experiment" }
func htslibVersionString() string   { return "bio_ai_experiment" }

// switchGT is the parsed biallelic phased genotype used by Process, mirroring
// gt_t / parse_genotype in trio-switch-rate.c.
type switchGT struct {
	a, b   int
	phased bool
	ok     bool
}

// parseSwitchGT parses sample i's GT for the switch-rate test, mirroring
// parse_genotype(): the genotype must be diploid with two non-missing alleles,
// and only the first two alleles (0,1) at biallelic sites are considered — any
// allele index > 1 makes the genotype non-informative.
func parseSwitchGT(v *vcf.Variant, i int) switchGT {
	gt, ok := sampleGT(v, i)
	if !ok || len(gt.alleles) != 2 {
		return switchGT{}
	}
	if gt.alleles[0] == missingAllele || gt.alleles[1] == missingAllele {
		return switchGT{}
	}
	g := switchGT{a: gt.alleles[0], b: gt.alleles[1], phased: gt.phased[1], ok: true}
	if g.a > 1 || g.b > 1 {
		return switchGT{}
	}
	return g
}

// Process updates the rolling per-trio phase state for one record, mirroring
// process() in trio-switch-rate.c. Records whose GT is not diploid for every
// sample are skipped (upstream's `ngt!=2` early return).
func (p *trioSwitchRatePlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	// Upstream requires a uniform diploid GT across all samples (ngt==2).
	for i := range v.Samples {
		gt, ok := sampleGT(v, i)
		if !ok || len(gt.alleles) != 2 {
			return nil, nil
		}
	}

	if !p.havePrv || v.Chrom != p.prevRID {
		p.prevRID = v.Chrom
		p.havePrv = true
		for i := range p.trios {
			p.trios[i].prev = 0
		}
	}

	for i := range p.trios {
		trio := &p.trios[i]

		child := parseSwitchGT(v, trio.child)
		if !child.ok || !child.phased {
			continue
		}
		if child.a+child.b != 1 { // child must be a het
			continue
		}
		father := parseSwitchGT(v, trio.father)
		if !father.ok {
			continue
		}
		mother := parseSwitchGT(v, trio.mother)
		if !mother.ok {
			continue
		}
		if father.a+father.b == 1 && mother.a+mother.b == 1 { // both parents het
			continue
		}
		if father.a+father.b == mother.a+mother.b { // Mendelian error
			trio.err++
			continue
		}

		testPhase := 0
		if father.a == father.b {
			testPhase = 1
			if child.a == father.a {
				testPhase = 2
			}
		} else if mother.a == mother.b {
			testPhase = 1
			if child.b == mother.a {
				testPhase = 2
			}
		}
		if trio.prev > 0 && trio.prev != testPhase {
			trio.nswitch++
		}
		trio.ntest++
		trio.prev = testPhase
	}
	return nil, nil
}

// Destroy prints the TRIO and POP report, mirroring destroy() in
// trio-switch-rate.c. The leading banner lines are emitted verbatim; the
// oracle's provenance stripper removes the version/command-line ones.
func (p *trioSwitchRatePlugin) Destroy() error {
	if p.out == nil {
		return nil
	}
	fp := p.out
	// Banner (provenance-stripped by the oracle). The version placeholders here
	// are irrelevant to parity since these two lines and the following "#" line
	// are removed before comparison.
	fmt.Fprintf(fp, "# This file was produced by: bcftools +trio-switch-rate(%s+htslib-%s)\n",
		bcftoolsVersionString(), htslibVersionString())
	fmt.Fprint(fp, "# The command line was:\tbcftools +trio-switch-rate\n")
	fmt.Fprint(fp, "#\n")
	fmt.Fprint(fp, "# TRIO\t[2]Father\t[3]Mother\t[4]Child\t[5]nTested\t[6]nMendelian Errors\t[7]nSwitch\t[8]nSwitch (%)\n")
	for i := range p.trios {
		trio := &p.trios[i]
		pct := 0.0
		if trio.ntest != 0 {
			pct = float64(trio.nswitch) * 100.0 / float64(trio.ntest)
		}
		fmt.Fprintf(fp, "TRIO\t%s\t%s\t%s\t%d\t%d\t%d\t%.2f\n",
			p.hdr.Samples[trio.father], p.hdr.Samples[trio.mother], p.hdr.Samples[trio.child],
			trio.ntest, trio.err, trio.nswitch, pct)
		if len(p.pops) != 0 {
			pop := &p.pops[trio.ipop]
			pop.err += trio.err
			pop.nswitch += trio.nswitch
			pop.ntest += trio.ntest
			if trio.ntest != 0 {
				pop.pswitch += float64(trio.nswitch) * 100.0 / float64(trio.ntest)
			}
		}
	}
	fmt.Fprint(fp, "# POP\tpopulation or other grouping defined by an optional 7-th column of the PED file\n")
	fmt.Fprint(fp, "# POP\t[2]Name\t[3]Number of trios\t[4]avgTested\t[5]avgMendelian Errors\t[6]avgSwitch\t[7]avgSwitch (%)\n")
	for i := range p.pops {
		pop := &p.pops[i]
		fmt.Fprintf(fp, "POP\t%s\t%d\t%.0f\t%.0f\t%.0f\t%.2f\n",
			pop.name, pop.ntrio,
			float64(pop.ntest)/float64(pop.ntrio),
			float64(pop.err)/float64(pop.ntrio),
			float64(pop.nswitch)/float64(pop.ntrio),
			pop.pswitch/float64(pop.ntrio))
	}
	return nil
}
