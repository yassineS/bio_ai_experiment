// Native port of the upstream `check-ploidy` plugin (plugins/check-ploidy.c).
// It reports, per sample, the contiguous chromosomal regions over which a
// consistent genotype ploidy is observed. The VCF/BCF output is suppressed; a
// tab-separated report is printed to stdout (a header from Init, one line per
// region emitted as runs close out, and the trailing runs at Destroy).
package bcftools

import (
	"fmt"
	"io"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() {
	registerNativePlugin("check-ploidy", func() NativePlugin { return &checkPloidyPlugin{} })
}

// checkPloidyDat tracks one sample's current ploidy run: the begin/end position
// (0-based, as upstream stores rec->pos) and the run's ploidy (0 == no run).
type checkPloidyDat struct {
	sample   string
	beg, end int
	ploidy   int
}

// checkPloidyPlugin implements the `check-ploidy` plugin. It carries per-sample
// run state across records and so is not parallel; it suppresses VCF output.
type checkPloidyPlugin struct {
	dat           []checkPloidyDat
	curChrom      string // chromosome of the current run (rid)
	haveChrom     bool   // whether curChrom has been set
	ignoreMissing bool
	out           io.Writer
}

// SuppressVCF reports true: check-ploidy emits only its textual report.
func (p *checkPloidyPlugin) SuppressVCF() bool { return true }

// SetStdout wires the host stdout writer the report is printed to.
func (p *checkPloidyPlugin) SetStdout(w io.Writer) { p.out = w }

// Name returns the plugin name.
func (p *checkPloidyPlugin) Name() string { return "check-ploidy" }

// About returns the one-line description, matching check-ploidy.c about().
func (p *checkPloidyPlugin) About() string {
	return "Check if ploidy of samples is consistent for all sites"
}

// Parallel reports false: the per-sample run state is updated serially.
func (p *checkPloidyPlugin) Parallel() bool { return false }

// Init parses -m/--use-missing, validates that GT is present in the header, and
// prints the report header line.
func (p *checkPloidyPlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	p.ignoreMissing = true
	for _, a := range args {
		switch a {
		case "-m", "--use-missing":
			p.ignoreMissing = false
		default:
			return nil, fmt.Errorf("check-ploidy: unsupported option %q", a)
		}
	}
	if !hasFormatHeader(hdr.MetaInfo, "GT") {
		return nil, fmt.Errorf("check-ploidy: GT field is not present")
	}
	p.dat = make([]checkPloidyDat, len(hdr.Samples))
	for i, s := range hdr.Samples {
		p.dat[i].sample = s
	}
	if p.out != nil {
		fmt.Fprint(p.out, "# [1]Sample\t[2]Chromosome\t[3]Region Start\t[4]Region End\t[5]Ploidy\n")
	}
	return hdr, nil
}

// Process updates the per-sample ploidy runs for one record and drops it. When
// the chromosome changes, the open runs from the previous chromosome are
// flushed before the current record's per-sample updates.
func (p *checkPloidyPlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	if !formatHasTag(v, "GT") {
		return nil, nil
	}

	// On a chromosome change, flush every open run under the previous
	// chromosome name and reset, matching the rid-change block upstream.
	if p.haveChrom && p.curChrom != v.Chrom {
		for i := range p.dat {
			dat := &p.dat[i]
			if dat.ploidy != 0 {
				p.printRun(dat, p.curChrom)
			}
			dat.ploidy = 0
		}
	}
	p.curChrom = v.Chrom
	p.haveChrom = true

	for i := range v.Samples {
		nal, missing := genotypePloidy(v, i, p.ignoreMissing)
		if nal == 0 || missing {
			continue // missing genotype
		}
		dat := &p.dat[i]
		if dat.ploidy == nal {
			dat.end = v.Pos
			continue
		}
		if dat.ploidy != 0 {
			p.printRun(dat, v.Chrom)
		}
		dat.ploidy = nal
		dat.beg = v.Pos
		dat.end = v.Pos
	}
	return nil, nil
}

// Destroy flushes any open per-sample runs under the last chromosome seen.
func (p *checkPloidyPlugin) Destroy() error {
	if p.out == nil {
		return nil
	}
	for i := range p.dat {
		dat := &p.dat[i]
		if dat.ploidy != 0 {
			p.printRun(dat, p.curChrom)
		}
		dat.ploidy = 0
	}
	return nil
}

// printRun emits one report line. Positions stored 0-based are printed 1-based,
// matching upstream's beg+1/end+1.
func (p *checkPloidyPlugin) printRun(dat *checkPloidyDat, chrom string) {
	fmt.Fprintf(p.out, "%s\t%s\t%d\t%d\t%d\n", dat.sample, chrom, dat.beg, dat.end, dat.ploidy)
}

// genotypePloidy returns (nal, missing) for sample i: nal is the number of
// alleles before the genotype's vector end (its ploidy), and missing reports a
// missing allele encountered while ignoreMissing is set. It mirrors the
// BRANCH_INT loop: scanning stops at the first vector end, and if ignoreMissing
// a missing allele short-circuits the sample as missing. The textual GT model
// stores exactly the called alleles, so there is no trailing vector end to
// account for; ploidy is the number of parsed alleles.
func genotypePloidy(v *vcf.Variant, i int, ignoreMissing bool) (nal int, missing bool) {
	gt, ok := sampleGT(v, i)
	if !ok {
		return 0, false
	}
	for _, a := range gt.alleles {
		if a == missingAllele && ignoreMissing {
			return nal, true
		}
		nal++
	}
	return nal, false
}
