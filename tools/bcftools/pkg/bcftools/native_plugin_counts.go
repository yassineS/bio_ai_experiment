// Native port of the upstream `counts` plugin (plugins/counts.c). It counts
// the number of samples, SNPs, indels, MNPs, "others", and total sites,
// suppresses all record output, and prints the totals to stdout at the end.
package bcftools

import (
	"fmt"
	"io"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() { registerNativePlugin("counts", func() NativePlugin { return &countsPlugin{} }) }

// countsPlugin implements the `counts` plugin. Because Destroy must report
// running totals, it accumulates state across records and so is NOT parallel.
type countsPlugin struct {
	nsamples              int
	nsnps, nindels, nmnps int
	nothers, nsites       int
	out                   io.Writer // host stdout; counts reports here
}

// SuppressVCF reports true: `+counts` emits no VCF/BCF output, only its
// textual summary on stdout (upstream init() returns 1).
func (p *countsPlugin) SuppressVCF() bool { return true }

// SetStdout wires the host stdout writer the summary is printed to.
func (p *countsPlugin) SetStdout(w io.Writer) { p.out = w }

// Name returns the plugin name.
func (p *countsPlugin) Name() string { return "counts" }

// About returns the one-line description, matching counts.c about().
func (p *countsPlugin) About() string {
	return "A minimal plugin which counts number of samples, SNPs, INDELs, MNPs and total number of sites."
}

// Parallel reports false: counts accumulates totals and must see records
// serially. (Order does not affect the totals, but the shared counters are
// not concurrency-safe.)
func (p *countsPlugin) Parallel() bool { return false }

// Init records the sample count. Upstream's init() returns 1 to suppress the
// VCF/BCF output entirely; the native pipeline mirrors this via SuppressVCF
// (no header, no records), and counts emits only its textual summary on stdout
// from Destroy — matching the observable upstream behaviour of `+counts`.
func (p *countsPlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	if len(args) > 0 {
		return nil, fmt.Errorf("counts: unexpected argument %q", args[0])
	}
	p.nsamples = len(hdr.Samples)
	return hdr, nil
}

// Process classifies the variant by type and drops it (returns no records).
func (p *countsPlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	t := variantTypeMask(v)
	if t&vtSNP != 0 {
		p.nsnps++
	}
	if t&vtINDEL != 0 {
		p.nindels++
	}
	if t&vtMNP != 0 {
		p.nmnps++
	}
	if t&vtOTHER != 0 {
		p.nothers++
	}
	p.nsites++
	return nil, nil
}

// Destroy prints the accumulated totals to stdout in upstream's exact layout.
func (p *countsPlugin) Destroy() error {
	w := p.out
	if w == nil {
		return nil
	}
	fmt.Fprintf(w, "Number of samples: %d\n", p.nsamples)
	fmt.Fprintf(w, "Number of SNPs:    %d\n", p.nsnps)
	fmt.Fprintf(w, "Number of INDELs:  %d\n", p.nindels)
	fmt.Fprintf(w, "Number of MNPs:    %d\n", p.nmnps)
	fmt.Fprintf(w, "Number of others:  %d\n", p.nothers)
	fmt.Fprintf(w, "Number of sites:   %d\n", p.nsites)
	return nil
}
