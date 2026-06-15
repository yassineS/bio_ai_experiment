// Native port of the upstream `allele-length` plugin (plugins/allele-length.c).
// It accumulates the frequency distribution of REF, ALT (first ALT), and
// REF+ALT lengths across all records, suppresses the VCF/BCF output, and prints
// a fixed 512-row table plus a trailing totals line at the end.
package bcftools

import (
	"fmt"
	"io"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() {
	registerNativePlugin("allele-length", func() NativePlugin { return &alleleLengthPlugin{} })
}

// alleleLengthMaxLen mirrors the upstream MAXLEN: lengths at or above it are
// clamped into the last bucket.
const alleleLengthMaxLen = 512

// alleleLengthPlugin implements the `allele-length` plugin. It accumulates
// totals across records and emits its summary in Destroy, so it is not
// parallel and suppresses the VCF/BCF output.
type alleleLengthPlugin struct {
	numvar, numxvar                       uint64
	reflen, altlen, refaltlen, xrefaltlen [alleleLengthMaxLen]uint64
	out                                   io.Writer
}

// SuppressVCF reports true: allele-length emits only its summary table.
func (p *alleleLengthPlugin) SuppressVCF() bool { return true }

// SetStdout wires the host stdout writer the table is printed to.
func (p *alleleLengthPlugin) SetStdout(w io.Writer) { p.out = w }

// Name returns the plugin name.
func (p *alleleLengthPlugin) Name() string { return "allele-length" }

// About returns the one-line description, matching allele-length.c about().
func (p *alleleLengthPlugin) About() string {
	return "Count the frequency of the length of REF, ALT and REF+ALT"
}

// Parallel reports false: the per-length counters are accumulated serially.
func (p *alleleLengthPlugin) Parallel() bool { return false }

// Init takes no plugin options.
func (p *alleleLengthPlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	if len(args) > 0 {
		return nil, fmt.Errorf("allele-length: unexpected argument %q", args[0])
	}
	return hdr, nil
}

// Process accumulates the REF/ALT/REF+ALT length buckets for one record and
// drops it. The first ALT allele is used, matching upstream's d.allele[1].
func (p *alleleLengthPlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	alt := ""
	if len(v.Alt) > 0 {
		alt = v.Alt[0]
	}
	rl := len(v.Ref)
	al := len(alt)
	ral := rl + al
	if rl >= alleleLengthMaxLen {
		rl = alleleLengthMaxLen - 1
	}
	if al >= alleleLengthMaxLen {
		al = alleleLengthMaxLen - 1
	}
	if ral >= alleleLengthMaxLen {
		ral = alleleLengthMaxLen - 1
	}
	p.reflen[rl]++
	p.altlen[al]++
	p.refaltlen[ral]++
	if containNonBase(v.Ref) || containNonBase(alt) {
		p.xrefaltlen[ral]++
		p.numxvar++
	}
	p.numvar++
	return nil, nil
}

// Destroy prints the accumulated length distribution in upstream's exact
// layout: a header row, one row per length 0..511, and a trailing totals line.
func (p *alleleLengthPlugin) Destroy() error {
	if p.out == nil {
		return nil
	}
	fmt.Fprintf(p.out, "LENGTH\tREF\tALT\tREF+ALT\tREF+ALT WITH NON-BASE NUCLEOTIDES\n")
	for i := 0; i < alleleLengthMaxLen; i++ {
		fmt.Fprintf(p.out, "%d\t%d\t%d\t%d\t%d\n", i, p.reflen[i], p.altlen[i], p.refaltlen[i], p.xrefaltlen[i])
	}
	fmt.Fprintf(p.out, "\t\t\t%d\t%d\n", p.numvar, p.numxvar)
	return nil
}

// containNonBase reports whether s contains any character other than the
// standard ACGT bases (case-insensitive), mirroring contain_non_base. An empty
// string contains no non-base characters.
func containNonBase(s string) bool {
	for i := 0; i < len(s); i++ {
		switch s[i] {
		case 'A', 'a', 'C', 'c', 'G', 'g', 'T', 't':
		default:
			return true
		}
	}
	return false
}
