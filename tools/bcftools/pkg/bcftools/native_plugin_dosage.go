// Native port of the upstream `dosage` plugin (plugins/dosage.c). It prints
// genotype dosage determined from FORMAT tags requested by the user (-t,
// default PL,GL,GT), suppresses the VCF/BCF output, and emits a TSV with one
// row per record: CHROM, POS, REF, ALT, and a per-sample dosage column.
package bcftools

import (
	"fmt"
	"io"
	"math"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() { registerNativePlugin("dosage", func() NativePlugin { return &dosagePlugin{} }) }

// dosagePlugin implements the `dosage` plugin. It accumulates no state across
// records but must print its header once before the per-record rows and must
// suppress the VCF/BCF output, so it is driven serially via the suppress path.
type dosagePlugin struct {
	hdr      *vcf.Header
	handlers []string // ordered subset of {"PL","GL","GT"} that have a usable header
	out      io.Writer
}

// SuppressVCF reports true: dosage emits only its TSV on stdout.
func (p *dosagePlugin) SuppressVCF() bool { return true }

// SetStdout wires the host stdout writer the TSV is printed to.
func (p *dosagePlugin) SetStdout(w io.Writer) { p.out = w }

// Name returns the plugin name.
func (p *dosagePlugin) Name() string { return "dosage" }

// About returns the one-line description, matching dosage.c about().
func (p *dosagePlugin) About() string {
	return "Prints genotype dosage determined from tags requested by the user."
}

// Parallel reports false: the header is printed once and rows are emitted in
// input order, so records are processed serially.
func (p *dosagePlugin) Parallel() bool { return false }

// Init parses -t/--tags, validates each tag against the header, and prints the
// TSV header line. The requested tags are kept in order; PL/GL are used only
// when the FORMAT line exists, while GT is always usable.
func (p *dosagePlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	tagsStr := "PL,GL,GT"
	for i := 0; i < len(args); i++ {
		a := args[i]
		switch a {
		case "-t", "--tags":
			if i+1 >= len(args) {
				return nil, fmt.Errorf("dosage: -t requires an argument")
			}
			i++
			tagsStr = args[i]
		default:
			if strings.HasPrefix(a, "-t") && len(a) > 2 {
				tagsStr = a[2:]
				continue
			}
			return nil, fmt.Errorf("dosage: unsupported option %q", a)
		}
	}
	p.hdr = hdr
	for _, t := range strings.Split(tagsStr, ",") {
		switch t {
		case "PL", "GL":
			if hasFormatHeader(hdr.MetaInfo, t) {
				p.handlers = append(p.handlers, t)
			}
		case "GT":
			p.handlers = append(p.handlers, "GT")
		default:
			return nil, fmt.Errorf("dosage: no handler for tag %q", t)
		}
	}

	var b strings.Builder
	b.WriteString("#[1]CHROM\t[2]POS\t[3]REF\t[4]ALT")
	for i, s := range hdr.Samples {
		fmt.Fprintf(&b, "\t[%d]%s", i+5, s)
	}
	b.WriteByte('\n')
	if p.out != nil {
		_, _ = io.WriteString(p.out, b.String())
	}
	return hdr, nil
}

// Process prints the dosage row for one record and drops it from the output.
func (p *dosagePlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	var b strings.Builder
	b.WriteString(v.Chrom)
	b.WriteByte('\t')
	b.WriteString(strconv.Itoa(v.Pos))
	b.WriteByte('\t')
	b.WriteString(v.Ref)
	nals := 1 + len(v.Alt)
	if nals == 1 {
		b.WriteString("\t.")
	} else {
		b.WriteByte('\t')
		b.WriteString(strings.Join(v.Alt, ","))
	}

	if nals == 1 {
		for range v.Samples {
			b.WriteString("\t0.0")
		}
	} else {
		printed := false
		for _, h := range p.handlers {
			if p.emit(&b, v, h, nals) {
				printed = true
				break
			}
		}
		if !printed {
			for range v.Samples {
				b.WriteString("\t-1.0")
			}
		}
	}
	b.WriteByte('\n')
	if p.out != nil {
		_, _ = io.WriteString(p.out, b.String())
	}
	return nil, nil
}

// Destroy releases resources (none held).
func (p *dosagePlugin) Destroy() error { return nil }

// emit appends one handler's dosage columns for every sample and reports
// whether the handler successfully produced output (the tag was present and
// well-formed for the diploid expectation), mirroring the C handlers' return
// convention (0 == success).
func (p *dosagePlugin) emit(b *strings.Builder, v *vcf.Variant, tag string, nals int) bool {
	switch tag {
	case "PL", "GL":
		return p.emitLikelihood(b, v, tag, nals)
	case "GT":
		return p.emitGT(b, v, nals)
	}
	return false
}

// emitLikelihood handles PL/GL: it converts the per-sample likelihoods to
// linear scale, normalises them, then accumulates per-allele dosage over the
// diploid genotype index ordering, exactly as calc_dosage_PL/GL do. The whole
// record is rejected (no output) if any sample lacks the tag or the value count
// is not the diploid n*(n+1)/2.
func (p *dosagePlugin) emitLikelihood(b *strings.Builder, v *vcf.Variant, tag string, nals int) bool {
	if !formatHasTag(v, tag) {
		return false
	}
	ndip := nals * (nals + 1) / 2
	rows := make([][]float32, len(v.Samples))
	for i := range v.Samples {
		raw, ok := v.Samples[i].Data[tag]
		if !ok {
			return false
		}
		parts := strings.Split(raw, ",")
		if len(parts) != ndip {
			return false
		}
		// Upstream computes vals/sum/dsg in C `float` precision, so the
		// pow/sum/divide/accumulate all run in 32 bits; mirroring that here is
		// what makes the printed dosages match byte-for-byte.
		dsg := make([]float32, nals)
		vals := make([]float32, ndip)
		broke := false
		var sum float32
		for j, s := range parts {
			f, err := strconv.ParseFloat(s, 64)
			if s == "." || err != nil {
				broke = true
				break
			}
			if tag == "PL" {
				vals[j] = float32(math.Pow(10, -0.1*f))
			} else {
				vals[j] = float32(math.Pow(10, f))
			}
			sum += vals[j]
		}
		if broke {
			for j := range dsg {
				dsg[j] = -1
			}
			rows[i] = dsg
			continue
		}
		if sum != 0 {
			for j := range vals {
				vals[j] /= sum
			}
		}
		vals[0] = 0
		l := 0
		for j := 0; j < nals; j++ {
			for k := 0; k <= j; k++ {
				dsg[j] += vals[l]
				dsg[k] += vals[l]
				l++
			}
		}
		rows[i] = dsg
	}
	for _, dsg := range rows {
		writeDosageRow32(b, dsg, nals)
	}
	return true
}

// emitGT handles GT: it counts allele occurrences per sample, yielding integer
// dosages, or -1 across all alleles for a fully missing genotype. Unlike the
// likelihood handlers this never rejects the record (calc_dosage_GT returns 0
// whenever GT is present, which it is here since the records are textual VCF).
func (p *dosagePlugin) emitGT(b *strings.Builder, v *vcf.Variant, nals int) bool {
	for i := range v.Samples {
		dsg := make([]float64, nals)
		gt, ok := sampleGT(v, i)
		ncalled := 0
		if ok {
			for _, a := range gt.alleles {
				if a == missingAllele {
					break
				}
				ncalled++
				if a >= 0 && a < nals {
					dsg[a]++
				}
			}
		}
		if ncalled == 0 {
			for j := range dsg {
				dsg[j] = -1
			}
		}
		writeDosageRowGT(b, dsg, nals)
	}
	return true
}

// writeDosageRow32 appends the per-ALT likelihood dosages (alleles 1..nals-1)
// for one sample as C printf("%f", ...) would: six digits after the decimal
// point, tab-separating the first column and comma-separating the rest.
func writeDosageRow32(b *strings.Builder, dsg []float32, nals int) {
	for j := 1; j < nals; j++ {
		if j == 1 {
			b.WriteByte('\t')
		} else {
			b.WriteByte(',')
		}
		b.WriteString(strconv.FormatFloat(float64(dsg[j]), 'f', 6, 64))
	}
}

// writeDosageRowGT appends the per-ALT GT dosages (alleles 1..nals-1) for one
// sample as C printf("%.1f", ...) would: one digit after the decimal point.
func writeDosageRowGT(b *strings.Builder, dsg []float64, nals int) {
	for j := 1; j < nals; j++ {
		if j == 1 {
			b.WriteByte('\t')
		} else {
			b.WriteByte(',')
		}
		b.WriteString(strconv.FormatFloat(dsg[j], 'f', 1, 64))
	}
}
