// Native port of the upstream `af-dist` plugin (plugins/af-dist.c). It collects
// allele-frequency deviation statistics and a genotype-probability distribution
// (assuming Hardy-Weinberg equilibrium) from an INFO allele-frequency tag and
// the per-sample genotypes, suppresses the VCF/BCF output, and prints its tables
// to stdout. Only non-reference genotypes with a known site allele frequency are
// considered; the probabilities are 2*AF*(1-AF) for the RA genotype and AF**2
// for AA, matching upstream exactly (including its 32-bit float arithmetic).
package bcftools

import (
	"fmt"
	"io"
	"math"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() { registerNativePlugin("af-dist", func() NativePlugin { return &afDistPlugin{} }) }

// afDistPlugin implements the `af-dist` plugin. It accumulates distribution
// counters across records and prints its tables once from Destroy, so it must
// run serially and suppresses the VCF/BCF output.
type afDistPlugin struct {
	hdr      *vcf.Header
	afTag    string
	devBins  *bins
	probBins *bins
	devDist  []uint64
	probDist []uint64

	listSet bool
	listMin float32
	listMax float32

	samples   bool
	nsmplProb []uint64
	smplProb  []float64

	out io.Writer
}

// SuppressVCF reports true: `+af-dist` emits no VCF/BCF output, only its
// textual tables on stdout (upstream init() returns 1).
func (p *afDistPlugin) SuppressVCF() bool { return true }

// SetStdout wires the host stdout writer the tables are printed to.
func (p *afDistPlugin) SetStdout(w io.Writer) { p.out = w }

// Name returns the plugin name.
func (p *afDistPlugin) Name() string { return "af-dist" }

// About returns the one-line description, matching af-dist.c about().
func (p *afDistPlugin) About() string { return "AF and GT probability distribution stats." }

// Parallel reports false: the counters are shared and the header tables are
// printed once, so records are processed serially.
func (p *afDistPlugin) Parallel() bool { return false }

// Init parses the plugin options (-t/--af-tag, -d/--dev-bins, -p/--prob-bins,
// -s/--samples, -l/--list), builds the two bin tables, and prints the leading
// comment header. It mirrors af-dist.c init(), including the printed command
// line and the optional per-genotype listing header.
func (p *afDistPlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	p.hdr = hdr
	p.afTag = "AF"
	p.listMin = -1
	devBinsDef := "0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1"
	probBinsDef := "0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1"

	for i := 0; i < len(args); i++ {
		a := args[i]
		next := func() (string, error) {
			if i+1 >= len(args) {
				return "", fmt.Errorf("af-dist: %s requires an argument", a)
			}
			i++
			return args[i], nil
		}
		switch a {
		case "-t", "--af-tag":
			v, err := next()
			if err != nil {
				return nil, err
			}
			p.afTag = v
		case "-d", "--dev-bins":
			v, err := next()
			if err != nil {
				return nil, err
			}
			devBinsDef = v
		case "-p", "--prob-bins":
			v, err := next()
			if err != nil {
				return nil, err
			}
			probBinsDef = v
		case "-s", "--samples":
			p.samples = true
		case "-l", "--list":
			v, err := next()
			if err != nil {
				return nil, err
			}
			comma := strings.IndexByte(v, ',')
			if comma < 0 {
				return nil, fmt.Errorf("af-dist: could not parse: --list %s", v)
			}
			lo, err1 := strconv.ParseFloat(v[:comma], 64)
			hi, err2 := strconv.ParseFloat(v[comma+1:], 64)
			if err1 != nil || err2 != nil {
				return nil, fmt.Errorf("af-dist: could not parse: --list %s", v)
			}
			p.listSet = true
			p.listMin = float32(lo)
			p.listMax = float32(hi)
		default:
			return nil, fmt.Errorf("af-dist: unsupported option %q", a)
		}
	}

	var err error
	if p.devBins, err = binInit(devBinsDef, 0, 1); err != nil {
		return nil, err
	}
	p.devDist = make([]uint64, p.devBins.size())
	if p.probBins, err = binInit(probBinsDef, 0, 1); err != nil {
		return nil, err
	}
	p.probDist = make([]uint64, p.probBins.size())

	if p.samples {
		p.nsmplProb = make([]uint64, len(hdr.Samples))
		p.smplProb = make([]float64, len(hdr.Samples))
	}

	if p.out != nil {
		// The version banner is provenance and is stripped by the oracle; the
		// command-line echo is reproduced for completeness, matching the layout
		// (it too is part of the stripped leading comment block in practice).
		fmt.Fprintf(p.out, "# This file was produced by: bcftools +af-dist\n")
		fmt.Fprintf(p.out, "# The command line was:\tbcftools +af-dist %s\n#\n", strings.Join(args, " "))
		if p.listSet {
			fmt.Fprintf(p.out, "# GT, genotypes with P(AF) in [%f,%f]; [2]Chromosome\t[3]Position[4]Sample\t[5]Genotype\t[6]AF-based probability\n",
				p.listMin, p.listMax)
		}
	}
	return hdr, nil
}

// Process accumulates the distribution counters for one record and drops it.
// It ports af-dist.c process(): it reads INFO/<af-tag>, computes the HWE
// genotype probabilities, bins the RA/AA probabilities, optionally lists the
// matching genotypes, and tallies the AF deviation. Records without the AF tag
// are skipped, as upstream does (process returns NULL).
func (p *afDistPlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	afStr, ok := v.Info[p.afTag]
	if !ok || afStr == "" {
		return nil, nil
	}
	// Only the first AF value is used (args->af[0]).
	if comma := strings.IndexByte(afStr, ','); comma >= 0 {
		afStr = afStr[:comma]
	}
	af64, err := strconv.ParseFloat(afStr, 64)
	if err != nil {
		return nil, nil
	}
	af := float32(af64)

	pRR := (1 - af) * (1 - af)
	pRA := 2 * af * (1 - af)
	pAA := af * af
	iRA := p.probBins.idx(pRA)
	iAA := p.probBins.idx(pAA)
	lRR := math.Log(float64(pRR))
	lRA := math.Log(float64(pRA))
	lAA := math.Log(float64(pAA))

	listRA := p.listSet && !(pRA < p.listMin || pRA > p.listMax)
	listAA := p.listSet && !(pAA < p.listMin || pAA > p.listMax)

	nals, nalt := 0, 0
	for i := range v.Samples {
		gt, ok := sampleGT(v, i)
		if !ok {
			continue
		}
		dosage := 0
		full := true
		j := 0
		for ; j < len(gt.alleles); j++ {
			a := gt.alleles[j]
			if a == missingAllele {
				full = false
				break
			}
			if a == 1 {
				dosage++
			}
		}
		if !full {
			continue
		}

		nals += j
		nalt += dosage

		lprob := lRR
		if dosage == 1 {
			lprob = lRA
			p.probDist[iRA]++
			if listRA && p.out != nil {
				fmt.Fprintf(p.out, "GT\t%s\t%d\t%s\t1\t%f\n", v.Chrom, v.Pos, p.hdr.Samples[i], pRA)
			}
		} else if dosage == 2 {
			lprob = lAA
			p.probDist[iAA]++
			if listAA && p.out != nil {
				fmt.Fprintf(p.out, "GT\t%s\t%d\t%s\t2\t%f\n", v.Chrom, v.Pos, p.hdr.Samples[i], pAA)
			}
		}
		if p.samples {
			// NB: upstream stores (not accumulates) lprob here, so only the last
			// record's value survives; this is reproduced faithfully.
			p.smplProb[i] = lprob
			p.nsmplProb[i]++
		}
	}

	if nals != 0 && (nalt != 0 || af != 0) {
		afDev := abs32(af - float32(nalt)/float32(nals))
		iAF := p.devBins.idx(afDev)
		p.devDist[iAF]++
	}
	return nil, nil
}

// Destroy prints the probability distribution, the AF-deviation distribution,
// and (with -s) the per-sample HWE log probability, in upstream's exact layout.
func (p *afDistPlugin) Destroy() error {
	w := p.out
	if w == nil {
		return nil
	}
	fmt.Fprintf(w, "# PROB_DIST, genotype probability distribution, assumes HWE\n")
	n := p.probBins.size()
	for i := 0; i < n-1; i++ {
		fmt.Fprintf(w, "PROB_DIST\t%f\t%f\t%d\n", p.probBins.value(i), p.probBins.value(i+1), p.probDist[i])
	}
	fmt.Fprintf(w, "# DEV_DIST, distribution of AF deviation, based on %s and INFO/AN, AC calculated on the fly\n", p.afTag)
	n = p.devBins.size()
	for i := 0; i < n-1; i++ {
		fmt.Fprintf(w, "DEV_DIST\t%f\t%f\t%d\n", p.devBins.value(i), p.devBins.value(i+1), p.devDist[i])
	}
	if p.samples {
		fmt.Fprintf(w, "# SMPL_PROB, per-sample HWE log probability (geometric mean) and the number of genotypes\n")
		for i, name := range p.hdr.Samples {
			val := 0.0
			if p.nsmplProb[i] != 0 {
				val = p.smplProb[i] / float64(p.nsmplProb[i])
			}
			fmt.Fprintf(w, "SMPL_PROB\t%s\t%s\t%d\n", name, formatCDoubleE(val), p.nsmplProb[i])
		}
	}
	return nil
}

// formatCDoubleE renders v the way C's printf("%e", ...) does: a mantissa with
// six fractional digits and a two-or-more-digit signed exponent (e.g.
// "-1.203973e+00"). Go's 'e' verb already matches this, including the minimum
// two-digit exponent.
func formatCDoubleE(v float64) string {
	return strconv.FormatFloat(v, 'e', 6, 64)
}
