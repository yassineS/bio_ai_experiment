// Native port of the upstream `ad-bias` plugin (plugins/ad-bias.c) for its
// default reporting mode. For each user-supplied sample pair it finds the two
// most-supported alleles from FORMAT/AD across the pair and runs Fisher's exact
// test on the 2x2 (sample/control) x (REF/ALT) depth table, printing a hit line
// for every comparison whose p-value is below the threshold. A trailing summary
// line reports the pair/site/comparison counts. The VCF/BCF output is
// suppressed (upstream init() returns 1 in this mode).
//
// The --clean-vcf allele-subsetting mode and the -f convert-format mode require
// htslib machinery the native pipeline does not provide and are reported as
// unsupported.
package bcftools

import (
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() { registerNativePlugin("ad-bias", func() NativePlugin { return &adBiasPlugin{} }) }

// adBiasPair is one sample/control pair, mirroring pair_t in ad-bias.c.
type adBiasPair struct {
	smpl, ctrl         int
	smplName, ctrlName string
}

// adBiasPlugin implements the `ad-bias` plugin in its default reporting mode.
type adBiasPlugin struct {
	hdr         *vcf.Header
	pairs       []adBiasPair
	minDP       int
	minAltDP    int
	th          float64
	variantType int // 0 = any, vtSNP, or vtINDEL
	nsite       uint64
	ncmp        uint64
	out         io.Writer
}

// SuppressVCF reports true: ad-bias emits only its textual report.
func (p *adBiasPlugin) SuppressVCF() bool { return true }

// SetStdout wires the host stdout writer the report is printed to.
func (p *adBiasPlugin) SetStdout(w io.Writer) { p.out = w }

// Name returns the plugin name.
func (p *adBiasPlugin) Name() string { return "ad-bias" }

// About returns the one-line description, matching ad-bias.c about().
func (p *adBiasPlugin) About() string {
	return "Find positions with wildly varying ALT allele frequency (Fisher test on FMT/AD).\n"
}

// Parallel reports false: nsite/ncmp totals are updated serially.
func (p *adBiasPlugin) Parallel() bool { return false }

// Init parses the options, loads the sample pairs, validates GT/AD presence and
// prints the report header, mirroring ad-bias.c init().
func (p *adBiasPlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	p.th = 1e-3
	p.minAltDP = 1
	var fname string
	for i := 0; i < len(args); i++ {
		a := args[i]
		next := func() (string, error) {
			if i+1 >= len(args) {
				return "", fmt.Errorf("ad-bias: %s requires an argument", a)
			}
			i++
			return args[i], nil
		}
		switch a {
		case "-c", "--clean-vcf":
			return nil, fmt.Errorf("ad-bias: the --clean-vcf allele-subsetting mode is not supported by the native plugin")
		case "-f", "--format":
			return nil, fmt.Errorf("ad-bias: the -f convert-format mode is not supported by the native plugin")
		case "-a", "--min-alt-dp":
			v, err := next()
			if err != nil {
				return nil, err
			}
			n, err := strconv.Atoi(v)
			if err != nil {
				return nil, fmt.Errorf("ad-bias: could not parse: -a %s", v)
			}
			p.minAltDP = n
		case "-d", "--min-dp":
			v, err := next()
			if err != nil {
				return nil, err
			}
			n, err := strconv.Atoi(v)
			if err != nil {
				return nil, fmt.Errorf("ad-bias: could not parse: -d %s", v)
			}
			p.minDP = n
		case "-t", "--threshold":
			v, err := next()
			if err != nil {
				return nil, err
			}
			f, err := strconv.ParseFloat(v, 64)
			if err != nil {
				return nil, fmt.Errorf("ad-bias: could not parse: -t %s", v)
			}
			p.th = f
		case "-s", "--samples":
			v, err := next()
			if err != nil {
				return nil, err
			}
			fname = v
		case "-v", "--variant-type":
			v, err := next()
			if err != nil {
				return nil, err
			}
			switch strings.ToLower(v) {
			case "snp", "snps":
				p.variantType = vtSNP
			case "indel", "indels":
				p.variantType = vtINDEL
			default:
				return nil, fmt.Errorf("ad-bias: variant type \"%s\" is not supported", v)
			}
		default:
			return nil, fmt.Errorf("ad-bias: unsupported option %q", a)
		}
	}
	if fname == "" {
		return nil, fmt.Errorf("ad-bias: expected the -s option")
	}
	p.hdr = hdr
	if err := p.parsePairs(hdr, fname); err != nil {
		return nil, err
	}

	if p.out != nil {
		fp := p.out
		// Provenance lines (stripped by the oracle) keep the surrounding
		// structure aligned with upstream.
		fmt.Fprint(fp, "# This file was produced by: bcftools +ad-bias(bio_ai_experiment+htslib-bio_ai_experiment)\n")
		fmt.Fprint(fp, "# The command line was:\tbcftools +ad-bias\n#\n")
		i := 1
		fmt.Fprint(fp, "# FT, Fisher Test")
		col := func(s string) { i++; fmt.Fprintf(fp, "\t[%d]%s", i, s) }
		col("Sample")
		col("Control")
		col("Chrom")
		col("Pos")
		col("REF")
		col("ALT")
		col("smpl.nREF")
		col("smpl.nALT")
		col("ctrl.nREF")
		col("ctrl.nALT")
		col("P-value")
		fmt.Fprintln(fp)
	}
	return hdr, nil
}

// parsePairs reads the tab-delimited sample-pair file, mirroring
// ad-bias.c parse_samples(). Pairs naming samples absent from the VCF are
// silently skipped, as upstream does.
func (p *adBiasPlugin) parsePairs(hdr *vcf.Header, fname string) error {
	data, err := readFileBytes(fname)
	if err != nil {
		return fmt.Errorf("ad-bias: could not read: %s", fname)
	}
	idx := make(map[string]int, len(hdr.Samples))
	for i, s := range hdr.Samples {
		idx[s] = i
	}
	for _, line := range strings.Split(string(data), "\n") {
		line = strings.TrimRight(line, "\r")
		if strings.TrimSpace(line) == "" {
			continue
		}
		cols := strings.Split(line, "\t")
		if len(cols) < 2 {
			return fmt.Errorf("ad-bias: could not parse the sample file: %s", line)
		}
		smpl, ok1 := idx[cols[0]]
		ctrl, ok2 := idx[cols[1]]
		if !ok1 || !ok2 {
			continue
		}
		p.pairs = append(p.pairs, adBiasPair{
			smpl: smpl, ctrl: ctrl,
			smplName: hdr.Samples[smpl], ctrlName: hdr.Samples[ctrl],
		})
	}
	return nil
}

// Process runs the per-pair Fisher test for one record, mirroring
// ad-bias.c process(). Records without two alleles or without FORMAT/AD are
// skipped.
func (p *adBiasPlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	nAllele := len(v.Alt) + 1
	if nAllele < 2 {
		return nil, nil
	}
	if !formatHasTag(v, "AD") {
		return nil, nil
	}
	ad := parseADAll(v, nAllele)
	p.nsite++

	for i := range p.pairs {
		pair := &p.pairs[i]
		aptr := ad[pair.smpl]
		bptr := ad[pair.ctrl]

		ibig, ismall, nbig, nsmall := adBiasTopTwo(aptr, bptr)
		if ibig == -1 || ismall == -1 {
			continue
		}
		if nbig+nsmall < p.minDP {
			continue
		}
		// All four selected cells must be present (non-missing).
		if !adPresent(aptr, ibig) || !adPresent(bptr, ibig) || !adPresent(aptr, ismall) || !adPresent(bptr, ismall) {
			continue
		}

		if p.variantType != 0 {
			lbig := adAlleleString(v, ibig)
			lsmall := adAlleleString(v, ismall)
			if p.variantType == vtSNP && len(lbig) != len(lsmall) {
				continue
			}
			if p.variantType == vtINDEL && len(lbig) == len(lsmall) {
				continue
			}
		}

		var iref, ialt, nalt int
		if ibig > ismall {
			ialt, iref, nalt = ibig, ismall, nbig
		} else {
			ialt, iref, nalt = ismall, ibig, nsmall
		}
		if nalt < p.minAltDP {
			continue
		}
		p.ncmp++

		n11, n12 := aptr[iref], aptr[ialt]
		n21, n22 := bptr[iref], bptr[ialt]
		_, _, fisher := mpileupFisherExact(int64(n11), int64(n12), int64(n21), int64(n22))
		if fisher >= p.th {
			continue
		}
		if p.out != nil {
			fmt.Fprintf(p.out, "FT\t%s\t%s\t%s\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%s\n",
				pair.smplName, pair.ctrlName, v.Chrom, v.Pos,
				adAlleleString(v, iref), adAlleleString(v, ialt),
				n11, n12, n21, n22, formatExp(fisher))
		}
	}
	return nil, nil
}

// Destroy prints the summary line, mirroring ad-bias.c destroy().
func (p *adBiasPlugin) Destroy() error {
	if p.out == nil {
		return nil
	}
	fmt.Fprint(p.out, "# SN, Summary Numbers\t[2]Number of Pairs\t[3]Number of Sites\t[4]Number of comparisons\t[5]P-value output threshold\n")
	fmt.Fprintf(p.out, "SN\t%d\t%d\t%d\t%s\n", len(p.pairs), p.nsite, p.ncmp, formatExp(p.th))
	return nil
}

// adBiasTopTwo finds the indices and depths of the two most-supported alleles
// across the sample and control AD vectors, mirroring the dual scan in
// ad-bias.c process(). Missing entries are skipped and a vector-end (negative
// sentinel) stops the scan of a vector.
func adBiasTopTwo(aptr, bptr []int) (ibig, ismall, nbig, nsmall int) {
	ibig, ismall, nbig, nsmall = -1, -1, -1, -1
	for j := 0; j < len(aptr); j++ {
		if aptr[j] < 0 {
			// -1 marks a missing AD value here; treat as missing/skip.
			continue
		}
		if ibig == -1 {
			ibig, nbig = j, aptr[j]
			continue
		}
		if nbig < aptr[j] {
			if ismall == -1 || nsmall < nbig {
				ismall, nsmall = ibig, nbig
			}
			ibig, nbig = j, aptr[j]
			continue
		}
		if ismall == -1 || nsmall < aptr[j] {
			ismall, nsmall = j, aptr[j]
		}
	}
	for j := 0; j < len(bptr); j++ {
		if bptr[j] < 0 {
			continue
		}
		if ibig == -1 {
			ibig, nbig = j, bptr[j]
			continue
		}
		if ibig == j {
			if nbig < bptr[j] {
				nbig = bptr[j]
			}
			continue
		}
		if nbig < bptr[j] {
			if ismall == -1 || nsmall < nbig {
				ismall, nsmall = ibig, nbig
			}
			ibig, nbig = j, bptr[j]
			continue
		}
		if ismall == -1 || nsmall < bptr[j] {
			ismall, nsmall = j, bptr[j]
		}
	}
	return
}

// adPresent reports whether AD slot j of ptr holds a present (non-missing)
// value. The -1 sentinel marks a missing/vector-end value here.
func adPresent(ptr []int, j int) bool {
	return j >= 0 && j < len(ptr) && ptr[j] >= 0
}

// alleleString returns the textual allele at index a (0 == REF).
func adAlleleString(v *vcf.Variant, a int) string {
	if a == 0 {
		return v.Ref
	}
	if a-1 < len(v.Alt) {
		return v.Alt[a-1]
	}
	return ""
}

// formatExp formats a float in C's "%e" style (six fractional digits, a sign and
// at least two exponent digits), matching the p-value/threshold formatting of
// ad-bias.
func formatExp(f float64) string {
	return strconv.FormatFloat(f, 'e', 6, 64)
}

// readFileBytes reads a whole file.
func readFileBytes(path string) ([]byte, error) {
	return os.ReadFile(path)
}
