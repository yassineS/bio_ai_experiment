// Native port of the upstream `parental-origin` plugin
// (plugins/parental-origin.c). It determines, for a trio over a CNV region, the
// parental origin of a deletion (-t del) or duplication (-t dup) by accumulating
// the log-probability of paternal vs maternal origin across informative SNP
// sites, then printing a tab-delimited summary (type, predicted_origin,
// quality, nmarkers) plus optional per-site DBG lines (-d).
//
// The numeric output is built from genotype likelihoods (gl = 10^(-0.1*PL),
// renormalised), the running sums log(ppat)/log(pmat), the final
// quality = 4.3429*|ppat-pmat|, and — for the dup branch — the regularized
// incomplete-beta binomial tails calc_binom_one_sided/two_sided. The binomial
// tails go through our faithful kfunc port (native_kfunc.go), which mirrors
// htslib's deterministic AS245 kf_lgamma + Lentz kf_betai rather than libm's
// lgamma, so the dup probabilities printed with %e match upstream byte-for-byte
// on linux/amd64. The plus/minus/times/divide/log accumulation likewise matches
// C exactly because the IEEE-754 operations are identical.
//
// parental-origin is a run()-style plugin (its options precede the input file
// with no `--` separator). It is registered as a fullPlugin so runNativePlugin
// hands it the whole invocation; it reads the input via ViewFile (with the
// plugin's own -r region, if any) exactly as the generic native pipeline does.
package bcftools

import (
	"bytes"
	"fmt"
	"io"
	"math"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() {
	registerNativePlugin("parental-origin", func() NativePlugin { return &parentalOriginPlugin{} })
}

// CNV type constants, matching CNV_DEL / CNV_DUP in parental-origin.c.
const (
	cnvDel = 0
	cnvDup = 1
)

// parentalOriginPlugin implements the `parental-origin` plugin end to end.
type parentalOriginPlugin struct{}

// Name returns the plugin name.
func (p *parentalOriginPlugin) Name() string { return "parental-origin" }

// About returns the one-line description, matching parental-origin.c about().
func (p *parentalOriginPlugin) About() string {
	return "Determine parental origin of a CNV region in a trio.\n"
}

// RunStyle reports that parental-origin is a run()-style plugin (options precede
// the input file, no `--` separator).
func (p *parentalOriginPlugin) RunStyle() bool { return true }

// FlagTakesValue reports whether one of parental-origin's value-taking flags
// consumes the following CLI token, used by the host to split the input-file
// positional out of the plugin options.
func (p *parentalOriginPlugin) FlagTakesValue(flag string) bool {
	switch flag {
	case "-b", "--min-binom-prob", "-e", "--exclude", "-i", "--include",
		"-p", "--pfm", "-r", "--region", "-t", "--type", "-v", "--verbosity":
		return true
	}
	return false
}

// Init satisfies NativePlugin; the real work runs in RunFull.
func (p *parentalOriginPlugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	return hdr, nil
}

// Process satisfies NativePlugin; never reached (RunFull owns the run).
func (p *parentalOriginPlugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	return []*vcf.Variant{v}, nil
}

// Destroy satisfies NativePlugin.
func (p *parentalOriginPlugin) Destroy() error { return nil }

// parentalOriginArgs holds the parsed plugin options.
type parentalOriginArgs struct {
	pfm        string
	cnvType    int
	debug      bool
	greedy     bool
	region     string
	minPBinom  float64
	filterExpr string
	filterExcl bool
	haveFilter bool
}

// RunFull executes parental-origin: parse options, read the input (with -r), run
// the per-record accumulation, and print the summary, matching run() in
// parental-origin.c.
func (p *parentalOriginPlugin) RunFull(opts PluginOptions, out io.Writer, stderr io.Writer) error {
	pa := parentalOriginArgs{minPBinom: 1e-2}
	args := opts.Args
	for i := 0; i < len(args); i++ {
		a := args[i]
		next := func() (string, error) {
			if i+1 >= len(args) {
				return "", fmt.Errorf("parental-origin: %s requires an argument", a)
			}
			i++
			return args[i], nil
		}
		switch a {
		case "-i", "--include", "-e", "--exclude":
			if pa.haveFilter {
				return fmt.Errorf("parental-origin: only one -i or -e expression can be given, and they cannot be combined")
			}
			v, err := next()
			if err != nil {
				return err
			}
			pa.filterExpr = v
			pa.filterExcl = a == "-e" || a == "--exclude"
			pa.haveFilter = true
		case "-p", "--pfm":
			v, err := next()
			if err != nil {
				return err
			}
			pa.pfm = v
		case "-r", "--region":
			v, err := next()
			if err != nil {
				return err
			}
			pa.region = v
		case "-t", "--type":
			v, err := next()
			if err != nil {
				return err
			}
			switch strings.ToLower(v) {
			case "dup":
				pa.cnvType = cnvDup
			case "del":
				pa.cnvType = cnvDel
			}
		case "-d", "--debug":
			pa.debug = true
		case "-g", "--greedy":
			pa.greedy = true
		case "-b", "--min-binom-prob":
			v, err := next()
			if err != nil {
				return err
			}
			f, perr := strconv.ParseFloat(v, 64)
			if perr != nil {
				return fmt.Errorf("parental-origin: could not parse: -b %s", v)
			}
			if f < 0 || f > 1 {
				return fmt.Errorf("parental-origin: expected value from the interval [0,1] with --min-binom-prob")
			}
			pa.minPBinom = f
		case "-v", "--verbosity":
			if _, err := next(); err != nil {
				return err
			}
		default:
			return fmt.Errorf("parental-origin: unsupported option %q", a)
		}
	}
	if pa.pfm == "" {
		return fmt.Errorf("parental-origin: missing the -p option")
	}

	// Read the input (transparently VCF/VCF.gz/BCF), applying the plugin's own
	// -r region exactly as upstream's bcf_sr_set_regions does.
	var regions []string
	if pa.region != "" {
		regions = []string{pa.region}
	}
	if len(opts.Regions) != 0 {
		regions = append(regions, opts.Regions...)
	}
	input := opts.InputFile
	if input == "" {
		input = "-"
	}
	var vcfText bytes.Buffer
	if _, err := ViewFile(input, &vcfText, ViewOptions{OutputFormat: OutputVCF, Regions: regions}, stderr); err != nil {
		return fmt.Errorf("parental-origin: reading input: %w", err)
	}
	r := vcf.NewReader(bytes.NewReader(vcfText.Bytes()))
	hdr, err := r.ReadHeader()
	if err != nil {
		return fmt.Errorf("parental-origin: %w", err)
	}

	// Resolve the trio (proband, father, mother). hts_readlist splits on commas;
	// exactly three names are required.
	names := strings.Split(pa.pfm, ",")
	if len(names) != 3 {
		return fmt.Errorf("parental-origin: expected three sample names with -t")
	}
	idx := sampleIndex(hdr)
	trio := [3]int{}
	for i, name := range names {
		j, ok := idx[name]
		if !ok {
			return fmt.Errorf("parental-origin: the sample is not present: %s", name)
		}
		trio[i] = j
	}

	// Required FORMAT tags.
	if !hasFormatHeader(hdr.MetaInfo, "PL") {
		return fmt.Errorf("parental-origin: the tag FORMAT/PL is not present in %s", input)
	}
	if !hasFormatHeader(hdr.MetaInfo, "AD") {
		return fmt.Errorf("parental-origin: the tag FORMAT/AD is not present in %s", input)
	}
	if !hasFormatHeader(hdr.MetaInfo, "GT") {
		return fmt.Errorf("parental-origin: the tag FORMAT/GT is not present in %s", input)
	}

	var filter *pluginFilter
	if pa.haveFilter {
		filter, err = newPluginFilterWithHeader(pa.filterExpr, pa.filterExcl, hdr)
		if err != nil {
			return fmt.Errorf("parental-origin: %w", err)
		}
	}

	st := &parentalOriginState{
		hdr:    hdr,
		trio:   trio,
		args:   pa,
		filter: filter,
		out:    out,
	}

	if pa.debug {
		if pa.cnvType == cnvDel {
			fmt.Fprintf(out, "# DBG: position; paternal probability; maternal probability; PLs of child, father, mother\n")
		} else {
			fmt.Fprintf(out, "# DBG: position; paternal probability; maternal probability; ADs of child, father, mother; PLs of child, father, mother\n")
		}
	}

	for {
		v, rerr := r.Read()
		if rerr == io.EOF {
			break
		}
		if rerr != nil {
			return fmt.Errorf("parental-origin: %w", rerr)
		}
		st.processRecord(v)
	}

	qual := 4.3429 * math.Abs(st.ppat-st.pmat)
	origin := "uncertain"
	if st.ppat > st.pmat {
		origin = "paternal"
	} else if st.ppat < st.pmat {
		origin = "maternal"
	}

	// Provenance banner: "# bcftools +parental-origin <opts...> <file>". The host
	// argv mirrors upstream's run()-style argv[0]=name.
	fmt.Fprintf(out, "# bcftools +%s", p.Name())
	for _, a := range opts.Args {
		fmt.Fprintf(out, " %s", a)
	}
	if input != "-" {
		fmt.Fprintf(out, " %s", input)
	}
	fmt.Fprintf(out, "\n")
	fmt.Fprintf(out, "# [1]type\t[2]predicted_origin\t[3]quality\t[4]nmarkers\n")
	cnvLabel := "del"
	if pa.cnvType == cnvDup {
		cnvLabel = "dup"
	}
	fmt.Fprintf(out, "%s\t%s\t%s\t%d\n", cnvLabel, origin, fmt.Sprintf("%f", qual), st.ntest)
	return nil
}

// parentalOriginState carries the running accumulators across records.
type parentalOriginState struct {
	hdr    *vcf.Header
	trio   [3]int
	args   parentalOriginArgs
	filter *pluginFilter
	out    io.Writer

	ppat, pmat float64
	ntest      int
}

// processRecord mirrors process_record() in parental-origin.c for a single
// record.
func (st *parentalOriginState) processRecord(v *vcf.Variant) {
	// Only biallelic SNPs.
	if len(v.Alt) != 1 {
		return
	}
	if variantTypeMask(v)&vtSNP == 0 {
		return
	}

	// Apply the -i/-e filter with the plugin's custom per-trio bookkeeping,
	// mirroring the FLT_INCLUDE / FLT_EXCLUDE branch in process_record().
	if st.filter != nil {
		siteMatch, raw, exclude := st.filter.rawSamples(v)
		if exclude {
			if siteMatch {
				if raw == nil {
					return
				}
				passSite := false
				for i := 0; i < 3; i++ {
					if raw[st.trio[i]] {
						raw[st.trio[i]] = false
					} else {
						raw[st.trio[i]] = true
						passSite = true
					}
				}
				if !passSite {
					return
				}
			} else {
				// pass_site false: keep, all trio members included.
				raw = nil
			}
			if raw != nil {
				for i := 0; i < 3; i++ {
					if !raw[st.trio[i]] {
						return
					}
				}
			}
		} else {
			if !siteMatch {
				return
			}
			if raw != nil {
				for i := 0; i < 3; i++ {
					if !raw[st.trio[i]] {
						return
					}
				}
			}
		}
	}

	pls, npl := parseFormatIntAll(v, "PL")
	if pls == nil || npl <= 0 {
		return
	}
	ads, nad := parseFormatIntAll(v, "AD")
	if ads == nil || nad <= 0 {
		// upstream prints a note and returns; we just skip (no parity-visible
		// output for our fixtures, and the AD header presence was checked).
		fmt.Fprintf(st.out, "The FORMAT/AD tag not present at %s:%d\n", v.Chrom, v.Pos)
		return
	}
	nAllele := len(v.Alt) + 1
	if npl != nAllele*(nAllele+1)/2 {
		fmt.Fprintf(st.out, "todo: not a diploid site at %s:%d: %d alleles, %d PLs\n", v.Chrom, v.Pos, nAllele, npl)
		return
	}

	// Gather per-member gl (normalised), dsg (alt-dosage from GT), ad[ref,alt].
	var gl [9]float64 // [child*3, father*3, mother*3]
	var dsg [3]int
	var ad [6]int // [child(2), father(2), mother(2)]
	rawPL := make([][3]int, 3)
	for i := 0; i < 3; i++ {
		src := pls[st.trio[i]]
		if len(src) < 3 {
			return
		}
		isum := 0
		sum := 0.0
		for j := 0; j < 3; j++ {
			if src[j] == intMissing {
				return
			}
			gl[3*i+j] = math.Pow(10, -0.1*float64(src[j]))
			sum += gl[3*i+j]
			isum += src[j]
			rawPL[i][j] = src[j]
		}
		if isum == 0 {
			return
		}
		for j := 0; j < 3; j++ {
			gl[3*i+j] /= sum
		}

		gt, ok := sampleGT(v, st.trio[i])
		if !ok || len(gt.alleles) == 0 {
			return
		}
		for _, a := range gt.alleles {
			if a == missingAllele {
				return
			}
			if a != 0 {
				dsg[i]++
			}
		}

		adSrc := ads[st.trio[i]]
		if len(adSrc) < 2 {
			return
		}
		ad[2*i] = adSrc[0]
		ad[2*i+1] = adSrc[1]
	}

	glP := gl[0:3]
	glF := gl[3:6]
	glM := gl[6:9]
	dsgP, dsgF, dsgM := dsg[0], dsg[1], dsg[2]
	adP := ad[0:2]

	if st.args.cnvType == cnvDel {
		if dsgP != 0 && dsgP != 2 {
			return
		}
		if dsgF == dsgM {
			return
		}
		if !st.args.greedy {
			if dsgF == 1 && dsgP == dsgM {
				return
			}
			if dsgM == 1 && dsgP == dsgF {
				return
			}
		}
		pmat := glP[0]*(0.5*glM[0]*glF[0]+2/3.*glM[0]*glF[1]+glM[0]*glF[2]+1/3.*glM[1]*glF[0]+0.5*glM[1]*glF[1]+glM[1]*glF[2]) +
			glP[2]*(0.5*glM[2]*glF[2]+2/3.*glM[2]*glF[1]+glM[2]*glF[0]+1/3.*glM[1]*glF[2]+0.5*glM[1]*glF[1]+glM[1]*glF[0])
		ppat := glP[0]*(0.5*glM[0]*glF[0]+2/3.*glM[1]*glF[0]+glM[2]*glF[0]+1/3.*glM[0]*glF[1]+0.5*glM[1]*glF[1]+glM[2]*glF[1]) +
			glP[2]*(0.5*glM[2]*glF[2]+2/3.*glM[1]*glF[2]+glM[0]*glF[2]+1/3.*glM[2]*glF[1]+0.5*glM[1]*glF[1]+glM[0]*glF[1])

		// NB: args->pmat accumulates log(ppat), args->ppat accumulates log(pmat)
		// — the swap is upstream's (probability of the deleted vs observed allele).
		st.pmat += math.Log(ppat)
		st.ppat += math.Log(pmat)
		st.ntest++

		if st.args.debug {
			fmt.Fprintf(st.out, "DBG\t%d\t%s\t%s\t", v.Pos, cExp(ppat), cExp(pmat))
			for i := 0; i < 3; i++ {
				for j := 0; j < 3; j++ {
					fmt.Fprintf(st.out, " %d", rawPL[i][j])
				}
				fmt.Fprintf(st.out, "\t")
			}
			fmt.Fprintf(st.out, "\n")
		}
	} else { // cnvDup
		if adP[0] == 0 || adP[1] == 0 {
			return
		}
		if adP[0] == adP[1] {
			return
		}
		if dsgP != 1 {
			return
		}
		if dsgF == dsgM {
			return
		}
		if st.args.minPBinom != 0 {
			if dsgF == 1 && ad[2] != 0 && ad[3] != 0 && calcBinomTwoSided(ad[2], ad[3], 0.5) < st.args.minPBinom {
				return
			}
			if dsgM == 1 && ad[4] != 0 && ad[5] != 0 && calcBinomTwoSided(ad[4], ad[5], 0.5) < st.args.minPBinom {
				return
			}
		}

		prra := glP[1] * calcBinomOneSided(adP[1], adP[0], 1/3., true)
		praa := glP[1] * calcBinomOneSided(adP[1], adP[0], 2/3., false)
		ppat := prra*(glM[1]*glF[0]+glM[2]*glF[0]+0.5*glM[1]*glF[1]+glM[2]*glF[1]) +
			praa*(glM[1]*glF[2]+glM[0]*glF[2]+0.5*glM[1]*glF[1]+glM[0]*glF[1])
		pmat := prra*(glM[0]*glF[1]+glM[0]*glF[2]+0.5*glM[1]*glF[1]+glM[1]*glF[2]) +
			praa*(glM[2]*glF[1]+glM[2]*glF[0]+0.5*glM[1]*glF[1]+glM[1]*glF[0])
		st.pmat += math.Log(pmat)
		st.ppat += math.Log(ppat)
		st.ntest++

		if st.args.debug {
			fmt.Fprintf(st.out, "DBG\t%d\t%s\t%s\t", v.Pos, cExp(ppat), cExp(pmat))
			for i := 0; i < 3; i++ {
				fmt.Fprintf(st.out, "%d %d\t", ad[2*i], ad[2*i+1])
			}
			for i := 0; i < 3; i++ {
				for j := 0; j < 3; j++ {
					fmt.Fprintf(st.out, " %d", rawPL[i][j])
				}
				fmt.Fprintf(st.out, "\t")
			}
			fmt.Fprintf(st.out, "\n")
		}
	}
}
