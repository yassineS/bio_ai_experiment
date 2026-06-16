// Native port of the upstream `trio-dnm3` plugin (plugins/trio-dnm3.c) for its
// NAIVE scoring model (`--use-NAIVE`, or any `--dnm-tag TAG:flag`). The NAIVE
// model determines de-novo mutations purely from FORMAT/GT by checking
// Mendelian inheritance: a site/trio is flagged when the child's genotype is
// incompatible with the parents'. It is the only trio-dnm3 mode with no
// floating-point dependence — the per-trio verdict is an integer table lookup
// (priors.denovo[fi][mi][ci]) over genotype indices, so the FORMAT/DNM=1 and
// FORMAT/VA (de-novo allele) annotations are byte-reproducible against upstream.
//
// The other models (DMM, ALM, DNG; the default is DMM) compute a phred/log
// de-novo score from a Dirichlet-multinomial / DeNovoGear likelihood over
// AD/PL/QS with libm pow/log/exp and kf_lgamma. Those primary outputs are
// libm-precision-dependent and are reported as a clean unsupported Init error
// rather than emitting a silently-divergent score.
//
// trio-dnm3 is a run()-style plugin (options precede the input file, no `--`).
// It is registered as a fullPlugin so runNativePlugin hands it the whole
// invocation; in NAIVE mode it streams the input through ViewFile and re-emits
// the records (with the DNM/VA tags added) in the requested -O container.
package bcftools

import (
	"fmt"
	"io"
	"math"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

func init() {
	registerNativePlugin("trio-dnm3", func() NativePlugin { return &trioDNM3Plugin{} })
}

// trio-dnm3 model constants, matching USE_* in trio-dnm3.c.
const (
	dnmUseDNG   = 1
	dnmUseALM   = 2
	dnmUseNaive = 3
	dnmUseDMM   = 4
)

// trio-dnm3 score-type constants (the subset NAIVE needs), matching the
// DNM_FLAG bit in trio-dnm3.c.
const dnmTypeFlag = 1

// trioDNM3Plugin implements the `trio-dnm3` plugin (NAIVE mode end to end; the
// float models are reported unsupported by RunFull).
type trioDNM3Plugin struct{}

// Name returns the plugin name.
func (p *trioDNM3Plugin) Name() string { return "trio-dnm3" }

// About returns the one-line description, matching trio-dnm3.c about().
func (p *trioDNM3Plugin) About() string {
	return "Screen variants for possible de-novo mutations in trios.\n"
}

// RunStyle reports that trio-dnm3 is a run()-style plugin.
func (p *trioDNM3Plugin) RunStyle() bool { return true }

// FlagTakesValue reports whether one of trio-dnm3's value-taking flags consumes
// the following CLI token, used by the host to split the input-file positional.
func (p *trioDNM3Plugin) FlagTakesValue(flag string) bool {
	switch flag {
	case "-p", "--pfm", "-P", "--ped", "-o", "--output", "-O", "--output-type",
		"-e", "--exclude", "-i", "--include", "-r", "--regions", "-R", "--regions-file",
		"-t", "--targets", "-T", "--targets-file", "--dnm-tag", "--va", "--vaf",
		"-u", "--use", "--mrate", "--pn", "--pns", "--max-QM", "--phi", "--noise-prior",
		"--np", "--ad", "--allelic-dropout", "--strand-bias", "--sb", "--min-vaf",
		"-X", "--chrX", "-m", "--min-score", "--regions-overlap", "--targets-overlap",
		"-v", "--verbosity":
		return true
	}
	return false
}

// Init satisfies NativePlugin; the real work runs in RunFull.
func (p *trioDNM3Plugin) Init(args []string, hdr *vcf.Header) (*vcf.Header, error) {
	return hdr, nil
}

// Process satisfies NativePlugin; never reached (RunFull owns the run).
func (p *trioDNM3Plugin) Process(v *vcf.Variant) ([]*vcf.Variant, error) {
	return []*vcf.Variant{v}, nil
}

// Destroy satisfies NativePlugin.
func (p *trioDNM3Plugin) Destroy() error { return nil }

// dnmTrio is a resolved trio (child, father, mother sample indices) plus the
// male-proband flag (chrX inheritance), mirroring trio_t in trio-dnm3.c.
type dnmTrio struct {
	child, father, mother int
	isMale                bool
}

// RunFull executes trio-dnm3: parse options, and in NAIVE mode stream the input
// re-emitting the DNM/VA annotations; in the float models report unsupported.
func (p *trioDNM3Plugin) RunFull(opts PluginOptions, out io.Writer, stderr io.Writer) error {
	var (
		pfm, pedFname             string
		dnmScoreTag               = ""
		dnmAlleleTag              = "VA"
		dnmVAFTag                 = "VAF"
		useModel                  = 0
		scoreType                 = 0
		strictlyNovel             bool
		chrXListStr               string
		filterExpr                string
		filterExclude, haveFilter bool
		outputFile                string
		outputType                = OutputVCF
		regions                   []string
		regionsFile               string

		// Float-model knobs, defaults from run() in trio-dnm3.c.
		useDNGPriors bool
		withPPL      bool
		withPAD      bool
		withCAD      bool
		forceAD      bool
		mrate        = 1e-8
		phi          = 1e3
		noisePrior   = 1e-3
		minVAF       = -1.0
		allelicDrop  = 0.0
		sbCoeff      = -1.0
		minScore     = math.Inf(-1)
		minQM        = phred2num(30)
		pnSNV        = pnoise{frac: 0.011, frac1: 0.045}
		pnIndel      pnoise
		pnIndelSet   bool
	)
	parsePN := func(s string, into *pnoise, isPNS bool) error {
		// FRAC[,NUM][:TYPE]
		typ := ""
		if c := strings.LastIndexByte(s, ':'); c >= 0 {
			typ = strings.ToLower(s[c+1:])
			s = s[:c]
		}
		frac, absv := 0.0, 0.0
		if c := strings.IndexByte(s, ','); c >= 0 {
			f, err := strconv.ParseFloat(s[:c], 64)
			if err != nil {
				return fmt.Errorf("trio-dnm3: could not parse --pn/--pns %q", s)
			}
			a, err := strconv.ParseFloat(s[c+1:], 64)
			if err != nil {
				return fmt.Errorf("trio-dnm3: could not parse --pn/--pns %q", s)
			}
			frac, absv = f, a
		} else {
			f, err := strconv.ParseFloat(s, 64)
			if err != nil {
				return fmt.Errorf("trio-dnm3: could not parse --pn/--pns %q", s)
			}
			frac = f
		}
		if frac < 0 || frac > 1 {
			return fmt.Errorf("trio-dnm3: expected value from the interval [0,1] for --pn/--pns")
		}
		set := func(p *pnoise) {
			if isPNS {
				p.frac1, p.abs1 = frac, absv
			} else {
				p.frac, p.abs = frac, absv
			}
		}
		switch typ {
		case "snv":
			set(&pnSNV)
		case "indel":
			set(&pnIndel)
			pnIndelSet = true
		case "":
			set(&pnSNV)
			set(&pnIndel)
			pnIndelSet = true
		default:
			return fmt.Errorf("trio-dnm3: could not parse --pn/--pns %q", typ)
		}
		_ = into
		return nil
	}

	args := opts.Args
	for i := 0; i < len(args); i++ {
		a := args[i]
		next := func() (string, error) {
			if i+1 >= len(args) {
				return "", fmt.Errorf("trio-dnm3: %s requires an argument", a)
			}
			i++
			return args[i], nil
		}
		switch a {
		case "-p", "--pfm":
			v, err := next()
			if err != nil {
				return err
			}
			pfm = v
		case "-P", "--ped":
			v, err := next()
			if err != nil {
				return err
			}
			pedFname = v
		case "--dnm-tag":
			v, err := next()
			if err != nil {
				return err
			}
			dnmScoreTag = v
		case "--va":
			v, err := next()
			if err != nil {
				return err
			}
			dnmAlleleTag = v
		case "-u", "--use":
			v, err := next()
			if err != nil {
				return err
			}
			// -u/--use accepts model names and the -u key=value knobs handled by
			// set_option() upstream; we support the model selectors and dng-priors.
			low := strings.ToLower(v)
			switch {
			case low == "naive":
				useModel = dnmUseNaive
			case low == "dmm":
				useModel = dnmUseDMM
			case low == "alm":
				useModel = dnmUseALM
			case low == "dng":
				useModel = dnmUseDNG
				useDNGPriors = true
			case low == "dng-priors":
				useDNGPriors = true
			case low == "ppl":
				withPPL = true
			case strings.HasPrefix(low, "mrate="):
				f, err := parseDNMFloat("-u mrate", v[len("mrate="):])
				if err != nil {
					return err
				}
				mrate = f
			default:
				return fmt.Errorf("trio-dnm3: the option \"-u %s\" is not recognised", v)
			}
		case "--use-NAIVE":
			useModel = dnmUseNaive
		case "--use-DMM":
			useModel = dnmUseDMM
		case "--use-ALM":
			useModel = dnmUseALM
		case "--use-DNG":
			useModel = dnmUseDNG
			useDNGPriors = true
		case "-n", "--strictly-novel":
			strictlyNovel = true
		case "-X", "--chrX":
			v, err := next()
			if err != nil {
				return err
			}
			chrXListStr = v
		case "-i", "--include", "-e", "--exclude":
			if haveFilter {
				return fmt.Errorf("trio-dnm3: only one -i or -e expression can be given, and they cannot be combined")
			}
			v, err := next()
			if err != nil {
				return err
			}
			filterExpr = v
			filterExclude = a == "-e" || a == "--exclude"
			haveFilter = true
		case "-o", "--output":
			v, err := next()
			if err != nil {
				return err
			}
			outputFile = v
		case "-O", "--output-type":
			v, err := next()
			if err != nil {
				return err
			}
			ot, oerr := parseDNMOutputType(v)
			if oerr != nil {
				return oerr
			}
			outputType = ot
		case "-r", "--regions":
			v, err := next()
			if err != nil {
				return err
			}
			regions = append(regions, v)
		case "-R", "--regions-file":
			v, err := next()
			if err != nil {
				return err
			}
			regionsFile = v
		case "--vaf":
			v, err := next()
			if err != nil {
				return err
			}
			dnmVAFTag = v
		case "--dng-priors":
			useDNGPriors = true
		case "--ppl", "--with-pPL", "--with-ppl":
			withPPL = true
		case "--with-pAD", "--with-pad":
			withPAD = true
		case "--with-cAD", "--with-cad":
			withCAD = true
		case "--force-AD":
			forceAD = true
		case "--no-version":
			// Provenance is stripped in tests; nothing to do here.
		case "--mrate":
			v, err := next()
			if err != nil {
				return err
			}
			if mrate, err = parseDNMFloat(a, v); err != nil {
				return err
			}
		case "--pn":
			v, err := next()
			if err != nil {
				return err
			}
			if err := parsePN(v, nil, false); err != nil {
				return err
			}
		case "--pns":
			v, err := next()
			if err != nil {
				return err
			}
			if err := parsePN(v, nil, true); err != nil {
				return err
			}
		case "--max-QM":
			v, err := next()
			if err != nil {
				return err
			}
			q, perr := parseDNMFloat(a, v)
			if perr != nil {
				return perr
			}
			if q < 0 {
				minQM = -phred2num(-q)
			} else {
				minQM = phred2num(q)
			}
			if math.Abs(minQM) > 1 {
				if minQM < 0 {
					minQM = -1
				} else {
					minQM = 1
				}
			}
		case "--phi":
			v, err := next()
			if err != nil {
				return err
			}
			if phi, err = parseDNMFloat(a, v); err != nil {
				return err
			}
			if phi <= 0 {
				return fmt.Errorf("trio-dnm3: expected positive values: --phi %s", v)
			}
		case "--noise-prior", "--np":
			v, err := next()
			if err != nil {
				return err
			}
			if noisePrior, err = parseDNMFloat(a, v); err != nil {
				return err
			}
			if math.Abs(noisePrior) > 1 {
				return fmt.Errorf("trio-dnm3: expected value [1,1] --noise-prior %s", v)
			}
		case "--ad", "--allelic-dropout":
			v, err := next()
			if err != nil {
				return err
			}
			if allelicDrop, err = parseDNMFloat(a, v); err != nil {
				return err
			}
		case "--strand-bias", "--sb":
			v, err := next()
			if err != nil {
				return err
			}
			if sbCoeff, err = parseDNMFloat(a, v); err != nil {
				return err
			}
		case "--min-vaf":
			v, err := next()
			if err != nil {
				return err
			}
			if minVAF, err = parseDNMFloat(a, v); err != nil {
				return err
			}
			if minVAF < 0 || minVAF > 0.5 {
				return fmt.Errorf("trio-dnm3: expected value [0,0.5]: --min-vaf %s", v)
			}
		case "-m", "--min-score":
			v, err := next()
			if err != nil {
				return err
			}
			if minScore, err = parseDNMFloat(a, v); err != nil {
				return err
			}
		case "--regions-overlap", "--targets-overlap", "-v", "--verbosity":
			// Value-taking options that do not affect the score; consume and ignore.
			if _, err := next(); err != nil {
				return err
			}
		case "-t", "--targets", "-T", "--targets-file":
			return fmt.Errorf("trio-dnm3: streaming targets (%s) are not supported by the native plugin; pre-slice with `bcftools view -t ... | bcftools +trio-dnm3`", a)
		default:
			// Attached getopt forms: -O<type>, --output-type=<type>, -o<file>.
			switch {
			case strings.HasPrefix(a, "-O") && len(a) > 2:
				ot, oerr := parseDNMOutputType(a[2:])
				if oerr != nil {
					return oerr
				}
				outputType = ot
			case strings.HasPrefix(a, "--output-type="):
				ot, oerr := parseDNMOutputType(a[len("--output-type="):])
				if oerr != nil {
					return oerr
				}
				outputType = ot
			case strings.HasPrefix(a, "-o") && len(a) > 2:
				outputFile = a[2:]
			case strings.HasPrefix(a, "--output="):
				outputFile = a[len("--output="):]
			default:
				return fmt.Errorf("trio-dnm3: unsupported option %q", a)
			}
		}
	}

	// Resolve the dnm score tag and type (mirrors init_data()): a bare "DNM"
	// defaults to :log (float); ":flag" implies NAIVE; NAIVE requires :flag.
	scoreTag := "DNM"
	floatScoreType := dnmScoreLog
	if dnmScoreTag != "" {
		if idx := strings.IndexByte(dnmScoreTag, ':'); idx >= 0 {
			if idx == 0 {
				return fmt.Errorf("trio-dnm3: could not parse --dnm-tag %s", dnmScoreTag)
			}
			scoreTag = dnmScoreTag[:idx]
			switch strings.ToLower(dnmScoreTag[idx+1:]) {
			case "flag":
				scoreType = dnmTypeFlag
			case "log":
				scoreType = 0
				floatScoreType = dnmScoreLog
			case "phred":
				scoreType = 0
				floatScoreType = dnmScorePhred
			case "prob":
				scoreType = 0
				floatScoreType = dnmScoreProb
			default:
				return fmt.Errorf("trio-dnm3: the type \"%s\" is not supported", dnmScoreTag[idx+1:])
			}
		} else {
			scoreTag = dnmScoreTag
		}
	}
	if scoreType == dnmTypeFlag {
		if useModel == 0 {
			useModel = dnmUseNaive
		} else if useModel != dnmUseNaive {
			return fmt.Errorf("trio-dnm3: the output type FLAG can be used only with --use-NAIVE")
		}
	}
	if useModel == dnmUseNaive {
		if scoreType == 0 {
			scoreType = dnmTypeFlag
		} else if scoreType != dnmTypeFlag {
			return fmt.Errorf("trio-dnm3: the output type FLAG is required with --use-NAIVE")
		}
	}
	if useModel == 0 {
		useModel = dnmUseDMM
	}
	if dnmScoreTag == "" && useModel == dnmUseNaive {
		scoreTag = "DNM"
	}

	if pedFname == "" && pfm == "" {
		return fmt.Errorf("trio-dnm3: missing the -p or -P option")
	}
	if pedFname != "" && pfm != "" {
		return fmt.Errorf("trio-dnm3: expected only -p or -P option, not both")
	}

	// Default-value resolution that depends on the model, mirroring init_data().
	if minVAF < 0 {
		if useModel == dnmUseDMM {
			minVAF = 0.2
		} else {
			minVAF = 0
		}
	}
	if sbCoeff < 0 {
		if useModel == dnmUseDMM {
			sbCoeff = 1e-2
		} else {
			sbCoeff = 0
		}
	}
	if useModel == dnmUseDNG && sbCoeff < 0 {
		sbCoeff = 0
	}
	sbCoeff = math.Abs(sbCoeff)
	if !pnIndelSet {
		pnIndel = pnoise{} // default --pns 0:indel --pn 0:indel
	}

	if outputFile != "" && outputFile != "-" {
		return fmt.Errorf("trio-dnm3: writing to a file (-o) is not supported by the native plugin; use stdout")
	}

	// Read input.
	if regionsFile != "" {
		regs, rerr := LoadRegionsFile(regionsFile)
		if rerr != nil {
			return rerr
		}
		regions = append(regions, regs...)
	}
	input := opts.InputFile
	if input == "" {
		input = "-"
	}
	hdr, variants, err := readPluginInput(PluginOptions{InputFile: input, Regions: regions}, stderr)
	if err != nil {
		return fmt.Errorf("trio-dnm3: %w", err)
	}

	// Header-level tag requirements, mirroring init_data().
	hasAD := hasFormatHeader(hdr.MetaInfo, "AD")
	hasSP := hasFormatHeader(hdr.MetaInfo, "SP")
	needPL := false
	needQS := false
	if useModel == dnmUseNaive {
		if !hasFormatHeader(hdr.MetaInfo, "GT") {
			return fmt.Errorf("trio-dnm3: the tag FORMAT/GT is not present in %s", input)
		}
	} else {
		needPL = true
		if withCAD && useModel == dnmUseDMM {
			needPL = false
		}
		if needPL && !hasFormatHeader(hdr.MetaInfo, "PL") {
			return fmt.Errorf("trio-dnm3: the tag FORMAT/PL is not present in %s", input)
		}
		if !hasAD && stderr != nil {
			fmt.Fprintf(stderr, "Warning: the tag FORMAT/AD is not present in %s, the output tag FORMAT/VAF will not be added\n", input)
		}
		if withPAD && !hasAD {
			return fmt.Errorf("trio-dnm3: no FORMAT/AD is present in %s, cannot run with --with-pAD", input)
		}
		if useModel == dnmUseALM && !withPPL && !withPAD {
			needQS = true
			if !hasFormatHeader(hdr.MetaInfo, "QS") {
				return fmt.Errorf("trio-dnm3: the FORMAT/QS tag is not present; add --with-pAD or --with-pPL, or generate QS with `bcftools mpileup -a AD,QS`")
			}
		}
		if useModel == dnmUseDMM {
			needQS = hasFormatHeader(hdr.MetaInfo, "QS")
			okData := true
			if !needQS && !hasAD {
				okData = false
			}
			if minQM >= 0 && !hasFormatHeader(hdr.MetaInfo, "QM") {
				okData = false
			}
			if !okData {
				return fmt.Errorf("trio-dnm3: the DMM method requires FORMAT/AD, QM and optionally FORMAT/SP; use --max-QM with a negative value to run without QM")
			}
		}
	}

	// Resolve trios.
	idx := sampleIndex(hdr)
	var trios []dnmTrio
	if pfm != "" {
		t, terr := parseDNMPFM(pfm, idx)
		if terr != nil {
			return fmt.Errorf("trio-dnm3: %w", terr)
		}
		trios = []dnmTrio{t}
	} else {
		trios, err = parseDNMPED(pedFname, idx)
		if err != nil {
			return fmt.Errorf("trio-dnm3: %w", err)
		}
		if len(trios) == 0 {
			return fmt.Errorf("trio-dnm3: no complete trio present")
		}
		if stderr != nil {
			plural := "s"
			if len(trios) == 1 {
				plural = ""
			}
			fmt.Fprintf(stderr, "Identified %d complete trio%s in the VCF file\n", len(trios), plural)
		}
	}

	var filter *pluginFilter
	if haveFilter {
		filter, err = newPluginFilterWithHeader(filterExpr, filterExclude, hdr)
		if err != nil {
			return fmt.Errorf("trio-dnm3: %w", err)
		}
	}

	// Build the chrX region matcher (default GRCh37 PAR-exclusion list).
	chrX := buildChrXMatcher(chrXListStr)

	// Float models (DMM/ALM/DNG): full pprob priors + the per-record float
	// pipeline. Scores are libm-tolerance-matched, not byte-exact (see the
	// proximity helper in numeric_parity_test.go).
	if useModel != dnmUseNaive {
		return p.runFloatModels(opts, out, hdr, variants, trios, filter, chrX, runFloatConfig{
			useModel: useModel, scoreType: floatScoreType, scoreTag: scoreTag,
			alleleTag: dnmAlleleTag, vafTag: dnmVAFTag, strictlyNovel: strictlyNovel,
			useDNGPriors: useDNGPriors, withPPL: withPPL, withPAD: withPAD, withCAD: withCAD,
			forceAD: forceAD, hasAD: hasAD, hasSP: hasSP, needPL: needPL, needQS: needQS,
			mrate: mrate, phi: phi, noisePrior: noisePrior, minVAF: minVAF, allelicDrop: allelicDrop,
			sbCoeff: sbCoeff, minScore: minScore, minQM: minQM, pnSNV: pnSNV, pnIndel: pnIndel,
			outputType: outputType,
		})
	}

	// Build the Mendelian priors tables.
	pa := newDNMPriors(strictlyNovel, autosomalPriors)
	pX := newDNMPriors(strictlyNovel, chrXPriors)
	pXX := newDNMPriors(strictlyNovel, chrXXPriors)

	// Augment the output header with the DNM and VA FORMAT lines (appended after
	// the existing meta, mirroring bcf_hdr_printf in init_data()).
	outHdr := &vcf.Header{Samples: hdr.Samples}
	outHdr.MetaInfo = append(outHdr.MetaInfo, hdr.MetaInfo...)
	outHdr.MetaInfo = append(outHdr.MetaInfo,
		fmt.Sprintf(`##FORMAT=<ID=%s,Number=1,Type=Integer,Description="De-novo mutation score given as 1 for Mendelian-incompatible genotypes">`, scoreTag))
	outHdr.MetaInfo = append(outHdr.MetaInfo,
		fmt.Sprintf(`##FORMAT=<ID=%s,Number=1,Type=Integer,Description="The de-novo allele">`, dnmAlleleTag))

	st := &dnmNaiveState{
		hdr:          hdr,
		trios:        trios,
		chrX:         chrX,
		priors:       pa,
		priorsX:      pX,
		priorsXX:     pXX,
		scoreTag:     scoreTag,
		dnmAlleleTag: dnmAlleleTag,
		filter:       filter,
		trioPass:     make([]bool, len(trios)),
	}

	results := make([]*vcf.Variant, 0, len(variants))
	for _, v := range variants {
		results = append(results, st.processRecord(v))
	}

	w, cleanup, err := openOutput(out, ViewOptions{OutputFormat: outputType, CompressLevel: opts.CompressLevel, Threads: opts.Threads}, outHdr)
	if err != nil {
		return err
	}
	if err := w.WriteHeader(); err != nil {
		cleanup()
		return err
	}
	for _, v := range results {
		if err := w.Write(v); err != nil {
			cleanup()
			return err
		}
	}
	if err := w.Flush(); err != nil {
		cleanup()
		return err
	}
	cleanup()
	return nil
}

// parseDNMOutputType maps the -O/--output-type value to a container, mirroring
// the leading-character switch in run(). A trailing compression digit is
// accepted and ignored (NAIVE output is integer-only; the level is irrelevant
// to parity once decoded).
func parseDNMOutputType(v string) (OutputFormat, error) {
	if v == "" {
		return OutputVCF, fmt.Errorf("trio-dnm3: the output type \"\" not recognised")
	}
	switch v[0] {
	case 'v':
		return OutputVCF, nil
	case 'z':
		return OutputVCFGz, nil
	case 'u':
		return OutputBCFUncompressed, nil
	case 'b':
		return OutputBCF, nil
	default:
		return OutputVCF, fmt.Errorf("trio-dnm3: the output type \"%s\" not recognised", v)
	}
}

// parseDNMPFM parses the -p/--pfm "C,F,M" value (child, father, mother), with
// the optional 1X:/2X: prefix on the child marking a male proband, mirroring the
// pfm branch of init_data() in trio-dnm3.c.
func parseDNMPFM(pfm string, idx map[string]int) (dnmTrio, error) {
	parts := strings.Split(pfm, ",")
	if len(parts) != 3 {
		return dnmTrio{}, fmt.Errorf("expected three sample names with -t")
	}
	var t dnmTrio
	childName := parts[0]
	if c, ok := idx[childName]; ok {
		t.child = c
	} else if len(childName) > 3 && strings.EqualFold(childName[:3], "1X:") {
		c, ok := idx[childName[3:]]
		if !ok {
			return dnmTrio{}, fmt.Errorf("the sample is not present: %s", childName[3:])
		}
		t.child = c
		t.isMale = true
	} else if len(childName) > 3 && strings.EqualFold(childName[:3], "2X:") {
		c, ok := idx[childName[3:]]
		if !ok {
			return dnmTrio{}, fmt.Errorf("the sample is not present: %s", childName[3:])
		}
		t.child = c
	} else {
		return dnmTrio{}, fmt.Errorf("the sample is not present: %s", childName)
	}
	f, ok := idx[parts[1]]
	if !ok {
		return dnmTrio{}, fmt.Errorf("the sample is not present: %s", parts[1])
	}
	t.father = f
	m, ok := idx[parts[2]]
	if !ok {
		return dnmTrio{}, fmt.Errorf("the sample is not present: %s", parts[2])
	}
	t.mother = m
	return t, nil
}

// parseDNMPED reads a PED file and returns the complete trios (father, mother,
// child all present in the header), sorted by minimum sample index and checked
// for duplicates, mirroring parse_ped() + cmp_trios() in trio-dnm3.c. The 5th
// (sex) column, when 1, marks a male proband.
func parseDNMPED(path string, idx map[string]int) ([]dnmTrio, error) {
	rows, err := parsePEDRows(path, idx)
	if err != nil {
		return nil, err
	}
	trios := make([]dnmTrio, 0, len(rows))
	for _, r := range rows {
		trios = append(trios, dnmTrio{child: r.child, father: r.father, mother: r.mother, isMale: r.sex == 1})
	}
	// Sort by minimum sample index (cmp_trios).
	for i := 1; i < len(trios); i++ {
		for j := i; j > 0 && dnmMinIdx(trios[j]) < dnmMinIdx(trios[j-1]); j-- {
			trios[j], trios[j-1] = trios[j-1], trios[j]
		}
	}
	for i := 1; i < len(trios); i++ {
		a, b := trios[i-1], trios[i]
		if a.child == b.child && a.father == b.father && a.mother == b.mother {
			return nil, fmt.Errorf("duplicate trio entries detected in the PED file: %s", path)
		}
	}
	return trios, nil
}

// dnmMinIdx returns the minimum of a trio's three sample indices (father,
// mother, child), the cmp_trios() sort key in trio-dnm3.c.
func dnmMinIdx(t dnmTrio) int {
	m := t.father
	if t.mother < m {
		m = t.mother
	}
	if t.child < m {
		m = t.child
	}
	return m
}

// dnmNaiveState carries the per-record NAIVE transform state.
type dnmNaiveState struct {
	hdr                    *vcf.Header
	trios                  []dnmTrio
	chrX                   *chrXMatcher
	priors                 *dnmPriors
	priorsX, priorsXX      *dnmPriors
	scoreTag, dnmAlleleTag string
	filter                 *pluginFilter // compiled -i/-e per-trio pre-filter, nil if none
	trioPass               []bool        // per-trio pass flag for the current record
}

// testFilters applies the -i/-e expression with trio-dnm3's per-trio
// bookkeeping for the NAIVE model; see dnmRunFilters.
func (st *dnmNaiveState) testFilters(v *vcf.Variant) bool {
	return dnmRunFilters(st.filter, st.trios, st.trioPass, v)
}

// dnmRunFilters applies the -i/-e expression with trio-dnm3's per-trio
// bookkeeping, mirroring test_filters() in trio-dnm3.c. It returns whether the
// site has at least one passing trio and fills trioPass. A nil filter passes
// every trio. Shared by the NAIVE and float model states.
func dnmRunFilters(filter *pluginFilter, trios []dnmTrio, trioPass []bool, v *vcf.Variant) bool {
	if filter == nil {
		for i := range trioPass {
			trioPass[i] = true
		}
		return true
	}
	siteMatch, raw, exclude := filter.rawSamples(v)
	if exclude {
		if siteMatch {
			if raw == nil {
				return false // -e mode, the expression failed at site level
			}
			passSite := false
			for i, tr := range trios {
				passTrio := true
				for _, idx := range [3]int{tr.child, tr.father, tr.mother} {
					if raw[idx] { // with -e one sample passing the expr fails the trio
						passTrio = false
						break
					}
				}
				trioPass[i] = passTrio
				if passTrio {
					passSite = true
				}
			}
			return passSite
		}
		for i := range trioPass {
			trioPass[i] = true
		}
		return true
	}
	if !siteMatch {
		return false
	}
	if raw != nil {
		passSite := false
		for i, tr := range trios {
			passTrio := true
			for _, idx := range [3]int{tr.child, tr.father, tr.mother} {
				if !raw[idx] {
					passTrio = false
					break
				}
			}
			trioPass[i] = passTrio
			if passTrio {
				passSite = true
			}
		}
		return passSite
	}
	for i := range trioPass {
		trioPass[i] = true
	}
	return true
}

// processRecord adds the DNM/VA annotations to a record when any trio shows a
// Mendelian-incompatible (de-novo) child genotype, mirroring
// process_record()/process_record_naive() in trio-dnm3.c. It returns the record
// (possibly mutated in place).
func (st *dnmNaiveState) processRecord(v *vcf.Variant) *vcf.Variant {
	nAllele := len(v.Alt) + 1
	// Skip ref-only sites (no ALT, or ALT is "."), matching n_allele==1 ||
	// bcf_get_variant_types(rec)==VCF_REF.
	if nAllele == 1 || (len(v.Alt) == 1 && (v.Alt[0] == "." || v.Alt[0] == "")) {
		return v
	}
	if !formatHasTag(v, "GT") {
		return v
	}
	// Apply the -i/-e per-trio filter (mirrors the test_filters() gate before
	// process_record_naive in process_record). A site with no passing trio is
	// emitted unchanged.
	if !st.testFilters(v) {
		return v
	}

	isChrX := st.chrX.overlaps(v.Chrom, v.Pos)

	type dnmOut struct {
		sample int
		score  int
		allele int
	}
	var outs []dnmOut
	writeDNM := false
	for ti, tr := range st.trios {
		if st.filter != nil && !st.trioPass[ti] {
			continue
		}
		ignoreFather := false
		priors := st.priors
		if isChrX {
			if tr.isMale {
				priors = st.priorsX
				ignoreFather = true
			} else {
				priors = st.priorsXX
			}
		}
		gts, ok := setTrioGTNaive(v, tr, nAllele, ignoreFather)
		if !ok {
			continue
		}
		fi := dnmSeq3[gts[1]] // father
		mi := dnmSeq3[gts[2]] // mother
		ci := dnmSeq3[gts[0]] // child
		if fi < 0 || mi < 0 || ci < 0 {
			continue
		}
		isDNM := priors.denovo[fi][mi][ci]
		allele := priors.denovoAllele[fi][mi][ci]
		outs = append(outs, dnmOut{sample: tr.child, score: isDNM, allele: allele})
		if isDNM != 0 {
			writeDNM = true
		}
	}
	if !writeDNM {
		return v
	}

	// Per-sample DNM/VA: only the trio children get a value; everyone else is
	// missing (".").
	scoreVals := make(map[int]string)
	alleleVals := make(map[int]string)
	for _, o := range outs {
		scoreVals[o.sample] = strconv.Itoa(o.score)
		alleleVals[o.sample] = strconv.Itoa(o.allele)
	}
	v.Format = append(v.Format, st.scoreTag, st.dnmAlleleTag)
	for s := range v.Samples {
		if sv, ok := scoreVals[s]; ok {
			v.Samples[s].Data[st.scoreTag] = sv
			v.Samples[s].Data[st.dnmAlleleTag] = alleleVals[s]
		}
		// Samples not in scoreVals are left without the tag; the writer emits ".".
	}
	return v
}

// setTrioGTNaive builds the per-member allele bitmask (1<<a)|(1<<b) for the
// trio, mirroring set_trio_GT()/set_trio_GT_many_alts(). It returns the three
// masks ordered [child, father, mother] and ok=false when any required member's
// GT is missing (or, for >4 alleles, more than four distinct alleles are used).
func setTrioGTNaive(v *vcf.Variant, tr dnmTrio, nAllele int, ignoreFather bool) (gts [3]int, ok bool) {
	order := [3]struct {
		idx      int
		isFather bool
		outSlot  int // 0=child,1=father,2=mother in gts (matches seq3 lookup below)
	}{
		{tr.father, true, 1},
		{tr.mother, false, 2},
		{tr.child, false, 0},
	}
	if nAllele <= 4 {
		for _, m := range order {
			mask, mok := dnmGTMask(v, m.idx)
			if !mok {
				if m.isFather && ignoreFather {
					// Missing father GT allowed for a male proband on chrX; the mask
					// is irrelevant (chrX priors ignore the father), so use a dummy
					// homozygous-ref mask to keep the seq3 lookup valid.
					gts[m.outSlot] = 1
					continue
				}
				return gts, false
			}
			gts[m.outSlot] = mask
		}
		return gts, true
	}
	// >4 alleles: remap used alleles to compact 0..3 indices, mirroring
	// set_trio_GT_many_alts. A trio using >4 distinct alleles is rejected.
	altIdx := make([]int, nAllele)
	for i := range altIdx {
		altIdx[i] = -1
	}
	nused := 0
	for _, m := range order {
		alleles, present := dnmGTAlleles(v, m.idx)
		if !present {
			if m.isFather && ignoreFather {
				gts[m.outSlot] = 1
				continue
			}
			return gts, false
		}
		mask := 0
		for _, ial := range alleles {
			if ial < 0 { // missing within GT
				if m.isFather && ignoreFather {
					ial = 0
				} else {
					return gts, false
				}
			}
			if altIdx[ial] == -1 {
				altIdx[ial] = nused
				nused++
				if nused > 4 {
					return gts, false
				}
			}
			mask |= 1 << altIdx[ial]
		}
		if mask == 0 {
			return gts, false
		}
		gts[m.outSlot] = mask
	}
	return gts, true
}

// dnmGTMask returns (1<<a)|(1<<b)|... for sample i's GT, or ok=false if the
// whole GT is missing. A "." within the GT makes the sample unusable here.
func dnmGTMask(v *vcf.Variant, i int) (int, bool) {
	gt, parsed := sampleGT(v, i)
	if !parsed || len(gt.alleles) == 0 {
		return 0, false
	}
	mask := 0
	for _, a := range gt.alleles {
		if a == missingAllele {
			return 0, false
		}
		mask |= 1 << a
	}
	if mask == 0 {
		return 0, false
	}
	return mask, true
}

// dnmGTAlleles returns the per-allele indices of sample i's GT (missingAllele
// for "."), and present=false if the sample has no GT at all.
func dnmGTAlleles(v *vcf.Variant, i int) ([]int, bool) {
	gt, parsed := sampleGT(v, i)
	if !parsed || len(gt.alleles) == 0 {
		return nil, false
	}
	return gt.alleles, true
}
