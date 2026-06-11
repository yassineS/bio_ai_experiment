// CLI runners for `bcftools cnv` and `bcftools csq`. Follows the
// project parity rule (docs/PARITY_ROADMAP.md "Definition of 1:1"):
// every documented upstream flag must parse cleanly. Flags whose
// underlying behaviour is deferred either no-op (when their default
// is no-op) or hard-reject with a roadmap pointer when explicitly set.
package main

import (
	"flag"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/tools/bcftools/pkg/bcftools"
)

const cnvUsage = `bcftools cnv - copy-number variation caller.

Usage:
  bcftools cnv [options] <in.vcf[.gz]|in.bcf>

Copy-number HMM over BAF/LRR with hidden states CN0/CN1/CN2/CN3, a
faithful port of upstream's vcfcnv.c. With -c the HMM runs the paired
16-state model. Output is the upstream summary.tab "RG" region rows on
stdout. Tracked in docs/PARITY_ROADMAP.md#bcftools.

I/O options:
  -c, --control-sample NAME      Control sample name.
  -s, --query-sample NAME        Query sample name.
  -o, --output-dir DIR           Output directory (required).
                                 The summary TSV is streamed to stdout regardless
                                 of this value; tracked in PARITY_ROADMAP.
  -p, --plot-threshold FLOAT     Plot threshold; this port emits no plots.
  -r, --regions LIST             Region list (post-filter).
  -R, --regions-file FILE        BED-like regions file.
      --regions-overlap 0|1|2    Accepted; always POS-in-region.
  -t, --targets LIST             Like -r but always a post-filter.
  -T, --targets-file FILE        BED-like targets file.
      --targets-overlap 0|1|2    Accepted; always POS-in-region.
  -v, --verbosity INT            Accepted; ignored.

HMM options:
  -a, --aberrant FLOAT[,FLOAT]   Fraction of aberrant cells (query,control).
  -b, --BAF-weight FLOAT         Relative BAF contribution.
  -d, --BAF-dev FLOAT[,FLOAT]    Expected BAF std-dev (query,control).
  -e, --err-prob FLOAT           Uniform error probability.
  -k, --LRR-dev FLOAT[,FLOAT]    Expected LRR std-dev (query,control).
  -l, --LRR-weight FLOAT         Relative LRR contribution.
  -L, --LRR-smooth-win INT       LRR moving-average window.
  -O, --optimize FLOAT           Estimate fraction of aberrant cells down to FLOAT.
  -P, --same-prob FLOAT          Prior of -s/-c being the same.
  -W, --baum-welch FLOAT         Baum-Welch convergence threshold (hidden upstream).
  -x, --xy-prob FLOAT            Transition probability.
      --AF-file FILE             Read allele frequencies from FILE
                                 (CHR<TAB>POS<TAB>REF,ALT<TAB>AF). Acts as a
                                 targets filter and sets per-site genotype freqs.

  -h, -?, --help                 Show this help.
      --version                  Show version.
`

func runCNV(args []string) int {
	fs := flag.NewFlagSet("bcftools cnv", flag.ContinueOnError)
	fs.SetOutput(io.Discard)

	var (
		controlSample  string
		querySample    string
		outputDir      string
		plotThreshold  float64
		regions        string
		regionsFile    string
		regionsOverlap int
		targets        string
		targetsFile    string
		targetsOverlap int
		verbosity      int

		aberrant     string
		bafWeight    float64
		bafDev       string
		errProb      float64
		lrrDev       string
		lrrWeight    float64
		lrrSmoothWin int
		optimize     float64
		sameProb     float64
		baumWelch    float64
		xyProb       float64
		afFile       string

		showHelp    bool
		showHelpAlt bool
		showVersion bool
	)
	cliflag.StringVar(fs, &controlSample, "c", "control-sample", "", "Control sample")
	cliflag.StringVar(fs, &querySample, "s", "query-sample", "", "Query sample")
	cliflag.StringVar(fs, &outputDir, "o", "output-dir", "", "Output dir")
	cliflag.Float64Var(fs, &plotThreshold, "p", "plot-threshold", 1e9, "Plot threshold")
	cliflag.StringVar(fs, &regions, "r", "regions", "", "Regions")
	cliflag.StringVar(fs, &regionsFile, "R", "regions-file", "", "Regions file")
	fs.IntVar(&regionsOverlap, "regions-overlap", 1, "")
	cliflag.StringVar(fs, &targets, "t", "targets", "", "Targets")
	cliflag.StringVar(fs, &targetsFile, "T", "targets-file", "", "Targets file")
	fs.IntVar(&targetsOverlap, "targets-overlap", 0, "")
	cliflag.IntVar(fs, &verbosity, "v", "verbosity", 0, "Verbosity")

	cliflag.StringVar(fs, &aberrant, "a", "aberrant", "1.0,1.0", "Aberrant fraction")
	cliflag.Float64Var(fs, &bafWeight, "b", "BAF-weight", 1.0, "BAF weight")
	cliflag.StringVar(fs, &bafDev, "d", "BAF-dev", "0.04,0.04", "BAF std dev")
	cliflag.Float64Var(fs, &errProb, "e", "err-prob", 1e-4, "Error prob")
	cliflag.StringVar(fs, &lrrDev, "k", "LRR-dev", "0.20,0.20", "LRR std dev")
	cliflag.Float64Var(fs, &lrrWeight, "l", "LRR-weight", 0.2, "LRR weight")
	cliflag.IntVar(fs, &lrrSmoothWin, "L", "LRR-smooth-win", 10, "LRR smoothing window")
	cliflag.Float64Var(fs, &optimize, "O", "optimize", 1.0, "Optimize fraction")
	cliflag.Float64Var(fs, &sameProb, "P", "same-prob", 0.5, "Same-prob prior")
	cliflag.Float64Var(fs, &baumWelch, "W", "baum-welch", 0, "Baum-Welch threshold (hidden upstream)")
	cliflag.Float64Var(fs, &xyProb, "x", "xy-prob", 1e-9, "Transition prob")
	fs.StringVar(&afFile, "AF-file", "", "AF tab file")

	fs.BoolVar(&showHelp, "h", false, "")
	fs.BoolVar(&showHelpAlt, "?", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVersion, "version", false, "")

	if err := parseFlags(fs, args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, cnvUsage)
		return 2
	}
	if showHelp || showHelpAlt {
		fmt.Print(cnvUsage)
		return 0
	}
	if showVersion {
		fmt.Println(version)
		return 0
	}
	// The HMM tuning knobs (bafWeight, errProb, lrrWeight, optimize,
	// sameProb, baumWelch, xyProb, lrrSmoothWin) feed CNVOptions below
	// and drive the HMM. --AF-file is passed through to CNVFile, which
	// uses it as a targets filter and per-site genotype-frequency source.
	// These remaining flags are accepted for CLI parity but have no
	// effect: --regions-overlap / --targets-overlap select the overlap
	// mode, and --verbosity is a logging level.
	_ = regionsOverlap
	_ = targetsOverlap
	_ = verbosity

	rest := fs.Args()
	if len(rest) == 0 {
		fmt.Fprintln(os.Stderr, "bcftools cnv: missing input file")
		fmt.Fprint(os.Stderr, cnvUsage)
		return 2
	}
	if outputDir == "" {
		fmt.Fprintln(os.Stderr, "bcftools cnv: -o/--output-dir is required (v1 streams to stdout but still requires the upstream flag for parity)")
		return 2
	}

	// Parse the FLOAT[,FLOAT] knobs (we only honour the first).
	bafDevVal, _, err := parseCNVPair(bafDev, "-d/--BAF-dev")
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}
	lrrDevVal, _, err := parseCNVPair(lrrDev, "-k/--LRR-dev")
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}
	aberrantQ, aberrantC, err := parseCNVPair(aberrant, "-a/--aberrant")
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}

	opts := bcftools.CNVOptions{
		QuerySample:     querySample,
		ControlSample:   controlSample,
		OutputDir:       outputDir,
		PlotThreshold:   plotThreshold,
		AberrantQuery:   aberrantQ,
		AberrantControl: aberrantC,
		BAFDev:          bafDevVal,
		LRRDev:          lrrDevVal,
		LRRWeight:       lrrWeight,
		LRRSmoothWin:    lrrSmoothWin,
		ErrProb:         errProb,
		XYProb:          xyProb,
		SameProb:        sameProb,
		BAFWeight:       bafWeight,
		Optimize:        optimize,
		BaumWelch:       baumWelch,
		AFFile:          afFile,
	}
	if regions != "" {
		opts.Regions = bcftools.SplitCommaList(regions)
	}
	if regionsFile != "" {
		regs, err := bcftools.LoadRegionsFile(regionsFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "bcftools cnv: %v\n", err)
			return 1
		}
		opts.Regions = append(opts.Regions, regs...)
		opts.RegionsFile = regionsFile
	}
	if targets != "" {
		opts.Targets = bcftools.SplitCommaList(targets)
	}
	if targetsFile != "" {
		regs, err := bcftools.LoadRegionsFile(targetsFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "bcftools cnv: %v\n", err)
			return 1
		}
		opts.Targets = append(opts.Targets, regs...)
		opts.TargetsFile = targetsFile
	}

	if _, err := bcftools.CNVFile(rest[0], os.Stdout, opts); err != nil {
		fmt.Fprintf(os.Stderr, "bcftools cnv: %v\n", err)
		return 1
	}
	return 0
}

// parseCNVPair parses an upstream FLOAT[,FLOAT] knob; the second value
// is optional and defaults to the first.
func parseCNVPair(s, name string) (float64, float64, error) {
	if s == "" {
		return 0, 0, nil
	}
	parts := strings.SplitN(s, ",", 2)
	a, err := strconv.ParseFloat(parts[0], 64)
	if err != nil {
		return 0, 0, fmt.Errorf("bcftools cnv: bad %s %q: %v", name, s, err)
	}
	b := a
	if len(parts) == 2 {
		b, err = strconv.ParseFloat(parts[1], 64)
		if err != nil {
			return 0, 0, fmt.Errorf("bcftools cnv: bad %s %q: %v", name, s, err)
		}
	}
	return a, b, nil
}

const csqUsage = `bcftools csq - predict variant consequences against a GFF.

Usage:
  bcftools csq [options] -f ref.fa -g anno.gff3 <in.vcf[.gz]|in.bcf>

The haplotype-aware consequence engine annotates an INFO/BCSQ tag and a
per-haplotype FORMAT/BCSQ bitmask (expand with the TBCSQ convert tag via
bcftools query). Output is VCF text (-O v), BGZF VCF (-O z), or BCF
(-O b|u). The remaining -l/--local-csq per-record caller is tracked in
docs/PARITY_ROADMAP.md#bcftools. The BCSQ tag has the form

  consequence|gene|transcript|biotype|strand|aa_change|dna_change

Required:
  -f, --fasta-ref FILE           Reference FASTA.
  -g, --gff-annot FILE           GFF3 annotation.

CSQ-specific options:
  -b, --brief-predictions        Brief amino-acid predictions (alias for -B 1).
  -B, --trim-protein-seq INT     Abbreviate amino-acid predictions to INT residues.
  -C, --genetic-code INT|l       NCBI translation table (0,1,2,3,5; "l" lists them).
  -c, --custom-tag STRING        INFO tag name (default BCSQ).
  -l, --local-csq                Local (non-haplotype-aware) calling.
  -n, --ncsq INT                 Maximum per-haplotype consequences per site [16].
  -p, --phase a|m|r|R|s          Accepted; v1 SNP classifier ignores phasing.
      --dump-gff FILE            Dump the parsed GFF (for debugging).
      --force                    Accepted; v1 has no sanity checks to skip.
      --unify-chr-names 0|LIST   Three comma-separated VCF,GFF,FAI prefixes,
                                 each '-' for none. "0" disables.

General options:
  -e, --exclude EXPR             Accepted; v1 ignores (every record is processed).
  -i, --include EXPR             Accepted; v1 ignores.
      --no-version               Accepted; v1 never appends a version line.
  -o, --output FILE              Output file (default stdout).
  -O, --output-type b|u|z|v      Output format: b (BCF), u (uncompressed BCF),
                                 z (BGZF VCF), v (VCF text). 't' is unsupported.
  -r, --regions LIST             Region list (post-filter in v1).
  -R, --regions-file FILE        BED-like regions file.
      --regions-overlap 0|1|2    Accepted; v1 ignores.
  -s, --samples -|LIST           Sample list. Accepted; v1 does not subset.
  -S, --samples-file FILE        Samples file.
  -t, --targets LIST             Like -r but always a post-filter.
  -T, --targets-file FILE        BED-like targets file.
      --targets-overlap 0|1|2    Accepted; v1 ignores.
  -@, --threads INT              Worker threads for parallel BGZF compression of -O z/-O b.
  -v, --verbose / --verbosity INT  Accepted; v1 ignores.
  -W, --write-index[=FMT]        Accepted; v1 never auto-indexes outputs.
  -q, --quiet                    Accepted; deprecated upstream.

  -h, -?, --help                 Show this help.
      --version                  Show version.
`

func runCSQ(args []string) int {
	fs := flag.NewFlagSet("bcftools csq", flag.ContinueOnError)
	fs.SetOutput(io.Discard)

	var (
		fastaRef       string
		gffAnnot       string
		brief          bool
		trimProtein    int
		geneticCode    string
		customTag      string
		localCSQ       bool
		ncsq           int
		phase          string
		dumpGFF        string
		force          bool
		unifyChrNames  string
		excludeExpr    string
		includeExpr    string
		noVersion      bool
		outputFile     string
		outputType     string
		regions        string
		regionsFile    string
		regionsOverlap int
		samples        string
		samplesFile    string
		targets        string
		targetsFile    string
		targetsOverlap int
		threads        int
		verbosity      int
		writeIndex     string
		quiet          bool
		showHelp       bool
		showHelpAlt    bool
		showVersion    bool
	)
	cliflag.StringVar(fs, &fastaRef, "f", "fasta-ref", "", "Reference FASTA (required)")
	cliflag.StringVar(fs, &gffAnnot, "g", "gff-annot", "", "GFF3 annotation (required)")
	cliflag.BoolVar(fs, &brief, "b", "brief-predictions", false, "Brief predictions (deprecated alias)")
	cliflag.IntVar(fs, &trimProtein, "B", "trim-protein-seq", 0, "Trim long predictions")
	cliflag.StringVar(fs, &geneticCode, "C", "genetic-code", "0", "Genetic code table")
	cliflag.StringVar(fs, &customTag, "c", "custom-tag", "BCSQ", "INFO tag")
	cliflag.BoolVar(fs, &localCSQ, "l", "local-csq", false, "Local CSQ")
	cliflag.IntVar(fs, &ncsq, "n", "ncsq", 15, "Max per-haplotype consequences")
	cliflag.StringVar(fs, &phase, "p", "phase", "r", "Phase handling")
	fs.StringVar(&dumpGFF, "dump-gff", "", "Dump parsed GFF")
	fs.BoolVar(&force, "force", false, "Skip sanity checks")
	fs.StringVar(&unifyChrNames, "unify-chr-names", "0", "Unify chrom names")
	cliflag.StringVar(fs, &excludeExpr, "e", "exclude", "", "Exclude expression")
	cliflag.StringVar(fs, &includeExpr, "i", "include", "", "Include expression")
	fs.BoolVar(&noVersion, "no-version", false, "Do not append version")
	cliflag.StringVar(fs, &outputFile, "o", "output", "", "Output file")
	cliflag.StringVar(fs, &outputType, "O", "output-type", "v", "Output type")
	cliflag.StringVar(fs, &regions, "r", "regions", "", "Regions")
	cliflag.StringVar(fs, &regionsFile, "R", "regions-file", "", "Regions file")
	fs.IntVar(&regionsOverlap, "regions-overlap", 1, "")
	cliflag.StringVar(fs, &samples, "s", "samples", "", "Samples list")
	cliflag.StringVar(fs, &samplesFile, "S", "samples-file", "", "Samples file")
	cliflag.StringVar(fs, &targets, "t", "targets", "", "Targets")
	cliflag.StringVar(fs, &targetsFile, "T", "targets-file", "", "Targets file")
	fs.IntVar(&targetsOverlap, "targets-overlap", 0, "")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "Worker threads for parallel BGZF compression")
	cliflag.IntVar(fs, &verbosity, "v", "verbosity", 1, "Verbosity (accepted, ignored)")
	fs.IntVar(&verbosity, "verbose", 1, "Verbose (alias for --verbosity)")
	cliflag.StringVar(fs, &writeIndex, "W", "write-index", "", "Auto-index outputs (accepted, ignored)")
	cliflag.BoolVar(fs, &quiet, "q", "quiet", false, "Quiet (deprecated)")
	fs.BoolVar(&showHelp, "h", false, "")
	fs.BoolVar(&showHelpAlt, "?", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVersion, "version", false, "")

	if err := parseFlags(fs, args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, csqUsage)
		return 2
	}
	if showHelp || showHelpAlt {
		fmt.Print(csqUsage)
		return 0
	}
	if showVersion {
		fmt.Println(version)
		return 0
	}
	// Hard-reject deferred flags.
	if deferred := checkCSQDeferred(checkCSQDeferredInputs{
		unifyChrNames: unifyChrNames,
		dumpGFF:       dumpGFF,
		outputType:    outputType,
		localCSQ:      localCSQ,
	}); deferred != "" {
		fmt.Fprintf(os.Stderr, "bcftools csq: %s is not implemented in v1; tracked in docs/PARITY_ROADMAP.md#bcftools\n", deferred)
		return 2
	}

	// -C/--genetic-code: "l"/"L" lists the supported tables and exits;
	// otherwise an integer selects an NCBI translation table.
	if geneticCode == "l" || geneticCode == "L" {
		fmt.Print(bcftools.GeneticCodeListing())
		return 0
	}
	geneticCodeID, err := strconv.Atoi(geneticCode)
	if err != nil {
		fmt.Fprintf(os.Stderr, "bcftools csq: could not parse -C/--genetic-code %q\n", geneticCode)
		return 2
	}
	if !bcftools.GeneticCodeKnown(geneticCodeID) {
		fmt.Fprintf(os.Stderr, "bcftools csq: -C/--genetic-code %d: no such table (supported: %s); more tables tracked in docs/PARITY_ROADMAP.md#bcftools\n", geneticCodeID, bcftools.GeneticCodeIDs())
		return 2
	}

	// -b/--brief-predictions is upstream's alias for -B 1; -B/--trim-protein-seq
	// takes an explicit length. Both feed args->brief_predictions, with the
	// later flag winning, so honour -B when given and fall back to -b.
	// Upstream rejects -B < 1, so do the same (0 stays "unset").
	if trimProtein < 0 {
		fmt.Fprintf(os.Stderr, "bcftools csq: could not parse -B/--trim-protein-seq %d\n", trimProtein)
		return 2
	}
	briefLen := trimProtein
	if briefLen == 0 && brief {
		briefLen = 1
	}

	if _, err := bcftools.ParseCSQPhase(phase); err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}

	// Silence "declared but not used".
	_ = excludeExpr
	_ = includeExpr
	_ = samples
	_ = samplesFile
	_ = regionsOverlap
	_ = targetsOverlap
	_ = verbosity
	_ = writeIndex
	_ = quiet
	_ = noVersion
	_ = force

	rest := fs.Args()
	if len(rest) == 0 {
		fmt.Fprintln(os.Stderr, "bcftools csq: missing input file")
		fmt.Fprint(os.Stderr, csqUsage)
		return 2
	}
	if fastaRef == "" {
		fmt.Fprintln(os.Stderr, "bcftools csq: -f/--fasta-ref is required")
		return 2
	}
	if gffAnnot == "" {
		fmt.Fprintln(os.Stderr, "bcftools csq: -g/--gff-annot is required")
		return 2
	}

	opts := bcftools.CSQOptions{
		FastaRef:       fastaRef,
		GFFAnnot:       gffAnnot,
		CustomTag:      customTag,
		LocalCSQ:       localCSQ,
		NCSQ:           ncsq,
		TrimProteinSeq: briefLen,
		GeneticCode:    geneticCodeID,
		Verbosity:      verbosity,
		Quiet:          quiet,
		Force:          force,
		NoVersion:      noVersion,
		IncludeExpr:    includeExpr,
		ExcludeExpr:    excludeExpr,
		DumpGFF:        dumpGFF,
		UnifyChrNames:  unifyChrNames,
		Threads:        threads,
	}
	if ph, err := bcftools.ParseCSQPhase(phase); err == nil {
		opts.Phase = ph
	}
	if of, err := bcftools.ParseOutputFormat(outputType); err == nil {
		opts.OutputFormat = of
	} else {
		fmt.Fprintf(os.Stderr, "bcftools csq: %v\n", err)
		return 2
	}
	if regions != "" {
		opts.Regions = bcftools.SplitCommaList(regions)
	}
	if regionsFile != "" {
		regs, err := bcftools.LoadRegionsFile(regionsFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "bcftools csq: %v\n", err)
			return 1
		}
		opts.Regions = append(opts.Regions, regs...)
		opts.RegionsFile = regionsFile
	}
	if targets != "" {
		opts.Targets = bcftools.SplitCommaList(targets)
	}
	if targetsFile != "" {
		regs, err := bcftools.LoadRegionsFile(targetsFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "bcftools csq: %v\n", err)
			return 1
		}
		opts.Targets = append(opts.Targets, regs...)
		opts.TargetsFile = targetsFile
	}
	if samples != "" {
		opts.Samples = bcftools.SplitCommaList(samples)
	}
	if samplesFile != "" {
		opts.SamplesFile = samplesFile
	}

	out, err := openOutFile(outputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "bcftools csq: %v\n", err)
		return 1
	}
	defer out.Close()

	if _, err := bcftools.CSQFile(rest[0], out, opts); err != nil {
		fmt.Fprintf(os.Stderr, "bcftools csq: %v\n", err)
		return 1
	}
	return 0
}

// checkCSQDeferredInputs is the surface of upstream flags that v1
// rejects (rather than silently accepts) because their behaviour is
// non-trivial and would be a misleading no-op.
type checkCSQDeferredInputs struct {
	unifyChrNames string
	dumpGFF       string
	outputType    string
	localCSQ      bool
}

func checkCSQDeferred(in checkCSQDeferredInputs) string {
	if in.localCSQ {
		// -l/--local-csq selects the per-record (non-haplotype-aware)
		// caller (upstream test_cds_local). The port only implements the
		// haplotype-aware path, so honouring -l would silently produce
		// haplotype-aware output under a flag that promises otherwise —
		// a misleading no-op. Reject until test_cds_local is ported.
		return "-l/--local-csq"
	}
	switch in.outputType {
	case "", "v", "z", "b", "u":
		// All four output formats are supported via openCSQOutput.
	default:
		return "-O " + in.outputType + " (expect v|z|b|u)"
	}
	// --unify-chr-names is honoured by the CSQ engine; no deferral needed.
	_ = in.unifyChrNames
	// --dump-gff is honoured by CSQFile / CSQOptions.DumpGFF; no deferral.
	_ = in.dumpGFF
	return ""
}
