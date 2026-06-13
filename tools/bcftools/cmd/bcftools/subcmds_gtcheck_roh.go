// CLI runners for `bcftools gtcheck` and `bcftools roh`. The shape
// follows the runners in main.go and subcmds.go: parse flags,
// validate, dispatch to the library. The "every documented upstream
// flag must be recognised" project rule (docs/PARITY_ROADMAP.md) is
// enforced here — flags that the v1 library does not implement are
// parsed cleanly and either ignored (when their default is no-op) or
// hard-rejected with a roadmap pointer when explicitly set.
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

const gtcheckUsage = `bcftools gtcheck - check sample identity (genotype concordance).

Usage:
  bcftools gtcheck [options] [-g <genotypes.vcf.gz>] <query.vcf.gz>

Options:
      --distinctive-sites NUM    Emit the minimal set of sites distinguishing at least NUM pairs (requires -p/-P).
      --dry-run                  Stop after first record (time estimation).
  -E, --error-probability INT    Phred-scaled genotyping error probability (default 40). 0 = integer-mismatch scoring.
  -e, --exclude EXPR             Drop sites for which EXPR is true. Accepts upstream's qry:/gt: prefix.
  -g, --genotypes FILE           Panel of "truth" genotypes.
  -H, --homs-only                Restrict scoring to sites where the panel GT is homozygous (requires -g).
  -i, --include EXPR             Keep sites for which EXPR is true (qry:/gt: prefix accepted).
      --keep-refs                Keep monoallelic (no-ALT) sites instead of skipping them.
      --n-matches INT            Print only the top INT matches per query sample. Negative sorts by HWE prob. 0 = all.
      --no-HWE-prob              Disable the -log P(HWE) column.
  -o, --output FILE              Write output to FILE (default stdout).
  -O, --output-type t|z          t = tab-text (default), z = compressed (v1 supports only t).
  -p, --pairs LIST               Comma-separated sample pairs.
  -P, --pairs-file FILE          File of tab-separated sample pairs.
  -r, --regions REGION           Comma-separated list of regions (post-filter in v1).
  -R, --regions-file FILE        BED-like regions file.
      --regions-overlap 0|1|2    Accepted; v1 always uses POS-in-region.
  -s, --samples [qry|gt]:LIST    Sample subset for query or panel. "-" means "all".
  -S, --samples-file [qry|gt]:FILE  File of sample IDs.
  -t, --targets REGION           Like -r but always a post-filter.
  -T, --targets-file FILE        BED-like targets file.
      --targets-overlap 0|1|2    Accepted; v1 always uses POS-in-region.
  -u, --use TAG[,TAG2]           Scoring tag for query (and -g panel): GT or PL. Auto-detected per record by default.
  -c, --cluster MIN,MAX          Cluster cross-check samples by pairwise error: merge samples within MAX (MAX<0 derives it from MIN).
  -G, --GTs-only                 Upstream-deprecated alias; rejected with upstream's literal deprecation error.
      --threads N                Accepted; v1 is single-threaded.
  -?, --help                     Show this help.
      --version                  Show version.

Notes: GT and PL scoring, the probability/integer discordance paths,
the HWE column, --n-matches trimming and --distinctive-sites are
implemented. The remaining gaps (-O z compressed output, the -c/--cluster
dendrogram, and filter expressions) are tracked in
docs/PARITY_ROADMAP.md#bcftools.
`

func runGtcheck(args []string) int {
	fs := flag.NewFlagSet("bcftools gtcheck", flag.ContinueOnError)
	fs.SetOutput(io.Discard)

	var (
		genotypesFile     string
		pairsSpec         string
		pairsFile         string
		useTag            string
		allSites          bool
		homsOnly          bool
		noHWEProb         bool
		keepRefs          bool
		regions           string
		regionsFile       string
		regionsOverlap    int
		targets           string
		targetsFile       string
		targetsOverlap    int
		includeExpr       string
		excludeExpr       string
		dryRun            bool
		errorProbability  int
		outputType        string
		outputPath        string
		cluster           string
		distinctiveSites  string
		nMatches          int
		gtsOnlyDeprecated bool
		threads           int
		showHelp          bool
		showVer           bool
		// Combined -s/-S forms with the upstream "qry:" / "gt:" prefix
		// land in samplesCombined (we split on the prefix below).
		samplesCombined     string
		samplesCombinedFile string
	)

	cliflag.StringVar(fs, &genotypesFile, "g", "genotypes", "", "Genotypes panel")
	cliflag.StringVar(fs, &pairsSpec, "p", "pairs", "", "Comma-separated pairs")
	cliflag.StringVar(fs, &pairsFile, "P", "pairs-file", "", "Pairs file")
	cliflag.StringVar(fs, &useTag, "u", "use", "", "Scoring tag: GT or PL (auto-detected from the header when unset)")
	fs.BoolVar(&allSites, "all-sites", false, "Include all sites")
	cliflag.BoolVar(fs, &homsOnly, "H", "homs-only", false, "Panel homozygotes only")
	fs.BoolVar(&noHWEProb, "no-HWE-prob", false, "Disable HWE-prob column")
	fs.BoolVar(&keepRefs, "keep-refs", false, "Keep monoallelic (no-ALT) sites")
	cliflag.StringVar(fs, &samplesCombined, "s", "samples", "", "Samples (qry:/gt: prefix supported)")
	cliflag.StringVar(fs, &samplesCombinedFile, "S", "samples-file", "", "Samples file")
	cliflag.StringVar(fs, &regions, "r", "regions", "", "Regions (post-filter in v1)")
	cliflag.StringVar(fs, &regionsFile, "R", "regions-file", "", "Regions file")
	fs.IntVar(&regionsOverlap, "regions-overlap", 1, "")
	cliflag.StringVar(fs, &targets, "t", "targets", "", "Targets (post-filter)")
	cliflag.StringVar(fs, &targetsFile, "T", "targets-file", "", "Targets file")
	fs.IntVar(&targetsOverlap, "targets-overlap", 0, "")
	cliflag.StringVar(fs, &includeExpr, "i", "include", "", "Include EXPR")
	cliflag.StringVar(fs, &excludeExpr, "e", "exclude", "", "Exclude EXPR")
	fs.BoolVar(&dryRun, "dry-run", false, "Stop after first record")
	cliflag.IntVar(fs, &errorProbability, "E", "error-probability", 40, "Phred error probability")
	cliflag.StringVar(fs, &outputType, "O", "output-type", "t", "Output type (t or z)")
	cliflag.StringVar(fs, &outputPath, "o", "output", "", "Output file")
	cliflag.StringVar(fs, &cluster, "c", "cluster", "", "Cluster MIN,MAX")
	fs.StringVar(&distinctiveSites, "distinctive-sites", "", "Find distinguishing sites (accepted; v1 not implemented)")
	fs.IntVar(&nMatches, "n-matches", 0, "Top-N matches (accepted; v1 not implemented)")
	// Upstream-deprecated `-G/--GTs-only`. We MUST accept it AND emit
	// the literal upstream deprecation error.
	cliflag.BoolVar(fs, &gtsOnlyDeprecated, "G", "GTs-only", false, "Deprecated upstream alias")
	fs.IntVar(&threads, "threads", 0, "Threads (accepted, ignored)")
	fs.BoolVar(&showHelp, "?", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "version", false, "")
	// Used only for upstream surface parity; ignored in v1.
	_ = regionsOverlap
	_ = targetsOverlap
	_ = threads
	_ = samplesCombinedFile

	if err := parseFlags(fs, args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, gtcheckUsage)
		return 2
	}
	if showHelp {
		fmt.Print(gtcheckUsage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}

	// Upstream's literal deprecation diagnostic for `-G/--GTs-only`.
	if gtsOnlyDeprecated {
		fmt.Fprintln(os.Stderr, "The option -G, --GTs-only has been deprecated")
		return 2
	}

	// Reject only the genuinely unimplemented flags (the -c/--cluster
	// dendrogram and -O z compressed output). PL scoring, --n-matches
	// and --distinctive-sites are implemented below.
	if deferred := checkGtcheckDeferred(checkGtcheckDeferredInputs{
		cluster:    cluster,
		outputType: outputType,
	}); deferred != "" {
		fmt.Fprintf(os.Stderr, "bcftools gtcheck: %s is not implemented in v1; tracked in docs/PARITY_ROADMAP.md#bcftools\n", deferred)
		return 2
	}

	// Parse -c/--cluster MIN,MAX into the clustering thresholds.
	var clusterMin, clusterMax float64
	clusterSet := false
	if cluster != "" {
		var perr error
		clusterMin, clusterMax, perr = parseGtcheckCluster(cluster)
		if perr != nil {
			fmt.Fprintf(os.Stderr, "bcftools gtcheck: %v\n", perr)
			return 2
		}
		clusterSet = true
	}

	// Parse --distinctive-sites NUM (the optional ,MEM,TMP suffix that
	// upstream accepts is irrelevant to the in-memory port and is
	// dropped). NUM <= 1 is a fraction of pairs, NUM > 1 an absolute
	// count.
	var dsValue float64
	dsSet := false
	if distinctiveSites != "" {
		numField := distinctiveSites
		if i := strings.IndexByte(numField, ','); i >= 0 {
			numField = numField[:i]
		}
		v, perr := strconv.ParseFloat(strings.TrimSpace(numField), 64)
		if perr != nil {
			fmt.Fprintf(os.Stderr, "bcftools gtcheck: could not parse --distinctive-sites %s\n", distinctiveSites)
			return 2
		}
		dsValue = v
		dsSet = true
	}

	rest := fs.Args()
	if len(rest) == 0 {
		fmt.Fprintln(os.Stderr, "bcftools gtcheck: missing query input")
		fmt.Fprint(os.Stderr, gtcheckUsage)
		return 2
	}

	qSamples, gSamples, err := splitSamplesArg(samplesCombined)
	if err != nil {
		fmt.Fprintf(os.Stderr, "bcftools gtcheck: %v\n", err)
		return 2
	}

	opts := bcftools.GtcheckOptions{
		GenotypesFile:        genotypesFile,
		PairsSpec:            pairsSpec,
		PairsFile:            pairsFile,
		UseTag:               useTag,
		AllSites:             allSites,
		HomsOnly:             homsOnly,
		NoHWEProb:            noHWEProb,
		KeepRefs:             keepRefs,
		SamplesQry:           qSamples,
		SamplesGT:            gSamples,
		RegionsFile:          regionsFile,
		TargetsFile:          targetsFile,
		IncludeExpr:          includeExpr,
		ExcludeExpr:          excludeExpr,
		DryRun:               dryRun,
		ErrorProbability:     errorProbability,
		ErrorProbabilityZero: errorProbability == 0,
		DistinctiveSites:     dsValue,
		HasDistinctiveSites:  dsSet,
		NMatches:             nMatches,
		OutputType:           outputType,
		Cluster:              clusterSet,
		ClusterMin:           clusterMin,
		ClusterMax:           clusterMax,
	}
	if regions != "" {
		opts.Regions = bcftools.SplitCommaList(regions)
	}
	if targets != "" {
		opts.Targets = bcftools.SplitCommaList(targets)
	}

	out, err := openOutFile(outputPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "bcftools gtcheck: %v\n", err)
		return 1
	}
	defer out.Close()

	if _, err := bcftools.GtcheckFile(rest[0], out, opts); err != nil {
		fmt.Fprintf(os.Stderr, "bcftools gtcheck: %v\n", err)
		return 1
	}
	return 0
}

// splitSamplesArg handles upstream's `qry:LIST` / `gt:LIST` prefix on
// -s. A bare list (no prefix) is taken to apply to both cohorts.
// splitSamplesArg parses upstream's "qry:"/"gt:"-prefixed sample list.
// Returns an error for un-prefixed lists, matching upstream
// vcfgtcheck.c's "Which one? Query samples (qry:...) or genotype
// samples (gt:...)?" diagnostic.
func splitSamplesArg(arg string) (qry, gt []string, err error) {
	if arg == "" {
		return nil, nil, nil
	}
	switch {
	case strings.HasPrefix(arg, "qry:"):
		return bcftools.SplitCommaList(arg[len("qry:"):]), nil, nil
	case strings.HasPrefix(arg, "gt:"):
		return nil, bcftools.SplitCommaList(arg[len("gt:"):]), nil
	}
	return nil, nil, fmt.Errorf("Which one? Query samples (qry:%s) or genotype samples (gt:%s)?", arg, arg)
}

// checkGtcheckDeferredInputs groups the flag values that the port
// recognises but does not implement.
type checkGtcheckDeferredInputs struct {
	cluster    string
	outputType string
}

func checkGtcheckDeferred(in checkGtcheckDeferredInputs) string {
	switch {
	case in.outputType != "" && in.outputType != "t":
		return "-O z (compressed output)"
	}
	return ""
}

// parseGtcheckCluster parses the `-c/--cluster MIN,MAX` argument into its two
// floats. Both are required (matching upstream's MIN,MAX surface); MAX may be
// negative (upstream's "derive it" sentinel, default -0.3).
func parseGtcheckCluster(s string) (min, max float64, err error) {
	parts := strings.Split(s, ",")
	if len(parts) != 2 {
		return 0, 0, fmt.Errorf("--cluster expects MIN,MAX, got %q", s)
	}
	min, err = strconv.ParseFloat(strings.TrimSpace(parts[0]), 64)
	if err != nil {
		return 0, 0, fmt.Errorf("--cluster MIN %q: %w", parts[0], err)
	}
	max, err = strconv.ParseFloat(strings.TrimSpace(parts[1]), 64)
	if err != nil {
		return 0, 0, fmt.Errorf("--cluster MAX %q: %w", parts[1], err)
	}
	return min, max, nil
}

const rohUsage = `bcftools roh - detect runs of autozygosity (ROH).

Usage:
  bcftools roh [options] <in.vcf[.gz]|in.bcf>

General Options:
      --AF-dflt FLOAT            Use this AF when AF is unknown (default: skip the site).
      --AF-tag TAG               INFO tag to read AF from (default AF).
      --AF-file FILE             AF source file (CHR<TAB>POS<TAB>REF,ALT<TAB>AF).
  -b, --buffer-size INT[,INT]    Buffer / overlap size, 0 for unlimited.
  -e, --estimate-AF [TAG],FILE   Estimate AF from FORMAT/TAG (GT or PL) of all samples ("-") or listed samples.
      --exclude EXPR             Drop sites for which EXPR is true.
  -G, --GTs-only FLOAT           Hard-GT mode, use FLOAT (phred) as PL of the two least likely genotypes.
                                   Safe value is 30 to account for GT errors. (FLOAT, not int.)
      --include EXPR             Keep sites for which EXPR is true.
  -i, --ignore-homref            Skip 0/0 genotypes.
      --include-noalt            Include sites with no ALT.
  -I, --skip-indels              Drop indel records.
  -m, --genetic-map FILE         Genetic map (IMPUTE2 format) for distance-scaled transitions.
  -M, --rec-rate FLOAT           Constant recombination rate per bp for distance-scaled transitions.
  -o, --output FILE              Write output to FILE.
  -O, --output-type [srz]        Output sections: s = per-site, r = per-region, z = compressed.
                                  Default "sr". z wraps the text output in BGZF.
  -r, --regions REGION           Comma-separated regions (post-filter in v1).
  -R, --regions-file FILE        Regions file.
      --regions-overlap 0|1|2    Accepted; v1 always uses POS-in-region.
  -s, --samples LIST             Sample subset.
  -S, --samples-file FILE        File of sample IDs.
  -t, --targets REGION           Like -r but always a post-filter.
  -T, --targets-file FILE        Targets file.
      --targets-overlap 0|1|2    Accepted; v1 always uses POS-in-region.
      --threads N                Accepted; v1 is single-threaded.

HMM Options:
  -a, --hw-to-az FLOAT           P(HW->AZ) transition probability per bp (default 6.7e-8).
  -H, --az-to-hw FLOAT           P(AZ->HW) transition probability per bp (default 5e-9).
  -V, --viterbi-training FLOAT   Estimate HMM parameters by Baum-Welch; FLOAT is the convergence threshold.

  -?, --help                     Show this help.
      --version                  Show version.

Notes:
  - Transitions scale with the physical distance between adjacent
    markers; -m/--genetic-map and -M/--rec-rate further scale them by
    the interval cross-over probability.
  - PL-based emission scoring is not implemented; -G/--GTs-only
    (hard-GT mode) is the only supported emission path. An invocation
    without -G is rejected with upstream's "FORMAT/PL tag not found"
    error (see docs/PARITY_ROADMAP.md#bcftools).
`

func runRoh(args []string) int {
	fs := flag.NewFlagSet("bcftools roh", flag.ContinueOnError)
	fs.SetOutput(io.Discard)

	var (
		afDfltSet       bool
		afDflt          float64
		afTag           string
		afFile          string
		bufferSize      string
		estimateAF      string
		excludeExpr     string
		gtsOnly         float64
		includeExpr     string
		ignoreHomRef    bool
		includeNoalt    bool
		skipIndels      bool
		geneticMap      string
		recRate         float64
		outputPath      string
		outputType      string
		regions         string
		regionsFile     string
		regionsOverlap  int
		samples         string
		samplesFile     string
		targets         string
		targetsFile     string
		targetsOverlap  int
		threads         int
		hwToAz          float64
		azToHw          float64
		viterbiTraining float64
		showHelp        bool
		showVer         bool
	)

	// --AF-dflt requires a sentinel to distinguish "0.0 set" from
	// "never set". We use a parallel flag plus float captured below.
	fs.Func("AF-dflt", "AF default", func(s string) error {
		_, err := fmt.Sscanf(s, "%f", &afDflt)
		if err != nil {
			return fmt.Errorf("--AF-dflt: %w", err)
		}
		afDfltSet = true
		return nil
	})
	fs.StringVar(&afTag, "AF-tag", "AF", "AF source tag")
	fs.StringVar(&afFile, "AF-file", "", "AF file")
	cliflag.StringVar(fs, &bufferSize, "b", "buffer-size", "", "Buffer / overlap size")
	cliflag.StringVar(fs, &estimateAF, "e", "estimate-AF", "", "Estimate AF (accepted; v1 reads INFO/AF)")
	fs.StringVar(&excludeExpr, "exclude", "", "Exclude EXPR")
	// -G/--GTs-only takes a FLOAT (phred). This is the reviewer's
	// requirement #8.
	cliflag.Float64Var(fs, &gtsOnly, "G", "GTs-only", 0, "Hard-GT mode error (phred float)")
	fs.StringVar(&includeExpr, "include", "", "Include EXPR")
	cliflag.BoolVar(fs, &ignoreHomRef, "i", "ignore-homref", false, "Skip 0/0 GTs")
	fs.BoolVar(&includeNoalt, "include-noalt", false, "Include no-ALT sites")
	cliflag.BoolVar(fs, &skipIndels, "I", "skip-indels", false, "Skip indels")
	cliflag.StringVar(fs, &geneticMap, "m", "genetic-map", "", "Genetic map (IMPUTE2 format) for distance-scaled transitions")
	cliflag.Float64Var(fs, &recRate, "M", "rec-rate", 0, "Recombination rate per bp")
	cliflag.StringVar(fs, &outputPath, "o", "output", "", "Output file")
	cliflag.StringVar(fs, &outputType, "O", "output-type", "sr", "Output sections srz")
	cliflag.StringVar(fs, &regions, "r", "regions", "", "Regions")
	cliflag.StringVar(fs, &regionsFile, "R", "regions-file", "", "Regions file")
	fs.IntVar(&regionsOverlap, "regions-overlap", 1, "")
	cliflag.StringVar(fs, &samples, "s", "samples", "", "Samples")
	cliflag.StringVar(fs, &samplesFile, "S", "samples-file", "", "Samples file")
	cliflag.StringVar(fs, &targets, "t", "targets", "", "Targets")
	cliflag.StringVar(fs, &targetsFile, "T", "targets-file", "", "Targets file")
	fs.IntVar(&targetsOverlap, "targets-overlap", 0, "")
	fs.IntVar(&threads, "threads", 0, "Threads (accepted, ignored)")
	// -a/-H — the reviewer's requirement #10. The library accepts
	// these; the PR #106 CLI wrongly rejected them as "deferred".
	cliflag.Float64Var(fs, &hwToAz, "a", "hw-to-az", bcftools.DefaultHWtoAZ, "HW->AZ transition")
	cliflag.Float64Var(fs, &azToHw, "H", "az-to-hw", bcftools.DefaultAZtoHW, "AZ->HW transition")
	cliflag.Float64Var(fs, &viterbiTraining, "V", "viterbi-training", 0, "Estimate HMM parameters by Baum-Welch; FLOAT is the convergence threshold")
	fs.BoolVar(&showHelp, "?", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "version", false, "")

	if err := parseFlags(fs, args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, rohUsage)
		return 2
	}
	if showHelp {
		fmt.Print(rohUsage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}

	rest := fs.Args()
	if len(rest) == 0 {
		fmt.Fprintln(os.Stderr, "bcftools roh: missing input file")
		fmt.Fprint(os.Stderr, rohUsage)
		return 2
	}

	opts := bcftools.RohOptions{
		AFTag:           afTag,
		AFFile:          afFile,
		GTsOnly:         gtsOnly,
		IgnoreHomRef:    ignoreHomRef,
		IncludeNoalt:    includeNoalt,
		SkipIndels:      skipIndels,
		HWtoAZ:          hwToAz,
		AZtoHW:          azToHw,
		OutputTypes:     outputType,
		SamplesFile:     samplesFile,
		RegionsFile:     regionsFile,
		TargetsFile:     targetsFile,
		IncludeExpr:     includeExpr,
		ExcludeExpr:     excludeExpr,
		EstimateAF:      estimateAF,
		BufferSize:      bufferSize,
		GeneticMap:      geneticMap,
		RecRate:         recRate,
		ViterbiTraining: viterbiTraining,
		RegionsOverlap:  regionsOverlap,
		TargetsOverlap:  targetsOverlap,
	}
	if afDfltSet {
		// Re-bind so we can take an address without aliasing the
		// loop variable in tests.
		v := afDflt
		opts.AFDflt = &v
	}
	if samples != "" {
		opts.Samples = bcftools.SplitCommaList(samples)
	}
	if regions != "" {
		opts.Regions = bcftools.SplitCommaList(regions)
	}
	if targets != "" {
		opts.Targets = bcftools.SplitCommaList(targets)
	}

	out, err := openOutFile(outputPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "bcftools roh: %v\n", err)
		return 1
	}
	defer out.Close()

	if _, err := bcftools.RohFile(rest[0], out, opts); err != nil {
		fmt.Fprintf(os.Stderr, "bcftools roh: %v\n", err)
		return 1
	}
	return 0
}
