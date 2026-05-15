// CLI runners for the bcftools subcommands added by the gtcheck / roh
// PR. The shape follows the runners in main.go and subcmds.go:
// parse flags, validate, dispatch to the library package. Every
// upstream flag is either implemented or accept-and-ignored / rejected
// with a `docs/PARITY_ROADMAP.md` pointer.
package main

import (
	"flag"
	"fmt"
	"io"
	"math"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/tools/bcftools/pkg/bcftools"
)

const gtcheckUsage = `bcftools gtcheck - genotype concordance between a query VCF and a panel.

Usage:
  bcftools gtcheck [options] <query.vcf[.gz]|query.bcf>

Required:
  -g, --genotypes FILE       Panel VCF/BCF.

Pair selection:
  -p, --pairs SPEC           Sample pairs to compare, "Q1,P1,Q2,P2" or "Q1:P1,Q2:P2".
                             Default: every query sample × every panel sample.
  -P, --pairs-file FILE      File of QUERY <tab|comma> PANEL lines.

Scoring:
  -u, --use {GT|PL|GL}       Scoring metric (default GT — only GT implemented in v1;
                             PL/GL rejected with a roadmap pointer).
  -G, --ngs-error FLOAT      Genotype-error rate (accepted; v1 ignores).

Filters:
  -r, --regions LIST         Region list (post-filter in v1).
  -R, --regions-file FILE    BED-like regions file.
  -t, --targets LIST         Like -r but always a post-filter.
  -T, --targets-file FILE    BED-like targets file.
  -i, --include EXPR         Accepted; v1 ignores per-record expressions.
  -e, --exclude EXPR         Accepted; v1 ignores per-record expressions.

I/O:
  -O, --output-type {t}      Output format (text TSV only in v1).
  -o, --output PATH          Output file (default stdout).
      --threads N            Accepted; v1 is single-threaded.
  -h, --help                 Show this help.
      --version              Show version.

Accepted-and-deferred upstream flags (parse cleanly; v1 rejects when set):
  --all-sites, --homs-only, --no-HWE-prob, --GTs-only,
  --pl-units, --dosage, --tags, --cluster, --normalize.
`

func runGtcheck(args []string) int {
	fs := flag.NewFlagSet("bcftools gtcheck", flag.ContinueOnError)
	fs.SetOutput(io.Discard)

	var (
		panelPath   string
		pairsSpec   string
		pairsFile   string
		useMode     string
		ngsError    float64
		regions     string
		regionsFile string
		targets     string
		targetsFile string
		includeExpr string
		excludeExpr string
		outputType  string
		outputPath  string
		threads     int
		showHelp    bool
		showVer     bool
		// Deferred upstream flags — parse cleanly, error on use.
		allSites  bool
		homsOnly  bool
		noHWEProb bool
		gtsOnly   string
		plUnits   string
		dosage    bool
		tagsFlag  string
		cluster   string
		normalize string
	)
	cliflag.StringVar(fs, &panelPath, "g", "genotypes", "", "Panel VCF/BCF")
	cliflag.StringVar(fs, &pairsSpec, "p", "pairs", "", "Sample pairs")
	cliflag.StringVar(fs, &pairsFile, "P", "pairs-file", "", "Sample-pairs file")
	cliflag.StringVar(fs, &useMode, "u", "use", "GT", "Scoring metric GT|PL|GL")
	cliflag.Float64Var(fs, &ngsError, "G", "ngs-error", 0, "NGS error rate (accepted, ignored)")
	cliflag.StringVar(fs, &regions, "r", "regions", "", "Region(s)")
	cliflag.StringVar(fs, &regionsFile, "R", "regions-file", "", "Regions file")
	cliflag.StringVar(fs, &targets, "t", "targets", "", "Targets")
	cliflag.StringVar(fs, &targetsFile, "T", "targets-file", "", "Targets file")
	cliflag.StringVar(fs, &includeExpr, "i", "include", "", "Include expression (accepted, ignored)")
	cliflag.StringVar(fs, &excludeExpr, "e", "exclude", "", "Exclude expression (accepted, ignored)")
	cliflag.StringVar(fs, &outputType, "O", "output-type", "t", "Output type (only 't' supported)")
	cliflag.StringVar(fs, &outputPath, "o", "output", "", "Output path")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "Threads (accepted, ignored)")
	fs.BoolVar(&allSites, "all-sites", false, "")
	fs.BoolVar(&homsOnly, "homs-only", false, "")
	fs.BoolVar(&noHWEProb, "no-HWE-prob", false, "")
	fs.StringVar(&gtsOnly, "GTs-only", "", "")
	fs.StringVar(&plUnits, "pl-units", "", "")
	fs.BoolVar(&dosage, "dosage", false, "")
	fs.StringVar(&tagsFlag, "tags", "", "")
	fs.StringVar(&cluster, "cluster", "", "")
	fs.StringVar(&normalize, "normalize", "", "")
	fs.BoolVar(&showHelp, "h", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "version", false, "")

	if err := fs.Parse(args); err != nil {
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

	// Reject deferred upstream paths with a clear roadmap pointer.
	if deferred := checkGtcheckDeferred(checkGtcheckDeferredInputs{
		allSites:  allSites,
		homsOnly:  homsOnly,
		noHWEProb: noHWEProb,
		gtsOnly:   gtsOnly,
		plUnits:   plUnits,
		dosage:    dosage,
		tagsFlag:  tagsFlag,
		cluster:   cluster,
		normalize: normalize,
	}); deferred != "" {
		fmt.Fprintf(os.Stderr, "bcftools gtcheck: %s not implemented in v1; tracked in docs/PARITY_ROADMAP.md#bcftools\n", deferred)
		return 2
	}
	// Output-type "t" is the only emit we support today; reject anything
	// else with a roadmap pointer rather than silently emitting TSV.
	if outputType != "" && outputType != "t" {
		fmt.Fprintf(os.Stderr, "bcftools gtcheck: -O %s not implemented in v1 (only 't' is supported); tracked in docs/PARITY_ROADMAP.md#bcftools\n", outputType)
		return 2
	}

	rest := fs.Args()
	if len(rest) == 0 {
		fmt.Fprintln(os.Stderr, "bcftools gtcheck: missing query file")
		fmt.Fprint(os.Stderr, gtcheckUsage)
		return 2
	}
	if panelPath == "" {
		fmt.Fprintln(os.Stderr, "bcftools gtcheck: -g/--genotypes is required")
		return 2
	}

	use, err := bcftools.ParseGtcheckUseMode(useMode)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}

	opts := bcftools.GtcheckOptions{
		PanelPath:   panelPath,
		PairsFile:   pairsFile,
		Use:         use,
		NGSError:    ngsError,
		RegionsFile: regionsFile,
		TargetsFile: targetsFile,
		IncludeExpr: includeExpr,
		ExcludeExpr: excludeExpr,
	}
	if pairsSpec != "" {
		pairs, err := bcftools.ParseGtcheckPairsSpec(pairsSpec)
		if err != nil {
			fmt.Fprintln(os.Stderr, err)
			return 2
		}
		opts.Pairs = append(opts.Pairs, pairs...)
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

type checkGtcheckDeferredInputs struct {
	allSites  bool
	homsOnly  bool
	noHWEProb bool
	gtsOnly   string
	plUnits   string
	dosage    bool
	tagsFlag  string
	cluster   string
	normalize string
}

func checkGtcheckDeferred(in checkGtcheckDeferredInputs) string {
	switch {
	case in.allSites:
		return "--all-sites"
	case in.homsOnly:
		return "--homs-only"
	case in.noHWEProb:
		return "--no-HWE-prob"
	case in.gtsOnly != "":
		return "--GTs-only"
	case in.plUnits != "":
		return "--pl-units"
	case in.dosage:
		return "--dosage"
	case in.tagsFlag != "":
		return "--tags"
	case in.cluster != "":
		return "--cluster"
	case in.normalize != "":
		return "--normalize"
	}
	return ""
}

const rohUsage = `bcftools roh - call runs of homozygosity via a two-state HMM.

Usage:
  bcftools roh [options] <in.vcf[.gz]|in.bcf>

Samples:
  -s, --samples SAMPLE       Single sample (or comma list).
  -S, --samples-file FILE    File of sample names, one per line.

Model:
  -G, --genotype-only INT    Use hard genotypes with the given phred error rate
                             (default 30 → 1e-3). PL-mode is accepted but v1
                             always uses GT.
      --AF-tag TAG           INFO tag carrying ALT-allele frequency (default AF).
      --AF-file FILE         External AF table (accepted; v1 reads AFTag/AFDflt).
      --AF-dflt FLOAT        Fallback allele frequency (default 0.4).

Filters:
  -r, --regions LIST         Region list (post-filter in v1).
  -R, --regions-file FILE    BED-like regions file.
  -t, --targets LIST         Like -r but always a post-filter.
  -T, --targets-file FILE    BED-like targets file.
  -i, --include EXPR         Accepted; v1 ignores per-record expressions.
  -e, --exclude EXPR         Accepted; v1 ignores per-record expressions.

Output:
  -O, --output-type {r|s|sr} Output mode (default r — regions).
  -o, --output PATH          Output file (default stdout).
      --threads N            Accepted; v1 is single-threaded.
  -h, --help                 Show this help.
      --version              Show version.

Accepted-and-deferred upstream flags (parse cleanly; v1 rejects when set):
  -M/--rec-rate, -V/--genetic-map, -a/--hw-to-az, -H/--az-to-hw (use direct probs instead),
  --buffer-size, --skip-indels, --include-noalt, --estimate-AF, -I/--skip-indels.
`

func runRoh(args []string) int {
	fs := flag.NewFlagSet("bcftools roh", flag.ContinueOnError)
	fs.SetOutput(io.Discard)

	var (
		samples     string
		samplesFile string
		genoErrInt  int
		afTag       string
		afFile      string
		afDflt      float64
		regions     string
		regionsFile string
		targets     string
		targetsFile string
		includeExpr string
		excludeExpr string
		outputType  string
		outputPath  string
		threads     int
		showHelp    bool
		showVer     bool
		// Deferred upstream knobs — parse cleanly, error on use.
		recRate    float64
		geneticMap string
		bufferSize string
		skipIndels bool
		includeNo  bool
		estimateAF bool
		hwToAzKnob float64
		azToHwKnob float64
	)
	cliflag.StringVar(fs, &samples, "s", "samples", "", "Sample(s) comma list")
	cliflag.StringVar(fs, &samplesFile, "S", "samples-file", "", "Samples file")
	cliflag.IntVar(fs, &genoErrInt, "G", "genotype-only", -1, "Phred genotype error (e.g. 30)")
	fs.StringVar(&afTag, "AF-tag", "AF", "")
	fs.StringVar(&afFile, "AF-file", "", "")
	fs.Float64Var(&afDflt, "AF-dflt", 0.4, "")
	cliflag.StringVar(fs, &regions, "r", "regions", "", "Region(s)")
	cliflag.StringVar(fs, &regionsFile, "R", "regions-file", "", "Regions file")
	cliflag.StringVar(fs, &targets, "t", "targets", "", "Targets")
	cliflag.StringVar(fs, &targetsFile, "T", "targets-file", "", "Targets file")
	cliflag.StringVar(fs, &includeExpr, "i", "include", "", "Include expression (accepted, ignored)")
	cliflag.StringVar(fs, &excludeExpr, "e", "exclude", "", "Exclude expression (accepted, ignored)")
	cliflag.StringVar(fs, &outputType, "O", "output-type", "r", "Output mode r|s|sr")
	cliflag.StringVar(fs, &outputPath, "o", "output", "", "Output path")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "Threads (accepted, ignored)")
	cliflag.Float64Var(fs, &recRate, "M", "rec-rate", 0, "")
	cliflag.StringVar(fs, &geneticMap, "V", "genetic-map", "", "")
	fs.StringVar(&bufferSize, "buffer-size", "", "")
	fs.BoolVar(&skipIndels, "skip-indels", false, "")
	fs.BoolVar(&includeNo, "include-noalt", false, "")
	fs.BoolVar(&estimateAF, "estimate-AF", false, "")
	cliflag.Float64Var(fs, &hwToAzKnob, "a", "hw-to-az", 0, "")
	cliflag.Float64Var(fs, &azToHwKnob, "H", "az-to-hw", 0, "")
	fs.BoolVar(&showHelp, "h", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "version", false, "")

	if err := fs.Parse(args); err != nil {
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

	// Reject deferred upstream paths with a clear roadmap pointer.
	if deferred := checkRohDeferred(checkRohDeferredInputs{
		recRate:    recRate,
		geneticMap: geneticMap,
		bufferSize: bufferSize,
		skipIndels: skipIndels,
		includeNo:  includeNo,
		estimateAF: estimateAF,
		afFile:     afFile,
		hwToAz:     hwToAzKnob,
		azToHw:     azToHwKnob,
	}); deferred != "" {
		fmt.Fprintf(os.Stderr, "bcftools roh: %s not implemented in v1; tracked in docs/PARITY_ROADMAP.md#bcftools\n", deferred)
		return 2
	}

	rest := fs.Args()
	if len(rest) == 0 {
		fmt.Fprintln(os.Stderr, "bcftools roh: missing input file")
		fmt.Fprint(os.Stderr, rohUsage)
		return 2
	}

	outMode, err := bcftools.ParseRohOutputMode(outputType)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}

	// Convert phred-encoded genotype-error to a linear rate. -G 30 → 1e-3.
	genoErr := 0.0
	if genoErrInt > 0 {
		genoErr = phredToProb(float64(genoErrInt))
	}

	opts := bcftools.RohOptions{
		SamplesFile:   samplesFile,
		GenotypeError: genoErr,
		AFTag:         afTag,
		AFFile:        afFile,
		AFDflt:        afDflt,
		RegionsFile:   regionsFile,
		TargetsFile:   targetsFile,
		IncludeExpr:   includeExpr,
		ExcludeExpr:   excludeExpr,
		Output:        outMode,
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

type checkRohDeferredInputs struct {
	recRate    float64
	geneticMap string
	bufferSize string
	skipIndels bool
	includeNo  bool
	estimateAF bool
	afFile     string
	hwToAz     float64
	azToHw     float64
}

func checkRohDeferred(in checkRohDeferredInputs) string {
	switch {
	case in.recRate != 0:
		return "-M/--rec-rate"
	case in.geneticMap != "":
		return "-V/--genetic-map"
	case in.bufferSize != "":
		return "--buffer-size"
	case in.skipIndels:
		return "--skip-indels"
	case in.includeNo:
		return "--include-noalt"
	case in.estimateAF:
		return "--estimate-AF"
	case in.afFile != "":
		return "--AF-file"
	case in.hwToAz != 0:
		return "-a/--hw-to-az"
	case in.azToHw != 0:
		return "-H/--az-to-hw"
	}
	return ""
}

// phredToProb converts a phred-encoded error to a linear probability.
// 0 -> 1, 10 -> 0.1, 20 -> 0.01, 30 -> 1e-3.
func phredToProb(q float64) float64 {
	if q <= 0 {
		return 1
	}
	return math.Pow(10, -q/10)
}
