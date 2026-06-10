// CLI runners for the bcftools subcommands added by the convert /
// mendelian PR. The shape follows the runners in main.go and
// subcmds.go: parse flags, validate, dispatch to the library package.
package main

import (
	"flag"
	"fmt"
	"io"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/tools/bcftools/pkg/bcftools"
)

const convertUsage = `bcftools convert - re-emit VCF/BCF in a different format.

Usage:
  bcftools convert [options] <in.vcf[.gz]|in.bcf>

Options:
  -O, --output-type {v|z|u|b}  Output format (v=VCF, z=VCF.gz, u=uncompressed BCF, b=compressed BCF).
  -o, --output PATH            Output file (default stdout).
  -s, --samples LIST           Restrict per-sample columns to these names (comma list).
  -S, --samples-file FILE      File with sample IDs (one per line).
      --force-samples          Allow --samples names that are missing from the input header.
  -r, --regions LIST           Region post-filter chr[:beg-end[,...]] (post-filter in v1).
  -R, --regions-file FILE      BED-like regions file.
  -t, --targets LIST           Like -r but always a post-filter.
  -T, --targets-file FILE      BED-like targets file (post-filter).
  -i, --include EXPR           Keep records matching expression.
  -e, --exclude EXPR           Drop records matching expression.
      --threads N              Accepted; v1 is single-threaded.
  -h, --help                   Show this help.
      --version                Show version.

IMPUTE2 HAP/legend conversion modes (implemented):
      --hapsample PREFIX           VCF/BCF -> .hap + .samples.
      --hapsample2vcf PREFIX       .hap + .samples -> VCF/BCF.
      --haplegendsample PREFIX     VCF/BCF -> .hap + .legend + .samples.
      --haplegendsample2vcf PREFIX .hap + .legend + .samples -> VCF/BCF.
      --haploid2diploid            Emit haploid genotypes as diploid homozygotes.

Accepted-and-deferred conversion modes (parse cleanly; v1 emits a
"not implemented" error pointing at docs/PARITY_ROADMAP.md if any
are actually set):
  --gvcf2vcf, --gvcf, -f/--fasta-ref, -g/--gensample,
  -G/--gensample2vcf, --3N6, --tag, --chrom, --keep-duplicates,
  --sex, --vcf-ids, --tsv2vcf, -c/--columns.

Accepted-and-ignored stubs (no-op against the round-trip flow):
  --regions-overlap, --targets-overlap, --no-version,
  -W/--write-index, --verbosity.
`

func runConvert(args []string) int {
	fs := flag.NewFlagSet("bcftools convert", flag.ContinueOnError)
	fs.SetOutput(io.Discard)

	var (
		outputType    string
		outputPath    string
		samples       string
		samplesFile   string
		forceSamples  bool
		regions       string
		regionsFile   string
		targets       string
		targetsFile   string
		includeExpr   string
		excludeExpr   string
		threads       int
		compressLevel int
		showHelp      bool
		showVer       bool
		// Accepted-and-deferred upstream flags (project parity rule:
		// every documented upstream flag must parse cleanly, even
		// when not implemented). Each that gets non-default values
		// triggers a one-line "not implemented" error with a roadmap
		// pointer. Tracked in docs/PARITY_ROADMAP.md.
		regionsOverlap      int
		targetsOverlap      int
		noVersion           bool
		writeIndex          string
		verbosity           int
		gvcf2vcf            bool
		fastaRef            string
		gvcfBlocks          string
		gensample           string
		gensample2vcf       string
		threeN6             bool
		tagFlag             string
		chromFlag           string
		keepDuplicates      bool
		sexFlag             string
		vcfIds              bool
		hapsample           string
		hapsample2vcf       string
		haploid2diploid     bool
		haplegendsample     string
		haplegendsample2vcf string
		tsv2vcf             string
		columnsFlag         string
	)
	cliflag.StringVar(fs, &outputType, "O", "output-type", "v", "Output type")
	cliflag.StringVar(fs, &outputPath, "o", "output", "", "Output path")
	cliflag.StringVar(fs, &samples, "s", "samples", "", "Samples list")
	cliflag.StringVar(fs, &samplesFile, "S", "samples-file", "", "Samples file")
	fs.BoolVar(&forceSamples, "force-samples", false, "Allow missing requested samples")
	cliflag.StringVar(fs, &regions, "r", "regions", "", "Region(s)")
	cliflag.StringVar(fs, &regionsFile, "R", "regions-file", "", "Regions file")
	cliflag.StringVar(fs, &targets, "t", "targets", "", "Targets")
	cliflag.StringVar(fs, &targetsFile, "T", "targets-file", "", "Targets file")
	cliflag.StringVar(fs, &includeExpr, "i", "include", "", "Include expression")
	cliflag.StringVar(fs, &excludeExpr, "e", "exclude", "", "Exclude expression")
	// --- accepted-and-deferred stubs ---
	// Long-form only for upstream short-flag letters that already
	// collide with `-h/--help` / `-v/--version` in our convention.
	fs.IntVar(&regionsOverlap, "regions-overlap", -1, "")
	fs.IntVar(&targetsOverlap, "targets-overlap", -1, "")
	fs.BoolVar(&noVersion, "no-version", false, "")
	cliflag.StringVar(fs, &writeIndex, "W", "write-index", "", "")
	fs.IntVar(&verbosity, "verbosity", 0, "")
	fs.BoolVar(&gvcf2vcf, "gvcf2vcf", false, "")
	cliflag.StringVar(fs, &fastaRef, "f", "fasta-ref", "", "")
	fs.StringVar(&gvcfBlocks, "gvcf", "", "")
	cliflag.StringVar(fs, &gensample, "g", "gensample", "", "")
	cliflag.StringVar(fs, &gensample2vcf, "G", "gensample2vcf", "", "")
	fs.BoolVar(&threeN6, "3N6", false, "")
	fs.StringVar(&tagFlag, "tag", "", "")
	fs.StringVar(&chromFlag, "chrom", "", "")
	fs.BoolVar(&keepDuplicates, "keep-duplicates", false, "")
	fs.StringVar(&sexFlag, "sex", "", "")
	fs.BoolVar(&vcfIds, "vcf-ids", false, "")
	fs.StringVar(&hapsample, "hapsample", "", "")
	fs.StringVar(&hapsample2vcf, "hapsample2vcf", "", "")
	fs.BoolVar(&haploid2diploid, "haploid2diploid", false, "")
	fs.StringVar(&haplegendsample, "haplegendsample", "", "")
	fs.StringVar(&haplegendsample2vcf, "haplegendsample2vcf", "", "")
	fs.StringVar(&tsv2vcf, "tsv2vcf", "", "")
	cliflag.StringVar(fs, &columnsFlag, "c", "columns", "", "")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "Threads (accepted, ignored)")
	cliflag.IntVar(fs, &compressLevel, "l", "compression-level", -1, "gzip level for -O z")
	fs.BoolVar(&showHelp, "h", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "version", false, "")

	if err := fs.Parse(args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, convertUsage)
		return 2
	}
	if showHelp {
		fmt.Print(convertUsage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}

	// Reject conversion modes that aren't implemented in v1 with a
	// clear error pointing at docs/PARITY_ROADMAP.md (rather than
	// silently accepting and producing the wrong output).
	if deferred := checkConvertDeferred(checkConvertDeferredInputs{
		gvcf2vcf:       gvcf2vcf,
		fastaRef:       fastaRef,
		gvcfBlocks:     gvcfBlocks,
		gensample:      gensample,
		gensample2vcf:  gensample2vcf,
		threeN6:        threeN6,
		tagFlag:        tagFlag,
		chromFlag:      chromFlag,
		keepDuplicates: keepDuplicates,
		sexFlag:        sexFlag,
		vcfIds:         vcfIds,
		tsv2vcf:        tsv2vcf,
		columnsFlag:    columnsFlag,
	}); deferred != "" {
		fmt.Fprintf(os.Stderr, "bcftools convert: %s not implemented in v1; tracked in docs/PARITY_ROADMAP.md#bcftools\n", deferred)
		return 2
	}

	// HAP/legend (IMPUTE2) conversion modes. These are mutually
	// exclusive with the round-trip flow and with each other; dispatch
	// to the dedicated entry points in convert_hap.go.
	if hapMode, err := runConvertHap(convertHapInputs{
		hapsample:           hapsample,
		hapsample2vcf:       hapsample2vcf,
		haplegendsample:     haplegendsample,
		haplegendsample2vcf: haplegendsample2vcf,
		haploid2diploid:     haploid2diploid,
		outputType:          outputType,
		outputPath:          outputPath,
		compressLevel:       compressLevel,
		samples:             samples,
		samplesFile:         samplesFile,
		forceSamples:        forceSamples,
		regions:             regions,
		regionsFile:         regionsFile,
		targets:             targets,
		targetsFile:         targetsFile,
		includeExpr:         includeExpr,
		excludeExpr:         excludeExpr,
		rest:                fs.Args(),
	}); hapMode {
		return err
	}
	// Silently-ignored stubs (matching upstream when the flag is a
	// no-op or an indexing optimisation the v1 doesn't need).
	_ = regionsOverlap
	_ = targetsOverlap
	_ = noVersion
	_ = writeIndex
	_ = verbosity

	rest := fs.Args()
	if len(rest) == 0 {
		fmt.Fprintln(os.Stderr, "bcftools convert: missing input file")
		fmt.Fprint(os.Stderr, convertUsage)
		return 2
	}

	format, err := bcftools.ParseOutputFormat(outputType)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}
	opts := bcftools.ConvertOptions{
		OutputFormat:  format,
		CompressLevel: compressLevel,
		ForceSamples:  forceSamples,
		IncludeExpr:   includeExpr,
		ExcludeExpr:   excludeExpr,
		SamplesFile:   samplesFile,
		RegionsFile:   regionsFile,
		TargetsFile:   targetsFile,
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
		fmt.Fprintf(os.Stderr, "bcftools convert: %v\n", err)
		return 1
	}
	defer out.Close()
	if _, err := bcftools.ConvertFile(rest[0], out, opts); err != nil {
		fmt.Fprintf(os.Stderr, "bcftools convert: %v\n", err)
		return 1
	}
	return 0
}

// convertHapInputs groups the flag values consumed by the HAP/legend
// (IMPUTE2) conversion modes of `bcftools convert`.
type convertHapInputs struct {
	hapsample           string
	hapsample2vcf       string
	haplegendsample     string
	haplegendsample2vcf string
	haploid2diploid     bool
	outputType          string
	outputPath          string
	compressLevel       int
	samples             string
	samplesFile         string
	forceSamples        bool
	regions             string
	regionsFile         string
	targets             string
	targetsFile         string
	includeExpr         string
	excludeExpr         string
	rest                []string
}

// runConvertHap dispatches the HAP/legend conversion modes. The first
// return value reports whether one of these modes was actually requested;
// when true, the second value is the process exit code. When false the
// caller falls through to the standard round-trip flow.
func runConvertHap(in convertHapInputs) (bool, int) {
	// Exactly one of the four file-naming modes may be set.
	set := 0
	if in.hapsample != "" {
		set++
	}
	if in.hapsample2vcf != "" {
		set++
	}
	if in.haplegendsample != "" {
		set++
	}
	if in.haplegendsample2vcf != "" {
		set++
	}
	if set == 0 {
		// --haploid2diploid alone is a modifier with no base mode; it is
		// silently ignored by the round-trip flow, matching upstream
		// where it only affects the hap exporters.
		return false, 0
	}
	if set > 1 {
		fmt.Fprintln(os.Stderr, "bcftools convert: only one of --hapsample/--hapsample2vcf/--haplegendsample/--haplegendsample2vcf may be given")
		return true, 2
	}

	hapOpts := bcftools.HapConvertOptions{
		Hap2Dip:      in.haploid2diploid,
		ForceSamples: in.forceSamples,
		SamplesFile:  in.samplesFile,
		RegionsFile:  in.regionsFile,
		TargetsFile:  in.targetsFile,
		IncludeExpr:  in.includeExpr,
		ExcludeExpr:  in.excludeExpr,
	}
	if in.samples != "" {
		hapOpts.Samples = bcftools.SplitCommaList(in.samples)
	}
	if in.regions != "" {
		hapOpts.Regions = bcftools.SplitCommaList(in.regions)
	}
	if in.targets != "" {
		hapOpts.Targets = bcftools.SplitCommaList(in.targets)
	}

	format, err := bcftools.ParseOutputFormat(in.outputType)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return true, 2
	}
	hapOpts.OutputFormat = format
	hapOpts.CompressLevel = in.compressLevel

	// VCF -> hap exporters need a positional input file.
	switch {
	case in.hapsample != "":
		if len(in.rest) == 0 {
			fmt.Fprintln(os.Stderr, "bcftools convert: missing input file")
			return true, 2
		}
		hapOpts.Prefix = in.hapsample
		if _, err := bcftools.VCFToHapSample(in.rest[0], hapOpts, os.Stderr); err != nil {
			fmt.Fprintf(os.Stderr, "bcftools convert: %v\n", err)
			return true, 1
		}
		return true, 0
	case in.haplegendsample != "":
		if len(in.rest) == 0 {
			fmt.Fprintln(os.Stderr, "bcftools convert: missing input file")
			return true, 2
		}
		hapOpts.Prefix = in.haplegendsample
		if _, err := bcftools.VCFToHapLegendSample(in.rest[0], hapOpts, os.Stderr); err != nil {
			fmt.Fprintf(os.Stderr, "bcftools convert: %v\n", err)
			return true, 1
		}
		return true, 0
	case in.hapsample2vcf != "":
		out, err := openOutFile(in.outputPath)
		if err != nil {
			fmt.Fprintf(os.Stderr, "bcftools convert: %v\n", err)
			return true, 1
		}
		defer out.Close()
		if _, err := bcftools.HapSampleToVCF(in.hapsample2vcf, out, hapOpts, os.Stderr); err != nil {
			fmt.Fprintf(os.Stderr, "bcftools convert: %v\n", err)
			return true, 1
		}
		return true, 0
	case in.haplegendsample2vcf != "":
		out, err := openOutFile(in.outputPath)
		if err != nil {
			fmt.Fprintf(os.Stderr, "bcftools convert: %v\n", err)
			return true, 1
		}
		defer out.Close()
		if _, err := bcftools.HapLegendSampleToVCF(in.haplegendsample2vcf, out, hapOpts, os.Stderr); err != nil {
			fmt.Fprintf(os.Stderr, "bcftools convert: %v\n", err)
			return true, 1
		}
		return true, 0
	}
	return false, 0
}

// checkConvertDeferredInputs groups the conversion-mode flag values that
// v1 recognises but does not implement. checkConvertDeferred returns
// the upstream flag name when any are set, or "" if none.
type checkConvertDeferredInputs struct {
	gvcf2vcf       bool
	fastaRef       string
	gvcfBlocks     string
	gensample      string
	gensample2vcf  string
	threeN6        bool
	tagFlag        string
	chromFlag      string
	keepDuplicates bool
	sexFlag        string
	vcfIds         bool
	tsv2vcf        string
	columnsFlag    string
}

func checkConvertDeferred(in checkConvertDeferredInputs) string {
	switch {
	case in.gvcf2vcf:
		return "--gvcf2vcf"
	case in.fastaRef != "":
		return "-f/--fasta-ref"
	case in.gvcfBlocks != "":
		return "--gvcf"
	case in.gensample != "":
		return "-g/--gensample"
	case in.gensample2vcf != "":
		return "-G/--gensample2vcf"
	case in.threeN6:
		return "--3N6"
	case in.tagFlag != "":
		return "--tag"
	case in.chromFlag != "":
		return "--chrom"
	case in.keepDuplicates:
		return "--keep-duplicates"
	case in.sexFlag != "":
		return "--sex"
	case in.vcfIds:
		return "--vcf-ids"
	case in.tsv2vcf != "":
		return "--tsv2vcf"
	case in.columnsFlag != "":
		return "-c/--columns"
	}
	return ""
}

const mendelianUsage = `bcftools mendelian - detect Mendelian-inconsistent genotypes.

Usage:
  bcftools mendelian [options] <in.vcf[.gz]|in.bcf>

Options:
  -t, --trio CHILD,FATHER,MOTHER  Single trio (may be supplied multiple times).
  -T, --trio-file FILE            File of CHILD,FATHER,MOTHER (or CHILD<TAB>FATHER<TAB>MOTHER) lines.
  -c, --count                     Emit a TSV trio-level summary instead of VCF (alias for -m c).
  -d, --delete                    Drop records with at least one Mendel error (alias for -m d).
  -m, --mode {a|c|x|d|+|g}        Output mode:
                                    a  annotate INFO/MERR (default)
                                    c  TSV summary
                                    x  X-chromosome aware (father haploid on chrX)
                                    d  delete records with errors
                                    +  retain everything (annotate-only synonym)
                                    g  good-only (same set as 'd' keeps)
      --rules FILE                Ploidy rules file (accepted; v1 honours only the chrX heuristic).
  -O, --output-type {v|z|u|b}     Output format (ignored under -c).
  -o, --output PATH               Output file (default stdout).
  -l, --compression-level N       gzip level for -O z output.
      --threads N                 Accepted; v1 is single-threaded.
  -h, --help                      Show this help.
      --version                   Show version.
`

func runMendelian(args []string) int {
	fs := flag.NewFlagSet("bcftools mendelian", flag.ContinueOnError)
	fs.SetOutput(io.Discard)

	var (
		trioFlag      multiString
		trioFile      string
		count         bool
		deleteFlag    bool
		modeFlag      string
		rules         string
		outputType    string
		outputPath    string
		compressLevel int
		threads       int
		showHelp      bool
		showVer       bool
	)
	fs.Var(&trioFlag, "t", "Trio CHILD,FATHER,MOTHER (may repeat)")
	fs.Var(&trioFlag, "trio", "Trio CHILD,FATHER,MOTHER (may repeat)")
	cliflag.StringVar(fs, &trioFile, "T", "trio-file", "", "Trio file")
	cliflag.BoolVar(fs, &count, "c", "count", false, "Summary mode")
	cliflag.BoolVar(fs, &deleteFlag, "d", "delete", false, "Delete bad records")
	cliflag.StringVar(fs, &modeFlag, "m", "mode", "", "Output mode (a|c|x|d|+|g)")
	fs.StringVar(&rules, "rules", "", "Ploidy rules file")
	cliflag.StringVar(fs, &outputType, "O", "output-type", "v", "Output type")
	cliflag.StringVar(fs, &outputPath, "o", "output", "", "Output path")
	cliflag.IntVar(fs, &compressLevel, "l", "compression-level", -1, "gzip level for -O z")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "Threads (accepted, ignored)")
	fs.BoolVar(&showHelp, "h", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "version", false, "")

	if err := fs.Parse(args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, mendelianUsage)
		return 2
	}
	if showHelp {
		fmt.Print(mendelianUsage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}

	rest := fs.Args()
	if len(rest) == 0 {
		fmt.Fprintln(os.Stderr, "bcftools mendelian: missing input file")
		fmt.Fprint(os.Stderr, mendelianUsage)
		return 2
	}

	mode, err := bcftools.ParseMendelianMode(modeFlag)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}
	format, err := bcftools.ParseOutputFormat(outputType)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}
	opts := bcftools.MendelianOptions{
		TrioFile:      trioFile,
		Mode:          mode,
		Count:         count,
		Delete:        deleteFlag,
		RulesFile:     rules,
		OutputFormat:  format,
		CompressLevel: compressLevel,
	}
	for _, s := range trioFlag {
		t, err := bcftools.ParseTrioFlag(s)
		if err != nil {
			fmt.Fprintln(os.Stderr, err)
			return 2
		}
		opts.Trios = append(opts.Trios, t)
	}

	out, err := openOutFile(outputPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "bcftools mendelian: %v\n", err)
		return 1
	}
	defer out.Close()
	if _, err := bcftools.MendelianFile(rest[0], out, opts); err != nil {
		fmt.Fprintf(os.Stderr, "bcftools mendelian: %v\n", err)
		return 1
	}
	return 0
}

// multiString implements flag.Value for repeated string flags like
// -t CHILD,FATHER,MOTHER. It accumulates each appearance instead of
// overwriting the prior value (which is what bare StringVar would do).
type multiString []string

func (m *multiString) String() string { return fmt.Sprint([]string(*m)) }

func (m *multiString) Set(v string) error {
	*m = append(*m, v)
	return nil
}
