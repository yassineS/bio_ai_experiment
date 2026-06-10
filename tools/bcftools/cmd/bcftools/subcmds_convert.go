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

GEN/sample conversion modes (Oxford .gen+.sample):
  -g, --gensample PREFIX|GEN,SAMPLE     VCF/BCF -> .gen+.sample.
  -G, --gensample2vcf PREFIX|GEN,SAMPLE .gen+.sample -> VCF/BCF.
      --tag {GT|PL|GP}                  FORMAT tag driving .gen probabilities (default GT).
      --3N6                             3*N+6 column .gen layout (leading CHROM column).
      --sex FILE                        Add a sex column to the .sample file ("ID<TAB>[MF]").
      --vcf-ids                         Use VCF IDs in the .gen ID column.
      --chrom                           Deprecated; errors (use --3N6).

Accepted-and-deferred conversion modes (parse cleanly; v1 emits a
"not implemented" error pointing at docs/PARITY_ROADMAP.md if any
are actually set):
  --gvcf2vcf, --gvcf, -f/--fasta-ref, --keep-duplicates,
  --hapsample, --hapsample2vcf, --haploid2diploid,
  --haplegendsample, --haplegendsample2vcf, --tsv2vcf, -c/--columns.

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

	// --chrom is deprecated upstream and simply errors, pointing the user at
	// --3N6. Mirror that behaviour exactly (including the exit status).
	if chromFlag != "" {
		fmt.Fprintln(os.Stderr, "The --chrom option has been deprecated, please use --3N6 instead")
		return 1
	}

	format, err := bcftools.ParseOutputFormat(outputType)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}

	// --tsv2vcf and --gvcf2vcf are implemented modes; dispatch them before
	// the deferred-mode gate (which still rejects gen/hap export).
	if tsv2vcf != "" {
		out, err := openOutFile(outputPath)
		if err != nil {
			fmt.Fprintf(os.Stderr, "bcftools convert: %v\n", err)
			return 1
		}
		defer out.Close()
		opts := bcftools.TSV2VCFOptions{
			FastaRef:       fastaRef,
			Columns:        columnsFlag,
			SamplesFile:    samplesFile,
			KeepDuplicates: keepDuplicates,
			OutputFormat:   format,
			CompressLevel:  compressLevel,
			NoVersion:      noVersion,
		}
		if samples != "" {
			opts.Samples = bcftools.SplitCommaList(samples)
		}
		if _, err := bcftools.TSV2VCFFile(tsv2vcf, out, opts); err != nil {
			fmt.Fprintf(os.Stderr, "bcftools convert: %v\n", err)
			return 1
		}
		return 0
	}
	if gvcf2vcf {
		rest := fs.Args()
		if len(rest) == 0 {
			fmt.Fprintln(os.Stderr, "bcftools convert: missing input file")
			fmt.Fprint(os.Stderr, convertUsage)
			return 2
		}
		out, err := openOutFile(outputPath)
		if err != nil {
			fmt.Fprintf(os.Stderr, "bcftools convert: %v\n", err)
			return 1
		}
		defer out.Close()
		opts := bcftools.GVCFToVCFOptions{
			FastaRef:      fastaRef,
			IncludeExpr:   includeExpr,
			ExcludeExpr:   excludeExpr,
			OutputFormat:  format,
			CompressLevel: compressLevel,
			NoVersion:     noVersion,
		}
		if _, err := bcftools.GVCFToVCFFile(rest[0], out, opts); err != nil {
			fmt.Fprintf(os.Stderr, "bcftools convert: %v\n", err)
			return 1
		}
		return 0
	}

	// Reject conversion modes that aren't implemented in v1 with a
	// clear error pointing at docs/PARITY_ROADMAP.md (rather than
	// silently accepting and producing the wrong output). The GEN/sample
	// family (-g/-G/--tag/--3N6/--sex/--vcf-ids) is handled below.
	if deferred := checkConvertDeferred(checkConvertDeferredInputs{
		fastaRef:            fastaRef,
		gvcfBlocks:          gvcfBlocks,
		keepDuplicates:      keepDuplicates,
		hapsample:           hapsample,
		hapsample2vcf:       hapsample2vcf,
		haploid2diploid:     haploid2diploid,
		haplegendsample:     haplegendsample,
		haplegendsample2vcf: haplegendsample2vcf,
	}); deferred != "" {
		fmt.Fprintf(os.Stderr, "bcftools convert: %s not implemented in v1; tracked in docs/PARITY_ROADMAP.md#bcftools\n", deferred)
		return 2
	}

	// GEN/sample conversion modes (Oxford .gen+.sample).
	if gensample != "" || gensample2vcf != "" {
		return runConvertGenSample(convertGenInputs{
			gensample:     gensample,
			gensample2vcf: gensample2vcf,
			tag:           tagFlag,
			threeN6:       threeN6,
			vcfIds:        vcfIds,
			sexFile:       sexFlag,
			keepDup:       keepDuplicates,
			noVersion:     noVersion,
			includeExpr:   includeExpr,
			excludeExpr:   excludeExpr,
			outputType:    outputType,
			compressLevel: compressLevel,
			outputPath:    outputPath,
			rest:          fs.Args(),
		})
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

// checkConvertDeferredInputs groups the conversion-mode flag values that
// v1 recognises but does not implement. checkConvertDeferred returns
// the upstream flag name when any are set, or "" if none.
type checkConvertDeferredInputs struct {
	fastaRef            string
	gvcfBlocks          string
	keepDuplicates      bool
	hapsample           string
	hapsample2vcf       string
	haploid2diploid     bool
	haplegendsample     string
	haplegendsample2vcf string
}

func checkConvertDeferred(in checkConvertDeferredInputs) string {
	switch {
	case in.fastaRef != "":
		// Bare -f/--fasta-ref (without --tsv2vcf or --gvcf2vcf, both of
		// which are dispatched earlier) implies an unimplemented gen/hap
		// export mode owned by a sibling PR.
		return "-f/--fasta-ref"
	case in.gvcfBlocks != "":
		return "--gvcf"
	case in.keepDuplicates:
		return "--keep-duplicates"
	case in.hapsample != "":
		return "--hapsample"
	case in.hapsample2vcf != "":
		return "--hapsample2vcf"
	case in.haploid2diploid:
		return "--haploid2diploid"
	case in.haplegendsample != "":
		return "--haplegendsample"
	case in.haplegendsample2vcf != "":
		return "--haplegendsample2vcf"
	}
	return ""
}

// convertGenInputs groups the parsed flag values that drive the GEN/sample
// conversion modes (-g/--gensample and -G/--gensample2vcf).
type convertGenInputs struct {
	gensample     string
	gensample2vcf string
	tag           string
	threeN6       bool
	vcfIds        bool
	sexFile       string
	keepDup       bool
	noVersion     bool
	includeExpr   string
	excludeExpr   string
	outputType    string
	compressLevel int
	outputPath    string
	rest          []string
}

// runConvertGenSample dispatches the VCF<->GEN/sample conversion. It returns
// the process exit code.
func runConvertGenSample(in convertGenInputs) int {
	opts := bcftools.GenSampleOptions{
		Tag:            in.tag,
		ThreeN6:        in.threeN6,
		VCFIDs:         in.vcfIds,
		SexFile:        in.sexFile,
		KeepDuplicates: in.keepDup,
		NoVersion:      in.noVersion,
		IncludeExpr:    in.includeExpr,
		ExcludeExpr:    in.excludeExpr,
		CompressLevel:  in.compressLevel,
	}

	// -G/--gensample2vcf: GEN+sample -> VCF/BCF.
	if in.gensample2vcf != "" {
		format, err := bcftools.ParseOutputFormat(in.outputType)
		if err != nil {
			fmt.Fprintln(os.Stderr, err)
			return 2
		}
		opts.OutputFormat = format
		out, err := openOutFile(in.outputPath)
		if err != nil {
			fmt.Fprintf(os.Stderr, "bcftools convert: %v\n", err)
			return 1
		}
		defer out.Close()
		if err := bcftools.GenSampleToVCFFile(in.gensample2vcf, out, opts); err != nil {
			fmt.Fprintf(os.Stderr, "bcftools convert: %v\n", err)
			return 1
		}
		return 0
	}

	// -g/--gensample: VCF/BCF -> GEN+sample.
	if len(in.rest) == 0 {
		fmt.Fprintln(os.Stderr, "bcftools convert: missing input file")
		return 2
	}
	if err := bcftools.VCFToGenSampleFile(in.rest[0], in.gensample, opts); err != nil {
		fmt.Fprintf(os.Stderr, "bcftools convert: %v\n", err)
		return 1
	}
	return 0
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
