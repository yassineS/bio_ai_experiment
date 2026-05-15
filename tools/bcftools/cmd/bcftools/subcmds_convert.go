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

Deferred output paths (tracked in docs/PARITY_ROADMAP.md):
  --gvcf2vcf, --tag2tag, --haplegendsample2vcf, --hapsample2vcf,
  --tsv2vcf, --gensample2vcf, --gvcf, and the PLINK / GEN / HAP
  family. The v1 port implements only the round-trip pass-through
  with sample / region filtering; emit a feature request if you need
  one of the deferred shapes.
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

const mendelianUsage = `bcftools mendelian - detect Mendelian-inconsistent genotypes.

Usage:
  bcftools mendelian [options] <in.vcf[.gz]|in.bcf>

Options:
  -t, --trio CHILD,FATHER,MOTHER  Single trio (may be supplied multiple times).
  -T, --trio-file FILE            File of CHILD,FATHER,MOTHER (or CHILD<TAB>FATHER<TAB>MOTHER) lines.
  -c, --count                     Emit a TSV trio-level summary instead of VCF (alias for -m c).
  -d, --delete                    Drop records with at least one Mendel error (alias for -m d).
  -m, --mode {a|c|x|d|+}          Output mode:
                                    a  annotate INFO/MERR (default)
                                    c  TSV summary
                                    x  X-chromosome aware (father haploid on chrX)
                                    d  delete records with errors
                                    +  retain everything (annotate-only synonym)
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
	cliflag.StringVar(fs, &modeFlag, "m", "mode", "", "Output mode (a|c|x|d|+)")
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
