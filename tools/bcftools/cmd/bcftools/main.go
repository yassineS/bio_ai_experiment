// Command bcftools is a pure-Go reimplementation of selected bcftools
// subcommands. This first slice ships the `view` subcommand together with
// the pkg/bioformats/bcf BCF decoder; other subcommands (query, stats, norm,
// concat, merge) and BCF writing will follow in subsequent PRs.
package main

import (
	"flag"
	"fmt"
	"io"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/tools/bcftools/pkg/bcftools"
)

const version = "0.1.0"

const rootUsage = `bcftools - pure-Go reimplementation.

Usage:
  bcftools <subcommand> [options]

Subcommands:
  view      Print, filter, or convert VCF/BCF records.
  norm      Left-align indels, split/join multiallelics, drop duplicates.
  index     Build a CSI (or .tbi) index for a BCF / VCF.gz file.
  help      Show this help (also via -? on subcommands).
  version   Show version.
`

func main() {
	if len(os.Args) < 2 {
		fmt.Fprint(os.Stderr, rootUsage)
		os.Exit(1)
	}
	switch os.Args[1] {
	case "view":
		os.Exit(runView(os.Args[2:]))
	case "norm":
		os.Exit(runNorm(os.Args[2:]))
	case "index":
		os.Exit(runIndex(os.Args[2:]))
	case "help", "--help":
		fmt.Print(rootUsage)
		return
	case "version", "-v", "--version":
		fmt.Println(version)
		return
	default:
		fmt.Fprintf(os.Stderr, "bcftools: unknown subcommand %q\n", os.Args[1])
		fmt.Fprint(os.Stderr, rootUsage)
		os.Exit(1)
	}
}

const indexUsage = `bcftools index - build a CSI (or .tbi) index for a BCF / VCF.gz file.

Usage:
  bcftools index [options] <in.bcf|in.vcf.gz>

Options:
  -c, --csi                Emit CSI (default; required for BCF).
  -t, --tbi                Emit .tbi (VCF.gz only).
      --csi-min-shift N    Use a non-default min_shift for the CSI bin scheme (default 14).
  -o, --output PATH        Output index path (default <in>.csi or <in>.tbi).
  -f, --force              Overwrite an existing index file.
  -@, --threads N          Accepted; v1 is single-threaded.
  -h, --help               Show this help.
      --version            Show version.
`

func runIndex(args []string) int {
	fs := flag.NewFlagSet("bcftools index", flag.ContinueOnError)
	fs.SetOutput(io.Discard)

	var (
		csiFlag     bool
		tbiFlag     bool
		minShift    int
		outPath     string
		force       bool
		threads     int
		showHelp    bool
		showVersion bool
	)
	cliflag.BoolVar(fs, &csiFlag, "c", "csi", false, "Emit CSI")
	cliflag.BoolVar(fs, &tbiFlag, "t", "tbi", false, "Emit TBI (VCF.gz only)")
	fs.IntVar(&minShift, "csi-min-shift", 0, "CSI min_shift")
	cliflag.StringVar(fs, &outPath, "o", "output", "", "Output index path")
	cliflag.BoolVar(fs, &force, "f", "force", false, "Overwrite existing index")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "Threads (accepted, ignored)")
	fs.BoolVar(&showHelp, "h", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVersion, "version", false, "")

	if err := fs.Parse(args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, indexUsage)
		return 2
	}
	if showHelp {
		fmt.Print(indexUsage)
		return 0
	}
	if showVersion {
		fmt.Println(version)
		return 0
	}

	rest := fs.Args()
	if len(rest) == 0 {
		fmt.Fprintln(os.Stderr, "bcftools index: missing input file")
		fmt.Fprint(os.Stderr, indexUsage)
		return 2
	}
	input := rest[0]

	format := bcftools.IndexCSI
	if tbiFlag && !csiFlag {
		format = bcftools.IndexTBI
	}

	opts := bcftools.IndexOptions{
		Format:     format,
		MinShift:   int32(minShift),
		OutputPath: outPath,
		Force:      force,
	}
	written, err := bcftools.BuildIndex(input, opts)
	if err != nil {
		fmt.Fprintf(os.Stderr, "bcftools index: %v\n", err)
		return 1
	}
	_ = written
	return 0
}

const viewUsage = `bcftools view - print, filter, or convert VCF/BCF records.

Usage:
  bcftools view [options] <in.vcf[.gz]|in.bcf> [region ...]

Options:
  -O, --output-type {v|z|u|b}     v=VCF (default), z=VCF.gz, u=uncompressed BCF, b=compressed BCF.
                                  (u and b are NOT YET IMPLEMENTED in this slice.)
  -o, --output PATH               Output file (default stdout).
  -h, --header-only               Emit only the header.
  -H, --no-header                 Drop the header.
  -G, --drop-genotypes            Drop FORMAT and per-sample columns.
  -c, --min-ac N                  Minimum allele count.
  -C, --max-ac N                  Maximum allele count.
  -q, --min-af FLOAT              Minimum allele frequency.
  -Q, --max-af FLOAT              Maximum allele frequency.
  -i, --include EXPR              Keep records matching expression.
  -e, --exclude EXPR              Drop records matching expression.
  -f, --apply-filters NAME[,..]   Keep only PASS or named filters.
  -r, --regions chr[:beg-end[,..]] Region(s) (uses .tbi when available).
  -R, --regions-file PATH         BED-like regions file.
  -t, --targets chr[:beg-end[,..]] Like -r but always a post-filter.
  -T, --targets-file PATH         BED-like targets file (post-filter).
  -s, --samples LIST              Restrict to these samples (comma list).
  -S, --samples-file PATH         File with sample IDs (one per line).
  -l, --compression-level N       gzip level for z output.
      --threads N                 Accepted; v1 is single-threaded.
  -?, --help                      Show this help.
      --version                   Show version.

Note on -h:
  Upstream bcftools uses -h for "header-only" rather than "help" — we follow
  that convention here. Use -? or --help for help.
`

func runView(args []string) int {
	fs := flag.NewFlagSet("bcftools view", flag.ContinueOnError)
	fs.SetOutput(io.Discard) // we print usage ourselves

	var (
		outputType    string
		outputPath    string
		headerOnly    bool
		noHeader      bool
		dropGT        bool
		minAC         int
		maxAC         int
		minAF         float64
		maxAF         float64
		includeExpr   string
		excludeExpr   string
		applyFilters  string
		regions       string
		regionsFile   string
		targets       string
		targetsFile   string
		samples       string
		samplesFile   string
		compressLevel int
		threads       int
		showHelp      bool
		showVer       bool
	)
	cliflag.StringVar(fs, &outputType, "O", "output-type", "v", "Output type")
	cliflag.StringVar(fs, &outputPath, "o", "output", "", "Output path")
	cliflag.BoolVar(fs, &headerOnly, "h", "header-only", false, "Header only")
	cliflag.BoolVar(fs, &noHeader, "H", "no-header", false, "Drop header")
	cliflag.BoolVar(fs, &dropGT, "G", "drop-genotypes", false, "Drop genotypes")
	cliflag.IntVar(fs, &minAC, "c", "min-ac", 0, "Min allele count")
	cliflag.IntVar(fs, &maxAC, "C", "max-ac", 0, "Max allele count")
	cliflag.Float64Var(fs, &minAF, "q", "min-af", 0, "Min allele frequency")
	cliflag.Float64Var(fs, &maxAF, "Q", "max-af", 0, "Max allele frequency")
	cliflag.StringVar(fs, &includeExpr, "i", "include", "", "Include expression")
	cliflag.StringVar(fs, &excludeExpr, "e", "exclude", "", "Exclude expression")
	cliflag.StringVar(fs, &applyFilters, "f", "apply-filters", "", "Filter list")
	cliflag.StringVar(fs, &regions, "r", "regions", "", "Region(s)")
	cliflag.StringVar(fs, &regionsFile, "R", "regions-file", "", "Regions file")
	cliflag.StringVar(fs, &targets, "t", "targets", "", "Targets (post-filter)")
	cliflag.StringVar(fs, &targetsFile, "T", "targets-file", "", "Targets file")
	cliflag.StringVar(fs, &samples, "s", "samples", "", "Samples")
	cliflag.StringVar(fs, &samplesFile, "S", "samples-file", "", "Samples file")
	cliflag.IntVar(fs, &compressLevel, "l", "compression-level", -1, "gzip level")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "Threads (accepted, ignored)")
	fs.BoolVar(&showHelp, "?", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "version", false, "")

	if err := fs.Parse(args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, viewUsage)
		return 2
	}
	if showHelp {
		fmt.Print(viewUsage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}

	positional := fs.Args()
	if len(positional) == 0 {
		fmt.Fprintln(os.Stderr, "bcftools view: missing input file")
		fmt.Fprint(os.Stderr, viewUsage)
		return 2
	}
	input := positional[0]
	cliRegions := positional[1:]

	format, err := bcftools.ParseOutputFormat(outputType)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}

	opts := bcftools.ViewOptions{
		OutputFormat:   format,
		HeaderOnly:     headerOnly,
		NoHeader:       noHeader,
		DropGenotypes:  dropGT,
		MinAlleleCount: minAC,
		MaxAlleleCount: maxAC,
		MinAlleleFreq:  minAF,
		MaxAlleleFreq:  maxAF,
		IncludeExpr:    includeExpr,
		ExcludeExpr:    excludeExpr,
		ApplyFilters:   bcftools.SplitCommaList(applyFilters),
		CompressLevel:  compressLevel,
	}
	if regions != "" {
		opts.Regions = bcftools.SplitCommaList(regions)
	}
	opts.Regions = append(opts.Regions, cliRegions...)
	if targets != "" {
		opts.Targets = bcftools.SplitCommaList(targets)
	}
	if regionsFile != "" {
		regs, err := bcftools.LoadRegionsFile(regionsFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "bcftools view: %v\n", err)
			return 1
		}
		opts.Regions = append(opts.Regions, regs...)
		opts.RegionsFile = regionsFile
	}
	if targetsFile != "" {
		regs, err := bcftools.LoadRegionsFile(targetsFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "bcftools view: %v\n", err)
			return 1
		}
		opts.Targets = append(opts.Targets, regs...)
		opts.TargetsFile = targetsFile
	}
	if samples != "" {
		opts.Samples = bcftools.SplitCommaList(samples)
	}
	if samplesFile != "" {
		names, err := bcftools.LoadSamplesFile(samplesFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "bcftools view: %v\n", err)
			return 1
		}
		opts.Samples = append(opts.Samples, names...)
		opts.SamplesFile = samplesFile
	}

	out, err := openOutFile(outputPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "bcftools view: %v\n", err)
		return 1
	}
	defer out.Close()

	if _, err := bcftools.ViewFile(input, out, opts, os.Stderr); err != nil {
		fmt.Fprintf(os.Stderr, "bcftools view: %v\n", err)
		return 1
	}
	return 0
}

// openOutFile mirrors the helper used by samtools view: stdout when path is
// empty or "-", otherwise create the file.
func openOutFile(path string) (io.WriteCloser, error) {
	if path == "" || path == "-" {
		return nopCloser{os.Stdout}, nil
	}
	return os.Create(path)
}

type nopCloser struct{ io.Writer }

func (nopCloser) Close() error { return nil }

const normUsage = `bcftools norm - left-align indels and normalize multiallelics.

Usage:
  bcftools norm [options] <in.vcf[.gz]|in.bcf>

Options:
  -f, --fasta-ref FASTA          Reference FASTA for left-alignment / REF check.
      --check-ref {e|w|s}        Action on REF/FASTA mismatch: e=error (default), w=warn, s=skip.
  -m, --multiallelics MODE       Split (-) or join (+) multiallelics. MODE = {-snps|-indels|-both|-any|+snps|+indels|+both|+any}.
  -d, --rm-dup MODE              Drop duplicates: snps|indels|both|all|none|exact.
  -a, --atomize                  Decompose complex variants into single-base events.
  -N, --do-not-normalize         Skip left-alignment (useful with -m alone).
  -s, --strict-filter            Apply -f filter list BEFORE splitting (default: after).
  -r, --regions chr[:beg-end]    Region(s) (post-filter on streaming input).
  -R, --regions-file PATH        BED-like regions file.
  -t, --targets chr[:beg-end]    Like -r but always a post-filter.
  -T, --targets-file PATH        BED-like targets file (post-filter).
  -f, --apply-filters NAMES      Keep only PASS / named filters.
  -O, --output-type {v|z|u|b}    Output format (b/u requires BCF writer).
  -o, --output PATH              Output file (default stdout).
  -l, --compression-level N      gzip level for -O z output.
      --threads N                Accepted; v1 is single-threaded.
  -h, --help                     Show this help.
      --version                  Show version.

Note:
  -f is overloaded by upstream bcftools for both "fasta-ref" and "apply-filters".
  We follow the same convention: when --fasta-ref is also set we accept --apply-filters
  via the long form for clarity. The short -f always means --fasta-ref to match upstream.
`

func runNorm(args []string) int {
	fs := flag.NewFlagSet("bcftools norm", flag.ContinueOnError)
	fs.SetOutput(io.Discard)

	var (
		fastaRef       string
		checkRef       string
		multiallelics  string
		rmDup          string
		atomize        bool
		doNotNormalize bool
		strictFilter   bool
		regions        string
		regionsFile    string
		targets        string
		targetsFile    string
		applyFilters   string
		outputType     string
		outputPath     string
		compressLevel  int
		threads        int
		showHelp       bool
		showVersion    bool
	)
	cliflag.StringVar(fs, &fastaRef, "f", "fasta-ref", "", "Reference FASTA")
	fs.StringVar(&checkRef, "check-ref", "e", "REF mismatch policy: e|w|s")
	cliflag.StringVar(fs, &multiallelics, "m", "multiallelics", "", "Split / join multiallelics")
	cliflag.StringVar(fs, &rmDup, "d", "rm-dup", "none", "Drop duplicate records")
	cliflag.BoolVar(fs, &atomize, "a", "atomize", false, "Atomize complex variants")
	cliflag.BoolVar(fs, &doNotNormalize, "N", "do-not-normalize", false, "Skip left-alignment")
	cliflag.BoolVar(fs, &strictFilter, "s", "strict-filter", false, "Apply -f filters before split")
	cliflag.StringVar(fs, &regions, "r", "regions", "", "Region(s)")
	cliflag.StringVar(fs, &regionsFile, "R", "regions-file", "", "Regions file")
	cliflag.StringVar(fs, &targets, "t", "targets", "", "Targets")
	cliflag.StringVar(fs, &targetsFile, "T", "targets-file", "", "Targets file")
	fs.StringVar(&applyFilters, "apply-filters", "", "Filter list (PASS,...)")
	cliflag.StringVar(fs, &outputType, "O", "output-type", "v", "Output type")
	cliflag.StringVar(fs, &outputPath, "o", "output", "", "Output path")
	cliflag.IntVar(fs, &compressLevel, "l", "compression-level", -1, "gzip level")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "Threads (accepted, ignored)")
	fs.BoolVar(&showHelp, "h", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVersion, "version", false, "")

	if err := fs.Parse(args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, normUsage)
		return 2
	}
	if showHelp {
		fmt.Print(normUsage)
		return 0
	}
	if showVersion {
		fmt.Println(version)
		return 0
	}

	rest := fs.Args()
	if len(rest) == 0 {
		fmt.Fprintln(os.Stderr, "bcftools norm: missing input file")
		fmt.Fprint(os.Stderr, normUsage)
		return 2
	}
	input := rest[0]

	checkRefMode, err := bcftools.ParseCheckRefMode(checkRef)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}
	mMode, err := bcftools.ParseMultiallelicMode(multiallelics)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}
	dupMode, err := bcftools.ParseRmDupMode(rmDup)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}
	format, err := bcftools.ParseOutputFormat(outputType)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}

	opts := bcftools.NormOptions{
		FastaRef:       fastaRef,
		CheckRef:       checkRefMode,
		Multiallelics:  mMode,
		RmDup:          dupMode,
		Atomize:        atomize,
		DoNotNormalize: doNotNormalize,
		StrictFilter:   strictFilter,
		ApplyFilters:   bcftools.SplitCommaList(applyFilters),
		OutputFormat:   format,
		CompressLevel:  compressLevel,
	}
	if regions != "" {
		opts.Regions = bcftools.SplitCommaList(regions)
	}
	if regionsFile != "" {
		regs, err := bcftools.LoadRegionsFile(regionsFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "bcftools norm: %v\n", err)
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
			fmt.Fprintf(os.Stderr, "bcftools norm: %v\n", err)
			return 1
		}
		opts.Targets = append(opts.Targets, regs...)
		opts.TargetsFile = targetsFile
	}

	out, err := openOutFile(outputPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "bcftools norm: %v\n", err)
		return 1
	}
	defer out.Close()

	if _, err := bcftools.NormFile(input, out, opts, os.Stderr); err != nil {
		fmt.Fprintf(os.Stderr, "bcftools norm: %v\n", err)
		return 1
	}
	return 0
}
