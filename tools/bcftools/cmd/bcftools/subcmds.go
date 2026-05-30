// CLI runners for the new bcftools subcommands added in the tail-closure
// PR: `merge`, `isec`, `sort`, `head`, `reheader`, and `annotate`. They
// follow the same parsing-and-dispatch shape as the runners in main.go.
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

const mergeUsage = `bcftools merge - combine multiple per-sample VCF/BCF files.

Usage:
  bcftools merge [options] <in1.vcf[.gz]|in1.bcf> [<in2> ...]

Options:
  -l, --file-list PATH       File of input paths (one per line, # comments).
  -m, --merge MODE           Collapse rule (none|snps|indels|both|all|id, default both).
  -r, --regions LIST         Region list chr[:beg-end[,...]]; v1 post-filter only.
  -R, --regions-file PATH    BED-like regions file.
  -O, --output-type {v|z|u|b}  Output format.
  -o, --output PATH          Output file (default stdout).
  -l, --compression-level N  gzip level for -O z output.
      --threads N            Accepted; v1 is single-threaded.
  -?, --help                 Show this help.
      --version              Show version.
`

func runMerge(args []string) int {
	fs := flag.NewFlagSet("bcftools merge", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	var (
		fileList      string
		mergeFlag     string
		regions       string
		regionsFile   string
		outputType    string
		outputPath    string
		compressLevel int
		threads       int
		showHelp      bool
		showVer       bool
	)
	cliflag.StringVar(fs, &fileList, "L", "file-list", "", "File of input paths")
	cliflag.StringVar(fs, &mergeFlag, "m", "merge", "both", "Collapse rule")
	cliflag.StringVar(fs, &regions, "r", "regions", "", "Region(s)")
	cliflag.StringVar(fs, &regionsFile, "R", "regions-file", "", "Regions file")
	cliflag.StringVar(fs, &outputType, "O", "output-type", "v", "Output type")
	cliflag.StringVar(fs, &outputPath, "o", "output", "", "Output path")
	cliflag.IntVar(fs, &compressLevel, "l", "compression-level", -1, "gzip level")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "Threads (accepted, ignored)")
	fs.BoolVar(&showHelp, "?", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "version", false, "")

	if err := parseFlags(fs, args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, mergeUsage)
		return 2
	}
	if showHelp {
		fmt.Print(mergeUsage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}
	paths := fs.Args()
	if len(paths) == 0 && fileList == "" {
		fmt.Fprintln(os.Stderr, "bcftools merge: missing input files")
		fmt.Fprint(os.Stderr, mergeUsage)
		return 2
	}
	mode, err := bcftools.ParseMergeMode(mergeFlag)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}
	format, err := bcftools.ParseOutputFormat(outputType)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}
	opts := bcftools.MergeOptions{
		FileList:      fileList,
		MergeMode:     mode,
		OutputFormat:  format,
		CompressLevel: compressLevel,
		RegionsFile:   regionsFile,
	}
	if regions != "" {
		opts.Regions = bcftools.SplitCommaList(regions)
	}

	out, err := openOutFile(outputPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "bcftools merge: %v\n", err)
		return 1
	}
	defer out.Close()
	if _, err := bcftools.MergeFiles(paths, out, opts); err != nil {
		fmt.Fprintf(os.Stderr, "bcftools merge: %v\n", err)
		return 1
	}
	return 0
}

const isecUsage = `bcftools isec - set operations on N VCFs.

Usage:
  bcftools isec [options] <in1.vcf[.gz]|in1.bcf> <in2> [<in3> ...]

Options:
  -n, --nfiles SPEC          Membership constraint:
                               +N  present in at least N inputs (default: all).
                               =N  present in exactly N inputs.
                               ~01 inclusion bitmask (one char per input).
  -c, --collapse MODE        Collapse rule (none|snps|indels|both|all|some|id; default none).
  -p, --prefix DIR           Per-input output dir (writes 0000.vcf, 0001.vcf, ... + sites.txt).
  -w, --write LIST           Comma-separated 1-based input indices to write to stdout.
  -O, --output-type {v|z}    Output format.
  -o, --output PATH          Output file (default stdout) — used when -p is not set.
      --threads N            Accepted; v1 is single-threaded.
  -?, --help                 Show this help.
      --version              Show version.
`

func runIsec(args []string) int {
	fs := flag.NewFlagSet("bcftools isec", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	var (
		nfiles      string
		collapse    string
		prefix      string
		writeList   string
		outputType  string
		outputPath  string
		threads     int
		showHelp    bool
		showVersion bool
	)
	cliflag.StringVar(fs, &nfiles, "n", "nfiles", "", "Membership constraint")
	cliflag.StringVar(fs, &collapse, "c", "collapse", "none", "Collapse mode")
	cliflag.StringVar(fs, &prefix, "p", "prefix", "", "Per-input output dir")
	cliflag.StringVar(fs, &writeList, "w", "write", "", "Inputs to dump to stdout")
	cliflag.StringVar(fs, &outputType, "O", "output-type", "v", "Output type")
	cliflag.StringVar(fs, &outputPath, "o", "output", "", "Output path")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "Threads (accepted, ignored)")
	fs.BoolVar(&showHelp, "?", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVersion, "version", false, "")

	if err := parseFlags(fs, args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, isecUsage)
		return 2
	}
	if showHelp {
		fmt.Print(isecUsage)
		return 0
	}
	if showVersion {
		fmt.Println(version)
		return 0
	}
	paths := fs.Args()
	if len(paths) < 2 {
		fmt.Fprintln(os.Stderr, "bcftools isec: need at least two input files")
		fmt.Fprint(os.Stderr, isecUsage)
		return 2
	}
	spec, err := bcftools.ParseNfilesSpec(nfiles)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}
	cmode, err := bcftools.ParseCollapseMode(collapse)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}
	format, err := bcftools.ParseOutputFormat(outputType)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}
	opts := bcftools.IsecOptions{
		Nfiles:       spec,
		Collapse:     cmode,
		Prefix:       prefix,
		OutputFormat: format,
	}
	if writeList != "" {
		for _, p := range strings.Split(writeList, ",") {
			n, err := strconv.Atoi(strings.TrimSpace(p))
			if err != nil {
				fmt.Fprintf(os.Stderr, "bcftools isec: bad -w %q: %v\n", writeList, err)
				return 2
			}
			opts.Write = append(opts.Write, n)
		}
	}

	out, err := openOutFile(outputPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "bcftools isec: %v\n", err)
		return 1
	}
	defer out.Close()
	if _, err := bcftools.IsecFiles(paths, out, opts); err != nil {
		fmt.Fprintf(os.Stderr, "bcftools isec: %v\n", err)
		return 1
	}
	return 0
}

const sortUsage = `bcftools sort - sort VCF/BCF by (CHROM, POS).

Usage:
  bcftools sort [options] <in.vcf[.gz]|in.bcf>

Options:
  -m, --max-mem MEM          Memory budget for the sort (default 768M, accepted but in-memory in v1).
  -T, --tmpdir DIR           Tmpdir for the external-merge step (accepted, unused in v1).
  -O, --output-type {v|z|u|b}  Output format.
  -o, --output PATH          Output file (default stdout).
  -l, --compression-level N  gzip level for -O z output.
      --threads N            Accepted; v1 is single-threaded.
  -?, --help                 Show this help.
      --version              Show version.
`

func runSort(args []string) int {
	fs := flag.NewFlagSet("bcftools sort", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	var (
		maxMem        string
		tmpDir        string
		outputType    string
		outputPath    string
		compressLevel int
		threads       int
		showHelp      bool
		showVer       bool
	)
	cliflag.StringVar(fs, &maxMem, "m", "max-mem", "768M", "Max RAM (accepted, no-op)")
	cliflag.StringVar(fs, &tmpDir, "T", "tmpdir", "", "Tmpdir (accepted, no-op)")
	cliflag.StringVar(fs, &outputType, "O", "output-type", "v", "Output type")
	cliflag.StringVar(fs, &outputPath, "o", "output", "", "Output path")
	cliflag.IntVar(fs, &compressLevel, "l", "compression-level", -1, "gzip level")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "Threads (accepted, ignored)")
	fs.BoolVar(&showHelp, "?", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "version", false, "")

	if err := parseFlags(fs, args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, sortUsage)
		return 2
	}
	if showHelp {
		fmt.Print(sortUsage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}
	rest := fs.Args()
	if len(rest) == 0 {
		fmt.Fprintln(os.Stderr, "bcftools sort: missing input file")
		fmt.Fprint(os.Stderr, sortUsage)
		return 2
	}
	format, err := bcftools.ParseOutputFormat(outputType)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}
	out, err := openOutFile(outputPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "bcftools sort: %v\n", err)
		return 1
	}
	defer out.Close()
	if _, err := bcftools.SortFile(rest[0], out, bcftools.SortOptions{
		OutputFormat:  format,
		CompressLevel: compressLevel,
		MaxMem:        maxMem,
		TmpDir:        tmpDir,
	}); err != nil {
		fmt.Fprintf(os.Stderr, "bcftools sort: %v\n", err)
		return 1
	}
	return 0
}

const headUsage = `bcftools head - print VCF/BCF header.

Usage:
  bcftools head [options] <in.vcf[.gz]|in.bcf>

Options:
  -h, --headers N            Print only the first N header lines.
  -n, --records N            Print the first N variant records after the header.
  -s, --samples              Print one sample-name per line and exit.
  -?, --help                 Show this help.
      --version              Show version.
`

func runHead(args []string) int {
	fs := flag.NewFlagSet("bcftools head", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	var (
		numLines    int
		numRecords  int
		samplesOnly bool
		showHelp    bool
		showVer     bool
	)
	cliflag.IntVar(fs, &numLines, "h", "headers", 0, "Number of header lines to print")
	cliflag.IntVar(fs, &numRecords, "n", "records", 0, "Number of variant records to print after the header")
	cliflag.BoolVar(fs, &samplesOnly, "s", "samples", false, "Print sample names only")
	fs.BoolVar(&showHelp, "?", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "version", false, "")

	if err := parseFlags(fs, args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, headUsage)
		return 2
	}
	if showHelp {
		fmt.Print(headUsage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}
	rest := fs.Args()
	if len(rest) == 0 {
		fmt.Fprintln(os.Stderr, "bcftools head: missing input file")
		fmt.Fprint(os.Stderr, headUsage)
		return 2
	}
	if err := bcftools.HeadFile(rest[0], os.Stdout, bcftools.HeadOptions{
		NumLines:    numLines,
		NumRecords:  numRecords,
		SamplesOnly: samplesOnly,
	}); err != nil {
		fmt.Fprintf(os.Stderr, "bcftools head: %v\n", err)
		return 1
	}
	return 0
}

const reheaderUsage = `bcftools reheader - replace VCF/BCF header in place.

Usage:
  bcftools reheader [options] <in.vcf[.gz]|in.bcf>

Options:
  -h, --header FILE          Replace entire header with the contents of FILE.
  -s, --samples FILE         Sample-rename file (one new name per line, OR
                             tab-separated old<TAB>new mapping).
  -f, --fai FILE             Rebuild ##contig lines from a samtools FAI.
  -O, --output-type {v|z|u|b}  Output format.
  -o, --output PATH          Output file (default stdout).
  -l, --compression-level N  gzip level for -O z output.
      --threads N            Accepted; v1 is single-threaded.
  -?, --help                 Show this help.
      --version              Show version.
`

func runReheader(args []string) int {
	fs := flag.NewFlagSet("bcftools reheader", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	var (
		headerFile    string
		samplesFile   string
		faiFile       string
		outputType    string
		outputPath    string
		compressLevel int
		threads       int
		showHelp      bool
		showVer       bool
	)
	cliflag.StringVar(fs, &headerFile, "H", "header", "", "Replacement header file")
	cliflag.StringVar(fs, &samplesFile, "s", "samples", "", "Sample rename file")
	cliflag.StringVar(fs, &faiFile, "f", "fai", "", "FAI for ##contig rebuild")
	cliflag.StringVar(fs, &outputType, "O", "output-type", "v", "Output type")
	cliflag.StringVar(fs, &outputPath, "o", "output", "", "Output path")
	cliflag.IntVar(fs, &compressLevel, "l", "compression-level", -1, "gzip level")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "Threads (accepted, ignored)")
	fs.BoolVar(&showHelp, "?", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "version", false, "")

	if err := parseFlags(fs, args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, reheaderUsage)
		return 2
	}
	if showHelp {
		fmt.Print(reheaderUsage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}
	rest := fs.Args()
	if len(rest) == 0 {
		fmt.Fprintln(os.Stderr, "bcftools reheader: missing input file")
		fmt.Fprint(os.Stderr, reheaderUsage)
		return 2
	}
	format, err := bcftools.ParseOutputFormat(outputType)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}
	out, err := openOutFile(outputPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "bcftools reheader: %v\n", err)
		return 1
	}
	defer out.Close()
	if _, err := bcftools.ReheaderFile(rest[0], out, bcftools.ReheaderOptions{
		HeaderFile:    headerFile,
		SamplesFile:   samplesFile,
		FaiFile:       faiFile,
		OutputFormat:  format,
		CompressLevel: compressLevel,
	}); err != nil {
		fmt.Fprintf(os.Stderr, "bcftools reheader: %v\n", err)
		return 1
	}
	return 0
}

const annotateUsage = `bcftools annotate - annotate INFO/ID/FILTER from a tab-indexed table.

Usage:
  bcftools annotate [options] <in.vcf[.gz]|in.bcf>

Options:
  -a, --annotations FILE     Annotation source: TAB-delimited (.tab.gz) or VCF.
  -c, --columns LIST         Comma list mapping annotation columns to record
                             fields. Examples: CHROM,POS,REF,ALT,INFO/AC,INFO/AN
  -h, --header-lines FILE    Inject these ##... lines into the output header.
  -x, --remove FIELD,...     Drop fields from the records (INFO/TAG, FILTER, ID).
  -r, --regions LIST         Region post-filter chr[:beg-end[,...]].
      --rename-chrs FILE     Two-column tab file (OLD<TAB>NEW) renaming CHROM.
  -O, --output-type {v|z|u|b}  Output format.
  -o, --output PATH          Output file (default stdout).
  -l, --compression-level N  gzip level for -O z output.
      --threads N            Accepted; v1 is single-threaded.
  -?, --help                 Show this help.
      --version              Show version.
`

func runAnnotate(args []string) int {
	fs := flag.NewFlagSet("bcftools annotate", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	var (
		annsFile      string
		columns       string
		headerLines   string
		remove        string
		regions       string
		renameChrs    string
		outputType    string
		outputPath    string
		compressLevel int
		threads       int
		setID         string
		showHelp      bool
		showVer       bool
	)
	cliflag.StringVar(fs, &annsFile, "a", "annotations", "", "Annotation source")
	cliflag.StringVar(fs, &columns, "c", "columns", "", "Column mapping")
	cliflag.StringVar(fs, &headerLines, "H", "header-lines", "", "Header lines file")
	cliflag.StringVar(fs, &remove, "x", "remove", "", "Fields to drop")
	cliflag.StringVar(fs, &regions, "r", "regions", "", "Region(s)")
	cliflag.StringVar(fs, &setID, "I", "set-id", "", "Set ID column using a query-like format string")
	fs.StringVar(&renameChrs, "rename-chrs", "", "Rename CHROM via two-col map")
	cliflag.StringVar(fs, &outputType, "O", "output-type", "v", "Output type")
	cliflag.StringVar(fs, &outputPath, "o", "output", "", "Output path")
	cliflag.IntVar(fs, &compressLevel, "l", "compression-level", -1, "gzip level")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "Threads (accepted, ignored)")
	fs.BoolVar(&showHelp, "?", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "version", false, "")

	if err := parseFlags(fs, args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, annotateUsage)
		return 2
	}
	if showHelp {
		fmt.Print(annotateUsage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}
	rest := fs.Args()
	if len(rest) == 0 {
		fmt.Fprintln(os.Stderr, "bcftools annotate: missing input file")
		fmt.Fprint(os.Stderr, annotateUsage)
		return 2
	}
	format, err := bcftools.ParseOutputFormat(outputType)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}
	opts := bcftools.AnnotateOptions{
		Annotations:    annsFile,
		Columns:        columns,
		HeaderLines:    headerLines,
		Remove:         remove,
		RegionsFile:    "",
		RenameChromMap: renameChrs,
		OutputFormat:   format,
		CompressLevel:  compressLevel,
		SetID:          setID,
	}
	if regions != "" {
		opts.Regions = bcftools.SplitCommaList(regions)
	}
	out, err := openOutFile(outputPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "bcftools annotate: %v\n", err)
		return 1
	}
	defer out.Close()
	if _, err := bcftools.AnnotateFile(rest[0], out, opts); err != nil {
		fmt.Fprintf(os.Stderr, "bcftools annotate: %v\n", err)
		return 1
	}
	return 0
}

// boolFlag is the optional interface the standard library's flag package
// uses to recognise boolean flags (those that do not consume a following
// value). We use it to distinguish value-taking short flags from boolean
// ones when normalising getopt-style attached values.
type boolFlag interface {
	IsBoolFlag() bool
}

// valueTakingShortFlags inspects fs and returns the set of registered
// single-character flag names that consume a value (i.e. are NOT boolean
// flags). These are the only short flags for which an attached value
// (`-Xvalue`) is meaningful in upstream getopt semantics.
func valueTakingShortFlags(fs *flag.FlagSet) map[byte]bool {
	set := make(map[byte]bool)
	fs.VisitAll(func(f *flag.Flag) {
		if len(f.Name) != 1 {
			return
		}
		if bf, ok := f.Value.(boolFlag); ok && bf.IsBoolFlag() {
			return
		}
		set[f.Name[0]] = true
	})
	return set
}

// normalizeShortFlags rewrites getopt-style attached short-flag values
// into the two-token form that Go's flag package accepts. Upstream
// bcftools is getopt-based and accepts a value attached directly to a
// single-letter flag (e.g. `-Ob`, `norm -m-`, `-m+`, `-mboth`); Go's
// flag package only accepts `-X value` or `-X=value`. For each argument
// of the form `-X...` where X is a registered value-taking short flag
// and extra characters follow X, it splits the token into `-X` and the
// remainder. So `-Ob` -> `-O b`, `-m-` -> `-m -`, `-mboth` -> `-m both`.
//
// It deliberately leaves untouched:
//   - long flags (`--foo`, `--foo=bar`),
//   - boolean short flags (which take no value),
//   - the `-X=value` form (already valid; passed through),
//   - a bare `-` (stdin/stdout),
//   - everything after a bare `--` (end-of-options marker).
func normalizeShortFlags(fs *flag.FlagSet, args []string) []string {
	values := valueTakingShortFlags(fs)
	out := make([]string, 0, len(args)+2)
	for i, a := range args {
		if a == "--" {
			out = append(out, args[i:]...)
			break
		}
		// Candidate: a single-dash flag with at least one char after
		// the flag letter, and not a long flag (`--`) or bare `-`.
		if len(a) > 2 && a[0] == '-' && a[1] != '-' {
			if values[a[1]] && a[2] != '=' {
				out = append(out, a[:2], a[2:])
				continue
			}
		}
		out = append(out, a)
	}
	return out
}

// parseFlags is the shared parse entry point for every bcftools
// subcommand. It (1) ensures `--no-version` is accepted (see
// registerNoVersionIfAbsent) and (2) normalises getopt-style attached
// short-flag values (see normalizeShortFlags) so all subcommands accept
// the upstream attached short-flag forms (`-Ob`, `-m-`, ...) and
// `--no-version` uniformly, before delegating to fs.Parse.
func parseFlags(fs *flag.FlagSet, args []string) error {
	registerNoVersionIfAbsent(fs)
	return fs.Parse(normalizeShortFlags(fs, args))
}

// registerNoVersionIfAbsent registers a no-op `--no-version` boolean flag
// on fs when one is not already present. Upstream bcftools accepts
// `--no-version` as a per-subcommand option that suppresses the
// `##bcftools_*Version`/`##bcftools_*Command` provenance header lines.
// Subcommands that already register `--no-version` (and wire it into
// their options) keep their own registration; this only fills the gap for
// subcommands such as `view` and `norm` which never emit a provenance
// line, so accepting the flag is a safe no-op.
func registerNoVersionIfAbsent(fs *flag.FlagSet) {
	if fs.Lookup("no-version") != nil {
		return
	}
	var noVersion bool
	fs.BoolVar(&noVersion, "no-version", false, "Accepted; no provenance line is emitted.")
}
