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
		Stderr:       os.Stderr,
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
		// Upstream's "Expected the -p option" message is printed raw
		// (no program prefix). Preserve byte-equality on stderr for it
		// — every other error still carries the "bcftools isec: " tag.
		if err.Error() == "Expected the -p option" {
			fmt.Fprintln(os.Stderr, err.Error())
		} else {
			fmt.Fprintf(os.Stderr, "bcftools isec: %v\n", err)
		}
		return 1
	}
	return 0
}

const sortUsage = `bcftools sort - sort VCF/BCF by (CHROM, POS).

Usage:
  bcftools sort [options] <in.vcf[.gz]|in.bcf>

Options:
  -m, --max-mem MEM          Memory budget for the sort (default 768M, accepted but in-memory in v1).
  -T, --temp-dir DIR         Tmpdir for the external-merge step (accepted, unused in v1).
  -O, --output-type {v|z|u|b}  Output format.
  -o, --output PATH          Output file (default stdout).
  -l, --compression-level N  gzip level for -O z output.
      --threads N            Accepted; v1 is single-threaded.
  -v, --verbosity INT        Verbosity level (accepted).
  -W, --write-index[=FMT]    Auto-index the output (csi|tbi) (accepted).
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
		verbosity     int
		writeIndex    string
		showHelp      bool
		showVer       bool
	)
	cliflag.StringVar(fs, &maxMem, "m", "max-mem", "768M", "Max RAM (accepted, no-op)")
	// -T accepts both upstream's `--temp-dir` and our legacy `--tmpdir`.
	cliflag.StringVar(fs, &tmpDir, "T", "temp-dir", "", "Tmpdir (accepted, no-op)")
	fs.StringVar(&tmpDir, "tmpdir", "", "Alias for --temp-dir (legacy)")
	cliflag.StringVar(fs, &outputType, "O", "output-type", "v", "Output type")
	cliflag.StringVar(fs, &outputPath, "o", "output", "", "Output path")
	cliflag.IntVar(fs, &compressLevel, "l", "compression-level", -1, "gzip level")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "Threads (accepted, ignored)")
	cliflag.IntVar(fs, &verbosity, "v", "verbosity", 0, "Verbosity level")
	cliflag.StringVar(fs, &writeIndex, "W", "write-index", "", "Auto-index output")
	fs.BoolVar(&showHelp, "?", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "version", false, "")

	args = preprocessOptionalArg(args, "-W", "csi")
	args = preprocessOptionalArg(args, "--write-index", "csi")

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
  -h, --headers INT      Display INT header lines [all].
  -n, --records INT      Display INT variant record lines [none].
  -s, --samples INT      Display INT records starting with the #CHROM line.
  -v, --verbosity INT    Verbosity level (accepted).
  -?, --help             Show this help.
      --version          Show version.
`

func runHead(args []string) int {
	fs := flag.NewFlagSet("bcftools head", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	var (
		numLines       int
		numRecords     int
		samplesRecords int
		verbosity      int
		showHelp       bool
		showVer        bool
	)
	cliflag.IntVar(fs, &numLines, "h", "headers", 0, "Number of header lines to print")
	cliflag.IntVar(fs, &numRecords, "n", "records", 0, "Number of variant records to print after the header")
	cliflag.IntVar(fs, &samplesRecords, "s", "samples", 0, "Print INT records starting with #CHROM line")
	cliflag.IntVar(fs, &verbosity, "v", "verbosity", 0, "Verbosity level")
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
		NumLines:       numLines,
		NumRecords:     numRecords,
		SamplesRecords: samplesRecords,
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
  -n, --samples-list LIST    New sample names as a comma-separated list.
  -N, --samples-file FILE    File with new sample names (one per line, OR
                             tab-separated old<TAB>new mapping). Upstream
                             alias -s / --samples is also accepted.
  -f, --fai FILE             Rebuild ##contig lines from a samtools FAI.
  -T, --temp-prefix PATH     Accepted; tmp-file template is unused in v1.
  -O, --output-type {v|z|u|b}  Output format.
  -o, --output PATH          Output file (default stdout).
  -l, --compression-level N  gzip level for -O z output.
      --threads N            Accepted; v1 is single-threaded.
  -v, --verbosity INT        Verbosity level (accepted).
  -?, --help                 Show this help.
      --version              Show version.
`

func runReheader(args []string) int {
	fs := flag.NewFlagSet("bcftools reheader", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	var (
		headerFile    string
		samplesFile   string
		samplesList   string
		faiFile       string
		tempPrefix    string
		outputType    string
		outputPath    string
		compressLevel int
		threads       int
		verbosity     int
		showHelp      bool
		showVer       bool
	)
	cliflag.StringVar(fs, &headerFile, "h", "header", "", "Replacement header file")
	cliflag.StringVar(fs, &samplesFile, "N", "samples-file", "", "New sample names file")
	cliflag.StringVar(fs, &samplesList, "n", "samples-list", "", "New sample names (comma list)")
	// `-s/--samples` is upstream's legacy alias for `-N/--samples-file`.
	fs.StringVar(&samplesFile, "s", "", "Alias for -N/--samples-file (legacy)")
	fs.StringVar(&samplesFile, "samples", "", "Alias for --samples-file (legacy)")
	cliflag.StringVar(fs, &faiFile, "f", "fai", "", "FAI for ##contig rebuild")
	cliflag.StringVar(fs, &tempPrefix, "T", "temp-prefix", "", "Tmp-file template (accepted, no-op)")
	cliflag.StringVar(fs, &outputType, "O", "output-type", "v", "Output type")
	cliflag.StringVar(fs, &outputPath, "o", "output", "", "Output path")
	cliflag.IntVar(fs, &compressLevel, "l", "compression-level", -1, "gzip level")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "Threads (accepted, ignored)")
	cliflag.IntVar(fs, &verbosity, "v", "verbosity", 0, "Verbosity level")
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
	var samplesListSlice []string
	if samplesList != "" {
		samplesListSlice = bcftools.SplitCommaList(samplesList)
	}
	if _, err := bcftools.ReheaderFile(rest[0], out, bcftools.ReheaderOptions{
		HeaderFile:    headerFile,
		SamplesFile:   samplesFile,
		SamplesList:   samplesListSlice,
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
  -C, --columns-file FILE    Read column names from FILE (one per row).
  -e, --exclude EXPR         Exclude sites for which the expression is true.
      --force                Continue past parsing errors (best effort).
  -H, --header-line STR      Literal ##... line appended to the header (repeatable).
  -h, --header-lines FILE    Inject these ##... lines into the output header.
  -i, --include EXPR         Include only sites for which the expression is true.
  -I, --set-id [+]FORMAT     Set ID via a bcftools query-like format string.
  -k, --keep-sites           Keep -i/-e-excluded sites in the output unmodified.
  -l, --merge-logic TAG:TYPE Merge logic for overlapping annotations (deferred).
  -m, --mark-sites [+-]TAG   Tag matched (+) or unmatched (-) sites with INFO/TAG.
      --min-overlap ANN:VCF  Required fractional overlap (deferred).
      --no-version           Do not stamp the provenance line on the header.
      --pair-logic STR       Variant-matching mode (deferred).
  -r, --regions LIST         Region post-filter chr[:beg-end[,...]].
  -R, --regions-file FILE    BED-like sidecar listing regions.
      --regions-overlap N    0|1|2 region inclusion rule (accepted; default 1).
      --rename-annots FILE   Rename annotations (deferred).
      --rename-chrs FILE     Two-column tab file (OLD<TAB>NEW) renaming CHROM.
  -s, --samples [^]LIST      Restrict (or exclude with ^) per-sample annotations.
  -S, --samples-file [^]FILE Same as -s but read from FILE.
      --single-overlaps      Single-overlap mode (deferred).
  -x, --remove FIELD,...     Drop fields from the records (INFO/TAG, FILTER, ID).
  -O, --output-type {v|z|u|b}  Output format.
  -o, --output PATH          Output file (default stdout).
      --threads N            Accepted; v1 is single-threaded.
  -v, --verbosity INT        Verbosity level (accepted, ignored).
  -W, --write-index[=FMT]    Automatically index the output (csi|tbi).
  -?, --help                 Show this help.
      --version              Show version.
`

// repeatedString is a flag.Value that accumulates each occurrence of a
// repeatable string flag (e.g. `-H ##line1 -H ##line2`).
type repeatedString struct{ values *[]string }

func (r repeatedString) String() string {
	if r.values == nil || len(*r.values) == 0 {
		return ""
	}
	return strings.Join(*r.values, ",")
}

func (r repeatedString) Set(s string) error {
	*r.values = append(*r.values, s)
	return nil
}

func runAnnotate(args []string) int {
	fs := flag.NewFlagSet("bcftools annotate", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	var (
		annsFile       string
		columns        string
		columnsFile    string
		headerLines    string
		headerLineList []string
		remove         string
		regions        string
		regionsFile    string
		regionsOverlap int
		renameChrs     string
		renameAnnots   string
		outputType     string
		outputPath     string
		compressLevel  int
		threads        int
		setID          string
		includeExpr    string
		excludeExpr    string
		keepSites      bool
		markSites      string
		mergeLogic     string
		minOverlap     string
		pairLogic      string
		singleOverlaps bool
		samples        string
		samplesFile    string
		force          bool
		writeIndex     string
		verbosity      int
		showHelp       bool
		showVer        bool
	)
	cliflag.StringVar(fs, &annsFile, "a", "annotations", "", "Annotation source")
	cliflag.StringVar(fs, &columns, "c", "columns", "", "Column mapping")
	cliflag.StringVar(fs, &columnsFile, "C", "columns-file", "", "Columns file")
	cliflag.StringVar(fs, &headerLines, "h", "header-lines", "", "Header lines file")
	// -H/--header-line is repeatable.
	fs.Var(repeatedString{&headerLineList}, "H", "Header line (repeatable)")
	fs.Var(repeatedString{&headerLineList}, "header-line", "Header line (repeatable)")
	cliflag.StringVar(fs, &remove, "x", "remove", "", "Fields to drop")
	cliflag.StringVar(fs, &regions, "r", "regions", "", "Region(s)")
	cliflag.StringVar(fs, &regionsFile, "R", "regions-file", "", "Regions file")
	fs.IntVar(&regionsOverlap, "regions-overlap", 1, "Region inclusion rule (0|1|2)")
	cliflag.StringVar(fs, &setID, "I", "set-id", "", "Set ID column using a query-like format string")
	cliflag.StringVar(fs, &includeExpr, "i", "include", "", "Include expression")
	cliflag.StringVar(fs, &excludeExpr, "e", "exclude", "", "Exclude expression")
	cliflag.BoolVar(fs, &keepSites, "k", "keep-sites", false, "Keep -i/-e-excluded sites unchanged")
	cliflag.StringVar(fs, &markSites, "m", "mark-sites", "", "Mark matched/unmatched sites with INFO/TAG")
	cliflag.StringVar(fs, &mergeLogic, "l", "merge-logic", "", "Merge logic (deferred)")
	fs.StringVar(&minOverlap, "min-overlap", "", "Min overlap (deferred)")
	fs.StringVar(&pairLogic, "pair-logic", "", "Pair logic (deferred)")
	fs.BoolVar(&singleOverlaps, "single-overlaps", false, "Single overlaps mode (deferred)")
	fs.StringVar(&renameAnnots, "rename-annots", "", "Rename annotations (deferred)")
	cliflag.StringVar(fs, &samples, "s", "samples", "", "Samples (^ to exclude)")
	cliflag.StringVar(fs, &samplesFile, "S", "samples-file", "", "Samples file (^ to exclude)")
	fs.BoolVar(&force, "force", false, "Continue past parsing errors")
	fs.StringVar(&renameChrs, "rename-chrs", "", "Rename CHROM via two-col map")
	cliflag.StringVar(fs, &outputType, "O", "output-type", "v", "Output type")
	cliflag.StringVar(fs, &outputPath, "o", "output", "", "Output path")
	fs.IntVar(&compressLevel, "compression-level", -1, "gzip level")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "Threads (accepted, ignored)")
	cliflag.StringVar(fs, &writeIndex, "W", "write-index", "", "Automatically index output (csi|tbi)")
	cliflag.IntVar(fs, &verbosity, "v", "verbosity", 0, "Verbosity level")
	fs.BoolVar(&showHelp, "?", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "version", false, "")
	registerNoVersionIfAbsent(fs)

	// -W/--write-index accepts a bare form; expand to the upstream default.
	args = preprocessOptionalArg(args, "-W", "csi")
	args = preprocessOptionalArg(args, "--write-index", "csi")

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
	// Deferred-but-accepted flags: rejection-parity diagnostics so the
	// user sees a clear signal rather than silently-wrong output.
	if mergeLogic != "" {
		fmt.Fprintln(os.Stderr, "bcftools annotate: --merge-logic is to be implemented, please open an issue on github")
		return 1
	}
	if minOverlap != "" {
		fmt.Fprintln(os.Stderr, "bcftools annotate: --min-overlap is to be implemented, please open an issue on github")
		return 1
	}
	if pairLogic != "" {
		fmt.Fprintln(os.Stderr, "bcftools annotate: --pair-logic is to be implemented, please open an issue on github")
		return 1
	}
	if singleOverlaps {
		fmt.Fprintln(os.Stderr, "bcftools annotate: --single-overlaps is to be implemented, please open an issue on github")
		return 1
	}
	if renameAnnots != "" {
		fmt.Fprintln(os.Stderr, "bcftools annotate: --rename-annots is to be implemented, please open an issue on github")
		return 1
	}
	if regionsOverlap < 0 || regionsOverlap > 2 {
		fmt.Fprintf(os.Stderr, "bcftools annotate: --regions-overlap must be 0, 1 or 2 (got %d)\n", regionsOverlap)
		return 2
	}
	if writeIndex != "" && writeIndex != "csi" && writeIndex != "tbi" {
		fmt.Fprintf(os.Stderr, "bcftools annotate: --write-index must be csi or tbi (got %q)\n", writeIndex)
		return 2
	}
	format, err := bcftools.ParseOutputFormat(outputType)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}

	// Parse `^`-prefixed sample lists.
	samplesExclude := false
	if strings.HasPrefix(samples, "^") {
		samples = samples[1:]
		samplesExclude = true
	}
	sf := samplesFile
	if strings.HasPrefix(sf, "^") {
		sf = sf[1:]
		samplesExclude = true
	}

	noVersionFlag := fs.Lookup("no-version")
	noVersion := noVersionFlag != nil && noVersionFlag.Value.String() == "true"

	opts := bcftools.AnnotateOptions{
		Annotations:    annsFile,
		Columns:        columns,
		ColumnsFile:    columnsFile,
		HeaderLines:    headerLines,
		HeaderLine:     headerLineList,
		Remove:         remove,
		RegionsFile:    regionsFile,
		RegionsOverlap: regionsOverlap,
		RenameChromMap: renameChrs,
		OutputFormat:   format,
		CompressLevel:  compressLevel,
		SetID:          setID,
		IncludeExpr:    includeExpr,
		ExcludeExpr:    excludeExpr,
		KeepSites:      keepSites,
		MarkSites:      markSites,
		SamplesExclude: samplesExclude,
		Force:          force,
		NoVersion:      noVersion,
		PGCommand:      strings.Join(os.Args, " "),
		WriteIndex:     writeIndex,
		Verbosity:      verbosity,
	}
	if regions != "" {
		opts.Regions = bcftools.SplitCommaList(regions)
	}
	if samples != "" {
		opts.Samples = bcftools.SplitCommaList(samples)
	}
	if sf != "" {
		names, err := bcftools.LoadSamplesFile(sf)
		if err != nil {
			fmt.Fprintf(os.Stderr, "bcftools annotate: %v\n", err)
			return 1
		}
		opts.Samples = append(opts.Samples, names...)
		opts.SamplesFile = sf
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
