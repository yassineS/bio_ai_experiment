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

// commaAccumFlag is a flag.Value that joins every occurrence of a repeated
// flag with commas, mirroring how upstream bcftools accumulates repeated
// `-l/--merge-logic` arguments into one comma-separated string.
type commaAccumFlag struct {
	value *string
}

func (f commaAccumFlag) String() string {
	if f.value == nil {
		return ""
	}
	return *f.value
}

func (f commaAccumFlag) Set(s string) error {
	if *f.value != "" {
		*f.value += ","
	}
	*f.value += s
	return nil
}

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
  -@, --threads N            Worker threads for parallel BGZF compression of
                             z/b output (>1 enables it).
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

	if err := fs.Parse(args); err != nil {
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
		Threads:       threads,
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
  -@, --threads N            Worker threads for parallel BGZF compression of
                             -O z / -O b output (default 0 = serial).
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
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "Worker threads for parallel BGZF compression")
	fs.BoolVar(&showHelp, "?", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVersion, "version", false, "")

	if err := fs.Parse(args); err != nil {
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
		Threads:      threads,
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
  -@, --threads N            Worker threads for parallel BGZF compression of
                             z/b output (>1 enables it).
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

	if err := fs.Parse(args); err != nil {
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
		Threads:       threads,
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
  -s, --samples              Print one sample-name per line and exit.
  -?, --help                 Show this help.
      --version              Show version.
`

func runHead(args []string) int {
	fs := flag.NewFlagSet("bcftools head", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	var (
		numLines    int
		samplesOnly bool
		showHelp    bool
		showVer     bool
	)
	cliflag.IntVar(fs, &numLines, "h", "headers", 0, "Number of header lines to print")
	cliflag.BoolVar(fs, &samplesOnly, "s", "samples", false, "Print sample names only")
	fs.BoolVar(&showHelp, "?", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "version", false, "")

	if err := fs.Parse(args); err != nil {
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
  -@, --threads N            Worker threads for parallel BGZF compression of
                             z/b output (>1 enables it).
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
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "Worker threads for parallel BGZF compression")
	fs.BoolVar(&showHelp, "?", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "version", false, "")

	if err := fs.Parse(args); err != nil {
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
		Threads:       threads,
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
      --rename-annots FILE   Rename annotations: TYPE/old<TAB>new (TYPE=INFO|FORMAT|FILTER).
  -I, --set-id [+]FORMAT     Set the ID column from a query-like macro string.
                             Macros: %CHROM,%POS,%REF,%ALT,%FIRST_ALT,%QUAL,
                             %FILTER,%TYPE,%END,%INFO/TAG. A leading + only sets
                             missing IDs.
      --merge-logic TAG:TYPE Multi-overlap merge logic for range tables
                             (first|append|append-missing|unique|sum|avg|min|max).
      --min-overlap ANN:VCF  Minimum overlap fraction (annotation:VCF, 0-1).
      --pair-logic STR       Allele pairing for VCF sources
                             (exact|some|all|snps|indels|both|id) [some].
      --single-overlaps      Apply only the first overlapping annotation row.
  -O, --output-type {v|z|u|b}  Output format.
  -o, --output PATH          Output file (default stdout).
  -l, --compression-level N  gzip level for -O z output.
  -@, --threads N            Worker threads for parallel BGZF compression of
                             z/b output (>1 enables it).
  -?, --help                 Show this help.
      --version              Show version.
`

func runAnnotate(args []string) int {
	fs := flag.NewFlagSet("bcftools annotate", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	var (
		annsFile       string
		columns        string
		headerLines    string
		remove         string
		regions        string
		renameChrs     string
		renameAnnots   string
		setID          string
		mergeLogic     string
		minOverlap     string
		pairLogic      string
		singleOverlaps bool
		outputType     string
		outputPath     string
		compressLevel  int
		threads        int
		showHelp       bool
		showVer        bool
	)
	cliflag.StringVar(fs, &annsFile, "a", "annotations", "", "Annotation source")
	cliflag.StringVar(fs, &columns, "c", "columns", "", "Column mapping")
	cliflag.StringVar(fs, &headerLines, "H", "header-lines", "", "Header lines file")
	cliflag.StringVar(fs, &remove, "x", "remove", "", "Fields to drop")
	cliflag.StringVar(fs, &regions, "r", "regions", "", "Region(s)")
	fs.StringVar(&renameChrs, "rename-chrs", "", "Rename CHROM via two-col map")
	fs.StringVar(&renameAnnots, "rename-annots", "", "Rename INFO/FORMAT/FILTER tags via map file")
	cliflag.StringVar(fs, &setID, "I", "set-id", "", "Set ID column from a macro string")
	fs.Var(commaAccumFlag{&mergeLogic}, "merge-logic", "Multi-overlap merge logic TAG:TYPE")
	fs.StringVar(&minOverlap, "min-overlap", "", "Minimum overlap fraction ANN:VCF")
	fs.StringVar(&pairLogic, "pair-logic", "", "Allele pairing mode for VCF sources")
	fs.BoolVar(&singleOverlaps, "single-overlaps", false, "Apply only the first overlapping row")
	cliflag.StringVar(fs, &outputType, "O", "output-type", "v", "Output type")
	cliflag.StringVar(fs, &outputPath, "o", "output", "", "Output path")
	cliflag.IntVar(fs, &compressLevel, "l", "compression-level", -1, "gzip level")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "Threads (accepted, ignored)")
	fs.BoolVar(&showHelp, "?", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "version", false, "")

	if err := fs.Parse(args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		io.WriteString(os.Stderr, annotateUsage)
		return 2
	}
	if showHelp {
		io.WriteString(os.Stdout, annotateUsage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}
	rest := fs.Args()
	if len(rest) == 0 {
		fmt.Fprintln(os.Stderr, "bcftools annotate: missing input file")
		io.WriteString(os.Stderr, annotateUsage)
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
		RenameAnnots:   renameAnnots,
		SetID:          setID,
		MergeLogic:     mergeLogic,
		MinOverlap:     minOverlap,
		PairLogic:      pairLogic,
		SingleOverlaps: singleOverlaps,
		OutputFormat:   format,
		CompressLevel:  compressLevel,
		Threads:        threads,
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
