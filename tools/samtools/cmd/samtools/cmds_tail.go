package main

// New subcommand runners landed in the "tail closure" wave: idxstats,
// quickcheck, dict, cat, reheader, addreplacerg, fixmate, merge, coverage,
// split. These each follow the same pattern as the existing runX
// functions in main.go: build a FlagSet via pkg/cliflag, parse, marshal
// into the library options struct, and call into pkg/samtools.

import (
	"flag"
	"fmt"
	"io"
	"os"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/tools/samtools/pkg/samtools"
)

// ----- idxstats ---------------------------------------------------------

const idxstatsUsage = `samtools idxstats - per-reference read counts.

Usage:
  samtools idxstats <in.sorted.bam>

Reads the BAI index next to the BAM (<in>.bam.bai) and emits one line
per reference: name<TAB>length<TAB>mapped<TAB>unmapped, plus a final
"*<TAB>0<TAB>0<TAB>n_no_coor" row covering unplaced reads.

Options:
  -h, --help     Show this help.
  -v, --version  Show version.
`

func runIdxstats(args []string) int {
	fs := flag.NewFlagSet("samtools idxstats", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	var showHelp, showVer bool
	fs.BoolVar(&showHelp, "h", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "v", false, "")
	fs.BoolVar(&showVer, "version", false, "")
	if err := fs.Parse(args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, idxstatsUsage)
		return 2
	}
	if showHelp {
		fmt.Print(idxstatsUsage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}
	if fs.NArg() == 0 {
		fmt.Fprint(os.Stderr, idxstatsUsage)
		return 2
	}
	if err := samtools.IdxstatsFile(fs.Arg(0), os.Stdout); err != nil {
		fmt.Fprintf(os.Stderr, "samtools idxstats: %v\n", err)
		return 1
	}
	return 0
}

// ----- quickcheck -------------------------------------------------------

const quickcheckUsage = `samtools quickcheck - fast BAM sanity check.

Usage:
  samtools quickcheck [options] <in.bam> [<in2.bam> ...]

Checks each input file for a BGZF leading header, the trailing 28-byte
BGZF EOF block, and a parseable BAM header. Returns 0 on success.

Options:
  -q                   Quiet (suppress per-failure path on stdout).
  -u                   Unmapped expected (allow header without @SQ lines).
  -v, --verbose        Per-file PASS/FAIL output.
  -h, --help           Show this help.
      --version        Show version.
`

func runQuickcheck(args []string) int {
	fs := flag.NewFlagSet("samtools quickcheck", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	var (
		quiet    bool
		unmapped bool
		verbose  bool
		showHelp bool
		showVer  bool
	)
	fs.BoolVar(&quiet, "q", false, "")
	fs.BoolVar(&unmapped, "u", false, "")
	cliflag.BoolVar(fs, &verbose, "v", "verbose", false, "")
	fs.BoolVar(&showHelp, "h", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "version", false, "")
	if err := fs.Parse(args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, quickcheckUsage)
		return 2
	}
	if showHelp {
		fmt.Print(quickcheckUsage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}
	if fs.NArg() == 0 {
		fmt.Fprint(os.Stderr, quickcheckUsage)
		return 2
	}
	out := io.Writer(os.Stdout)
	if quiet {
		out = io.Discard
	}
	_, err := samtools.Quickcheck(fs.Args(), samtools.QuickcheckOptions{
		Verbose:          verbose,
		UnmappedExpected: unmapped,
	}, out)
	if err != nil {
		return 1
	}
	return 0
}

// ----- dict -------------------------------------------------------------

const dictUsage = `samtools dict - emit a sequence dictionary from FASTA.

Usage:
  samtools dict [options] <ref.fa>

Options:
  -a, --assembly NAME    Populate the @SQ AS: tag.
  -s, --species SPECIES  Populate the @SQ SP: tag.
  -u, --uri URI          Populate the @SQ UR: tag (default file://<path>).
  -A, --alias            Emit AN: alias tags from header tokens.
  -H                     Suppress the @HD header line.
  -o, --output PATH      Output path (default stdout).
  -h, --help             Show this help.
  -v, --version          Show version.
`

func runDict(args []string) int {
	fs := flag.NewFlagSet("samtools dict", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	var (
		assembly string
		species  string
		uri      string
		alias    bool
		noHdr    bool
		outPath  string
		showHelp bool
		showVer  bool
	)
	cliflag.StringVar(fs, &assembly, "a", "assembly", "", "")
	cliflag.StringVar(fs, &species, "s", "species", "", "")
	cliflag.StringVar(fs, &uri, "u", "uri", "", "")
	cliflag.BoolVar(fs, &alias, "A", "alias", false, "")
	fs.BoolVar(&noHdr, "H", false, "")
	cliflag.StringVar(fs, &outPath, "o", "output", "", "")
	fs.BoolVar(&showHelp, "h", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "v", false, "")
	fs.BoolVar(&showVer, "version", false, "")

	if err := fs.Parse(args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, dictUsage)
		return 2
	}
	if showHelp {
		fmt.Print(dictUsage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}
	if fs.NArg() == 0 {
		fmt.Fprint(os.Stderr, dictUsage)
		return 2
	}
	out, err := openOut(outPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "samtools dict: %v\n", err)
		return 1
	}
	defer out.Close()
	if err := samtools.DictFile(fs.Arg(0), out, samtools.DictOptions{
		Assembly:        assembly,
		Species:         species,
		URI:             uri,
		AliasFromHeader: alias,
		NoHeader:        noHdr,
	}); err != nil {
		fmt.Fprintf(os.Stderr, "samtools dict: %v\n", err)
		return 1
	}
	return 0
}

// ----- cat --------------------------------------------------------------

const catUsage = `samtools cat - concatenate BAMs (no re-sort).

Usage:
  samtools cat [options] <in1.bam> [<in2.bam> ...]

Options:
  -h FILE             Header override (SAM text).
  -o, --output PATH   Output BAM path (default stdout).
  -b FILE-LIST        File of BAM paths, one per line.
      --threads N     Accepted; ignored.
      --help          Show this help.
      --version       Show version.
`

func runCat(args []string) int {
	fs := flag.NewFlagSet("samtools cat", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	var (
		hdrPath  string
		outPath  string
		fofn     string
		threads  int
		showHelp bool
		showVer  bool
	)
	fs.StringVar(&hdrPath, "h", "", "")
	cliflag.StringVar(fs, &outPath, "o", "output", "", "")
	fs.StringVar(&fofn, "b", "", "")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "version", false, "")

	if err := fs.Parse(args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, catUsage)
		return 2
	}
	if showHelp {
		fmt.Print(catUsage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}
	paths := append([]string{}, fs.Args()...)
	if fofn != "" {
		extra, err := samtools.LoadFOFN(fofn)
		if err != nil {
			fmt.Fprintf(os.Stderr, "samtools cat: %v\n", err)
			return 1
		}
		paths = append(paths, extra...)
	}
	if len(paths) == 0 {
		fmt.Fprint(os.Stderr, catUsage)
		return 2
	}
	out, err := openOut(outPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "samtools cat: %v\n", err)
		return 1
	}
	defer out.Close()
	if err := samtools.CatFiles(paths, out, samtools.CatOptions{
		HeaderOverride: hdrPath,
		Threads:        threads,
	}); err != nil {
		fmt.Fprintf(os.Stderr, "samtools cat: %v\n", err)
		return 1
	}
	return 0
}

// ----- reheader ---------------------------------------------------------

const reheaderUsage = `samtools reheader - replace BAM header.

Usage:
  samtools reheader [options] <new-header.sam> <in.bam>
  samtools reheader [options] -i <in.bam>

Options:
  -i           In-place rewrite (atomic rename).
  -c CMD       Filter the existing header through "sh -c CMD" instead of
               replacing wholesale.
  -o PATH      Output path (default stdout when not -i).
      --no-PG  Accepted; v1 never injects @PG.
  -h, --help   Show this help.
      --version
               Show version.
`

func runReheader(args []string) int {
	fs := flag.NewFlagSet("samtools reheader", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	var (
		inPlace  bool
		cmd      string
		outPath  string
		noPG     bool
		showHelp bool
		showVer  bool
	)
	fs.BoolVar(&inPlace, "i", false, "")
	fs.StringVar(&cmd, "c", "", "")
	fs.StringVar(&outPath, "o", "", "")
	fs.BoolVar(&noPG, "no-PG", false, "")
	fs.BoolVar(&showHelp, "h", false, "")
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

	positional := fs.Args()
	opts := samtools.ReheaderOptions{Command: cmd, InPlace: inPlace, NoPG: noPG}
	var bamPath string
	switch {
	case cmd != "" && len(positional) == 1:
		bamPath = positional[0]
	case cmd != "" && inPlace && len(positional) == 1:
		bamPath = positional[0]
	case len(positional) == 2:
		opts.HeaderPath = positional[0]
		bamPath = positional[1]
	case inPlace && len(positional) == 2:
		opts.HeaderPath = positional[0]
		bamPath = positional[1]
	default:
		fmt.Fprint(os.Stderr, reheaderUsage)
		return 2
	}

	if inPlace {
		if err := samtools.ReheaderFile(bamPath, "", opts); err != nil {
			fmt.Fprintf(os.Stderr, "samtools reheader: %v\n", err)
			return 1
		}
		return 0
	}
	out, err := openOut(outPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "samtools reheader: %v\n", err)
		return 1
	}
	defer out.Close()
	in, err := os.Open(bamPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "samtools reheader: %v\n", err)
		return 1
	}
	defer in.Close()
	if err := samtools.Reheader(in, out, opts); err != nil {
		fmt.Fprintf(os.Stderr, "samtools reheader: %v\n", err)
		return 1
	}
	return 0
}

// ----- addreplacerg -----------------------------------------------------

const addReplaceRGUsage = `samtools addreplacerg - add/replace @RG line and per-record RG tag.

Usage:
  samtools addreplacerg [options] <in.bam>

Options:
  -r, --rg-line STRING   Full @RG line (e.g. "ID:rgX\tSM:s1").
  -R, --rg-id ID         Existing RG ID from the header to apply.
  -m, --mode MODE        "orphan_only" (default) or "overwrite_all".
  -o, --output PATH      Output BAM (default stdout).
  -w, --no-PG            Accepted; v1 never injects @PG.
  -h, --help             Show this help.
  -v, --version          Show version.
`

func runAddReplaceRG(args []string) int {
	fs := flag.NewFlagSet("samtools addreplacerg", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	var (
		rgLine   string
		rgID     string
		mode     string
		outPath  string
		noPG     bool
		showHelp bool
		showVer  bool
	)
	cliflag.StringVar(fs, &rgLine, "r", "rg-line", "", "")
	cliflag.StringVar(fs, &rgID, "R", "rg-id", "", "")
	cliflag.StringVar(fs, &mode, "m", "mode", "orphan_only", "")
	cliflag.StringVar(fs, &outPath, "o", "output", "", "")
	cliflag.BoolVar(fs, &noPG, "w", "no-PG", false, "")
	fs.BoolVar(&showHelp, "h", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "v", false, "")
	fs.BoolVar(&showVer, "version", false, "")

	if err := fs.Parse(args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, addReplaceRGUsage)
		return 2
	}
	if showHelp {
		fmt.Print(addReplaceRGUsage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}
	if fs.NArg() == 0 {
		fmt.Fprint(os.Stderr, addReplaceRGUsage)
		return 2
	}
	var rgMode samtools.AddReplaceRGMode
	switch strings.ToLower(mode) {
	case "orphan_only", "":
		rgMode = samtools.AddReplaceRGOrphanOnly
	case "overwrite_all":
		rgMode = samtools.AddReplaceRGOverwriteAll
	default:
		fmt.Fprintf(os.Stderr, "samtools addreplacerg: bad -m %q (orphan_only|overwrite_all)\n", mode)
		return 2
	}
	in, err := iohelper.OpenReader(fs.Arg(0))
	if err != nil {
		fmt.Fprintf(os.Stderr, "samtools addreplacerg: %v\n", err)
		return 1
	}
	defer in.Close()
	out, err := openOut(outPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "samtools addreplacerg: %v\n", err)
		return 1
	}
	defer out.Close()
	if err := samtools.AddReplaceRG(in, out, samtools.AddReplaceRGOptions{
		RGLine: rgLine,
		RGID:   rgID,
		Mode:   rgMode,
		NoPG:   noPG,
	}); err != nil {
		fmt.Fprintf(os.Stderr, "samtools addreplacerg: %v\n", err)
		return 1
	}
	return 0
}

// ----- fixmate ----------------------------------------------------------

const fixmateUsage = `samtools fixmate - fill in mate-read fields.

Usage:
  samtools fixmate [options] <in.namesorted.bam> [<out.bam>]

Options:
  -m              Add MQ (mate MAPQ) and ms (mate score) tags.
  -c              Add MC (mate CIGAR) tag.
  -r              Remove unmapped reads (and pairs entirely unmapped).
  -p, --no-PG     Accepted; v1 never injects @PG.
      --threads N Accepted; ignored.
  -h, --help      Show this help.
  -v, --version   Show version.
`

func runFixmate(args []string) int {
	fs := flag.NewFlagSet("samtools fixmate", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	var (
		addMS    bool
		addMC    bool
		rmUnmap  bool
		noPG     bool
		threads  int
		showHelp bool
		showVer  bool
	)
	fs.BoolVar(&addMS, "m", false, "")
	fs.BoolVar(&addMC, "c", false, "")
	fs.BoolVar(&rmUnmap, "r", false, "")
	cliflag.BoolVar(fs, &noPG, "p", "no-PG", false, "")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "")
	fs.BoolVar(&showHelp, "h", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "v", false, "")
	fs.BoolVar(&showVer, "version", false, "")

	if err := fs.Parse(args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, fixmateUsage)
		return 2
	}
	if showHelp {
		fmt.Print(fixmateUsage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}
	if fs.NArg() == 0 {
		fmt.Fprint(os.Stderr, fixmateUsage)
		return 2
	}
	in, err := iohelper.OpenReader(fs.Arg(0))
	if err != nil {
		fmt.Fprintf(os.Stderr, "samtools fixmate: %v\n", err)
		return 1
	}
	defer in.Close()
	outPath := ""
	if fs.NArg() > 1 {
		outPath = fs.Arg(1)
	}
	out, err := openOut(outPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "samtools fixmate: %v\n", err)
		return 1
	}
	defer out.Close()
	if err := samtools.Fixmate(in, out, samtools.FixmateOptions{
		AddMateScore:   addMS,
		AddMateCigar:   addMC,
		RemoveUnmapped: rmUnmap,
		NoPG:           noPG,
		Threads:        threads,
	}); err != nil {
		fmt.Fprintf(os.Stderr, "samtools fixmate: %v\n", err)
		return 1
	}
	return 0
}

// ----- merge ------------------------------------------------------------

const mergeUsage = `samtools merge - merge sorted BAMs.

Usage:
  samtools merge [options] <out.bam> <in1.bam> [<in2.bam> ...]

Options:
  -n                K-way merge of name-sorted inputs (lexicographic).
  -b FILE-LIST      File of input BAM paths, one per line.
  -h FILE           Header override (SAM text).
  -r RG             Force every record's RG tag to this RG-line's ID.
  -c                Collapse identical @PG chains.
  -p                Preserve every @PG line.
  -l N              Output BGZF deflate level (0..9).
  -@, --threads N   Accepted; ignored.
  -h, --help        Show this help.
  -v, --version     Show version.
`

func runMerge(args []string) int {
	fs := flag.NewFlagSet("samtools merge", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	var (
		byName     bool
		fofn       string
		hdrPath    string
		forceRG    string
		collapsePG bool
		preservePG bool
		compLevel  int
		threads    int
		showHelp   bool
		showVer    bool
	)
	fs.BoolVar(&byName, "n", false, "")
	fs.StringVar(&fofn, "b", "", "")
	fs.StringVar(&hdrPath, "header", "", "")
	// -h is overloaded: "header file" (per upstream) AND "help" — we
	// resolve to header on the merge command and treat `--help` as help.
	hdrShort := fs.String("hf", "", "alias for -h header file (rare)")
	_ = hdrShort
	fs.StringVar(&forceRG, "r", "", "")
	fs.BoolVar(&collapsePG, "c", false, "")
	fs.BoolVar(&preservePG, "p", false, "")
	fs.IntVar(&compLevel, "l", -1, "")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "v", false, "")
	fs.BoolVar(&showVer, "version", false, "")

	// Pre-scan args: if the user passes "-h FILE" treat it as header
	// override (matches upstream); --help still goes to help.
	rawArgs := args
	for i := 0; i < len(rawArgs); i++ {
		if rawArgs[i] == "-h" && i+1 < len(rawArgs) {
			hdrPath = rawArgs[i+1]
			rawArgs = append(rawArgs[:i], rawArgs[i+2:]...)
			i--
		} else if rawArgs[i] == "--help" {
			fmt.Print(mergeUsage)
			return 0
		}
	}

	if err := fs.Parse(rawArgs); err != nil {
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

	positional := fs.Args()
	if len(positional) == 0 {
		fmt.Fprint(os.Stderr, mergeUsage)
		return 2
	}
	outPath := positional[0]
	inputs := append([]string{}, positional[1:]...)
	if fofn != "" {
		extra, err := samtools.LoadFOFN(fofn)
		if err != nil {
			fmt.Fprintf(os.Stderr, "samtools merge: %v\n", err)
			return 1
		}
		inputs = append(inputs, extra...)
	}
	if len(inputs) == 0 {
		fmt.Fprintln(os.Stderr, "samtools merge: no input files")
		return 2
	}
	out, err := openOut(outPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "samtools merge: %v\n", err)
		return 1
	}
	defer out.Close()
	if err := samtools.MergeFiles(inputs, out, samtools.MergeOptions{
		ByName:         byName,
		HeaderOverride: hdrPath,
		ForceRGLine:    forceRG,
		CompressLevel:  compLevel,
		CollapsePG:     collapsePG,
		PreservePG:     preservePG,
		Threads:        threads,
	}); err != nil {
		fmt.Fprintf(os.Stderr, "samtools merge: %v\n", err)
		return 1
	}
	return 0
}

// ----- coverage ---------------------------------------------------------

const coverageUsage = `samtools coverage - per-region coverage summary.

Usage:
  samtools coverage [options] <in.bam>

Tabular output (default):
  #rname startpos endpos numreads covbases coverage meandepth baseq mapq

Options:
  -r, --region SPEC      Restrict to region (repeatable).
  -q, --min-mapq N       Skip records with MAPQ below N.
  -Q, --min-baseq N      Skip bases with quality below N.
  -f, --include-flags N  Required flag bits.
  -F, --exclude-flags N  Excluded flag bits (default 0xF04).
  -o, --output PATH      Output path (default stdout).
  -H                     Suppress the column-header line.
  -A                     ASCII-histogram mode (NOT YET IMPLEMENTED).
  -h, --help             Show this help.
  -v, --version          Show version.
`

func runCoverage(args []string) int {
	fs := flag.NewFlagSet("samtools coverage", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	var (
		regions  multiString
		minMAPQ  int
		minBaseQ int
		incFlags int
		excFlags int
		outPath  string
		noHdr    bool
		hist     bool
		showHelp bool
		showVer  bool
	)
	fs.Var(&regions, "r", "")
	fs.Var(&regions, "region", "")
	cliflag.IntVar(fs, &minMAPQ, "q", "min-mapq", 0, "")
	cliflag.IntVar(fs, &minBaseQ, "Q", "min-baseq", 0, "")
	cliflag.IntVar(fs, &incFlags, "f", "include-flags", 0, "")
	cliflag.IntVar(fs, &excFlags, "F", "exclude-flags", 0, "")
	cliflag.StringVar(fs, &outPath, "o", "output", "", "")
	fs.BoolVar(&noHdr, "H", false, "")
	fs.BoolVar(&hist, "A", false, "")
	fs.BoolVar(&showHelp, "h", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "v", false, "")
	fs.BoolVar(&showVer, "version", false, "")

	if err := fs.Parse(args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, coverageUsage)
		return 2
	}
	if showHelp {
		fmt.Print(coverageUsage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}
	if fs.NArg() == 0 {
		fmt.Fprint(os.Stderr, coverageUsage)
		return 2
	}
	if hist {
		fmt.Fprintln(os.Stderr, "samtools coverage: -A (ASCII histogram) is not yet implemented; tabular mode below")
	}
	in, err := iohelper.OpenReader(fs.Arg(0))
	if err != nil {
		fmt.Fprintf(os.Stderr, "samtools coverage: %v\n", err)
		return 1
	}
	defer in.Close()
	out, err := openOut(outPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "samtools coverage: %v\n", err)
		return 1
	}
	defer out.Close()
	if err := samtools.Coverage(in, out, samtools.CoverageOptions{
		Regions:      []string(regions),
		MinMAPQ:      uint8(minMAPQ),
		MinBaseQ:     uint8(minBaseQ),
		IncludeFlags: uint16(incFlags),
		ExcludeFlags: uint16(excFlags),
		Histogram:    hist,
		HeaderOff:    noHdr,
	}); err != nil {
		fmt.Fprintf(os.Stderr, "samtools coverage: %v\n", err)
		return 1
	}
	return 0
}

// ----- split ------------------------------------------------------------

const splitUsage = `samtools split - split BAM by @RG.

Usage:
  samtools split [options] <in.bam>

Options:
  -f PATTERN     Output filename pattern. Tokens: %! (RG ID), %* (input
                 basename), %. (input extension). Default "%*_%!.bam".
  -u FILE        Output for unidentified reads.
  -d, --no-PG    Accepted; v1 never injects @PG.
  -h, --help     Show this help.
  -v, --version  Show version.
`

func runSplit(args []string) int {
	fs := flag.NewFlagSet("samtools split", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	var (
		pattern  string
		unident  string
		noPG     bool
		showHelp bool
		showVer  bool
	)
	fs.StringVar(&pattern, "f", "", "")
	fs.StringVar(&unident, "u", "", "")
	cliflag.BoolVar(fs, &noPG, "d", "no-PG", false, "")
	fs.BoolVar(&showHelp, "h", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "v", false, "")
	fs.BoolVar(&showVer, "version", false, "")

	if err := fs.Parse(args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, splitUsage)
		return 2
	}
	if showHelp {
		fmt.Print(splitUsage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}
	if fs.NArg() == 0 {
		fmt.Fprint(os.Stderr, splitUsage)
		return 2
	}
	if err := samtools.SplitFile(fs.Arg(0), samtools.SplitOptions{
		Pattern:      pattern,
		Unidentified: unident,
		NoPG:         noPG,
	}); err != nil {
		fmt.Fprintf(os.Stderr, "samtools split: %v\n", err)
		return 1
	}
	return 0
}

// ----- markdup ----------------------------------------------------------

const markdupUsage = `samtools markdup - mark / remove PCR duplicates.

Usage:
  samtools markdup [options] <in.bam> <out.bam>

Two-pass algorithm: pass 1 builds per-fragment buckets keyed by
(refID, unclipped coord, mate refID, mate unclipped coord, orientation);
pass 2 re-streams the input and marks all but the highest-scoring record
in each bucket with the 0x400 (duplicate) flag.

Options:
  -r, --remove-dups       Drop duplicates from output (vs just marking).
  -d, --max-dist N        Optical-dup distance. v1 accepts the flag but
                          does NOT implement optical-dup detection.
  -s, --mode {t|s|tp}     Key mode: template (default), sequence, or
                          template+position (folded into template).
  -T, --tmpdir PATH       Accepted; v1 keeps state in memory.
  -l, --max-len N         Max read length considered (default 300).
      --include-flags N   Require ALL bits set.
      --exclude-flags N   Drop records with ANY bit set.
  -c, --clear-tags        Clear pre-existing dup tags (do/dt/mc).
  -t, --add-tag           Write the 'do' aux tag on flagged duplicates.
  -@, --threads N         Accepted; v1 is single-threaded.
  -o, --output FILE       Output path (default stdout).
      --no-PG             Suppress @PG injection (we never inject anyway).
  -h, --help              Show this help.
  -v, --version           Show version.
`

func runMarkdup(args []string) int {
	fs := flag.NewFlagSet("samtools markdup", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	var (
		removeDups bool
		maxDist    int
		modeStr    string
		tmpDir     string
		maxLen     int
		includeF   int
		excludeF   int
		clearTags  bool
		addTag     bool
		threads    int
		outPath    string
		noPG       bool
		showHelp   bool
		showVer    bool
	)
	cliflag.BoolVar(fs, &removeDups, "r", "remove-dups", false, "")
	cliflag.IntVar(fs, &maxDist, "d", "max-dist", 0, "")
	cliflag.StringVar(fs, &modeStr, "s", "mode", "t", "")
	cliflag.StringVar(fs, &tmpDir, "T", "tmpdir", "", "")
	cliflag.IntVar(fs, &maxLen, "l", "max-len", 300, "")
	fs.IntVar(&includeF, "include-flags", 0, "")
	fs.IntVar(&excludeF, "exclude-flags", 0, "")
	cliflag.BoolVar(fs, &clearTags, "c", "clear-tags", false, "")
	cliflag.BoolVar(fs, &addTag, "t", "add-tag", false, "")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "")
	cliflag.StringVar(fs, &outPath, "o", "output", "", "")
	fs.BoolVar(&noPG, "no-PG", false, "")
	fs.BoolVar(&showHelp, "h", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "v", false, "")
	fs.BoolVar(&showVer, "version", false, "")
	if err := fs.Parse(args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, markdupUsage)
		return 2
	}
	if showHelp {
		fmt.Print(markdupUsage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}
	if fs.NArg() == 0 {
		fmt.Fprint(os.Stderr, markdupUsage)
		return 2
	}
	if maxDist != 0 {
		fmt.Fprintln(os.Stderr, "samtools markdup: warning: optical-dup detection (-d) is not yet implemented; PCR dups only")
	}
	var mode samtools.MarkdupMode
	switch modeStr {
	case "t", "":
		mode = samtools.MarkdupModeTemplate
	case "s":
		mode = samtools.MarkdupModeSequence
	case "tp":
		mode = samtools.MarkdupModeTemplatePos
	default:
		fmt.Fprintf(os.Stderr, "samtools markdup: unknown mode %q\n", modeStr)
		return 2
	}
	inPath := fs.Arg(0)
	if fs.NArg() > 1 && outPath == "" {
		outPath = fs.Arg(1)
	}
	opener := func() (io.ReadCloser, error) {
		return iohelper.OpenReader(inPath)
	}
	out, err := openOut(outPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "samtools markdup: %v\n", err)
		return 1
	}
	defer out.Close()
	if _, err := samtools.Markdup(opener, out, samtools.MarkdupOptions{
		RemoveDups:   removeDups,
		MaxDist:      maxDist,
		Mode:         mode,
		TmpDir:       tmpDir,
		MaxLen:       maxLen,
		IncludeFlags: uint16(includeF),
		ExcludeFlags: uint16(excludeF),
		ClearTags:    clearTags,
		AddTag:       addTag,
		Threads:      threads,
		NoPG:         noPG,
	}); err != nil {
		fmt.Fprintf(os.Stderr, "samtools markdup: %v\n", err)
		return 1
	}
	return 0
}

// ----- stats ------------------------------------------------------------

const statsUsage = `samtools stats - exhaustive per-file alignment statistics.

Usage:
  samtools stats [options] <in.bam> [region ...]

Emits an upstream-compatible text report. v1 ships the most-used sections:
SN (Summary Numbers), RL (read-length), MAPQ (MAPQ distribution),
IS (insert sizes), FFQ/LFQ (per-cycle qualities), GCF/GCL (GC-fraction
histograms), GCC (per-cycle GC). Other sections (COV/COV2/GCD/OXC/...)
are skipped with a documented reason; see PARITY_VALIDATION.md.

Options:
  -r, --ref-seq FASTA      Reference FASTA (accepted; sections that need
                           reference bases are skipped without it).
  -c, --coverage MIN[,MAX[,STEP]]
                           Coverage histogram bins (parsed but COV is
                           omitted in v1).
  -l, --required-flag N    Require ALL bits set.
  -F, --filtering-flag N   Drop records with ANY bit set.
  -d, --max-depth N        Cap depth used in coverage.
  -q, --min-mapq N         Skip records with MAPQ < N.
      --remove-dups        Skip duplicate-flagged records.
      --remove-overlaps    Accept (no-op in v1).
  -i, --insert-size N      Max insert size for the IS section (default 8000).
  -x, --sparse             Omit sections that would emit only zero lines.
  -t, --target-regions BED Restrict stats to this BED (skipped in v1).
  -@, --threads N          Accepted; v1 is single-threaded.
  -o, --output FILE        Output path (default stdout).
  -h, --help               Show this help.
  -v, --version            Show version.
`

func runStats(args []string) int {
	fs := flag.NewFlagSet("samtools stats", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	var (
		refSeq     string
		coverage   string
		reqFlag    int
		filtFlag   int
		maxDepth   int
		minMapQ    int
		removeDups bool
		removeOvl  bool
		insertSize int
		sparse     bool
		targetBED  string
		threads    int
		outPath    string
		showHelp   bool
		showVer    bool
	)
	cliflag.StringVar(fs, &refSeq, "r", "ref-seq", "", "")
	cliflag.StringVar(fs, &coverage, "c", "coverage", "", "")
	cliflag.IntVar(fs, &reqFlag, "l", "required-flag", 0, "")
	cliflag.IntVar(fs, &filtFlag, "F", "filtering-flag", 0, "")
	cliflag.IntVar(fs, &maxDepth, "d", "max-depth", 0, "")
	cliflag.IntVar(fs, &minMapQ, "q", "min-mapq", 0, "")
	fs.BoolVar(&removeDups, "remove-dups", false, "")
	fs.BoolVar(&removeOvl, "remove-overlaps", false, "")
	cliflag.IntVar(fs, &insertSize, "i", "insert-size", 8000, "")
	cliflag.BoolVar(fs, &sparse, "x", "sparse", false, "")
	cliflag.StringVar(fs, &targetBED, "t", "target-regions", "", "")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "")
	cliflag.StringVar(fs, &outPath, "o", "output", "", "")
	fs.BoolVar(&showHelp, "h", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "v", false, "")
	fs.BoolVar(&showVer, "version", false, "")
	if err := fs.Parse(args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, statsUsage)
		return 2
	}
	if showHelp {
		fmt.Print(statsUsage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}
	if fs.NArg() == 0 {
		fmt.Fprint(os.Stderr, statsUsage)
		return 2
	}
	in, err := iohelper.OpenReader(fs.Arg(0))
	if err != nil {
		fmt.Fprintf(os.Stderr, "samtools stats: %v\n", err)
		return 1
	}
	defer in.Close()
	out, err := openOut(outPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "samtools stats: %v\n", err)
		return 1
	}
	defer out.Close()
	if err := samtools.Stats(in, out, samtools.StatsOptions{
		RefSeq:         refSeq,
		Coverage:       coverage,
		RequiredFlag:   uint16(reqFlag),
		FilteringFlag:  uint16(filtFlag),
		MaxDepth:       maxDepth,
		MinMAPQ:        uint8(minMapQ),
		RemoveDups:     removeDups,
		RemoveOverlaps: removeOvl,
		MaxInsertSize:  insertSize,
		Sparse:         sparse,
		TargetBED:      targetBED,
		Threads:        threads,
	}); err != nil {
		fmt.Fprintf(os.Stderr, "samtools stats: %v\n", err)
		return 1
	}
	return 0
}

// ----- calmd -----------------------------------------------------------

const calmdUsage = `samtools calmd - compute MD + NM aux tags.

Usage:
  samtools calmd [options] <in.bam|in.sam> <ref.fa>

Walks each record's CIGAR against the reference FASTA, fills in or
updates the MD:Z and NM:i auxiliary tags, and writes the records back
out. Unmapped records pass through unchanged.

Options:
  -e             Replace MATCH bases in SEQ with '='.
  -b             Output BAM (default text SAM).
  -u             Uncompressed BAM out (implies -b).
  -S             Input is SAM (auto-detected — accepted, no-op).
  -A             Accept reads with mismatch/quality issues (no-op in v1).
  -r             Compute BQ tag from BAQ (accepted; BAQ recompute is
                 deferred — see docs/PARITY_ROADMAP.md#samtools).
  -E             Extended BAQ mode (accepted; deferred).
  -Q             Quiet: suppress per-record "different MD/NM" warnings.
  -@, --threads N  Accepted; v1 is single-threaded.
  -o, --output PATH  Output path (default stdout).
  -h, --help     Show this help.
  -v, --version  Show version.
`

func runCalmd(args []string) int {
	fs := flag.NewFlagSet("samtools calmd", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	var (
		useEqual bool
		outBAM   bool
		uncomp   bool
		sInFmt   bool
		adjustA  bool
		realnR   bool
		extBAQ   bool
		quiet    bool
		threads  int
		outPath  string
		showHelp bool
		showVer  bool
	)
	fs.BoolVar(&useEqual, "e", false, "")
	fs.BoolVar(&outBAM, "b", false, "")
	fs.BoolVar(&uncomp, "u", false, "")
	fs.BoolVar(&sInFmt, "S", false, "")
	fs.BoolVar(&adjustA, "A", false, "")
	fs.BoolVar(&realnR, "r", false, "")
	fs.BoolVar(&extBAQ, "E", false, "")
	fs.BoolVar(&quiet, "Q", false, "")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "")
	cliflag.StringVar(fs, &outPath, "o", "output", "", "")
	// Accept-and-ignore stubs for upstream parity (see
	// docs/PARITY_ROADMAP.md#samtools-calmd-deferred). Behaviour is
	// deferred; flag parse must not hard-error.
	var (
		clearMDNM bool
		maxNM     int
		capQ      int
		dropTag   string
		binQual   int
		noPG      bool
		hashQNM   bool
	)
	fs.BoolVar(&clearMDNM, "N", false, "")
	fs.IntVar(&maxNM, "n", 0, "")
	fs.IntVar(&capQ, "C", 0, "")
	fs.StringVar(&dropTag, "d", "", "")
	fs.IntVar(&binQual, "q", 0, "")
	fs.BoolVar(&noPG, "no-PG", false, "")
	fs.BoolVar(&hashQNM, "hash-qnm", false, "")
	fs.BoolVar(&showHelp, "h", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "v", false, "")
	fs.BoolVar(&showVer, "version", false, "")

	if err := fs.Parse(args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, calmdUsage)
		return 2
	}
	if showHelp {
		fmt.Print(calmdUsage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}
	if fs.NArg() < 2 {
		fmt.Fprintln(os.Stderr, "samtools calmd: need <in.bam> <ref.fa>")
		fmt.Fprint(os.Stderr, calmdUsage)
		return 2
	}
	inPath := fs.Arg(0)
	refPath := fs.Arg(1)
	_ = threads // accepted, ignored
	_ = sInFmt  // we auto-detect
	_ = clearMDNM
	_ = maxNM
	_ = capQ
	_ = dropTag
	_ = binQual
	_ = noPG
	_ = hashQNM
	out, err := openOut(outPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "samtools calmd: %v\n", err)
		return 1
	}
	defer out.Close()
	opts := samtools.CalmdOptions{
		UseEqual:     useEqual,
		OutputBAM:    outBAM,
		Uncompressed: uncomp,
		ExtendedBAQ:  extBAQ,
		AdjustCapQ:   adjustA,
		RealignBAQ:   realnR,
		Quiet:        quiet,
	}
	if err := samtools.CalmdFile(inPath, out, refPath, opts, os.Stderr); err != nil {
		fmt.Fprintf(os.Stderr, "samtools calmd: %v\n", err)
		return 1
	}
	return 0
}

// ----- import ----------------------------------------------------------

const importUsage = `samtools import - convert FASTQ to BAM/SAM.

Usage:
  samtools import [options] [file.fastq ...]

Positional arguments (where supplied without -0/-1/-2/-s):
  one file  → single / interleaved input
  two files → R1 then R2

Options:
  -0 FILE             Unpaired reads.
  -1 FILE             Read-1 input (paired output).
  -2 FILE             Read-2 input (paired output).
  -s FILE             Paired input from one file (with /1 /2 in QNAME).
  -r STRING           Full @RG line (or "ID:rgX\tSM:s1" without leading @RG).
  -R STRING           Build @RG with just this ID (shorthand for -r ID:STRING).
  -N, --name2         Keep /1 /2 suffix in QNAME (default strips it).
  -T TAGS             Aux-tag list. "*" = all, "" = none, "BC,QT" = explicit.
  --order TAG         Per-record counter aux. "TAG" = int, "TAG:N" = zero-pad N.
  -o FILE             Output path (default stdout).
  -u                  Uncompressed BAM out.
  -b                  Force BAM output (default; SAM only when -o ends in .sam).
      --no-PG         Accepted; v1 never injects @PG.
  -h, --help          Show this help.
      --version       Show version.
`

func runImport(args []string) int {
	fs := flag.NewFlagSet("samtools import", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	var (
		r0Path   string
		r1Path   string
		r2Path   string
		sPath    string
		rgID     string
		rgLine   string
		auxTags  string
		orderTag string
		outPath  string
		name2    bool
		outBAM   bool
		uncomp   bool
		noPG     bool
		showHelp bool
		showVer  bool
	)
	fs.StringVar(&r0Path, "0", "", "")
	fs.StringVar(&r1Path, "1", "", "")
	fs.StringVar(&r2Path, "2", "", "")
	fs.StringVar(&sPath, "s", "", "")
	fs.StringVar(&rgLine, "r", "", "")
	fs.StringVar(&rgID, "R", "", "")
	fs.StringVar(&auxTags, "T", "", "")
	fs.StringVar(&orderTag, "order", "", "")
	cliflag.StringVar(fs, &outPath, "o", "output", "", "")
	cliflag.BoolVar(fs, &name2, "N", "name2", false, "")
	fs.BoolVar(&outBAM, "b", true, "")
	fs.BoolVar(&uncomp, "u", false, "")
	fs.BoolVar(&noPG, "no-PG", false, "")
	// Accept-and-ignore stubs for upstream parity (see
	// docs/PARITY_ROADMAP.md#samtools-import-deferred).
	var (
		i1Path     string
		i2Path     string
		casavaForm bool
		barcodeTag string
		qualityTag string
		outputFmt  string
		threads    int
	)
	fs.StringVar(&i1Path, "i1", "", "")
	fs.StringVar(&i2Path, "i2", "", "")
	fs.BoolVar(&casavaForm, "i", false, "")
	fs.StringVar(&barcodeTag, "barcode-tag", "", "")
	fs.StringVar(&qualityTag, "quality-tag", "", "")
	cliflag.StringVar(fs, &outputFmt, "O", "output-fmt", "", "")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "")
	fs.BoolVar(&showHelp, "h", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "version", false, "")
	if err := fs.Parse(args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, importUsage)
		return 2
	}
	if showHelp {
		fmt.Print(importUsage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}
	out, err := openOut(outPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "samtools import: %v\n", err)
		return 1
	}
	defer out.Close()
	// outBAM defaults to true so plain `samtools import` produces BAM
	// (upstream behaviour). A `.sam` filename override flips us back to
	// text SAM, matching upstream's mode autodetection.
	if strings.HasSuffix(strings.ToLower(outPath), ".sam") {
		outBAM = false
	}
	opts := samtools.FastqImportOptions{
		SinglePath:      sPath,
		UnpairedPath:    r0Path,
		Read1Path:       r1Path,
		Read2Path:       r2Path,
		ReadGroup:       rgID,
		ReadGroupLine:   rgLine,
		AuxTags:         auxTags,
		OrderTag:        orderTag,
		StripPairSuffix: !name2,
		OutputBAM:       outBAM,
		Uncompressed:    uncomp,
		NoPG:            noPG,
	}
	_ = i1Path
	_ = i2Path
	_ = casavaForm
	_ = barcodeTag
	_ = qualityTag
	_ = outputFmt
	_ = threads
	if _, err := samtools.FastqImportFiles(fs.Args(), out, opts); err != nil {
		fmt.Fprintf(os.Stderr, "samtools import: %v\n", err)
		return 1
	}
	return 0
}

// ----- phase -----------------------------------------------------------

const phaseUsage = `samtools phase - phase haplotypes from heterozygous SNPs.

Usage:
  samtools phase [options] <in.bam|in.sam>

Walks reads against an in-memory pileup, calls heterozygous SNPs
(positions with at least two alleles, each backed by ≥ 2 reads), then
chains each het across overlapping reads using a greedy majority-vote
solver. Emits a 4-column TSV per het site:

  PS<TAB>chrom<TAB>pos<TAB>{0|1|2}

where 0 = ambiguous (no consistent cluster), 1 = hap1, 2 = hap2.
Positions are 1-based to match SAM POS.

Options:
  -k INT       Block-merge window: max number of unphased hets between
               two phased hets before a new block starts. (Default 13.)
  -b STR       Output prefix for per-haplotype BAMs. Accepted; v1 does
               not yet split the BAM (the TSV stream is always emitted
               to -o / stdout). Tracked in PARITY_ROADMAP.md.
  -q INT       Min MAPQ. (Default 13.)
  -Q INT       Min base quality. (Default 13.)
  -D INT       Max depth observed per position. (Default 256.)
  -F           Use the full read (ignore soft-clip trimming). Accepted;
               v1 always uses the aligned slice.
  -A           Mark drop in the dropped output. Accepted; v1 has no
               per-read chimera output.
  -e           Use empirical-Bayes prior. Accepted-and-ignored
               (upstream MCMC fallback is deferred — see roadmap).
  -l INT       Block-merge length cap. Accepted-and-ignored.
  -o, --output PATH  Output TSV path (default stdout).
  -h, --help   Show this help.
      --version  Show version.

MCMC fallback for chimera repair and tied junctions is deferred; see
docs/PARITY_ROADMAP.md.
`

func runPhase(args []string) int {
	fs := flag.NewFlagSet("samtools phase", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	var (
		blockK    int
		outPrefix string
		minMAPQ   int
		minBaseQ  int
		maxDepth  int
		fullRead  bool
		dropAmbig bool
		outPath   string
		showHelp  bool
		showVer   bool
	)
	fs.IntVar(&blockK, "k", samtools.DefaultPhaseBlockWindow, "")
	fs.StringVar(&outPrefix, "b", "", "")
	fs.IntVar(&minMAPQ, "q", samtools.DefaultPhaseMinMAPQ, "")
	fs.IntVar(&minBaseQ, "Q", samtools.DefaultPhaseMinBaseQ, "")
	fs.IntVar(&maxDepth, "D", samtools.DefaultPhaseMaxDepth, "")
	fs.BoolVar(&fullRead, "F", false, "")
	fs.BoolVar(&dropAmbig, "A", false, "")
	// Upstream phase.c:631 declares `-e` (use empirical-Bayes prior) and
	// `-l INT` (block-merge length cap). Both are accepted-and-ignored
	// for CLI parity per docs/PARITY_ROADMAP.md "phase MCMC" deferral.
	var (
		upstreamE bool
		upstreamL int
	)
	fs.BoolVar(&upstreamE, "e", false, "")
	fs.IntVar(&upstreamL, "l", 0, "")
	cliflag.StringVar(fs, &outPath, "o", "output", "", "")
	fs.BoolVar(&showHelp, "h", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "version", false, "")

	if err := fs.Parse(args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, phaseUsage)
		return 2
	}
	if showHelp {
		fmt.Print(phaseUsage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}
	_ = upstreamE
	_ = upstreamL
	if fs.NArg() == 0 {
		fmt.Fprint(os.Stderr, phaseUsage)
		return 2
	}
	in, err := iohelper.OpenReader(fs.Arg(0))
	if err != nil {
		fmt.Fprintf(os.Stderr, "samtools phase: %v\n", err)
		return 1
	}
	defer in.Close()
	out, err := openOut(outPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "samtools phase: %v\n", err)
		return 1
	}
	defer out.Close()
	if _, err := samtools.Phase(in, out, samtools.PhaseOptions{
		BlockWindow:   blockK,
		MinMAPQ:       uint8(minMAPQ),
		MinBaseQ:      uint8(minBaseQ),
		MaxDepth:      maxDepth,
		FullRead:      fullRead,
		DropAmbiguous: dropAmbig,
		OutputPrefix:  outPrefix,
	}); err != nil {
		fmt.Fprintf(os.Stderr, "samtools phase: %v\n", err)
		return 1
	}
	return 0
}

// ----- targetcut -------------------------------------------------------

const targetcutUsage = `samtools targetcut - emit a FASTA slice of each aligned read.

Usage:
  samtools targetcut [options] <in.bam>

For every mapped primary record, writes a FASTA entry containing the
read sequence that aligns to the reference (soft-clipped flanks
removed; insertions retained because they consume query bases;
deletions / refskips / padding contribute nothing).

Options:
  -Q INT             Min base quality. Bases with Phred quality below
                     the cutoff are dropped. (Default 13.)
  -o, --output FILE  Output FASTA (default stdout). Accepts "-".
  -h, --help         Show this help.
      --version      Show version.

Upstream-only flags (accepted-and-ignored — upstream samtools
targetcut is an HMM-consensus tool over fosmid pools; our v1 ships
the simpler aligned-slice-to-FASTA mode). The accepted upstream
flag letters are: -i INT (state-transition penalty), -f FILE (ref
FASTA), -0/-1/-2 INT (HMM emission scores). See
docs/PARITY_ROADMAP.md for the documented scope reduction.
`

func runTargetcut(args []string) int {
	fs := flag.NewFlagSet("samtools targetcut", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	var (
		minBaseQ int
		outPath  string
		// Upstream `cut_target.c` declares these flag letters with
		// different meanings (HMM-consensus tool). We accept them as
		// no-ops so a user copy-pasting an upstream invocation gets a
		// recognised parse rather than a silent flag-letter collision.
		upstreamI  int
		upstreamF  string
		upstreamE0 int
		upstreamE1 int
		upstreamE2 int
		showHelp   bool
		showVer    bool
	)
	fs.IntVar(&minBaseQ, "Q", int(samtools.DefaultTargetcutMinBaseQ), "")
	cliflag.StringVar(fs, &outPath, "o", "output", "", "")
	fs.IntVar(&upstreamI, "i", 0, "")
	fs.StringVar(&upstreamF, "f", "", "")
	fs.IntVar(&upstreamE0, "0", 0, "")
	fs.IntVar(&upstreamE1, "1", 0, "")
	fs.IntVar(&upstreamE2, "2", 0, "")
	fs.BoolVar(&showHelp, "h", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "version", false, "")

	if err := fs.Parse(args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, targetcutUsage)
		return 2
	}
	if showHelp {
		fmt.Print(targetcutUsage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}
	// Discard the upstream-only flag values (recognised for CLI parity).
	_ = upstreamI
	_ = upstreamF
	_ = upstreamE0
	_ = upstreamE1
	_ = upstreamE2
	inPath := "-"
	if fs.NArg() > 0 {
		inPath = fs.Arg(0)
	}
	in, err := iohelper.OpenReader(inPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "samtools targetcut: %v\n", err)
		return 1
	}
	defer in.Close()
	out, err := openOut(outPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "samtools targetcut: %v\n", err)
		return 1
	}
	defer out.Close()
	if _, err := samtools.Targetcut(in, out, samtools.TargetcutOptions{
		MinBaseQ: uint8(minBaseQ),
	}); err != nil {
		fmt.Fprintf(os.Stderr, "samtools targetcut: %v\n", err)
		return 1
	}
	return 0
}

// ----- consensus --------------------------------------------------------

const consensusUsage = `samtools consensus - call a per-position consensus base.

Usage:
  samtools consensus [options] <in.bam>

Options:
  -f, --format FMT          Output format: fasta (default), fastq, pileup.
  -m, --mode MODE           Algorithm: simple (default) or bayesian
                            (accepted but not yet implemented; falls
                            back to simple with a stderr warning).
  -a, --all                 Emit zero-coverage positions as N (FASTA/FASTQ);
                            in pileup mode emit a row for every position
                            in the contig.
  -r, --region REG          Restrict to a "chr[:start-end]" region; may
                            be repeated.
  -l, --positions FILE      BED file of positions to keep.
  -d, --min-depth INT       Minimum unfiltered depth (default 1).
      --min-fraction FLOAT  Minimum proportion of the dominant base over
                            total score (default 0.6).
  -c, --call-fract FLOAT    Minimum proportion (best+second) / total
                            required to make any call (default 0.75).
  -H, --het-fract FLOAT     Minimum proportion of second-best/best to
                            trigger an IUPAC het call (default 0.5);
                            requires -A.
  -A, --ambig               Emit IUPAC ambiguity codes for hets.
  -q, --min-MQ INT          Minimum MAPQ (default 0).
  -Q, --min-BQ INT          Minimum base quality (default 0).
      --show-del yes|no     Render deletions as '*' in FASTA/FASTQ output
                            (default no). Pileup always shows them.
      --show-ins yes|no     Include inserted bases in FASTA/FASTQ output
                            (default yes). Pileup ignores this in v1.
      --line-len INT        Wrap FASTA/FASTQ lines at this many bases
                            (default 70).
      --ignore-overlaps     Accepted; not implemented in v1 (we don't
                            deduplicate mate-pair overlaps).
  -o, --output FILE         Output file (default stdout).
  -@, --threads INT         Accepted; single-threaded in v1.
  -h, --help                Show this help.
  -v, --version             Show version.

Notes:
  - v1 implements the "simple" mode only (majority-vote by qual-weighted
    sum). The "bayesian" mode is accepted on the CLI and falls back to
    simple with a stderr warning; the Gap5-derived posterior caller is
    tracked in docs/PARITY_ROADMAP.md#samtools.
  - In pileup mode each row carries chrom\tpos\tnth\tdepth\tcall\tcq\tseq\tqual,
    matching upstream's basic_pileup format. nth is always 0 in v1
    (insertion-column rows are not yet supported).
`

func runConsensus(args []string) int {
	fs := flag.NewFlagSet("samtools consensus", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	var (
		formatStr string
		modeStr   string
		allPos    bool
		regions   multiString
		bedPath   string
		minDepth  int
		minFrac   float64
		callFract float64
		hetFract  float64
		ambig     bool
		minMAPQ   int
		minBaseQ  int
		showDel   string
		showIns   string
		lineLen   int
		ignoreOvl bool
		outPath   string
		threads   int
		showHelp  bool
		showVer   bool
	)
	cliflag.StringVar(fs, &formatStr, "f", "format", "fasta", "Output format")
	cliflag.StringVar(fs, &modeStr, "m", "mode", "simple", "Consensus mode")
	cliflag.BoolVar(fs, &allPos, "a", "all", false, "Emit zero-coverage positions")
	fs.Var(&regions, "r", "")
	fs.Var(&regions, "region", "")
	cliflag.StringVar(fs, &bedPath, "l", "positions", "", "BED of positions")
	cliflag.IntVar(fs, &minDepth, "d", "min-depth", 1, "Min depth")
	cliflag.Float64Var(fs, &minFrac, "", "min-fraction", 0.6, "Min dominant-base fraction")
	cliflag.Float64Var(fs, &callFract, "c", "call-fract", 0.75, "Min call fraction")
	cliflag.Float64Var(fs, &hetFract, "H", "het-fract", 0.5, "Min het fraction")
	cliflag.BoolVar(fs, &ambig, "A", "ambig", false, "Emit IUPAC ambig codes")
	cliflag.IntVar(fs, &minMAPQ, "q", "min-MQ", 0, "Min MAPQ")
	cliflag.IntVar(fs, &minBaseQ, "Q", "min-BQ", 0, "Min base quality")
	cliflag.StringVar(fs, &showDel, "", "show-del", "no", "Show deletions in FASTA/FASTQ")
	cliflag.StringVar(fs, &showIns, "", "show-ins", "yes", "Include insertions in FASTA/FASTQ")
	cliflag.IntVar(fs, &lineLen, "", "line-len", 70, "FASTA/FASTQ line wrap")
	cliflag.BoolVar(fs, &ignoreOvl, "", "ignore-overlaps", false, "Accepted; not implemented")
	cliflag.StringVar(fs, &outPath, "o", "output", "", "Output path")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "Threads (accepted, ignored)")
	fs.BoolVar(&showHelp, "h", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "v", false, "")
	fs.BoolVar(&showVer, "version", false, "")

	if err := fs.Parse(args); err != nil {
		if err == flag.ErrHelp {
			fmt.Print(consensusUsage)
			return 0
		}
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, consensusUsage)
		return 2
	}
	if showHelp {
		fmt.Print(consensusUsage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}
	if fs.NArg() == 0 {
		fmt.Fprintln(os.Stderr, "samtools consensus: missing input file")
		fmt.Fprint(os.Stderr, consensusUsage)
		return 2
	}
	format, err := samtools.ParseConsensusFormat(formatStr)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}
	mode, err := samtools.ParseConsensusMode(modeStr)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}
	opts := samtools.ConsensusOptions{
		Input:                fs.Arg(0),
		Format:               format,
		Mode:                 mode,
		AllPositions:         allPos,
		Regions:              []string(regions),
		BEDPath:              bedPath,
		MinDepth:             minDepth,
		MinConsensusFraction: minFrac,
		MinCallFraction:      callFract,
		MinHetFraction:       hetFract,
		AmbigCodes:           ambig,
		MinMAPQ:              uint8(minMAPQ),
		MinBaseQ:             uint8(minBaseQ),
		LineLen:              lineLen,
		IgnoreOverlaps:       ignoreOvl,
		Output:               outPath,
		Threads:              threads,
	}
	opts.ShowDel = parseYesNo(showDel, false)
	opts.NoShowIns = !parseYesNo(showIns, true)

	out, err := openOut(outPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "samtools consensus: %v\n", err)
		return 1
	}
	defer out.Close()
	if err := samtools.ConsensusFile(opts, out, os.Stderr); err != nil {
		fmt.Fprintf(os.Stderr, "samtools consensus: %v\n", err)
		return 1
	}
	return 0
}

// parseYesNo accepts upstream samtools' "yes"/"no"/"on"/"off"/"true"/"false"/"1"/"0"
// triplets, returning def for empty.
func parseYesNo(v string, def bool) bool {
	switch strings.ToLower(strings.TrimSpace(v)) {
	case "":
		return def
	case "yes", "y", "on", "true", "1":
		return true
	case "no", "n", "off", "false", "0":
		return false
	}
	return def
}
