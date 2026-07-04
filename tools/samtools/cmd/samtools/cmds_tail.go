package main

// New subcommand runners landed in the "tail closure" wave: idxstats,
// quickcheck, dict, cat, reheader, addreplacerg, fixmate, merge, coverage,
// split. These each follow the same pattern as the existing runX
// functions in main.go: build a FlagSet via pkg/cliflag, parse, marshal
// into the library options struct, and call into pkg/samtools.

import (
	"bytes"
	"flag"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
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
	// Upstream idxstats (bam_index.c getopt "@:X") accepts -@ (threads) and
	// -X (explicit index-file argument). The index path is read single-threaded
	// (it does not decode the BAM body); -@ drives block-parallel BGZF decode for
	// the no-index fallback scan. -X is accepted for compatibility.
	var (
		idxThreads   int
		idxCustomIdx bool
	)
	cliflag.IntVar(fs, &idxThreads, "@", "threads", 0, "BGZF input decode worker count (no-index scan)")
	fs.BoolVar(&idxCustomIdx, "X", false, "")
	if err := cliflag.Parse(fs, args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, idxstatsUsage)
		return 2
	}
	_ = idxCustomIdx
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
	if err := samtools.IdxstatsFile(fs.Arg(0), os.Stdout, idxThreads); err != nil {
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
	// Route through cliflag.Parse so POSIX getopt-style short-flag bundling
	// (-qu == -q -u) works the way upstream quickcheck's getopt parser
	// (bam_quickcheck.c getopt "vqu") accepts it.
	if err := cliflag.Parse(fs, args); err != nil {
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
  -A, --alias            Add an AN: tag of alternate names by adding or
                         removing a leading "chr" prefix (with chrM/chrMT/MT
                         mitochondrial synonyms).
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
	// Upstream dict (dict.c getopt "?AhHa:l:s:u:o:") accepts -l FILE, a file
	// of alternative sequence names. This port does not emit AN: aliases from
	// such a file yet, so -l is an accepted value-taking no-op kept for
	// compatibility (and so bundled clusters parse).
	var altNames string
	cliflag.StringVar(fs, &altNames, "l", "alt-names", "", "Alt-names file (accepted no-op)")

	if err := cliflag.Parse(fs, args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, dictUsage)
		return 2
	}
	_ = altNames
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
	// Upstream cat (bam_cat.c getopt "h:o:b:r:p:qf") accepts a few flags this
	// port does not act on. Register them as accepted stubs so legacy command
	// lines and bundled clusters parse:
	//   -r REG    restrict output to a region (accepted no-op; v1 concatenates
	//             whole files)
	//   -p FILE   add a partial-cat @PG header (accepted no-op)
	//   -q        suppress @PG-merge warnings (accepted no-op)
	//   -f        fast mode: concatenate compressed blocks verbatim (accepted)
	var (
		catRegion string
		catPart   string
		catQuiet  bool
		catFast   bool
	)
	fs.StringVar(&catRegion, "r", "", "")
	fs.StringVar(&catPart, "p", "", "")
	fs.BoolVar(&catQuiet, "q", false, "")
	fs.BoolVar(&catFast, "f", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "version", false, "")

	if err := cliflag.Parse(fs, args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, catUsage)
		return 2
	}
	_ = catRegion
	_ = catPart
	_ = catQuiet
	_ = catFast
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
	// Upstream reheader (bam_reheader.c getopt "hiPc:") spells no-PG as the
	// short flag -P; register it alongside our long --no-PG so bundled
	// clusters (e.g. -iP) parse.
	cliflag.BoolVar(fs, &noPG, "P", "no-PG", false, "")
	fs.BoolVar(&showHelp, "h", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "version", false, "")

	if err := cliflag.Parse(fs, args); err != nil {
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
  -m, --mode MODE        "overwrite_all" (default) or "orphan_only".
  -o, --output PATH      Output BAM (default stdout).
  -w                     Overwrite an existing @RG line with the same ID.
      --no-PG            Accepted; v1 never injects @PG.
  -h, --help             Show this help.
  -v, --version          Show version.
`

func runAddReplaceRG(args []string) int {
	fs := flag.NewFlagSet("samtools addreplacerg", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	var (
		rgLine     string
		rgID       string
		mode       string
		outPath    string
		overwriteW bool
		noPG       bool
		showHelp   bool
		showVer    bool
	)
	cliflag.StringVar(fs, &rgLine, "r", "rg-line", "", "")
	cliflag.StringVar(fs, &rgID, "R", "rg-id", "", "")
	// Upstream default is overwrite_all (bam_addrprg.c retval->mode =
	// overwrite_all), not orphan_only.
	cliflag.StringVar(fs, &mode, "m", "mode", "overwrite_all", "")
	cliflag.StringVar(fs, &outPath, "o", "output", "", "")
	// -w (overwrite existing header @RG) is short-only; --no-PG is long-only.
	fs.BoolVar(&overwriteW, "w", false, "")
	fs.BoolVar(&noPG, "no-PG", false, "")
	fs.BoolVar(&showHelp, "h", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "v", false, "")
	fs.BoolVar(&showVer, "version", false, "")
	// Upstream addreplacerg (bam_addrprg.c getopt "r:R:m:o:O:h@:uw") also
	// accepts -O (output format), -@ (threads), and -u (uncompressed BAM).
	// This port emits BAM single-threaded, so these are accepted no-ops kept
	// for compatibility (and so bundled clusters parse).
	var (
		arOutFmt  string
		arThreads int
		arUncomp  bool
	)
	cliflag.StringVar(fs, &arOutFmt, "O", "output-fmt", "", "Output format (accepted)")
	cliflag.IntVar(fs, &arThreads, "@", "threads", 0, "Threads (accepted, ignored)")
	fs.BoolVar(&arUncomp, "u", false, "")

	if err := cliflag.Parse(fs, args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, addReplaceRGUsage)
		return 2
	}
	_ = arOutFmt
	_ = arThreads
	_ = arUncomp
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
	case "overwrite_all", "":
		rgMode = samtools.AddReplaceRGOverwriteAll
	case "orphan_only":
		rgMode = samtools.AddReplaceRGOrphanOnly
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
		RGLine:            rgLine,
		RGID:              rgID,
		Mode:              rgMode,
		OverwriteHeaderRG: overwriteW,
		NoPG:              noPG,
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
	// Upstream fixmate (bam_mate.c getopt "rpcmMO:@:uz:") exposes a few more
	// flags. Register them as accepted stubs so legacy command lines and
	// bundled clusters parse:
	//   -M        add base-modification (MM/ML) mate handling (accepted no-op)
	//   -O FMT    output format (accepted; this port emits BAM)
	//   -u        uncompressed BAM output (accepted no-op)
	//   -z FLAGS  record-sanitisation flags (accepted no-op)
	var (
		fmBaseMods bool
		fmOutFmt   string
		fmUncomp   bool
		fmSanitize string
	)
	fs.BoolVar(&fmBaseMods, "M", false, "")
	cliflag.StringVar(fs, &fmOutFmt, "O", "output-fmt", "", "Output format (accepted)")
	fs.BoolVar(&fmUncomp, "u", false, "")
	fs.StringVar(&fmSanitize, "z", "", "")

	if err := cliflag.Parse(fs, args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, fixmateUsage)
		return 2
	}
	_ = fmBaseMods
	_ = fmOutFmt
	_ = fmUncomp
	_ = fmSanitize
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
  -c                Combine @RG headers with colliding IDs (default: rename
                    a colliding ID to a distinct "<id>-<8 hex>" form).
  -p                Combine @PG headers with colliding IDs.
  -s SEED           Seed for the @RG-collision rename PRNG (default: wall clock).
  -R STR            Merge only reads overlapping region STR (e.g. chr20:1-2000000);
                    inputs are queried through their sibling .bai/.csi index.
  -l N              Output BGZF deflate level (0..9).
  -@, --threads N   Accepted; ignored.
  -h, --help        Show this help.
  -v, --version     Show version.
`

func runMerge(args []string) int {
	fs := flag.NewFlagSet("samtools merge", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	var (
		byName    bool
		fofn      string
		hdrPath   string
		forceRG   string
		combineRG bool
		combinePG bool
		compLevel int
		threads   int
		showHelp  bool
		showVer   bool
	)
	fs.BoolVar(&byName, "n", false, "")
	fs.StringVar(&fofn, "b", "", "")
	// Upstream merge spells the header override as the value-taking short flag
	// -h FILE (bam_sort.c getopt "...h:..."), so register -h directly as a
	// string flag; --header is our long alias and --help is the separate help
	// switch below. Routing through cliflag.Parse means a bundled value flag
	// like `-h FILE` is handled by getopt's "next token is the value" rule
	// without a manual pre-scan.
	fs.StringVar(&hdrPath, "h", "", "")
	fs.StringVar(&hdrPath, "header", "", "")
	fs.StringVar(&forceRG, "r", "", "")
	fs.BoolVar(&combineRG, "c", false, "")
	fs.BoolVar(&combinePG, "p", false, "")
	fs.IntVar(&compLevel, "l", -1, "")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "v", false, "")
	fs.BoolVar(&showVer, "version", false, "")
	// Upstream merge (bam_sort.c getopt "h:nNru1R:o:f@:l:cps:b:O:t:XL:")
	// exposes several flags this port does not act on. Register them as
	// accepted stubs so legacy command lines and bundled clusters parse:
	//   -N        merge name-sorted with natural (numeric) ordering (accepted;
	//             this port's name merge is lexicographic — see -n)
	//   -u        uncompressed BAM output (accepted no-op)
	//   -1        fastest compression (level 1) (accepted no-op)
	//   -f        force overwrite of the output (accepted; we always write it)
	//   -o FILE   output file (accepted; v1 takes the output as the first
	//             positional, matching the historical merge CLI)
	//   -s SEED   subsample seed (accepted no-op)
	//   -O FMT    output format (accepted; this port emits BAM)
	//   -t TAG    merge by an aux tag's sort order (accepted no-op)
	//   -X        expect explicit index-file arguments (accepted no-op)
	//   -L BED    restrict to a BED of regions (accepted no-op)
	var (
		mgNatural bool
		mgUncomp  bool
		mgLevel1  bool
		mgForce   bool
		mgRegion  string
		mgOutFile string
		mgSeed    int
		mgOutFmt  string
		mgTag     string
		mgCustom  bool
		mgBed     string
	)
	fs.BoolVar(&mgNatural, "N", false, "")
	fs.BoolVar(&mgUncomp, "u", false, "")
	fs.BoolVar(&mgLevel1, "1", false, "")
	fs.BoolVar(&mgForce, "f", false, "")
	fs.StringVar(&mgRegion, "R", "", "")
	fs.StringVar(&mgOutFile, "o", "", "")
	fs.IntVar(&mgSeed, "s", 0, "")
	cliflag.StringVar(fs, &mgOutFmt, "O", "output-fmt", "", "Output format (accepted)")
	fs.StringVar(&mgTag, "t", "", "")
	fs.BoolVar(&mgCustom, "X", false, "")
	fs.StringVar(&mgBed, "L", "", "")

	if err := cliflag.Parse(fs, args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, mergeUsage)
		return 2
	}
	_ = mgNatural
	_ = mgUncomp
	_ = mgLevel1
	_ = mgForce
	_ = mgOutFmt
	_ = mgTag
	_ = mgCustom
	_ = mgBed
	// -s SEED seeds the @RG-collision disambiguation PRNG. Track whether it was
	// given so the merge can match upstream's seeded rename (and otherwise fall
	// back to a wall-clock seed, as upstream does).
	seedSet := false
	fs.Visit(func(f *flag.Flag) {
		if f.Name == "s" {
			seedSet = true
		}
	})
	if showHelp {
		fmt.Print(mergeUsage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}

	positional := fs.Args()
	// Upstream merge accepts the output either as the first positional
	// (legacy form: `merge out.bam in1.bam ...`) or via -o FILE. When -o is
	// given every positional is an input; otherwise the first positional is
	// the output.
	var outPath string
	var inputs []string
	if mgOutFile != "" {
		outPath = mgOutFile
		inputs = append([]string{}, positional...)
	} else {
		if len(positional) == 0 {
			fmt.Fprint(os.Stderr, mergeUsage)
			return 2
		}
		outPath = positional[0]
		inputs = append([]string{}, positional[1:]...)
	}
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
		CombineRG:      combineRG,
		CombinePG:      combinePG,
		RandomSeed:     int64(mgSeed),
		SeedSet:        seedSet,
		Region:         mgRegion,
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
  -l, --min-read-len N   Ignore reads shorter than N bp.
  -q, --min-mapq N       Skip records with MAPQ below N.
  -Q, --min-baseq N      Skip bases with quality below N.
  -f, --include-flags N  Required flag bits.
  -F, --exclude-flags N  Excluded flag bits (default 0xF04).
      --min-depth N      Minimum depth to count a position as covered [1].
  -d, --depth N          Accepted; depth cap (no-op — the Go port does not
                         clamp depth).
  -o, --output PATH      Output path (default stdout).
  -H                     Suppress the column-header line.
  -m, --histogram        Show the per-reference histogram instead of the table.
  -D, --plot-depth       Plot summed depth per bin (implies histogram).
  -A, --ascii            Render the histogram with ASCII (".", ":") glyphs
                         instead of UTF-8 block characters (implies histogram).
  -w, --n-bins N         Histogram bin count (implies histogram) [40].
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
	cliflag.Var(fs, &regions, "r", "region", "")
	cliflag.IntVar(fs, &minMAPQ, "q", "min-mapq", 0, "")
	cliflag.IntVar(fs, &minBaseQ, "Q", "min-baseq", 0, "")
	cliflag.IntVar(fs, &incFlags, "f", "include-flags", 0, "")
	cliflag.IntVar(fs, &excFlags, "F", "exclude-flags", 0, "")
	cliflag.StringVar(fs, &outPath, "o", "output", "", "")
	fs.BoolVar(&noHdr, "H", false, "")
	fs.BoolVar(&showHelp, "h", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "v", false, "")
	fs.BoolVar(&showVer, "version", false, "")
	// Histogram mode (upstream coverage.c getopt "Ao:l:q:Q:hHw:r:b:md:D").
	// Any of -A / -m / -D / -w turns on histogram output. -A selects the
	// ASCII glyph ramp; -D plots summed depth; -w sets the bin count.
	var (
		covAscii    bool
		covHistM    bool
		covHistD    bool
		covBins     int
		covMinLen   int
		covMinDepth int
	)
	cliflag.BoolVar(fs, &covAscii, "A", "ascii", false, "")
	cliflag.BoolVar(fs, &covHistM, "m", "histogram", false, "")
	cliflag.BoolVar(fs, &covHistD, "D", "plot-depth", false, "")
	cliflag.IntVar(fs, &covBins, "w", "n-bins", 0, "")
	cliflag.IntVar(fs, &covMinLen, "l", "min-read-len", 0, "")
	fs.IntVar(&covMinDepth, "min-depth", 0, "")
	// Upstream extras the Go port accepts but does not act on:
	//   -b FILE   file of input BAM paths (the Go port takes one input)
	//   -d N      maximum coverage depth (the Go delta-map path does not clamp)
	var (
		covFileList string
		covMaxDepth int
	)
	fs.StringVar(&covFileList, "b", "", "")
	cliflag.IntVar(fs, &covMaxDepth, "d", "depth", 0, "")

	if err := cliflag.Parse(fs, args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, coverageUsage)
		return 2
	}
	_ = covFileList
	_ = covMaxDepth
	// -A / -m / -D / -w each enable histogram mode upstream.
	if covAscii || covHistM || covHistD || covBins > 0 {
		hist = true
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
		MinReadLen:   covMinLen,
		MinDepth:     covMinDepth,
		Histogram:    hist,
		PlotDepth:    covHistD,
		ASCII:        covAscii,
		NBins:        covBins,
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
  -v             Verbose progress (accepted; no-op). Matches upstream split,
                 where -v is verbose rather than version.
  -M, -p, -@ N   Accepted upstream compatibility stubs (no-op).
  -h, --help     Show this help.
      --version  Show version.
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
	fs.BoolVar(&showVer, "version", false, "")
	// Upstream split (bam_split.c getopt "vf:h:u:d:M:p:@:") exposes a few
	// flags this port does not act on. Register them as accepted stubs so
	// legacy command lines and bundled clusters parse:
	//   -v        verbose progress (accepted no-op)
	//   -M N      maximum number of split files (accepted no-op)
	//   -p N      per-file read budget (accepted no-op)
	//   -@ N      threads (accepted no-op)
	// (-v is upstream's verbose switch here, not version; this port keeps
	// --version for the version banner.)
	var (
		splitVerbose  bool
		splitMaxSplit int
		splitPerFile  int
		splitThreads  int
	)
	fs.BoolVar(&splitVerbose, "v", false, "")
	fs.IntVar(&splitMaxSplit, "M", 0, "")
	fs.IntVar(&splitPerFile, "p", 0, "")
	cliflag.IntVar(fs, &splitThreads, "@", "threads", 0, "Threads (accepted, ignored)")

	if err := cliflag.Parse(fs, args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, splitUsage)
		return 2
	}
	_ = splitVerbose
	_ = splitMaxSplit
	_ = splitPerFile
	_ = splitThreads
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
  -d, --max-dist N        Optical-duplicate distance: flag a duplicate as
                          optical (dt:Z:SQ) when its Illumina read-name tile
                          coordinates are within N of the kept original.
  -s, --stats             Print duplicate statistics (to stderr / -f file).
  -f, --stats-file FILE   Write the statistics to FILE (implies -s).
  -S                      Also mark supplementary/secondary duplicates.
  -m, --mode {t|s|tp}     Key mode: template (default), sequence, or
                          template+position (folded into template).
  -T, --tmpdir PATH       Accepted; v1 keeps state in memory.
  -l, --max-len N         Max read length considered (default 300).
      --include-flags N   Require ALL bits set.
      --exclude-flags N   Drop records with ANY bit set.
  -c, --clear-tags        Clear pre-existing dup tags (do/dt/mc).
  -t, --add-tag           Write the 'do' aux tag on flagged duplicates.
  -@, --threads N         BGZF compression worker count for BAM output.
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
		doStats    bool
		statsFile  string
		supp       bool
	)
	cliflag.BoolVar(fs, &removeDups, "r", "remove-dups", false, "")
	cliflag.IntVar(fs, &maxDist, "d", "max-dist", 0, "")
	// -s/--stats prints the duplicate-statistics summary (upstream's -s).
	cliflag.BoolVar(fs, &doStats, "s", "stats", false, "")
	cliflag.StringVar(fs, &statsFile, "f", "stats-file", "", "")
	fs.BoolVar(&supp, "S", false, "")
	// -m/--mode is the key-mode flag (upstream's -m). The mode list is
	// {t|s|tp}; default template.
	cliflag.StringVar(fs, &modeStr, "m", "mode", "t", "")
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
	// Upstream markdup (bam_markdup.c getopt "rsl:StT:O:@:f:d:cm:u") extras the
	// Go port accepts but does not act on:
	//   -O FMT    output format (this port emits BAM)
	//   -u        uncompressed BAM output (accepted no-op)
	var (
		mdOutFmt string
		mdUncomp bool
	)
	cliflag.StringVar(fs, &mdOutFmt, "O", "output-fmt", "", "Output format (accepted)")
	fs.BoolVar(&mdUncomp, "u", false, "")
	if err := cliflag.Parse(fs, args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, markdupUsage)
		return 2
	}
	_ = mdOutFmt
	_ = mdUncomp
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
	// -f implies -s. The CLI owns the destination so it can prepend the
	// upstream-style "COMMAND:" line (which echoes the argv) to both the
	// stderr report and the -f stats file; the library only formats the
	// counter block.
	wantStats := doStats || statsFile != ""
	var statsBuf bytes.Buffer
	res, err := samtools.Markdup(opener, out, samtools.MarkdupOptions{
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
		Supp:         supp,
		Stats:        wantStats,
		StatsWriter:  &statsBuf,
		NoPG:         noPG,
	})
	if err != nil {
		fmt.Fprintf(os.Stderr, "samtools markdup: %v\n", err)
		return 1
	}
	if wantStats {
		// Upstream prints "COMMAND: <argv>" then the counter block then a
		// trailing blank line.
		command := fmt.Sprintf("COMMAND: samtools markdup %s\n", strings.Join(args, " "))
		report := append([]byte(command), statsBuf.Bytes()...)
		report = append(report, '\n')
		if statsFile != "" {
			if werr := os.WriteFile(statsFile, report, 0o644); werr != nil {
				fmt.Fprintf(os.Stderr, "samtools markdup: cannot write stats to %s: %v\n", statsFile, werr)
				return 1
			}
		} else {
			_, _ = os.Stderr.Write(report)
		}
	}
	_ = res
	return 0
}

// ----- stats ------------------------------------------------------------

const statsUsage = `samtools stats - exhaustive per-file alignment statistics.

Usage:
  samtools stats [options] <in.bam> [region ...]

Emits an upstream-compatible text report covering every section: SN
(Summary Numbers), RL/FRL/LRL (read-length), MAPQ, IS (insert sizes),
FFQ/LFQ (per-cycle qualities), GCF/GCL (GC-fraction histograms), GCC/GCT
(per-cycle GC), FBC/FTC/LBC/LTC (per-fragment ACGT), the per-barcode
ACGT-content and quality tables (<tag>C/<tag>Q for the BC/CR/OX/RX tag
pairs, emitted automatically when those aux tags are present), IC/ID
(indels), COV (coverage distribution), GCD (GC-depth distribution), MPC
(mismatches per cycle, with --ref-seq) and RFS (reference statistics,
with --ref-stats).

Options:
  -r, --ref-seq FASTA      Reference FASTA. When given, GCD derives GC content
                           from the reference instead of the read sequences,
                           and the MPC mismatches-per-cycle section is emitted.
      --ref-stats          Emit the RFS reference-statistics section.
      --ref-stats-chunk N  Reference-fetch chunk width in megabytes for RFS
                           (default 1); affects only how the FASTA is read.
  -c, --coverage MIN[,MAX[,STEP]]
                           Coverage histogram bins (default 1,1000,1).
      --GC-depth N         GC-depth bin width in reference bases
                           (default 20000).
  -l, --required-flag N    Require ALL bits set.
  -F, --filtering-flag N   Drop records with ANY bit set.
  -d, --max-depth N        Cap depth used in coverage.
      --min-mapq N         Skip records with MAPQ < N.
  -q, --trim-quality N     BWA-style 3'-end quality-trim threshold; feeds the
                           "bases trimmed" SN counter (default 0 = disabled).
      --remove-dups        Skip duplicate-flagged records.
      --remove-overlaps    Accept (no-op in v1).
  -i, --insert-size N      Max insert size for the IS section (default 8000).
  -x, --sparse             Suppress IS rows that have no insertions.
  -t, --target-regions F   Restrict stats to a target-regions file. Each line
                           is "seq-name beg end", 1-based inclusive (NOT BED).
  -g, --cov-threshold N    Coverage threshold for the target-genome coverage
                           SN line (default 0; requires -t or a region).
  -@, --threads N          Accepted; v1 is single-threaded.
  -o, --output FILE        Output path (default stdout).
  -h, --help               Show this help.
  -v, --version            Show version.

Positional region arguments (chr, chr:beg-end) after the input restrict the
statistics to reads overlapping those regions, like -t. They require an
indexed input (.bai/.csi). -t and positional regions are mutually exclusive.
`

// openStatsInput opens the stats input file. With threads >= 2 it opens the raw
// (still-compressed) bytes so samtools.Stats can decode BGZF in parallel;
// otherwise it uses the standard decompressing opener. The decoded records are
// identical for both paths.
func openStatsInput(path string, threads int) (io.ReadCloser, error) {
	if threads >= 2 {
		return iohelper.OpenRaw(path)
	}
	return iohelper.OpenReader(path)
}

// statsInputHasIndex reports whether a coordinate index (.bai or .csi) sits
// alongside the stats input file. It mirrors htslib's sam_index_load search
// order: "<file>.bai"/"<file>.csi" and the extension-replaced "<base>.bai"/
// "<base>.csi". Stats only needs to know an index exists to satisfy upstream's
// "Random alignment retrieval only works for indexed files" precondition for
// command-line region arguments; the actual restriction is a streaming
// isInRegions overlap filter, so we never read the index bytes.
func statsInputHasIndex(path string) bool {
	candidates := []string{path + ".bai", path + ".csi"}
	if ext := filepath.Ext(path); ext != "" {
		base := path[:len(path)-len(ext)]
		candidates = append(candidates, base+".bai", base+".csi")
	}
	for _, c := range candidates {
		if fi, err := os.Stat(c); err == nil && !fi.IsDir() {
			return true
		}
	}
	return false
}

func runStats(args []string) int {
	fs := flag.NewFlagSet("samtools stats", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	var (
		refSeq      string
		coverage    string
		gcdBinSize  int
		reqFlag     int
		filtFlag    int
		maxDepth    int
		minMapQ     int
		trimQuality int
		covThresh   int
		removeDups  bool
		removeOvl   bool
		insertSize  int
		sparse      bool
		targetBED   string
		refStats    bool
		refStatsChk int
		threads     int
		outPath     string
		showHelp    bool
		showVer     bool
	)
	cliflag.StringVar(fs, &refSeq, "r", "ref-seq", "", "")
	cliflag.StringVar(fs, &coverage, "c", "coverage", "", "")
	fs.IntVar(&gcdBinSize, "GC-depth", 0, "")
	cliflag.IntVar(fs, &reqFlag, "l", "required-flag", 0, "")
	cliflag.IntVar(fs, &filtFlag, "F", "filtering-flag", 0, "")
	cliflag.IntVar(fs, &maxDepth, "d", "max-depth", 0, "")
	cliflag.IntVar(fs, &minMapQ, "", "min-mapq", 0, "")
	cliflag.IntVar(fs, &trimQuality, "q", "trim-quality", 0, "")
	cliflag.IntVar(fs, &covThresh, "g", "cov-threshold", 0, "")
	fs.BoolVar(&removeDups, "remove-dups", false, "")
	fs.BoolVar(&removeOvl, "remove-overlaps", false, "")
	cliflag.IntVar(fs, &insertSize, "i", "insert-size", 8000, "")
	cliflag.BoolVar(fs, &sparse, "x", "sparse", false, "")
	cliflag.StringVar(fs, &targetBED, "t", "target-regions", "", "")
	fs.BoolVar(&refStats, "ref-stats", false, "")
	fs.IntVar(&refStatsChk, "ref-stats-chunk", 0, "")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "")
	cliflag.StringVar(fs, &outPath, "o", "output", "", "")
	fs.BoolVar(&showHelp, "h", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "v", false, "")
	fs.BoolVar(&showVer, "version", false, "")
	// Upstream stats (stats.c getopt
	// "?hdsXxpr:c:l:i:t:m:q:f:F:g:I:S:P:@:") exposes a number of additional
	// short flags. Register them as accepted stubs so legacy command lines and
	// bundled clusters parse:
	//   -p        remove overlaps (upstream short for our --remove-overlaps)
	//   -s        deprecated/no-op upstream
	//   -X        expect an explicit index-file argument (accepted no-op)
	//   -m FLOAT  insert-size main-bulk fraction (accepted no-op)
	//   -I STR    restrict to a single read-group ID (accepted no-op)
	//   -S STR    split output by an aux tag (accepted no-op)
	//   -P STR    prefix for -S split output (accepted no-op)
	// (-d, -l and -f are intentionally this port's max-depth / required-flag /
	// included flags rather than upstream's remove-dups / read-len / required
	// spellings; see the per-flag usage above.)
	var (
		stRemoveOvlP bool
		stDeprecS    bool
		stCustomIdx  bool
		stIsizeBulk  float64
		stGroupID    string
		stSplitTag   string
		stSplitPfx   string
	)
	fs.BoolVar(&stRemoveOvlP, "p", false, "")
	fs.BoolVar(&stDeprecS, "s", false, "")
	fs.BoolVar(&stCustomIdx, "X", false, "")
	fs.Float64Var(&stIsizeBulk, "m", 0, "")
	fs.StringVar(&stGroupID, "I", "", "")
	fs.StringVar(&stSplitTag, "S", "", "")
	fs.StringVar(&stSplitPfx, "P", "", "")
	if err := cliflag.Parse(fs, args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, statsUsage)
		return 2
	}
	// -p is upstream's short form of --remove-overlaps; fold it in.
	if stRemoveOvlP {
		removeOvl = true
	}
	_ = stDeprecS
	_ = stIsizeBulk
	_ = stGroupID
	_ = stSplitTag
	_ = stSplitPfx
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
	// Positional region arguments after the input file (e.g.
	// `samtools stats in.bam chr1:100-200 chr2`) restrict the statistics to
	// reads overlapping those regions. With -X (customized-index-file) the
	// first positional after the input is the index path (consumed here as a
	// no-op locator) and the regions begin one slot later, mirroring upstream
	// stats.c's `if (has_index_file) bam_idx_fname = argv[optind++]`.
	statsArgs := fs.Args()[1:]
	if stCustomIdx && len(statsArgs) > 0 {
		statsArgs = statsArgs[1:]
	}
	statsRegions := statsArgs
	if covThresh > 0 && targetBED == "" && len(statsRegions) == 0 {
		fmt.Fprintln(os.Stderr, "samtools stats: coverage percentage calculation requires a list of target regions")
		return 2
	}
	if len(statsRegions) > 0 {
		// Upstream requires an index for command-line regions
		// (sam_itr_regarray): "Random alignment retrieval only works for
		// indexed files." Mirror that precondition.
		if !stCustomIdx && !statsInputHasIndex(fs.Arg(0)) {
			fmt.Fprintln(os.Stderr, "samtools stats: random alignment retrieval only works for indexed files")
			return 1
		}
	}
	// Resolve the effective inflate worker count once: a default (no -@) opts
	// into parallel BGZF decode across cores. When that resolves to >= 2 the
	// file is opened raw (undecompressed) so samtools.Stats can inflate the
	// BGZF blocks in parallel; otherwise the standard decompressing opener is
	// used. The decoded records — and thus the report — are identical either
	// way; only decode throughput changes.
	effThreads := samtools.ReadDecodeThreads(threads)
	in, err := openStatsInput(fs.Arg(0), effThreads)
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
		GcdBinSize:     gcdBinSize,
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
		Regions:        statsRegions,
		RefStats:       refStats,
		RefStatsChunk:  refStatsChk,
		TrimQuality:    trimQuality,
		CovThreshold:   covThresh,
		Threads:        effThreads,
		Argv:           os.Args,
		Version:        fmt.Sprintf("%s+htslib-%s", version, version),
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
  -A             With -r, cap base qualities by BAQ (write ZQ tag).
  -r             Compute the BQ:Z BAQ tag (or, with -A, cap baseQ by BAQ).
  -E             Extended BAQ mode (used with -r).
  -C INT         Cap MAPQ using mismatch quality, threshold INT (>10).
  -d             Drop all aux tags except RG.
  -q             Reduce base-quality resolution (qual/10*10+7 for qual>=3).
  -n INT         Mask matching bases of reads whose NM >= INT.
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
	var (
		dropTag bool
		binQual bool
		maxNM   int
	)
	fs.BoolVar(&useEqual, "e", false, "")
	fs.BoolVar(&outBAM, "b", false, "")
	fs.BoolVar(&uncomp, "u", false, "")
	fs.BoolVar(&sInFmt, "S", false, "")
	fs.BoolVar(&adjustA, "A", false, "")
	fs.BoolVar(&realnR, "r", false, "")
	fs.BoolVar(&extBAQ, "E", false, "")
	fs.BoolVar(&dropTag, "d", false, "")
	fs.BoolVar(&binQual, "q", false, "")
	fs.IntVar(&maxNM, "n", 0, "")
	fs.BoolVar(&quiet, "Q", false, "")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "")
	cliflag.StringVar(fs, &outPath, "o", "output", "", "")
	// Accept-and-ignore stubs for upstream parity (see
	// docs/PARITY_ROADMAP.md#samtools-calmd-deferred). Behaviour is
	// deferred; flag parse must not hard-error.
	var (
		clearMDNM bool
		capQ      int
		noPG      bool
		hashQNM   bool
	)
	fs.BoolVar(&clearMDNM, "N", false, "")
	fs.IntVar(&capQ, "C", 0, "")
	fs.BoolVar(&noPG, "no-PG", false, "")
	fs.BoolVar(&hashQNM, "hash-qnm", false, "")
	fs.BoolVar(&showHelp, "h", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "v", false, "")
	fs.BoolVar(&showVer, "version", false, "")

	// Route through cliflag.Parse so POSIX getopt-style short-flag bundling
	// works the way upstream calmd's getopt parser (bam_md.c getopt
	// "EqQreuNhbSC:n:Ad@:") accepts it. All upstream short flags are already
	// registered above. (Our -h is help, whereas upstream -h is the HASH_QNM
	// switch; the latter remains available as --hash-qnm.)
	if err := cliflag.Parse(fs, args); err != nil {
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
		CapMapQ:      capQ,
		Quiet:        quiet,
		DropTags:     dropTag,
		BinQual:      binQual,
		MaxNM:        maxNM,
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
	// Route through cliflag.Parse so POSIX getopt-style short-flag bundling
	// works the way upstream import's getopt parser (bam_import.c getopt
	// "1:2:s:0:bhiT:r:R:o:O:u@:N") accepts it. Upstream's --i1/--i2 index
	// inputs are long-only options; cliflag.Parse therefore reads a single
	// dash "-i1 FILE" as the getopt cluster "-i -1 FILE" (casava + read1),
	// exactly as upstream does, while "--i1 FILE" still names the index input.
	if err := cliflag.Parse(fs, args); err != nil {
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
  -b STR       Output prefix for per-haplotype BAMs. When set, three
               BAM files are written: <prefix>.0.bam, <prefix>.1.bam,
               and <prefix>.chimera.bam, alongside the TSV stream on
               stdout/-o.
  -q INT       Min het Phred-LOD: a pileup column is treated as a
               variant (het) site only when its gl2cns LOD reaches
               this. (Default 37.)
  -Q, --min-BQ INT
               Min base quality in het calling. (Default 13.)
  -D INT       Max depth observed per position. (Default 256.)
  -F, --no-fix-chimera
               Do not attempt to fix chimeras (disable the per-read
               chimera-repair pass; the read goes to its majority
               haplotype bucket regardless of split evidence).
  -A           In -b mode, drop reads with ambiguous phase: route them
               to <prefix>.chimera.bam rather than keeping them in
               their majority bucket.
  -l FILE      List of sites to phase (CHROM<TAB>POS, 1-based).
  -e           With -l, phase ONLY the listed sites (exclusive): sites
               outside the list are dropped even when their LOD passes.
  -o, --output PATH  Output TSV path (default stdout).
  -h, --help   Show this help.
      --version  Show version.

The chimera-repair flip-point search is the deterministic Go port of
upstream's fragphase scoring loop. Ambiguous reads in -b mode are
routed using a seeded math/rand source (seed 1); RNG byte-parity with
upstream's drand48() is not pursued. See docs/PARITY_ROADMAP.md.
`

func runPhase(args []string) int {
	fs := flag.NewFlagSet("samtools phase", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	var (
		blockK       int
		outPrefix    string
		minVarLOD    int
		minBaseQ     int
		maxDepth     int
		noFixChimera bool
		dropAmbig    bool
		outPath      string
		showHelp     bool
		showVer      bool
	)
	fs.IntVar(&blockK, "k", samtools.DefaultPhaseBlockWindow, "")
	fs.StringVar(&outPrefix, "b", "", "")
	// Upstream phase.c getopt string is "Q:eFq:k:b:l:D:A": -q is the
	// minimum het Phred-LOD (g.min_varLOD, default 37) and -Q is the
	// minimum base quality (g.min_baseQ, default 13). phase has NO
	// MAPQ CLI flag — MAPQ only enters via min(baseQ, mapq) during het
	// detection and the per-read core.qual==0 skip in fragment build.
	fs.IntVar(&minVarLOD, "q", samtools.DefaultPhaseMinVarLOD, "")
	fs.IntVar(&minBaseQ, "Q", samtools.DefaultPhaseMinBaseQ, "")
	fs.IntVar(&maxDepth, "D", samtools.DefaultPhaseMaxDepth, "")
	cliflag.BoolVar(fs, &noFixChimera, "F", "no-fix-chimera", false, "")
	fs.BoolVar(&dropAmbig, "A", false, "")
	// Upstream phase.c: `-l FILE` loads a site list and `-e` makes
	// the list exclusive (positions outside the list are dropped).
	// Both are commented out of upstream's usage block (phase.c
	// usage()), but the flag table still consumes them and the
	// loadpos path remains live (phase.c:526-559).
	var (
		listFile      string
		listExclusive bool
	)
	fs.BoolVar(&listExclusive, "e", false, "")
	fs.StringVar(&listFile, "l", "", "")
	// --no-PG: upstream uses this to suppress @PG injection in the
	// per-haplotype BAMs written under -b. Our port never injects @PG
	// (those BAMs use a verbatim copy of the input header), so the flag
	// is accepted-and-ignored. The TSV stream itself never carries @PG
	// either, so the flag is a no-op for the byte-parity comparison.
	var phaseNoPG bool
	fs.BoolVar(&phaseNoPG, "no-PG", false, "")
	cliflag.StringVar(fs, &outPath, "o", "output", "", "")
	fs.BoolVar(&showHelp, "h", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "version", false, "")

	// Route through cliflag.Parse so POSIX getopt-style short-flag bundling
	// works the way upstream phase's getopt parser (phase.c getopt
	// "Q:eFq:k:b:l:D:A") accepts it. All upstream short flags are registered
	// above.
	if err := cliflag.Parse(fs, args); err != nil {
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
	_ = phaseNoPG
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
		BlockWindow:    blockK,
		MinVarLOD:      minVarLOD,
		MinBaseQ:       uint8(minBaseQ),
		MaxDepth:       maxDepth,
		NoFixChimera:   noFixChimera,
		DropAmbiguous:  dropAmbig,
		OutputPrefix:   outPrefix,
		UpstreamSchema: true,
		SiteListPath:   listFile,
		ListExclusive:  listExclusive,
	}); err != nil {
		fmt.Fprintf(os.Stderr, "samtools phase: %v\n", err)
		return 1
	}
	return 0
}

// ----- targetcut -------------------------------------------------------

const targetcutUsage = `samtools targetcut - HMM consensus over a pileup.

Usage:
  samtools targetcut [options] <in.bam>

Faithful port of upstream samtools cut_target.c: builds a per-position
consensus from the SAM/BAM pileup using the MAQ revised error model,
runs a 2-state Viterbi over the per-chrom consensus track to segment
"covered, callable" regions away from "no info or uninformative"
regions, and emits one consensus SAM record per identified region.

The output line for each region has shape:

  <chrom>:<start>-<end>  0  <chrom>  <start>  60  <len>M  *  0  0  <seq>  <qual>

where <seq> is the per-position consensus base (ACGT or 'N' when no
read provided usable evidence) and <qual> is the per-position
consensus quality (Phred+33-encoded, after upstream's >>2 shift).

Options:
  -Q INT             Per-base quality cutoff (default 13).
  -i INT             HMM entry penalty, i.e. magnitude of the 0->1
                     state transition penalty (default 14000).
  -0 INT             HMM emission score in state 1 for "no info"
                     positions (default -4).
  -1 INT             HMM emission score in state 1 for "depth but no
                     callable base" positions (default 1).
  -2 INT             HMM emission score in state 1 for callable-base
                     positions (default 6).
  -f FILE            Reference FASTA. When supplied, every per-record
                     SEQ is run through BAQ (sam_prob_realn, apply+
                     extend) before its bases enter the pileup,
                     matching upstream cut_target.c::read_aln. No
                     effect with --simple.
      --simple       Emit the v1 aligned-slice FASTA per-read mode
                     instead of the HMM consensus (legacy behaviour,
                     retained for backward compatibility). With
                     --simple, -i / -0 / -1 / -2 are ignored.
  -o, --output FILE  Output file (default stdout). Accepts "-".
  -h, --help         Show this help.
      --version      Show version.
`

func runTargetcut(args []string) int {
	fs := flag.NewFlagSet("samtools targetcut", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	var (
		minBaseQ     int
		entryPenalty int
		emNoInfo     int
		emDepth      int
		emCallable   int
		fastaRef     string
		simpleMode   bool
		outPath      string
		showHelp     bool
		showVer      bool
	)
	fs.IntVar(&minBaseQ, "Q", int(samtools.DefaultTargetcutMinBaseQ), "")
	fs.IntVar(&entryPenalty, "i", samtools.DefaultTargetcutEntryPenalty, "")
	fs.IntVar(&emNoInfo, "0", samtools.DefaultTargetcutEmissionNoInfo, "")
	fs.IntVar(&emDepth, "1", samtools.DefaultTargetcutEmissionDepth, "")
	fs.IntVar(&emCallable, "2", samtools.DefaultTargetcutEmissionCallable, "")
	fs.StringVar(&fastaRef, "f", "", "")
	fs.BoolVar(&simpleMode, "simple", false, "")
	cliflag.StringVar(fs, &outPath, "o", "output", "", "")
	fs.BoolVar(&showHelp, "h", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "version", false, "")

	// Route through cliflag.Parse so POSIX getopt-style short-flag bundling
	// works the way upstream cut_target's getopt parser (cut_target.c getopt
	// "f:Q:i:o:0:1:2:") accepts it. All upstream short flags are registered
	// above.
	if err := cliflag.Parse(fs, args); err != nil {
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
		MinBaseQ:         uint8(minBaseQ),
		EntryPenalty:     entryPenalty,
		EmissionNoInfo:   emNoInfo,
		EmissionDepth:    emDepth,
		EmissionCallable: emCallable,
		SimpleMode:       simpleMode,
		FastaRef:         fastaRef,
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
  -r, --region REG          Limit to "chr[:start-end]" region; may repeat.
  -f, --format FMT          Output format: fasta (default), fastq, pileup.
  -l, --line-len INT        Wrap FASTA/FASTQ lines at INT (default 70).
  -o, --output FILE         Output file (default stdout).
  -m, --mode STR            Algorithm: simple, bayesian, bayesian_r,
                            bayesian_p, bayesian_m, bayesian_116
                            (default bayesian, i.e. bayesian_r/RECALL).
  -a                        Output all bases; repeat (-aa) to also emit
                            contigs with no reads.
      --rf, --incl-flags N  Require ALL these flag bits set (accepted; v1 has
                            a fixed include set).
      --ff, --excl-flags N  Drop reads with ANY of these flag bits set
                            (default UNMAP,SECONDARY,QCFAIL,DUP).
      --min-MQ INT          Skip reads with MAPQ below INT (default 0).
                            No short alias upstream.
      --min-BQ INT          Skip bases with quality below INT (default 0).
      --show-del yes|no     Show deletions as '*' (default no).
                            Honoured in pileup mode too.
      --show-ins yes|no     Include insertions in FASTA/FASTQ (default yes).
      --mark-ins            Prepend '+' before inserted base/qual (default off).
  -A, --ambig               Emit IUPAC ambiguity codes for hets.
  -d, --min-depth INT       Minimum depth (default 1).
      --het-only            Restrict output to heterozygous-called positions;
                            homozygous/no-call positions become 'N' (FASTA/
                            FASTQ) or are omitted (pileup). NOTE: this fixes an
                            upstream dead-option bug — samtools parses
                            --het-only but never acts on it.
      --ref-qual INT        QUAL for reference-filled bases in FASTQ output
                            (default 0; used with -T/--reference).
      --default-qual INT    Default qual when a base has none (accepted; v1
                            uses the per-base qual unchanged).
  -Z, --block-size INT      Chromosome block size (accepted; v1 is single-pass).
      --input-fmt-option OPT[=VAL]
                            Input-format option (accepted; ignored).
      --verbosity INT       Verbosity level (accepted; ignored).
  -O, --output-fmt FMT      Output format (accepted; v1 honours -f).
      --write-index         Write index alongside output (accepted; v1
                            emits text).

For simple mode (the v1 fallback):
  -q, --use-qual            Use base quality in the score (default off).
      --no-use-qual         Force frequency-only counting (the default).
  -c, --call-fract FLOAT    Minimum (best+second)/total score required to
                            make a call (default 0.75).
  -H, --het-fract FLOAT     Minimum second/best score for a het call (default 0.5).

For the bayesian mode (the default):
  -C, --cutoff INT          Bayesian cutoff quality (default 10).
      --adj-qual            Modify quality with local minima (default on).
      --no-adj-qual         Disable adj-qual.
      --use-MQ              Use MAPQ in calculation (default on).
      --no-use-MQ           Disable use-MQ.
      --adj-MQ              Modify MAPQ by local NM (default on).
      --no-adj-MQ           Disable adj-MQ.
      --NM-halo INT         Window for NM count in adj-MQ (default 50).
      --SC-cost INT         Soft-clip cost per base (default 60).
      --scale-MQ FLOAT      Scale MAPQ (default 1.00).
      --low-MQ INT          Floor MAPQ (default 1).
      --high-MQ INT         Cap MAPQ (default 60).
      --P-het FLOAT         Het-site probability.
      --P-indel FLOAT       Indel-site probability.
      --het-scale FLOAT     Het SNP probability multiplier.
  -p, --homopoly-fix        Spread low-qual bases at homopolymer ends.
      --homopoly-score FLOAT  Quality fraction adjustment for -p.
      --homopoly-redux FLOAT  Quality reduction for -p (default 0.01).
  -t, --qual-calibration STR  Quality calibration file or :preset
                            (accepted; the FLAT identity table is used).
  -X, --config STR          Predefined config (accepted; not applied).

Global options:
  -T, --reference FILE      Reference FASTA. No-coverage / gap positions are
                            filled from the reference instead of 'N' (pileup,
                            FASTA and FASTQ; FASTQ gap quality is --ref-qual).
  -@, --threads INT         Threads (accepted; v1 is single-threaded).
      --ignore-overlaps     Accepted; v1 does not deduplicate mate overlaps.
  -h, --help                Show this help.
  -v, --version             Show version.

Notes:
  - The default mode is bayesian (RECALL); the Gap5 posterior caller,
    the NM-halo MAPQ adjustment, and all four bayesian sub-modes are
    implemented.
  - Frequency-only counting is the default for simple mode (upstream
    use_qual=0). Pass -q/--use-qual to weight by per-base quality.
  - -t/--qual-calibration and -X/--config are accepted but apply the
    FLAT identity calibration only. -T/--reference no-coverage fill is
    implemented (pileup/FASTA/FASTQ, byte-validated vs upstream).
    Tracked in docs/PARITY_ROADMAP.md#samtools.
`

// consensusBayesianFlags collects the bayesian-only flag values so the
// CLI can accept them all without bloating the local var declarations.
type consensusBayesianFlags struct {
	cutoff        int
	adjQual       bool
	noAdjQual     bool
	useMQ         bool
	noUseMQ       bool
	adjMQ         bool
	noAdjMQ       bool
	nmHalo        int
	scCost        int
	scaleMQ       float64
	lowMQ         int
	highMQ        int
	pHet          float64
	pIndel        float64
	hetScale      float64
	homopolyFix   bool
	homopolyScore float64
	homopolyRedux float64
	qualCal       string
	config        string
}

// countFlag is a boolean-style flag that counts repeats, used for the
// consensus -a/-aa option (upstream's repeatable `-a`).
type countFlag int

func (c *countFlag) String() string   { return strconv.Itoa(int(*c)) }
func (c *countFlag) Set(string) error { *c++; return nil }
func (c *countFlag) IsBoolFlag() bool { return true }

// permuteFlagArgs reorders args so every recognised option flag (and its
// value, if any) precedes the positional arguments, letting Go's flag
// package — which otherwise stops at the first non-flag token — accept
// flags given after the input file. It returns the reordered slice and
// true; if anything looks unrecognised it returns (nil, false) so the
// caller falls back to plain parsing and the normal error path.
func permuteFlagArgs(fs *flag.FlagSet, args []string) ([]string, bool) {
	// Collect the flag names and which take a value.
	isBool := map[string]bool{}
	known := map[string]bool{}
	fs.VisitAll(func(f *flag.Flag) {
		known[f.Name] = true
		if bf, ok := f.Value.(interface{ IsBoolFlag() bool }); ok && bf.IsBoolFlag() {
			isBool[f.Name] = true
		}
	})
	known["h"], known["help"], known["v"], known["version"] = true, true, true, true
	isBool["h"], isBool["help"], isBool["v"], isBool["version"] = true, true, true, true

	var flags, positional []string
	for i := 0; i < len(args); i++ {
		a := args[i]
		if a == "--" {
			positional = append(positional, args[i+1:]...)
			break
		}
		if len(a) < 2 || a[0] != '-' {
			positional = append(positional, a)
			continue
		}
		name := strings.TrimLeft(a, "-")
		eq := strings.IndexByte(name, '=')
		if eq >= 0 {
			name = name[:eq]
		}
		if !known[name] {
			if a[1] != '-' && len(name) > 1 {
				short := name[:1]
				// Glued short-flag form "-C0": a known non-bool flag
				// immediately followed by its value.
				if known[short] && !isBool[short] {
					flags = append(flags, "-"+short, name[1:])
					continue
				}
				// Repeated bool short flag "-aa": expand to "-a -a".
				if known[short] && isBool[short] {
					allSame := true
					for k := 0; k < len(name); k++ {
						if name[k] != short[0] {
							allSame = false
							break
						}
					}
					if allSame {
						for k := 0; k < len(name); k++ {
							flags = append(flags, "-"+short)
						}
						continue
					}
				}
			}
			// Unknown flag: bail and let plain parsing report it.
			return nil, false
		}
		flags = append(flags, a)
		if eq < 0 && !isBool[name] && i+1 < len(args) {
			i++
			flags = append(flags, args[i])
		}
	}
	return append(flags, positional...), true
}

func runConsensus(args []string) int {
	fs := flag.NewFlagSet("samtools consensus", flag.ContinueOnError)
	fs.SetOutput(io.Discard)

	var (
		// Output / format / region.
		formatStr string
		modeStr   string
		allPos    countFlag
		regions   multiString
		outPath   string
		lineLen   int

		// Simple-mode scoring.
		callFract float64
		hetFract  float64
		minDepth  int
		ambig     bool
		useQual   bool
		noUseQual bool

		// Filtering.
		inclFlags string
		exclFlags string
		minMQ     int
		minBQ     int

		// Insertion / deletion display.
		showDel string
		showIns string
		markIns bool
		hetOnly bool

		// Reference / global.
		refFasta   string
		refQual    int
		defaultQ   int
		blockSize  int
		threads    int
		ignoreOvl  bool
		inFmtOpt   multiString
		verbosity  int
		outputFmt  string
		writeIndex bool

		// Bayesian-only (accepted, fall back to simple).
		bay consensusBayesianFlags

		showHelp bool
		showVer  bool
	)
	// Default values are wired up here to match upstream's
	// consensus_opts initialisers (bam_consensus.c:2981+).
	cliflag.StringVar(fs, &formatStr, "f", "format", "fasta", "Output format")
	cliflag.StringVar(fs, &modeStr, "m", "mode", "bayesian", "Consensus mode")
	fs.Var(&allPos, "a", "Output all bases (repeat for all contigs)")
	cliflag.Var(fs, &regions, "r", "region", "")
	cliflag.StringVar(fs, &outPath, "o", "output", "", "Output path")
	cliflag.IntVar(fs, &lineLen, "l", "line-len", 70, "Line wrap")

	cliflag.Float64Var(fs, &callFract, "c", "call-fract", 0.75, "Min call fraction")
	cliflag.Float64Var(fs, &hetFract, "H", "het-fract", 0.5, "Min het fraction")
	cliflag.IntVar(fs, &minDepth, "d", "min-depth", 1, "Min depth")
	cliflag.BoolVar(fs, &ambig, "A", "ambig", false, "Emit IUPAC ambig codes")
	cliflag.BoolVar(fs, &useQual, "q", "use-qual", false, "Use base quality")
	fs.BoolVar(&noUseQual, "no-use-qual", false, "")

	// Upstream spells these --rf/--incl-flags and --ff/--excl-flags.
	cliflag.StringVar(fs, &inclFlags, "", "incl-flags", "", "Required flag bits")
	cliflag.StringVar(fs, &inclFlags, "", "rf", "", "Required flag bits (alias)")
	cliflag.StringVar(fs, &exclFlags, "", "excl-flags", "", "Excluded flag bits")
	cliflag.StringVar(fs, &exclFlags, "", "ff", "", "Excluded flag bits (alias)")
	cliflag.IntVar(fs, &minMQ, "", "min-MQ", 0, "Min MAPQ")
	cliflag.IntVar(fs, &minBQ, "", "min-BQ", 0, "Min base quality")

	cliflag.StringVar(fs, &showDel, "", "show-del", "no", "Show deletions")
	cliflag.StringVar(fs, &showIns, "", "show-ins", "yes", "Include insertions")
	cliflag.BoolVar(fs, &markIns, "", "mark-ins", false, "Mark inserted bases with '+'")
	cliflag.BoolVar(fs, &hetOnly, "", "het-only", false, "Restrict output to heterozygous positions (fixes upstream dead-option bug)")

	cliflag.StringVar(fs, &refFasta, "T", "reference", "", "Reference FASTA")
	cliflag.IntVar(fs, &refQual, "", "ref-qual", 0, "Reference qual")
	cliflag.IntVar(fs, &defaultQ, "", "default-qual", 10, "Default qual")
	cliflag.IntVar(fs, &blockSize, "Z", "block-size", 500000, "Block size")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "Threads")
	cliflag.BoolVar(fs, &ignoreOvl, "", "ignore-overlaps", false, "Ignore overlapping mates (accepted; not implemented)")
	fs.Var(&inFmtOpt, "input-fmt-option", "")
	cliflag.IntVar(fs, &verbosity, "", "verbosity", 0, "Verbosity")
	cliflag.StringVar(fs, &outputFmt, "O", "output-fmt", "", "Output format (accepted; v1 honours -f instead)")
	cliflag.BoolVar(fs, &writeIndex, "", "write-index", false, "Write index alongside output (accepted; v1 emits text)")

	// Bayesian-only flag landing pads.
	cliflag.IntVar(fs, &bay.cutoff, "C", "cutoff", 10, "Bayesian cutoff (accepted)")
	cliflag.BoolVar(fs, &bay.adjQual, "", "adj-qual", false, "Bayesian adj-qual on (accepted)")
	cliflag.BoolVar(fs, &bay.noAdjQual, "", "no-adj-qual", false, "Bayesian adj-qual off (accepted)")
	cliflag.BoolVar(fs, &bay.useMQ, "", "use-MQ", false, "Bayesian use-MQ on (accepted)")
	cliflag.BoolVar(fs, &bay.noUseMQ, "", "no-use-MQ", false, "Bayesian use-MQ off (accepted)")
	cliflag.BoolVar(fs, &bay.adjMQ, "", "adj-MQ", false, "Bayesian adj-MQ on (accepted)")
	cliflag.BoolVar(fs, &bay.noAdjMQ, "", "no-adj-MQ", false, "Bayesian adj-MQ off (accepted)")
	cliflag.IntVar(fs, &bay.nmHalo, "", "NM-halo", 50, "Bayesian NM-halo (accepted)")
	cliflag.IntVar(fs, &bay.scCost, "", "SC-cost", 60, "Bayesian SC-cost (accepted)")
	cliflag.Float64Var(fs, &bay.scaleMQ, "", "scale-MQ", 1.0, "Bayesian scale-MQ (accepted)")
	cliflag.IntVar(fs, &bay.lowMQ, "", "low-MQ", 1, "Bayesian low-MQ (accepted)")
	cliflag.IntVar(fs, &bay.highMQ, "", "high-MQ", 60, "Bayesian high-MQ (accepted)")
	cliflag.Float64Var(fs, &bay.pHet, "", "P-het", 1.0e-3, "Bayesian P-het (accepted)")
	cliflag.Float64Var(fs, &bay.pIndel, "", "P-indel", 2.0e-4, "Bayesian P-indel")
	cliflag.Float64Var(fs, &bay.hetScale, "", "het-scale", 1.0, "Bayesian het-scale (accepted)")
	cliflag.BoolVar(fs, &bay.homopolyFix, "p", "homopoly-fix", false, "Homopolymer fix (accepted)")
	cliflag.Float64Var(fs, &bay.homopolyScore, "", "homopoly-score", 0.0, "Homopolymer score (accepted)")
	cliflag.Float64Var(fs, &bay.homopolyRedux, "", "homopoly-redux", 0.01, "Homopolymer redux (accepted)")
	cliflag.StringVar(fs, &bay.qualCal, "t", "qual-calibration", "", "Quality calibration (accepted)")
	cliflag.StringVar(fs, &bay.config, "X", "config", "", "Predefined config (accepted)")

	// Upstream consensus' getopt string ("...r:5f:...") reserves a bare -5
	// flag that its switch never handles (a no-op). Register it so a bundled
	// cluster that includes it still parses.
	var consDash5 bool
	fs.BoolVar(&consDash5, "5", false, "")

	fs.BoolVar(&showHelp, "h", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "v", false, "")
	fs.BoolVar(&showVer, "version", false, "")
	_ = consDash5

	// Two-stage argument handling for upstream parity:
	//   1. cliflag.Normalize expands POSIX getopt short-flag clusters
	//      (`-aA`, `-C0`, `-aa`) into canonical one-flag-per-token form, the
	//      way upstream consensus' getopt parser (bam_consensus.c getopt
	//      "@:qd:c:H:r:5f:C:aAl:o:m:pt:X:T:Z:") reads them.
	//   2. permuteFlagArgs then floats the option flags ahead of the
	//      positional input, because upstream's tests put the input BAM
	//      before its flags (e.g. `consensus in.bam -m bayesian -C 0`) and
	//      Go's flag package otherwise stops at the first non-flag token.
	normalized, nerr := cliflag.Normalize(fs, args)
	if nerr != nil {
		fmt.Fprintln(os.Stderr, nerr)
		fmt.Fprint(os.Stderr, consensusUsage)
		return 2
	}
	args = normalized
	if perm, ok := permuteFlagArgs(fs, args); ok {
		args = perm
	}
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

	// Record which flags were explicitly set so the bayesian on/off
	// toggles (--adj-qual / --no-adj-qual etc.) and the homopoly-redux
	// default can be distinguished from their zero values.
	explicit := map[string]bool{}
	fs.Visit(func(f *flag.Flag) { explicit[f.Name] = true })

	// The accepted-but-not-implemented knobs (incl-flags, excl-flags,
	// block-size, input-fmt-option, verbosity, and the -t/-X bayesian config
	// knobs) live in the symbol table because cliflag took their address; we
	// deliberately don't route them into ConsensusOptions yet. They're listed
	// in consensusUsage and tracked in docs/PARITY_ROADMAP.md. -T/--reference
	// and --ref-qual ARE routed (see opts.Reference / opts.RefQual below).
	_ = inclFlags
	_ = exclFlags
	_ = blockSize
	_ = inFmtOpt
	_ = verbosity
	_ = outputFmt
	_ = writeIndex

	format, ferr := samtools.ParseConsensusFormat(formatStr)
	if ferr != nil {
		fmt.Fprintln(os.Stderr, ferr)
		return 2
	}
	mode, merr := samtools.ParseConsensusMode(modeStr)
	if merr != nil {
		fmt.Fprintln(os.Stderr, merr)
		return 2
	}

	// --no-use-qual wins over -q/--use-qual when both are given (last
	// flag wins is what upstream does; flag package gives us either
	// "both seen" or "only one"; if both, we default to the explicit
	// "no" since that's the safer fallback for accidental composition).
	if noUseQual {
		useQual = false
	}

	opts := samtools.ConsensusOptions{
		Input:           fs.Arg(0),
		Format:          format,
		Mode:            mode,
		AllPositions:    allPos >= 1,
		AllContigs:      allPos >= 2,
		Regions:         []string(regions),
		MinDepth:        minDepth,
		MinCallFraction: callFract,
		MinHetFraction:  hetFract,
		AmbigCodes:      ambig,
		UseQual:         useQual,
		MinMAPQ:         uint8(minMQ),
		MinBaseQ:        uint8(minBQ),
		LineLen:         lineLen,
		IgnoreOverlaps:  ignoreOvl,
		Output:          outPath,
		Threads:         threads,
		Reference:       refFasta,
		RefQual:         refQual,
	}
	opts.ShowDel = parseYesNo(showDel, false)
	opts.NoShowIns = !parseYesNo(showIns, true)
	opts.MarkIns = markIns
	opts.HetOnly = hetOnly

	// Bayesian-mode knobs. The sub-mode comes from the -m string; the
	// remaining knobs route straight through, with the on/off toggles
	// resolved against which flags were explicitly given.
	opts.SetBayesianMode(modeStr)
	opts.ConsCutoff = bay.cutoff
	opts.ConsCutoffSet = explicit["C"] || explicit["cutoff"]
	opts.PHet = bay.pHet
	opts.PIndel = bay.pIndel
	opts.HetScale = bay.hetScale
	opts.NMHalo = bay.nmHalo
	opts.SCCost = bay.scCost
	opts.ScaleMQual = bay.scaleMQ
	opts.LowMQual = bay.lowMQ
	opts.HighMQual = bay.highMQ
	opts.DefaultQual = defaultQ
	// adj-qual: default on; --no-adj-qual disables, --adj-qual forces on.
	opts.AdjQual = true
	if explicit["no-adj-qual"] {
		opts.AdjQual = false
	}
	if explicit["adj-qual"] {
		opts.AdjQual = true
	}
	opts.AdjQualSet = explicit["no-adj-qual"] || explicit["adj-qual"]
	// use-MQ: default on; --no-use-MQ disables.
	opts.UseMQual = true
	if explicit["no-use-MQ"] {
		opts.UseMQual = false
	}
	if explicit["use-MQ"] {
		opts.UseMQual = true
	}
	opts.UseMQualSet = explicit["no-use-MQ"] || explicit["use-MQ"]
	// adj-MQ: default on; --no-adj-MQ disables.
	opts.NMAdjust = true
	if explicit["no-adj-MQ"] {
		opts.NMAdjust = false
	}
	if explicit["adj-MQ"] {
		opts.NMAdjust = true
	}
	opts.NMAdjustSet = explicit["no-adj-MQ"] || explicit["adj-MQ"]
	// homopoly-fix: -p sets the P_HOMOPOLY (0.5) multiplier;
	// --homopoly-score overrides it with an explicit value.
	if bay.homopolyFix {
		opts.HomopolyFix = 0.5
	}
	if explicit["homopoly-score"] {
		opts.HomopolyFix = bay.homopolyScore
	}
	opts.HomopolyRedux = bay.homopolyRedux
	opts.HomopolyReduxSet = explicit["homopoly-redux"]

	out, oerr := openOut(outPath)
	if oerr != nil {
		fmt.Fprintf(os.Stderr, "samtools consensus: %v\n", oerr)
		return 1
	}
	defer out.Close()
	if err := samtools.ConsensusFile(opts, out, os.Stderr); err != nil {
		fmt.Fprintf(os.Stderr, "samtools consensus: %v\n", err)
		return 1
	}
	return 0
}

// parseYesNo accepts upstream samtools' "yes"/"no"/"on"/"off"/
// "true"/"false"/"1"/"0" forms (case-insensitive). Empty returns def.
// Used by --show-del / --show-ins.
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
