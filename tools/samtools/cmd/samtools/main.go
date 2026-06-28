// Command samtools is a pure-Go reimplementation of selected samtools
// subcommands. This first slice ships `view` and `flagstat`; other
// subcommands (sort, index, depth, fastq, mpileup) will follow in later PRs.
package main

import (
	"flag"
	"fmt"
	"io"
	"os"
	"runtime/debug"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
	"github.com/yassineS/bio_ai_experiment/tools/samtools/pkg/samtools"
)

const version = "0.1.0"

const rootUsage = `samtools - pure-Go reimplementation.

Usage:
  samtools <subcommand> [options]

Subcommands:
  view          Print, filter, or convert SAM/BAM records.
  sort          Sort alignments by coordinate, name, or tag.
  index         Build a BAI index for a coordinate-sorted BAM.
  flagstat      Print flag statistics for a SAM/BAM file.
  depth         Per-position depth across one or more BAMs.
  fastq         Convert SAM/BAM to FASTQ.
  bam2fq        Alias for fastq.
  mpileup       Per-position pileup across one or more BAMs.
  idxstats      Per-reference read counts from an indexed BAM.
  quickcheck    Fast format sanity check.
  dict          Emit a sequence dictionary from a FASTA file.
  faidx         Index or extract regions from a FASTA file.
  fqidx         Index or extract regions from a FASTQ file.
  cat           Concatenate BAMs without re-sorting.
  reheader      Replace the header of a BAM in place.
  addreplacerg  Add or replace @RG line and per-record RG tag.
  fixmate       Fill in mate-read fields on a name-sorted BAM.
  merge         Merge sorted BAMs.
  coverage      Per-region coverage summary table.
  split         Split a BAM by @RG ID into per-RG files.
  markdup       Mark or remove PCR duplicate records.
  stats         Print exhaustive alignment statistics.
  calmd         Compute MD + NM aux tags by walking CIGAR vs reference.
  import        Convert FASTQ to BAM (single / paired / interleaved).
  phase         Phase haplotypes from heterozygous SNPs.
  targetcut     Emit FASTA slice of each aligned read.
  consensus     Call per-position consensus base (FASTA/FASTQ/pileup).
  tview         Text/HTML alignment viewer (-d T | -d H).
  help          Show this help.
  version       Show version.
`

func main() {
	if len(os.Args) < 2 {
		fmt.Fprint(os.Stderr, rootUsage)
		os.Exit(1)
	}
	switch os.Args[1] {
	case "view":
		os.Exit(runView(os.Args[2:]))
	case "sort":
		os.Exit(runSort(os.Args[2:]))
	case "index":
		os.Exit(runIndex(os.Args[2:]))
	case "flagstat":
		os.Exit(runFlagstat(os.Args[2:]))
	case "depth":
		os.Exit(runDepth(os.Args[2:]))
	case "fastq", "bam2fq":
		os.Exit(runFastq(os.Args[2:]))
	case "mpileup":
		os.Exit(runMpileup(os.Args[2:]))
	case "idxstats":
		os.Exit(runIdxstats(os.Args[2:]))
	case "quickcheck":
		os.Exit(runQuickcheck(os.Args[2:]))
	case "dict":
		os.Exit(runDict(os.Args[2:]))
	case "faidx":
		os.Exit(runFaidx(os.Args[2:]))
	case "fqidx":
		os.Exit(runFqidx(os.Args[2:]))
	case "cat":
		os.Exit(runCat(os.Args[2:]))
	case "reheader":
		os.Exit(runReheader(os.Args[2:]))
	case "addreplacerg":
		os.Exit(runAddReplaceRG(os.Args[2:]))
	case "fixmate":
		os.Exit(runFixmate(os.Args[2:]))
	case "merge":
		os.Exit(runMerge(os.Args[2:]))
	case "coverage":
		os.Exit(runCoverage(os.Args[2:]))
	case "split":
		os.Exit(runSplit(os.Args[2:]))
	case "markdup":
		os.Exit(runMarkdup(os.Args[2:]))
	case "stats":
		os.Exit(runStats(os.Args[2:]))
	case "calmd":
		os.Exit(runCalmd(os.Args[2:]))
	case "import":
		os.Exit(runImport(os.Args[2:]))
	case "phase":
		os.Exit(runPhase(os.Args[2:]))
	case "targetcut":
		os.Exit(runTargetcut(os.Args[2:]))
	case "consensus":
		os.Exit(runConsensus(os.Args[2:]))
	case "tview":
		os.Exit(runTview(os.Args[2:]))
	case "help", "-h", "--help":
		fmt.Print(rootUsage)
		return
	case "version", "-v", "--version":
		fmt.Println(version)
		return
	default:
		fmt.Fprintf(os.Stderr, "samtools: unknown subcommand %q\n", os.Args[1])
		fmt.Fprint(os.Stderr, rootUsage)
		os.Exit(1)
	}
}

const viewUsage = `samtools view - print, filter, or convert SAM/BAM records.

Usage:
  samtools view [options] <in.bam|in.sam> [region ...]

Options:
  -b, --bam                   Output BAM (default text SAM).
  -u, --uncompressed          Output uncompressed (level-0) BAM (implies -b).
  -C, --cram                  Output CRAM (reference-free, v3.0).
  -O, --output-fmt FMT        Force output format ('sam', 'bam' or 'cram').
      --output-fmt-option OPT CRAM output tuning, KEY=VALUE form. Repeatable.
                              Supported: qbin=8|4|2|none — apply lossy
                              Illumina quality-score binning (8-/4-/2-level)
                              to CRAM output. Ignored for SAM/BAM output.
  -h, --with-header           Include the header in SAM output.
  -H, --header-only           Print the header only.
  -c, --count                 Print only the count of matching records.
  -f, --include-flags <int>   Require ALL bits set.
  -F, --exclude-flags <int>   Drop records with ANY bit set.
  -G, --exclude-flags-all <int>  Drop only when ALL bits set.
  -q, --min-mapq <int>        Minimum MAPQ.
  -r, --read-group <string>   Keep records with this RG.
  -R, --read-groups-file <f>  File of RG IDs (one per line).
  -L, --regions-file <f>      BED of regions to keep (linear scan).
  -M, --use-multi-region-iterator
                              Accepted; we always do the full intersection.
  -d, --tag STR[:VAL]         Keep records with aux tag STR (and optional VAL).
                              May be repeated; multiple values OR within
                              the same tag.
  -D, --tag-file STR:FILE     Keep records with aux tag STR matching one of
                              the values listed in FILE (one per line).
  -N, --qname-file <f>        Keep records whose QNAME is listed in FILE.
                              Leading "^" inverts: -N ^FILE drops
                              records whose QNAME is in FILE.
  -s, --subsample <float>     Keep fraction (or "<seed>.<frac>").
  -o, --output <file>         Output file (default stdout).
  -T, --reference <fasta>     Reference FASTA for decoding reference-backed
                              CRAM input (needs a sibling .fai). Ignored
                              for SAM/BAM. CRAM input is auto-detected.
  -@, --threads <int>         BGZF compression worker count for BAM output.
      --no-PG                 Suppress @PG line emission.
      --help                  Show this help.
      --version               Show version.

Region queries (chr:start-end) use a sibling <input>.bai index when one
exists; otherwise samtools view falls back to a linear scan with a
warning to stderr. The -L/--regions-file form always performs a linear
scan over the input and keeps records whose [Pos, Pos+refLen) range
intersects any BED interval on the record's reference.
`

// inputIsCRAM reports whether the file at path begins with the CRAM magic
// number ("CRAM", 0x43 0x52 0x41 0x4D). It opens a SEPARATE handle, reads the
// first four bytes, and closes it — so it never consumes bytes from the read
// stream the decode path opens through its own normal route. Any error (the
// path is a pipe, is shorter than four bytes, or cannot be opened) is treated
// as "not a detectable CRAM file": the caller then simply skips the memory
// limit rather than failing, and the real open reports any genuine error.
func inputIsCRAM(path string) bool {
	f, err := os.Open(path)
	if err != nil {
		return false
	}
	defer f.Close()
	var magic [4]byte
	if _, err := io.ReadFull(f, magic[:]); err != nil {
		return false
	}
	return magic[0] == 'C' && magic[1] == 'R' && magic[2] == 'A' && magic[3] == 'M'
}

func runView(args []string) int {
	fs := flag.NewFlagSet("samtools view", flag.ContinueOnError)
	fs.SetOutput(io.Discard) // we print usage ourselves

	var (
		outBAM      bool
		outCRAM     bool
		outFmt      string
		outFmtOpts  multiString
		withHdr     bool
		hdrOnly     bool
		countOnly   bool
		incFlags    int
		excFlags    int
		excFlagsG   int
		minMAPQ     int
		rg          string
		rgFile      string
		regFile     string
		multiRegion bool
		customIdx   bool
		tagSpecs    multiString
		tagFiles    multiString
		qnameFile   string
		subsample   string
		outFile     string
		refFile     string
		threads     int
		noPG        bool
		showHelp    bool
		showVer     bool
	)
	cliflag.BoolVar(fs, &outBAM, "b", "bam", false, "Output BAM")
	cliflag.BoolVar(fs, &outCRAM, "C", "cram", false, "Output CRAM")
	var outUncomp bool
	cliflag.BoolVar(fs, &outUncomp, "u", "uncompressed", false, "Uncompressed BAM output (implies -b)")
	// -S is upstream samtools' legacy "input is SAM" switch (sam_view.c:999
	// `case 'S': break;`). Modern samtools auto-detects the input format, so
	// the flag is a no-op kept for backward compatibility — and so the common
	// bundled invocation `samtools view -bS in.sam` parses.
	var ignoreSAMInput bool
	cliflag.BoolVar(fs, &ignoreSAMInput, "S", "", false, "Ignored: input format is auto-detected (legacy)")
	cliflag.StringVar(fs, &outFmt, "O", "output-fmt", "", "Output format (sam|bam|cram)")
	fs.Var(&outFmtOpts, "output-fmt-option", "")
	cliflag.BoolVar(fs, &withHdr, "h", "with-header", false, "Include header")
	cliflag.BoolVar(fs, &hdrOnly, "H", "header-only", false, "Header only")
	cliflag.BoolVar(fs, &countOnly, "c", "count", false, "Count records")
	cliflag.IntVar(fs, &incFlags, "f", "include-flags", 0, "Required flags")
	cliflag.IntVar(fs, &excFlags, "F", "exclude-flags", 0, "Forbidden flags")
	cliflag.IntVar(fs, &excFlagsG, "G", "exclude-flags-all", 0, "Forbidden flags (all)")
	cliflag.IntVar(fs, &minMAPQ, "q", "min-mapq", 0, "Minimum MAPQ")
	cliflag.StringVar(fs, &rg, "r", "read-group", "", "Read-group filter")
	cliflag.StringVar(fs, &rgFile, "R", "read-groups-file", "", "File of read-group IDs")
	cliflag.StringVar(fs, &regFile, "L", "regions-file", "", "BED of regions")
	cliflag.BoolVar(fs, &multiRegion, "M", "use-multi-region-iterator", false, "Deduplicate overlapping regions (each record at most once, coordinate order)")
	cliflag.BoolVar(fs, &customIdx, "X", "customized-index", false, "Expect an explicit index-file argument after <in.bam>")
	// Upstream samtools also spells the long form `--use-index`. Accept it
	// as an alias for parity.
	fs.BoolVar(&multiRegion, "use-index", false, "")
	cliflag.Var(fs, &tagSpecs, "d", "tag", "")
	cliflag.Var(fs, &tagFiles, "D", "tag-file", "")
	cliflag.StringVar(fs, &qnameFile, "N", "qname-file", "", "File of QNAMEs to keep")
	cliflag.StringVar(fs, &subsample, "s", "subsample", "", "Subsample fraction")
	cliflag.StringVar(fs, &outFile, "o", "output", "", "Output file")
	cliflag.StringVar(fs, &refFile, "T", "reference", "", "Reference FASTA for CRAM input")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "BGZF compression worker count for BAM output")
	cliflag.BoolVar(fs, &noPG, "", "no-PG", false, "Suppress @PG line")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "version", false, "")

	// Route through cliflag.Parse so POSIX getopt-style short-flag bundling
	// (-bS == -b -S) and value concatenation (-q20 == -q 20) work the way
	// upstream samtools' getopt parser accepts them.
	if err := cliflag.Parse(fs, args); err != nil {
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
		fmt.Fprintln(os.Stderr, "samtools view: missing input file")
		fmt.Fprint(os.Stderr, viewUsage)
		return 2
	}
	input := positional[0]

	// CRAM→BAM/SAM decode keeps a per-slice working set (the decompressed
	// data-series blocks plus the slice's reconstructed records) of ~50 MB on
	// whole-chromosome inputs. Under the default GOGC the collector lets the
	// heap grow to ~2× that working set before collecting and is slow to return
	// freed pages to the OS, so the ru_maxrss high-water mark reaches ~110 MB —
	// ~5.5× upstream's exceptionally lean ~20 MB C decoder. A soft memory limit
	// pulls the resident set down toward the working-set floor without trimming
	// GOGC (which would tax CPU on every cycle): below the cap the GC stays lazy,
	// so wall time is unaffected; near the cap it scavenges idle pages back to
	// the OS. 48 MiB is the measured knee on real GIAB (chr20): it brings the
	// worst-of-N peak RSS from ~70 MB (64 MiB cap, ≈3.5×) down to ~57–61 MB
	// (≈2.9×) with no wall regression (worst-of-7 wall ~10.0 s vs ~9.8 s
	// uncapped, within run-to-run jitter). Tighter caps stop lowering RSS — it
	// plateaus at ~50–55 MB, the genuine per-slice working set — and start
	// thrashing the GC: 36 MiB costs ~+4 % wall, 32 MiB ~+11 %, 28 MiB ~+26 %,
	// 24 MiB ~+120 %, all for no further RSS gain. Reaching ≤2× would require
	// streaming the slice's series blocks rather than decompressing them all up
	// front (the codec-scratch restructure tracked separately as task #52). The
	// limit is applied ONLY for a real CRAM file input (BAM/SAM view is
	// unaffected) and skipped for stdin; SetMemoryLimit returns the previous
	// (no-limit) value, restored on return so this stays scoped to the view
	// command.
	if input != "-" && inputIsCRAM(input) {
		prevLimit := debug.SetMemoryLimit(48 << 20)
		defer debug.SetMemoryLimit(prevLimit)
	}

	var indexPath string
	var regions []string
	if customIdx {
		// samtools view -X: the index file is supplied as an explicit
		// positional argument immediately after <in.bam>, with any
		// regions following it (sam_view.c:1441-1450).
		if len(positional) < 2 {
			fmt.Fprintln(os.Stderr, "samtools view: -X requires an index-file argument after the input")
			fmt.Fprint(os.Stderr, viewUsage)
			return 2
		}
		indexPath = positional[1]
		regions = positional[2:]
	} else {
		regions = positional[1:]
	}

	// Resolve the output format from -b/-C and -O/--output-fmt. -O takes
	// precedence over the boolean shorthands; an unknown -O value is an
	// error. CRAM (-C / -O cram) wins over BAM when both are requested.
	switch strings.ToLower(outFmt) {
	case "":
		// No explicit -O; the booleans below decide.
	case "sam":
		outBAM, outCRAM = false, false
	case "bam":
		outBAM, outCRAM = true, false
	case "cram":
		outCRAM = true
	default:
		fmt.Fprintf(os.Stderr, "samtools view: unknown output format %q (want sam, bam or cram)\n", outFmt)
		return 2
	}
	if outCRAM {
		outBAM = false
	}
	// -u (uncompressed) implies BAM output, matching upstream samtools view,
	// where -u sets the output format to BAM at compression level 0. An
	// explicit CRAM request still wins; -u then only affects the (BAM) path.
	if outUncomp && !outCRAM {
		outBAM = true
	}

	// --output-fmt-option carries KEY=VALUE CRAM tuning knobs. Only qbin
	// (lossy quality-score binning) is recognised today; an unknown key is
	// rejected so a typo is not silently ignored.
	cramQBin, perr := parseOutputFmtOptions(outFmtOpts)
	if perr != nil {
		fmt.Fprintln(os.Stderr, perr)
		return 2
	}

	opts := samtools.ViewOptions{
		OutputBAM:          outBAM,
		OutputCRAM:         outCRAM,
		Uncompressed:       outUncomp,
		WithHeader:         withHdr,
		HeaderOnly:         hdrOnly,
		Count:              countOnly,
		IncludeFlags:       uint16(incFlags),
		ExcludeFlags:       uint16(excFlags),
		ExcludeFlagsAll:    uint16(excFlagsG),
		UseExcludeAll:      excFlagsG != 0,
		MinMAPQ:            uint8(minMAPQ),
		ReadGroup:          rg,
		Regions:            append([]string{}, regions...),
		RegionsEnabled:     len(regions) > 0 || regFile != "",
		BedPath:            regFile,
		MultiRegion:        multiRegion,
		NoPG:               noPG,
		Reference:          refFile,
		CRAMQualityBinning: cramQBin,
		IndexPath:          indexPath,
		Threads:            threads,
	}

	// Honour the output file extension when no format was given: a .bam
	// name implies BAM, a .cram name implies CRAM. An explicit -O is
	// authoritative and suppresses the inference, matching upstream
	// samtools (-O sam -o foo.cram still writes SAM).
	if outFmt == "" && !opts.OutputBAM && !opts.OutputCRAM && outFile != "" {
		switch {
		case strings.HasSuffix(outFile, ".cram"):
			opts.OutputCRAM = true
		case strings.HasSuffix(outFile, ".bam"):
			opts.OutputBAM = true
		}
	}

	if rgFile != "" {
		set, err := samtools.LoadReadGroupsFile(rgFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "samtools view: %v\n", err)
			return 1
		}
		opts.ReadGroupSet = set
	}

	// -d / -D tag-value filters compose with AND across distinct tags but
	// upstream rejects mixing tags, so all -d/-D must reference the same
	// tag. MergeTagFilter enforces that and unions the value sets.
	for _, spec := range tagSpecs {
		f, perr := samtools.ParseTagFilterSpec(spec)
		if perr != nil {
			fmt.Fprintln(os.Stderr, perr)
			return 2
		}
		merged, merr := samtools.MergeTagFilter(opts.TagFilters, f)
		if merr != nil {
			fmt.Fprintln(os.Stderr, merr)
			return 2
		}
		opts.TagFilters = merged
	}
	for _, spec := range tagFiles {
		f, perr := samtools.ParseTagFileSpec(spec)
		if perr != nil {
			fmt.Fprintf(os.Stderr, "samtools view: %v\n", perr)
			return 1
		}
		merged, merr := samtools.MergeTagFilter(opts.TagFilters, f)
		if merr != nil {
			fmt.Fprintln(os.Stderr, merr)
			return 2
		}
		opts.TagFilters = merged
	}

	if qnameFile != "" {
		// Upstream sam_view.c:347-352: a leading `^` inverts the
		// keep-set into a drop-set. We strip the prefix here and let
		// View flip via opts.QNameInvert.
		path := qnameFile
		if strings.HasPrefix(path, "^") {
			opts.QNameInvert = true
			path = path[1:]
		}
		set, qerr := samtools.LoadLinesFile(path)
		if qerr != nil {
			fmt.Fprintf(os.Stderr, "samtools view: %v\n", qerr)
			return 1
		}
		opts.QNameSet = set
	}

	if subsample != "" {
		frac, seed, err := samtools.ParseSubsample(subsample)
		if err != nil {
			fmt.Fprintln(os.Stderr, err)
			return 2
		}
		opts.Subsample = frac
		opts.SubsampleSeed = seed
	}

	out, err := openOut(outFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "samtools view: %v\n", err)
		return 1
	}
	defer out.Close()

	// When the input is a real file path, prefer ViewFile so we can use a
	// sibling .bai for indexed region queries. The streaming View path is
	// used for stdin / `-`.
	if len(opts.Regions) > 0 && input != "-" {
		if _, err := samtools.ViewFile(input, out, opts, os.Stderr); err != nil {
			fmt.Fprintf(os.Stderr, "samtools view: %v\n", err)
			return 1
		}
		return 0
	}

	in, err := iohelper.OpenReader(input)
	if err != nil {
		fmt.Fprintf(os.Stderr, "samtools view: %v\n", err)
		return 1
	}
	defer in.Close()

	if _, err := samtools.View(in, out, opts); err != nil {
		fmt.Fprintf(os.Stderr, "samtools view: %v\n", err)
		return 1
	}
	return 0
}

const flagstatUsage = `samtools flagstat - print flag statistics for SAM/BAM.

Usage:
  samtools flagstat <in.bam|in.sam>
`

func runFlagstat(args []string) int {
	fs := flag.NewFlagSet("samtools flagstat", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	var showHelp, showVer bool
	fs.BoolVar(&showHelp, "h", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "version", false, "")
	// Upstream flagstat (bam_stat.c getopt "@:O:") accepts -@ (threads) and
	// -O (output format). -@ now drives block-parallel BGZF *input* decode for
	// BAM (the report is identical for any thread count); -O is accepted but the
	// default text report is the only format this port emits.
	var (
		fsThreads int
		fsOutFmt  string
	)
	cliflag.IntVar(fs, &fsThreads, "@", "threads", 0, "BGZF input decode worker count")
	cliflag.StringVar(fs, &fsOutFmt, "O", "output-fmt", "", "Output format (accepted; default text)")
	if err := cliflag.Parse(fs, args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, flagstatUsage)
		return 2
	}
	_ = fsOutFmt
	if showHelp {
		fmt.Print(flagstatUsage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}
	if fs.NArg() == 0 {
		fmt.Fprint(os.Stderr, flagstatUsage)
		return 2
	}
	if err := samtools.FlagstatFile(fs.Arg(0), os.Stdout, fsThreads); err != nil {
		fmt.Fprintf(os.Stderr, "samtools flagstat: %v\n", err)
		return 1
	}
	return 0
}

const sortUsage = `samtools sort - sort SAM/BAM records by coordinate, name, or tag.

Usage:
  samtools sort [options] <in.bam|in.sam>

Options:
  -o, --output PATH           Output file (default stdout, BAM unless extension says SAM).
  -O, --output-fmt FMT        Force output format ('bam' or 'sam').
  -n, --by-name               Sort by read name (natural numeric, upstream default).
  -N, --by-natural-name       Sort by read name (lexicographic; legacy alias kept).
  -t, --by-tag TAG            Sort by an Aux tag value (e.g. "NM").
  -m, --max-mem N[K|M|G]      Per-shard memory budget (default 768M).
  -T, --tmpdir PREFIX         Temporary-file prefix.
  -l, --compress-level N      Output BGZF deflate level 0..9.
  -@, --threads N             BGZF compression worker count for BAM/shard output.
      --no-PG                 Suppress @PG injection (v1 never injects).
  -h, --help                  Show this help.
  -v, --version               Show version.
`

func runSort(args []string) int {
	fs := flag.NewFlagSet("samtools sort", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	var (
		outFile   string
		outFmt    string
		byName    bool
		byNatural bool
		byTag     string
		maxMem    string
		tmpdir    string
		compLevel int
		threads   int
		noPG      bool
		showHelp  bool
		showVer   bool
	)
	cliflag.StringVar(fs, &outFile, "o", "output", "", "Output path")
	cliflag.StringVar(fs, &outFmt, "O", "output-fmt", "", "Output format (bam|sam)")
	cliflag.BoolVar(fs, &byName, "n", "by-name", false, "Sort by QName (natural numeric, upstream default)")
	cliflag.BoolVar(fs, &byNatural, "N", "by-natural-name", false, "Sort by QName (lexicographic)")
	cliflag.StringVar(fs, &byTag, "t", "by-tag", "", "Sort by aux tag")
	cliflag.StringVar(fs, &maxMem, "m", "max-mem", "", "Per-shard memory budget")
	cliflag.StringVar(fs, &tmpdir, "T", "tmpdir", "", "Temp file prefix")
	cliflag.IntVar(fs, &compLevel, "l", "compress-level", -1, "BGZF deflate level")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "BGZF compression worker count for BAM/shard output")
	cliflag.BoolVar(fs, &noPG, "", "no-PG", false, "No @PG injection")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "version", false, "")
	// Upstream sort (bam_sort.c getopt "l:m:nNo:O:T:@:t:MI:K:uRw:H") also
	// accepts several minimiser/minhash sort-order knobs and an uncompressed
	// switch. The MinHash ordering itself is out of this port's scope; the
	// flags are registered as accepted stubs so legacy command lines — and
	// bundled clusters that include them — still parse:
	//   -M        sort by minimiser/minhash order (accepted; coordinate sort)
	//   -I FILE   minimiser reference (accepted no-op)
	//   -K N      minimiser k-mer size (accepted no-op)
	//   -w N      minimiser window (accepted no-op)
	//   -R        disable reverse-complement minimiser probing (accepted no-op)
	//   -H        do not squash homopolymers in minimiser order (accepted no-op)
	var (
		sortMinHash  bool
		sortMinRef   string
		sortKmer     int
		sortWindow   int
		sortTryRev   bool
		sortNoSquash bool
		sortUncomp   bool
	)
	cliflag.BoolVar(fs, &sortMinHash, "M", "", false, "Minimiser sort order (accepted; coordinate sort)")
	cliflag.StringVar(fs, &sortMinRef, "I", "", "", "Minimiser reference (accepted no-op)")
	cliflag.IntVar(fs, &sortKmer, "K", "", 0, "Minimiser k-mer (accepted no-op)")
	cliflag.IntVar(fs, &sortWindow, "w", "", 0, "Minimiser window (accepted no-op)")
	cliflag.BoolVar(fs, &sortTryRev, "R", "", false, "Disable reverse minimiser probe (accepted no-op)")
	cliflag.BoolVar(fs, &sortNoSquash, "H", "", false, "No homopolymer squash (accepted no-op)")
	// -u: upstream sets compression level 0 (uncompressed BAM). Force BAM
	// output at the lowest level so this matches the upstream effect.
	cliflag.BoolVar(fs, &sortUncomp, "u", "", false, "Uncompressed BAM output")

	if err := cliflag.Parse(fs, args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, sortUsage)
		return 2
	}
	if sortUncomp {
		compLevel = 0
		if outFmt == "" {
			outFmt = "bam"
		}
	}
	_ = sortMinHash
	_ = sortMinRef
	_ = sortKmer
	_ = sortWindow
	_ = sortTryRev
	_ = sortNoSquash
	if showHelp {
		fmt.Print(sortUsage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}
	if fs.NArg() == 0 {
		fmt.Fprintln(os.Stderr, "samtools sort: missing input file")
		fmt.Fprint(os.Stderr, sortUsage)
		return 2
	}

	order := samtools.SortCoordinate
	// CLI flag mapping matches upstream samtools sort:
	//   -n  -> natural numeric name sort (the upstream default for name-sort)
	//   -N  -> plain lexicographic name sort (sets natural_sort=0 upstream)
	// Our library exposes both modes; the CLI just wires them up to match.
	switch {
	case byTag != "":
		order = samtools.SortByTag
	case byNatural:
		order = samtools.SortByName
	case byName:
		order = samtools.SortByNameNatural
	}
	mem, err := samtools.ParseMemBudget(maxMem)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}

	outputSAM := false
	outputBAM := false
	switch strings.ToLower(outFmt) {
	case "sam":
		outputSAM = true
	case "bam":
		outputBAM = true
	case "":
		// Infer from output file extension. Stdout defaults to BAM (matches
		// upstream samtools when format is unspecified).
		if strings.HasSuffix(strings.ToLower(outFile), ".sam") {
			outputSAM = true
		} else {
			outputBAM = true
		}
	default:
		fmt.Fprintf(os.Stderr, "samtools sort: unknown -O format %q\n", outFmt)
		return 2
	}

	opts := samtools.SortOptions{
		Order:         order,
		Tag:           byTag,
		MaxMemBytes:   mem,
		TmpPrefix:     tmpdir,
		CompressLevel: compLevel,
		OutputBAM:     outputBAM,
		OutputSAM:     outputSAM,
		Threads:       threads,
		NoPG:          noPG,
	}

	// With -@ >= 2 the input is opened raw (still BGZF-framed) so samtools.Sort
	// can inflate its blocks in parallel; otherwise the standard decompressing
	// opener is used. The decoded records — and thus the sorted output — are
	// identical either way; only input-decode throughput changes.
	var in io.ReadCloser
	if threads >= 2 {
		in, err = iohelper.OpenRaw(fs.Arg(0))
	} else {
		in, err = iohelper.OpenReader(fs.Arg(0))
	}
	if err != nil {
		fmt.Fprintf(os.Stderr, "samtools sort: %v\n", err)
		return 1
	}
	defer in.Close()
	out, err := openOut(outFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "samtools sort: %v\n", err)
		return 1
	}
	defer out.Close()
	if err := samtools.Sort(in, out, opts); err != nil {
		fmt.Fprintf(os.Stderr, "samtools sort: %v\n", err)
		return 1
	}
	return 0
}

const indexUsage = `samtools index - build an index for a coordinate-sorted BAM or CRAM.

Usage:
  samtools index [options] <in.sorted.bam|in.cram>

Options:
  -b, --bai                   Produce a .bai index (BAM input; default).
  -c, --csi                   Produce a .csi index (required for reference
                              sequences longer than ~512 Mbp).
  -m, --min-shift N           CSI bin-hierarchy min_shift (default 14;
                              only used with -c). Accepts --csi-min-shift
                              as an alias for upstream compatibility.
  -o, --output PATH           Index output path. Default is <in>.bai (or
                              <in>.csi with -c) for BAM input, or
                              <in>.crai for CRAM input.
  -@, --threads N             Accepted; v1 is single-threaded.
  -h, --help                  Show this help.
  -v, --version               Show version.

The index kind is chosen from the input format and options: a CRAM file
is given a .crai index, a BAM file a .bai index (or a .csi index with -c).
`

func runIndex(args []string) int {
	fs := flag.NewFlagSet("samtools index", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	var (
		wantBAI     bool
		wantCSI     bool
		csiMinShift int
		outFile     string
		threads     int
		showHelp    bool
		showVer     bool
	)
	cliflag.BoolVar(fs, &wantBAI, "b", "bai", false, "Emit .bai")
	cliflag.BoolVar(fs, &wantCSI, "c", "csi", false, "Emit .csi")
	cliflag.IntVar(fs, &csiMinShift, "m", "min-shift", 14, "CSI bin-hierarchy min_shift")
	// --csi-min-shift is accepted as an alias for upstream compatibility.
	fs.IntVar(&csiMinShift, "csi-min-shift", 14, "")
	cliflag.StringVar(fs, &outFile, "o", "output", "", "Output path")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "Threads (accepted, ignored)")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "version", false, "")
	// Upstream index (bam_index.c getopt "bcm:Mo:@:") accepts -M to index
	// multiple input files at once. This port indexes a single file; -M is
	// accepted as a no-op for compatibility so bundled clusters parse.
	var indexMulti bool
	cliflag.BoolVar(fs, &indexMulti, "M", "", false, "Multiple inputs (accepted; single-file in v1)")

	if err := cliflag.Parse(fs, args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, indexUsage)
		return 2
	}
	_ = indexMulti
	if showHelp {
		fmt.Print(indexUsage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}
	if fs.NArg() == 0 {
		fmt.Fprintln(os.Stderr, "samtools index: missing input file")
		fmt.Fprint(os.Stderr, indexUsage)
		return 2
	}
	opts := samtools.IndexOptions{
		SelectCSI:   wantCSI,
		CSIMinShift: csiMinShift,
		Threads:     threads,
	}
	if err := samtools.IndexFile(fs.Arg(0), outFile, opts); err != nil {
		fmt.Fprintf(os.Stderr, "samtools index: %v\n", err)
		return 1
	}
	return 0
}

func openOut(path string) (io.WriteCloser, error) {
	if path == "" || path == "-" {
		return &nopCloser{os.Stdout}, nil
	}
	return os.Create(path)
}

type nopCloser struct {
	io.Writer
}

func (nopCloser) Close() error { return nil }

const depthUsage = `samtools depth - print per-position depth across one or more BAMs.

Usage:
  samtools depth [options] <in1.bam> [<in2.bam> ...]

Options:
  -a, --all                 Emit positions with 0 depth too (within covered
                            ranges). Repeatable: -aa (or -a -a) extends to every
                            position of every reference, like -A.
  -A, --all-trans           Emit every position of every reference.
  -r, --region chr[:S-E]    Limit to region (chr name only or range).
  -b, --bed FILE            Limit to BED regions.
  -q, --min-BQ N            Skip bases with base quality below N (default 0).
  -Q, --min-MQ N            Skip reads with mapping quality below N (default 0).
  -l, --min-readlen N       Skip reads shorter than N query bases (default 0).
  -f, --include-flags N     Require ALL these flag bits set (default 0).
  -F, --exclude-flags N     Drop reads with ANY of these flag bits set (default 0x4).
  -d, --max-depth N         Cap reported depth (0 = no cap).
  -o, --output PATH         Output path (default stdout).
  -@, --threads N           Accepted; single-threaded in v1.
  -h, --help                Show this help.
  -v, --version             Show version.
`

func runDepth(args []string) int {
	fs := flag.NewFlagSet("samtools depth", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	var (
		allCount countFlag
		allLong  bool
		allTrans bool
		regions  multiString
		bedPath  string
		minMAPQ  int
		minBaseQ int
		minReadL int
		incFlags int
		excFlags int
		maxDepth int
		outPath  string
		threads  int
		showHelp bool
		showVer  bool
	)
	// Upstream depth treats -a as repeatable (bam2depth.c `opt.all_pos++`):
	// one -a emits zero-depth positions inside covered regions; a second
	// (-aa, or -a -a) extends to every reference position. Registering the
	// short -a as a count flag lets cliflag.Parse expand the fused `-aa` to
	// `-a -a` with the same meaning. The long --all spelling keeps the
	// single-step behaviour; --all-trans / -A still force the all-references
	// mode.
	cliflag.Var(fs, &allCount, "a", "", "Emit zero-depth positions (repeat for all reference positions)")
	cliflag.BoolVar(fs, &allLong, "", "all", false, "Emit zero-depth positions")
	cliflag.BoolVar(fs, &allTrans, "A", "all-trans", false, "Emit every reference position")
	cliflag.Var(fs, &regions, "r", "region", "")
	cliflag.StringVar(fs, &bedPath, "b", "bed", "", "BED of regions")
	// Upstream depth maps -q/--min-BQ to the BASE-quality threshold and
	// -Q/--min-MQ to the MAPPING-quality threshold (bam2depth.c: case 'q' ->
	// opt.min_qual, the per-base floor; case 'Q' -> opt.min_mqual, the read
	// MAPQ floor). An earlier port had these reversed. Register the upstream
	// long spellings (--min-BQ / --min-MQ) plus the legacy --min-baseq /
	// --min-mapq aliases so older command lines still parse to the right knob.
	cliflag.IntVar(fs, &minBaseQ, "q", "min-BQ", 0, "Min base quality")
	cliflag.IntVar(fs, &minMAPQ, "Q", "min-MQ", 0, "Min MAPQ")
	cliflag.IntVar(fs, &minBaseQ, "", "min-baseq", 0, "Min base quality (alias of --min-BQ)")
	cliflag.IntVar(fs, &minMAPQ, "", "min-mapq", 0, "Min MAPQ (alias of --min-MQ)")
	cliflag.IntVar(fs, &minReadL, "l", "min-readlen", 0, "Min read length")
	cliflag.IntVar(fs, &incFlags, "f", "include-flags", 0, "Required flags")
	cliflag.IntVar(fs, &excFlags, "F", "exclude-flags", int(samtools.DefaultDepthExcludeFlags), "Excluded flags")
	cliflag.IntVar(fs, &maxDepth, "d", "max-depth", 0, "Max depth cap")
	cliflag.StringVar(fs, &outPath, "o", "output", "", "Output path")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "Threads (accepted, ignored)")
	fs.BoolVar(&showHelp, "h", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "v", false, "")
	fs.BoolVar(&showVer, "version", false, "")
	// Upstream depth (bam2depth.c getopt "@:q:Q:JHd:m:l:g:G:o:ar:Xf:b:s")
	// exposes a few short flags this port does not act on. Register them as
	// accepted stubs (matching upstream where each is a no-op or unsupported
	// here) so legacy command lines and bundled clusters parse:
	//   -J  include deletions in the depth (upstream default is to skip them;
	//       this port already counts coverage spans, so accepted no-op)
	//   -s  remove overlapping mate bases (accepted no-op in v1)
	//   -X  expect an explicit index-file argument (accepted no-op)
	//   -H  emit a column header line (accepted no-op)
	var (
		depthInclDel   bool
		depthRmOvl     bool
		depthCustomIdx bool
		depthHeader    bool
	)
	fs.BoolVar(&depthInclDel, "J", false, "")
	fs.BoolVar(&depthRmOvl, "s", false, "")
	fs.BoolVar(&depthCustomIdx, "X", false, "")
	fs.BoolVar(&depthHeader, "H", false, "")

	if err := cliflag.Parse(fs, args); err != nil {
		if err == flag.ErrHelp {
			fmt.Print(depthUsage)
			return 0
		}
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, depthUsage)
		return 2
	}
	_ = depthInclDel
	_ = depthRmOvl
	_ = depthCustomIdx
	_ = depthHeader
	if showHelp {
		fmt.Print(depthUsage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}
	if fs.NArg() == 0 {
		fmt.Fprintln(os.Stderr, "samtools depth: missing input file")
		fmt.Fprint(os.Stderr, depthUsage)
		return 2
	}

	// Resolve the -a count, --all long flag, and -A/--all-trans into the two
	// booleans the engine uses: any -a (or --all) → all positions inside
	// covered regions; a second -a (i.e. -aa) or -A/--all-trans → every
	// reference position.
	allPos := allCount >= 1 || allLong
	if allCount >= 2 {
		allTrans = true
	}

	// Resolve the effective inflate worker count once: a default (no -@)
	// opts into parallel BGZF decode across cores. The same resolved value
	// gates the raw open below and is handed to samtools.Depth.
	effThreads := samtools.ReadDecodeThreads(threads)

	opts := samtools.DepthOptions{
		AllPositions:      allPos,
		AllTransPositions: allTrans,
		Regions:           []string(regions),
		BedPath:           bedPath,
		MinMAPQ:           uint8(minMAPQ),
		MinBaseQ:          uint8(minBaseQ),
		MinReadLen:        minReadL,
		IncludeFlags:      uint16(incFlags),
		ExcludeFlags:      uint16(excFlags),
		MaxDepth:          maxDepth,
		Threads:           effThreads,
	}

	paths := fs.Args()
	readers := make([]io.Reader, 0, len(paths))
	closers := make([]io.Closer, 0, len(paths))
	for _, path := range paths {
		// When parallel decode is in effect (resolved threads >= 2) open the
		// file raw so samtools.Depth can inflate the BGZF blocks in parallel;
		// otherwise use the decompressing opener. The per-position depth output
		// is identical for both paths.
		var (
			r   io.ReadCloser
			err error
		)
		if effThreads >= 2 {
			r, err = iohelper.OpenRaw(path)
		} else {
			r, err = iohelper.OpenReader(path)
		}
		if err != nil {
			for _, c := range closers {
				_ = c.Close()
			}
			fmt.Fprintf(os.Stderr, "samtools depth: %v\n", err)
			return 1
		}
		readers = append(readers, r)
		closers = append(closers, r)
	}
	out, err := openOut(outPath)
	if err != nil {
		for _, c := range closers {
			_ = c.Close()
		}
		fmt.Fprintf(os.Stderr, "samtools depth: %v\n", err)
		return 1
	}
	defer out.Close()

	// With an -r region query against on-disk inputs, DepthFile seeks straight
	// to the indexed BGZF chunks (a sibling .csi / .bai) so only the region's
	// blocks are inflated, falling back to the linear scan over the readers we
	// already opened when no index is usable. The streaming Depth handles every
	// other case (no region, stdin, -b BED, -A). DepthFile owns its own seek
	// handles; the readers above remain the fallback source it reuses.
	if err := samtools.DepthFile(paths, readers, out, opts); err != nil {
		for _, c := range closers {
			_ = c.Close()
		}
		fmt.Fprintf(os.Stderr, "samtools depth: %v\n", err)
		return 1
	}
	for _, c := range closers {
		_ = c.Close()
	}
	return 0
}

// multiString implements flag.Value for collecting repeated string flags
// (used by `-r`/`--region`, which may appear multiple times).
type multiString []string

func (m *multiString) String() string     { return strings.Join(*m, ",") }
func (m *multiString) Set(v string) error { *m = append(*m, v); return nil }

// parseOutputFmtOptions interprets the KEY=VALUE strings collected from
// `--output-fmt-option`. It returns the raw qbin value (passed on to
// alnio.ParseQualityBinning, which validates it) and an error for any
// malformed entry or unrecognised key. A later entry for the same key
// overrides an earlier one, matching upstream samtools' last-wins
// behaviour.
func parseOutputFmtOptions(opts multiString) (qbin string, err error) {
	for _, o := range opts {
		eq := strings.IndexByte(o, '=')
		if eq < 1 {
			return "", fmt.Errorf("samtools view: malformed --output-fmt-option %q (want KEY=VALUE)", o)
		}
		key, val := o[:eq], o[eq+1:]
		switch key {
		case "qbin", "quality-binning":
			qbin = val
		default:
			return "", fmt.Errorf("samtools view: unknown --output-fmt-option key %q (supported: qbin)", key)
		}
	}
	return qbin, nil
}

const fastqUsage = `samtools fastq - convert SAM/BAM to FASTQ.

Usage:
  samtools fastq [options] <in.bam|in.sam>

Options:
  -1, --read1 FILE          Output for read1 (mate 1).
  -2, --read2 FILE          Output for read2 (mate 2).
  -0, --read-orphan FILE    Reads where 0x40/0x80 are both set or both unset.
  -s, --singleton FILE      Output for unpaired reads.
  -o, --output FILE         Default sink (interleaved if -1/-2 unset).
  -N, --output-name         Always append /1 or /2 to read names.
  -n, --no-suffix           Never append /1 /2.
  -f, --include-flags N     Required flag bits.
  -F, --exclude-flags N     Excluded flag bits (default 0x900).
  -G, --exclude-flags-all N Drop only when ALL bits match.
  -T, --add-tags TAGS       Comma-separated aux tags to append.
  -t, --no-CO               Accepted for compatibility (no-op).
  -c, --compress-level N    Gzip level for .gz outputs.
  -O, --use-qq              Use OQ aux tag for quality when present.
      --threads N           Accepted; single-threaded.
  -h, --help                Show this help.
  -v, --version             Show version.

Paired output (-1/-2) requires name-sorted input. Coordinate-sorted input
falls back to writing every record through to -o (or singletons).
`

func runFastq(args []string) int {
	fs := flag.NewFlagSet("samtools fastq", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	var (
		r1Path    string
		r2Path    string
		orphPath  string
		singPath  string
		outPath   string
		alwaysSfx bool
		noSuffix  bool
		incFlags  int
		excFlags  int
		excFlagsG int
		addTags   string
		noCO      bool
		compLevel string
		useOQ     bool
		threads   int
		showHelp  bool
		showVer   bool
	)
	cliflag.StringVar(fs, &r1Path, "1", "read1", "", "Read1 output")
	cliflag.StringVar(fs, &r2Path, "2", "read2", "", "Read2 output")
	cliflag.StringVar(fs, &orphPath, "0", "read-orphan", "", "Orphan output")
	cliflag.StringVar(fs, &singPath, "s", "singleton", "", "Singleton output")
	cliflag.StringVar(fs, &outPath, "o", "output", "", "Default output")
	cliflag.BoolVar(fs, &alwaysSfx, "N", "output-name", false, "Always add /1 /2 suffix")
	cliflag.BoolVar(fs, &noSuffix, "n", "no-suffix", false, "Never add /1 /2 suffix")
	cliflag.IntVar(fs, &incFlags, "f", "include-flags", 0, "Required flags")
	cliflag.IntVar(fs, &excFlags, "F", "exclude-flags", int(sam.FlagSecondary|sam.FlagSupplementary), "Excluded flags")
	cliflag.IntVar(fs, &excFlagsG, "G", "exclude-flags-all", 0, "Drop if all of these set")
	cliflag.StringVar(fs, &addTags, "T", "add-tags", "", "Aux tags to append (CSV)")
	cliflag.BoolVar(fs, &noCO, "t", "no-CO", false, "No @CO emission (no-op)")
	cliflag.StringVar(fs, &compLevel, "c", "compress-level", "", "Gzip level for .gz outputs")
	cliflag.BoolVar(fs, &useOQ, "O", "use-qq", false, "Use OQ aux tag for quality")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "Threads (accepted, ignored)")
	fs.BoolVar(&showHelp, "h", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "v", false, "")
	fs.BoolVar(&showVer, "version", false, "")
	// Upstream fastq (bam_fastq.c getopt "0:1:2:o:f:F:G:niNOs:c:tT:v:@:d:D:")
	// exposes a few short flags this port does not act on. Register them as
	// accepted stubs so legacy command lines and bundled clusters parse:
	//   -i        emit the Illumina-style index/casava tag (accepted no-op)
	//   -d STR[:VAL]  keep records with this aux tag (accepted no-op)
	//   -D STR:FILE   keep records with aux tag in FILE (accepted no-op)
	// (-v is intentionally this port's version switch rather than upstream's
	// default-quality value; -@ is the canonical threads short flag.)
	var (
		fqIllumina bool
		fqTagSpec  string
		fqTagFile  string
	)
	fs.BoolVar(&fqIllumina, "i", false, "")
	fs.StringVar(&fqTagSpec, "d", "", "")
	fs.StringVar(&fqTagFile, "D", "", "")

	if err := cliflag.Parse(fs, args); err != nil {
		if err == flag.ErrHelp {
			fmt.Print(fastqUsage)
			return 0
		}
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, fastqUsage)
		return 2
	}
	_ = fqIllumina
	_ = fqTagSpec
	_ = fqTagFile
	if showHelp {
		fmt.Print(fastqUsage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}
	if fs.NArg() == 0 {
		fmt.Fprintln(os.Stderr, "samtools fastq: missing input file")
		fmt.Fprint(os.Stderr, fastqUsage)
		return 2
	}
	if r1Path == "" && r2Path == "" && outPath == "" && singPath == "" && orphPath == "" {
		// Default to interleaved stdout.
		outPath = "-"
	}
	level, err := samtools.ParseCompressLevel(compLevel)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}
	opts := samtools.FastqOptions{
		Read1Path:       r1Path,
		Read2Path:       r2Path,
		OrphanPath:      orphPath,
		SingletonPath:   singPath,
		OutputPath:      outPath,
		AlwaysAddSuffix: alwaysSfx,
		NoSuffix:        noSuffix,
		IncludeFlags:    uint16(incFlags),
		ExcludeFlags:    uint16(excFlags),
		ExcludeFlagsAll: uint16(excFlagsG),
		UseExcludeAll:   excFlagsG != 0,
		AddTags:         samtools.ParseAddTags(addTags),
		CompressLevel:   level,
		UseOQ:           useOQ,
		NoCO:            noCO,
		Threads:         threads,
	}
	in, err := iohelper.OpenReader(fs.Arg(0))
	if err != nil {
		fmt.Fprintf(os.Stderr, "samtools fastq: %v\n", err)
		return 1
	}
	defer in.Close()
	counts, err := samtools.Fastq(in, opts)
	if err != nil {
		fmt.Fprintf(os.Stderr, "samtools fastq: %v\n", err)
		return 1
	}
	if counts.PairedCoordinateWarn {
		fmt.Fprintln(os.Stderr, "samtools fastq: paired output (-1/-2) requires name-sorted input; coordinate-sorted input falls back to -o/singleton")
	}
	return 0
}

const mpileupUsage = `samtools mpileup - per-position pileup across one or more BAMs.

Usage:
  samtools mpileup [options] <in1.bam> [<in2.bam> ...]

Options:
  -f, --fasta-ref FASTA      Reference FASTA (random-access via .fai).
  -l, --positions FILE       BED or 2-col positions file restricting output.
  -r, --region chr:start-end Restrict to region (uses .bai when available).
  -b, --bam-list FILE        File of BAM paths, one per line.
  -q, --min-mapq INT         Min MAPQ (default 0).
  -Q, --min-baseq INT        Min base quality (default 13).
  -d, --max-depth INT        Max reads per position (default 8000).
  -A, --count-orphans        Include reads with unmapped mates / anomalous pairs.
  -x, --ignore-overlaps      Discard overlapping mate-pair bases.
  -E, --redo-baq             Recompute BAQ on the fly, ignoring any existing BQ tag.
  -B, --no-BAQ               Disable BAQ (per-Base Alignment Quality) realignment.
  -a, --all-positions        Emit zero-depth positions inside covered regions.
      --all-positions-all-chroms
                             Emit every reference position (-aa). Pass -a twice or this flag.
  -s, --output-mapq          Append MAPQs column.
  -O, --output-BP            Append per-read positions column.
  -o, --output PATH          Output file (default stdout).
  -u, --uncompressed-bcf     Uncompressed BCF genotype-likelihood output.
  -g, --bcf                  BCF genotype-likelihood output (FORMAT/PL).
  -@, --threads N            Accepted; single-threaded in v1.
  -h, --help                 Show this help.
  -v, --version              Show version.

Notes:
  - Multi-BAM input emits parallel "depth bases quals" triples per BAM.
  - The bases column encodes . / , (forward/reverse strand match), upper/lower
    base for mismatch, * for deletion/refskip placeholder, +<len><seq> for
    insertions after the base, -<len><seq> for deletions starting after this
    position, ^<charq> for read start (charq = mapq + 33), $ for read end.
  - -g/-u request BCF (and uncompressed-BCF) genotype-likelihood output
    instead of the text pileup. Upstream samtools removed this path (1.10+)
    in favour of "bcftools mpileup"; this port restores it by delegating to
    the shared bcftools mpileup engine, emitting per-site FORMAT/PL with the
    <*> unseen allele. -f/--fasta-ref is required. Full indel calling is the
    deferred remainder tracked in docs/PARITY_ROADMAP.md#samtools.
  - -E (--redo-baq) is honoured by the -g/-u path; in the text path it is
    deferred per docs/PARITY_ROADMAP.md#samtools.
`

func runMpileup(args []string) int {
	// mpileup over a whole chromosome keeps a lean but non-trivial live set
	// (the contig reference plus a bounded sliding pileup-event matrix). Under
	// the default GOGC=100 the collector lets the heap grow to ~2x that live set
	// before collecting, so the resident high-water mark — what ru_maxrss
	// reports — overshoots well past the working set. Rather than trimming GOGC
	// (which makes the collector run aggressively on *every* cycle, taxing CPU
	// even when there is ample headroom), we install a soft memory limit that
	// bounds total Go-managed memory — the quantity ru_maxrss tracks, including
	// the off-heap BGZF read buffers. Below the cap the GC stays lazy (default
	// GOGC) so wall time is unaffected; it only gets more frequent as the heap
	// approaches the ceiling. The limit is soft, so resident memory can overshoot
	// it transiently by a few MB (off-heap BGZF buffers, runtime overhead, a
	// jittery high-water sample); 120 MiB leaves margin so even the worst-of-five
	// peak on whole-chromosome -r stays comfortably within ~2x of upstream's
	// (~74 MB), while keeping ~40 MB of headroom over the live working set so the
	// GC stays lazy and wall time is unaffected — all measured on real GIAB data.
	// SetMemoryLimit returns the previous limit (the no-limit math.MaxInt64
	// default when unset), restored on return so this stays scoped to the mpileup
	// command.
	prevLimit := debug.SetMemoryLimit(120 << 20)
	defer debug.SetMemoryLimit(prevLimit)

	fs := flag.NewFlagSet("samtools mpileup", flag.ContinueOnError)
	fs.SetOutput(io.Discard)
	var (
		fastaRef  string
		posFile   string
		regions   multiString
		bamList   string
		minMAPQ   int
		minBaseQ  int
		maxDepth  int
		orphans   bool
		ignoreOvl bool
		redoBAQ   bool
		noBAQ     bool
		allPosCnt countFlag
		allChrom  bool
		outMapq   bool
		outBP     bool
		outPath   string
		bcf       bool
		ubcf      bool
		threads   int
		showHelp  bool
		showVer   bool
	)
	cliflag.StringVar(fs, &fastaRef, "f", "fasta-ref", "", "Reference FASTA")
	cliflag.StringVar(fs, &posFile, "l", "positions", "", "Positions/BED file")
	cliflag.Var(fs, &regions, "r", "region", "")
	cliflag.StringVar(fs, &bamList, "b", "bam-list", "", "BAM list file")
	cliflag.IntVar(fs, &minMAPQ, "q", "min-mapq", 0, "Min MAPQ")
	cliflag.IntVar(fs, &minBaseQ, "Q", "min-baseq", samtools.DefaultMpileupMinBaseQ, "Min base quality")
	cliflag.IntVar(fs, &maxDepth, "d", "max-depth", samtools.DefaultMpileupMaxDepth, "Max depth")
	cliflag.BoolVar(fs, &orphans, "A", "count-orphans", false, "Include orphan reads")
	cliflag.BoolVar(fs, &ignoreOvl, "x", "ignore-overlaps", false, "Discard overlapping mates")
	cliflag.BoolVar(fs, &redoBAQ, "E", "redo-baq", false, "Recompute BAQ, ignoring existing BQ tags")
	cliflag.BoolVar(fs, &noBAQ, "B", "no-BAQ", false, "Disable BAQ realignment")
	// -a is repeatable upstream (mplp.all++): one -a emits zero-depth
	// positions inside covered regions; a second (-aa, or -a -a) extends to
	// every reference position. A count flag lets cliflag.Parse expand the
	// bundled `-aa` to `-a -a` and have it mean the same thing — no special
	// pre-pass needed. The long --all-positions-all-chroms remains a
	// shortcut for "all twice".
	cliflag.Var(fs, &allPosCnt, "a", "all-positions", "Emit zero-depth positions (repeat for all reference positions)")
	cliflag.BoolVar(fs, &allChrom, "", "all-positions-all-chroms", false, "Emit every reference position (-aa)")
	cliflag.BoolVar(fs, &outMapq, "s", "output-mapq", false, "Append MAPQs column")
	cliflag.BoolVar(fs, &outBP, "O", "output-BP", false, "Append per-read positions column")
	cliflag.StringVar(fs, &outPath, "o", "output", "", "Output path")
	cliflag.BoolVar(fs, &ubcf, "u", "uncompressed-bcf", false, "Uncompressed BCF genotype-likelihood output")
	cliflag.BoolVar(fs, &bcf, "g", "bcf", false, "BCF genotype-likelihood output")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "Threads")
	// Upstream mpileup (bam_plcmd.c getopt "Af:r:l:q:Q:RC:Bd:b:o:EG:6OsxXaM")
	// has several short flags this port does not act on. Register them as
	// accepted stubs so legacy command lines and bundled clusters parse:
	//   -R        ignore read groups (accepted no-op)
	//   -6        input quals are Illumina-1.3+ encoded (accepted no-op)
	//   -M        print base modifications column (accepted no-op)
	//   -X        expect an explicit index-file argument (accepted no-op)
	//   -C N      mapping-quality adjustment coefficient (accepted no-op)
	//   -G FILE   exclude read groups listed in FILE (accepted no-op)
	var (
		mpIgnoreRG   bool
		mpIllumina   bool
		mpPrintMods  bool
		mpCustomIdx  bool
		mpCapQ       int
		mpExclRGFile string
	)
	fs.BoolVar(&mpIgnoreRG, "R", false, "")
	fs.BoolVar(&mpIllumina, "6", false, "")
	fs.BoolVar(&mpPrintMods, "M", false, "")
	fs.BoolVar(&mpCustomIdx, "X", false, "")
	fs.IntVar(&mpCapQ, "C", 0, "")
	fs.StringVar(&mpExclRGFile, "G", "", "")
	fs.BoolVar(&showHelp, "h", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "v", false, "")
	fs.BoolVar(&showVer, "version", false, "")

	// Route through cliflag.Parse so POSIX getopt-style short-flag bundling
	// works the way upstream samtools' getopt parser accepts them — including
	// the fused `-aa` ("all positions, all chromosomes"), which the -a count
	// flag above turns into two increments.
	if err := cliflag.Parse(fs, args); err != nil {
		if err == flag.ErrHelp {
			fmt.Print(mpileupUsage)
			return 0
		}
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, mpileupUsage)
		return 2
	}
	_ = mpIgnoreRG
	_ = mpIllumina
	_ = mpPrintMods
	_ = mpCustomIdx
	_ = mpCapQ
	_ = mpExclRGFile
	// Resolve the repeatable -a count into the two booleans the engine uses:
	// one -a → all positions in covered range; two or more (or
	// --all-positions-all-chroms) → every reference position.
	allPos := allPosCnt >= 1
	if allPosCnt >= 2 {
		allChrom = true
	}
	if showHelp {
		fmt.Print(mpileupUsage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}
	if fs.NArg() == 0 && bamList == "" {
		fmt.Fprintln(os.Stderr, "samtools mpileup: missing input file")
		fmt.Fprint(os.Stderr, mpileupUsage)
		return 2
	}

	// -g/-u: BCF / uncompressed-BCF genotype-likelihood output. Upstream
	// removed this from `samtools mpileup` (1.10+) in favour of `bcftools
	// mpileup`; we restore it by delegating to the shared bcftools mpileup
	// engine (the htslib bam2bcf port). See tools/samtools/pkg/samtools/
	// mpileup_bcf.go. The samtools text-path defaults for -Q (13) and -d
	// (8000) do not apply to the genotype-likelihood caller; only honour
	// these flags when the user set them explicitly, letting the bcftools
	// engine fill its own upstream defaults (min-BQ 1, max-depth 250)
	// otherwise.
	if bcf || ubcf {
		bopts := samtools.MpileupBCFOptions{
			Inputs:         fs.Args(),
			FastaRef:       fastaRef,
			BamList:        bamList,
			Regions:        []string(regions),
			PositionsFile:  posFile,
			MinMAPQ:        uint8(minMAPQ),
			CountOrphans:   orphans,
			IgnoreOverlaps: ignoreOvl,
			NoBAQ:          noBAQ,
			RedoBAQ:        redoBAQ,
			Uncompressed:   ubcf,
			Output:         outPath,
			Threads:        threads,
		}
		fs.Visit(func(f *flag.Flag) {
			switch f.Name {
			case "Q", "min-baseq":
				bopts.MinBaseQ = uint8(minBaseQ)
			case "d", "max-depth":
				bopts.MaxDepth = maxDepth
			}
		})
		out, err := openOut(outPath)
		if err != nil {
			fmt.Fprintf(os.Stderr, "samtools mpileup: %v\n", err)
			return 1
		}
		defer out.Close()
		if err := samtools.MpileupBCF(bopts, out); err != nil {
			fmt.Fprintf(os.Stderr, "%v\n", err)
			return 1
		}
		return 0
	}
	opts := samtools.MpileupOptions{
		Inputs:                fs.Args(),
		FastaRef:              fastaRef,
		Regions:               []string(regions),
		PositionsFile:         posFile,
		MinMAPQ:               uint8(minMAPQ),
		MinBaseQ:              uint8(minBaseQ),
		MaxDepth:              maxDepth,
		CountOrphans:          orphans,
		IgnoreOverlaps:        ignoreOvl,
		AllPositions:          allPos,
		AllPositionsAllChroms: allChrom,
		OutputMapQ:            outMapq,
		OutputBP:              outBP,
		NoBAQ:                 noBAQ,
		RedoBAQ:               redoBAQ,
		Output:                outPath,
		Threads:               threads,
		BamList:               bamList,
	}

	out, err := openOut(outPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "samtools mpileup: %v\n", err)
		return 1
	}
	defer out.Close()
	if err := samtools.MpileupFile(opts, out); err != nil {
		fmt.Fprintf(os.Stderr, "samtools mpileup: %v\n", err)
		return 1
	}
	return 0
}
