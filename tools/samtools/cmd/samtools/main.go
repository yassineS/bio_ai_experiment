// Command samtools is a pure-Go reimplementation of selected samtools
// subcommands. This first slice ships `view` and `flagstat`; other
// subcommands (sort, index, depth, fastq, mpileup) will follow in later PRs.
package main

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

const version = "0.1.0"

const rootUsage = `samtools - pure-Go reimplementation (first slice).

Usage:
  samtools <subcommand> [options]

Subcommands:
  view       Print, filter, or convert SAM/BAM records.
  flagstat   Print flag statistics for a SAM/BAM file.
  help       Show this help.
  version    Show version.
`

func main() {
	if len(os.Args) < 2 {
		fmt.Fprint(os.Stderr, rootUsage)
		os.Exit(1)
	}
	switch os.Args[1] {
	case "view":
		os.Exit(runView(os.Args[2:]))
	case "flagstat":
		os.Exit(runFlagstat(os.Args[2:]))
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
  -h, --with-header           Include the header in SAM output.
  -H, --header-only           Print the header only.
  -c, --count                 Print only the count of matching records.
  -f, --include-flags <int>   Require ALL bits set.
  -F, --exclude-flags <int>   Drop records with ANY bit set.
  -G, --exclude-flags-all <int>  Drop only when ALL bits set.
  -q, --min-mapq <int>        Minimum MAPQ.
  -r, --read-group <string>   Keep records with this RG.
  -R, --read-groups-file <f>  File of RG IDs (one per line).
  -L, --regions-file <f>      BED of regions (deferred — see notes).
  -s, --subsample <float>     Keep fraction (or "<seed>.<frac>").
  -o, --output <file>         Output file (default stdout).
  -T, --reference <fasta>     Accepted; CRAM is not supported in v1.
  -@, --threads <int>         Accepted; single-threaded in v1.
      --no-PG                 Suppress @PG line emission.
      --help                  Show this help.
      --version               Show version.

Region queries (chr:start-end and -L) require BAI indexing, which is not
yet implemented. Providing a region will cause an explicit error.
`

func runView(args []string) int {
	fs := flag.NewFlagSet("samtools view", flag.ContinueOnError)
	fs.SetOutput(io.Discard) // we print usage ourselves

	var (
		outBAM    bool
		withHdr   bool
		hdrOnly   bool
		countOnly bool
		incFlags  int
		excFlags  int
		excFlagsG int
		minMAPQ   int
		rg        string
		rgFile    string
		regFile   string
		subsample string
		outFile   string
		refFile   string
		threads   int
		noPG      bool
		showHelp  bool
		showVer   bool
	)
	cliflag.BoolVar(fs, &outBAM, "b", "bam", false, "Output BAM")
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
	cliflag.StringVar(fs, &subsample, "s", "subsample", "", "Subsample fraction")
	cliflag.StringVar(fs, &outFile, "o", "output", "", "Output file")
	cliflag.StringVar(fs, &refFile, "T", "reference", "", "Reference FASTA (CRAM unsupported)")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "Threads (accepted, ignored)")
	cliflag.BoolVar(fs, &noPG, "", "no-PG", false, "Suppress @PG line")
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
		fmt.Fprintln(os.Stderr, "samtools view: missing input file")
		fmt.Fprint(os.Stderr, viewUsage)
		return 2
	}
	input := positional[0]
	regions := positional[1:]

	opts := samtools.ViewOptions{
		OutputBAM:       outBAM,
		WithHeader:      withHdr,
		HeaderOnly:      hdrOnly,
		Count:           countOnly,
		IncludeFlags:    uint16(incFlags),
		ExcludeFlags:    uint16(excFlags),
		ExcludeFlagsAll: uint16(excFlagsG),
		UseExcludeAll:   excFlagsG != 0,
		MinMAPQ:         uint8(minMAPQ),
		ReadGroup:       rg,
		Regions:         append([]string{}, regions...),
		RegionsEnabled:  len(regions) > 0 || regFile != "",
		NoPG:            noPG,
	}

	// Honour the .bam output extension even without explicit -b.
	if !opts.OutputBAM && outFile != "" && strings.HasSuffix(outFile, ".bam") {
		opts.OutputBAM = true
	}

	if rgFile != "" {
		set, err := samtools.LoadReadGroupsFile(rgFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "samtools view: %v\n", err)
			return 1
		}
		opts.ReadGroupSet = set
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

	in, err := iohelper.OpenReader(input)
	if err != nil {
		fmt.Fprintf(os.Stderr, "samtools view: %v\n", err)
		return 1
	}
	defer in.Close()

	out, err := openOut(outFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "samtools view: %v\n", err)
		return 1
	}
	defer out.Close()

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
	if err := fs.Parse(args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, flagstatUsage)
		return 2
	}
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
	in, err := iohelper.OpenReader(fs.Arg(0))
	if err != nil {
		fmt.Fprintf(os.Stderr, "samtools flagstat: %v\n", err)
		return 1
	}
	defer in.Close()
	if err := samtools.Flagstat(in, os.Stdout); err != nil {
		fmt.Fprintf(os.Stderr, "samtools flagstat: %v\n", err)
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
