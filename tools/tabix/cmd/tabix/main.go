// Command tabix builds a `.tbi` index for a bgzipped tab-delimited file
// (VCF, BED, GFF, SAM, custom) and supports random-access region queries
// of the indexed file.
//
// Usage:
//
//	tabix [options] file.gz                  # build index → file.gz.tbi
//	tabix [options] file.gz REGION [REGION...]   # query records by region
//
// Use --preset {vcf|bed|gff|sam} to select the standard column layouts; or
// specify -s/-b/-e/-S/-c/-0 manually.
package main

import (
	"bufio"
	"flag"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/tabix"
)

const version = "1.0.0"

const usage = `tabix - Generic random-access index for bgzipped TAB-delimited files.

Usage:
  tabix [options] file.gz                       # build index
  tabix [options] file.gz REGION [REGION...]    # query

A REGION is CHROM, CHROM:START-END, or CHROM:START (all 1-based, inclusive).

Options:
  -p, --preset {gff|bed|sam|vcf}   Select a standard column layout.
  -s, --seq-col N                  1-based column of sequence (chrom) name.
  -b, --begin-col N                1-based column of begin position.
  -e, --end-col N                  1-based column of end (0 = use begin).
  -S, --skip-lines N               Header lines to skip (default 0).
  -c, --meta-char CHAR             Comment-line prefix (default '#').
  -0, --zero-based                 0-based half-open coordinates (BED-style).
  -f, --force                      Overwrite an existing .tbi index.
  -R, --regions FILE               Read regions from a BED-like file (index-jump).
  -T, --targets FILE               Stream records overlapping intervals in FILE
                                   (strict post-filter, distinct from -R).
  -r, --reheader FILE              Replace the header with the contents of FILE.
  -l, --list-chroms                Print chromosome names from the index.
  -h, --print-header               Also emit header lines when querying.
      --only-header                Emit only the header from the file.
  -D                               Do not save the index (only relevant for build).
      --help                       Show this help and exit.
  -v, --version                    Show version and exit.
`

func main() {
	os.Exit(run(os.Args[1:], os.Stdin, os.Stdout, os.Stderr))
}

type opts struct {
	preset      string
	seqCol      int
	begCol      int
	endCol      int
	skip        int
	metaCh      string
	zeroBased   bool
	force       bool
	regionsFile string
	targetsFile string
	reheader    string
	listChroms  bool
	printHdr    bool
	onlyHdr     bool
	noSaveIdx   bool
	csiOutput   bool
	csiMinShift int
	showHelp    bool
	showVersion bool
}

func run(args []string, _ io.Reader, stdout, stderr io.Writer) int {
	fs := flag.NewFlagSet("tabix", flag.ContinueOnError)
	fs.SetOutput(stderr)
	o := opts{}

	cliflag.StringVar(fs, &o.preset, "p", "preset", "", "preset {gff|bed|sam|vcf}")
	cliflag.IntVar(fs, &o.seqCol, "s", "seq-col", 0, "sequence column (1-based)")
	cliflag.IntVar(fs, &o.begCol, "b", "begin-col", 0, "begin column (1-based)")
	cliflag.IntVar(fs, &o.endCol, "e", "end-col", 0, "end column (1-based, 0 to use begin)")
	cliflag.IntVar(fs, &o.skip, "S", "skip-lines", 0, "header lines to skip")
	cliflag.StringVar(fs, &o.metaCh, "c", "meta-char", "#", "comment-line prefix")
	cliflag.BoolVar(fs, &o.zeroBased, "0", "zero-based", false, "0-based half-open coordinates")
	cliflag.BoolVar(fs, &o.force, "f", "force", false, "overwrite existing index")
	cliflag.StringVar(fs, &o.regionsFile, "R", "regions", "", "BED-like regions file")
	cliflag.StringVar(fs, &o.targetsFile, "T", "targets", "", "strict overlap post-filter from targets file")
	cliflag.StringVar(fs, &o.reheader, "r", "reheader", "", "replace the header with the contents of FILE")
	cliflag.BoolVar(fs, &o.listChroms, "l", "list-chroms", false, "list chromosomes in the index")
	cliflag.BoolVar(fs, &o.printHdr, "h", "print-header", false, "emit header lines")
	fs.BoolVar(&o.onlyHdr, "only-header", false, "emit only the header")
	fs.BoolVar(&o.noSaveIdx, "D", false, "do not save the index")
	fs.BoolVar(&o.csiOutput, "csi", false, "emit a CSI index instead of .tbi")
	fs.IntVar(&o.csiMinShift, "csi-min-shift", 0, "CSI min_shift parameter (default 14)")
	fs.BoolVar(&o.showHelp, "help", false, "show help")
	cliflag.BoolVar(fs, &o.showVersion, "v", "version", false, "show version")

	fs.Usage = func() { fmt.Fprint(stderr, usage) }
	if err := fs.Parse(args); err != nil {
		return 2
	}
	if o.showHelp {
		fmt.Fprint(stdout, usage)
		return 0
	}
	if o.showVersion {
		fmt.Fprintf(stdout, "tabix %s (Go implementation)\n", version)
		return 0
	}

	rest := fs.Args()
	if len(rest) < 1 {
		fmt.Fprintln(stderr, "tabix: a data file argument is required")
		return 2
	}
	dataPath := rest[0]
	regions := rest[1:]

	// Resolve configuration: preset overrides individual columns when one
	// is given, but individual flags then override the preset's values.
	cfg, err := resolveConfig(o)
	if err != nil {
		fmt.Fprintf(stderr, "tabix: %v\n", err)
		return 2
	}

	switch {
	case o.listChroms:
		return runListChroms(dataPath, stdout, stderr)
	case len(regions) > 0 || o.regionsFile != "" || o.targetsFile != "" || o.onlyHdr:
		// Query mode: any of positional regions, -R, -T, or --only-header
		// (mirrors upstream tabix's main dispatch condition).
		return runQuery(dataPath, regions, o, cfg, stdout, stderr)
	case o.reheader != "":
		// Reheader mode: replace the header and re-emit a bgzipped stream.
		return runReheader(dataPath, o, cfg, stdout, stderr)
	default:
		// Build mode.
		return runBuild(dataPath, cfg, o.force, o.noSaveIdx, stderr)
	}
}

// runReheader replaces the leading header of the bgzipped file at dataPath
// with the contents of o.reheader and writes a fresh bgzipped stream to
// stdout, mirroring upstream tabix's `--reheader` behavior.
func runReheader(dataPath string, o opts, cfg tabix.Config, stdout, stderr io.Writer) int {
	if err := tabix.Reheader(dataPath, o.reheader, byte(cfg.Meta), stdout); err != nil {
		fmt.Fprintf(stderr, "tabix: reheader: %v\n", err)
		return 1
	}
	return 0
}

func resolveConfig(o opts) (tabix.Config, error) {
	var cfg tabix.Config
	if o.preset != "" {
		p, err := tabix.PresetConfig(o.preset)
		if err != nil {
			return cfg, err
		}
		cfg = p
	} else {
		cfg.Meta = '#'
	}
	if o.seqCol > 0 {
		cfg.ColSeq = int32(o.seqCol)
	}
	if o.begCol > 0 {
		cfg.ColBeg = int32(o.begCol)
	}
	if o.endCol > 0 {
		cfg.ColEnd = int32(o.endCol)
	}
	if o.skip > 0 {
		cfg.Skip = int32(o.skip)
	}
	if o.metaCh != "" {
		cfg.Meta = int32(o.metaCh[0])
	}
	if o.zeroBased {
		cfg.Format |= tabix.FlagZeroBased
	}
	return cfg, nil
}

func runBuild(dataPath string, cfg tabix.Config, force, noSave bool, stderr io.Writer) int {
	if cfg.ColSeq < 1 || cfg.ColBeg < 1 {
		fmt.Fprintln(stderr, "tabix: missing --preset or -s/-b column specification")
		return 2
	}
	tbiPath := dataPath + ".tbi"
	if !force && !noSave {
		if _, err := os.Stat(tbiPath); err == nil {
			fmt.Fprintf(stderr, "tabix: %s already exists; use -f to overwrite\n", tbiPath)
			return 1
		}
	}
	idx, err := tabix.Build(dataPath, cfg)
	if err != nil {
		fmt.Fprintf(stderr, "tabix: build failed: %v\n", err)
		return 1
	}
	if noSave {
		return 0
	}
	if err := idx.WriteFile(tbiPath); err != nil {
		fmt.Fprintf(stderr, "tabix: write %s failed: %v\n", tbiPath, err)
		return 1
	}
	return 0
}

func runListChroms(dataPath string, stdout, stderr io.Writer) int {
	tbiPath := dataPath + ".tbi"
	idx, err := tabix.ReadFile(tbiPath)
	if err != nil {
		fmt.Fprintf(stderr, "tabix: %v\n", err)
		return 1
	}
	for _, n := range idx.Names {
		fmt.Fprintln(stdout, n)
	}
	return 0
}

type region struct {
	chrom    string
	beg, end int // 0-based half-open
}

func parseRegion(s string) (region, error) {
	// CHROM, CHROM:START, CHROM:START-END
	idx := strings.IndexByte(s, ':')
	if idx < 0 {
		return region{chrom: s, beg: 0, end: 1 << 30}, nil
	}
	chrom := s[:idx]
	rest := s[idx+1:]
	var begStr, endStr string
	if d := strings.IndexByte(rest, '-'); d >= 0 {
		begStr = rest[:d]
		endStr = rest[d+1:]
	} else {
		begStr = rest
		endStr = rest
	}
	beg, err := strconv.Atoi(begStr)
	if err != nil {
		return region{}, fmt.Errorf("invalid start in %q", s)
	}
	end := beg
	if endStr != "" {
		end, err = strconv.Atoi(endStr)
		if err != nil {
			return region{}, fmt.Errorf("invalid end in %q", s)
		}
	}
	// CLI uses 1-based, inclusive; convert.
	return region{chrom: chrom, beg: beg - 1, end: end}, nil
}

func runQuery(dataPath string, regions []string, o opts, _ tabix.Config, stdout, stderr io.Writer) int {
	tbiPath := dataPath + ".tbi"
	idx, err := tabix.ReadFile(tbiPath)
	if err != nil {
		fmt.Fprintf(stderr, "tabix: %v\n", err)
		return 1
	}

	if o.printHdr || o.onlyHdr {
		if err := emitHeader(dataPath, idx.Config, stdout); err != nil {
			fmt.Fprintf(stderr, "tabix: %v\n", err)
			return 1
		}
	}
	if o.onlyHdr {
		return 0
	}

	// Load the -T/--targets strict overlap filter, if requested.
	var targets *tabix.Targets
	if o.targetsFile != "" {
		targets, err = tabix.LoadTargets(o.targetsFile)
		if err != nil {
			fmt.Fprintf(stderr, "tabix: targets: %v\n", err)
			return 1
		}
	}

	// Collect query regions from --regions FILE and positional args.
	var rs []region
	if o.regionsFile != "" {
		f, err := os.Open(o.regionsFile)
		if err != nil {
			fmt.Fprintf(stderr, "tabix: %v\n", err)
			return 1
		}
		scan := bufio.NewScanner(f)
		for scan.Scan() {
			line := strings.TrimSpace(scan.Text())
			if line == "" || strings.HasPrefix(line, "#") {
				continue
			}
			fields := strings.Fields(line)
			if len(fields) >= 3 {
				beg, _ := strconv.Atoi(fields[1])
				end, _ := strconv.Atoi(fields[2])
				rs = append(rs, region{chrom: fields[0], beg: beg, end: end})
			} else if len(fields) == 1 {
				rs = append(rs, region{chrom: fields[0], beg: 0, end: 1 << 30})
			}
		}
		f.Close()
	}
	for _, r := range regions {
		reg, err := parseRegion(r)
		if err != nil {
			fmt.Fprintf(stderr, "tabix: %v\n", err)
			return 2
		}
		rs = append(rs, reg)
	}

	// When -T is given without any explicit region, upstream streams the
	// whole file (the "." region), letting the targets filter do the work.
	if len(rs) == 0 && targets != nil {
		for _, chrom := range idx.Chroms() {
			rs = append(rs, region{chrom: chrom, beg: 0, end: 1 << 30})
		}
	}

	for _, r := range rs {
		records, err := idx.QueryRecords(dataPath, r.chrom, r.beg, r.end)
		if err != nil {
			fmt.Fprintf(stderr, "tabix: query %s: %v\n", r.chrom, err)
			return 1
		}
		for _, rec := range records {
			if targets != nil && !targets.Overlaps(r.chrom, rec.Beg, rec.End) {
				continue
			}
			stdout.Write(rec.Line)
			fmt.Fprintln(stdout)
		}
	}
	return 0
}

// emitHeader streams every line at the top of dataPath that begins with the
// configured comment character to stdout, stopping at the first data line.
func emitHeader(dataPath string, cfg tabix.Config, stdout io.Writer) error {
	f, err := os.Open(dataPath)
	if err != nil {
		return err
	}
	defer f.Close()
	br, err := bgzip.NewReader(f)
	if err != nil {
		return err
	}
	defer br.Close()
	scan := bufio.NewScanner(br)
	scan.Buffer(make([]byte, 0, 1<<16), 1<<24)
	for scan.Scan() {
		line := scan.Bytes()
		if len(line) == 0 || line[0] != byte(cfg.Meta) {
			break
		}
		stdout.Write(line)
		fmt.Fprintln(stdout)
	}
	return scan.Err()
}
