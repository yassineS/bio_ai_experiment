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
  -R, --regions FILE               Read regions from a BED-like file.
  -T, --targets FILE               Stream and post-filter to records inside FILE.
  -r, --reheader FILE              Replace the header with the content of FILE.
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
	cliflag.StringVar(fs, &o.targetsFile, "T", "targets", "", "filter records by targets file")
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
	case len(regions) == 0 && o.regionsFile == "" && o.targetsFile == "" && !o.onlyHdr:
		// Build mode.
		return runBuild(dataPath, cfg, o.force, o.noSaveIdx, stderr)
	default:
		return runQuery(dataPath, regions, o, cfg, stdout, stderr)
	}
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

	// --targets builds a post-filter applied to every emitted record. Unlike
	// --regions, the targets file does NOT contribute query regions: the file
	// is streamed (default: whole file) and only records that strictly overlap
	// a target interval are emitted. This mirrors upstream tabix.c, where the
	// targets file is loaded into a regidx and each record passing the index
	// query is checked with regidx_overlap(reg_idx, seq, curr_beg, curr_end-1).
	var tgt *targetSet
	if o.targetsFile != "" {
		t, err := loadTargets(o.targetsFile)
		if err != nil {
			fmt.Fprintf(stderr, "tabix: %v\n", err)
			return 1
		}
		tgt = t
	}

	// When no query regions were supplied (neither --regions nor positional),
	// scan the whole file. Upstream represents this as the region ".", which
	// expands to every sequence in the index.
	if len(rs) == 0 {
		for _, name := range idx.Names {
			rs = append(rs, region{chrom: name, beg: 0, end: 1 << 30})
		}
	}

	for _, r := range rs {
		records, err := idx.QueryBytes(dataPath, r.chrom, r.beg, r.end)
		if err != nil {
			fmt.Fprintf(stderr, "tabix: query %s: %v\n", r.chrom, err)
			return 1
		}
		for _, rec := range records {
			if tgt != nil && !tgt.overlapsRecord(rec, idx.Config) {
				continue
			}
			stdout.Write(rec)
			fmt.Fprintln(stdout)
		}
	}
	return 0
}

// targetInterval is one [beg, end) half-open 0-based target region.
type targetInterval struct {
	beg, end int
}

// targetSet holds the parsed --targets intervals grouped by sequence name.
type targetSet struct {
	byChrom map[string][]targetInterval
}

// loadTargets parses a targets file into a targetSet. The coordinate
// convention matches htslib's regidx: files whose name ends in .bed/.bed.gz/
// .bed.bgz are read as 0-based half-open BED (cols: chrom, start, end);
// every other file is read as 1-based inclusive "tab" format (chrom, beg,
// end), which is internally converted to 0-based half-open. A line with only
// a chromosome name covers the whole sequence.
func loadTargets(path string) (*targetSet, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()

	isBED := strings.HasSuffix(path, ".bed") ||
		strings.HasSuffix(path, ".bed.gz") ||
		strings.HasSuffix(path, ".bed.bgz")

	ts := &targetSet{byChrom: make(map[string][]targetInterval)}
	scan := bufio.NewScanner(f)
	scan.Buffer(make([]byte, 0, 1<<16), 1<<24)
	for scan.Scan() {
		line := strings.TrimSpace(scan.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Fields(line)
		chrom := fields[0]
		if len(fields) < 2 {
			ts.byChrom[chrom] = append(ts.byChrom[chrom], targetInterval{0, 1 << 30})
			continue
		}
		beg, err := strconv.Atoi(fields[1])
		if err != nil {
			return nil, fmt.Errorf("targets %s: bad begin %q", path, fields[1])
		}
		var end int
		if isBED {
			// BED: 0-based half-open. end defaults to beg+1 if absent.
			if len(fields) >= 3 {
				end, err = strconv.Atoi(fields[2])
				if err != nil {
					return nil, fmt.Errorf("targets %s: bad end %q", path, fields[2])
				}
			} else {
				end = beg + 1
			}
		} else {
			// Tab: 1-based inclusive -> 0-based half-open.
			beg--
			if len(fields) >= 3 {
				e, err := strconv.Atoi(fields[2])
				if err == nil {
					end = e // inclusive 1-based end == exclusive 0-based end
				} else {
					end = beg + 1
				}
			} else {
				end = beg + 1
			}
		}
		if end <= beg {
			end = beg + 1
		}
		ts.byChrom[chrom] = append(ts.byChrom[chrom], targetInterval{beg, end})
	}
	return ts, scan.Err()
}

// overlapsRecord reports whether rec strictly overlaps any target interval.
// The record's [beg, end) interval is derived exactly as the tabix index
// does (via cfg's columns and coordinate flags), so the comparison matches
// upstream's regidx_overlap(reg_idx, seq, curr_beg, curr_end-1) post-filter.
func (ts *targetSet) overlapsRecord(rec []byte, cfg tabix.Config) bool {
	chrom, beg, end, ok := recordInterval(rec, cfg)
	if !ok {
		return false
	}
	for _, iv := range ts.byChrom[chrom] {
		// Half-open overlap: [beg,end) intersects [iv.beg,iv.end).
		if beg < iv.end && iv.beg < end {
			return true
		}
	}
	return false
}

// recordInterval extracts the sequence name and 0-based half-open [beg, end)
// span of a single record line, replicating the tabix index's column parsing
// (htslib's get_intv / tbx_parse1).
func recordInterval(rec []byte, cfg tabix.Config) (chrom string, beg, end int, ok bool) {
	fields := strings.Split(string(rec), "\t")
	if int(cfg.ColSeq) < 1 || int(cfg.ColSeq) > len(fields) || int(cfg.ColBeg) < 1 || int(cfg.ColBeg) > len(fields) {
		return "", 0, 0, false
	}
	chrom = fields[cfg.ColSeq-1]
	bv, err := strconv.Atoi(fields[cfg.ColBeg-1])
	if err != nil {
		return "", 0, 0, false
	}
	if cfg.Format&tabix.FlagZeroBased == 0 {
		bv-- // 1-based inclusive -> 0-based half-open begin
	}
	beg = bv
	end = beg + 1
	switch {
	case cfg.ColEnd > 0:
		if int(cfg.ColEnd) <= len(fields) {
			if ev, err := strconv.Atoi(fields[cfg.ColEnd-1]); err == nil {
				end = ev
			}
		}
	case cfg.Format&0xFFFF == tabix.FormatVCF:
		if len(fields) >= 4 {
			end = beg + len(fields[3])
		}
	}
	if end <= beg {
		end = beg + 1
	}
	if beg < 0 {
		beg = 0
	}
	return chrom, beg, end, true
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
