// Command mosdepth is a pure-Go re-implementation of brentp's mosdepth
// per-base / per-region depth-of-coverage tool. See
// tools/mosdepth/README.md for usage and the per-flag parity matrix.
package main

import (
	"flag"
	"fmt"
	"os"
	"strconv"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/tools/mosdepth/pkg/mosdepth"
)

const version = "0.1.0"

const usage = `mosdepth - pure-Go per-base/per-region BAM depth.

Usage:
  mosdepth [options] <prefix> <in.bam>

Outputs:
  <prefix>.mosdepth.global.dist.txt   cumulative depth distribution.
  <prefix>.mosdepth.summary.txt       per-chrom + total summary.
  <prefix>.per-base.bed.gz [.csi]     per-base depth (omitted with --by, --no-per-base, or --d4).
  <prefix>.per-base.d4                per-base depth in D4 binary format (when --d4 is set).
  <prefix>.regions.bed.gz [.csi]      per-region depth (when --by is set).
  <prefix>.thresholds.bed.gz [.csi]   threshold proportions (when --thresholds set).

Options:
  -t, --threads INT       accepted; v1 is single-threaded.
  -b, --by FILE_OR_INT    BED of regions or fixed integer window size.
  -Q, --mapq INT          minimum MAPQ (default 0).
  -F, --flag INT          exclude reads with ANY of these flag bits (default 1796).
  -i, --include-flag INT  keep only reads with ALL of these flag bits.
  -x, --fast-mode         skip CIGAR walking (faster, slightly inaccurate near indels).
  -n, --no-per-base       suppress the per-base output file.
  -T, --thresholds LIST   comma list of integer thresholds, e.g. 1,5,10,30.
  -c, --chrom STRING      restrict to one chromosome.
      --d4                write per-base depth to <prefix>.per-base.d4 (D4 format) instead of BED.
  -d                      port-only short alias for --d4.
  -R, --read-groups LIST  comma list of allowed RG ids, or "OPS:X,Y" for the OPS aux tag.
  -r                      port-only lowercase alias for -R/--read-groups.
  -l, --min-frag-len INT  minimum absolute TLEN.
  -u, --max-frag-len INT  maximum absolute TLEN.
  -f, --fasta FILE        FASTA reference for CRAM input (accepted; CRAM not yet supported, ignored).
  -a, --fragment-mode     full-fragment coverage (upstream flag; not yet implemented — rejected).
  -q, --quantize SEGS     quantized output (upstream flag; not yet implemented — rejected).
  -m, --use-median        per-region median (upstream flag; not yet implemented — rejected).
  -h, --help              show this help.
  -v, --version           print version and exit.

Short-flag bundling (docopt parity): single-char short flags may be
clustered like upstream's docopt parser, e.g. "-nx" == "-n -x" and
"-Q20" == "-Q 20".

Deviations from upstream mosdepth (Nim):
  - D4 output (--d4) is byte-identical to the upstream mosdepth
    binary: a real D4 framefile with a 7-bit-packed primary table
    (SimpleRange{0,128}), validated against the real mosdepth_d4 binary.
    Upstream exposes only the long --d4; this port additionally accepts
    a -d short alias.
  - Upstream's short read-groups flag is -R; this port also accepts a
    lowercase -r alias.
  - Threads is accepted for compatibility; the v1 engine is
    single-threaded.
  - -a/--fragment-mode, -q/--quantize and -m/--use-median are parsed for
    CLI parity but not yet implemented; supplying them is rejected
    (exit 2) rather than silently ignored.
`

type runOptions struct {
	threads     int
	by          string
	mapq        int
	flag        int
	includeFlag int
	fastMode    bool
	noPerBase   bool
	noPerBase2  bool
	thresholds  string
	chrom       string
	d4          bool
	readGroups  string
	readGroupsR string
	minFragLen  int
	maxFragLen  int
	fasta       string
	fragmentLen bool
	quantize    string
	useMedian   bool
	showHelp    bool
	showVersion bool
}

func parseFlags(args []string) (*runOptions, []string, error) {
	fs := flag.NewFlagSet("mosdepth", flag.ContinueOnError)
	fs.SetOutput(os.Stderr)
	fs.Usage = func() { fmt.Fprint(os.Stderr, usage) }
	opts := &runOptions{flag: int(mosdepth.DefaultExcludeFlag)}
	cliflag.IntVar(fs, &opts.threads, "t", "threads", 1, "accepted; single-threaded")
	cliflag.StringVar(fs, &opts.by, "b", "by", "", "BED file or window size")
	cliflag.IntVar(fs, &opts.mapq, "Q", "mapq", 0, "min MAPQ")
	cliflag.IntVar(fs, &opts.flag, "F", "flag", int(mosdepth.DefaultExcludeFlag), "exclude flag bits")
	cliflag.IntVar(fs, &opts.includeFlag, "i", "include-flag", 0, "include flag bits")
	cliflag.BoolVar(fs, &opts.fastMode, "x", "fast-mode", false, "skip CIGAR walking")
	cliflag.BoolVar(fs, &opts.noPerBase, "n", "no-per-base", false, "suppress per-base output")
	// Also accept long form alone for users that don't know the short alias.
	fs.BoolVar(&opts.noPerBase2, "no-per-base-only", false, "suppress per-base output (deprecated alias)")
	cliflag.StringVar(fs, &opts.thresholds, "T", "thresholds", "", "thresholds list")
	cliflag.StringVar(fs, &opts.chrom, "c", "chrom", "", "restrict to one chromosome")
	// --d4 is the upstream-canonical name (upstream has no short form for
	// it); -d is a port-only convenience alias retained for backward
	// compatibility with earlier releases of this port. Both target the
	// same boolean.
	cliflag.BoolVar(fs, &opts.d4, "d", "d4", false, "D4 output")
	// -R/--read-groups is the upstream spelling (capital R, docopt). -r is
	// a port-only lowercase alias retained for backward compatibility;
	// both feed the same value (last-set wins, resolved below).
	cliflag.StringVar(fs, &opts.readGroupsR, "R", "read-groups", "", "comma list of RG IDs (or OPS:...)")
	fs.StringVar(&opts.readGroups, "r", "", "comma list of RG IDs (port-only lowercase alias for -R)")
	cliflag.IntVar(fs, &opts.minFragLen, "l", "min-frag-len", 0, "min |TLEN|")
	cliflag.IntVar(fs, &opts.maxFragLen, "u", "max-frag-len", 0, "max |TLEN|")
	// Upstream flags this port parses for CLI parity but does not yet
	// implement. -f/--fasta only matters for CRAM input (not supported
	// yet) so it is accepted and ignored; -a/--fragment-mode,
	// -q/--quantize and -m/--use-median change the numerical output, so
	// supplying them is rejected in run() rather than silently ignored.
	cliflag.StringVar(fs, &opts.fasta, "f", "fasta", "", "FASTA reference for CRAM (accepted; CRAM not yet supported)")
	cliflag.BoolVar(fs, &opts.fragmentLen, "a", "fragment-mode", false, "count full-fragment coverage (not yet implemented)")
	cliflag.StringVar(fs, &opts.quantize, "q", "quantize", "", "quantized output segments (not yet implemented)")
	cliflag.BoolVar(fs, &opts.useMedian, "m", "use-median", false, "use per-region median (not yet implemented)")
	cliflag.BoolVar(fs, &opts.showHelp, "h", "help", false, "help")
	cliflag.BoolVar(fs, &opts.showVersion, "v", "version", false, "version")
	if err := cliflag.Parse(fs, args); err != nil {
		return nil, nil, err
	}
	// Resolve the -r / -R read-group aliases: -R (upstream) wins when both
	// are set; otherwise whichever is non-empty.
	if opts.readGroupsR != "" {
		opts.readGroups = opts.readGroupsR
	}
	return opts, fs.Args(), nil
}

func run(args []string) int {
	opts, positional, err := parseFlags(args)
	if err != nil {
		return 2
	}
	if opts.showHelp {
		fmt.Print(usage)
		return 0
	}
	if opts.showVersion {
		fmt.Println(version)
		return 0
	}
	// Reject upstream flags that this port parses for CLI parity but does
	// not yet implement, when they would change the output. Silently
	// ignoring them would produce results that disagree with upstream.
	if opts.fragmentLen {
		fmt.Fprintln(os.Stderr, "mosdepth: -a/--fragment-mode is not yet implemented in this port")
		return 2
	}
	if opts.quantize != "" {
		fmt.Fprintln(os.Stderr, "mosdepth: -q/--quantize is not yet implemented in this port")
		return 2
	}
	if opts.useMedian {
		fmt.Fprintln(os.Stderr, "mosdepth: -m/--use-median is not yet implemented in this port")
		return 2
	}
	if len(positional) != 2 {
		fmt.Fprint(os.Stderr, usage)
		fmt.Fprintln(os.Stderr, "mosdepth: expected exactly two positional arguments: <prefix> <in.bam>")
		return 2
	}
	prefix := positional[0]
	bam := positional[1]

	winSize, bedPath, err := mosdepth.ParseByArg(opts.by)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 1
	}
	th, err := parseThresholdsFlag(opts.thresholds)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 1
	}
	mo := mosdepth.Options{
		Prefix:      prefix,
		ByBED:       bedPath,
		ByWindow:    winSize,
		MinMAPQ:     uint8Clamp(opts.mapq),
		ExcludeFlag: uint16(opts.flag),
		IncludeFlag: uint16(opts.includeFlag),
		FastMode:    opts.fastMode,
		NoPerBase:   opts.noPerBase || opts.noPerBase2,
		Thresholds:  th,
		Chrom:       opts.chrom,
		D4Output:    opts.d4,
		ReadGroups:  mosdepth.ParseReadGroups(opts.readGroups),
		MinFragLen:  opts.minFragLen,
		MaxFragLen:  opts.maxFragLen,
		Threads:     opts.threads,
	}
	if err := mosdepth.OpenAndRun(bam, mo); err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 1
	}
	return 0
}

func parseThresholdsFlag(spec string) ([]int, error) {
	if spec == "" {
		return nil, nil
	}
	out := []int{}
	for _, part := range splitCSV(spec) {
		v, err := strconv.Atoi(part)
		if err != nil {
			return nil, fmt.Errorf("mosdepth: bad --thresholds value %q: %w", part, err)
		}
		if v < 0 {
			return nil, fmt.Errorf("mosdepth: threshold must be >= 0, got %d", v)
		}
		out = append(out, v)
	}
	return out, nil
}

func splitCSV(s string) []string {
	out := []string{}
	cur := ""
	for i := 0; i < len(s); i++ {
		if s[i] == ',' {
			if cur != "" {
				out = append(out, cur)
				cur = ""
			}
			continue
		}
		cur += string(s[i])
	}
	if cur != "" {
		out = append(out, cur)
	}
	return out
}

func uint8Clamp(v int) uint8 {
	if v < 0 {
		return 0
	}
	if v > 255 {
		return 255
	}
	return uint8(v)
}

func main() {
	os.Exit(run(os.Args[1:]))
}
