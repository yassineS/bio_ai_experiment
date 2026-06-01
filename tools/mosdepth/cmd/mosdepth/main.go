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
  <prefix>.per-base.bed.gz [.tbi]     per-base depth (omitted with --by or --no-per-base).
  <prefix>.regions.bed.gz [.tbi]      per-region depth (when --by is set).
  <prefix>.thresholds.bed.gz [.tbi]   threshold proportions (when --thresholds set).

Options:
  -t, --threads INT       accepted; v1 is single-threaded.
  -b, --by FILE_OR_INT    BED of regions or fixed integer window size.
  -Q, --mapq INT          minimum MAPQ (default 0).
  -F, --flag INT          exclude reads with ANY of these flag bits (default 1796).
  -i, --include-flag INT  keep only reads with ALL of these flag bits.
  -x, --fast-mode         skip CIGAR walking (faster, slightly inaccurate near indels).
      --no-per-base       suppress the per-base output file.
  -n, --no-per-base       alias for --no-per-base.
  -T, --thresholds LIST   comma list of integer thresholds, e.g. 1,5,10,30.
  -c, --chrom STRING      restrict to one chromosome.
  -d, --d4                request D4 output (not yet implemented).
  -r, --read-groups LIST  comma list of allowed RG ids, or "OPS:X,Y" for the OPS aux tag.
  -l, --min-frag-len INT  minimum absolute TLEN.
  -u, --max-frag-len INT  maximum absolute TLEN.
  -a, --fragment-mode     score each pair as a single fragment span.
  -m, --use-median        report region MEDIAN depth (with --by) instead of mean.
  -q, --quantize CUTOFFS  colon-separated cutoffs (e.g. 0:1:1000) for the
                          quantized.bed.gz output. Bin labels default to
                          "cutoffs[i]:cutoffs[i+1]" (e.g. "0:1", "1:1000");
                          override per bin via the MOSDEPTH_Q{i}
                          environment variables. The implicit
                          below-first-cutoff and open-ended top bins are
                          omitted from output to match upstream.
  -h, --help              show this help.
  -v, --version           print version and exit.

Deviations from upstream mosdepth (Nim):
  - Indexes are tabix .tbi (not .csi). Consumers that read TBI work
    transparently; the underlying chunk/bin layout is the same.
  - D4 output (-d/--d4) is rejected with a clear error in v1.
  - Threads is accepted for compatibility; the v1 engine is
    single-threaded.
`

type runOptions struct {
	threads      int
	by           string
	mapq         int
	flag         int
	includeFlag  int
	fastMode     bool
	noPerBase    bool
	noPerBase2   bool
	thresholds   string
	chrom        string
	d4           bool
	readGroups   string
	minFragLen   int
	maxFragLen   int
	fragmentMode bool
	useMedian    bool
	quantize     string
	showHelp     bool
	showVersion  bool
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
	cliflag.BoolVar(fs, &opts.d4, "d", "d4", false, "D4 output")
	cliflag.StringVar(fs, &opts.readGroups, "r", "read-groups", "", "comma list of RG IDs (or OPS:...)")
	cliflag.IntVar(fs, &opts.minFragLen, "l", "min-frag-len", 0, "min |TLEN|")
	cliflag.IntVar(fs, &opts.maxFragLen, "u", "max-frag-len", 0, "max |TLEN|")
	cliflag.BoolVar(fs, &opts.fragmentMode, "a", "fragment-mode", false, "score each pair as a fragment span")
	cliflag.BoolVar(fs, &opts.useMedian, "m", "use-median", false, "report region median depth instead of mean")
	cliflag.StringVar(fs, &opts.quantize, "q", "quantize", "", "colon-separated cutoffs (e.g. 0:1:1000) for quantized.bed.gz")
	cliflag.BoolVar(fs, &opts.showHelp, "h", "help", false, "help")
	cliflag.BoolVar(fs, &opts.showVersion, "v", "version", false, "version")
	if err := fs.Parse(args); err != nil {
		return nil, nil, err
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
	qz, err := mosdepth.ParseQuantize(opts.quantize)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 1
	}
	mo := mosdepth.Options{
		Prefix:       prefix,
		ByBED:        bedPath,
		ByWindow:     winSize,
		MinMAPQ:      uint8Clamp(opts.mapq),
		ExcludeFlag:  uint16(opts.flag),
		IncludeFlag:  uint16(opts.includeFlag),
		FastMode:     opts.fastMode,
		NoPerBase:    opts.noPerBase || opts.noPerBase2,
		Thresholds:   th,
		Chrom:        opts.chrom,
		D4Output:     opts.d4,
		ReadGroups:   mosdepth.ParseReadGroups(opts.readGroups),
		MinFragLen:   opts.minFragLen,
		MaxFragLen:   opts.maxFragLen,
		FragmentMode: opts.fragmentMode,
		UseMedian:    opts.useMedian,
		Quantize:     qz,
		Threads:      opts.threads,
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
