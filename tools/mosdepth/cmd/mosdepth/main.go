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

const usage = `mosdepth - pure-Go per-base/per-region BAM/CRAM depth.

Usage:
  mosdepth [options] <prefix> <in.bam|in.cram>

Outputs:
  <prefix>.mosdepth.global.dist.txt   cumulative depth distribution.
  <prefix>.mosdepth.summary.txt       per-chrom + total summary.
  <prefix>.per-base.bed.gz [.csi]     per-base depth (omitted with --by, --no-per-base, or --d4).
  <prefix>.per-base.d4                per-base depth in D4 binary format (when --d4 is set).
  <prefix>.regions.bed.gz [.csi]      per-region depth (when --by is set).
  <prefix>.quantized.bed.gz [.csi]    quantized depth segments (when --quantize set).
  <prefix>.thresholds.bed.gz [.csi]   threshold proportions (when --thresholds set).

Options:
  -t, --threads INT       BGZF/BAM decompression threads (output is identical for any value).
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
  -f, --fasta FILE        FASTA reference for decoding CRAM input (--reference alias; honours REF_CACHE).
  -a, --fragment-mode     count coverage across the whole fragment (proper pairs only); excludes -x.
  -q, --quantize SEGS     ':'-separated depth bins, e.g. 0:1:4:; writes <prefix>.quantized.bed.gz.
  -m, --use-median        report the per-region MEDIAN depth (with --by) instead of the mean.
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
  - -t/--threads spreads BGZF block decompression across N goroutines;
    the decoded stream and every output file are byte-identical for any
    thread count (it only affects throughput).
  - -m/--use-median changes ONLY the regions.bed.gz depth column to the
    histogram median (matching upstream's depthstat.CountStat); the
    summary, distribution, thresholds, and per-base outputs are unaffected.
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
	reference   string
	fragmentLen bool
	quantize    string
	useMedian   bool
	showHelp    bool
	showVersion bool
}

// reorderArgs permutes an argv so that all option tokens (and the values they
// consume) precede the positional arguments, mirroring the GNU getopt /
// docopt permutation that upstream mosdepth relies on. Go's flag package stops
// parsing options at the first non-flag argument, so without this pre-pass an
// interspersed command line like "mosdepth PREFIX in.bam --chrom MT" would
// treat "--chrom" and "MT" as a (rejected) third and fourth positional. After
// reordering, the existing parser sees the canonical flags-first form.
//
// The fs argument is introspected (via cliflag.Normalize for short-flag
// clusters and the per-flag value/bool classification below) so a value-taking
// flag's argument — even one that begins with '-', such as "-q -5" — is moved
// alongside its flag rather than mistaken for a positional. Everything after a
// "--" terminator is treated as positional verbatim. When the argv is already
// in flags-first order the output is identical to the input, so the function
// is a no-op for non-interspersed command lines.
func reorderArgs(fs *flag.FlagSet, args []string) ([]string, error) {
	// First expand any short-flag clusters ("-nx", "-Q20") into canonical
	// one-token-per-flag form and split out inline values. Normalize also
	// honours the "--" terminator and leaves long options untouched.
	norm, err := cliflag.Normalize(fs, args)
	if err != nil {
		return nil, err
	}
	options := make([]string, 0, len(norm))
	positionals := make([]string, 0, len(norm))
	for i := 0; i < len(norm); i++ {
		arg := norm[i]
		if arg == "--" {
			// End of options: the terminator itself and every remaining
			// token are positional. Keep the terminator with the options so
			// flag.Parse still stops option scanning at the right place.
			options = append(options, arg)
			positionals = append(positionals, norm[i+1:]...)
			break
		}
		if len(arg) >= 2 && arg[0] == '-' && arg != "-" {
			// An option token. Determine whether it carries an inline value
			// ("--chrom=MT" or, post-Normalize, "-c" "MT") or consumes the
			// following argument as its value.
			options = append(options, arg)
			name, hasInline := flagName(arg)
			if !hasInline && !isBoolFlag(fs, name) && i+1 < len(norm) {
				// Value-taking flag with no inline value: the next token is
				// its option-argument and must travel with it, even if it
				// begins with '-'.
				options = append(options, norm[i+1])
				i++
			}
			continue
		}
		// A bare "-" or any other token is a positional argument.
		positionals = append(positionals, arg)
	}
	return append(options, positionals...), nil
}

// flagName extracts the registered flag name from an option token and reports
// whether the token already carries an inline value. For "--chrom=MT" it
// returns ("chrom", true); for "--chrom" or "-c" it returns the bare name and
// false. Both single-dash short flags and double-dash long flags are handled.
func flagName(arg string) (name string, hasInline bool) {
	s := arg
	if len(s) >= 2 && s[0] == '-' && s[1] == '-' {
		s = s[2:]
	} else if len(s) >= 1 && s[0] == '-' {
		s = s[1:]
	}
	if eq := indexByte(s, '='); eq >= 0 {
		return s[:eq], true
	}
	return s, false
}

// isBoolFlag reports whether the flag registered under name on fs is a boolean
// switch (one that takes no value). It returns false for unknown or
// value-taking flags so reorderArgs treats their following argument as a
// consumed option-value rather than a positional.
func isBoolFlag(fs *flag.FlagSet, name string) bool {
	f := fs.Lookup(name)
	if f == nil {
		return false
	}
	if bf, ok := f.Value.(interface{ IsBoolFlag() bool }); ok {
		return bf.IsBoolFlag()
	}
	return false
}

// indexByte returns the index of the first occurrence of b in s, or -1.
func indexByte(s string, b byte) int {
	for i := 0; i < len(s); i++ {
		if s[i] == b {
			return i
		}
	}
	return -1
}

func parseFlags(args []string) (*runOptions, []string, error) {
	fs := flag.NewFlagSet("mosdepth", flag.ContinueOnError)
	fs.SetOutput(os.Stderr)
	fs.Usage = func() { fmt.Fprint(os.Stderr, usage) }
	opts := &runOptions{flag: int(mosdepth.DefaultExcludeFlag)}
	cliflag.IntVar(fs, &opts.threads, "t", "threads", 1, "BGZF decompression threads")
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
	// -f/--fasta names the FASTA reference used to decode reference-backed
	// CRAM input; it is ignored for BAM and SAM (which carry sequence
	// inline). --reference is a samtools-style long alias for the same value.
	// -a/--fragment-mode, -q/--quantize and -m/--use-median are all
	// implemented and wired into mosdepth.Options.
	cliflag.StringVar(fs, &opts.fasta, "f", "fasta", "", "FASTA reference for CRAM input")
	fs.StringVar(&opts.reference, "reference", "", "FASTA reference for CRAM input (samtools-style alias for -f/--fasta)")
	cliflag.BoolVar(fs, &opts.fragmentLen, "a", "fragment-mode", false, "count full-fragment coverage (proper pairs only)")
	cliflag.StringVar(fs, &opts.quantize, "q", "quantize", "", "quantized output segments, e.g. 0:1:4:")
	cliflag.BoolVar(fs, &opts.useMedian, "m", "use-median", false, "report per-region median depth instead of mean")
	cliflag.BoolVar(fs, &opts.showHelp, "h", "help", false, "help")
	cliflag.BoolVar(fs, &opts.showVersion, "v", "version", false, "version")
	// Permute options ahead of the two positionals so flags interspersed among
	// or after the positionals (the natural docopt order upstream accepts, e.g.
	// "mosdepth PREFIX in.bam --chrom MT") parse identically to the flags-first
	// form. reorderArgs already performs the short-flag cluster expansion that
	// cliflag.Parse would, so the FlagSet is parsed directly afterwards.
	reordered, err := reorderArgs(fs, args)
	if err != nil {
		return nil, nil, err
	}
	if err := fs.Parse(reordered); err != nil {
		return nil, nil, err
	}
	// Resolve the -r / -R read-group aliases: -R (upstream) wins when both
	// are set; otherwise whichever is non-empty.
	if opts.readGroupsR != "" {
		opts.readGroups = opts.readGroupsR
	}
	// Resolve the -f/--fasta and --reference aliases for the CRAM decode
	// reference: --reference (samtools-style) wins when both are set;
	// otherwise whichever is non-empty.
	if opts.reference != "" {
		opts.fasta = opts.reference
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
	// Upstream's hot loop lets --fragment-mode take precedence over
	// --fast-mode; the two are conceptually exclusive (fragment coverage vs.
	// per-read fast coverage). Reject the combination so users get a clear
	// error instead of a silently-fragment result.
	if opts.fragmentLen && opts.fastMode {
		fmt.Fprintln(os.Stderr, "mosdepth: -a/--fragment-mode and -x/--fast-mode are mutually exclusive")
		return 2
	}
	// A maximum fragment length below the minimum keeps no reads; upstream
	// rejects this with exit 2 before opening the input. Mirror that here so
	// the CLI never produces empty output for an impossible filter.
	if opts.minFragLen > 0 && opts.maxFragLen > 0 && opts.maxFragLen < opts.minFragLen {
		fmt.Fprintln(os.Stderr, mosdepth.ErrFragLenBounds)
		return 2
	}
	quant, err := mosdepth.ParseQuantize(opts.quantize)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
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
		FragmentMode: opts.fragmentLen,
		Quantize:     quant,
		Threads:      opts.threads,
		UseMedian:    opts.useMedian,
		Fasta:        opts.fasta,
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
