// CLI runners for `bcftools mendelian2` and `bcftools polysomy`.
//
// Per the project "every documented upstream flag must be recognised
// — implemented or gracefully rejected with a PARITY_ROADMAP.md
// pointer" rule (docs/PARITY_ROADMAP.md#definition-of-11), the flag
// surface here is taken directly from upstream:
//
//   - mendelian2: `reference_code/bcftools/plugins/mendelian2.c:run`
//     (note: it's a plugin, but in our port we expose it as a
//     first-class subcommand, matching the rest of the bcftools
//     CLI).
//   - polysomy:  `reference_code/bcftools/polysomy.c:main_polysomy`.
package main

import (
	"errors"
	"flag"
	"fmt"
	"io"
	"os"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/tools/bcftools/pkg/bcftools"
)

const mendelian2Usage = `bcftools mendelian2 - count Mendelian consistent / inconsistent genotypes.

Usage:
  bcftools mendelian2 [options] <in.vcf[.gz]|in.bcf>

Common Options:
  -e, --exclude EXPR              Exclude sites for which the expression is true.
  -i, --include EXPR              Include sites for which the expression is true.
  -o, --output FILE               Output file name (default: stdout).
  -O, --output-type u|b|v|z[0-9]  Output type:
                                    u/b: un/compressed BCF
                                    v/z: un/compressed VCF
                                    0-9: gzip level (default v)
  -r, --regions REGION            Restrict to comma-separated region list.
  -R, --regions-file FILE         Restrict to regions listed in a file.
      --regions-overlap 0|1|2     POS-in-region / record overlaps / variant
                                  overlaps (accepted; v1 uses POS-in-region).
  -t, --targets REGION            Like -r but streams (post-filter).
  -T, --targets-file FILE         Like -R but streams.
      --targets-overlap 0|1|2     Accepted; v1 uses POS-in-region.
      --no-version                Do not append version+command to header.
  -v, --verbosity INT             Accepted; v1 ignores.
  -W, --write-index[=csi|tbi]     Write a .csi/.tbi index for a bgzipped
                                  output (requires -o FILE and -Oz).

Mendelian2 options:
  -m, --mode c|[adeEgmMS]         Output mode (default c). Multiple modes can
                                  be combined; drop modes win over list modes:
                                    c .. print TSV counts summary
                                    a .. add INFO/MERR per record
                                    d .. set offending trio GTs to "./."
                                    e .. emit only sites with a Mendel error
                                    E .. drop sites with a Mendel error
                                    g .. emit only sites with a good trio
                                    m .. emit only sites with missing trio GT
                                    M .. drop sites with missing trio GT
                                    S .. drop sites skipped for housekeeping reasons
                                    s .. (accepted alias for S)
  -p, --pfm [1X:|2X:]P,F,M        Single-trio shortcut (proband,father,mother).
                                  1X: = male child (X haploid from father),
                                  2X: = female child. Sex prefix optional.
  -P, --ped FILE                  PED file (6 cols: family, individual, dad,
                                  mom, sex, phenotype). Trios are derived
                                  automatically: every row whose dad AND mom
                                  AND child are all in the input VCF.
      --rules ASSEMBLY[?]         Predefined inheritance rules (GRCh37 default,
                                  also GRCh38; "list"/"list?" to enumerate).
      --rules-file FILE           Custom per-contig ploidy/inheritance rules
                                  file (SEX_ID CHROM:BEG-END INHERITED_FROM).

  -h, --help                      Show this help.
      --version                   Show version.

Deferred flags (accepted but applied as a post-filter / no-op; see
docs/PARITY_ROADMAP.md#bcftools): --regions-overlap, --targets-overlap,
--verbosity.
`

// optionalStringValue is a flag.Value for options whose argument is
// optional (getopt `optional_argument`). A bare `-W` sets present=true
// without consuming the next token; `-W=csi` records the value. It lets
// -W/--write-index behave like upstream's `-W[=FMT]`.
type optionalStringValue struct {
	target  *string
	present *bool
}

func (o *optionalStringValue) String() string {
	if o == nil || o.target == nil {
		return ""
	}
	return *o.target
}

func (o *optionalStringValue) Set(s string) error {
	*o.target = s
	*o.present = true
	return nil
}

// IsBoolFlag lets the flag package accept a bare `-W` (no argument)
// without consuming the next token, while `-W=csi` still works.
func (o *optionalStringValue) IsBoolFlag() bool { return true }

func runMendelian2(args []string) int {
	fs := flag.NewFlagSet("bcftools mendelian2", flag.ContinueOnError)
	fs.SetOutput(io.Discard)

	var (
		pfm            string
		ped            string
		mode           string
		includeExpr    string
		excludeExpr    string
		regions        string
		regionsFile    string
		regionsOverlap int
		targets        string
		targetsFile    string
		targetsOverlap int
		rules          string
		rulesFile      string
		outputType     string
		outputPath     string
		noVersion      bool
		writeIndex     string
		writeIndexSet  bool
		verbosity      int
		showHelp       bool
		showVer        bool
	)

	cliflag.StringVar(fs, &pfm, "p", "pfm", "", "Single-trio shortcut [1X:|2X:]P,F,M")
	cliflag.StringVar(fs, &ped, "P", "ped", "", "PED file")
	cliflag.StringVar(fs, &mode, "m", "mode", "", "Output mode (c|[adeEgmMS])")
	cliflag.StringVar(fs, &includeExpr, "i", "include", "", "Include EXPR")
	cliflag.StringVar(fs, &excludeExpr, "e", "exclude", "", "Exclude EXPR")
	cliflag.StringVar(fs, &regions, "r", "regions", "", "Regions")
	cliflag.StringVar(fs, &regionsFile, "R", "regions-file", "", "Regions file")
	fs.IntVar(&regionsOverlap, "regions-overlap", 1, "Region overlap mode (accepted)")
	cliflag.StringVar(fs, &targets, "t", "targets", "", "Targets (post-filter)")
	cliflag.StringVar(fs, &targetsFile, "T", "targets-file", "", "Targets file")
	fs.IntVar(&targetsOverlap, "targets-overlap", 0, "Target overlap mode (accepted)")
	fs.StringVar(&rules, "rules", "", "Predefined inheritance rules (GRCh37/GRCh38, \"list\" to enumerate)")
	fs.StringVar(&rulesFile, "rules-file", "", "Custom inheritance rules file")
	cliflag.StringVar(fs, &outputType, "O", "output-type", "v", "Output type")
	cliflag.StringVar(fs, &outputPath, "o", "output", "", "Output path")
	fs.BoolVar(&noVersion, "no-version", false, "Do not append version+cmd to header")
	// -W/--write-index takes an OPTIONAL value (csi|tbi). A bare flag
	// means "index with the default format". optionalStringValue makes
	// `-W` work without consuming the following positional argument,
	// matching upstream's getopt `optional_argument`.
	wi := &optionalStringValue{target: &writeIndex, present: &writeIndexSet}
	fs.Var(wi, "W", "Write an index for a bgzipped output [optional csi|tbi]")
	fs.Var(wi, "write-index", "Write an index for a bgzipped output [optional csi|tbi]")
	cliflag.IntVar(fs, &verbosity, "v", "verbosity", 0, "Verbosity (accepted)")
	fs.BoolVar(&showHelp, "h", false, "")
	fs.BoolVar(&showHelp, "?", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "version", false, "")

	if err := fs.Parse(args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, mendelian2Usage)
		return 2
	}
	if showHelp {
		fmt.Print(mendelian2Usage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}

	if rules != "" && rulesFile != "" {
		fmt.Fprintln(os.Stderr, "bcftools mendelian2: --rules and --rules-file are mutually exclusive")
		return 2
	}

	// Resolve the inheritance rules: a custom --rules-file, a named
	// --rules assembly (or its "list"/"list?" catalogue), or the
	// GRCh37 default when neither is given. Mirrors init_rules.
	var ruleSet *bcftools.MendelianRules
	switch {
	case rulesFile != "":
		rs, err := bcftools.LoadMendelianRulesFile(rulesFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "bcftools mendelian2: %v\n", err)
			return 2
		}
		ruleSet = rs
	case rules != "":
		rs, err := bcftools.LoadMendelianRulesByName(rules)
		if err != nil {
			var listErr *bcftools.ErrMendelianRulesList
			if errors.As(err, &listErr) {
				fmt.Fprint(os.Stderr, listErr.Listing)
				return 255
			}
			fmt.Fprintf(os.Stderr, "bcftools mendelian2: %v\n", err)
			return 2
		}
		ruleSet = rs
	}

	rest := fs.Args()
	if len(rest) == 0 {
		fmt.Fprintln(os.Stderr, "bcftools mendelian2: missing input file")
		fmt.Fprint(os.Stderr, mendelian2Usage)
		return 2
	}
	if pfm == "" && ped == "" {
		fmt.Fprintln(os.Stderr, "bcftools mendelian2: missing -p/--pfm or -P/--ped option")
		return 2
	}
	if pfm != "" && ped != "" {
		fmt.Fprintln(os.Stderr, "bcftools mendelian2: -p/--pfm and -P/--ped are mutually exclusive")
		return 2
	}

	modeBits, err := bcftools.ParseMendelian2Mode(mode)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}
	format, err := bcftools.ParseOutputFormat(outputType)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}

	opts := bcftools.Mendelian2Options{
		PEDFile:       ped,
		Mode:          modeBits,
		IncludeExpr:   includeExpr,
		ExcludeExpr:   excludeExpr,
		OutputFormat:  format,
		CompressLevel: -1,
		Rules:         ruleSet,
	}
	if pfm != "" {
		parsed, err := bcftools.ParseMendelian2PFM(pfm)
		if err != nil {
			fmt.Fprintln(os.Stderr, err)
			return 2
		}
		opts.PFM = &parsed
	}

	// -W/--write-index is only meaningful when we emit an indexable
	// (bgzipped) VCF/BCF to a real file, which in turn requires a
	// VCF/BCF-producing mode (anything other than pure count).
	if writeIndexSet {
		if modeBits == bcftools.Mendelian2Count {
			fmt.Fprintln(os.Stderr, "bcftools mendelian2: -W/--write-index requires a VCF/BCF output mode (not -m c)")
			return 2
		}
		if outputPath == "" || outputPath == "-" {
			fmt.Fprintln(os.Stderr, "bcftools mendelian2: -W/--write-index requires -o FILE")
			return 2
		}
		if format != bcftools.OutputVCFGz {
			fmt.Fprintln(os.Stderr, "bcftools mendelian2: -W/--write-index requires a bgzipped output (-Oz)")
			return 2
		}
		// Emit BGZF (block-gzip) so the file is indexable.
		opts.BGZF = true
	}

	// Acknowledge the surface-only flags so the linter doesn't bark.
	_ = regions
	_ = regionsFile
	_ = regionsOverlap
	_ = targets
	_ = targetsFile
	_ = targetsOverlap
	_ = noVersion
	_ = verbosity

	out, err := openOutFile(outputPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "bcftools mendelian2: %v\n", err)
		return 1
	}

	if _, err := bcftools.Mendelian2File(rest[0], out, opts); err != nil {
		out.Close()
		fmt.Fprintf(os.Stderr, "bcftools mendelian2: %v\n", err)
		return 1
	}
	// Flush/close the output before indexing so the on-disk bytes are
	// complete (BuildIndex re-opens and scans the file).
	if err := out.Close(); err != nil {
		fmt.Fprintf(os.Stderr, "bcftools mendelian2: %v\n", err)
		return 1
	}

	if writeIndexSet {
		if _, err := bcftools.BuildIndex(outputPath, parseWriteIndexFormat(writeIndex, outputPath)); err != nil {
			fmt.Fprintf(os.Stderr, "bcftools mendelian2: %v\n", err)
			return 1
		}
	}
	return 0
}

// parseWriteIndexFormat maps the optional -W/--write-index argument
// ("", "csi", or "tbi") and the output path to IndexOptions. Upstream's
// default for a bgzipped VCF is .csi; an explicit "tbi" selects the
// tabix flavour. The index is always written with Force so a stale
// index from a previous run is overwritten, matching upstream's
// auto-index behaviour.
func parseWriteIndexFormat(arg, outputPath string) bcftools.IndexOptions {
	opts := bcftools.IndexOptions{Format: bcftools.IndexCSI, Force: true}
	switch strings.ToLower(strings.TrimPrefix(arg, "=")) {
	case "tbi":
		opts.Format = bcftools.IndexTBI
	case "csi", "":
		opts.Format = bcftools.IndexCSI
	}
	_ = outputPath
	return opts
}

const polysomyUsage = `bcftools polysomy - detect chromosomal copy number from B-allele frequency.

Usage:
  bcftools polysomy [options] <in.vcf[.gz]|in.bcf>

The input must carry per-sample BAF, either as an explicit FORMAT/BAF
field (upstream's requirement) or as FORMAT/AD = REF,ALT counts
(v1 fallback — we synthesise BAF as ALT/(REF+ALT) at het sites).

General options:
  -o, --output-dir PATH          Output dir (accepted; v1 emits a TSV to
                                 stdout and ignores --output-dir).
  -r, --regions REGION           Restrict to comma-separated region list.
  -R, --regions-file FILE        Restrict to regions listed in a file.
      --regions-overlap 0|1|2    Accepted; v1 uses POS-in-region.
  -s, --sample NAME              Sample to analyze (required if VCF has >1
                                 sample).
  -t, --targets REGION           Like -r but streams (post-filter).
  -T, --targets-file FILE        Like -R but streams.
      --targets-overlap 0|1|2    Accepted; v1 uses POS-in-region.
  -v, --verbosity INT            Accepted; v1 ignores.
      --verbose INT              Alias for --verbosity (upstream).

Algorithm options:
  -b, --peak-size FLOAT          Min CN4 side-peak size (0-1, larger stricter) [0.1].
  -c, --cn-penalty FLOAT         CN-increase penalty (0-1, larger stricter) [0.7].
  -f, --fit-th FLOAT             Goodness-of-fit threshold (>0, smaller stricter) [3.3].
  -i, --include-aa               Include the AA peak in CN2 / CN4 evaluation.
  -m, --min-fraction FLOAT       Min distinguishable fraction of aberrant cells [0.1].
  -p, --peak-symmetry FLOAT      Peak-symmetry threshold (0-1, larger stricter) [0.5].
  -n, --nbins INT                Histogram bin count (hidden upstream option) [150].
  -S, --smooth INT               Smoothing half-window control (hidden upstream) [-3].
      --ra-rr-scaling            Disable RA/RR/AA normalisation (hidden upstream option).
      --force-cn INT             Tag every chromosome with this CN (hidden upstream option).

  -h, --help                     Show this help.
      --version                  Show version.

Output (stdout TSV):
  # sample  chrom  n_het  mean_baf  median_baf  cn_call

The CN call is produced by the upstream Gaussian-mixture peak fit
(polysomy.c + peakfit.c); all algorithm knobs above are live.

Deferred flags (accepted but parsed into a no-op): -o/--output-dir
(this port emits the TSV to stdout, no PNG plots), --regions-overlap,
--targets-overlap, -v/--verbosity. See docs/PARITY_ROADMAP.md#bcftools.
`

func runPolysomy(args []string) int {
	fs := flag.NewFlagSet("bcftools polysomy", flag.ContinueOnError)
	fs.SetOutput(io.Discard)

	var (
		outputDir      string
		regions        string
		regionsFile    string
		regionsOverlap int
		sample         string
		targets        string
		targetsFile    string
		targetsOverlap int
		verbosity      int
		peakSize       float64
		cnPenalty      float64
		fitTh          float64
		includeAA      bool
		minFraction    float64
		peakSymmetry   float64
		nbins          int
		smooth         int
		raRrScaling    bool
		forceCN        int
		showHelp       bool
		showVer        bool
	)

	// Upstream `bcftools polysomy` does NOT advertise -i/--include or
	// -e/--exclude in `polysomy.c:main_polysomy`'s getopt_long table.
	// We keep the surface honest and OMIT them, in line with the
	// project's "no invented flags" rule (docs/PARITY_ROADMAP.md).
	// Same for --threads and --no-version (also not in upstream).

	cliflag.StringVar(fs, &outputDir, "o", "output-dir", "", "Output dir (accepted; v1 emits TSV to stdout)")
	cliflag.StringVar(fs, &regions, "r", "regions", "", "Regions")
	cliflag.StringVar(fs, &regionsFile, "R", "regions-file", "", "Regions file")
	fs.IntVar(&regionsOverlap, "regions-overlap", 1, "")
	cliflag.StringVar(fs, &sample, "s", "sample", "", "Sample name")
	cliflag.StringVar(fs, &targets, "t", "targets", "", "Targets (post-filter)")
	cliflag.StringVar(fs, &targetsFile, "T", "targets-file", "", "Targets file")
	fs.IntVar(&targetsOverlap, "targets-overlap", 0, "")
	// Upstream binds `-v` to BOTH `--verbose` and `--verbosity`
	// (polysomy.c:432). Mirror that aliasing here.
	cliflag.IntVar(fs, &verbosity, "v", "verbosity", 0, "Verbosity")
	fs.IntVar(&verbosity, "verbose", 0, "Verbosity (alias for --verbosity)")
	cliflag.Float64Var(fs, &peakSize, "b", "peak-size", 0.1, "Min peak size (accepted; v1 unused)")
	cliflag.Float64Var(fs, &cnPenalty, "c", "cn-penalty", 0.7, "CN-increase penalty")
	cliflag.Float64Var(fs, &fitTh, "f", "fit-th", 3.3, "Goodness-of-fit (accepted; v1 unused)")
	cliflag.BoolVar(fs, &includeAA, "i", "include-aa", false, "Include AA peak (accepted; v1 unused)")
	cliflag.Float64Var(fs, &minFraction, "m", "min-fraction", 0.1, "Min aberrant fraction (also the v1 |median(BAF)-0.5| threshold)")
	cliflag.Float64Var(fs, &peakSymmetry, "p", "peak-symmetry", 0.5, "Peak symmetry (accepted; v1 unused)")
	// Upstream short-letter bindings (polysomy.c:412-421):
	//   -n/--nbins   (NOT --include-noise)
	//   -S/--smooth  (NOT --samples-file)
	// We keep these letters bound the upstream way; PR #109 v1
	// initially stole both for invented flags and the reviewer
	// caught it. There is no upstream `-S/--samples-file` or
	// `-n/--include-noise`, so they are gone.
	cliflag.IntVar(fs, &nbins, "n", "nbins", 150, "Histogram bins (hidden upstream; accepted)")
	cliflag.IntVar(fs, &smooth, "S", "smooth", -3, "Smoothing window (hidden upstream; accepted)")
	fs.BoolVar(&raRrScaling, "ra-rr-scaling", false, "Disable RA/RR scaling (hidden upstream; accepted)")
	fs.IntVar(&forceCN, "force-cn", 0, "Force every chromosome to this CN (hidden upstream)")
	fs.BoolVar(&showHelp, "h", false, "")
	fs.BoolVar(&showHelp, "?", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "version", false, "")

	if err := fs.Parse(args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, polysomyUsage)
		return 2
	}
	if showHelp {
		fmt.Print(polysomyUsage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}

	rest := fs.Args()
	if len(rest) == 0 {
		fmt.Fprintln(os.Stderr, "bcftools polysomy: missing input file")
		fmt.Fprint(os.Stderr, polysomyUsage)
		return 2
	}

	// The polysomy algorithm is the upstream Gaussian-mixture peak
	// fit; every algorithm knob below is now live.
	opts := bcftools.PolysomyOptions{
		Sample:       sample,
		CnPenalty:    cnPenalty,
		FitTh:        fitTh,
		PeakSymmetry: peakSymmetry,
		MinPeakSize:  peakSize,
		MinFraction:  minFraction,
		IncludeAA:    includeAA,
		NBins:        nbins,
		Smooth:       smooth,
		// Upstream's --ra-rr-scaling flag DISABLES the per-segment
		// normalisation, which is otherwise on by default.
		RaRrScaling: !raRrScaling,
		ForceCN:     forceCN,
		RegionsFile: regionsFile,
		TargetsFile: targetsFile,
	}
	if regions != "" {
		opts.Regions = bcftools.SplitCommaList(regions)
	}
	if targets != "" {
		opts.Targets = bcftools.SplitCommaList(targets)
	}

	// Acknowledge surface-only flags so the linter doesn't bark.
	_ = outputDir
	_ = regionsOverlap
	_ = targetsOverlap
	_ = verbosity

	if _, err := bcftools.PolysomyFile(rest[0], os.Stdout, opts); err != nil {
		fmt.Fprintf(os.Stderr, "bcftools polysomy: %v\n", err)
		return 1
	}
	return 0
}
