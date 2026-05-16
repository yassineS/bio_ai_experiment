// CLI runner for `bcftools mpileup`. Follows the project parity rule
// (docs/PARITY_ROADMAP.md "Definition of 1:1"): every documented
// upstream flag must parse cleanly. Flags whose underlying behaviour
// is deferred either no-op (when their default is no-op) or
// hard-reject with a roadmap pointer when explicitly set.
//
// Upstream reference: reference_code/bcftools/mpileup.c — see the
// `lopts[]` and `getopt_long` block around line 1411.
package main

import (
	"compress/gzip"
	"flag"
	"fmt"
	"io"
	"os"
	"strconv"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/tools/bcftools/pkg/bcftools"
)

const mpileupUsage = `bcftools mpileup - per-position genotype likelihoods from BAM.

Usage:
  bcftools mpileup [options] -f ref.fa <in.bam> [<in2.bam> ...]

V1 SIMPLIFICATION: SNP-only, uniform-error genotype-likelihood model.
Indels, BAQ recalibration, MAQ-style likelihoods, and the full upstream
indel realigner are deferred. See docs/PARITY_ROADMAP.md#bcftools for
the option-tail gap list.

Input:
  -b, --bam-list FILE            File of BAM paths (one per line).
  -f, --fasta-ref FILE           Reference FASTA (required).
  -G, --read-groups FILE         Read-group file (accepted; v1 ignores).
      --ignore-RG                Ignore @RG tags (long-only upstream).

Output:
  -o, --output FILE              Output path (default stdout).
  -O, --output-type v|z|u|b      VCF (default), VCF.gz, BCF (u/b
                                 not yet implemented).
      --no-version               Skip version line in the header.
      --threads INT              Accepted; v1 is single-threaded.

Read filtering:
  -A, --count-orphans            Keep anomalous read pairs.
  -d, --max-depth INT            Max per-file depth (default 250).
  -q, --min-MQ INT               Skip reads with mapping quality < INT.
  -Q, --min-BQ INT               Skip bases with quality < INT (default 13).
      --max-bq INT               Cap base quality at INT (default 60).
      --delta-BQ INT             Accepted; v1 ignores.
  -x, --ignore-overlaps          Disable read-pair overlap detection.
  -C, --adjust-MQ INT            Accepted; v1 ignores.
  -r, --regions LIST             Region(s) (chr:beg-end[,...]).
  -R, --regions-file FILE        BED-like regions file.
  -t, --targets LIST             Like -r but always a post-filter.
  -T, --targets-file FILE        BED-like targets file.
  -s, --samples LIST             Sample list (BAM SM-tag filter).
  -S, --samples-file FILE        File of sample names.
      --skip-any-unset MASK      Accepted; v1 ignores.
      --skip-all-unset MASK      Accepted; v1 ignores.
      --skip-any-set MASK        Accepted; v1 ignores.
      --skip-all-set MASK        Accepted; v1 ignores.
      --ls MASK                  Accepted; v1 ignores.

Likelihood model (accepted; v1 uses uniform-error model):
  -B, --no-BAQ                   Disable BAQ (v1 default; no-op).
  -D, --full-BAQ                 Accepted; v1 ignores.
  -E, --redo-BAQ                 NOT IMPLEMENTED (hard-reject).
  -P, --platforms LIST           Accepted; v1 ignores.
  -p, --per-sample-mF            Accepted; v1 ignores.
  -6, --illumina1.3+             Accepted; v1 ignores.
  -X, --config STR               Predefined config (e.g. "1.12", "2.1").
                                 Accepted; v1 ignores.
      --seed INT                 Random seed for subsampling.
                                 Accepted; v1 ignores.
      --ambig-reads / --ar STR   Accepted; v1 ignores.

Indel model (accepted; v1 emits no indel records):
  -I, --skip-indels              Always on in v1.
  -L, --max-idepth INT           Accepted; v1 ignores.
  -m, --min-ireads INT           Accepted; v1 ignores.
  -h, --tandem-qual INT          Accepted; v1 ignores.
  -e, --ext-prob INT             Accepted; v1 ignores.
  -F, --gap-frac FLOAT           Accepted; v1 ignores.
  -M, --max-read-len INT         Accepted; v1 ignores.
      --indel-bias FLOAT         Accepted; v1 ignores.
      --indel-size INT           Accepted; v1 ignores.
      --indels-cns               Accepted; v1 ignores.
      --indels-2.0               Accepted; v1 ignores.
      --no-indels-cns            Accepted; v1 ignores.
      --open-prob INT            Accepted; v1 ignores.
      --del-bias FLOAT           Accepted; v1 ignores.
      --poly-mqual               Accepted; v1 ignores.
      --no-poly-mqual            Accepted; v1 ignores.
      --score-vs-ref FLOAT       Accepted; v1 ignores.
      --seqq-offset INT          Accepted; v1 ignores.

Output annotation:
  -a, --annotate LIST            FORMAT/INFO tags to add (accepted;
                                 v1 emits default DP, I16, PL).
  -g, --gvcf BLOCK               gVCF blocking (accepted; v1 emits
                                 one record per variant site).
      --no-reference             Accepted; v1 ignores.
  -W, --write-index[=FMT]        Accepted; v1 never auto-indexes.

  -v, --verbosity INT            Accepted; v1 ignores.
  -?, --help                     Show this help.
      --version                  Show version.

Note on -h:
  Upstream bcftools mpileup uses -h for "tandem-qual" rather than "help" —
  we follow that convention. Use -? or --help for help.
`

// mpileupFlags holds every flag value, including the accepted-but-
// deferred ones. The package-level helpers downstream consume only the
// supported subset; the rest stays here so the CLI surface is parity-
// clean (see PARITY_ROADMAP "Definition of 1:1").
type mpileupFlags struct {
	bamList        string
	fastaRef       string
	readGroups     string
	ignoreRG       bool
	outputPath     string
	outputType     string
	noVersion      bool
	threads        int
	countOrphans   bool
	maxDepth       int
	minMQ          int
	minBQ          int
	maxBQ          int
	deltaBQ        int
	ignoreOverlaps bool
	adjustMQ       int
	regions        string
	regionsFile    string
	targets        string
	targetsFile    string
	samples        string
	samplesFile    string
	flagAnyUnset   string
	flagAllUnset   string
	flagAnySet     string
	flagAllSet     string
	flagLS         string
	noBAQ          bool
	fullBAQ        bool
	redoBAQ        bool
	platforms      string
	perSampleMF    bool
	illumina13     bool
	config         string
	seed           int64
	ambigReads     string
	skipIndels     bool
	maxIDepth      int
	minIReads      int
	tandemQual     int
	extProb        int
	gapFrac        float64
	maxReadLen     int
	indelBias      float64
	indelSize      int
	indelsCNS      bool
	indels20       bool
	noIndelsCNS    bool
	openProb       int
	delBias        float64
	polyMQual      bool
	noPolyMQual    bool
	scoreVsRef     float64
	seqQOffset     int
	annotate       string
	gvcf           string
	noReference    bool
	writeIndex     string
	verbosity      int
	showHelp       bool
	showHelpAlt    bool
	showVersion    bool
}

// registerMpileupFlags binds every upstream-supported flag to its
// storage slot on fs. Keep one block per flag for readability; this
// table is the source of truth for the parity-locked CLI surface.
func registerMpileupFlags(fs *flag.FlagSet, mf *mpileupFlags) {
	cliflag.StringVar(fs, &mf.bamList, "b", "bam-list", "", "BAM list file")
	cliflag.StringVar(fs, &mf.fastaRef, "f", "fasta-ref", "", "Reference FASTA")
	cliflag.StringVar(fs, &mf.readGroups, "G", "read-groups", "", "Read-group filter file")
	// Upstream `--ignore-RG` (mpileup.c:1423) is long-only. The earlier
	// `-Z` short binding was an invention; reviewer caught it. Keep
	// the long form and the `--ignore-rg` lowercase alias.
	fs.BoolVar(&mf.ignoreRG, "ignore-RG", false, "Ignore @RG tags")
	fs.BoolVar(&mf.ignoreRG, "ignore-rg", false, "")

	cliflag.StringVar(fs, &mf.outputPath, "o", "output", "", "Output file")
	cliflag.StringVar(fs, &mf.outputType, "O", "output-type", "v", "Output type")
	fs.BoolVar(&mf.noVersion, "no-version", false, "Skip version line")
	fs.IntVar(&mf.threads, "threads", 0, "Threads (accepted, ignored)")

	cliflag.BoolVar(fs, &mf.countOrphans, "A", "count-orphans", false, "Keep orphan reads")
	cliflag.IntVar(fs, &mf.maxDepth, "d", "max-depth", 250, "Max depth per BAM")
	cliflag.IntVar(fs, &mf.minMQ, "q", "min-MQ", 0, "Minimum MAPQ")
	fs.IntVar(&mf.minMQ, "min-mq", 0, "")
	cliflag.IntVar(fs, &mf.minBQ, "Q", "min-BQ", 13, "Minimum base quality")
	fs.IntVar(&mf.minBQ, "min-bq", 13, "")
	fs.IntVar(&mf.maxBQ, "max-bq", 60, "Max base quality cap")
	fs.IntVar(&mf.maxBQ, "max-BQ", 60, "")
	fs.IntVar(&mf.deltaBQ, "delta-BQ", 0, "Delta base quality (accepted, ignored)")
	cliflag.BoolVar(fs, &mf.ignoreOverlaps, "x", "ignore-overlaps", false, "Disable overlap detection")
	cliflag.IntVar(fs, &mf.adjustMQ, "C", "adjust-MQ", 0, "Adjust MAPQ for excess mismatches")
	fs.IntVar(&mf.adjustMQ, "adjust-mq", 0, "")

	cliflag.StringVar(fs, &mf.regions, "r", "regions", "", "Region(s)")
	fs.StringVar(&mf.regions, "region", "", "")
	cliflag.StringVar(fs, &mf.regionsFile, "R", "regions-file", "", "Regions file")
	cliflag.StringVar(fs, &mf.targets, "t", "targets", "", "Targets (post-filter)")
	cliflag.StringVar(fs, &mf.targetsFile, "T", "targets-file", "", "Targets file")
	cliflag.StringVar(fs, &mf.samples, "s", "samples", "", "Samples list")
	cliflag.StringVar(fs, &mf.samplesFile, "S", "samples-file", "", "Samples file")

	fs.StringVar(&mf.flagAnyUnset, "skip-any-unset", "", "Skip reads with any flag unset")
	fs.StringVar(&mf.flagAnyUnset, "nu", "", "")
	fs.StringVar(&mf.flagAllUnset, "skip-all-unset", "", "Skip reads with all flags unset")
	fs.StringVar(&mf.flagAllUnset, "lu", "", "")
	fs.StringVar(&mf.flagAllUnset, "rf", "", "")
	fs.StringVar(&mf.flagAnySet, "skip-any-set", "", "Skip reads with any flag set")
	fs.StringVar(&mf.flagAnySet, "ns", "", "")
	fs.StringVar(&mf.flagAnySet, "ff", "", "")
	fs.StringVar(&mf.flagAllSet, "skip-all-set", "", "Skip reads with all flags set")
	fs.StringVar(&mf.flagLS, "ls", "", "Reads-listed filter")

	cliflag.BoolVar(fs, &mf.noBAQ, "B", "no-BAQ", false, "Disable BAQ")
	fs.BoolVar(&mf.noBAQ, "no-baq", false, "")
	cliflag.BoolVar(fs, &mf.fullBAQ, "D", "full-BAQ", false, "Full BAQ")
	fs.BoolVar(&mf.fullBAQ, "full-baq", false, "")
	cliflag.BoolVar(fs, &mf.redoBAQ, "E", "redo-BAQ", false, "Redo BAQ")
	fs.BoolVar(&mf.redoBAQ, "redo-baq", false, "")
	cliflag.StringVar(fs, &mf.platforms, "P", "platforms", "", "Platforms")
	cliflag.BoolVar(fs, &mf.perSampleMF, "p", "per-sample-mF", false, "Per-sample mF")
	fs.BoolVar(&mf.perSampleMF, "per-sample-mf", false, "")
	cliflag.BoolVar(fs, &mf.illumina13, "6", "illumina1.3+", false, "Illumina 1.3+ encoding")
	cliflag.StringVar(fs, &mf.config, "X", "config", "", "Predefined config")
	fs.Int64Var(&mf.seed, "seed", 0, "Random seed (accepted, ignored)")
	fs.StringVar(&mf.ambigReads, "ambig-reads", "", "Ambiguous-read handling")
	fs.StringVar(&mf.ambigReads, "ar", "", "")

	cliflag.BoolVar(fs, &mf.skipIndels, "I", "skip-indels", false, "Skip indels")
	cliflag.IntVar(fs, &mf.maxIDepth, "L", "max-idepth", 0, "Max indel depth")
	cliflag.IntVar(fs, &mf.minIReads, "m", "min-ireads", 0, "Min indel reads")
	cliflag.IntVar(fs, &mf.tandemQual, "h", "tandem-qual", 0, "Tandem-quality penalty")
	cliflag.IntVar(fs, &mf.extProb, "e", "ext-prob", 0, "Gap extension probability")
	cliflag.Float64Var(fs, &mf.gapFrac, "F", "gap-frac", 0, "Gap-open fraction")
	cliflag.IntVar(fs, &mf.maxReadLen, "M", "max-read-len", 0, "Max read length")
	fs.Float64Var(&mf.indelBias, "indel-bias", 0, "Indel bias (accepted, ignored)")
	fs.IntVar(&mf.indelSize, "indel-size", 0, "Indel size (accepted, ignored)")
	fs.BoolVar(&mf.indelsCNS, "indels-cns", false, "Indels via consensus (accepted, ignored)")
	fs.BoolVar(&mf.indels20, "indels-2.0", false, "Indels 2.0 (accepted, ignored)")
	fs.BoolVar(&mf.noIndelsCNS, "no-indels-cns", false, "No indels-cns (accepted, ignored)")
	fs.IntVar(&mf.openProb, "open-prob", 0, "Gap-open probability (accepted, ignored)")
	fs.Float64Var(&mf.delBias, "del-bias", 0, "Del bias (accepted, ignored, hidden upstream)")
	fs.BoolVar(&mf.polyMQual, "poly-mqual", false, "Poly MQ (accepted, ignored)")
	fs.BoolVar(&mf.noPolyMQual, "no-poly-mqual", false, "No-poly MQ (accepted, ignored)")
	fs.Float64Var(&mf.scoreVsRef, "score-vs-ref", 0, "Score-vs-ref (accepted, ignored)")
	fs.IntVar(&mf.seqQOffset, "seqq-offset", 0, "SeqQ offset (accepted, ignored)")

	cliflag.StringVar(fs, &mf.annotate, "a", "annotate", "", "FORMAT/INFO tags to add")
	cliflag.StringVar(fs, &mf.gvcf, "g", "gvcf", "", "gVCF block spec")
	fs.BoolVar(&mf.noReference, "no-reference", false, "No reference (accepted, ignored)")
	cliflag.StringVar(fs, &mf.writeIndex, "W", "write-index", "", "Auto-index output (accepted, ignored)")
	cliflag.IntVar(fs, &mf.verbosity, "v", "verbosity", 0, "Verbosity (accepted, ignored)")
	fs.BoolVar(&mf.showHelp, "?", false, "")
	fs.BoolVar(&mf.showHelp, "help", false, "")
	fs.BoolVar(&mf.showVersion, "version", false, "")
}

// runMpileup is the bcftools mpileup entry point.
func runMpileup(args []string) int {
	fs := flag.NewFlagSet("bcftools mpileup", flag.ContinueOnError)
	fs.SetOutput(io.Discard)

	mf := &mpileupFlags{}
	registerMpileupFlags(fs, mf)

	if err := fs.Parse(args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, mpileupUsage)
		return 2
	}
	if mf.showHelp || mf.showHelpAlt {
		fmt.Print(mpileupUsage)
		return 0
	}
	if mf.showVersion {
		fmt.Println(version)
		return 0
	}

	if deferred := checkMpileupDeferred(mf); deferred != "" {
		fmt.Fprintf(os.Stderr, "bcftools mpileup: %s is not implemented in v1; tracked in docs/PARITY_ROADMAP.md#bcftools\n", deferred)
		return 2
	}

	rest := fs.Args()
	if len(rest) == 0 && mf.bamList == "" {
		fmt.Fprintln(os.Stderr, "bcftools mpileup: at least one BAM input is required (positional or via -b/--bam-list)")
		fmt.Fprint(os.Stderr, mpileupUsage)
		return 2
	}
	if mf.fastaRef == "" {
		fmt.Fprintln(os.Stderr, "bcftools mpileup: -f/--fasta-ref is required")
		return 2
	}

	format, err := bcftools.ParseOutputFormat(mf.outputType)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}

	opts := bcftools.MpileupOptions{
		Inputs:         rest,
		FastaRef:       mf.fastaRef,
		BamList:        mf.bamList,
		MaxDepth:       mf.maxDepth,
		MinMQ:          uint8(clampU8(mf.minMQ)),
		MinBQ:          uint8(clampU8(mf.minBQ)),
		MaxBQ:          uint8(clampU8(mf.maxBQ)),
		CountOrphans:   mf.countOrphans,
		IgnoreOverlaps: mf.ignoreOverlaps,
		NoBAQ:          mf.noBAQ,
		RedoBAQ:        mf.redoBAQ,
		FullBAQ:        mf.fullBAQ,
		AdjustMQ:       mf.adjustMQ,
		Annotate:       mf.annotate,
		ReadGroups:     mf.readGroups,
		IgnoreRG:       mf.ignoreRG,
		Platforms:      mf.platforms,
		Config:         mf.config,
		PerSampleMF:    mf.perSampleMF,
		Seed:           mf.seed,
		TandemQual:     mf.tandemQual,
		ExtProb:        mf.extProb,
		GapFrac:        mf.gapFrac,
		OpenProb:       mf.openProb,
		IndelBias:      mf.indelBias,
		IndelSize:      mf.indelSize,
		MinIReads:      mf.minIReads,
		MaxIDepth:      mf.maxIDepth,
		ARProb:         0,
		AmbigReads:     mf.ambigReads,
		MaxReadLen:     mf.maxReadLen,
		DelBias:        mf.delBias,
		PolyMQual:      mf.polyMQual,
		ScoreVsRef:     mf.scoreVsRef,
		SeqQOffset:     mf.seqQOffset,
		SkipIndels:     mf.skipIndels,
		IndelsCNS:      mf.indelsCNS,
		NoIndelsCNS:    mf.noIndelsCNS,
		GVCFBlock:      mf.gvcf,
		NoReference:    mf.noReference,
		OutputFormat:   format,
		Output:         mf.outputPath,
		Threads:        mf.threads,
		NoVersion:      mf.noVersion,
		Verbosity:      mf.verbosity,
		FlagIncl:       mf.flagAllUnset,
		FlagExcl:       mf.flagAnySet,
		FlagAny:        mf.flagAnyUnset,
		FlagLS:         mf.flagLS,
	}
	if mf.regions != "" {
		opts.Regions = bcftools.SplitCommaList(mf.regions)
	}
	if mf.regionsFile != "" {
		opts.RegionsFile = mf.regionsFile
	}
	if mf.targets != "" {
		opts.Targets = bcftools.SplitCommaList(mf.targets)
	}
	if mf.targetsFile != "" {
		opts.TargetsFile = mf.targetsFile
	}
	if mf.samples != "" {
		opts.Samples = bcftools.SplitCommaList(mf.samples)
	}
	if mf.samplesFile != "" {
		opts.SamplesFile = mf.samplesFile
	}

	out, err := openOutFile(mf.outputPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "bcftools mpileup: %v\n", err)
		return 1
	}
	defer out.Close()

	// `-O z` wraps the destination in a gzip.Writer. Previously the
	// validator accepted "z" but neither the library nor the CLI
	// actually compressed, silently emitting plain VCF — reviewer
	// caught it.
	var writer io.WriteCloser = out
	if mf.outputType == "z" {
		gz := gzip.NewWriter(out)
		// Close the gzip layer FIRST (flushes trailer) before the
		// underlying file is closed by the deferred out.Close().
		defer gz.Close()
		writer = gz
	}

	if err := bcftools.MpileupFile(opts, writer); err != nil {
		fmt.Fprintf(os.Stderr, "bcftools mpileup: %v\n", err)
		return 1
	}
	return 0
}

// checkMpileupDeferred returns a non-empty string when the caller has
// explicitly set a flag whose underlying behaviour is deferred to v2.
// Returning "" means "all-OK, the v1 model accepts these settings".
//
// The full list of accepted-but-inert flags lives in PARITY_ROADMAP
// under the bcftools mpileup option-tail gap section; we only
// hard-reject the small set whose set-but-ignored behaviour would
// be silently wrong.
func checkMpileupDeferred(mf *mpileupFlags) string {
	if mf.redoBAQ {
		return "-E/--redo-BAQ"
	}
	switch mf.outputType {
	case "", "v", "z":
		// OK.
	case "u", "b":
		return "-O " + mf.outputType + " (BCF output)"
	default:
		return "-O " + mf.outputType
	}
	return ""
}

// clampU8 truncates n into the uint8 range [0,255], returning 0 for
// negative values. Used for the -q/-Q/--max-bq flags which the CLI
// parses as ints to preserve upstream's "negative is invalid" error
// messages but the library wants as uint8.
func clampU8(n int) int {
	if n < 0 {
		return 0
	}
	if n > 255 {
		return 255
	}
	return n
}

// Ensure strconv is referenced for the eventual --gvcf parser; for v1
// we ignore the block argument but keep the import so a follow-up that
// honours block syntax (e.g. "5,15,30") compiles cleanly.
var _ = strconv.Atoi
