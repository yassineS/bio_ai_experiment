// Command bcftools is a pure-Go reimplementation of selected bcftools
// subcommands. Today it ships `view`, `index`, `stats`, `query`, `concat`,
// `norm`, `call`, `merge`, `isec`, `sort`, `head`, `reheader`,
// `annotate`, `convert`, `mendelian`, `mendelian2`, `gtcheck`, `roh`,
// `filter`, `consensus`, `polysomy`, `cnv`, `csq`, `mpileup`, and the
// `plugin` subprocess plugin system (also reachable as `bcftools +<name>`).
package main

import (
	"flag"
	"fmt"
	"io"
	"os"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/tools/bcftools/pkg/bcftools"
)

const version = "0.1.0"

const rootUsage = `bcftools - pure-Go reimplementation.

Usage:
  bcftools <subcommand> [options]

Subcommands:
  view      Print, filter, or convert VCF/BCF records.
  query     Format-string output for VCF/BCF records.
  concat    Concatenate VCF/BCF files.
  norm      Left-align indels, split/join multiallelics, drop duplicates.
  call      Variant calling from per-position genotype likelihoods.
  merge     Combine multiple per-sample VCFs into a multi-sample VCF.
  isec      Set operations on N VCFs (intersection / union / complement).
  sort      Sort VCF/BCF by (CHROM, POS) following the header contig order.
  head      Print just the header (or a slice of it).
  reheader  Replace header in-place (sample renames / FAI contigs).
  annotate  Annotate INFO from a tab-indexed table or VCF.
  convert   Re-emit VCF/BCF in a different format (with sample/region filters).
  mendelian Detect Mendelian-inconsistent genotypes given a PED-style trio.
  mendelian2 Newer Mendelian-inheritance checker with PED-file ingestion.
  gtcheck   Check sample identity (genotype concordance).
  roh       Detect runs of autozygosity (ROH) via a 2-state HMM.
  filter    Soft-filter records by include / exclude expression.
  consensus Apply VCF variants to a reference FASTA.
  polysomy  Detect chromosomal copy number from B-allele frequency.
  cnv       Copy-number variation caller (v1: heuristic CN-call).
  csq       Predict variant consequences against a GFF (v1: SNPs only).
  mpileup   Per-position genotype likelihoods from BAM (v1: SNPs only).
  index     Build a CSI (or .tbi) index for a BCF / VCF.gz file.
  stats     Produce summary statistics from VCF/BCF (plot-vcfstats compatible).
  plugin    Run a user plugin (subprocess); also reachable as +<name>.
  help      Show this help (also via -? on subcommands).
  version   Show version.
`

func main() {
	if len(os.Args) < 2 {
		fmt.Fprint(os.Stderr, rootUsage)
		os.Exit(1)
	}
	// `bcftools +<name> ...` is shorthand for `bcftools plugin <name> ...`.
	if len(os.Args[1]) > 1 && os.Args[1][0] == '+' {
		os.Exit(runPlugin(os.Args[2:], os.Args[1][1:]))
	}
	switch os.Args[1] {
	case "view":
		os.Exit(runView(os.Args[2:]))
	case "query":
		os.Exit(runQuery(os.Args[2:]))
	case "concat":
		os.Exit(runConcat(os.Args[2:]))
	case "norm":
		os.Exit(runNorm(os.Args[2:]))
	case "call":
		os.Exit(runCall(os.Args[2:]))
	case "index":
		os.Exit(runIndex(os.Args[2:]))
	case "stats":
		os.Exit(runStatsCmd(os.Args[2:]))
	case "merge":
		os.Exit(runMerge(os.Args[2:]))
	case "isec":
		os.Exit(runIsec(os.Args[2:]))
	case "sort":
		os.Exit(runSort(os.Args[2:]))
	case "head":
		os.Exit(runHead(os.Args[2:]))
	case "reheader":
		os.Exit(runReheader(os.Args[2:]))
	case "annotate":
		os.Exit(runAnnotate(os.Args[2:]))
	case "convert":
		os.Exit(runConvert(os.Args[2:]))
	case "mendelian":
		os.Exit(runMendelian(os.Args[2:]))
	case "mendelian2":
		os.Exit(runMendelian2(os.Args[2:]))
	case "gtcheck":
		os.Exit(runGtcheck(os.Args[2:]))
	case "roh":
		os.Exit(runRoh(os.Args[2:]))
	case "filter":
		os.Exit(runFilter(os.Args[2:]))
	case "consensus":
		os.Exit(runConsensus(os.Args[2:]))
	case "polysomy":
		os.Exit(runPolysomy(os.Args[2:]))
	case "cnv":
		os.Exit(runCNV(os.Args[2:]))
	case "csq":
		os.Exit(runCSQ(os.Args[2:]))
	case "mpileup":
		os.Exit(runMpileup(os.Args[2:]))
	case "plugin":
		os.Exit(runPlugin(os.Args[2:], ""))
	case "help", "--help":
		fmt.Print(rootUsage)
		return
	case "version", "-v", "--version":
		fmt.Println(version)
		return
	default:
		fmt.Fprintf(os.Stderr, "bcftools: unknown subcommand %q\n", os.Args[1])
		fmt.Fprint(os.Stderr, rootUsage)
		os.Exit(1)
	}
}

const indexUsage = `bcftools index - build a CSI (or .tbi) index for a BCF / VCF.gz file.

Usage:
  bcftools index [options] <in.bcf|in.vcf.gz>

Options:
  -c, --csi                Emit CSI (default; required for BCF).
  -t, --tbi                Emit .tbi (VCF.gz only).
      --csi-min-shift N    Use a non-default min_shift for the CSI bin scheme (default 14).
  -o, --output PATH        Output index path (default <in>.csi or <in>.tbi).
  -f, --force              Overwrite an existing index file.
  -@, --threads N          Accepted; v1 is single-threaded.
  -h, --help               Show this help.
      --version            Show version.
`

func runIndex(args []string) int {
	fs := flag.NewFlagSet("bcftools index", flag.ContinueOnError)
	fs.SetOutput(io.Discard)

	var (
		csiFlag     bool
		tbiFlag     bool
		minShift    int
		outPath     string
		force       bool
		threads     int
		showHelp    bool
		showVersion bool
	)
	cliflag.BoolVar(fs, &csiFlag, "c", "csi", false, "Emit CSI")
	cliflag.BoolVar(fs, &tbiFlag, "t", "tbi", false, "Emit TBI (VCF.gz only)")
	fs.IntVar(&minShift, "csi-min-shift", 0, "CSI min_shift")
	cliflag.StringVar(fs, &outPath, "o", "output", "", "Output index path")
	cliflag.BoolVar(fs, &force, "f", "force", false, "Overwrite existing index")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "Threads (accepted, ignored)")
	fs.BoolVar(&showHelp, "h", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVersion, "version", false, "")

	if err := parseFlags(fs, args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, indexUsage)
		return 2
	}
	if showHelp {
		fmt.Print(indexUsage)
		return 0
	}
	if showVersion {
		fmt.Println(version)
		return 0
	}

	rest := fs.Args()
	if len(rest) == 0 {
		fmt.Fprintln(os.Stderr, "bcftools index: missing input file")
		fmt.Fprint(os.Stderr, indexUsage)
		return 2
	}
	input := rest[0]

	format := bcftools.IndexCSI
	if tbiFlag && !csiFlag {
		format = bcftools.IndexTBI
	}

	opts := bcftools.IndexOptions{
		Format:     format,
		MinShift:   int32(minShift),
		OutputPath: outPath,
		Force:      force,
	}
	written, err := bcftools.BuildIndex(input, opts)
	if err != nil {
		fmt.Fprintf(os.Stderr, "bcftools index: %v\n", err)
		return 1
	}
	_ = written
	return 0
}

const viewUsage = `bcftools view - print, filter, or convert VCF/BCF records.

Usage:
  bcftools view [options] <in.vcf[.gz]|in.bcf> [region ...]

Options:
  -O, --output-type {v|z|u|b}     v=VCF (default), z=VCF.gz, u=uncompressed BCF, b=compressed BCF.
  -o, --output PATH               Output file (default stdout).
  -h, --header-only               Emit only the header.
  -H, --no-header                 Drop the header.
      --with-header               Print both header and records (default).
      --no-version                Do not append provenance lines.
  -G, --drop-genotypes            Drop FORMAT and per-sample columns.
  -A, --trim-unseen-allele [N]    Remove <*>/<NON_REF> (-A) or at all sites (-AA).
  -a, --trim-alt-alleles          Trim ALTs not seen in any remaining sample's GT.
  -c, --min-ac N                  Minimum allele count.
  -C, --max-ac N                  Maximum allele count.
  -q, --min-af FLOAT              Minimum allele frequency.
  -Q, --max-af FLOAT              Maximum allele frequency.
  -g, --genotype [^]hom|het|miss  Keep sites with at least one matching GT.
  -i, --include EXPR              Keep records matching expression.
  -e, --exclude EXPR              Drop records matching expression.
  -f, --apply-filters NAME[,..]   Keep only PASS or named filters.
  -k, --known                     Keep sites with non-'.' ID.
  -n, --novel                     Keep sites with '.' ID.
  -m, --min-alleles N             Minimum REF+ALT count.
  -M, --max-alleles N             Maximum REF+ALT count.
  -p, --phased                    Keep sites where all called GTs are phased.
  -P, --exclude-phased            Drop such sites.
  -u, --uncalled                  Keep sites without a called genotype.
  -U, --exclude-uncalled          Drop such sites.
  -x, --private                   (Deferred — rejection-parity.)
  -X, --exclude-private           (Deferred — rejection-parity.)
  -r, --regions chr[:beg-end[,..]] Region(s) (uses .tbi when available).
  -R, --regions-file PATH         BED-like regions file.
      --regions-overlap 0|1|2     Region inclusion rule [1].
  -t, --targets chr[:beg-end[,..]] Like -r but always a post-filter.
  -T, --targets-file PATH         BED-like targets file (post-filter).
      --targets-overlap 0|1|2     Target inclusion rule [0].
  -s, --samples LIST              Restrict to these samples (prefix with '^' to exclude).
  -S, --samples-file PATH         File with sample IDs (prefix with '^' to exclude).
      --force-samples             Warn about unknown subset samples instead of failing.
  -I, --no-update                 Do not (re)calculate INFO/AC and INFO/AN.
  -v, --types LIST                Include comma-list of variant types.
  -V, --exclude-types LIST        Exclude variant types.
  -W, --write-index[=FMT]         Auto-index output (csi|tbi).
      --verbosity INT             Verbosity level.
  -l, --compression-level N       gzip level for z output.
  -@, --threads N                 Number of BGZF compression threads.
  -?, --help                      Show this help.
      --version                   Show version.

Note on -h:
  Upstream bcftools uses -h for "header-only" rather than "help" — we follow
  that convention here. Use -? or --help for help.
`

func runView(args []string) int {
	fs := flag.NewFlagSet("bcftools view", flag.ContinueOnError)
	fs.SetOutput(io.Discard) // we print usage ourselves

	var (
		outputType      string
		outputPath      string
		headerOnly      bool
		noHeader        bool
		withHeader      bool
		dropGT          bool
		minAC           int
		maxAC           int
		minAF           float64
		maxAF           float64
		includeExpr     string
		excludeExpr     string
		applyFilters    string
		regions         string
		regionsFile     string
		regionsOverlap  int
		targets         string
		targetsFile     string
		targetsOverlap  int
		samples         string
		samplesFile     string
		forceSamples    bool
		compressLevel   int
		threads         int
		includeTypes    string
		excludeTypes    string
		noUpdate        bool
		trimAlts        bool
		trimUnseen      int
		genotypeFilter  string
		knownOnly       bool
		novelOnly       bool
		minAlleles      int
		maxAlleles      int
		phasedOnly      bool
		excludePhased   bool
		uncalledOnly    bool
		excludeUncalled bool
		privateOnly     bool
		excludePrivate  bool
		writeIndex      string
		verbosity       int
		showHelp        bool
		showVer         bool
	)
	cliflag.StringVar(fs, &outputType, "O", "output-type", "v", "Output type")
	cliflag.StringVar(fs, &outputPath, "o", "output", "", "Output path")
	cliflag.BoolVar(fs, &headerOnly, "h", "header-only", false, "Header only")
	cliflag.BoolVar(fs, &noHeader, "H", "no-header", false, "Drop header")
	fs.BoolVar(&withHeader, "with-header", false, "Print both header and records (default)")
	cliflag.BoolVar(fs, &dropGT, "G", "drop-genotypes", false, "Drop genotypes")
	cliflag.IntVar(fs, &minAC, "c", "min-ac", 0, "Min allele count")
	cliflag.IntVar(fs, &maxAC, "C", "max-ac", 0, "Max allele count")
	cliflag.Float64Var(fs, &minAF, "q", "min-af", 0, "Min allele frequency")
	cliflag.Float64Var(fs, &maxAF, "Q", "max-af", 0, "Max allele frequency")
	cliflag.StringVar(fs, &includeExpr, "i", "include", "", "Include expression")
	cliflag.StringVar(fs, &excludeExpr, "e", "exclude", "", "Exclude expression")
	cliflag.StringVar(fs, &applyFilters, "f", "apply-filters", "", "Filter list")
	cliflag.StringVar(fs, &regions, "r", "regions", "", "Region(s)")
	cliflag.StringVar(fs, &regionsFile, "R", "regions-file", "", "Regions file")
	fs.IntVar(&regionsOverlap, "regions-overlap", 1, "Region inclusion rule (0|1|2)")
	cliflag.StringVar(fs, &targets, "t", "targets", "", "Targets (post-filter)")
	cliflag.StringVar(fs, &targetsFile, "T", "targets-file", "", "Targets file")
	fs.IntVar(&targetsOverlap, "targets-overlap", 0, "Target inclusion rule (0|1|2)")
	cliflag.StringVar(fs, &samples, "s", "samples", "", "Samples (prefix with '^' to exclude)")
	cliflag.StringVar(fs, &samplesFile, "S", "samples-file", "", "Samples file (prefix with '^' to exclude)")
	fs.BoolVar(&forceSamples, "force-samples", false, "Only warn about unknown subset samples")
	cliflag.IntVar(fs, &compressLevel, "l", "compression-level", -1, "gzip level")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "Threads (accepted, ignored)")
	cliflag.StringVar(fs, &includeTypes, "v", "types", "", "Include only listed variant types (snps,indels,mnps,bnd,other)")
	cliflag.StringVar(fs, &excludeTypes, "V", "exclude-types", "", "Exclude listed variant types")
	cliflag.BoolVar(fs, &noUpdate, "I", "no-update", false, "Do not (re)calculate INFO fields for the subset (currently INFO/AC and INFO/AN)")
	cliflag.BoolVar(fs, &trimAlts, "a", "trim-alt-alleles", false, "Trim ALT alleles not seen in genotypes")
	cliflag.IntVar(fs, &trimUnseen, "A", "trim-unseen-allele", 0, "Remove '<*>'/'<NON_REF>' (1) or at all sites (2)")
	cliflag.StringVar(fs, &genotypeFilter, "g", "genotype", "", "Require [^]hom|het|miss genotype on any sample")
	cliflag.BoolVar(fs, &knownOnly, "k", "known", false, "Keep only sites with non-'.' ID")
	cliflag.BoolVar(fs, &novelOnly, "n", "novel", false, "Keep only sites with '.' ID")
	cliflag.IntVar(fs, &minAlleles, "m", "min-alleles", 0, "Minimum number of alleles (REF+ALT)")
	cliflag.IntVar(fs, &maxAlleles, "M", "max-alleles", 0, "Maximum number of alleles (REF+ALT)")
	cliflag.BoolVar(fs, &phasedOnly, "p", "phased", false, "Select sites where all called genotypes are phased")
	cliflag.BoolVar(fs, &excludePhased, "P", "exclude-phased", false, "Exclude sites where all called genotypes are phased")
	cliflag.BoolVar(fs, &uncalledOnly, "u", "uncalled", false, "Select sites without a called genotype")
	cliflag.BoolVar(fs, &excludeUncalled, "U", "exclude-uncalled", false, "Exclude sites without a called genotype")
	cliflag.BoolVar(fs, &privateOnly, "x", "private", false, "Select sites where non-ref alleles are private to subset samples (deferred)")
	cliflag.BoolVar(fs, &excludePrivate, "X", "exclude-private", false, "Exclude sites where non-ref alleles are private to subset samples (deferred)")
	cliflag.StringVar(fs, &writeIndex, "W", "write-index", "", "Automatically index output (csi|tbi)")
	fs.IntVar(&verbosity, "verbosity", 0, "Verbosity level")
	fs.BoolVar(&showHelp, "?", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "version", false, "")

	// -W/--write-index accepts a bare form (no value); expand to upstream's default.
	args = preprocessOptionalArg(args, "-W", "csi")
	args = preprocessOptionalArg(args, "--write-index", "csi")
	// -A is upstream's counting flag: bare `-A` ≡ 1 (record-level trim),
	// `-AA` ≡ 2 (trim everywhere). Collapse to `-A=N`.
	args = collapseRepeatedShortFlag(args, "-A", "trim-unseen-allele")
	args = preprocessFlagOrBare(args, "-A", "1")
	args = preprocessFlagOrBare(args, "--trim-unseen-allele", "1")

	if err := parseFlags(fs, args); err != nil {
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
		fmt.Fprintln(os.Stderr, "bcftools view: missing input file")
		fmt.Fprint(os.Stderr, viewUsage)
		return 2
	}
	input := positional[0]
	cliRegions := positional[1:]

	format, err := bcftools.ParseOutputFormat(outputType)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}

	if includeTypes != "" && excludeTypes != "" {
		fmt.Fprintln(os.Stderr, "bcftools view: only one of -v/--types or -V/--exclude-types can be given")
		return 2
	}
	if knownOnly && novelOnly {
		fmt.Fprintln(os.Stderr, "bcftools view: only one of -k/--known and -n/--novel can be given")
		return 2
	}
	if phasedOnly && excludePhased {
		fmt.Fprintln(os.Stderr, "bcftools view: only one of -p/--phased and -P/--exclude-phased can be given")
		return 2
	}
	if uncalledOnly && excludeUncalled {
		fmt.Fprintln(os.Stderr, "bcftools view: only one of -u/--uncalled and -U/--exclude-uncalled can be given")
		return 2
	}
	if regionsOverlap < 0 || regionsOverlap > 2 {
		fmt.Fprintf(os.Stderr, "bcftools view: --regions-overlap must be 0, 1 or 2 (got %d)\n", regionsOverlap)
		return 2
	}
	if targetsOverlap < 0 || targetsOverlap > 2 {
		fmt.Fprintf(os.Stderr, "bcftools view: --targets-overlap must be 0, 1 or 2 (got %d)\n", targetsOverlap)
		return 2
	}
	if writeIndex != "" && writeIndex != "csi" && writeIndex != "tbi" {
		fmt.Fprintf(os.Stderr, "bcftools view: --write-index must be csi or tbi (got %q)\n", writeIndex)
		return 2
	}
	if privateOnly || excludePrivate {
		// Rejection-parity: requires sample-subset bookkeeping that's
		// not yet wired into the view loop. Queued for follow-up.
		fmt.Fprintln(os.Stderr, "bcftools view: -x/-X --private/--exclude-private is to be implemented, please open an issue on github")
		return 1
	}
	if withHeader && headerOnly {
		fmt.Fprintln(os.Stderr, "bcftools view: --with-header conflicts with --header-only")
		return 2
	}
	// `-A` / `--trim-unseen-allele`: 0=off, 1=record-level, 2=trim everywhere.
	if trimUnseen < 0 || trimUnseen > 2 {
		fmt.Fprintf(os.Stderr, "bcftools view: -A may appear at most twice (got count=%d)\n", trimUnseen)
		return 2
	}

	noVersionFlag := fs.Lookup("no-version")
	noVersion := noVersionFlag != nil && noVersionFlag.Value.String() == "true"

	// Parse `^`-prefixed sample lists.
	samplesExclude := false
	if strings.HasPrefix(samples, "^") {
		samples = samples[1:]
		samplesExclude = true
	}
	samplesFileExclude := false
	sf := samplesFile
	if strings.HasPrefix(sf, "^") {
		sf = sf[1:]
		samplesFileExclude = true
	}

	opts := bcftools.ViewOptions{
		OutputFormat:        format,
		HeaderOnly:          headerOnly,
		NoHeader:            noHeader,
		DropGenotypes:       dropGT,
		MinAlleleCount:      minAC,
		MaxAlleleCount:      maxAC,
		MinAlleleFreq:       minAF,
		MaxAlleleFreq:       maxAF,
		IncludeExpr:         includeExpr,
		ExcludeExpr:         excludeExpr,
		ApplyFilters:        bcftools.SplitCommaList(applyFilters),
		CompressLevel:       compressLevel,
		IncludeTypes:        bcftools.SplitCommaList(includeTypes),
		ExcludeTypes:        bcftools.SplitCommaList(excludeTypes),
		NoUpdateINFO:        noUpdate,
		Threads:             threads,
		TrimAltAlleles:      trimAlts,
		TrimUnseenAllele:    trimUnseen >= 1,
		TrimUnseenAlleleAll: trimUnseen >= 2,
		ForceSamples:        forceSamples,
		SamplesExclude:      samplesExclude || samplesFileExclude,
		GenotypeFilter:      genotypeFilter,
		Known:               knownOnly,
		Novel:               novelOnly,
		MinAlleles:          minAlleles,
		MaxAlleles:          maxAlleles,
		Phased:              phasedOnly,
		ExcludePhased:       excludePhased,
		Uncalled:            uncalledOnly,
		ExcludeUncalled:     excludeUncalled,
		RegionsOverlap:      regionsOverlap,
		TargetsOverlap:      targetsOverlap,
		WriteIndex:          writeIndex,
		Verbosity:           verbosity,
		NoVersion:           noVersion,
		PGCommand:           strings.Join(os.Args, " "),
	}
	if regions != "" {
		opts.Regions = bcftools.SplitCommaList(regions)
	}
	opts.Regions = append(opts.Regions, cliRegions...)
	if targets != "" {
		opts.Targets = bcftools.SplitCommaList(targets)
	}
	if regionsFile != "" {
		regs, err := bcftools.LoadRegionsFile(regionsFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "bcftools view: %v\n", err)
			return 1
		}
		opts.Regions = append(opts.Regions, regs...)
		opts.RegionsFile = regionsFile
	}
	if targetsFile != "" {
		regs, err := bcftools.LoadRegionsFile(targetsFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "bcftools view: %v\n", err)
			return 1
		}
		opts.Targets = append(opts.Targets, regs...)
		opts.TargetsFile = targetsFile
	}
	if samples != "" {
		opts.Samples = bcftools.SplitCommaList(samples)
	}
	if sf != "" {
		names, err := bcftools.LoadSamplesFile(sf)
		if err != nil {
			fmt.Fprintf(os.Stderr, "bcftools view: %v\n", err)
			return 1
		}
		opts.Samples = append(opts.Samples, names...)
		opts.SamplesFile = sf
	}

	out, err := openOutFile(outputPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "bcftools view: %v\n", err)
		return 1
	}
	defer out.Close()

	if _, err := bcftools.ViewFile(input, out, opts, os.Stderr); err != nil {
		fmt.Fprintf(os.Stderr, "bcftools view: %v\n", err)
		return 1
	}
	return 0
}

// openOutFile mirrors the helper used by samtools view: stdout when path is
// empty or "-", otherwise create the file.
func openOutFile(path string) (io.WriteCloser, error) {
	if path == "" || path == "-" {
		return nopCloser{os.Stdout}, nil
	}
	return os.Create(path)
}

type nopCloser struct{ io.Writer }

func (nopCloser) Close() error { return nil }

const statsUsage = `bcftools stats - produce summary statistics from VCF/BCF.

Usage:
  bcftools stats [options] <in.vcf[.gz]|in.bcf>

Options:
  -s, --samples LIST              Restrict to these samples (comma list).
  -S, --samples-file PATH         File with sample IDs (one per line).
  -r, --regions chr[:beg-end[,..]] Region(s) — applied as a post-filter (no index).
  -R, --regions-file PATH         BED-like regions file.
  -t, --targets chr[:beg-end[,..]] Like -r but always a post-filter.
  -T, --targets-file PATH         BED-like targets file.
  -i, --include EXPR              Keep records matching expression.
  -e, --exclude EXPR              Drop records matching expression.
  -f, --apply-filters NAME[,..]   Keep only PASS or named filters.
  -d, --depth MIN,MAX,STEP        Depth-distribution bins (default 0,500,1).
  -a, --af-bins LIST              Allele-frequency bin edges (default 0,0.1,...,0.9,0.99,1.0).
  -c, --collapse {none|snps|indels|both|all|some|id}
                                  Multi-allelic site collapse rule. (Accepted; v1 always treats each ALT separately.)
  -1, --1st-allele-only           Count only the 1st ALT allele.
      --af-tag TAG                INFO tag to read AF from (default: compute from GT).
  -o, --output PATH               Output file (default stdout).
      --threads N                 Accepted; v1 is single-threaded.
  -h, --help                      Show this help.
      --version                   Show version.
`

func runStatsCmd(args []string) int {
	fs := flag.NewFlagSet("bcftools stats", flag.ContinueOnError)
	fs.SetOutput(io.Discard)

	var (
		samples         string
		samplesFile     string
		regions         string
		regionsFile     string
		targets         string
		targetsFile     string
		includeExpr     string
		excludeExpr     string
		applyFilters    string
		depthSpec       string
		afBinsSpec      string
		collapse        string
		firstAlleleOnly bool
		afTag           string
		outputPath      string
		threads         int
		showHelp        bool
		showVersion     bool
	)
	cliflag.StringVar(fs, &samples, "s", "samples", "", "Samples")
	cliflag.StringVar(fs, &samplesFile, "S", "samples-file", "", "Samples file")
	cliflag.StringVar(fs, &regions, "r", "regions", "", "Region(s)")
	cliflag.StringVar(fs, &regionsFile, "R", "regions-file", "", "Regions file")
	cliflag.StringVar(fs, &targets, "t", "targets", "", "Targets (post-filter)")
	cliflag.StringVar(fs, &targetsFile, "T", "targets-file", "", "Targets file")
	cliflag.StringVar(fs, &includeExpr, "i", "include", "", "Include expression")
	cliflag.StringVar(fs, &excludeExpr, "e", "exclude", "", "Exclude expression")
	cliflag.StringVar(fs, &applyFilters, "f", "apply-filters", "", "Filter list")
	cliflag.StringVar(fs, &depthSpec, "d", "depth", "", "Depth spec MIN,MAX,STEP")
	cliflag.StringVar(fs, &afBinsSpec, "a", "af-bins", "", "AF bin edges")
	cliflag.StringVar(fs, &collapse, "c", "collapse", "none", "Collapse rule")
	cliflag.BoolVar(fs, &firstAlleleOnly, "1", "1st-allele-only", false, "1st ALT only")
	fs.StringVar(&afTag, "af-tag", "", "INFO AF tag")
	cliflag.StringVar(fs, &outputPath, "o", "output", "", "Output path")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "Threads (accepted, ignored)")
	fs.BoolVar(&showHelp, "h", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVersion, "version", false, "")

	if err := parseFlags(fs, args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, statsUsage)
		return 2
	}
	if showHelp {
		fmt.Print(statsUsage)
		return 0
	}
	if showVersion {
		fmt.Println(version)
		return 0
	}

	rest := fs.Args()
	if len(rest) == 0 {
		fmt.Fprintln(os.Stderr, "bcftools stats: missing input file")
		fmt.Fprint(os.Stderr, statsUsage)
		return 2
	}
	input := rest[0]

	opts := bcftools.StatsOptions{
		IncludeExpr:     includeExpr,
		ExcludeExpr:     excludeExpr,
		ApplyFilters:    bcftools.SplitCommaList(applyFilters),
		Collapse:        collapse,
		FirstAlleleOnly: firstAlleleOnly,
		AFTag:           afTag,
		InputFile:       input,
	}
	// Upstream gates the per-sample sections (PSC/PSI/HWE) and the DP
	// genotype histogram on whether -s/-S was supplied at all, regardless
	// of value (`-s -` selects every sample). Detect that via fs.Visit so
	// even an empty `-s ''` counts as "given".
	fs.Visit(func(f *flag.Flag) {
		switch f.Name {
		case "s", "samples", "S", "samples-file":
			opts.SamplesGiven = true
		}
	})
	if samples != "" {
		opts.Samples = bcftools.SplitCommaList(samples)
	}
	if samplesFile != "" {
		names, err := bcftools.LoadSamplesFile(samplesFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "bcftools stats: %v\n", err)
			return 1
		}
		opts.Samples = append(opts.Samples, names...)
		opts.SamplesFile = samplesFile
	}
	if regions != "" {
		opts.Regions = bcftools.SplitCommaList(regions)
	}
	if regionsFile != "" {
		regs, err := bcftools.LoadRegionsFile(regionsFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "bcftools stats: %v\n", err)
			return 1
		}
		opts.Regions = append(opts.Regions, regs...)
		opts.RegionsFile = regionsFile
	}
	if targets != "" {
		opts.Targets = bcftools.SplitCommaList(targets)
	}
	if targetsFile != "" {
		regs, err := bcftools.LoadRegionsFile(targetsFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "bcftools stats: %v\n", err)
			return 1
		}
		opts.Targets = append(opts.Targets, regs...)
		opts.TargetsFile = targetsFile
	}
	if depthSpec != "" {
		min, max, step, err := bcftools.ParseDepthSpec(depthSpec)
		if err != nil {
			fmt.Fprintln(os.Stderr, err)
			return 2
		}
		opts.DepthMin, opts.DepthMax, opts.DepthStep = min, max, step
	}
	if afBinsSpec != "" {
		bins, err := bcftools.ParseAFBins(afBinsSpec)
		if err != nil {
			fmt.Fprintln(os.Stderr, err)
			return 2
		}
		opts.AFBins = bins
	}

	out, err := openOutFile(outputPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "bcftools stats: %v\n", err)
		return 1
	}
	defer out.Close()

	if _, err := bcftools.StatsFile(input, out, opts); err != nil {
		fmt.Fprintf(os.Stderr, "bcftools stats: %v\n", err)
		return 1
	}
	return 0
}

const queryUsage = `bcftools query - format-string output for VCF/BCF records.

Usage:
  bcftools query [options] <in.vcf[.gz]|in.bcf>

Options:
  -f, --format STRING        Format string. Supported tokens:
                               %CHROM %POS %REF %ALT %QUAL %ID %FILTER
                               %TYPE %TGT %GT %INFO/<TAG> %FMT/<TAG>
                               [%TOKEN ...] for per-sample expansion
                               \n \t literal newline / tab
  -H, --print-header         Emit a header row derived from the format string.
  -l, --list-samples         Print one sample name per line and exit.
  -s, --samples LIST         Restrict per-sample expansion to these names.
  -S, --samples-file PATH    File with sample IDs (one per line).
  -r, --regions LIST         Region list (chr:beg-end[,...]) — uses .tbi/.csi.
  -R, --regions-file PATH    BED-like regions file.
  -t, --targets LIST         Like -r but always a post-filter.
  -T, --targets-file PATH    BED-like targets file (post-filter).
  -i, --include EXPR         Keep records matching expression.
  -e, --exclude EXPR         Drop records matching expression.
  -F, --apply-filters NAMES  Comma list of FILTER names to keep.
  -o, --output PATH          Output file (default stdout).
      --threads N            Accepted; v1 is single-threaded.
  -?, --help                 Show this help.
      --version              Show version.
`

func runQuery(args []string) int {
	fs := flag.NewFlagSet("bcftools query", flag.ContinueOnError)
	fs.SetOutput(io.Discard)

	var (
		format         string
		printHeader    bool
		listSamples    bool
		samples        string
		samplesFile    string
		regions        string
		regionsFile    string
		regionsOverlap int
		targets        string
		targetsFile    string
		targetsOverlap int
		includeExpr    string
		excludeExpr    string
		applyFilters   string
		printFiltered  string
		disableNewline bool
		allowUndefTags bool
		forceSamples   bool
		vcfList        string
		verbosity      int
		outputPath     string
		threads        int
		showHelp       bool
		showVer        bool
	)
	cliflag.StringVar(fs, &format, "f", "format", "", "Format string")
	cliflag.BoolVar(fs, &printHeader, "H", "print-header", false, "Print header row")
	cliflag.BoolVar(fs, &listSamples, "l", "list-samples", false, "List samples and exit")
	cliflag.StringVar(fs, &samples, "s", "samples", "", "Sample list")
	cliflag.StringVar(fs, &samplesFile, "S", "samples-file", "", "Samples file")
	cliflag.StringVar(fs, &regions, "r", "regions", "", "Regions")
	cliflag.StringVar(fs, &regionsFile, "R", "regions-file", "", "Regions file")
	fs.IntVar(&regionsOverlap, "regions-overlap", 1, "")
	cliflag.StringVar(fs, &targets, "t", "targets", "", "Targets (post-filter)")
	cliflag.StringVar(fs, &targetsFile, "T", "targets-file", "", "Targets file")
	fs.IntVar(&targetsOverlap, "targets-overlap", 1, "")
	cliflag.StringVar(fs, &includeExpr, "i", "include", "", "Include expression")
	cliflag.StringVar(fs, &excludeExpr, "e", "exclude", "", "Exclude expression")
	cliflag.StringVar(fs, &applyFilters, "F", "apply-filters", "", "FILTER name list to keep")
	fs.StringVar(&printFiltered, "print-filtered", "", "Print STR for filter-failing records")
	cliflag.BoolVar(fs, &disableNewline, "N", "disable-automatic-newline", false, "Suppress implicit newline")
	cliflag.BoolVar(fs, &allowUndefTags, "u", "allow-undef-tags", false, "Allow undefined tag references")
	fs.BoolVar(&forceSamples, "force-samples", false, "Continue past missing sample names")
	cliflag.StringVar(fs, &vcfList, "v", "vcf-list", "", "File of VCF paths")
	fs.IntVar(&verbosity, "verbosity", 0, "Verbosity (accepted, ignored)")
	cliflag.StringVar(fs, &outputPath, "o", "output", "", "Output path")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "Threads (accepted, ignored)")
	fs.BoolVar(&showHelp, "?", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "version", false, "")

	if err := parseFlags(fs, args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		os.Stderr.WriteString(queryUsage)
		return 2
	}
	if showHelp {
		os.Stdout.WriteString(queryUsage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}
	rest := fs.Args()
	if len(rest) == 0 {
		fmt.Fprintln(os.Stderr, "bcftools query: missing input file")
		os.Stderr.WriteString(queryUsage)
		return 2
	}
	if !listSamples && format == "" {
		fmt.Fprintln(os.Stderr, "bcftools query: -f/--format is required (use -l for sample list)")
		return 2
	}

	opts := bcftools.QueryOptions{
		Format:             format,
		PrintHeader:        printHeader,
		ListSamples:        listSamples,
		IncludeExpr:        includeExpr,
		ExcludeExpr:        excludeExpr,
		ApplyFilters:       bcftools.SplitCommaList(applyFilters),
		SamplesFile:        samplesFile,
		RegionsFile:        regionsFile,
		TargetsFile:        targetsFile,
		PrintFiltered:      printFiltered,
		DisableAutoNewline: disableNewline,
		AllowUndefTags:     allowUndefTags,
		ForceSamples:       forceSamples,
		VCFList:            vcfList,
		RegionsOverlap:     regionsOverlap,
		TargetsOverlap:     targetsOverlap,
		Verbosity:          verbosity,
	}
	if samples != "" {
		opts.Samples = bcftools.SplitCommaList(samples)
	}
	if samplesFile != "" {
		names, err := bcftools.LoadSamplesFile(samplesFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "bcftools query: %v\n", err)
			return 1
		}
		opts.Samples = append(opts.Samples, names...)
	}
	if regions != "" {
		opts.Regions = bcftools.SplitCommaList(regions)
	}
	if regionsFile != "" {
		regs, err := bcftools.LoadRegionsFile(regionsFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "bcftools query: %v\n", err)
			return 1
		}
		opts.Regions = append(opts.Regions, regs...)
	}
	if targets != "" {
		opts.Targets = bcftools.SplitCommaList(targets)
	}
	if targetsFile != "" {
		regs, err := bcftools.LoadRegionsFile(targetsFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "bcftools query: %v\n", err)
			return 1
		}
		opts.Targets = append(opts.Targets, regs...)
	}

	out, err := openOutFile(outputPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "bcftools query: %v\n", err)
		return 1
	}
	defer out.Close()

	if _, err := bcftools.QueryFile(rest[0], out, opts, os.Stderr); err != nil {
		fmt.Fprintf(os.Stderr, "bcftools query: %v\n", err)
		return 1
	}
	return 0
}

const concatUsage = `bcftools concat - concatenate VCF/BCF files.

Usage:
  bcftools concat [options] <in1.vcf[.gz]|in1.bcf> [<in2> ...]

Options:
  -a, --allow-overlaps       Sort-merge across inputs (default is plain concat).
  -D, --remove-duplicates    Drop adjacent duplicate records.
  -f, --file-list PATH       File of input paths (one per line).
  -O, --output-type {v|z|u|b}  Output format. (u/b need a BCF writer.)
  -o, --output PATH          Output file (default stdout).
  -q, --min-PQ INT           Accepted but no-op in v1.
  -l, --ligate               Accepted but no-op in v1 (imputation chunks).
      --threads N            Accepted; v1 is single-threaded.
      --compression-level N  gzip level for -O z output.
  -?, --help                 Show this help.
      --version              Show version.
`

func runConcat(args []string) int {
	fs := flag.NewFlagSet("bcftools concat", flag.ContinueOnError)
	fs.SetOutput(io.Discard)

	var (
		allowOverlaps    bool
		removeDuplicates bool
		fileList         string
		outputType       string
		outputPath       string
		minPQ            int
		ligate           bool
		compressLevel    int
		threads          int
		showHelp         bool
		showVer          bool
	)
	cliflag.BoolVar(fs, &allowOverlaps, "a", "allow-overlaps", false, "Sort-merge across inputs")
	cliflag.BoolVar(fs, &removeDuplicates, "D", "remove-duplicates", false, "Drop adjacent duplicate records")
	cliflag.StringVar(fs, &fileList, "f", "file-list", "", "File of input paths")
	cliflag.StringVar(fs, &outputType, "O", "output-type", "v", "Output type")
	cliflag.StringVar(fs, &outputPath, "o", "output", "", "Output path")
	cliflag.IntVar(fs, &minPQ, "q", "min-PQ", 0, "Minimum PQ (accepted, ignored)")
	cliflag.BoolVar(fs, &ligate, "l", "ligate", false, "Ligate (accepted, ignored)")
	fs.IntVar(&compressLevel, "compression-level", -1, "")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "Threads (accepted, ignored)")
	fs.BoolVar(&showHelp, "?", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "version", false, "")

	if err := parseFlags(fs, args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, concatUsage)
		return 2
	}
	if showHelp {
		fmt.Print(concatUsage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}

	paths := fs.Args()
	if len(paths) == 0 && fileList == "" {
		fmt.Fprintln(os.Stderr, "bcftools concat: missing input files")
		fmt.Fprint(os.Stderr, concatUsage)
		return 2
	}
	format, err := bcftools.ParseOutputFormat(outputType)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}

	out, err := openOutFile(outputPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "bcftools concat: %v\n", err)
		return 1
	}
	defer out.Close()

	opts := bcftools.ConcatOptions{
		OutputFormat:     format,
		AllowOverlaps:    allowOverlaps,
		RemoveDuplicates: removeDuplicates,
		FileList:         fileList,
		MinPQ:            minPQ,
		Ligate:           ligate,
		CompressLevel:    compressLevel,
	}
	if _, err := bcftools.ConcatFiles(paths, out, opts); err != nil {
		fmt.Fprintf(os.Stderr, "bcftools concat: %v\n", err)
		return 1
	}
	return 0
}

const normUsage = `bcftools norm - left-align indels and normalize multiallelics.

Usage:
  bcftools norm [options] <in.vcf[.gz]|in.bcf>

Options:
  -a, --atomize                  Decompose complex variants (e.g. MNVs become consecutive SNVs).
      --atom-overlaps '*'|.      Star allele (default) or missing for overlapping ALTs.
  -c, --check-ref {e|w|x|s}      REF mismatch policy: e=error (default), w=warn, x=exclude, s=set.
  -D, --remove-duplicates        Remove duplicate lines of the same type.
  -d, --rm-dup MODE              Drop duplicates: snps|indels|both|all|none|exact.
  -e, --exclude EXPR             Skip records matching the filter expression.
  -f, --fasta-ref FASTA          Reference FASTA for left-alignment / REF check.
      --force                    Continue past malformed records (experimental).
  -g, --gff-annot FILE           HGVS 3'-rule right-alignment (deferred — see roadmap).
  -i, --include EXPR             Keep only records matching the filter expression.
      --keep-sum TAG[,TAG...]    Keep INFO/<TAG> vector sum constant when splitting.
  -m, --multiallelics MODE       Split (-) or join (+) multiallelics. MODE = {-snps|-indels|-both|-any|+snps|+indels|+both|+any}.
      --multi-overlaps 0|.       Reference (0) or missing (.) fill when splitting [0].
      --no-version               Omit the bcftools provenance lines from the header.
  -N, --no-realign [NUM]         Skip left-alignment (or only when event > NUM bp).
      --do-not-normalize         Alias for --no-realign (legacy).
      --old-rec-tag STR          Annotate modified records with INFO/<STR>=<orig>.
  -o, --output PATH              Output file (default stdout).
  -O, --output-type {v|z|u|b}    Output format (b/u requires BCF writer).
  -r, --regions chr[:beg-end]    Region(s) (post-filter on streaming input).
  -R, --regions-file PATH        BED-like regions file.
      --regions-overlap 0|1|2    Region inclusion rule [1].
  -s, --strict-filter            When merging (-m+), merged site is PASS iff all are PASS.
  -S, --sort METHOD              Sort order: chr_pos (default) | lex.
  -t, --targets chr[:beg-end]    Like -r but always a post-filter.
  -T, --targets-file PATH        BED-like targets file (post-filter).
      --targets-overlap 0|1|2    Target inclusion rule [0].
      --threads N                Accepted; v1 is single-threaded.
  -v, --verbosity INT            Verbosity level (accepted; no extra output yet).
  -w, --site-win INT             Buffer (bp) for re-sorting realigned records [1000].
  -W, --write-index[=FMT]        Automatically index the output files [off].
  -l, --compression-level N      gzip level for -O z output.
  -h, --help                     Show this help.
      --version                  Show version.

Notes:
  -f is overloaded by upstream bcftools for both "fasta-ref" and "apply-filters".
  We follow the same convention: when --fasta-ref is also set we accept --apply-filters
  via the long form for clarity. The short -f always means --fasta-ref to match upstream.
  --gff-annot is rejected verbatim; see docs/PARITY_ROADMAP.md.
`

func runNorm(args []string) int {
	fs := flag.NewFlagSet("bcftools norm", flag.ContinueOnError)
	fs.SetOutput(io.Discard)

	var (
		fastaRef       string
		checkRef       string
		multiallelics  string
		rmDup          string
		removeDups     bool
		atomize        bool
		atomOverlaps   string
		multiOverlaps  string
		doNotNormalize bool
		noRealign      string
		strictFilter   bool
		force          bool
		gffAnnot       string
		includeExpr    string
		excludeExpr    string
		keepSum        string
		oldRecTag      string
		regions        string
		regionsFile    string
		regionsOverlap int
		targets        string
		targetsFile    string
		targetsOverlap int
		sortMethod     string
		siteWin        int
		applyFilters   string
		outputType     string
		outputPath     string
		compressLevel  int
		threads        int
		writeIndex     string
		verbosity      int
		showHelp       bool
		showVersion    bool
	)
	cliflag.StringVar(fs, &fastaRef, "f", "fasta-ref", "", "Reference FASTA")
	cliflag.StringVar(fs, &checkRef, "c", "check-ref", "e", "REF mismatch policy: e|w|x|s")
	cliflag.StringVar(fs, &multiallelics, "m", "multiallelics", "", "Split / join multiallelics")
	cliflag.BoolVar(fs, &removeDups, "D", "remove-duplicates", false, "Remove duplicate lines of the same type")
	cliflag.StringVar(fs, &rmDup, "d", "rm-dup", "none", "Drop duplicate records")
	cliflag.BoolVar(fs, &atomize, "a", "atomize", false, "Atomize complex variants")
	fs.StringVar(&atomOverlaps, "atom-overlaps", "*", "Symbol for overlapping atomized ALTs: '*' or '.'")
	fs.StringVar(&multiOverlaps, "multi-overlaps", "0", "Fill for missing-when-splitting: 0 (REF) or '.'")
	// --do-not-normalize is the legacy alias kept for callers; -N is the
	// upstream --no-realign with optional NUM.
	fs.BoolVar(&doNotNormalize, "do-not-normalize", false, "Alias for --no-realign")
	cliflag.StringVar(fs, &noRealign, "N", "no-realign", "", "Skip left-alignment (or only when event > NUM bp)")
	cliflag.BoolVar(fs, &strictFilter, "s", "strict-filter", false, "Apply -f filters before split")
	fs.BoolVar(&force, "force", false, "Continue past malformed records")
	cliflag.StringVar(fs, &gffAnnot, "g", "gff-annot", "", "HGVS 3'-rule right-alignment (deferred)")
	cliflag.StringVar(fs, &includeExpr, "i", "include", "", "Keep records matching expression")
	cliflag.StringVar(fs, &excludeExpr, "e", "exclude", "", "Drop records matching expression")
	fs.StringVar(&keepSum, "keep-sum", "", "Keep INFO/<TAG> vector sum constant when splitting")
	fs.StringVar(&oldRecTag, "old-rec-tag", "", "Annotate modified records with INFO/<STR>=<orig>")
	cliflag.StringVar(fs, &regions, "r", "regions", "", "Region(s)")
	cliflag.StringVar(fs, &regionsFile, "R", "regions-file", "", "Regions file")
	fs.IntVar(&regionsOverlap, "regions-overlap", 1, "Region inclusion rule (0|1|2)")
	cliflag.StringVar(fs, &targets, "t", "targets", "", "Targets")
	cliflag.StringVar(fs, &targetsFile, "T", "targets-file", "", "Targets file")
	fs.IntVar(&targetsOverlap, "targets-overlap", 0, "Target inclusion rule (0|1|2)")
	cliflag.StringVar(fs, &sortMethod, "S", "sort", "chr_pos", "Sort method: chr_pos|lex")
	cliflag.IntVar(fs, &siteWin, "w", "site-win", 1000, "Buffer (bp) for re-sorting realigned records")
	fs.StringVar(&applyFilters, "apply-filters", "", "Filter list (PASS,...)")
	cliflag.StringVar(fs, &outputType, "O", "output-type", "v", "Output type")
	cliflag.StringVar(fs, &outputPath, "o", "output", "", "Output path")
	cliflag.IntVar(fs, &compressLevel, "l", "compression-level", -1, "gzip level")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "Threads (accepted, ignored)")
	cliflag.StringVar(fs, &writeIndex, "W", "write-index", "", "Automatically index output (csi|tbi)")
	cliflag.IntVar(fs, &verbosity, "v", "verbosity", 0, "Verbosity level")
	fs.BoolVar(&showHelp, "h", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVersion, "version", false, "")
	registerNoVersionIfAbsent(fs)

	// -W/--write-index accepts a bare form; expand to the upstream default.
	args = preprocessOptionalArg(args, "-W", "csi")
	args = preprocessOptionalArg(args, "--write-index", "csi")
	// -N/--no-realign similarly accepts an optional NUM. Bare (or
	// followed by another flag) means "skip realignment for all"; the
	// attached form `-N1000` reaches us as `-N 1000` after
	// normalizeShortFlags so we expand only when the next token is a
	// flag or absent.
	args = preprocessFlagOrBare(args, "-N", "0")
	args = preprocessFlagOrBare(args, "--no-realign", "0")

	if err := parseFlags(fs, args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, normUsage)
		return 2
	}
	if showHelp {
		fmt.Print(normUsage)
		return 0
	}
	if showVersion {
		fmt.Println(version)
		return 0
	}

	rest := fs.Args()
	if len(rest) == 0 {
		fmt.Fprintln(os.Stderr, "bcftools norm: missing input file")
		fmt.Fprint(os.Stderr, normUsage)
		return 2
	}
	input := rest[0]

	if gffAnnot != "" {
		// Rejection-parity: upstream's -g/--gff-annot drives HGVS 3'-rule
		// right-alignment via a GFF3 transcript map. Porting is queued.
		fmt.Fprintln(os.Stderr, "bcftools norm: --gff-annot is to be implemented, please open an issue on github")
		return 1
	}

	checkRefMode, err := bcftools.ParseCheckRefMode(checkRef)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}
	mMode, err := bcftools.ParseMultiallelicMode(multiallelics)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}
	dupMode, err := bcftools.ParseRmDupMode(rmDup)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}
	if removeDups {
		// Upstream -D ≡ -d exact (drop byte-identical duplicates).
		dupMode = bcftools.RmDupExact
	}
	format, err := bcftools.ParseOutputFormat(outputType)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}
	// --atom-overlaps validation
	atomChar := byte('*')
	switch atomOverlaps {
	case "", "*":
		atomChar = '*'
	case ".":
		atomChar = '.'
	default:
		fmt.Fprintf(os.Stderr, "bcftools norm: --atom-overlaps must be '*' or '.', got %q\n", atomOverlaps)
		return 2
	}
	multiChar := byte('0')
	switch multiOverlaps {
	case "", "0":
		multiChar = '0'
	case ".":
		multiChar = '.'
	default:
		fmt.Fprintf(os.Stderr, "bcftools norm: --multi-overlaps must be '0' or '.', got %q\n", multiOverlaps)
		return 2
	}
	if regionsOverlap < 0 || regionsOverlap > 2 {
		fmt.Fprintf(os.Stderr, "bcftools norm: --regions-overlap must be 0, 1 or 2 (got %d)\n", regionsOverlap)
		return 2
	}
	if targetsOverlap < 0 || targetsOverlap > 2 {
		fmt.Fprintf(os.Stderr, "bcftools norm: --targets-overlap must be 0, 1 or 2 (got %d)\n", targetsOverlap)
		return 2
	}
	if writeIndex != "" && writeIndex != "csi" && writeIndex != "tbi" {
		fmt.Fprintf(os.Stderr, "bcftools norm: --write-index must be csi or tbi (got %q)\n", writeIndex)
		return 2
	}

	// Resolve -N / --no-realign.
	var (
		noRealignBool bool
		noRealignMax  int
	)
	if noRealign != "" {
		noRealignBool = true
		// "0" (sentinel for "bare -N") means skip for all events.
		if n, err := strconv.Atoi(noRealign); err == nil {
			noRealignMax = n
		} else {
			fmt.Fprintf(os.Stderr, "bcftools norm: --no-realign value %q must be an integer\n", noRealign)
			return 2
		}
	}

	noVersionFlag := fs.Lookup("no-version")
	noVersion := noVersionFlag != nil && noVersionFlag.Value.String() == "true"

	opts := bcftools.NormOptions{
		FastaRef:        fastaRef,
		CheckRef:        checkRefMode,
		Multiallelics:   mMode,
		RmDup:           dupMode,
		Atomize:         atomize,
		AtomOverlaps:    atomChar,
		MultiOverlaps:   multiChar,
		DoNotNormalize:  doNotNormalize,
		NoRealign:       noRealignBool,
		NoRealignMaxLen: noRealignMax,
		StrictFilter:    strictFilter,
		Force:           force,
		GffAnnot:        gffAnnot,
		IncludeExpr:     includeExpr,
		ExcludeExpr:     excludeExpr,
		KeepSum:         bcftools.SplitCommaList(keepSum),
		OldRecTag:       oldRecTag,
		RegionsOverlap:  regionsOverlap,
		TargetsOverlap:  targetsOverlap,
		Sort:            sortMethod,
		SiteWin:         siteWin,
		NoVersion:       noVersion,
		PGCommand:       strings.Join(os.Args, " "),
		WriteIndex:      writeIndex,
		Verbosity:       verbosity,
		ApplyFilters:    bcftools.SplitCommaList(applyFilters),
		OutputFormat:    format,
		CompressLevel:   compressLevel,
	}
	if regions != "" {
		opts.Regions = bcftools.SplitCommaList(regions)
	}
	if regionsFile != "" {
		regs, err := bcftools.LoadRegionsFile(regionsFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "bcftools norm: %v\n", err)
			return 1
		}
		opts.Regions = append(opts.Regions, regs...)
		opts.RegionsFile = regionsFile
	}
	if targets != "" {
		opts.Targets = bcftools.SplitCommaList(targets)
	}
	if targetsFile != "" {
		regs, err := bcftools.LoadRegionsFile(targetsFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "bcftools norm: %v\n", err)
			return 1
		}
		opts.Targets = append(opts.Targets, regs...)
		opts.TargetsFile = targetsFile
	}

	out, err := openOutFile(outputPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "bcftools norm: %v\n", err)
		return 1
	}
	defer out.Close()

	if _, err := bcftools.NormFile(input, out, opts, os.Stderr); err != nil {
		fmt.Fprintf(os.Stderr, "bcftools norm: %v\n", err)
		return 1
	}
	return 0
}

const callUsage = `bcftools call - variant calling from per-position genotype likelihoods.

Usage:
  bcftools call [options] <in.vcf[.gz]|in.bcf>

The input is the BCF / VCF produced by ` + "`samtools mpileup -g`" + ` (or any
other tool that emits FORMAT/PL per position per sample).

Calling model:
  -c, --consensus-caller         Old (Li 2011) consensus caller.
  -m, --multiallelic-caller      Multi-allelic caller (v1: biallelic-only;
                                 falls back to consensus on multi-allelic
                                 sites — see docs/PARITY_ROADMAP.md).

Filters:
  -A, --keep-alts                Emit every declared ALT even without support.
  -v, --variants-only            Drop all-reference sites.
  -P, --prior FLOAT              Mutation rate prior (default 1.1e-3).
  -p, --pval-threshold FLOAT     Variant-posterior threshold (default 0.5).
      --ploidy {1|2|GRCh37|GRCh38}  Ploidy spec (default 2). GRCh37/38 use the
                                 per-region, per-sex tables from upstream; the
                                 default sex is F.
  -X, --chromosome-X             Legacy alias for --ploidy 1.
  -g, --gvcf INT[,INT...]        Group non-variant sites into gVCF blocks by
                                 minimum per-sample DP (one block per bin).
                                 Rejects --variants-only.

I/O:
  -O, --output-type {v|z|u|b}    Output format (b/u require BCF writer).
  -o, --output PATH              Output file (default stdout).
  -r, --regions chr[:beg-end[,..]] Region(s) (post-filter in v1).
  -R, --regions-file PATH        BED-like regions file.
  -t, --targets chr[:beg-end[,..]] Like -r but always a post-filter.
  -T, --targets-file PATH        BED-like targets file (post-filter).
  -s, --samples LIST             Restrict to these samples.
  -S, --samples-file PATH        File of sample IDs (one per line).
      --threads N                Accepted; v1 is single-threaded.
  -?, --help                     Show this help.
      --version                  Show version.
`

// collapseRepeatedShortFlag finds runs of `-X` (consecutive occurrences
// of the same short flag) and rewrites them as `--longName=N` where N is
// the run length. Mirrors upstream's "counting" short-flag convention
// (e.g. `view -AA` ≡ 2, `view -A` ≡ 1). Long-flag occurrences are not
// rewritten. Only operates on the same flag character — `-AB` is left
// alone.
func collapseRepeatedShortFlag(args []string, shortFlag, longName string) []string {
	if len(shortFlag) != 2 || shortFlag[0] != '-' {
		return args
	}
	out := make([]string, 0, len(args))
	i := 0
	for i < len(args) {
		a := args[i]
		// Match `-AAA` style: leading '-', followed by 1+ repetitions of
		// the flag char. Avoid touching long flags or compound `-aB`.
		if len(a) >= 2 && a[0] == '-' && a[1] != '-' && a[1] == shortFlag[1] {
			n := 0
			for n < len(a)-1 && a[1+n] == shortFlag[1] {
				n++
			}
			if n == len(a)-1 {
				out = append(out, "--"+longName+"="+strconv.Itoa(n))
				i++
				continue
			}
		}
		out = append(out, a)
		i++
	}
	return out
}

// preprocessFlagOrBare expands `flagName` to `flagName=defaultVal`
// whenever it is bare (last arg) OR immediately followed by another
// flag (`-X` or `--Y`) OR a non-integer non-flag token. Only the
// attached form (`-N1000`, normalized to `-N 1000` upstream of this
// pass) leaves an integer in the next slot — that is the optional NUM.
// Mirrors upstream's "-N1000, no space" rule.
func preprocessFlagOrBare(args []string, flagName, defaultVal string) []string {
	out := make([]string, 0, len(args))
	for i := 0; i < len(args); i++ {
		a := args[i]
		if a == flagName {
			if i+1 >= len(args) {
				out = append(out, flagName+"="+defaultVal)
				continue
			}
			next := args[i+1]
			isFlag := len(next) > 0 && next[0] == '-' && next != "-"
			_, atoiErr := strconv.Atoi(next)
			isInt := atoiErr == nil
			if isFlag || !isInt {
				out = append(out, flagName+"="+defaultVal)
				continue
			}
		}
		out = append(out, a)
	}
	return out
}

// preprocessOptionalArg expands `--flag` (or `-W`) without a value into
// `--flag=defaultVal` so the Go flag package doesn't treat the next
// positional as the flag's value. `--flag=X` and `--flag X` are
// untouched.
func preprocessOptionalArg(args []string, flagName, defaultVal string) []string {
	out := make([]string, 0, len(args))
	for i := 0; i < len(args); i++ {
		a := args[i]
		if a == flagName {
			// Peek next arg: if it looks like a value (no leading
			// '-'), assume the user meant `--flag VALUE` and pass
			// through; otherwise expand to `--flag=defaultVal`.
			if i+1 < len(args) {
				next := args[i+1]
				if next == "" || next[0] != '-' {
					// Check if next looks like a positional file (e.g. has
					// '/' or ends in .vcf / .bcf / .gz) — heuristic.
					if looksLikePositional(next) {
						out = append(out, flagName+"="+defaultVal)
						continue
					}
				}
			} else {
				// Trailing bare flag.
				out = append(out, flagName+"="+defaultVal)
				continue
			}
		}
		out = append(out, a)
	}
	return out
}

// looksLikePositional reports whether s looks more like an input file
// path than an option value (heuristic: contains '/', '.', or ends in
// common suffixes).
func looksLikePositional(s string) bool {
	if strings.ContainsAny(s, "/") {
		return true
	}
	for _, suf := range []string{".vcf", ".bcf", ".gz", ".bgz", ".tbi", ".csi"} {
		if strings.HasSuffix(s, suf) {
			return true
		}
	}
	return false
}

// priorANFromSpec returns the first field of "AN,AC".
func priorANFromSpec(s string) string {
	if s == "" {
		return ""
	}
	if i := strings.IndexByte(s, ','); i >= 0 {
		return strings.TrimSpace(s[:i])
	}
	return strings.TrimSpace(s)
}

// priorACFromSpec returns the second field of "AN,AC".
func priorACFromSpec(s string) string {
	if i := strings.IndexByte(s, ','); i >= 0 {
		return strings.TrimSpace(s[i+1:])
	}
	return ""
}

func runCall(args []string) int {
	fs := flag.NewFlagSet("bcftools call", flag.ContinueOnError)
	fs.SetOutput(io.Discard)

	var (
		consensus      bool
		multiallelic   bool
		keepAlts       bool
		variantsOnly   bool
		prior          float64
		pval           float64
		ploidy         string
		ploidyFile     string
		haploidX       bool
		outputType     string
		outputPath     string
		regions        string
		regionsFile    string
		regionsOverlap int
		targets        string
		targetsFile    string
		samples        string
		samplesFile    string
		threads        int
		gvcf           string
		constrain      string
		groupSamples   string
		groupSmplTag   string
		noVersion      bool
		keepUnseen     bool
		keepMaskedRef  bool
		skipVariants   string
		annotate       string
		writeIndex     string
		insertMissed   bool
		priorFreqs     string
		novelRate      string
		verbosity      int
		showHelp       bool
		showVer        bool
	)
	cliflag.BoolVar(fs, &consensus, "c", "consensus-caller", false, "Use consensus caller")
	cliflag.BoolVar(fs, &multiallelic, "m", "multiallelic-caller", false, "Use multi-allelic caller")
	cliflag.BoolVar(fs, &keepAlts, "A", "keep-alts", false, "Keep all ALT alleles")
	cliflag.BoolVar(fs, &variantsOnly, "v", "variants-only", false, "Drop all-reference sites")
	cliflag.Float64Var(fs, &prior, "P", "prior", 1.1e-3, "Mutation rate prior")
	cliflag.Float64Var(fs, &pval, "p", "pval-threshold", 0.5, "P-value threshold")
	fs.StringVar(&ploidy, "ploidy", "2", "Ploidy spec")
	fs.StringVar(&ploidyFile, "ploidy-file", "", "Ploidy file (CHROM,FROM,TO,SEX,PLOIDY)")
	cliflag.BoolVar(fs, &haploidX, "X", "chromosome-X", false, "Treat samples as haploid")
	cliflag.StringVar(fs, &outputType, "O", "output-type", "v", "Output type")
	cliflag.StringVar(fs, &outputPath, "o", "output", "", "Output path")
	cliflag.StringVar(fs, &regions, "r", "regions", "", "Region(s)")
	cliflag.StringVar(fs, &regionsFile, "R", "regions-file", "", "Regions file")
	cliflag.IntVar(fs, &regionsOverlap, "", "regions-overlap", 1, "Region overlap semantic (0|1|2)")
	cliflag.StringVar(fs, &targets, "t", "targets", "", "Targets (post-filter)")
	cliflag.StringVar(fs, &targetsFile, "T", "targets-file", "", "Targets file")
	cliflag.StringVar(fs, &samples, "s", "samples", "", "Samples list")
	cliflag.StringVar(fs, &samplesFile, "S", "samples-file", "", "Samples file")
	cliflag.IntVar(fs, &threads, "@", "threads", 0, "Threads (accepted, ignored)")
	cliflag.StringVar(fs, &gvcf, "g", "gvcf", "", "Group non-variant sites into gVCF blocks by minimum per-sample DP (INT[,INT...])")
	cliflag.StringVar(fs, &constrain, "C", "constrain", "", "Constrain calling to one of: alleles, trio (requires -T)")
	cliflag.StringVar(fs, &groupSamples, "G", "group-samples", "", "Sample-group file for per-pool calling (`-` for per-sample groups)")
	fs.StringVar(&groupSmplTag, "group-samples-tag", "", "FORMAT tag for -G (default auto-detect QS or AD)")
	fs.BoolVar(&noVersion, "no-version", false, "Do not append version/command to header")
	cliflag.BoolVar(fs, &keepUnseen, "*", "keep-unseen-allele", false, "Keep the <*> / <NON_REF> allele")
	cliflag.BoolVar(fs, &keepMaskedRef, "M", "keep-masked-ref", false, "Keep sites whose REF base is N")
	cliflag.StringVar(fs, &skipVariants, "V", "skip-variants", "", "Skip records of type 'indels' or 'snps'")
	cliflag.StringVar(fs, &annotate, "a", "annotate", "", "Optional output tags (comma-separated; '?' to list)")
	fs.StringVar(&writeIndex, "write-index", "", "Index the output (csi|tbi)")
	fs.BoolVar(&insertMissed, "i", false, "Insert records for sites in -T not seen by mpileup")
	fs.BoolVar(&insertMissed, "insert-missed", false, "Insert records for sites in -T not seen by mpileup")
	cliflag.StringVar(fs, &priorFreqs, "F", "prior-freqs", "", "AN,AC INFO tags providing prior allele frequencies")
	cliflag.StringVar(fs, &novelRate, "n", "novel-rate", "", "Novel-rate spec for constrained trio calling")
	fs.IntVar(&verbosity, "verbosity", 0, "Verbosity level (accepted, ignored)")
	fs.BoolVar(&showHelp, "?", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "version", false, "")

	// Pre-process bare `--write-index` / `-W` (no value) into the
	// explicit `--write-index=csi` form so Go's flag package
	// doesn't swallow the next positional as its argument.
	args = preprocessOptionalArg(args, "--write-index", "csi")
	args = preprocessOptionalArg(args, "-W", "csi")
	if err := parseFlags(fs, args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, callUsage)
		return 2
	}
	if showHelp {
		fmt.Print(callUsage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}
	rest := fs.Args()
	if len(rest) == 0 {
		fmt.Fprintln(os.Stderr, "bcftools call: missing input file")
		fmt.Fprint(os.Stderr, callUsage)
		return 2
	}

	model := bcftools.CallModelNone
	switch {
	case consensus && multiallelic:
		fmt.Fprintln(os.Stderr, "bcftools call: choose either -c/--consensus-caller or -m/--multiallelic-caller, not both")
		return 2
	case consensus:
		model = bcftools.CallModelConsensus
	case multiallelic:
		model = bcftools.CallModelMultiallelic
	default:
		fmt.Fprintln(os.Stderr, "bcftools call: a caller must be selected (-c or -m)")
		return 2
	}

	ploidySpec, ploidyText, err := bcftools.ParsePloidySpec(ploidy)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}
	if haploidX {
		ploidySpec = bcftools.PloidyHaploid
		ploidyText = "1"
	}

	format, err := bcftools.ParseOutputFormat(outputType)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}

	opts := bcftools.CallOptions{
		Model:           model,
		KeepAlts:        keepAlts,
		VariantsOnly:    variantsOnly,
		Prior:           prior,
		PvalThreshold:   pval,
		Ploidy:          ploidySpec,
		PloidySpec:      ploidyText,
		PloidyFile:      ploidyFile,
		OutputFormat:    format,
		GVCFSpec:        gvcf,
		NoVersion:       noVersion,
		PGCommand:       strings.Join(os.Args, " "),
		KeepUnseen:      keepUnseen,
		KeepMaskedRef:   keepMaskedRef,
		SkipVariants:    skipVariants,
		RegionsOverlap:  regionsOverlap,
		GroupSamplesTag: groupSmplTag,
		Annotate:        annotate,
		WriteIndex:      writeIndex,
		InsertMissed:    insertMissed,
		PriorAN:         priorANFromSpec(priorFreqs),
		PriorAC:         priorACFromSpec(priorFreqs),
		NovelRate:       novelRate,
		Verbosity:       verbosity,
	}
	switch constrain {
	case "":
		// no constraint
	case "alleles":
		opts.Constrain = bcftools.CallConstrainAlleles
		opts.ConstrainSites = targetsFile
		if opts.ConstrainSites == "" {
			fmt.Fprintln(os.Stderr, "bcftools call: -C alleles requires -T sites_file")
			return 2
		}
	case "trio":
		opts.Constrain = bcftools.CallConstrainTrio
	default:
		fmt.Fprintf(os.Stderr, "bcftools call: unknown -C value %q (expect: alleles, trio)\n", constrain)
		return 2
	}
	if groupSamples != "" {
		opts.GroupSamplesFile = groupSamples
	}
	if regions != "" {
		opts.Regions = bcftools.SplitCommaList(regions)
	}
	if regionsFile != "" {
		regs, err := bcftools.LoadRegionsFile(regionsFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "bcftools call: %v\n", err)
			return 1
		}
		opts.Regions = append(opts.Regions, regs...)
		opts.RegionsFile = regionsFile
	}
	if targets != "" {
		opts.Targets = bcftools.SplitCommaList(targets)
	}
	if targetsFile != "" && constrain != "alleles" {
		regs, err := bcftools.LoadRegionsFile(targetsFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "bcftools call: %v\n", err)
			return 1
		}
		opts.Targets = append(opts.Targets, regs...)
		opts.TargetsFile = targetsFile
	}
	if samples != "" {
		opts.Samples = bcftools.SplitCommaList(samples)
	}
	if samplesFile != "" {
		names, err := bcftools.LoadSamplesFile(samplesFile)
		if err != nil {
			fmt.Fprintf(os.Stderr, "bcftools call: %v\n", err)
			return 1
		}
		opts.Samples = append(opts.Samples, names...)
		opts.SamplesFile = samplesFile
	}

	out, err := openOutFile(outputPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "bcftools call: %v\n", err)
		return 1
	}

	if _, err := bcftools.CallFile(rest[0], out, opts, os.Stderr); err != nil {
		_ = out.Close()
		fmt.Fprintf(os.Stderr, "bcftools call: %v\n", err)
		return 1
	}
	if err := out.Close(); err != nil {
		fmt.Fprintf(os.Stderr, "bcftools call: %v\n", err)
		return 1
	}
	// -W/--write-index: after the output is fully written, index
	// it via the in-tree bcftools.BuildIndex machinery. The
	// optional FMT suffix selects "csi" (default) or "tbi".
	if writeIndex != "" && outputPath != "" {
		fmt2 := bcftools.IndexCSI
		switch strings.ToLower(writeIndex) {
		case "", "csi":
			fmt2 = bcftools.IndexCSI
		case "tbi":
			fmt2 = bcftools.IndexTBI
		default:
			fmt.Fprintf(os.Stderr, "bcftools call: unknown --write-index format %q (expect csi|tbi)\n", writeIndex)
			return 2
		}
		if _, err := bcftools.BuildIndex(outputPath, bcftools.IndexOptions{Format: fmt2}); err != nil {
			fmt.Fprintf(os.Stderr, "bcftools call: --write-index: %v\n", err)
			return 1
		}
	}
	return 0
}
