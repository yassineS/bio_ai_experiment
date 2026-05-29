// CLI runners for `bcftools filter` and `bcftools consensus`. The
// shape follows the runners in main.go and the earlier subcmds_*.go
// files: parse flags, validate, dispatch to the library. The
// project parity rule (docs/PARITY_ROADMAP.md "Definition of 1:1")
// is enforced here — every documented upstream flag must parse
// cleanly, and flags whose underlying behaviour is deferred either
// no-op (when their default is no-op) or hard-reject with a roadmap
// pointer when explicitly set.
package main

import (
	"flag"
	"fmt"
	"io"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/cliflag"
	"github.com/yassineS/bio_ai_experiment/tools/bcftools/pkg/bcftools"
)

const filterUsage = `bcftools filter - soft-filter records by expression.

Usage:
  bcftools filter [options] <in.vcf[.gz]|in.bcf>

Soft-filter mode differs from "view -i/-e": records that fail the
expression are KEPT in the output but tagged in the FILTER column with
the -s/--soft-filter NAME (rather than being dropped).

Options:
  -e, --exclude EXPR             Soft-filter sites where EXPR is true.
  -g, --SnpGap INT[:TYPE]        Filter SNPs within INT bp of an indel.
                                 The :TYPE qualifier (indel|mnp|bnd|other|overlap)
                                 is accepted; v1 treats every neighbour as "indel".
  -G, --IndelGap INT             Filter clusters of indels separated by INT or fewer bp.
  -i, --include EXPR             Soft-filter sites where EXPR is false.
      --mask [^]REGION           Soft-filter the named region. Accepted; v1 not
                                 implemented (see docs/PARITY_ROADMAP.md#bcftools).
  -M, --mask-file [^]FILE        Soft-filter regions in a BED-like file. Accepted; v1
                                 not implemented (see docs/PARITY_ROADMAP.md#bcftools).
      --mask-overlap 0|1|2       Mask-overlap mode. Accepted; v1 ignores.
  -m, --mode +|x|+x              "+": append to existing FILTER (default: replace);
                                 "x": reset FILTER to PASS on passing sites.
      --no-version               Do not append a ##bcftools_filterCommand= header line.
  -o, --output FILE              Write output to FILE (default stdout).
  -O, --output-type {v|z|u|b}    v=VCF (default), z=VCF.gz, u=uncompressed BCF, b=BCF.
  -r, --regions LIST             Region post-filter (uses .tbi/.csi when available
                                 in v1 this is always a post-filter).
  -R, --regions-file FILE        BED-like regions file.
      --regions-overlap 0|1|2    Accepted; v1 always uses POS-in-region.
  -s, --soft-filter STRING       FILTER tag for failing records (or "+" for an
                                 auto-picked unique "FilterN" name).
  -S, --set-GTs .|0              Replace GTs of failing-record samples with missing (.)
                                 or homozygous reference (0).
  -t, --targets LIST             Like -r but always a post-filter.
  -T, --targets-file FILE        BED-like targets file (post-filter).
      --targets-overlap 0|1|2    Accepted; v1 always uses POS-in-region.
      --threads N                Accepted; v1 is single-threaded.
  -v, --verbosity INT            Accepted; v1 ignores.
  -W, --write-index[=FMT]        Accepted; v1 does not auto-index outputs.
  -?, --help                     Show this help.
      --version                  Show version.

Notes:
  - --mask / --mask-file replace failing records' FILTER with the soft
    filter name from the per-region annotation. v1 doesn't load the
    mask file, so passing -M/--mask-file errors out with a roadmap
    pointer.
  - -g/--SnpGap accepts the upstream :TYPE qualifier (indel|mnp|bnd|other|overlap)
    but always treats every nearby indel as an "indel" type in v1.
`

func runFilter(args []string) int {
	fs := flag.NewFlagSet("bcftools filter", flag.ContinueOnError)
	fs.SetOutput(io.Discard)

	var (
		excludeExpr    string
		snpGapSpec     string
		indelGap       int
		includeExpr    string
		maskRegion     string
		maskFile       string
		maskOverlap    int
		modeFlag       string
		noVersion      bool
		outputPath     string
		outputType     string
		regions        string
		regionsFile    string
		regionsOverlap int
		softFilter     string
		setGTs         string
		targets        string
		targetsFile    string
		targetsOverlap int
		threads        int
		verbosity      int
		writeIndex     string
		showHelp       bool
		showVer        bool
	)

	cliflag.StringVar(fs, &excludeExpr, "e", "exclude", "", "Exclude expression")
	cliflag.StringVar(fs, &snpGapSpec, "g", "SnpGap", "", "SNP gap (INT[:TYPE])")
	cliflag.IntVar(fs, &indelGap, "G", "IndelGap", 0, "Indel gap (INT)")
	cliflag.StringVar(fs, &includeExpr, "i", "include", "", "Include expression")
	fs.StringVar(&maskRegion, "mask", "", "Mask region (accepted; v1 not implemented)")
	cliflag.StringVar(fs, &maskFile, "M", "mask-file", "", "Mask file (accepted; v1 not implemented)")
	fs.IntVar(&maskOverlap, "mask-overlap", 1, "")
	cliflag.StringVar(fs, &modeFlag, "m", "mode", "", "Mode (+|x|+x)")
	fs.BoolVar(&noVersion, "no-version", false, "Do not append ##bcftools_filterCommand")
	cliflag.StringVar(fs, &outputPath, "o", "output", "", "Output file")
	cliflag.StringVar(fs, &outputType, "O", "output-type", "v", "Output type")
	cliflag.StringVar(fs, &regions, "r", "regions", "", "Regions")
	cliflag.StringVar(fs, &regionsFile, "R", "regions-file", "", "Regions file")
	fs.IntVar(&regionsOverlap, "regions-overlap", 1, "")
	cliflag.StringVar(fs, &softFilter, "s", "soft-filter", "", "Soft-filter name")
	cliflag.StringVar(fs, &setGTs, "S", "set-GTs", "", "Replace failing GTs (.|0)")
	cliflag.StringVar(fs, &targets, "t", "targets", "", "Targets (post-filter)")
	cliflag.StringVar(fs, &targetsFile, "T", "targets-file", "", "Targets file")
	fs.IntVar(&targetsOverlap, "targets-overlap", 0, "")
	fs.IntVar(&threads, "threads", 0, "Threads (accepted, ignored)")
	cliflag.IntVar(fs, &verbosity, "v", "verbosity", 0, "Verbosity (accepted, ignored)")
	cliflag.StringVar(fs, &writeIndex, "W", "write-index", "", "Auto-index outputs (accepted, ignored)")
	fs.BoolVar(&showHelp, "?", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "version", false, "")
	// Reference once to silence "declared but not used" / "key part of
	// the upstream surface but no v1 behaviour" notes.
	_ = maskOverlap
	_ = regionsOverlap
	_ = targetsOverlap
	_ = threads
	_ = verbosity
	_ = writeIndex

	if err := parseFlags(fs, args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, filterUsage)
		return 2
	}
	if showHelp {
		fmt.Print(filterUsage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}

	if deferred := checkFilterDeferred(checkFilterDeferredInputs{
		maskRegion: maskRegion,
		maskFile:   maskFile,
	}); deferred != "" {
		fmt.Fprintf(os.Stderr, "bcftools filter: %s is not implemented in v1; tracked in docs/PARITY_ROADMAP.md#bcftools\n", deferred)
		return 2
	}

	mode, err := bcftools.ParseFilterMode(modeFlag)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}
	setGTsMode, err := bcftools.ParseSetGTsMode(setGTs)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}
	format, err := bcftools.ParseOutputFormat(outputType)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}

	snpGap, err := parseSnpGap(snpGapSpec)
	if err != nil {
		fmt.Fprintln(os.Stderr, err)
		return 2
	}

	rest := fs.Args()
	if len(rest) == 0 {
		fmt.Fprintln(os.Stderr, "bcftools filter: missing input file")
		fmt.Fprint(os.Stderr, filterUsage)
		return 2
	}

	if setGTsMode != bcftools.SetGTsOff && softFilter == "" {
		fmt.Fprintln(os.Stderr, "bcftools filter: -S/--set-GTs requires -s/--soft-filter")
		return 2
	}
	if mode&bcftools.FilterModeAdd != 0 && softFilter == "" {
		// -m + on its own with no -s makes no sense.
		fmt.Fprintln(os.Stderr, "bcftools filter: -m + requires -s/--soft-filter")
		return 2
	}

	opts := bcftools.VCFFilterOptions{
		OutputFormat: format,
		IncludeExpr:  includeExpr,
		ExcludeExpr:  excludeExpr,
		SoftFilter:   softFilter,
		Mode:         mode,
		SetGTs:       setGTsMode,
		SnpGap:       snpGap,
		IndelGap:     indelGap,
		RegionsFile:  regionsFile,
		TargetsFile:  targetsFile,
		NoVersion:    noVersion,
	}
	if regions != "" {
		opts.Regions = bcftools.SplitCommaList(regions)
	}
	if targets != "" {
		opts.Targets = bcftools.SplitCommaList(targets)
	}

	out, err := openOutFile(outputPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "bcftools filter: %v\n", err)
		return 1
	}
	defer out.Close()

	if _, err := bcftools.VCFFilterFile(rest[0], out, opts); err != nil {
		fmt.Fprintf(os.Stderr, "bcftools filter: %v\n", err)
		return 1
	}
	return 0
}

type checkFilterDeferredInputs struct {
	maskRegion string
	maskFile   string
}

func checkFilterDeferred(in checkFilterDeferredInputs) string {
	switch {
	case in.maskRegion != "":
		return "--mask"
	case in.maskFile != "":
		return "-M/--mask-file"
	}
	return ""
}

// parseSnpGap accepts the upstream "INT" or "INT:TYPE[,TYPE,...]" form
// and returns the integer part. The :TYPE qualifier is parsed-and-accepted
// but always treated as "indel" in v1 (the upstream default).
func parseSnpGap(s string) (int, error) {
	if s == "" {
		return 0, nil
	}
	// Strip any ":TYPE" qualifier.
	for i := 0; i < len(s); i++ {
		if s[i] == ':' {
			s = s[:i]
			break
		}
	}
	n := 0
	for i := 0; i < len(s); i++ {
		if s[i] < '0' || s[i] > '9' {
			return 0, fmt.Errorf("bcftools filter: bad --SnpGap %q", s)
		}
		n = n*10 + int(s[i]-'0')
	}
	return n, nil
}

const consensusUsage = `bcftools consensus - apply variants from a VCF to a reference FASTA.

Usage:
  bcftools consensus [options] <in.vcf[.gz]|in.bcf>

Options:
  -c, --chain FILE               Write a chain file for liftover. Accepted;
                                 v1 not implemented (see docs/PARITY_ROADMAP.md#bcftools).
  -a, --absent CHAR              Replace positions absent from VCF with CHAR.
  -e, --exclude EXPR             Exclude sites where EXPR is true.
  -f, --fasta-ref FILE           Reference sequence in FASTA format (required).
  -H, --haplotype WHICH          Choose which allele to use from FORMAT/GT:
                                   N: 1-based allele index, regardless of phasing.
                                   R: REF allele in het genotypes.
                                   A: ALT allele in het genotypes.
                                   I: IUPAC code (see also -I).
                                   LR,LA: longer allele, REF/ALT on ties.
                                   SR,SA: shorter allele, REF/ALT on ties.
                                   NpIu: phased index / unphased IUPAC (accepted;
                                         v1 not implemented).
  -i, --include EXPR             Include sites where EXPR is true.
  -I, --iupac-codes              Encode het genotypes as IUPAC ambiguity codes.
      --mark-del CHAR            Replace deleted bases with CHAR (not removed).
      --mark-ins uc|lc|CHAR      Highlight inserted bases.
      --mark-snv uc|lc|CHAR      Highlight substituted bases.
  -m, --mask FILE                BED of regions to replace per --mask-with.
      --mask-with CHAR|uc|lc     Replacement spec for -m mask (default 'N').
  -M, --missing CHAR             Emit CHAR for missing genotypes "./." (default skip).
  -o, --output FILE              Write output to FILE (default stdout).
  -p, --prefix STRING            Prefix to add to each output sequence name.
      --regions-overlap 0|1|2    Accepted; v1 ignores (no index-backed regions).
  -s, --samples LIST             Restrict to these samples (only the first is honoured
                                 in v1; multi-sample apply is a follow-up).
  -S, --samples-file FILE        File of samples to include.
  -v, --verbosity INT            Accepted; v1 ignores.
  -?, --help                     Show this help.
      --version                  Show version.

Notes:
  - v1 applies SNPs, simple insertions (REF=A, ALT=AC), and simple
    deletions (REF=AC, ALT=A). Complex MNPs and structural variants
    are tracked in docs/PARITY_ROADMAP.md#bcftools.
  - Overlapping variants on the same contig: first wins (left-to-right).
`

func runConsensus(args []string) int {
	fs := flag.NewFlagSet("bcftools consensus", flag.ContinueOnError)
	fs.SetOutput(io.Discard)

	var (
		chainFile      string
		absent         string
		excludeExpr    string
		fastaRef       string
		haplotype      string
		includeExpr    string
		iupacCodes     bool
		markDel        string
		markIns        string
		markSnv        string
		maskFile       string
		maskWith       string
		missing        string
		outputPath     string
		prefix         string
		regionsOverlap int
		samples        string
		samplesFile    string
		verbosity      int
		showHelp       bool
		showVer        bool
	)

	cliflag.StringVar(fs, &chainFile, "c", "chain", "", "Liftover chain file (accepted; v1 not implemented)")
	cliflag.StringVar(fs, &absent, "a", "absent", "", "Replace absent positions with CHAR")
	cliflag.StringVar(fs, &excludeExpr, "e", "exclude", "", "Exclude expression")
	cliflag.StringVar(fs, &fastaRef, "f", "fasta-ref", "", "Reference FASTA")
	cliflag.StringVar(fs, &haplotype, "H", "haplotype", "", "Haplotype selector")
	cliflag.StringVar(fs, &includeExpr, "i", "include", "", "Include expression")
	cliflag.BoolVar(fs, &iupacCodes, "I", "iupac-codes", false, "IUPAC codes for het GTs")
	fs.StringVar(&markDel, "mark-del", "", "Replace deleted bases with CHAR")
	fs.StringVar(&markIns, "mark-ins", "", "Highlight insertions (uc|lc|CHAR)")
	fs.StringVar(&markSnv, "mark-snv", "", "Highlight SNVs (uc|lc|CHAR)")
	cliflag.StringVar(fs, &maskFile, "m", "mask", "", "BED of regions to replace")
	fs.StringVar(&maskWith, "mask-with", "", "Replacement spec for -m mask (default N)")
	cliflag.StringVar(fs, &missing, "M", "missing", "", "Emit CHAR for missing GTs")
	cliflag.StringVar(fs, &outputPath, "o", "output", "", "Output file")
	cliflag.StringVar(fs, &prefix, "p", "prefix", "", "Sequence name prefix")
	fs.IntVar(&regionsOverlap, "regions-overlap", 1, "")
	cliflag.StringVar(fs, &samples, "s", "samples", "", "Sample list")
	// Upstream consensus.c:1245 also accepts `--sample` (singular)
	// as an alias for `-s/--samples`.
	fs.StringVar(&samples, "sample", "", "Sample (alias for --samples)")
	cliflag.StringVar(fs, &samplesFile, "S", "samples-file", "", "Samples file")
	cliflag.IntVar(fs, &verbosity, "v", "verbosity", 0, "Verbosity (accepted, ignored)")
	fs.BoolVar(&showHelp, "?", false, "")
	fs.BoolVar(&showHelp, "help", false, "")
	fs.BoolVar(&showVer, "version", false, "")
	_ = regionsOverlap
	_ = verbosity

	if err := parseFlags(fs, args); err != nil {
		fmt.Fprintln(os.Stderr, err)
		fmt.Fprint(os.Stderr, consensusUsage)
		return 2
	}
	if showHelp {
		fmt.Print(consensusUsage)
		return 0
	}
	if showVer {
		fmt.Println(version)
		return 0
	}

	if deferred := checkConsensusDeferred(checkConsensusDeferredInputs{
		chainFile: chainFile,
		haplotype: haplotype,
	}); deferred != "" {
		fmt.Fprintf(os.Stderr, "bcftools consensus: %s is not implemented in v1; tracked in docs/PARITY_ROADMAP.md#bcftools\n", deferred)
		return 2
	}

	rest := fs.Args()
	if len(rest) == 0 {
		fmt.Fprintln(os.Stderr, "bcftools consensus: missing input file")
		fmt.Fprint(os.Stderr, consensusUsage)
		return 2
	}
	if fastaRef == "" {
		fmt.Fprintln(os.Stderr, "bcftools consensus: -f/--fasta-ref is required")
		return 2
	}

	hapSel, hapIdx, err := bcftools.ParseHaplotypeSelector(haplotype)
	if err != nil {
		fmt.Fprintln(os.Stderr, "bcftools consensus:", err)
		return 2
	}

	markDelSpec, err := parseSingleCharMark(markDel, "--mark-del")
	if err != nil {
		fmt.Fprintln(os.Stderr, "bcftools consensus:", err)
		return 2
	}
	markInsSpec, err := bcftools.ParseMarkSpec(markIns)
	if err != nil {
		fmt.Fprintln(os.Stderr, "bcftools consensus: --mark-ins:", err)
		return 2
	}
	markSnvSpec, err := bcftools.ParseMarkSpec(markSnv)
	if err != nil {
		fmt.Fprintln(os.Stderr, "bcftools consensus: --mark-snv:", err)
		return 2
	}
	maskWithSpec, err := bcftools.ParseMarkSpec(maskWith)
	if err != nil {
		fmt.Fprintln(os.Stderr, "bcftools consensus: --mask-with:", err)
		return 2
	}

	// Resolve sample list. -s "-" is upstream's "ignore samples and use
	// REF/ALT" mode.
	sampleName := ""
	if samplesFile != "" {
		names, err := bcftools.LoadSamplesFile(samplesFile)
		if err != nil {
			fmt.Fprintln(os.Stderr, "bcftools consensus:", err)
			return 2
		}
		if len(names) > 0 {
			sampleName = names[0]
		}
	}
	if samples != "" && samples != "-" {
		ss := bcftools.SplitCommaList(samples)
		if len(ss) > 0 {
			sampleName = ss[0]
		}
	}

	var missingByte byte
	if missing != "" {
		if len(missing) != 1 {
			fmt.Fprintf(os.Stderr, "bcftools consensus: -M/--missing must be a single character (got %q)\n", missing)
			return 2
		}
		missingByte = missing[0]
	}
	var absentByte byte
	if absent != "" {
		if len(absent) != 1 {
			fmt.Fprintf(os.Stderr, "bcftools consensus: -a/--absent must be a single character (got %q)\n", absent)
			return 2
		}
		absentByte = absent[0]
	}

	opts := bcftools.ConsensusOptions{
		Sample:         sampleName,
		Haplotype:      hapSel,
		HaplotypeIndex: hapIdx,
		MaskWith:       maskWithSpec,
		MaskBED:        maskFile,
		MarkIns:        markInsSpec,
		MarkSnv:        markSnvSpec,
		MarkDel:        markDelSpec,
		Missing:        missingByte,
		Absent:         absentByte,
		Prefix:         prefix,
		IUPACCodes:     iupacCodes,
		IncludeExpr:    includeExpr,
		ExcludeExpr:    excludeExpr,
	}
	if maskFile != "" {
		regs, err := bcftools.LoadMaskBED(maskFile)
		if err != nil {
			fmt.Fprintln(os.Stderr, "bcftools consensus:", err)
			return 2
		}
		opts.Mask = regs
	}

	out, err := openOutFile(outputPath)
	if err != nil {
		fmt.Fprintf(os.Stderr, "bcftools consensus: %v\n", err)
		return 1
	}
	defer out.Close()

	if _, err := bcftools.ConsensusFile(rest[0], fastaRef, out, opts); err != nil {
		fmt.Fprintf(os.Stderr, "bcftools consensus: %v\n", err)
		return 1
	}
	return 0
}

type checkConsensusDeferredInputs struct {
	chainFile string
	haplotype string
}

func checkConsensusDeferred(in checkConsensusDeferredInputs) string {
	switch {
	case in.chainFile != "":
		return "-c/--chain (liftover chain output)"
	case isPhasedIUPAC(in.haplotype):
		return "-H NpIu (phased-index / unphased-IUPAC)"
	}
	return ""
}

// isPhasedIUPAC matches upstream's "NpIu" haplotype encoding ("phased
// index, unphased IUPAC"). Example: "2pIu".
func isPhasedIUPAC(s string) bool {
	if len(s) < 3 {
		return false
	}
	// Must contain "pI" between an integer prefix and a "u" suffix.
	for i := 1; i < len(s)-1; i++ {
		if s[i] == 'p' && i+2 < len(s)+1 && len(s)-i-1 >= 2 && s[i+1] == 'I' && s[i+2] == 'u' && i+3 == len(s) {
			// Prefix must be a positive integer.
			for k := 0; k < i; k++ {
				if s[k] < '0' || s[k] > '9' {
					return false
				}
			}
			return i > 0
		}
	}
	return false
}

// parseSingleCharMark validates --mark-del (single character only).
func parseSingleCharMark(s, flag string) (bcftools.MarkSpec, error) {
	if s == "" {
		return bcftools.MarkSpec{}, nil
	}
	if len(s) != 1 {
		return bcftools.MarkSpec{}, fmt.Errorf("%s must be a single character (got %q)", flag, s)
	}
	return bcftools.MarkSpec{Mode: bcftools.MarkChar, Char: s[0]}, nil
}
