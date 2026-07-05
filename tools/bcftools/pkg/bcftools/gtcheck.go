package bcftools

// bcftools gtcheck — sample-identity / contamination check.
//
// Faithful port of reference_code/bcftools/vcfgtcheck.c. The output
// matches upstream's "#DCv2" discordance table byte-for-byte (modulo the
// non-reproducible "# This file was produced by ..." provenance header
// and the "INFO\tTime required ..." timing line, which depend on the
// command line, working directory and wall-clock time).
//
// Scoring model (vcfgtcheck.c process_line):
//
//   - Each sample's genotype at a biallelic site is reduced to a dosage
//     bitmask: 1<<dose, so dose 0 -> 1, dose 1 -> 2, dose 2 -> 4. PL
//     genotypes can be ambiguous (multiple minima), giving a bitmask with
//     several bits set. A missing genotype is 0.
//   - Two samples are concordant at a site when (qry_dsg & gt_dsg) != 0.
//
//   - With -E/--error-probability != 0 (DEFAULT is 40): the discordance
//     is a probability score. Each dosage/PL is turned into a vector of
//     three negative-log probabilities (one per genotype 0/0, 0/1, 1/1)
//     and the per-site contribution is min_G(qry_prob[G] + gt_prob[G]).
//     The column is printed in "%e" scientific notation.
//   - With -E 0 (only reachable via --distinctive-sites today): the
//     discordance is the integer count of discordant sites, printed "%u".
//
//   - The Average -log P(HWE) column (unless --no-HWE-prob) accumulates,
//     over MATCHING sites only, -log of the HWE genotype probability for
//     the matched dosage, using the site allele frequency from INFO/AC,AN
//     when present and otherwise counted from FORMAT/GT across all
//     samples of the record. The reported value is hwe_prob / nmatch.
//
// Tag selection mirrors upstream's -u/--use: GT or PL, auto-detected per
// record (GT preferred, falling back to PL when GT is absent, or the
// reverse when PL is requested). Diploid GT/PL only.
//
// Pair enumeration:
//   - -p/-P explicit pairs: emitted in the given order (qry, gt).
//   - cross-check (no -g): the lower triangle, query = sample[i],
//     gt = sample[j] for j < i, i.e. rows (S2,S1), (S3,S1), (S3,S2)...
//   - paired (-g): every query sample against every panel sample.
//
// --distinctive-sites finds the minimal set of sites that distinguish at
// least NUM sample pairs and emits the "#DS" block after the table.

import (
	"bufio"
	"fmt"
	"io"
	"math"
	"os"
	"sort"
	"strconv"
	"strings"

	bgzf "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// GtcheckOptions controls the behaviour of Gtcheck / GtcheckFile.
type GtcheckOptions struct {
	// GenotypesFile is the -g panel of "truth" genotypes. When empty,
	// gtcheck runs in cross-check mode: every query sample versus
	// every other query sample (lower triangle).
	GenotypesFile string

	// PairsSpec is the upstream -p comma list. Layout depends on
	// whether GenotypesFile is set: qry,gt[,qry,gt...] when -g is
	// given, qry,qry[,qry,qry...] otherwise.
	PairsSpec string
	// PairsFile is the -P sidecar of pairs, one tab-separated pair
	// per line. Blank lines and '#' comments are ignored.
	PairsFile string

	// UseTag picks the scoring tag, upstream -u/--use. One or two
	// comma-separated tags ("GT", "PL"): the first applies to the
	// query, the second (default = first) to the -g panel. Empty
	// means "auto" (GT preferred, PL fallback).
	UseTag string

	// AllSites is accepted for surface parity; upstream itself rejects
	// -a as "to be implemented", so the field is currently unused.
	AllSites bool
	// HomsOnly restricts scoring to sites where the panel genotype is
	// homozygous (-H/--homs-only).
	HomsOnly bool
	// NoHWEProb disables the HWE-probability column (set to 0.0).
	NoHWEProb bool
	// KeepRefs disables the monoallelic-site skip (upstream --keep-refs).
	KeepRefs bool

	// SamplesQry / SamplesGT restrict the cohorts. Both are post-filtered
	// against the input headers.
	SamplesQry []string
	SamplesGT  []string

	// Regions / Targets restrict processing to CHROM[:beg-end] (applied as
	// a post-filter). IncludeExpr / ExcludeExpr are upstream's -i/-e filter
	// expressions: they drop non-matching sites before scoring and accept
	// the `qry:` / `gt:` scope prefix (an unprefixed expression applies to
	// both the query and the genotype reader).
	Regions     []string
	RegionsFile string
	Targets     []string
	TargetsFile string
	IncludeExpr string
	ExcludeExpr string

	// DryRun mirrors upstream `--dry-run`: process exactly one record
	// then return (used for time-estimation). No table is written.
	DryRun bool

	// ErrorProbability is the phred-scaled -E. The effective default is
	// 40 (matching upstream): a zero value here means "unset, use 40".
	// To force the integer-mismatch path set ErrorProbabilityZero (or
	// use --distinctive-sites, which forces it). Negative is treated as
	// 0.
	ErrorProbability int
	// ErrorProbabilityZero forces -E 0 (integer mismatch scoring) even
	// though the ErrorProbability zero value otherwise means "default
	// 40". The CLI sets this when the user passes -E 0 explicitly.
	ErrorProbabilityZero bool

	// DistinctiveSites enables the --distinctive-sites block. The value
	// is the NUM field: a fraction (<=1) of pairs or an absolute count
	// (>1). Requires -p/-P. Setting it forces ErrorProbability to 0.
	DistinctiveSites float64
	// HasDistinctiveSites distinguishes "--distinctive-sites 0" (which
	// upstream rejects) from "not set".
	HasDistinctiveSites bool

	// NMatches is upstream --n-matches: print only the top INT matches
	// per query sample, sorted by average discordance. A negative value
	// sorts by HWE probability instead. 0 means unlimited.
	NMatches int

	// OutputType: "t" (default tab-delimited text) or "z" (BGZF-compressed
	// text), mirroring upstream vcfgtcheck.c's -O t|z. The same report
	// bytes are written either way; "z" wraps them in a BGZF stream.
	OutputType string

	// CompressLevel is the BGZF deflate level for OutputType "z" (the
	// optional digit suffix in -O z<N>). -1 selects the default level.
	CompressLevel int

	// Cluster enables `-c/--cluster MIN,MAX`: group the (cross-check) query
	// samples into clusters of putatively-identical individuals by their
	// pairwise discordance. ClusterMax is the maximum within-cluster average
	// error rate (samples no further apart than this are merged); ClusterMin
	// is the minimum between-cluster error expected to separate distinct
	// individuals (reported, and used as the threshold when ClusterMax is
	// negative — upstream's default is [0.23,-0.3], where the negative MAX
	// means "derive it"). Upstream leaves -c as an unimplemented error stub,
	// so the clustering and its output are this port's own design.
	Cluster    bool
	ClusterMin float64
	ClusterMax float64
}

// GtcheckPair captures one row of the #DCv2 output. Discordance is held
// both as the floating-point probability score (DiscScore) and, when the
// integer path is active, as the count (DiscCount); IsInteger selects
// which one is emitted.
type GtcheckPair struct {
	QuerySample     string
	GenotypedSample string
	DiscScore       float64
	DiscCount       int
	IsInteger       bool
	AvgLogPHWE      float64
	NumSites        int
	NumMatching     int

	// queryIdx/gtIdx are the indices into the chosen sample cohorts,
	// retained so --n-matches and cross-check ordering can be applied.
	queryIdx int
	gtIdx    int
}

// GtcheckResult is the full pair table.
type GtcheckResult struct {
	Pairs []GtcheckPair
}

// GtcheckFile is the file-aware entry point used by the CLI.
func GtcheckFile(queryPath string, out io.Writer, opts GtcheckOptions) (GtcheckResult, error) {
	if opts.GenotypesFile != "" {
		gIn, err := openVariantReader(opts.GenotypesFile)
		if err != nil {
			return GtcheckResult{}, fmt.Errorf("bcftools gtcheck: open -g %s: %w", opts.GenotypesFile, err)
		}
		defer gIn.Close()
		qIn, err := openVariantReader(queryPath)
		if err != nil {
			return GtcheckResult{}, fmt.Errorf("bcftools gtcheck: open %s: %w", queryPath, err)
		}
		defer qIn.Close()
		return GtcheckPaired(qIn, gIn, out, opts)
	}
	in, err := openVariantReader(queryPath)
	if err != nil {
		return GtcheckResult{}, fmt.Errorf("bcftools gtcheck: open %s: %w", queryPath, err)
	}
	defer in.Close()
	return Gtcheck(in, out, opts)
}

// openVariantReader opens a VCF/BCF/gzipped-VCF input.
func openVariantReader(path string) (io.ReadCloser, error) {
	return iohelper.OpenReader(path)
}

// Gtcheck is the cross-check entry point: every query sample vs every
// other query sample (lower triangle).
func Gtcheck(in io.Reader, out io.Writer, opts GtcheckOptions) (GtcheckResult, error) {
	hdr, vars, err := readAllVariants(in)
	if err != nil {
		return GtcheckResult{}, fmt.Errorf("bcftools gtcheck: %w", err)
	}
	qInc, qExc, _, _, ferr := compileGtcheckFilters(opts)
	if ferr != nil {
		return GtcheckResult{}, ferr
	}
	// In cross-check mode there is only the query reader, so just the
	// qry-scoped filter applies (a bare, unprefixed expression is qry-scoped
	// here too). Each dropped site is one skipped position.
	vars, dropped := applyGtcheckFilter(vars, qInc, qExc)
	return runGtcheck(hdr, vars, hdr, vars, out, opts, true, uint32(len(dropped)))
}

// GtcheckPaired is the -g panel entry point.
func GtcheckPaired(qIn, gIn io.Reader, out io.Writer, opts GtcheckOptions) (GtcheckResult, error) {
	hdrQ, varsQ, err := readAllVariants(qIn)
	if err != nil {
		return GtcheckResult{}, fmt.Errorf("bcftools gtcheck: query: %w", err)
	}
	hdrG, varsG, err := readAllVariants(gIn)
	if err != nil {
		return GtcheckResult{}, fmt.Errorf("bcftools gtcheck: panel: %w", err)
	}
	qInc, qExc, gInc, gExc, ferr := compileGtcheckFilters(opts)
	if ferr != nil {
		return GtcheckResult{}, ferr
	}
	var dropQ, dropG []*vcf.Variant
	varsQ, dropQ = applyGtcheckFilter(varsQ, qInc, qExc)
	varsG, dropG = applyGtcheckFilter(varsG, gInc, gExc)
	// Upstream counts one skipped site per merged position at which either
	// reader's record was filtered out (vcfgtcheck.c:1131-1133), so the skip
	// count is the number of distinct CHROM:POS positions dropped from either
	// file rather than the raw sum.
	skipped := countDistinctPositions(dropQ, dropG)
	return runGtcheck(hdrQ, varsQ, hdrG, varsG, out, opts, false, skipped)
}

// countDistinctPositions counts the distinct CHROM:POS positions across the
// supplied variant slices.
func countDistinctPositions(groups ...[]*vcf.Variant) uint32 {
	seen := make(map[string]struct{})
	for _, g := range groups {
		for _, v := range g {
			seen[v.Chrom+"\x00"+strconv.Itoa(v.Pos)] = struct{}{}
		}
	}
	return uint32(len(seen))
}

// gtcheckFilterScope splits an upstream gtcheck filter expression into its
// scope ("qry", "gt", or "" for an unprefixed expression that applies to
// both readers) and the bare expression, mirroring vcfgtcheck.c's handling of
// the `qry:` / `gt:` prefix on -i/-e.
func gtcheckFilterScope(expr string) (scope, bare string) {
	switch {
	case strings.HasPrefix(expr, "qry:"):
		return "qry", expr[len("qry:"):]
	case strings.HasPrefix(expr, "gt:"):
		return "gt", expr[len("gt:"):]
	default:
		return "", expr
	}
}

// compileGtcheckFilters compiles the -i/--include and -e/--exclude expressions
// into per-reader (query and genotype) include/exclude Filters, honouring the
// `qry:` / `gt:` scope prefix. An unprefixed expression applies to both
// readers, matching upstream (vcfgtcheck.c sets both qry_filter and gt_filter
// for a bare expression).
func compileGtcheckFilters(opts GtcheckOptions) (qInc, qExc, gInc, gExc *Filter, err error) {
	if opts.IncludeExpr != "" {
		scope, bare := gtcheckFilterScope(opts.IncludeExpr)
		f, cerr := CompileFilter(bare)
		if cerr != nil {
			return nil, nil, nil, nil, fmt.Errorf("bcftools gtcheck: --include: %w", cerr)
		}
		if scope == "qry" || scope == "" {
			qInc = f
		}
		if scope == "gt" || scope == "" {
			gInc = f
		}
	}
	if opts.ExcludeExpr != "" {
		scope, bare := gtcheckFilterScope(opts.ExcludeExpr)
		f, cerr := CompileFilter(bare)
		if cerr != nil {
			return nil, nil, nil, nil, fmt.Errorf("bcftools gtcheck: --exclude: %w", cerr)
		}
		if scope == "qry" || scope == "" {
			qExc = f
		}
		if scope == "gt" || scope == "" {
			gExc = f
		}
	}
	return qInc, qExc, gInc, gExc, nil
}

// applyGtcheckFilter splits vars into the variants that pass the include
// filter (if any) and fail the exclude filter (if any) — kept — and those that
// do not — dropped — mirroring upstream's per-site filter_test gate that skips
// non-matching sites before scoring.
func applyGtcheckFilter(vars []*vcf.Variant, inc, exc *Filter) (kept, dropped []*vcf.Variant) {
	if inc == nil && exc == nil {
		return vars, nil
	}
	for _, v := range vars {
		if (inc == nil || inc.Eval(v)) && (exc == nil || !exc.Eval(v)) {
			kept = append(kept, v)
		} else {
			dropped = append(dropped, v)
		}
	}
	return kept, dropped
}

// tagMode is a per-sample-set tag selection: GT, PL, or auto.
type tagMode int

const (
	tagAuto tagMode = -1
	tagPL   tagMode = 0
	tagGT   tagMode = 1
)

// parseUseTag parses upstream's -u/--use spec into (query, panel) modes.
func parseUseTag(spec string) (qry, gt tagMode, err error) {
	s := strings.TrimSpace(spec)
	if s == "" {
		return tagAuto, tagAuto, nil
	}
	parts := strings.Split(s, ",")
	if len(parts) > 2 {
		return 0, 0, fmt.Errorf("bcftools gtcheck: failed to parse --use %q", spec)
	}
	parse1 := func(t string) (tagMode, error) {
		switch strings.ToUpper(strings.TrimSpace(t)) {
		case "GT":
			return tagGT, nil
		case "PL":
			return tagPL, nil
		default:
			return 0, fmt.Errorf("bcftools gtcheck: failed to parse --use %q; only GT and PL are supported", spec)
		}
	}
	qry, err = parse1(parts[0])
	if err != nil {
		return 0, 0, err
	}
	if len(parts) == 2 {
		gt, err = parse1(parts[1])
		if err != nil {
			return 0, 0, err
		}
	} else {
		gt = qry
	}
	return qry, gt, nil
}

// gtcheckState carries the scoring lookup tables for one run.
type gtcheckState struct {
	opts     GtcheckOptions
	useErr   bool // gt_err != 0: probability scoring
	calcHWE  bool
	dsg2prob map[uint8][3]float64
}

// effectiveErrorProbability resolves the phred-scaled -E that scoring
// should use. The zero value of ErrorProbability means "default 40";
// ErrorProbabilityZero forces 0 (the integer path). A negative value is
// clamped to 0.
func effectiveErrorProbability(opts GtcheckOptions) int {
	if opts.ErrorProbabilityZero {
		return 0
	}
	e := opts.ErrorProbability
	if e == 0 {
		return 40
	}
	if e < 0 {
		return 0
	}
	return e
}

// newState builds the scoring lookup tables.
func newState(opts GtcheckOptions) *gtcheckState {
	st := &gtcheckState{opts: opts}
	gtErr := effectiveErrorProbability(opts)
	st.calcHWE = !opts.NoHWEProb
	if gtErr != 0 {
		st.useErr = true
		eprob := math.Pow(10, -0.1*float64(gtErr))
		l := -math.Log(eprob)
		inf := math.Inf(1)
		// Index by the single-bit dosage mask (1,2,4). Mirrors
		// upstream args->dsg2prob.
		st.dsg2prob = map[uint8][3]float64{
			0: {inf, inf, inf},
			1: {0, l, 2 * l},
			2: {l, 0, l},
			4: {2 * l, l, 0},
		}
	}
	return st
}

// runGtcheck is the shared scoring path.
func runGtcheck(
	hdrQ *vcf.Header, varsQ []*vcf.Variant,
	hdrG *vcf.Header, varsG []*vcf.Variant,
	out io.Writer, opts GtcheckOptions, crossCheck bool,
	filterSkipped uint32,
) (GtcheckResult, error) {
	qryMode, gtMode, err := parseUseTag(opts.UseTag)
	if err != nil {
		return GtcheckResult{}, err
	}

	// Resolve the auto (-u unset) tag selection from the headers, exactly
	// once, mirroring upstream init_data(): the query prefers PL (falls
	// back to GT), the panel prefers GT (falls back to PL). In
	// cross-check mode the panel tag follows the query tag.
	qryUseGT, err := resolveHeaderTag(qryMode, varsQ, true, "query")
	if err != nil {
		return GtcheckResult{}, err
	}
	var gtUseGT bool
	if crossCheck {
		gtUseGT = qryUseGT
	} else {
		gtUseGT, err = resolveHeaderTag(gtMode, varsG, false, "panel")
		if err != nil {
			return GtcheckResult{}, err
		}
	}

	if opts.RegionsFile != "" {
		extra, err := LoadRegionsFile(opts.RegionsFile)
		if err != nil {
			return GtcheckResult{}, fmt.Errorf("bcftools gtcheck: regions-file: %w", err)
		}
		opts.Regions = append(opts.Regions, extra...)
	}
	if opts.TargetsFile != "" {
		extra, err := LoadRegionsFile(opts.TargetsFile)
		if err != nil {
			return GtcheckResult{}, fmt.Errorf("bcftools gtcheck: targets-file: %w", err)
		}
		opts.Targets = append(opts.Targets, extra...)
	}

	qSamples := selectSamples(hdrQ.Samples, opts.SamplesQry)
	gSamples := qSamples
	if !crossCheck {
		gSamples = selectSamples(hdrG.Samples, opts.SamplesGT)
	}
	if len(qSamples) == 0 || len(gSamples) == 0 {
		return GtcheckResult{}, fmt.Errorf("bcftools gtcheck: no samples to compare (qry=%d gt=%d)", len(qSamples), len(gSamples))
	}

	usePairs := opts.PairsSpec != "" || opts.PairsFile != ""
	if opts.HasDistinctiveSites && !usePairs {
		return GtcheckResult{}, fmt.Errorf("bcftools gtcheck: the experimental option --distinctive-sites requires -p/-P")
	}
	if opts.HomsOnly && crossCheck {
		return GtcheckResult{}, fmt.Errorf("bcftools gtcheck: the option --homs-only requires --genotypes")
	}

	// --distinctive-sites forces the integer path (gt_err = 0).
	if opts.HasDistinctiveSites {
		opts.ErrorProbabilityZero = true
	}
	st := newState(opts)

	pairs, err := buildPairList(qSamples, gSamples, opts, crossCheck)
	if err != nil {
		return GtcheckResult{}, err
	}
	if len(pairs) == 0 {
		// A plain cross-check over a single query sample resolves to zero
		// output pairs: upstream's report loop uses `ngt = cross_check ? i`,
		// which is 0 for the sole sample, so no DCv2 data row is written.
		// Upstream still exits 0 and still scores every site into the INFO
		// stats block (process_line increments ncmp/nused per site,
		// independent of the pair count), so we must NOT bail out here —
		// fall through to the normal site-scoring loop and emit the empty
		// data table. The error path is preserved only when -p/-P explicit
		// pairs or -g genotypes were requested but resolved to zero pairs;
		// that remains a genuine user error.
		if usePairs || !crossCheck {
			return GtcheckResult{}, fmt.Errorf("bcftools gtcheck: no sample pairs to compare")
		}
	}

	// Validate the resolved --distinctive-sites threshold, mirroring
	// upstream diff_sites_init: a NUM that resolves to <= 0 sites is an
	// error.
	if opts.HasDistinctiveSites {
		if n := resolveDistinctiveNsites(opts.DistinctiveSites, len(pairs)); n <= 0 {
			return GtcheckResult{}, fmt.Errorf("bcftools gtcheck: the value for --distinctive-sites was set too low: %d", n)
		}
	}

	panelByPos := indexByPos(varsG)
	qSampleIdx := indexSamples(hdrQ)
	gSampleIdx := indexSamples(hdrG)

	for i := range pairs {
		qi, ok := qSampleIdx[pairs[i].QuerySample]
		if !ok {
			return GtcheckResult{}, fmt.Errorf("bcftools gtcheck: query sample %q not in header", pairs[i].QuerySample)
		}
		gi, ok2 := gSampleIdx[pairs[i].GenotypedSample]
		if !ok2 {
			return GtcheckResult{}, fmt.Errorf("bcftools gtcheck: genotyped sample %q not in header", pairs[i].GenotypedSample)
		}
		pairs[i].queryIdx = qi
		pairs[i].gtIdx = gi
		pairs[i].IsInteger = !st.useErr
	}

	// Upstream sorts explicit -p/-P pairs by (query index, gt index)
	// before scoring (qsort cmp_pair). The auto-generated cross-check /
	// paired lists are already in that order.
	if usePairs {
		sort.SliceStable(pairs, func(a, b int) bool {
			if pairs[a].queryIdx != pairs[b].queryIdx {
				return pairs[a].queryIdx < pairs[b].queryIdx
			}
			return pairs[a].gtIdx < pairs[b].gtIdx
		})
	}

	stats := &gtcheckStats{skipFilter: filterSkipped}
	ds := newDistinctiveCollector(opts, len(pairs))

	sites := 0
	for _, qv := range varsQ {
		if len(opts.Regions) > 0 && !regionMatches(qv, opts.Regions) {
			continue
		}
		if len(opts.Targets) > 0 && !regionMatches(qv, opts.Targets) {
			continue
		}
		gv := qv
		if !crossCheck {
			gv = panelByPos[posKey(qv)]
			if gv == nil {
				stats.skipNoMatch++
				continue
			}
		}

		// Skip multi-allelic sites, mirroring upstream is_input_okay()'s
		// `n_allele>2` check (vcfgtcheck.c:1162): if any reader's record
		// carries more than one ALT allele the site is counted in
		// sites-skipped-multiallelic and skipped (upstream exits 0 and
		// merely advises `bcftools norm -m -`, it does not error). This
		// check runs BEFORE the monoallelic check to match upstream order.
		if isMultiAllelic(qv) || (!crossCheck && isMultiAllelic(gv)) {
			stats.skipNotBA++
			continue
		}

		// Skip monoallelic (no-ALT) sites unless --keep-refs, mirroring
		// upstream is_input_okay(): a site is skipped when every reader's
		// record has no ALT allele.
		if !opts.KeepRefs && isMonoallelic(qv) && (crossCheck || isMonoallelic(gv)) {
			stats.skipMono++
			continue
		}

		// Per-record tag with the upstream set_data() fallback: use the
		// globally resolved tag, but flip to the other tag if this
		// record lacks it.
		qUseGT, qOK := recordTag(qv, qryUseGT)
		if !qOK {
			stats.skipNoData++
			continue
		}
		gUseGT := qUseGT
		gOK := true
		if !crossCheck {
			gUseGT, gOK = recordTag(gv, gtUseGT)
		}
		if !gOK {
			stats.skipNoData++
			continue
		}

		sites++
		stats.compared++
		stats.bumpUsed(qUseGT, gUseGT)

		var hweDsg [8]float64
		if st.calcHWE {
			hweDsg = computeHWEDsg(siteAF(gv))
		}

		ds.resetSite()
		for pi := range pairs {
			p := &pairs[pi]
			gDsg, gProb := sampleDsgProb(gv, p.gtIdx, gUseGT, st)
			if gDsg == 0 {
				continue
			}
			if opts.HomsOnly && gDsg&5 == 0 {
				continue
			}
			qDsg, qProb := sampleDsgProb(qv, p.queryIdx, qUseGT, st)
			if qDsg == 0 {
				continue
			}

			if st.useErr {
				min := qProb[0] + gProb[0]
				if v := qProb[1] + gProb[1]; v < min {
					min = v
				}
				if v := qProb[2] + gProb[2]; v < min {
					min = v
				}
				p.DiscScore += min
			}
			match := qDsg & gDsg
			if match == 0 {
				if !st.useErr {
					p.DiscCount++
				}
				ds.markDiff(pi)
			} else if st.calcHWE {
				p.AvgLogPHWE += hweDsg[match]
				p.NumMatching++
			}
			p.NumSites++
		}
		ds.pushSite(qv.Chrom, qv.Pos)

		if opts.DryRun {
			break
		}
	}

	if opts.DryRun {
		// Upstream writes no table on --dry-run.
		return GtcheckResult{Pairs: pairs}, nil
	}

	// Finalise averaged HWE.
	result := GtcheckResult{Pairs: pairs}
	for i := range result.Pairs {
		p := &result.Pairs[i]
		if st.calcHWE && p.NumMatching > 0 {
			p.AvgLogPHWE = p.AvgLogPHWE / float64(p.NumMatching)
		} else {
			p.AvgLogPHWE = 0
		}
		if !st.calcHWE {
			p.NumMatching = 0
		}
	}

	if err := emitGtcheckOutput(out, result, stats, st, opts, crossCheck, qSamples, gSamples, ds); err != nil {
		return result, err
	}
	return result, nil
}

// emitGtcheckOutput writes the gtcheck report to out, honouring -O z / -O t.
//
// -O z wraps the identical report bytes in a BGZF stream, exactly as upstream
// opens its single output handle via bgzf_open(.., "wg") (vcfgtcheck.c:445) and
// writes the same text through it. -O t writes plain bytes (upstream's "wu"
// handle, an uncompressed passthrough). The framed result decodes
// byte-identically to the text output.
func emitGtcheckOutput(
	out io.Writer, result GtcheckResult, stats *gtcheckStats, st *gtcheckState,
	opts GtcheckOptions, crossCheck bool, qSamples, gSamples []string, ds *distinctiveCollector,
) error {
	if opts.OutputType == "z" {
		level := opts.CompressLevel
		if level < 0 {
			level = bgzf.DefaultCompression
		}
		bz, berr := bgzf.NewWriterLevel(out, level)
		if berr != nil {
			return berr
		}
		if werr := writeGtcheckBody(bz, result, stats, st, opts, crossCheck, qSamples, gSamples, ds); werr != nil {
			bz.Close()
			return werr
		}
		return bz.Close()
	}
	return writeGtcheckBody(out, result, stats, st, opts, crossCheck, qSamples, gSamples, ds)
}

// writeGtcheckBody writes the gtcheck report, the optional distinctive-sites
// block, and the optional cluster block to out, in upstream's order. It is
// shared by the plain-text (-O t) and BGZF (-O z) output paths so both emit
// byte-identical content.
func writeGtcheckBody(
	out io.Writer, result GtcheckResult, stats *gtcheckStats, st *gtcheckState,
	opts GtcheckOptions, crossCheck bool, qSamples, gSamples []string, ds *distinctiveCollector,
) error {
	if err := writeGtcheckReport(out, result, stats, st, opts, crossCheck, qSamples, gSamples); err != nil {
		return err
	}
	if opts.HasDistinctiveSites {
		if err := ds.report(out); err != nil {
			return err
		}
	}
	if opts.Cluster && crossCheck {
		if err := writeGtcheckClusters(out, result, opts); err != nil {
			return err
		}
	}
	return nil
}

// pairAvgError returns the per-pair average discordance (an error rate in
// [0,1]): the integer mismatch count or the floating-point discordance score,
// divided by the number of compared sites. Pairs with no compared sites are
// treated as maximally discordant (1) so they never cluster together.
func pairAvgError(p GtcheckPair) float64 {
	if p.NumSites <= 0 {
		return 1
	}
	if p.IsInteger {
		return float64(p.DiscCount) / float64(p.NumSites)
	}
	return p.DiscScore / float64(p.NumSites)
}

// writeGtcheckClusters groups the cross-check query samples into clusters of
// putatively-identical individuals by single-linkage over the pairwise average
// error: two samples join the same cluster when their error is at most the
// effective max threshold (opts.ClusterMax, or opts.ClusterMin when ClusterMax
// is negative — upstream's "derive it" sentinel). The CLUSTER section lists
// each cluster's members and its mean within-cluster error. Upstream's -c is an
// unimplemented error stub, so this format is this port's own design.
func writeGtcheckClusters(out io.Writer, result GtcheckResult, opts GtcheckOptions) error {
	thr := opts.ClusterMax
	if thr < 0 {
		thr = opts.ClusterMin
	}

	// Collect the sample set (cross-check pairs are query×query) and index it.
	idx := map[string]int{}
	var samples []string
	add := func(name string) {
		if _, ok := idx[name]; !ok {
			idx[name] = len(samples)
			samples = append(samples, name)
		}
	}
	for _, p := range result.Pairs {
		add(p.QuerySample)
		add(p.GenotypedSample)
	}

	// Union-find single-linkage: merge any pair within the threshold.
	parent := make([]int, len(samples))
	for i := range parent {
		parent[i] = i
	}
	var find func(int) int
	find = func(x int) int {
		for parent[x] != x {
			parent[x] = parent[parent[x]]
			x = parent[x]
		}
		return x
	}
	union := func(a, b int) {
		ra, rb := find(a), find(b)
		if ra != rb {
			parent[ra] = rb
		}
	}
	for _, p := range result.Pairs {
		if pairAvgError(p) <= thr {
			union(idx[p.QuerySample], idx[p.GenotypedSample])
		}
	}

	// Group samples by cluster root, preserving first-seen order.
	members := map[int][]string{}
	var roots []int
	for i, name := range samples {
		r := find(i)
		if _, ok := members[r]; !ok {
			roots = append(roots, r)
		}
		members[r] = append(members[r], name)
	}

	// Per-cluster mean within-cluster error (over the pairs whose both ends
	// fall in the cluster), for reporting.
	clusterOf := func(name string) int { return find(idx[name]) }
	sumErr := map[int]float64{}
	cntErr := map[int]int{}
	for _, p := range result.Pairs {
		if clusterOf(p.QuerySample) == clusterOf(p.GenotypedSample) {
			r := clusterOf(p.QuerySample)
			sumErr[r] += pairAvgError(p)
			cntErr[r]++
		}
	}

	bw := bufio.NewWriter(out)
	defer bw.Flush()
	fmt.Fprintf(bw, "# Clustering of query samples by pairwise discordance (-c %g,%g; threshold %g).\n",
		opts.ClusterMin, opts.ClusterMax, thr)
	fmt.Fprintln(bw, "# CLUSTER, [2]Cluster, [3]N samples, [4]Mean within-cluster error, [5]Samples")
	for i, r := range roots {
		mean := 0.0
		if cntErr[r] > 0 {
			mean = sumErr[r] / float64(cntErr[r])
		}
		fmt.Fprintf(bw, "CLUSTER\t%d\t%d\t%.6g\t%s\n",
			i+1, len(members[r]), mean, strings.Join(members[r], ","))
	}
	return nil
}

// resolveHeaderTag picks GT or PL for an entire cohort, mirroring
// upstream init_data(). It returns true for GT, false for PL. An explicit
// mode (from -u) is honoured (erroring if that tag is absent from every
// record). The auto mode prefers PL for the query cohort and GT for the
// panel cohort, falling back to the other tag.
func resolveHeaderTag(mode tagMode, vars []*vcf.Variant, isQuery bool, label string) (bool, error) {
	hasGT := cohortHasTag(vars, "GT")
	hasPL := cohortHasTag(vars, "PL")
	switch mode {
	case tagGT:
		if !hasGT {
			return false, fmt.Errorf("bcftools gtcheck: the GT tag is not present in the %s", label)
		}
		return true, nil
	case tagPL:
		if !hasPL {
			return false, fmt.Errorf("bcftools gtcheck: the PL tag is not present in the %s", label)
		}
		return false, nil
	default: // auto
		if isQuery {
			if hasPL {
				return false, nil
			}
			if hasGT {
				return true, nil
			}
		} else {
			if hasGT {
				return true, nil
			}
			if hasPL {
				return false, nil
			}
		}
		return false, fmt.Errorf("bcftools gtcheck: neither PL nor GT tag is present in the %s", label)
	}
}

// recordTag applies upstream's per-record set_data() fallback: prefer the
// cohort-resolved tag, but flip to the other if this record lacks it.
// Returns (useGT, ok); ok is false when neither tag is present.
func recordTag(v *vcf.Variant, preferGT bool) (useGT bool, ok bool) {
	if preferGT {
		if recordHasTag(v, "GT") {
			return true, true
		}
		if recordHasTag(v, "PL") {
			return false, true
		}
		return false, false
	}
	if recordHasTag(v, "PL") {
		return false, true
	}
	if recordHasTag(v, "GT") {
		return true, true
	}
	return false, false
}

func recordHasTag(v *vcf.Variant, tag string) bool {
	for i := range v.Samples {
		if _, ok := v.Samples[i].Data[tag]; ok {
			return true
		}
	}
	return false
}

// cohortHasTag reports whether any record carries the FORMAT tag. Upstream
// checks the header; the htsgo VCF model does not retain a typed FORMAT
// header set, so we probe the records (equivalent for well-formed input).
func cohortHasTag(vars []*vcf.Variant, tag string) bool {
	for _, v := range vars {
		if recordHasTag(v, tag) {
			return true
		}
	}
	return false
}

// sampleDsgProb returns the dosage bitmask and the three negative-log
// probabilities (only meaningful when st.useErr) for a sample.
func sampleDsgProb(v *vcf.Variant, idx int, useGT bool, st *gtcheckState) (uint8, [3]float64) {
	var prob [3]float64
	if v == nil || idx < 0 || idx >= len(v.Samples) {
		return 0, prob
	}
	if useGT {
		dsg := gtToDsg(v.Samples[idx].Data["GT"])
		if dsg != 0 && st.useErr {
			prob = st.dsg2prob[dsg]
		}
		return dsg, prob
	}
	pls, ok := parsePL(v.Samples[idx].Data["PL"])
	if !ok {
		return 0, prob
	}
	dsg := plToDsg(pls)
	if dsg != 0 && st.useErr {
		prob = plToProb(pls)
	}
	return dsg, prob
}

// gtToDsg converts a diploid FORMAT/GT string to a single-bit dosage
// bitmask (1<<dose). Returns 0 for missing or non-diploid.
func gtToDsg(gt string) uint8 {
	if gt == "" {
		return 0
	}
	gt = strings.ReplaceAll(gt, "|", "/")
	parts := strings.Split(gt, "/")
	if len(parts) != 2 {
		return 0
	}
	dose := 0
	for _, p := range parts {
		switch p {
		case ".":
			return 0
		case "0":
			// no-op
		case "1":
			dose++
		default:
			return 0
		}
	}
	return 1 << uint(dose)
}

// parsePL parses a diploid (3-value) FORMAT/PL string. Returns false if
// missing or not exactly three values.
func parsePL(s string) ([3]int, bool) {
	var out [3]int
	if s == "" || s == "." {
		return out, false
	}
	parts := strings.Split(s, ",")
	if len(parts) != 3 {
		return out, false
	}
	for i, p := range parts {
		if p == "." {
			return out, false
		}
		n, err := strconv.Atoi(p)
		if err != nil {
			return out, false
		}
		out[i] = n
	}
	return out, true
}

// plToDsg mirrors upstream pl_to_dsg: the dosage bitmask is the set of
// genotypes attaining the minimum PL.
func plToDsg(pl [3]int) uint8 {
	min := pl[0]
	if pl[1] < min {
		min = pl[1]
	}
	if pl[2] < min {
		min = pl[2]
	}
	var dsg uint8
	if pl[0] == min {
		dsg |= 1
	}
	if pl[1] == min {
		dsg |= 2
	}
	if pl[2] == min {
		dsg |= 4
	}
	return dsg
}

// plToProb mirrors upstream pl_to_prob: convert PLs to linear
// probabilities, normalise, then take the negative log.
func plToProb(pl [3]int) [3]float64 {
	var prob [3]float64
	for i := 0; i < 3; i++ {
		p := pl[i]
		if p < 0 || p >= 255 {
			p = 255
		}
		prob[i] = math.Pow(10, -0.1*float64(p))
	}
	sum := prob[0] + prob[1] + prob[2]
	for i := 0; i < 3; i++ {
		prob[i] = -math.Log(prob[i] / sum)
	}
	return prob
}

// computeHWEDsg builds the hwe_dsg lookup keyed by the dosage-match
// bitmask, mirroring upstream's per-record construction.
func computeHWEDsg(af float64) [8]float64 {
	var hwe [3]float64
	hwe[0] = -math.Log((1 - af) * (1 - af))
	hwe[1] = -math.Log(2 * af * (1 - af))
	hwe[2] = -math.Log(af * af)
	var out [8]float64
	out[0] = 0
	inf := math.Inf(1)
	for i := 1; i < 8; i++ {
		out[i] = inf
		for k := 0; k < 3; k++ {
			if (1<<uint(k))&i != 0 && out[i] > hwe[k] {
				out[i] = hwe[k]
			}
		}
	}
	return out
}

// siteAF computes the ALT allele frequency for a record, preferring
// INFO/AC,AN and falling back to counting FORMAT/GT over all samples.
// Mirrors bcf_calc_ac(BCF_UN_INFO|BCF_UN_FMT). A site with no observed
// allele uses the upstream sentinel of 1e-6.
func siteAF(v *vcf.Variant) float64 {
	if v == nil {
		return 1e-6
	}
	ac0, ac1, ok := acFromInfo(v)
	if !ok {
		ac0, ac1 = acFromGT(v)
	}
	if len(v.Alt) < 1 || v.Alt[0] == "." || v.Alt[0] == "" {
		ac1 = 0
	}
	if ac0+ac1 == 0 {
		return 1e-6
	}
	return float64(ac1) / float64(ac0+ac1)
}

// acFromInfo reads INFO/AC and INFO/AN. For a biallelic record AC is the
// ALT count and AN-AC the REF count.
func acFromInfo(v *vcf.Variant) (ac0, ac1 int, ok bool) {
	acStr, hasAC := v.Info["AC"]
	anStr, hasAN := v.Info["AN"]
	if !hasAC || !hasAN {
		return 0, 0, false
	}
	acField := acStr
	if i := strings.IndexByte(acField, ','); i >= 0 {
		acField = acField[:i]
	}
	ac1, err := strconv.Atoi(strings.TrimSpace(acField))
	if err != nil {
		return 0, 0, false
	}
	an, err := strconv.Atoi(strings.TrimSpace(anStr))
	if err != nil {
		return 0, 0, false
	}
	ac0 = an - ac1
	if ac0 < 0 {
		ac0 = 0
	}
	return ac0, ac1, true
}

// acFromGT counts REF (ac0) and ALT (ac1) alleles from FORMAT/GT across
// all samples of the record.
func acFromGT(v *vcf.Variant) (ac0, ac1 int) {
	for i := range v.Samples {
		gt := v.Samples[i].Data["GT"]
		if gt == "" {
			continue
		}
		gt = strings.ReplaceAll(gt, "|", "/")
		for _, a := range strings.Split(gt, "/") {
			switch a {
			case "0":
				ac0++
			case "1":
				ac1++
			}
		}
	}
	return ac0, ac1
}

// validateUseTag is retained for callers/tests that only want to confirm
// a tag spec parses.
func validateUseTag(tag string) error {
	_, _, err := parseUseTag(tag)
	return err
}

// isMultiAllelic reports whether the record carries more than one ALT
// allele, mirroring upstream is_input_okay()'s `n_allele>2` check
// (vcfgtcheck.c:1162; n_allele counts REF+ALTs, so >2 means >1 ALT).
// Such sites are skipped and counted in sites-skipped-multiallelic
// rather than treated as a fatal error.
func isMultiAllelic(v *vcf.Variant) bool {
	if v == nil {
		return false
	}
	if len(v.Alt) > 1 {
		return true
	}
	if len(v.Alt) == 1 && strings.Contains(v.Alt[0], ",") {
		return true
	}
	return false
}

// isMonoallelic reports whether the record has no usable ALT allele
// (ALT is absent, ".", or "<NON_REF>"-style), i.e. bcf_get_variant_types
// would classify it as VCF_REF.
func isMonoallelic(v *vcf.Variant) bool {
	if v == nil {
		return true
	}
	if len(v.Alt) == 0 {
		return true
	}
	for _, a := range v.Alt {
		if a != "" && a != "." {
			return false
		}
	}
	return true
}

func posKey(v *vcf.Variant) string {
	return v.Chrom + "\x00" + posStr(v.Pos)
}

func indexByPos(vs []*vcf.Variant) map[string]*vcf.Variant {
	out := make(map[string]*vcf.Variant, len(vs))
	for _, v := range vs {
		out[posKey(v)] = v
	}
	return out
}

func indexSamples(hdr *vcf.Header) map[string]int {
	out := make(map[string]int, len(hdr.Samples))
	for i, s := range hdr.Samples {
		out[s] = i
	}
	return out
}

// selectSamples narrows hdrSamples to those listed in want, preserving
// the input header's order. An empty want means "all".
func selectSamples(hdrSamples, want []string) []string {
	if len(want) == 0 {
		return append([]string{}, hdrSamples...)
	}
	w := make(map[string]struct{}, len(want))
	for _, s := range want {
		w[s] = struct{}{}
	}
	out := make([]string, 0, len(want))
	for _, s := range hdrSamples {
		if _, ok := w[s]; ok {
			out = append(out, s)
		}
	}
	return out
}

// buildPairList assembles the (query, genotype) pair list. With -p/-P the
// explicit list is used. Otherwise cross-check uses the lower triangle
// (query = larger-indexed sample) and the paired (-g) mode uses every
// query vs every panel sample.
func buildPairList(qSamples, gSamples []string, opts GtcheckOptions, crossCheck bool) ([]GtcheckPair, error) {
	if opts.PairsSpec != "" || opts.PairsFile != "" {
		parts, err := loadPairs(opts.PairsSpec, opts.PairsFile)
		if err != nil {
			return nil, err
		}
		if len(parts)%2 != 0 {
			return nil, fmt.Errorf("bcftools gtcheck: pairs list must have an even count, got %d", len(parts))
		}
		pairs := make([]GtcheckPair, 0, len(parts)/2)
		for i := 0; i < len(parts); i += 2 {
			pairs = append(pairs, GtcheckPair{
				QuerySample:     parts[i],
				GenotypedSample: parts[i+1],
			})
		}
		return pairs, nil
	}
	if crossCheck {
		// Lower triangle: query = sample[i], gt = sample[j], j < i.
		pairs := make([]GtcheckPair, 0, len(qSamples)*(len(qSamples)-1)/2)
		for i := 0; i < len(qSamples); i++ {
			for j := 0; j < i; j++ {
				pairs = append(pairs, GtcheckPair{
					QuerySample:     qSamples[i],
					GenotypedSample: qSamples[j],
				})
			}
		}
		return pairs, nil
	}
	pairs := make([]GtcheckPair, 0, len(qSamples)*len(gSamples))
	for _, q := range qSamples {
		for _, g := range gSamples {
			pairs = append(pairs, GtcheckPair{
				QuerySample:     q,
				GenotypedSample: g,
			})
		}
	}
	return pairs, nil
}

// loadPairs merges -p (comma list) and -P (file) into one flat slice.
func loadPairs(spec, path string) ([]string, error) {
	var out []string
	if spec != "" {
		for _, p := range strings.Split(spec, ",") {
			p = strings.TrimSpace(p)
			if p != "" {
				out = append(out, p)
			}
		}
	}
	if path != "" {
		f, err := os.Open(path)
		if err != nil {
			return nil, fmt.Errorf("bcftools gtcheck: -P %s: %w", path, err)
		}
		defer f.Close()
		sc := bufio.NewScanner(f)
		sc.Buffer(make([]byte, 0, 64<<10), 1<<20)
		for sc.Scan() {
			line := strings.TrimSpace(sc.Text())
			if line == "" || strings.HasPrefix(line, "#") {
				continue
			}
			var fields []string
			if strings.Contains(line, "\t") {
				fields = strings.Split(line, "\t")
			} else if strings.Contains(line, ",") {
				fields = strings.Split(line, ",")
			} else {
				fields = strings.Fields(line)
			}
			for _, p := range fields {
				p = strings.TrimSpace(p)
				if p != "" {
					out = append(out, p)
				}
			}
		}
		if err := sc.Err(); err != nil {
			return nil, err
		}
	}
	return out, nil
}

// regionMatches returns true if the variant CHROM[:beg-end] is in any of
// the given specs. An empty spec list returns true (no filter).
func regionMatches(v *vcf.Variant, specs []string) bool {
	if len(specs) == 0 {
		return true
	}
	for _, s := range specs {
		chrom, beg, end, err := parseRegionSpec(s)
		if err != nil {
			continue
		}
		if v.Chrom != chrom {
			continue
		}
		if beg <= 0 && end <= 0 {
			return true
		}
		if v.Pos >= beg && (end <= 0 || v.Pos <= end) {
			return true
		}
	}
	return false
}

func parseRegionSpec(s string) (chrom string, beg, end int, err error) {
	if !strings.Contains(s, ":") {
		return s, 0, 0, nil
	}
	colon := strings.LastIndex(s, ":")
	chrom = s[:colon]
	rest := s[colon+1:]
	if dash := strings.Index(rest, "-"); dash >= 0 {
		_, e := fmt.Sscanf(rest[:dash], "%d", &beg)
		if e != nil {
			return "", 0, 0, e
		}
		_, e = fmt.Sscanf(rest[dash+1:], "%d", &end)
		if e != nil {
			return "", 0, 0, e
		}
		return chrom, beg, end, nil
	}
	_, e := fmt.Sscanf(rest, "%d", &beg)
	if e != nil {
		return "", 0, 0, e
	}
	return chrom, beg, beg, nil
}

// posStr is a small helper for position formatting.
func posStr(p int) string {
	return strconv.Itoa(p)
}
