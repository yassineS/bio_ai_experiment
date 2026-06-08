package bcftools

// bcftools gtcheck — sample-identity / contamination check.
//
// Upstream reference: reference_code/bcftools/vcfgtcheck.c. The output is
// a multi-section report:
//
//  1. An INFO block of site counters (sites-compared, sites-skipped-*,
//     sites-used-*-vs-*).
//  2. A "# DCv2, discordance version 2:" comment block describing each
//     column.
//  3. The "#DCv2" header row followed by one DCv2 data row per
//     (query, genotyped) sample pair.
//
// The discordance score is the default error-probability model
// (gt_err=40, calc_hwe_prob=1): genotypes are converted to per-dosage
// negative-log probabilities given a phred-scaled allele-read-error
// probability, and the pair score accumulates the minimum joint
// negative-log probability across the three diploid genotypes at each
// overlapping site. The HWE column accumulates the per-site negative-log
// HWE probability of the matching dosage, averaged over matching sites.
//
// Out of scope for v1, tracked in docs/PARITY_ROADMAP.md#bcftools:
//   - `-u PL` (likelihood-based scoring from FORMAT/PL),
//   - `--cluster`, `--distinctive-sites`, `--n-matches`,
//   - `-E/--error-probability` other than the default 40,
//   - `-z` compressed output.
//
// These are accepted at the CLI for surface-parity and rejected with a
// clear pointer to the roadmap (see subcmds_gtcheck_roh.go).

import (
	"bufio"
	"fmt"
	"io"
	"math"
	"math/bits"
	"os"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// gtcheckDefaultErr is the default phred-scaled allele-read error
// probability (upstream args->gt_err = 40).
const gtcheckDefaultErr = 40

// GtcheckOptions controls the behaviour of Gtcheck / GtcheckFile.
type GtcheckOptions struct {
	// GenotypesFile is the -g panel of "truth" genotypes. When empty,
	// gtcheck runs in cross-check mode: every query sample versus
	// every other query sample, sub-diagonal half (i>j).
	GenotypesFile string

	// PairsSpec is the upstream -p comma list. Layout depends on
	// whether GenotypesFile is set: qry,gt[,qry,gt...] when -g is
	// given, qry,gt[,qry,gt...] otherwise (pairs are always (qry,gt)).
	PairsSpec string
	// PairsFile is the -P sidecar of pairs, one tab-separated pair
	// per line. Blank lines and '#' comments are ignored.
	PairsFile string

	// UseTag picks the scoring tag: only "GT" is implemented in v1.
	// "PL" is parsed and rejected with a PARITY_ROADMAP pointer.
	UseTag string

	// AllSites includes every site, even those where the query GT is
	// missing for every sample. Off by default (matches upstream).
	AllSites bool
	// HomsOnly restricts scoring to sites where the panel genotype is
	// homozygous (-H/--homs-only).
	HomsOnly bool
	// NoHWEProb disables the HWE-probability column (set to 0.0).
	NoHWEProb bool

	// SamplesQry / SamplesGT restrict the cohorts. Both are post-filtered
	// against the input headers.
	SamplesQry []string
	SamplesGT  []string

	// Regions / Targets / Include / Exclude are accepted but applied
	// only as a post-filter on CHROM[:beg-end] in v1.
	Regions     []string
	RegionsFile string
	Targets     []string
	TargetsFile string
	IncludeExpr string
	ExcludeExpr string

	// DryRun mirrors upstream `--dry-run`: process exactly one record
	// then return (used for time-estimation).
	DryRun bool

	// ErrorProbability is the phred-scaled -E. Defaults to 40 when 0.
	ErrorProbability int

	// OutputType: "t" (default), "z" (bgzip-compressed text), or a
	// digit selecting bgzip compression level (matches upstream
	// vcfgtcheck.c -O parsing).
	OutputType string
	// NMatches caps the top-N matches printed per query sample in
	// cross-check mode (mirrors upstream `--n-matches`). Default 0
	// = no cap. A negative value triggers HWE-based sort upstream
	// — we capture only the magnitude and a sort flag.
	NMatches  int
	SortByHWE bool
	// DistinctiveSites is upstream `--distinctive-sites`: when set,
	// emit a DS section listing the smallest set of sites whose
	// genotype mismatches collectively distinguish the requested
	// number of pairs. Values in (0,1] are treated as a fraction
	// of total pairs; >1 are absolute counts. Mirrors
	// vcfgtcheck.c::diff_sites_init.
	DistinctiveSites float64
}

// distSiteRow is the per-site row collected during cross-check
// when --distinctive-sites is active.
type distSiteRow struct {
	ndiff int
	chrom string
	pos   int
	bits  []uint64
}

// GtcheckPair captures one DCv2 data row.
type GtcheckPair struct {
	QuerySample     string
	GenotypedSample string
	// Discordance is the error-probability discordance score
	// (smaller = better match).
	Discordance float64
	// AvgLogPHWE is the average negative-log HWE probability over the
	// matching sites for this pair.
	AvgLogPHWE float64
	// NumSites is the number of sites compared for this pair.
	NumSites int
	// NumMatching is the number of matching genotypes for this pair.
	NumMatching int
}

// GtcheckCounters mirrors the upstream INFO section.
type GtcheckCounters struct {
	SitesCompared       int
	SkippedNoMatch      int
	SkippedMultiallelic int
	SkippedMonoallelic  int
	SkippedNoData       int
	SkippedGTNotDiploid int
	SkippedPLNotDiploid int
	SkippedFilterExpr   int
	UsedPLvsPL          int
	UsedPLvsGT          int
	UsedGTvsPL          int
	UsedGTvsGT          int
}

// GtcheckResult is the full report: counters plus the pair table.
type GtcheckResult struct {
	Counters GtcheckCounters
	Pairs    []GtcheckPair
}

// GtcheckFile is the file-aware entry point used by the CLI. It opens
// the query (and optional panel) through iohelper, scores each pair,
// and writes the multi-section report to out.
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
// other query sample, sub-diagonal half (i>j using the input header's
// sample ordering).
func Gtcheck(in io.Reader, out io.Writer, opts GtcheckOptions) (GtcheckResult, error) {
	hdr, vars, err := readAllVariants(in)
	if err != nil {
		return GtcheckResult{}, fmt.Errorf("bcftools gtcheck: %w", err)
	}
	return runGtcheck(hdr, vars, hdr, vars, out, opts, true)
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
	return runGtcheck(hdrQ, varsQ, hdrG, varsG, out, opts, false)
}

// gtErrProbs builds the dsg-bitmask → per-genotype negative-log
// probability table for the given phred-scaled error probability,
// mirroring upstream's dsg2prob construction.
func gtErrProbs(phred int) map[int][3]float64 {
	eprob := math.Pow(10, -0.1*float64(phred))
	le := -math.Log(eprob)
	return map[int][3]float64{
		1: {0, le, 2 * le}, // dsg=00: P(00|0)=1, P(01|0)=e, P(11|0)=e^2
		2: {le, 0, le},     // dsg=01: P(00|1)=e, P(01|1)=1, P(11|1)=e
		4: {2 * le, le, 0}, // dsg=11: P(00|2)=e^2, P(01|2)=e, P(11|2)=1
	}
}

// gtcheckSample holds the per-sample dosage bitmask and the error-model
// probabilities for a single site.
type gtcheckSample struct {
	dsg  int        // 0 = missing, else 1<<dosage
	prob [3]float64 // negative-log P(G) for G in {00,01,11}
}

// runGtcheck is the shared scoring path. crossCheck=true means the
// query and panel are the same cohort and we emit only the sub-diagonal
// (i>j) half.
func runGtcheck(
	hdrQ *vcf.Header, varsQ []*vcf.Variant,
	hdrG *vcf.Header, varsG []*vcf.Variant,
	out io.Writer, opts GtcheckOptions, crossCheck bool,
) (GtcheckResult, error) {
	if err := validateUseTag(opts.UseTag); err != nil {
		return GtcheckResult{}, err
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

	phred := opts.ErrorProbability
	if phred == 0 {
		phred = gtcheckDefaultErr
	}
	dsg2prob := gtErrProbs(phred)

	qSamples := selectSamples(hdrQ.Samples, opts.SamplesQry)
	gSamples := qSamples
	if !crossCheck {
		gSamples = selectSamples(hdrG.Samples, opts.SamplesGT)
	}
	if len(qSamples) == 0 || len(gSamples) == 0 {
		return GtcheckResult{}, fmt.Errorf("bcftools gtcheck: no samples to compare (qry=%d gt=%d)", len(qSamples), len(gSamples))
	}

	pairs, explicit, err := buildPairs(qSamples, gSamples, opts, crossCheck)
	if err != nil {
		return GtcheckResult{}, err
	}

	// Build (CHROM,POS) -> *Variant lookup for the panel so cross-file
	// joining is O(1) per query record.
	panelByPos := indexByPos(varsG)
	qSampleIdx := indexSamples(hdrQ)
	gSampleIdx := indexSamples(hdrG)

	type acc struct {
		pdiff  float64
		match  int
		sites  int
		hweAcc float64
	}
	accums := make(map[[2]string]*acc, len(pairs))
	for _, p := range pairs {
		key := [2]string{p.QuerySample, p.GenotypedSample}
		accums[key] = &acc{}
	}
	// --distinctive-sites bookkeeping: per-site index of which
	// pairs were discordant. Allocated only when the option is
	// active (mirrors upstream diff_sites_init). Upstream requires
	// -p/-P paired mode; mirror that rejection.
	var distSites []distSiteRow
	distActive := opts.DistinctiveSites > 0
	if distActive && opts.PairsSpec == "" && opts.PairsFile == "" {
		return GtcheckResult{}, fmt.Errorf("bcftools gtcheck: the experimental option --distinctive-sites requires -p/-P")
	}
	// Resolve the target count: when 0 < f <= 1 treat as a
	// fraction of npairs; otherwise an absolute count.
	distTarget := 0
	if distActive {
		if opts.DistinctiveSites <= 1 {
			distTarget = int(opts.DistinctiveSites * float64(len(pairs)))
		} else {
			distTarget = int(opts.DistinctiveSites)
		}
		if distTarget <= 0 {
			return GtcheckResult{}, fmt.Errorf("bcftools gtcheck: --distinctive-sites value too low: %d", distTarget)
		}
		if distTarget > len(pairs) {
			distTarget = len(pairs)
		}
	}

	pairIndex := map[[2]string]int{}
	for i, p := range pairs {
		pairIndex[[2]string{p.QuerySample, p.GenotypedSample}] = i
	}

	var counters GtcheckCounters
	for _, qv := range varsQ {
		// Apply region / targets post-filters independently.
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
				counters.SkippedNoMatch++
				continue
			}
		}

		// Multi-allelic / mono-allelic site checks (mirrors
		// is_input_okay). n_allele>2 → skip; all REF → skip.
		if nAltAlleles(qv) > 1 || (!crossCheck && nAltAlleles(gv) > 1) {
			counters.SkippedMultiallelic++
			continue
		}
		if nAltAlleles(qv) == 0 && (crossCheck || nAltAlleles(gv) == 0) {
			counters.SkippedMonoallelic++
			continue
		}

		counters.SitesCompared++
		counters.UsedGTvsGT++

		// Per-site allele frequency for the HWE column.
		af := siteAF(gv)
		var hwe [3]float64
		if !opts.NoHWEProb {
			hwe[0] = -math.Log((1 - af) * (1 - af))
			hwe[1] = -math.Log(2 * af * (1 - af))
			hwe[2] = -math.Log(af * af)
		}

		// Pre-compute each sample's dosage + probs once per site.
		qData := make(map[string]gtcheckSample, len(qSamples))
		for _, s := range qSamples {
			qData[s] = sampleDosage(qv, qSampleIdx[s], dsg2prob)
		}
		var gData map[string]gtcheckSample
		if crossCheck {
			gData = qData
		} else {
			gData = make(map[string]gtcheckSample, len(gSamples))
			for _, s := range gSamples {
				gData[s] = sampleDosage(gv, gSampleIdx[s], dsg2prob)
			}
		}

		// Per-record bitset of discordant pairs (for
		// --distinctive-sites). Reset per record; populated
		// inside the pair loop.
		var diffBits []uint64
		if distActive {
			diffBits = make([]uint64, (len(pairs)+63)/64)
		}
		ndiffPerSite := 0

		for _, p := range pairs {
			a := accums[[2]string{p.QuerySample, p.GenotypedSample}]
			gs := gData[p.GenotypedSample]
			if gs.dsg == 0 {
				continue // missing panel value
			}
			if opts.HomsOnly && gs.dsg&5 == 0 {
				continue // not a hom
			}
			qs := qData[p.QuerySample]
			if qs.dsg == 0 {
				continue // missing query value
			}

			min := qs.prob[0] + gs.prob[0]
			if v := qs.prob[1] + gs.prob[1]; v < min {
				min = v
			}
			if v := qs.prob[2] + gs.prob[2]; v < min {
				min = v
			}
			a.pdiff += min

			matchedBits := qs.dsg & gs.dsg
			if !opts.NoHWEProb {
				if matchedBits != 0 {
					hd := math.Inf(1)
					for k := 0; k < 3; k++ {
						if (1<<k)&matchedBits != 0 && hwe[k] < hd {
							hd = hwe[k]
						}
					}
					a.hweAcc += hd
					a.match++
				}
			} else if matchedBits != 0 {
				a.match++
			}
			a.sites++

			// Discordance for --distinctive-sites: the pair is
			// distinguished at this site iff the dosages do
			// not share any bit.
			if distActive && matchedBits == 0 {
				idx := pairIndex[[2]string{p.QuerySample, p.GenotypedSample}]
				diffBits[idx>>6] |= 1 << (uint(idx) & 63)
				ndiffPerSite++
			}
		}
		if distActive && ndiffPerSite > 0 {
			distSites = append(distSites, distSiteRow{
				ndiff: ndiffPerSite,
				chrom: qv.Chrom,
				pos:   qv.Pos,
				bits:  diffBits,
			})
		}
		if opts.DryRun {
			break
		}
	}

	if !opts.AllSites && counters.SitesCompared == 0 && !explicit {
		// Upstream still prints an (empty) report; we mirror that by
		// continuing to the writer rather than erroring.
	}

	result := GtcheckResult{Counters: counters, Pairs: make([]GtcheckPair, 0, len(pairs))}
	for _, p := range pairs {
		a := accums[[2]string{p.QuerySample, p.GenotypedSample}]
		avgHWE := 0.0
		if a.match > 0 && !opts.NoHWEProb {
			avgHWE = a.hweAcc / float64(a.match)
		}
		result.Pairs = append(result.Pairs, GtcheckPair{
			QuerySample:     p.QuerySample,
			GenotypedSample: p.GenotypedSample,
			Discordance:     a.pdiff,
			AvgLogPHWE:      avgHWE,
			NumSites:        a.sites,
			NumMatching:     a.match,
		})
	}
	// --n-matches: per-query top-N filter (mirrors upstream
	// vcfgtcheck.c:980-1027). Only active in cross-check mode
	// (incompatible with -p/-P).
	if opts.NMatches > 0 && crossCheck && opts.PairsSpec == "" && opts.PairsFile == "" {
		result.Pairs = topNMatchesByQuery(result.Pairs, opts.NMatches, opts.SortByHWE)
	}
	if err := writeGtcheckReport(out, result); err != nil {
		return result, err
	}
	// --distinctive-sites: after the DCv2 body, emit the DS block
	// (greedy site selection covering distTarget pairs per block).
	if distActive && len(distSites) > 0 {
		if err := writeDistinctiveSites(out, distSites, len(pairs), distTarget); err != nil {
			return result, err
		}
	}
	return result, nil
}

// writeDistinctiveSites is the Go port of
// vcfgtcheck.c::report_distinctive_sites (line 832-871). It sorts the
// per-site discordant-pair bitsets by ndiff (descending), then greedily
// walks them assigning new pairs to a block until distTarget pairs are
// distinguished; at that point a new block starts. Emits one `DS` row
// per site that adds at least one new pair.
func writeDistinctiveSites(out io.Writer, sites []distSiteRow, npairs, distTarget int) error {
	// Sort by ndiff descending; tie-break by (chrom, pos) for
	// determinism (upstream uses a random tag — we use the natural
	// order so the output is reproducible across runs).
	sort.SliceStable(sites, func(i, j int) bool {
		if sites[i].ndiff != sites[j].ndiff {
			return sites[i].ndiff > sites[j].ndiff
		}
		if sites[i].chrom != sites[j].chrom {
			return sites[i].chrom < sites[j].chrom
		}
		return sites[i].pos < sites[j].pos
	})
	bw := bufio.NewWriter(out)
	bw.WriteString("# DS, distinctive sites:\n")
	bw.WriteString("#     - chromosome\n")
	bw.WriteString("#     - position\n")
	bw.WriteString("#     - cumulative number of pairs distinguished by this block\n")
	bw.WriteString("#     - block id\n")
	bw.WriteString("#DS\t[2]Chromosome\t[3]Position\t[4]Cumulative number of distinct pairs\t[5]Block id\n")
	blk := make([]uint64, (npairs+63)/64)
	ndiffTot := 0
	iblock := 0
	ndiffMin := distTarget
	if ndiffMin > npairs {
		ndiffMin = npairs
	}
	for _, s := range sites {
		ndiffNew := 0
		for w, mask := range s.bits {
			if mask == 0 {
				continue
			}
			// Bits set in this site but not yet in the block.
			fresh := mask & ^blk[w]
			if fresh == 0 {
				continue
			}
			blk[w] |= fresh
			ndiffNew += bits.OnesCount64(fresh)
		}
		if ndiffNew == 0 {
			continue
		}
		ndiffTot += ndiffNew
		fmt.Fprintf(bw, "DS\t%s\t%d\t%d\t%d\n", s.chrom, s.pos, ndiffTot, iblock)
		if ndiffTot < ndiffMin {
			continue
		}
		iblock++
		ndiffTot = 0
		for i := range blk {
			blk[i] = 0
		}
	}
	return bw.Flush()
}

// topNMatchesByQuery groups pairs by query sample, sorts each
// group by score (discordance ascending, or -avgHWE ascending when
// sortByHWE is true), and keeps the top N per group. Output order:
// (query first appearance, then top-N within each query).
func topNMatchesByQuery(pairs []GtcheckPair, n int, sortByHWE bool) []GtcheckPair {
	if n <= 0 || len(pairs) == 0 {
		return pairs
	}
	byQry := map[string][]GtcheckPair{}
	order := []string{}
	for _, p := range pairs {
		if _, ok := byQry[p.QuerySample]; !ok {
			order = append(order, p.QuerySample)
		}
		byQry[p.QuerySample] = append(byQry[p.QuerySample], p)
	}
	out := make([]GtcheckPair, 0, len(order)*n)
	for _, q := range order {
		group := byQry[q]
		sort.SliceStable(group, func(i, j int) bool {
			if sortByHWE {
				a := scoreByHWE(group[i])
				b := scoreByHWE(group[j])
				return a < b
			}
			a := scoreByDiscordance(group[i])
			b := scoreByDiscordance(group[j])
			return a < b
		})
		if len(group) > n {
			group = group[:n]
		}
		out = append(out, group...)
	}
	return out
}

func scoreByDiscordance(p GtcheckPair) float64 {
	if p.NumSites == 0 {
		return 0
	}
	return p.Discordance / float64(p.NumSites)
}

func scoreByHWE(p GtcheckPair) float64 {
	if p.NumMatching == 0 {
		return 0
	}
	return -p.AvgLogPHWE
}

func validateUseTag(tag string) error {
	t := strings.ToUpper(strings.TrimSpace(tag))
	if t == "" || t == "GT" {
		return nil
	}
	if t == "PL" || strings.HasPrefix(t, "PL,") || strings.HasSuffix(t, ",PL") {
		return fmt.Errorf("bcftools gtcheck: PL mode not implemented in v1; tracked in docs/PARITY_ROADMAP.md#bcftools")
	}
	return fmt.Errorf("bcftools gtcheck: -u %q not recognised; only GT is implemented (PL deferred, see docs/PARITY_ROADMAP.md#bcftools)", tag)
}

// nAltAlleles returns the number of real ALT alleles (excluding the "."
// / empty placeholder that denotes a monoallelic REF-only site).
func nAltAlleles(v *vcf.Variant) int {
	if v == nil {
		return 0
	}
	n := 0
	for _, a := range v.Alt {
		if a == "" || a == "." {
			continue
		}
		n++
	}
	return n
}

// siteAF computes the alternate-allele frequency used for the HWE
// column. It prefers INFO/AC and INFO/AN (matching bcf_calc_ac with
// BCF_UN_INFO) and falls back to counting FORMAT/GT alleles. A site with
// no observed alternate alleles uses the upstream 1e-6 sentinel.
func siteAF(v *vcf.Variant) float64 {
	if v == nil {
		return 1e-6
	}
	acStr := firstField(v.Info["AC"])
	anStr := v.Info["AN"]
	if acStr != "" && anStr != "" {
		ac, err1 := strconv.Atoi(acStr)
		an, err2 := strconv.Atoi(anStr)
		if err1 == nil && err2 == nil && an > 0 {
			if ac == 0 {
				return 1e-6
			}
			return float64(ac) / float64(an)
		}
	}
	// Fall back to counting GT alleles.
	var an, ac int
	for i := range v.Samples {
		gt, ok := v.Samples[i].Data["GT"]
		if !ok {
			continue
		}
		for _, p := range strings.FieldsFunc(gt, func(r rune) bool { return r == '/' || r == '|' }) {
			switch p {
			case ".":
				// missing allele, not counted
			case "0":
				an++
			default:
				an++
				ac++
			}
		}
	}
	if an == 0 || ac == 0 {
		return 1e-6
	}
	return float64(ac) / float64(an)
}

// firstField returns the substring before the first comma (the value
// for the first ALT allele in a Number=A INFO field).
func firstField(s string) string {
	if i := strings.IndexByte(s, ','); i >= 0 {
		return s[:i]
	}
	return s
}

// sampleDosage converts a sample's diploid FORMAT/GT into a gtcheckSample
// (dosage bitmask + per-genotype negative-log probabilities). A missing,
// haploid, or non-biallelic GT yields dsg=0 (skip).
func sampleDosage(v *vcf.Variant, idx int, dsg2prob map[int][3]float64) gtcheckSample {
	if v == nil || idx < 0 || idx >= len(v.Samples) {
		return gtcheckSample{}
	}
	gt, ok := v.Samples[idx].Data["GT"]
	if !ok || gt == "" {
		return gtcheckSample{}
	}
	gt = strings.ReplaceAll(gt, "|", "/")
	parts := strings.Split(gt, "/")
	if len(parts) != 2 {
		return gtcheckSample{} // only diploid GT supported
	}
	dose := 0
	for _, p := range parts {
		switch p {
		case ".":
			return gtcheckSample{} // missing
		case "0":
			// no-op
		case "1":
			dose++
		default:
			// allele index >1 (multi-allelic) — site already
			// rejected upstream of this call.
			return gtcheckSample{}
		}
	}
	dsg := 1 << dose
	return gtcheckSample{dsg: dsg, prob: dsg2prob[dsg]}
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

// buildPairs assembles the (query, genotype) pair list in the upstream
// output order. When -p/-P is set, pairs are taken verbatim (sorted by
// (iqry,igt) header index). Otherwise: cross-check emits the
// sub-diagonal (samples[i], samples[j]) for j<i; paired mode emits every
// (qry,gt) combination. The bool return is true when pairs came from
// -p/-P.
func buildPairs(qSamples, gSamples []string, opts GtcheckOptions, crossCheck bool) ([]GtcheckPair, bool, error) {
	if opts.PairsSpec != "" || opts.PairsFile != "" {
		parts, err := loadPairs(opts.PairsSpec, opts.PairsFile)
		if err != nil {
			return nil, true, err
		}
		if len(parts)%2 != 0 {
			return nil, true, fmt.Errorf("bcftools gtcheck: pairs list must have an even count, got %d", len(parts))
		}
		pairs := make([]GtcheckPair, 0, len(parts)/2)
		for i := 0; i < len(parts); i += 2 {
			pairs = append(pairs, GtcheckPair{
				QuerySample:     parts[i],
				GenotypedSample: parts[i+1],
			})
		}
		return pairs, true, nil
	}
	if crossCheck {
		// Sub-diagonal i>j in header order: query=samples[i],
		// genotyped=samples[j].
		pairs := make([]GtcheckPair, 0, len(qSamples)*(len(qSamples)-1)/2)
		for i := 0; i < len(qSamples); i++ {
			for j := 0; j < i; j++ {
				pairs = append(pairs, GtcheckPair{
					QuerySample:     qSamples[i],
					GenotypedSample: qSamples[j],
				})
			}
		}
		return pairs, false, nil
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
	return pairs, false, nil
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

// regionMatches returns true if the variant CHROM[:beg-end] is in any
// of the given specs. An empty spec list returns true (no filter).
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

// writeGtcheckReport emits the full upstream multi-section report:
// the INFO counter block, the DCv2 comment block, the DCv2 header row,
// and one DCv2 data row per pair. Floats use C `%e` formatting, which
// Go's `%e` reproduces byte-for-byte for these values.
func writeGtcheckReport(out io.Writer, r GtcheckResult) error {
	w := bufio.NewWriter(out)
	c := r.Counters
	fmt.Fprintf(w, "INFO\tsites-compared\t%d\n", c.SitesCompared)
	fmt.Fprintf(w, "INFO\tsites-skipped-no-match\t%d\n", c.SkippedNoMatch)
	fmt.Fprintf(w, "INFO\tsites-skipped-multiallelic\t%d\n", c.SkippedMultiallelic)
	fmt.Fprintf(w, "INFO\tsites-skipped-monoallelic\t%d\n", c.SkippedMonoallelic)
	fmt.Fprintf(w, "INFO\tsites-skipped-no-data\t%d\n", c.SkippedNoData)
	fmt.Fprintf(w, "INFO\tsites-skipped-GT-not-diploid\t%d\n", c.SkippedGTNotDiploid)
	fmt.Fprintf(w, "INFO\tsites-skipped-PL-not-diploid\t%d\n", c.SkippedPLNotDiploid)
	fmt.Fprintf(w, "INFO\tsites-skipped-filtering-expression\t%d\n", c.SkippedFilterExpr)
	fmt.Fprintf(w, "INFO\tsites-used-PL-vs-PL\t%d\n", c.UsedPLvsPL)
	fmt.Fprintf(w, "INFO\tsites-used-PL-vs-GT\t%d\n", c.UsedPLvsGT)
	fmt.Fprintf(w, "INFO\tsites-used-GT-vs-PL\t%d\n", c.UsedGTvsPL)
	fmt.Fprintf(w, "INFO\tsites-used-GT-vs-GT\t%d\n", c.UsedGTvsGT)
	w.WriteString("# DCv2, discordance version 2:\n")
	w.WriteString("#     - Query sample\n")
	w.WriteString("#     - Genotyped sample\n")
	w.WriteString("#     - Discordance, given either as an abstract score or number of mismatches, see the options -E/-u\n")
	w.WriteString("#       in man page for details. Note that samples with high missingness have fewer sites compared,\n")
	w.WriteString("#       which results in lower overall discordance. Therefore it is advisable to use the average score\n")
	w.WriteString("#       per site rather than the absolute value, i.e. divide the value by the number of sites compared\n")
	w.WriteString("#       (smaller value = better match)\n")
	w.WriteString("#     - Average negative log of HWE probability at matching sites, attempts to quantify the following\n")
	w.WriteString("#       intuition: rare genotype matches are more informative than common genotype matches, hence two\n")
	w.WriteString("#       samples with similar discordance can be further stratified by the HWE score (bigger value = better\n")
	w.WriteString("#       match, the observed concordance was less likely to occur by chance)\n")
	w.WriteString("#     - Number of sites compared for this pair of samples (bigger = more informative)\n")
	w.WriteString("#     - Number of matching genotypes\n")
	w.WriteString("#DCv2\t[2]Query Sample\t[3]Genotyped Sample\t[4]Discordance\t[5]Average -log P(HWE)\t[6]Number of sites compared\t[7]Number of matching genotypes\n")
	for _, p := range r.Pairs {
		fmt.Fprintf(w, "DCv2\t%s\t%s\t%e\t%e\t%d\t%d\n",
			p.QuerySample, p.GenotypedSample, p.Discordance, p.AvgLogPHWE, p.NumSites, p.NumMatching)
	}
	return w.Flush()
}

// posStr is a small helper to avoid pulling strconv just for one fmt.
func posStr(p int) string {
	return strconv.Itoa(p)
}
