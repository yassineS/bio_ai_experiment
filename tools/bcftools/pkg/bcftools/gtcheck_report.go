package bcftools

// Output formatting for bcftools gtcheck, mirroring vcfgtcheck.c report().
//
// The emitted bytes match upstream exactly EXCEPT for two
// non-reproducible lines that upstream writes and we intentionally omit:
//
//   - "# This file was produced by bcftools (...), the command line ..."
//     and the following provenance lines (depend on argv + cwd);
//   - "INFO\tTime required to process one record .. <seconds>" (wall
//     clock).
//
// Everything from the "INFO\tsites-compared" stats block onward is
// byte-for-byte identical.

import (
	"bufio"
	"fmt"
	"io"
	"sort"
)

// gtcheckStats accumulates the per-run INFO counters that upstream emits
// at the top of report().
type gtcheckStats struct {
	compared    uint32
	skipNoMatch uint32
	skipNotBA   uint32
	skipMono    uint32
	skipNoData  uint32
	skipDipGT   uint32
	skipDipPL   uint32
	skipFilter  uint32
	// used[qryGT][gtGT] counts which tag combination was used.
	used [2][2]uint32
}

// bumpUsed records the tag combination used for one site. The indices
// follow upstream's nused[qry_use_GT][gt_use_GT].
func (s *gtcheckStats) bumpUsed(qryGT, gtGT bool) {
	qi, gi := 0, 0
	if qryGT {
		qi = 1
	}
	if gtGT {
		gi = 1
	}
	s.used[qi][gi]++
}

// dcv2Header is the literal block upstream writes before the data rows.
const dcv2Header = `# DCv2, discordance version 2:
#     - Query sample
#     - Genotyped sample
#     - Discordance, given either as an abstract score or number of mismatches, see the options -E/-u
#       in man page for details. Note that samples with high missingness have fewer sites compared,
#       which results in lower overall discordance. Therefore it is advisable to use the average score
#       per site rather than the absolute value, i.e. divide the value by the number of sites compared
#       (smaller value = better match)
#     - Average negative log of HWE probability at matching sites, attempts to quantify the following
#       intuition: rare genotype matches are more informative than common genotype matches, hence two
#       samples with similar discordance can be further stratified by the HWE score (bigger value = better
#       match, the observed concordance was less likely to occur by chance)
#     - Number of sites compared for this pair of samples (bigger = more informative)
#     - Number of matching genotypes
#DCv2	[2]Query Sample	[3]Genotyped Sample	[4]Discordance	[5]Average -log P(HWE)	[6]Number of sites compared	[7]Number of matching genotypes
`

// writeGtcheckReport writes the INFO stats block, the DCv2 comment
// header and the discordance rows.
func writeGtcheckReport(
	out io.Writer, r GtcheckResult, stats *gtcheckStats, st *gtcheckState,
	opts GtcheckOptions, crossCheck bool, qSamples, gSamples []string,
) error {
	w := bufio.NewWriter(out)

	fmt.Fprintf(w, "INFO\tsites-compared\t%d\n", stats.compared)
	fmt.Fprintf(w, "INFO\tsites-skipped-no-match\t%d\n", stats.skipNoMatch)
	fmt.Fprintf(w, "INFO\tsites-skipped-multiallelic\t%d\n", stats.skipNotBA)
	fmt.Fprintf(w, "INFO\tsites-skipped-monoallelic\t%d\n", stats.skipMono)
	fmt.Fprintf(w, "INFO\tsites-skipped-no-data\t%d\n", stats.skipNoData)
	fmt.Fprintf(w, "INFO\tsites-skipped-GT-not-diploid\t%d\n", stats.skipDipGT)
	fmt.Fprintf(w, "INFO\tsites-skipped-PL-not-diploid\t%d\n", stats.skipDipPL)
	fmt.Fprintf(w, "INFO\tsites-skipped-filtering-expression\t%d\n", stats.skipFilter)
	fmt.Fprintf(w, "INFO\tsites-used-PL-vs-PL\t%d\n", stats.used[0][0])
	fmt.Fprintf(w, "INFO\tsites-used-PL-vs-GT\t%d\n", stats.used[0][1])
	fmt.Fprintf(w, "INFO\tsites-used-GT-vs-PL\t%d\n", stats.used[1][0])
	fmt.Fprintf(w, "INFO\tsites-used-GT-vs-GT\t%d\n", stats.used[1][1])
	if _, err := w.WriteString(dcv2Header); err != nil {
		return err
	}

	rows := selectReportRows(r, opts, crossCheck, qSamples, gSamples)
	for _, p := range rows {
		if err := writeDCv2Row(w, p, st.calcHWE); err != nil {
			return err
		}
	}
	return w.Flush()
}

// writeDCv2Row writes one DCv2 data line. The Discordance column is "%u"
// in the integer path and "%e" otherwise; the HWE column is always "%e".
func writeDCv2Row(w *bufio.Writer, p GtcheckPair, calcHWE bool) error {
	hwe := 0.0
	nmatch := 0
	if calcHWE {
		hwe = p.AvgLogPHWE
		nmatch = p.NumMatching
	}
	if p.IsInteger {
		_, err := fmt.Fprintf(w, "DCv2\t%s\t%s\t%d\t%s\t%d\t%d\n",
			p.QuerySample, p.GenotypedSample, p.DiscCount, cExp(hwe), p.NumSites, nmatch)
		return err
	}
	_, err := fmt.Fprintf(w, "DCv2\t%s\t%s\t%s\t%s\t%d\t%d\n",
		p.QuerySample, p.GenotypedSample, cExp(p.DiscScore), cExp(hwe), p.NumSites, nmatch)
	return err
}

// selectReportRows applies --n-matches trimming, mirroring upstream
// report(). With NMatches == 0, explicit -p/-P pairs, or a cohort no
// bigger than |NMatches|, every pair is reported in input order.
// Otherwise the top |NMatches| matches per query sample are kept, sorted
// by average discordance ascending (or, when NMatches is negative, by
// average HWE probability descending).
func selectReportRows(r GtcheckResult, opts GtcheckOptions, crossCheck bool, qSamples, gSamples []string) []GtcheckPair {
	usePairs := opts.PairsSpec != "" || opts.PairsFile != ""
	ntop := opts.NMatches
	sortByHWE := ntop < 0
	if ntop < 0 {
		ntop = -ntop
	}

	// No trimming when explicit pairs, ntop == 0, or the cohort is no
	// bigger than ntop (matches upstream's `trim` computation).
	trim := ntop
	switch {
	case usePairs || ntop == 0:
		trim = 0
	case crossCheck:
		if len(qSamples) <= ntop {
			trim = 0
		}
	default:
		if len(gSamples) <= ntop {
			trim = 0
		}
	}
	if trim == 0 {
		return r.Pairs
	}

	if crossCheck {
		return crossCheckTopN(r.Pairs, qSamples, trim, sortByHWE)
	}
	return pairedTopN(r.Pairs, qSamples, trim, sortByHWE)
}

// pairedTopN keeps, for each query sample (in cohort order), the top-N
// genotyped samples by score, mirroring upstream's non-cross-check trim.
func pairedTopN(pairs []GtcheckPair, qSamples []string, trim int, sortByHWE bool) []GtcheckPair {
	groups := map[string][]GtcheckPair{}
	for _, p := range pairs {
		groups[p.QuerySample] = append(groups[p.QuerySample], p)
	}
	out := make([]GtcheckPair, 0, len(qSamples)*trim)
	for _, q := range qSamples {
		g := groups[q]
		sortByScore(g, sortByHWE)
		n := trim
		if n > len(g) {
			n = len(g)
		}
		out = append(out, g[:n]...)
	}
	return out
}

// crossCheckTopN reproduces upstream's cross-check trim: for each sample
// i (in cohort order) it considers ALL other samples (using the stored
// half-triangle pair regardless of orientation), sorts them by score,
// takes the top-N, and emits only the rows where i is the query (the
// larger-indexed sample), i.e. where i ranks after its partner.
func crossCheckTopN(pairs []GtcheckPair, qSamples []string, trim int, sortByHWE bool) []GtcheckPair {
	n := len(qSamples)
	idxOf := map[string]int{}
	for i, s := range qSamples {
		idxOf[s] = i
	}
	// byIJ[min*n+max] = the stored pair (query=max, gt=min).
	byIJ := make(map[int]GtcheckPair, len(pairs))
	for _, p := range pairs {
		qi, gi := idxOf[p.QuerySample], idxOf[p.GenotypedSample]
		lo, hi := gi, qi
		if lo > hi {
			lo, hi = hi, lo
		}
		byIJ[lo*n+hi] = p
	}

	out := make([]GtcheckPair, 0, n*trim)
	for i := 0; i < n; i++ {
		type cand struct {
			ism int
			p   GtcheckPair
		}
		cands := make([]cand, 0, n-1)
		for j := 0; j < n; j++ {
			if j == i {
				continue
			}
			lo, hi := i, j
			if lo > hi {
				lo, hi = hi, lo
			}
			p, ok := byIJ[lo*n+hi]
			if !ok {
				continue
			}
			cands = append(cands, cand{ism: j, p: p})
		}
		sort.SliceStable(cands, func(a, b int) bool {
			return rowSortKey(cands[a].p, sortByHWE) < rowSortKey(cands[b].p, sortByHWE)
		})
		// Upstream walks exactly the first `trim` sorted candidates and
		// emits those where i is the larger index (the query side).
		limit := trim
		if limit > len(cands) {
			limit = len(cands)
		}
		for k := 0; k < limit; k++ {
			if i <= cands[k].ism {
				continue
			}
			out = append(out, cands[k].p)
		}
	}
	return out
}

// sortByScore stably sorts pairs by upstream's --n-matches sort key.
func sortByScore(g []GtcheckPair, sortByHWE bool) {
	sort.SliceStable(g, func(i, j int) bool {
		return rowSortKey(g[i], sortByHWE) < rowSortKey(g[j], sortByHWE)
	})
}

// cExp formats f the way C's printf("%e") does: six fractional digits
// and a sign plus at-least-two-digit exponent. Go's "%e" verb already
// matches this representation, so this is a thin wrapper that documents
// the parity requirement.
func cExp(f float64) string {
	return fmt.Sprintf("%e", f)
}

// rowSortKey is the value upstream sorts on in the --n-matches path.
func rowSortKey(p GtcheckPair, sortByHWE bool) float64 {
	if sortByHWE {
		if p.NumMatching > 0 {
			return -p.AvgLogPHWE
		}
		return 0
	}
	if p.NumSites > 0 {
		if p.IsInteger {
			return float64(p.DiscCount) / float64(p.NumSites)
		}
		return p.DiscScore / float64(p.NumSites)
	}
	return 0
}
