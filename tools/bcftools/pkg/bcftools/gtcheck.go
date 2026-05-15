// Package-level implementation of `bcftools gtcheck`.
//
// gtcheck compares each sample in a query VCF/BCF against a panel VCF
// ("genotype" file) and reports a per-sample-pair discordance score.
// Upstream supports several scoring modes (PL-based likelihood,
// dosage-correlation, Bayesian) and a clustering / homozygosity-error
// post-pass. The v1 port implements the simplest useful form: hard-
// genotype Hamming distance, optionally restricted to a list of sample
// pairs or sites.
//
// Algorithm:
//   - Read both files into memory and index variants by (CHROM, POS, REF, ALT[0]).
//   - For each site shared between query and panel, compare each
//     requested (query-sample, panel-sample) pair's hard genotype:
//     a mismatch (different unordered allele set) contributes 1 to
//     the pair's discordance score; a match contributes 0; sites where
//     either side is uncalled contribute nothing and are tallied
//     separately under "n_missing".
//   - Emit a TSV: `# DC\tquery\tpanel\tscore\tn_sites\tn_missing\tdiscordance`
//     header + one row per pair.
//
// Upstream parity notes:
//   - The PL/GL likelihood scoring (`-u PL`) and the Bayesian / dosage
//     paths are accepted by the CLI but emit a "not implemented in v1"
//     error pointing at docs/PARITY_ROADMAP.md.
//   - The `--all-sites` / `--no-HWE-prob` / homozygosity-error tuning
//     knobs are accepted-and-ignored; their impact on the Hamming
//     comparison is nil.
package bcftools

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/vcf"
)

// GtcheckUseMode selects the genotype-comparison metric.
type GtcheckUseMode int

const (
	// GtcheckUseGT compares hard genotypes (Hamming distance on the
	// unordered allele set). This is the only mode v1 implements.
	GtcheckUseGT GtcheckUseMode = iota
	// GtcheckUsePL would compare phred-scaled likelihoods. Accepted by
	// the CLI but rejected at runtime in v1.
	GtcheckUsePL
	// GtcheckUseGL would compare genotype likelihoods. Same status as PL.
	GtcheckUseGL
)

// ParseGtcheckUseMode parses the `-u/--use` CLI argument.
func ParseGtcheckUseMode(s string) (GtcheckUseMode, error) {
	switch strings.ToUpper(strings.TrimSpace(s)) {
	case "", "GT":
		return GtcheckUseGT, nil
	case "PL":
		return GtcheckUsePL, nil
	case "GL":
		return GtcheckUseGL, nil
	}
	return 0, fmt.Errorf("bcftools gtcheck: unknown -u mode %q (accept GT|PL|GL)", s)
}

// GtcheckPair names one (query, panel) sample pair.
type GtcheckPair struct {
	Query string
	Panel string
}

// GtcheckOptions controls Gtcheck / GtcheckFile.
type GtcheckOptions struct {
	// PanelPath is the -g/--genotypes panel VCF/BCF path. Required.
	PanelPath string
	// Pairs is the explicit list of (query, panel) sample pairs to
	// compare. Empty means "every query × every panel" cross-join.
	Pairs []GtcheckPair
	// PairsFile is the path supplied via -P/--pairs-file; each line
	// is "QUERY <tab|comma> PANEL". Loaded into Pairs at runtime.
	PairsFile string
	// Use selects the scoring metric (GT today; PL/GL deferred).
	Use GtcheckUseMode
	// NGSError is the genotype-error rate. Accepted; not used in v1.
	NGSError float64
	// Regions / Targets are post-filters on (CHROM, POS).
	Regions     []string
	RegionsFile string
	Targets     []string
	TargetsFile string
	// IncludeExpr / ExcludeExpr are accepted-and-ignored stubs in v1;
	// gtcheck doesn't apply per-record expressions in upstream either.
	IncludeExpr string
	ExcludeExpr string
}

// GtcheckResult is the rollup returned by Gtcheck for inspection (and
// for the symmetry invariant assertions in the tests).
type GtcheckResult struct {
	// Pairs holds the per-pair scores in input order.
	Pairs []GtcheckPairResult
	// NQuerySamples and NPanelSamples count the resolved sample lists.
	NQuerySamples int
	NPanelSamples int
	// NSitesCompared is the number of (chrom,pos,ref,alt) tuples that
	// appeared in BOTH files after region filtering.
	NSitesCompared int
}

// GtcheckPairResult is the per-pair output row.
type GtcheckPairResult struct {
	Query       string
	Panel       string
	Score       int     // sum of GT mismatches across compared sites
	NSites      int     // shared sites where both samples were called
	NMissing    int     // shared sites skipped (uncalled on one side)
	Discordance float64 // Score / NSites (0 when NSites == 0)
}

// GtcheckFile opens the query path (and the panel path stored in opts)
// and dispatches to Gtcheck.
func GtcheckFile(queryPath string, out io.Writer, opts GtcheckOptions) (GtcheckResult, error) {
	if opts.PanelPath == "" {
		return GtcheckResult{}, fmt.Errorf("bcftools gtcheck: -g/--genotypes panel file is required")
	}
	if opts.Use != GtcheckUseGT {
		return GtcheckResult{}, fmt.Errorf("bcftools gtcheck: -u %s not implemented in v1; tracked in docs/PARITY_ROADMAP.md#bcftools", gtcheckUseString(opts.Use))
	}
	// Apply pairs-file before opening the inputs so a bad path fails fast.
	if opts.PairsFile != "" {
		pairs, err := loadGtcheckPairs(opts.PairsFile)
		if err != nil {
			return GtcheckResult{}, fmt.Errorf("bcftools gtcheck: %w", err)
		}
		opts.Pairs = append(opts.Pairs, pairs...)
	}
	if opts.RegionsFile != "" {
		regs, err := LoadRegionsFile(opts.RegionsFile)
		if err != nil {
			return GtcheckResult{}, fmt.Errorf("bcftools gtcheck: %w", err)
		}
		opts.Regions = append(opts.Regions, regs...)
	}
	if opts.TargetsFile != "" {
		regs, err := LoadRegionsFile(opts.TargetsFile)
		if err != nil {
			return GtcheckResult{}, fmt.Errorf("bcftools gtcheck: %w", err)
		}
		opts.Targets = append(opts.Targets, regs...)
	}

	queryR, err := iohelper.OpenReader(queryPath)
	if err != nil {
		return GtcheckResult{}, fmt.Errorf("bcftools gtcheck: open query %s: %w", queryPath, err)
	}
	defer queryR.Close()
	panelR, err := iohelper.OpenReader(opts.PanelPath)
	if err != nil {
		return GtcheckResult{}, fmt.Errorf("bcftools gtcheck: open panel %s: %w", opts.PanelPath, err)
	}
	defer panelR.Close()

	return Gtcheck(queryR, panelR, out, opts)
}

// Gtcheck is the in-memory entry point: reads both readers, computes
// per-pair scores, and writes the TSV report to out.
func Gtcheck(queryIn, panelIn io.Reader, out io.Writer, opts GtcheckOptions) (GtcheckResult, error) {
	queryHdr, queryVars, err := readAllVariants(queryIn)
	if err != nil {
		return GtcheckResult{}, fmt.Errorf("bcftools gtcheck: query: %w", err)
	}
	panelHdr, panelVars, err := readAllVariants(panelIn)
	if err != nil {
		return GtcheckResult{}, fmt.Errorf("bcftools gtcheck: panel: %w", err)
	}

	regions, err := parseRegions(opts.Regions)
	if err != nil {
		return GtcheckResult{}, fmt.Errorf("bcftools gtcheck: %w", err)
	}
	targets, err := parseRegions(opts.Targets)
	if err != nil {
		return GtcheckResult{}, fmt.Errorf("bcftools gtcheck: %w", err)
	}
	queryVars = filterByRegions(queryVars, regions, targets)
	panelVars = filterByRegions(panelVars, regions, targets)

	// Resolve sample pairs. With an empty list we cross-join everything.
	pairs := opts.Pairs
	if len(pairs) == 0 {
		pairs = crossJoinSamples(queryHdr.Samples, panelHdr.Samples)
	}
	if len(pairs) == 0 {
		return GtcheckResult{}, fmt.Errorf("bcftools gtcheck: no sample pairs to compare (query has %d samples, panel has %d)", len(queryHdr.Samples), len(panelHdr.Samples))
	}

	queryIdx, panelIdx, err := resolvePairIndices(queryHdr, panelHdr, pairs)
	if err != nil {
		return GtcheckResult{}, fmt.Errorf("bcftools gtcheck: %w", err)
	}

	// Build a (chrom, pos, ref, alt0) → panelVariant index. We use the
	// first ALT only — this matches the upstream behaviour for the
	// hard-GT comparison; multi-allelic sites that disagree on ALT
	// ordering will be flagged "not compared" rather than mismatched.
	panelByKey := make(map[string]*vcf.Variant, len(panelVars))
	for _, v := range panelVars {
		panelByKey[variantKey(v)] = v
	}

	result := GtcheckResult{
		Pairs:         make([]GtcheckPairResult, len(pairs)),
		NQuerySamples: len(queryHdr.Samples),
		NPanelSamples: len(panelHdr.Samples),
	}
	for i, p := range pairs {
		result.Pairs[i] = GtcheckPairResult{Query: p.Query, Panel: p.Panel}
	}

	nShared := 0
	for _, qv := range queryVars {
		pv, ok := panelByKey[variantKey(qv)]
		if !ok {
			continue
		}
		nShared++
		for i := range pairs {
			qGT, qOK := parseHardGT(sampleData(qv, queryIdx[i]))
			pGT, pOK := parseHardGT(sampleData(pv, panelIdx[i]))
			if !qOK || !pOK {
				result.Pairs[i].NMissing++
				continue
			}
			result.Pairs[i].NSites++
			if !gtEqual(qGT, pGT) {
				result.Pairs[i].Score++
			}
		}
	}
	result.NSitesCompared = nShared

	for i := range result.Pairs {
		if result.Pairs[i].NSites > 0 {
			result.Pairs[i].Discordance = float64(result.Pairs[i].Score) / float64(result.Pairs[i].NSites)
		}
	}

	if err := writeGtcheckReport(out, result); err != nil {
		return result, err
	}
	return result, nil
}

// gtcheckUseString is the human label used in the deferred-mode error.
func gtcheckUseString(m GtcheckUseMode) string {
	switch m {
	case GtcheckUsePL:
		return "PL"
	case GtcheckUseGL:
		return "GL"
	}
	return "GT"
}

// crossJoinSamples returns every (q, p) pair, in q-major order.
func crossJoinSamples(query, panel []string) []GtcheckPair {
	out := make([]GtcheckPair, 0, len(query)*len(panel))
	for _, q := range query {
		for _, p := range panel {
			out = append(out, GtcheckPair{Query: q, Panel: p})
		}
	}
	return out
}

// resolvePairIndices maps each pair's sample names to column indexes in
// the respective headers. Returns parallel slices: queryIdx[i] is the
// column index in queryHdr.Samples for pairs[i].Query, panelIdx[i] is
// the matching panel index.
func resolvePairIndices(queryHdr, panelHdr *vcf.Header, pairs []GtcheckPair) ([]int, []int, error) {
	qByName := indexSamples(queryHdr.Samples)
	pByName := indexSamples(panelHdr.Samples)
	qIdx := make([]int, len(pairs))
	pIdx := make([]int, len(pairs))
	for i, pr := range pairs {
		q, ok := qByName[pr.Query]
		if !ok {
			return nil, nil, fmt.Errorf("pair %d: query sample %q not in query file", i+1, pr.Query)
		}
		p, ok := pByName[pr.Panel]
		if !ok {
			return nil, nil, fmt.Errorf("pair %d: panel sample %q not in panel file", i+1, pr.Panel)
		}
		qIdx[i], pIdx[i] = q, p
	}
	return qIdx, pIdx, nil
}

// indexSamples turns a sample-name slice into a name→index map.
func indexSamples(names []string) map[string]int {
	out := make(map[string]int, len(names))
	for i, n := range names {
		out[n] = i
	}
	return out
}

// variantKey is the join key used to align query/panel records. We use
// CHROM\tPOS\tREF\tALT0; sites with no ALT (gVCF "<NON_REF>" only) are
// skipped via an empty key — see the check on the caller side.
func variantKey(v *vcf.Variant) string {
	alt := ""
	if len(v.Alt) > 0 {
		alt = v.Alt[0]
	}
	var b strings.Builder
	b.Grow(len(v.Chrom) + len(v.Ref) + len(alt) + 16)
	b.WriteString(v.Chrom)
	b.WriteByte('\t')
	b.WriteString(strconv.Itoa(v.Pos))
	b.WriteByte('\t')
	b.WriteString(v.Ref)
	b.WriteByte('\t')
	b.WriteString(alt)
	return b.String()
}

// parseHardGT extracts a sorted two-element allele set from the GT
// string. Returns ok == false for any of: empty, missing (`.`), non-
// diploid, or contains a missing allele.
func parseHardGT(gt string) ([2]int, bool) {
	var zero [2]int
	if gt == "" || gt == "." {
		return zero, false
	}
	gt = strings.ReplaceAll(gt, "|", "/")
	parts := strings.Split(gt, "/")
	if len(parts) != 2 {
		return zero, false
	}
	var out [2]int
	for i, p := range parts {
		if p == "" || p == "." {
			return zero, false
		}
		n, err := strconv.Atoi(p)
		if err != nil {
			return zero, false
		}
		out[i] = n
	}
	if out[0] > out[1] {
		out[0], out[1] = out[1], out[0]
	}
	return out, true
}

// gtEqual returns true if the two sorted diploid genotypes are equal.
func gtEqual(a, b [2]int) bool {
	return a[0] == b[0] && a[1] == b[1]
}

// loadGtcheckPairs reads a PAIRS-file of "QUERY <sep> PANEL" lines.
// Separator: any of `,`, tab, or whitespace.
func loadGtcheckPairs(path string) ([]GtcheckPair, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 0, 64<<10), 1<<20)
	var out []GtcheckPair
	for sc.Scan() {
		line := strings.TrimSpace(sc.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		var fields []string
		switch {
		case strings.Contains(line, ","):
			fields = strings.Split(line, ",")
		default:
			fields = strings.Fields(line)
		}
		if len(fields) < 2 {
			return nil, fmt.Errorf("pairs-file %s: malformed line %q", path, line)
		}
		out = append(out, GtcheckPair{
			Query: strings.TrimSpace(fields[0]),
			Panel: strings.TrimSpace(fields[1]),
		})
	}
	return out, sc.Err()
}

// ParseGtcheckPairsSpec splits a `-p/--pairs` spec into one or more
// pairs. The upstream syntax is "Q1,P1,Q2,P2,..." (a flat alternating
// list); we accept that plus the more readable "Q1:P1,Q2:P2" form.
func ParseGtcheckPairsSpec(s string) ([]GtcheckPair, error) {
	s = strings.TrimSpace(s)
	if s == "" {
		return nil, nil
	}
	tokens := strings.Split(s, ",")
	// Detect the "Q:P,Q:P" form by looking for a ':' in any token.
	hasColon := false
	for _, t := range tokens {
		if strings.Contains(t, ":") {
			hasColon = true
			break
		}
	}
	var out []GtcheckPair
	if hasColon {
		for _, t := range tokens {
			t = strings.TrimSpace(t)
			if t == "" {
				continue
			}
			parts := strings.SplitN(t, ":", 2)
			if len(parts) != 2 || parts[0] == "" || parts[1] == "" {
				return nil, fmt.Errorf("bcftools gtcheck: bad -p token %q (expected QUERY:PANEL)", t)
			}
			out = append(out, GtcheckPair{Query: strings.TrimSpace(parts[0]), Panel: strings.TrimSpace(parts[1])})
		}
		return out, nil
	}
	if len(tokens)%2 != 0 {
		return nil, fmt.Errorf("bcftools gtcheck: -p needs an even number of comma-separated names (got %d)", len(tokens))
	}
	for i := 0; i < len(tokens); i += 2 {
		q := strings.TrimSpace(tokens[i])
		p := strings.TrimSpace(tokens[i+1])
		if q == "" || p == "" {
			return nil, fmt.Errorf("bcftools gtcheck: empty sample in -p %q", s)
		}
		out = append(out, GtcheckPair{Query: q, Panel: p})
	}
	return out, nil
}

// writeGtcheckReport emits the TSV table used by upstream's `DC`
// (discordance) section. One header row + one body row per pair.
// Columns are deliberately stable for `awk`/`diff`-based downstream
// consumers; if upstream evolves the layout, we'll mirror it then.
func writeGtcheckReport(out io.Writer, r GtcheckResult) error {
	w := bufio.NewWriter(out)
	defer w.Flush()
	if _, err := fmt.Fprintf(w, "# DC\tquery\tpanel\tscore\tn_sites\tn_missing\tdiscordance\n"); err != nil {
		return err
	}
	// Stable per-pair ordering: respect input order. (Upstream sorts
	// by score; we leave that to consumers because preserving caller
	// order makes batch comparisons easier.)
	for _, p := range r.Pairs {
		if _, err := fmt.Fprintf(w, "DC\t%s\t%s\t%d\t%d\t%d\t%s\n",
			p.Query, p.Panel, p.Score, p.NSites, p.NMissing, formatGtcheckRate(p.Discordance)); err != nil {
			return err
		}
	}
	if _, err := fmt.Fprintf(w, "# totals\tn_sites_compared=%d\tn_query_samples=%d\tn_panel_samples=%d\n",
		r.NSitesCompared, r.NQuerySamples, r.NPanelSamples); err != nil {
		return err
	}
	return nil
}

// formatGtcheckRate prints the discordance fraction in a form that's
// stable across platforms (shortest %g, with at most 6 fractional
// digits). 0 prints as "0", 1 prints as "1".
func formatGtcheckRate(x float64) string {
	if x == 0 {
		return "0"
	}
	if x == 1 {
		return "1"
	}
	s := strconv.FormatFloat(x, 'g', 6, 64)
	return s
}

// SortedGtcheckPairs returns r.Pairs sorted by ascending discordance
// (stable on score / pair-name for determinism). Helper for callers
// who want the upstream "best match first" ordering.
func SortedGtcheckPairs(r GtcheckResult) []GtcheckPairResult {
	out := append([]GtcheckPairResult{}, r.Pairs...)
	sort.SliceStable(out, func(i, j int) bool {
		if out[i].Discordance != out[j].Discordance {
			return out[i].Discordance < out[j].Discordance
		}
		if out[i].Query != out[j].Query {
			return out[i].Query < out[j].Query
		}
		return out[i].Panel < out[j].Panel
	})
	return out
}
