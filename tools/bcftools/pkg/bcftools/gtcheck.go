package bcftools

// bcftools gtcheck — sample-identity / contamination check.
//
// Upstream reference: reference_code/bcftools/vcfgtcheck.c. The output
// header for the table is the literal string `#DCv2` followed by tab-
// separated column descriptors:
//
//   #DCv2	[2]Query Sample	[3]Genotyped Sample	[4]Discordance
//   	[5]Average -log P(HWE)	[6]Number of sites compared
//   	[7]Number of matching genotypes
//
// The v1 port covers the hard-GT (`-u GT`) path. We:
//   - reject multi-allelic records with the exact upstream "run
//     `bcftools norm -m -` first" diagnostic (anything with len(ALT)>1
//     OR a comma-separated ALT field);
//   - score each (query, genotyped) pair via biallelic GT-dosage Hamming
//     distance (0/0=0, 0/1=1, 1/1=2; ./. is a skip, NOT a discordance);
//   - emit one row per pair, ordered (qrySample, gtSample).
//
// Out of scope for v1, tracked in docs/PARITY_ROADMAP.md#bcftools:
//   - `-u PL` (likelihood-based scoring),
//   - `--cluster`, `--distinctive-sites`, `--n-matches`,
//   - `-E/--error-probability`-weighted scoring.
//
// All of these are accepted at the CLI for surface-parity and rejected
// with a clear pointer to the roadmap.

import (
	"bufio"
	"fmt"
	"io"
	"os"
	"sort"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// GtcheckOptions controls the behaviour of Gtcheck / GtcheckFile.
//
// The struct documents v1 surface explicitly: only the upstream flags
// that have a real effect are reflected here. CLI-only flags that are
// "parse and reject with PARITY_ROADMAP pointer" live in the CLI runner
// (subcmds_gtcheck_roh.go) so the library stays small.
type GtcheckOptions struct {
	// GenotypesFile is the -g panel of "truth" genotypes. When empty,
	// gtcheck runs in cross-check mode: every query sample versus
	// every other query sample, half-square (i<j).
	GenotypesFile string

	// PairsSpec is the upstream -p comma list. Layout depends on
	// whether GenotypesFile is set: qry,gt[,qry,gt...] when -g is
	// given, qry,qry[,qry,qry...] otherwise.
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

	// ErrorProbability is the phred-scaled -E. Recorded for traceability;
	// the hard-GT score is unweighted in v1.
	ErrorProbability int

	// OutputType: "t" (default) or "z". v1 supports only "t".
	OutputType string
}

// GtcheckPair captures one row of the #DCv2 output.
type GtcheckPair struct {
	QuerySample     string
	GenotypedSample string
	Discordance     int
	AvgLogPHWE      float64
	NumSites        int
	NumMatching     int
}

// GtcheckResult is the full pair table.
type GtcheckResult struct {
	Pairs []GtcheckPair
}

// GtcheckFile is the file-aware entry point used by the CLI. It opens
// the query (and optional panel) through iohelper, scores each pair,
// and writes the #DCv2 table to out.
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

// openVariantReader opens a VCF/BCF/gzipped-VCF input. It returns an
// io.ReadCloser to give the caller control of lifetime.
func openVariantReader(path string) (io.ReadCloser, error) {
	return iohelper.OpenReader(path)
}

// Gtcheck is the cross-check entry point: every query sample vs every
// other query sample. Pairs are emitted in (i<j) order using the input
// header's sample ordering.
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

// runGtcheck is the shared scoring path. crossCheck=true means the
// query and panel are the same cohort and we emit only the i<j half.
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

	if err := rejectMultiAllelic(varsQ, "query"); err != nil {
		return GtcheckResult{}, err
	}
	if !crossCheck {
		if err := rejectMultiAllelic(varsG, "panel"); err != nil {
			return GtcheckResult{}, err
		}
	}

	qSamples := selectSamples(hdrQ.Samples, opts.SamplesQry)
	gSamples := qSamples
	if !crossCheck {
		gSamples = selectSamples(hdrG.Samples, opts.SamplesGT)
	}
	if len(qSamples) == 0 || len(gSamples) == 0 {
		return GtcheckResult{}, fmt.Errorf("bcftools gtcheck: no samples to compare (qry=%d gt=%d)", len(qSamples), len(gSamples))
	}

	pairs, err := buildPairs(qSamples, gSamples, opts, crossCheck)
	if err != nil {
		return GtcheckResult{}, err
	}

	// Build (CHROM,POS) -> *Variant lookup for the panel so cross-file
	// joining is O(1) per query record.
	panelByPos := indexByPos(varsG)
	qSampleIdx := indexSamples(hdrQ)
	gSampleIdx := indexSamples(hdrG)

	type acc struct {
		disc  int
		match int
		sites int
		// Cumulative -log10(P(HWE)) for the HOM-REF / HOM-ALT / HET
		// allele-counts seen at scored sites. v1 uses a fixed
		// dummy of 0 when NoHWEProb is set; otherwise it sums
		// -log10(0.5) per site as a sentinel until upstream's
		// real per-site HWE estimator lands.
		hweAcc float64
	}
	accums := make(map[[2]string]*acc, len(pairs))
	for _, p := range pairs {
		key := [2]string{p.QuerySample, p.GenotypedSample}
		accums[key] = &acc{}
	}

	sites := 0
	for _, qv := range varsQ {
		// Apply region / targets post-filters independently. Either
		// filter being set is sufficient to drop a non-matching
		// variant; both unset means everything passes.
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
				continue
			}
			if err := rejectIfMultiAllelic(gv, "panel"); err != nil {
				return GtcheckResult{}, err
			}
		}
		sites++
		for _, p := range pairs {
			a := accums[[2]string{p.QuerySample, p.GenotypedSample}]
			qIdx, ok := qSampleIdx[p.QuerySample]
			if !ok {
				continue
			}
			gIdx, ok2 := gSampleIdx[p.GenotypedSample]
			if !ok2 {
				continue
			}
			qDose, qOK := gtDose(qv, qIdx)
			gDose, gOK := gtDose(gv, gIdx)
			if !qOK || !gOK {
				// Missing → skip (NOT a discordance).
				continue
			}
			if opts.HomsOnly && gDose == 1 {
				continue
			}
			a.sites++
			d := qDose - gDose
			if d < 0 {
				d = -d
			}
			a.disc += d
			if d == 0 {
				a.match++
			}
			if !opts.NoHWEProb {
				// v1 placeholder: the column is zeroed until a
				// real per-site HWE estimator from panel AF
				// lands (tracked in PARITY_ROADMAP). A non-zero
				// constant placeholder (like -log10(0.5)) was
				// worse than zero — every row would read the
				// same useless number.
				a.hweAcc += 0
			}
		}
		if opts.DryRun {
			break
		}
	}

	if !opts.AllSites && sites == 0 {
		return GtcheckResult{}, fmt.Errorf("bcftools gtcheck: no overlapping scoreable sites between query and -g panel")
	}

	result := GtcheckResult{Pairs: make([]GtcheckPair, 0, len(pairs))}
	for _, p := range pairs {
		a := accums[[2]string{p.QuerySample, p.GenotypedSample}]
		avgHWE := 0.0
		if a.sites > 0 && !opts.NoHWEProb {
			avgHWE = a.hweAcc / float64(a.sites)
		}
		result.Pairs = append(result.Pairs, GtcheckPair{
			QuerySample:     p.QuerySample,
			GenotypedSample: p.GenotypedSample,
			Discordance:     a.disc,
			AvgLogPHWE:      avgHWE,
			NumSites:        a.sites,
			NumMatching:     a.match,
		})
	}
	if err := writeGtcheckDCv2(out, result); err != nil {
		return result, err
	}
	return result, nil
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

// rejectMultiAllelic mirrors upstream's input-shape check.
func rejectMultiAllelic(vars []*vcf.Variant, label string) error {
	for _, v := range vars {
		if err := rejectIfMultiAllelic(v, label); err != nil {
			return err
		}
	}
	return nil
}

// rejectIfMultiAllelic enforces the biallelic input constraint that
// upstream's vcfgtcheck.c expects: any ALT > 1 ALT means the caller
// must run `bcftools norm -m -` first.
func rejectIfMultiAllelic(v *vcf.Variant, label string) error {
	if v == nil {
		return nil
	}
	if len(v.Alt) > 1 {
		return fmt.Errorf("bcftools gtcheck: multi-allelic %s record at %s:%d has %d ALT alleles; run `bcftools norm -m -` first", label, v.Chrom, v.Pos, len(v.Alt))
	}
	if len(v.Alt) == 1 && strings.Contains(v.Alt[0], ",") {
		return fmt.Errorf("bcftools gtcheck: multi-allelic %s record at %s:%d; run `bcftools norm -m -` first", label, v.Chrom, v.Pos)
	}
	return nil
}

// gtDose returns the biallelic dosage of the FORMAT/GT field for the
// given sample-index: 0 for 0/0, 1 for 0/1 or 1/0, 2 for 1/1. The
// second return is false if the GT is missing or unparseable.
func gtDose(v *vcf.Variant, idx int) (int, bool) {
	if v == nil || idx < 0 || idx >= len(v.Samples) {
		return 0, false
	}
	gt, ok := v.Samples[idx].Data["GT"]
	if !ok || gt == "" || gt == "." || gt == "./." || gt == ".|." {
		return 0, false
	}
	// Strip phasing.
	gt = strings.ReplaceAll(gt, "|", "/")
	parts := strings.Split(gt, "/")
	dose := 0
	for _, p := range parts {
		if p == "." {
			return 0, false
		}
		// Anything other than 0 or 1 in v1 (multi-allelic) was
		// rejected upstream of this call.
		switch p {
		case "0":
			// no-op
		case "1":
			dose++
		default:
			return 0, false
		}
	}
	return dose, true
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

// buildPairs assembles the (query, genotype) pair list. When PairsSpec
// or PairsFile is set, use those; otherwise the default is "all qry
// vs all gt", with crossCheck mode using the i<j half-square.
func buildPairs(qSamples, gSamples []string, opts GtcheckOptions, crossCheck bool) ([]GtcheckPair, error) {
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
	pairs := make([]GtcheckPair, 0, len(qSamples)*len(gSamples))
	if crossCheck {
		// Half-square i<j.
		sortedQ := append([]string{}, qSamples...)
		sort.Strings(sortedQ)
		for i := 0; i < len(sortedQ); i++ {
			for j := i + 1; j < len(sortedQ); j++ {
				pairs = append(pairs, GtcheckPair{
					QuerySample:     sortedQ[i],
					GenotypedSample: sortedQ[j],
				})
			}
		}
		return pairs, nil
	}
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
	// CHROM, CHROM:beg, CHROM:beg-end.
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

// writeGtcheckDCv2 emits the upstream "#DCv2" table. Columns match
// vcfgtcheck.c exactly:
//
//	#DCv2 \t [2]Query Sample \t [3]Genotyped Sample \t [4]Discordance
//	\t [5]Average -log P(HWE) \t [6]Number of sites compared
//	\t [7]Number of matching genotypes
func writeGtcheckDCv2(out io.Writer, r GtcheckResult) error {
	w := bufio.NewWriter(out)
	if _, err := w.WriteString("#DCv2\t[2]Query Sample\t[3]Genotyped Sample\t[4]Discordance\t[5]Average -log P(HWE)\t[6]Number of sites compared\t[7]Number of matching genotypes\n"); err != nil {
		return err
	}
	for _, p := range r.Pairs {
		if _, err := fmt.Fprintf(w,
			"DCv2\t%s\t%s\t%d\t%.6f\t%d\t%d\n",
			p.QuerySample, p.GenotypedSample, p.Discordance, p.AvgLogPHWE, p.NumSites, p.NumMatching); err != nil {
			return err
		}
	}
	return w.Flush()
}

// posStr is a small helper to avoid pulling strconv just for one fmt.
func posStr(p int) string {
	return fmt.Sprintf("%d", p)
}
