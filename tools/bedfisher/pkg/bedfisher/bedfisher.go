// Package bedfisher computes Fisher's exact test of overlap enrichment
// between two BED files over a genome, mirroring `bedtools fisher`.
//
// The algorithm:
//
//  1. Read A and B (optionally merging A with -m); compute raw counts and
//     summed lengths (queryUnion, dbUnion).
//  2. Sweep A against B (chromosome by chromosome) to count overlap pairs.
//     Honour the -f / -F (fraction filters), -r (reciprocal), -s / -S
//     (strand) flags exactly as upstream does for the overlap predicate.
//  3. Compute the n22 (genome-bg) count via the same heuristic as upstream:
//     bMean = (1 + queryUnion/queryCounts) + (1 + dbUnion/dbCounts);
//     n22_full = max(n11+n12+n21, genomeSize/bMean).
//  4. Run Fisher's exact test (two-tailed) using the same hypergeometric
//     accumulator as upstream's kfunc.cpp, producing left, right, two-tail
//     p-values and the odds-ratio.
//
// The output format matches upstream byte-for-byte for the cases we test
// against (see tools/bedfisher/testdata/parity/).
package bedfisher

import (
	"bufio"
	"errors"
	"fmt"
	"io"
	"math"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bed"
)

// Options configures the Fisher computation.
type Options struct {
	// FractionA: minimum fraction of A required to overlap B (-f).
	FractionA float64
	// FractionB: minimum fraction of B required to overlap A (-F).
	FractionB float64
	// Reciprocal: require -f to apply to both A and B (-r / --reciprocal).
	Reciprocal bool
	// SameStrand limits overlaps to same-strand pairs (-s).
	SameStrand bool
	// OppositeStrand limits overlaps to opposite-strand pairs (-S).
	OppositeStrand bool
	// MergeInputs pre-merges overlapping records before the sweep (-m). Upstream
	// applies this to BOTH the query (A) and database (B) files — its
	// FileRecordMergeMgr merges every input file when -m is given, not just A.
	MergeInputs bool
}

// Result holds the contingency table and p-values.
type Result struct {
	QueryCount  int   // number of intervals in A (after optional merge)
	DBCount     int   // number of intervals in B
	OverlapPair int   // n11: overlapping pairs
	N22Full     int64 // estimated total possible intervals
	N12         int64
	N21         int64
	N22         int64

	Left    float64
	Right   float64
	TwoTail float64
	Ratio   float64 // NaN means -nan in output (typical 0/0 cases)
}

// Run reads A, B and the genome file, runs the test and writes the upstream
// report to w.
func Run(a, b, g io.Reader, w io.Writer, opts Options) (*Result, error) {
	if opts.SameStrand && opts.OppositeStrand {
		return nil, errors.New("-s and -S are mutually exclusive")
	}
	if opts.FractionA < 0 || opts.FractionA > 1 {
		return nil, fmt.Errorf("-f must be in [0,1]")
	}
	if opts.FractionB < 0 || opts.FractionB > 1 {
		return nil, fmt.Errorf("-F must be in [0,1]")
	}

	genomeSize, err := loadGenomeSize(g)
	if err != nil {
		return nil, fmt.Errorf("reading -g genome: %w", err)
	}

	aRecs, err := readAllBED(a)
	if err != nil {
		return nil, fmt.Errorf("reading A: %w", err)
	}
	bRecs, err := readAllBED(b)
	if err != nil {
		return nil, fmt.Errorf("reading B: %w", err)
	}
	if opts.MergeInputs {
		// -m merges every input file (both A and B), matching upstream's
		// FileRecordMergeMgr which is enabled for all files when -m is set.
		aRecs = mergeBED(aRecs)
		bRecs = mergeBED(bRecs)
	}

	// Sum raw counts and total lengths.
	qCount := len(aRecs)
	dCount := len(bRecs)
	var qUnion, dUnion int64
	for _, r := range aRecs {
		qUnion += int64(r.ChromEnd - r.ChromStart)
	}
	for _, r := range bRecs {
		dUnion += int64(r.ChromEnd - r.ChromStart)
	}

	// Count overlap pairs.
	overlapPairs := countOverlapPairs(aRecs, bRecs, opts)

	// Heuristic for the n22 background size, matching upstream:
	//   qMean = 1 + qUnion/qCount; dMean = 1 + dUnion/dCount; bMean = qMean+dMean.
	//   n22_full = max(n11+n12+n21, genomeSize/bMean).
	var n11 = int64(overlapPairs)
	var n12 = int64(qCount) - n11
	if n12 < 0 {
		n12 = 0
	}
	var n21 = int64(dCount) - n11
	if n21 < 0 {
		n21 = 0
	}

	var n22Full int64
	if qCount > 0 && dCount > 0 {
		qMean := 1.0 + float64(qUnion)/float64(qCount)
		dMean := 1.0 + float64(dUnion)/float64(dCount)
		bMean := qMean + dMean
		est := int64(float64(genomeSize) / bMean) // truncation matches C long-long cast
		filled := n11 + n12 + n21
		if est < filled {
			est = filled
		}
		n22Full = est
	} else {
		n22Full = n11 + n12 + n21
	}
	n22 := n22Full - n12 - n21 - n11
	if n22 < 0 {
		n22 = 0
	}

	left, right, two := ktFisherExact(n11, n12, n21, n22)
	ratio := math.NaN()
	if n12 != 0 && n21 != 0 && n22 != 0 {
		ratio = (float64(n11) / float64(n12)) / (float64(n21) / float64(n22))
	} else if n12 != 0 && n22 != 0 && n21 == 0 {
		ratio = math.Inf(1)
	}

	res := &Result{
		QueryCount: qCount, DBCount: dCount, OverlapPair: overlapPairs,
		N22Full: n22Full, N12: n12, N21: n21, N22: n22,
		Left: left, Right: right, TwoTail: two, Ratio: ratio,
	}

	if err := writeReport(w, res); err != nil {
		return nil, err
	}
	return res, nil
}

func writeReport(w io.Writer, r *Result) error {
	bw := bufio.NewWriter(w)
	defer bw.Flush()
	fmt.Fprintf(bw, "# Number of query intervals: %d\n", r.QueryCount)
	fmt.Fprintf(bw, "# Number of db intervals: %d\n", r.DBCount)
	fmt.Fprintf(bw, "# Number of overlaps: %d\n", r.OverlapPair)
	fmt.Fprintf(bw, "# Number of possible intervals (estimated): %d\n", r.N22Full)
	fmt.Fprintf(bw, "# phyper(%d - 1, %d, %d - %d, %d, lower.tail=F)\n",
		r.OverlapPair, r.QueryCount, r.N22Full, r.QueryCount, r.DBCount)
	fmt.Fprintln(bw, "# Contingency Table Of Counts")
	fmt.Fprintln(bw, "#_________________________________________")
	// Upstream uses printf("%-12s", "..."). The header row left-pads " in -b" /
	// "not in -b" within a 12-wide field; we replicate that.
	fmt.Fprintf(bw, "#           | %-12s | %-12s |\n", " in -b", "not in -b")
	fmt.Fprintf(bw, "#     in -a | %-12d | %-12d |\n", r.OverlapPair, r.N12)
	fmt.Fprintf(bw, "# not in -a | %-12d | %-12d |\n", r.N21, r.N22)
	fmt.Fprintln(bw, "#_________________________________________")
	fmt.Fprintln(bw, "# p-values for fisher's exact test")
	fmt.Fprintln(bw, "left\tright\ttwo-tail\tratio")
	leftS := format5g(r.Left)
	rightS := format5g(r.Right)
	twoS := format5g(r.TwoTail)
	if math.IsNaN(r.Ratio) {
		fmt.Fprintf(bw, "%s\t%s\t%s\t-nan\n", leftS, rightS, twoS)
	} else if math.IsInf(r.Ratio, 1) {
		fmt.Fprintf(bw, "%s\t%s\t%s\tinf\n", leftS, rightS, twoS)
	} else if math.IsInf(r.Ratio, -1) {
		fmt.Fprintf(bw, "%s\t%s\t%s\t-inf\n", leftS, rightS, twoS)
	} else {
		fmt.Fprintf(bw, "%s\t%s\t%s\t%.3f\n", leftS, rightS, twoS, r.Ratio)
	}
	return nil
}

// format5g renders a probability using `%.5g`, matching upstream's printf
// format for left/right/two-tail values.
func format5g(v float64) string {
	if math.IsNaN(v) {
		return "nan"
	}
	if math.IsInf(v, 1) {
		return "inf"
	}
	if math.IsInf(v, -1) {
		return "-inf"
	}
	return strconv.FormatFloat(v, 'g', 5, 64)
}

// countOverlapPairs returns the number of (A, B) record pairs that overlap
// under the filter options. Mirrors upstream's chromsweep accounting: each
// B that overlaps an A contributes one to the count.
func countOverlapPairs(aRecs, bRecs []*bed.Record, opts Options) int {
	// Build a per-chrom sorted slice of B records.
	byChrom := map[string][]*bed.Record{}
	for _, r := range bRecs {
		byChrom[r.Chrom] = append(byChrom[r.Chrom], r)
	}
	for k := range byChrom {
		sort.Slice(byChrom[k], func(i, j int) bool {
			return byChrom[k][i].ChromStart < byChrom[k][j].ChromStart
		})
	}

	// B records are sorted by ChromStart. A naive binary search on ChromEnd is
	// INVALID here because ChromEnd is not monotonic over a start-sorted slice —
	// a long B can begin before A yet extend well past A.Start and must still be
	// counted. (That bug under-counted overlaps, skewing the contingency table.)
	// Instead we binary-search the exact, monotonic upper bound on ChromStart
	// (the first B that starts at or beyond A.End cannot overlap, nor can any
	// later B) and test the end coordinate per record, exactly matching
	// upstream chromsweep's intersects() predicate.
	total := 0
	for _, a := range aRecs {
		bs := byChrom[a.Chrom]
		if len(bs) == 0 {
			continue
		}
		hi := sort.Search(len(bs), func(i int) bool { return bs[i].ChromStart >= a.ChromEnd })
		for j := 0; j < hi; j++ {
			b := bs[j]
			// Guaranteed b.ChromStart < a.ChromEnd by hi; reject the records
			// whose end is at or before A.Start (no positive-length overlap).
			if b.ChromEnd <= a.ChromStart {
				continue
			}
			if !strandOK(a, b, opts) {
				continue
			}
			start := a.ChromStart
			if b.ChromStart > start {
				start = b.ChromStart
			}
			end := a.ChromEnd
			if b.ChromEnd < end {
				end = b.ChromEnd
			}
			if end <= start {
				continue
			}
			ov := end - start
			if !fractionOK(a, b, ov, opts) {
				continue
			}
			total++
		}
	}
	return total
}

func strandOK(a, b *bed.Record, opts Options) bool {
	switch {
	case opts.SameStrand:
		if a.Strand == "" || a.Strand == "." || b.Strand == "" || b.Strand == "." {
			return false
		}
		return a.Strand == b.Strand
	case opts.OppositeStrand:
		if a.Strand == "" || a.Strand == "." || b.Strand == "" || b.Strand == "." {
			return false
		}
		return (a.Strand == "+" && b.Strand == "-") || (a.Strand == "-" && b.Strand == "+")
	}
	return true
}

func fractionOK(a, b *bed.Record, overlap int, opts Options) bool {
	// Reciprocal: require -f to apply to BOTH sides at the same threshold.
	fA := opts.FractionA
	fB := opts.FractionB
	if opts.Reciprocal && fA > 0 {
		fB = fA
	}
	if fA > 0 {
		lenA := a.ChromEnd - a.ChromStart
		if lenA == 0 || float64(overlap)/float64(lenA) < fA {
			return false
		}
	}
	if fB > 0 {
		lenB := b.ChromEnd - b.ChromStart
		if lenB == 0 || float64(overlap)/float64(lenB) < fB {
			return false
		}
	}
	return true
}

// mergeBED merges overlapping records on each chromosome, returning a flat
// sorted list. Strand/name are discarded — same behaviour as `bedtools merge`
// run without any keep-info options.
func mergeBED(in []*bed.Record) []*bed.Record {
	if len(in) == 0 {
		return nil
	}
	// Sort by chrom, start.
	cp := make([]*bed.Record, len(in))
	copy(cp, in)
	sort.Slice(cp, func(i, j int) bool {
		if cp[i].Chrom != cp[j].Chrom {
			return cp[i].Chrom < cp[j].Chrom
		}
		return cp[i].ChromStart < cp[j].ChromStart
	})

	var out []*bed.Record
	var cur *bed.Record
	for _, r := range cp {
		if cur == nil || r.Chrom != cur.Chrom || r.ChromStart > cur.ChromEnd {
			if cur != nil {
				out = append(out, cur)
			}
			c := *r
			cur = &c
			continue
		}
		if r.ChromEnd > cur.ChromEnd {
			cur.ChromEnd = r.ChromEnd
		}
	}
	if cur != nil {
		out = append(out, cur)
	}
	return out
}

// inputFormat tags the autodetected format of an input file.
type inputFormat int

const (
	formatUnknown inputFormat = iota
	formatBED
	formatVCF
	formatGFF
)

// isInteger reports whether s is a base-10 integer (optional leading sign),
// matching upstream's isInteger used by BedFile::parseLine.
func isInteger(s string) bool {
	if s == "" {
		return false
	}
	_, err := strconv.Atoi(s)
	return err == nil
}

// detectFormat classifies a tokenized data line as BED, VCF, or GFF using the
// same precedence as upstream's BedFile::parseLine: BED first (cols 2,3
// integer), then VCF (col 2 integer and ≥8 cols), then GFF (8 or 9 cols with
// cols 4,5 integer). Returns formatUnknown when none match.
func detectFormat(fields []string) inputFormat {
	n := len(fields)
	if n < 3 {
		return formatUnknown
	}
	if isInteger(fields[1]) && isInteger(fields[2]) {
		return formatBED
	}
	if isInteger(fields[1]) && n >= 8 {
		return formatVCF
	}
	if (n == 8 || n == 9) && isInteger(fields[3]) && isInteger(fields[4]) {
		return formatGFF
	}
	return formatUnknown
}

// parseTyped converts a tokenized data line into a bed.Record according to the
// resolved format, mirroring upstream's parseBedLine / parseVcfLine /
// parseGffLine coordinate conventions:
//
//   - BED: start = col2, end = col3 (already 0-based half-open); strand col6.
//   - VCF: start = POS-1, end = start + len(REF) (col4).
//   - GFF: start = col4-1, end = col5 (1-based inclusive → 0-based
//     half-open); strand col7.
func parseTyped(fields []string, format inputFormat) (*bed.Record, error) {
	switch format {
	case formatBED:
		start, err := strconv.Atoi(fields[1])
		if err != nil {
			return nil, fmt.Errorf("invalid BED start %q: %w", fields[1], err)
		}
		end, err := strconv.Atoi(fields[2])
		if err != nil {
			return nil, fmt.Errorf("invalid BED end %q: %w", fields[2], err)
		}
		rec := &bed.Record{Chrom: fields[0], ChromStart: start, ChromEnd: end}
		if len(fields) >= 6 {
			rec.Strand = fields[5]
		}
		return rec, nil
	case formatVCF:
		pos, err := strconv.Atoi(fields[1])
		if err != nil {
			return nil, fmt.Errorf("invalid VCF POS %q: %w", fields[1], err)
		}
		start := pos - 1
		end := start + len(fields[3])
		if start < 0 || end < start {
			return nil, fmt.Errorf("malformed VCF entry: start=%d end=%d", start, end)
		}
		return &bed.Record{Chrom: fields[0], ChromStart: start, ChromEnd: end}, nil
	case formatGFF:
		gstart, err := strconv.Atoi(fields[3])
		if err != nil {
			return nil, fmt.Errorf("invalid GFF start %q: %w", fields[3], err)
		}
		gend, err := strconv.Atoi(fields[4])
		if err != nil {
			return nil, fmt.Errorf("invalid GFF end %q: %w", fields[4], err)
		}
		start := gstart - 1
		end := gend
		if start < 0 || end < start {
			return nil, fmt.Errorf("malformed GFF entry: start=%d end=%d", start, end)
		}
		rec := &bed.Record{Chrom: fields[0], ChromStart: start, ChromEnd: end}
		if len(fields) >= 7 {
			rec.Strand = fields[6]
		}
		return rec, nil
	default:
		return nil, fmt.Errorf("unknown input format")
	}
}

// readAllBED reads every data line from r, autodetecting BED/VCF/GFF from the
// first non-header record and applying that format to the rest of the file
// (matching upstream, which locks the file type on the first data record).
// Header lines (`#`, `track`, `browser`) and blanks are skipped.
func readAllBED(r io.Reader) ([]*bed.Record, error) {
	sc := bufio.NewScanner(r)
	// GFF attribute columns can be very long; allow lines up to 16 MiB.
	sc.Buffer(make([]byte, 64*1024), 16*1024*1024)
	var out []*bed.Record
	format := formatUnknown
	for sc.Scan() {
		line := strings.TrimRight(sc.Text(), "\r")
		trimmed := strings.TrimSpace(line)
		if trimmed == "" || strings.HasPrefix(trimmed, "#") ||
			strings.HasPrefix(trimmed, "track") || strings.HasPrefix(trimmed, "browser") {
			continue
		}
		fields := strings.Split(line, "\t")
		if format == formatUnknown {
			format = detectFormat(fields)
			if format == formatUnknown {
				return nil, fmt.Errorf("could not detect input format from line: %q", line)
			}
		}
		rec, err := parseTyped(fields, format)
		if err != nil {
			return nil, err
		}
		out = append(out, rec)
	}
	return out, sc.Err()
}

// loadGenomeSize parses a 2-column chrom-sizes file (chrom\tsize) and returns
// the sum of sizes.
func loadGenomeSize(r io.Reader) (int64, error) {
	sc := bufio.NewScanner(r)
	var total int64
	for sc.Scan() {
		line := strings.TrimSpace(sc.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) < 2 {
			return 0, fmt.Errorf("genome file: expected 2+ fields, got %q", line)
		}
		n, err := strconv.ParseInt(fields[1], 10, 64)
		if err != nil {
			return 0, fmt.Errorf("genome file: bad size %q: %w", fields[1], err)
		}
		total += n
	}
	return total, sc.Err()
}

// ----- Fisher's exact test -----
//
// Ported from upstream's kfunc.cpp (htslib-style implementation).

// lgamma uses math.Lgamma; kfLgamma is a thin wrapper for parity.
func kfLgamma(z float64) float64 {
	v, _ := math.Lgamma(z)
	return v
}

// lbinom returns log(n choose k).
func lbinom(n, k int64) float64 {
	if k == 0 || n == k {
		return 0
	}
	return kfLgamma(float64(n+1)) - kfLgamma(float64(k+1)) - kfLgamma(float64(n-k+1))
}

// hypergeo returns the hypergeometric probability for table {n11; n1_; n_1; n}.
func hypergeo(n11, n1Row, n1Col, n int64) float64 {
	return math.Exp(lbinom(n1Row, n11) + lbinom(n-n1Row, n1Col-n11) - lbinom(n, n1Col))
}

// hgaccT holds the incremental state for hypergeo_acc.
type hgaccT struct {
	n11, n1Row, n1Col, n int64
	p                    float64
}

// hypergeoAcc mirrors upstream's incremental updater. When called with
// (n11, 0, 0, 0) and a non-empty aux, it only updates n11 incrementally.
func hypergeoAcc(n11, n1Row, n1Col, n int64, aux *hgaccT) float64 {
	if n1Row != 0 || n1Col != 0 || n != 0 {
		aux.n11 = n11
		aux.n1Row = n1Row
		aux.n1Col = n1Col
		aux.n = n
	} else { // only n11 changed
		// Replicate upstream's gating condition exactly.
		if n11%11 != 0 && n11+aux.n-aux.n1Row-aux.n1Col != 0 {
			if n11 == aux.n11+1 { // increment
				aux.p *= (float64(aux.n1Row-aux.n11) / float64(n11)) *
					(float64(aux.n1Col-aux.n11) / float64(n11+aux.n-aux.n1Row-aux.n1Col))
				aux.n11 = n11
				return aux.p
			}
			if n11 == aux.n11-1 { // decrement
				aux.p *= (float64(aux.n11) / float64(aux.n1Row-n11)) *
					(float64(aux.n11+aux.n-aux.n1Row-aux.n1Col) / float64(aux.n1Col-n11))
				aux.n11 = n11
				return aux.p
			}
		}
		aux.n11 = n11
	}
	aux.p = hypergeo(aux.n11, aux.n1Row, aux.n1Col, aux.n)
	return aux.p
}

// ktFisherExact returns (left, right, two-tail) p-values, mirroring htslib.
func ktFisherExact(n11, n12, n21, n22 int64) (float64, float64, float64) {
	n1Row := n11 + n12
	n1Col := n11 + n21
	n := n11 + n12 + n21 + n22

	max := n1Col
	if n1Row < max {
		max = n1Row
	}
	min := n1Row + n1Col - n
	if min < 0 {
		min = 0
	}

	if min == max {
		return 1.0, 1.0, 1.0
	}

	var aux hgaccT
	q := hypergeoAcc(n11, n1Row, n1Col, n, &aux)
	if q == 0.0 {
		// Underflow case: pick the right tail to be 0 if the mode sits to
		// the right of n11, else the left tail.
		if int64(n11)*(n+2) < (n1Col+1)*(n1Row+1) {
			return 0.0, 1.0, 0.0
		}
		return 1.0, 0.0, 0.0
	}

	// Left tail.
	p := hypergeoAcc(min, 0, 0, 0, &aux)
	var left, right float64
	i := min + 1
	for p < 0.99999999*q && i <= max {
		left += p
		p = hypergeoAcc(i, 0, 0, 0, &aux)
		i++
	}
	i--
	if p < 1.00000001*q {
		left += p
	} else {
		i--
	}

	// Right tail.
	p = hypergeoAcc(max, 0, 0, 0, &aux)
	j := max - 1
	for p < 0.99999999*q && j >= 0 {
		right += p
		p = hypergeoAcc(j, 0, 0, 0, &aux)
		j--
	}
	j++
	if p < 1.00000001*q {
		right += p
	} else {
		j++
	}

	two := left + right
	if two > 1 {
		two = 1
	}

	// Adjust left/right exactly the way upstream does so that left+right-q == 1
	// for the side that contains n11.
	if absI(i-n11) < absI(j-n11) {
		right = 1.0 - left + q
	} else {
		left = 1.0 - right + q
	}
	return left, right, two
}

func absI(x int64) int64 {
	if x < 0 {
		return -x
	}
	return x
}
