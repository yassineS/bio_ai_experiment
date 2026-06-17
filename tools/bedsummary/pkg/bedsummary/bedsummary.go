// Package bedsummary reports per-chromosome summary statistics for the
// intervals in a BED/GFF/VCF file, mirroring the bedtools `summary` subcommand.
//
// For every chromosome listed in the genome (-g) file — in genome-file order —
// it reports the chromosome length, number of intervals, total interval bp, the
// chromosome's fraction of the genome, its fraction of all intervals and of all
// interval bp, and the min/max/mean interval length. A final "all" row
// aggregates over the whole input. The output is a 10-column TSV whose header,
// column set, ordering, fixed-9 decimal precision, and per-row trailing tab
// match upstream byte-for-byte.
package bedsummary

import (
	"bufio"
	"fmt"
	"io"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bed"
)

// Options controls Run behaviour.
type Options struct {
	// NoHeader, when true, suppresses the column-header line. (Not an upstream
	// flag; provided for callers that want the raw rows.)
	NoHeader bool
}

// genome holds the parsed -g chrom-sizes file: per-chrom sizes plus the order
// in which chromosomes appear in the file (which determines output order).
type genome struct {
	sizes map[string]int64
	order []string
	total int64
}

// chromStats holds the accumulated interval data for one chromosome.
type chromStats struct {
	count   int64
	totalBP int64
	minLen  int64
	maxLen  int64
}

// ParseGenome reads a 2-column (chrom\tsize) genome file, preserving the
// chromosome order and summing the total genome size. Blank lines and lines
// beginning with '#' are skipped.
func ParseGenome(r io.Reader) (*genome, error) {
	g := &genome{sizes: map[string]int64{}}
	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 64*1024), 16*1024*1024)
	for sc.Scan() {
		line := strings.TrimRight(sc.Text(), "\r")
		if strings.TrimSpace(line) == "" || strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) < 2 {
			return nil, fmt.Errorf("genome file: expected 2 columns, got %q", line)
		}
		size, err := strconv.ParseInt(fields[1], 10, 64)
		if err != nil {
			return nil, fmt.Errorf("genome file: bad size %q: %w", fields[1], err)
		}
		if _, seen := g.sizes[fields[0]]; !seen {
			g.order = append(g.order, fields[0])
		}
		g.sizes[fields[0]] = size
		g.total += size
	}
	if err := sc.Err(); err != nil {
		return nil, err
	}
	return g, nil
}

// Run reads BED records from r, computes the summary against the genome g, and
// writes the upstream-format report to w.
func Run(r io.Reader, g *genome, w io.Writer, opts Options) error {
	if g == nil {
		return fmt.Errorf("a genome (-g) file is required")
	}

	stats := map[string]*chromStats{}
	var totalIntervals, totalLength int64

	br := bed.NewReader(r)
	for {
		rec, err := br.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return fmt.Errorf("reading BED: %w", err)
		}
		if _, ok := g.sizes[rec.Chrom]; !ok {
			// Upstream aborts when an interval's chromosome is absent from the
			// genome file (it cannot compute that chromosome's genome fraction).
			return fmt.Errorf("requested chromosome %s does not exist in the genome file. Exiting.", rec.Chrom)
		}
		length := int64(rec.ChromEnd - rec.ChromStart)
		s := stats[rec.Chrom]
		if s == nil {
			s = &chromStats{minLen: length, maxLen: length}
			stats[rec.Chrom] = s
		}
		s.count++
		s.totalBP += length
		if length < s.minLen {
			s.minLen = length
		}
		if length > s.maxLen {
			s.maxLen = length
		}
		totalIntervals++
		totalLength += length
	}

	bw := bufio.NewWriter(w)
	defer bw.Flush()

	if !opts.NoHeader {
		if _, err := fmt.Fprintln(bw, "chrom\tchrom_length\tnum_ivls\ttotal_ivl_bp\t"+
			"chrom_frac_genome\tfrac_all_ivls\tfrac_all_bp\tmin\tmax\tmean"); err != nil {
			return err
		}
	}

	var overallMin int64 = -1 // sentinel; set on first data chromosome
	var overallMax int64

	for _, chrom := range g.order {
		if chrom == "" {
			continue
		}
		chromSize := g.sizes[chrom]
		pctGenome := float64(chromSize) / float64(g.total)

		s := stats[chrom]
		if s == nil || s.count == 0 {
			// No intervals on this chromosome: upstream emits a default row
			// with -1 for min/max/mean and 0 for the count/bp/fractions, with
			// NO trailing tab.
			if _, err := fmt.Fprintf(bw, "%s\t%d\t0\t0\t%s\t0.000000000\t0.000000000\t-1\t-1\t-1\n",
				chrom, chromSize, frac9(pctGenome)); err != nil {
				return err
			}
			continue
		}

		fracAllIvls := float64(s.count) / float64(totalIntervals)
		fracAllBP := float64(s.totalBP) / float64(totalLength)
		mean := float64(s.totalBP) / float64(s.count)

		if overallMin == -1 || s.minLen < overallMin {
			overallMin = s.minLen
		}
		if s.maxLen > overallMax {
			overallMax = s.maxLen
		}

		// Data rows END WITH a trailing tab, matching upstream.
		if _, err := fmt.Fprintf(bw, "%s\t%d\t%d\t%d\t%s\t%s\t%s\t%d\t%d\t%s\t\n",
			chrom, chromSize, s.count, s.totalBP,
			frac9(pctGenome), frac9(fracAllIvls), frac9(fracAllBP),
			s.minLen, s.maxLen, frac9(mean)); err != nil {
			return err
		}
	}

	// The "all" aggregate row. Upstream prints the literal "1.0" for the three
	// fraction columns and has no trailing tab.
	var allMean string
	if totalIntervals > 0 {
		allMean = frac9(float64(totalLength) / float64(totalIntervals))
	} else {
		allMean = frac9(0)
	}
	if overallMin == -1 {
		overallMin = 0
	}
	if _, err := fmt.Fprintf(bw, "all\t%d\t%d\t%d\t1.0\t1.0\t1.0\t%d\t%d\t%s\n",
		g.total, totalIntervals, totalLength, overallMin, overallMax, allMean); err != nil {
		return err
	}
	return nil
}

// frac9 formats a floating-point value with fixed 9-decimal precision, matching
// upstream's `std::fixed << std::setprecision(9)`.
func frac9(v float64) string {
	return strconv.FormatFloat(v, 'f', 9, 64)
}
