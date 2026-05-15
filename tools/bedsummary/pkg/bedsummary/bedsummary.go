// Package bedsummary computes per-chromosome summary statistics for a BED
// file: number of intervals, total length covered, and per-interval length
// min / max / mean / median.
//
// This is the Go port of the bedtools `summary` subcommand. Output is a
// 7-column TSV: chrom, num_intervals, total_length, min_length, max_length,
// mean_length, median_length. A trailing "all" row aggregates over all
// chromosomes.
package bedsummary

import (
	"bufio"
	"fmt"
	"io"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/bed"
)

// Options controls Run behaviour. Currently only output flags exist; future
// upstream additions (e.g. genome-fraction) can be added here.
type Options struct {
	// NoHeader, when true, suppresses the column-header line.
	NoHeader bool
	// SkipAll, when true, suppresses the trailing "all" aggregate row.
	SkipAll bool
}

// ChromSummary holds the stats for one chromosome (or the "all" aggregate).
type ChromSummary struct {
	Chrom       string
	Count       int
	TotalLength int
	MinLength   int
	MaxLength   int
	MeanLength  float64
	MedianLen   float64
}

// Compute reads BED records from r and returns per-chromosome summaries
// in the original-input chromosome order, plus an "all" aggregate as the
// final element (unless opts.SkipAll is set).
func Compute(r io.Reader, opts Options) ([]ChromSummary, error) {
	br := bed.NewReader(r)
	lengthsByChrom := map[string][]int{}
	var chromOrder []string

	for {
		rec, err := br.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, fmt.Errorf("reading BED: %w", err)
		}
		if _, seen := lengthsByChrom[rec.Chrom]; !seen {
			chromOrder = append(chromOrder, rec.Chrom)
		}
		lengthsByChrom[rec.Chrom] = append(lengthsByChrom[rec.Chrom], rec.ChromEnd-rec.ChromStart)
	}

	out := make([]ChromSummary, 0, len(chromOrder)+1)
	allLengths := make([]int, 0)
	for _, c := range chromOrder {
		lens := lengthsByChrom[c]
		out = append(out, summarise(c, lens))
		allLengths = append(allLengths, lens...)
	}
	if !opts.SkipAll && len(allLengths) > 0 {
		out = append(out, summarise("all", allLengths))
	}
	return out, nil
}

// summarise builds a ChromSummary from a slice of interval lengths.
func summarise(chrom string, lens []int) ChromSummary {
	sorted := append([]int(nil), lens...)
	sort.Ints(sorted)

	total := 0
	for _, l := range sorted {
		total += l
	}
	n := len(sorted)
	var med float64
	if n%2 == 1 {
		med = float64(sorted[n/2])
	} else {
		med = float64(sorted[n/2-1]+sorted[n/2]) / 2.0
	}
	mean := float64(total) / float64(n)
	return ChromSummary{
		Chrom:       chrom,
		Count:       n,
		TotalLength: total,
		MinLength:   sorted[0],
		MaxLength:   sorted[n-1],
		MeanLength:  mean,
		MedianLen:   med,
	}
}

// Run is the streaming entry point used by the CLI: read BED from r, write
// the formatted summary table to w.
func Run(r io.Reader, w io.Writer, opts Options) error {
	rows, err := Compute(r, opts)
	if err != nil {
		return err
	}
	bw := bufio.NewWriter(w)
	defer bw.Flush()

	if !opts.NoHeader {
		if _, err := fmt.Fprintln(bw, strings.Join([]string{
			"chrom", "num_ivls", "total_ivl_bp",
			"min_ivl_bp", "max_ivl_bp", "mean_ivl_bp", "median_ivl_bp",
		}, "\t")); err != nil {
			return err
		}
	}
	for _, row := range rows {
		if _, err := fmt.Fprintf(bw, "%s\t%d\t%d\t%d\t%d\t%s\t%s\n",
			row.Chrom, row.Count, row.TotalLength,
			row.MinLength, row.MaxLength,
			formatNum(row.MeanLength), formatNum(row.MedianLen),
		); err != nil {
			return err
		}
	}
	return nil
}

// formatNum prints integer-valued floats without a fractional part, otherwise
// uses three-digit-decimal precision (matches what we emit for fractions
// elsewhere in the project).
func formatNum(v float64) string {
	if v == float64(int64(v)) {
		return strconv.FormatInt(int64(v), 10)
	}
	return strconv.FormatFloat(v, 'f', 3, 64)
}
