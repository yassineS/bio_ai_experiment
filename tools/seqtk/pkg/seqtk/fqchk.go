// fqchk.go: implementation of `seqtk fqchk`.
//
// Upstream source: reference_code/seqtk/seqtk.c::stk_fqchk and the helper
// stk_fqchk's static fqc_aux (v1.5-r133, lines 1841-1942).
//
// Behaviour: walks one FASTQ stream, accumulating per-position base
// composition (A / C / G / T / N) and per-position quality histograms,
// then emits a TSV report on stdout.
//
// Output format (matches upstream byte-for-byte):
//
//	min_len: <min>; max_len: <max>; avg_len: <%.2f>; <K> distinct quality values
//	POS\t#bases\t%A\t%C\t%G\t%T\t%N\tavgQ\terrQ\t...     <- header
//	ALL\t...                                              <- aggregate row
//	1\t...                                                <- per-position rows
//	2\t...
//	...
//
// Per-row columns after "errQ" depend on the -q threshold:
//
//   - qthres <= 0  -> one "%Qk" column per distinct observed quality k
//   - qthres >  0  -> exactly two columns: "%low" (q < qthres) and
//     "%high" (q >= qthres)
//
// Notes:
//   - "errQ" is upstream's "-4.343 * log((psum + 1e-6) / (sum + 1e-6))",
//     where psum is sum of perr[q] over all bases at this position, with
//     perr[0..3] all forced to 0.5 (matching upstream's seqtk.c:1894-95).
//   - Quality is decoded as PHRED+33 by upstream (offset 33; not
//     configurable). Quality values are clamped to [0, 93].
//   - Records with zero quality length are skipped (matching upstream's
//     `if (seq->qual.l == 0) continue`).
//
// Upstream getopt surface (seqtk.c:1879): "q:"
//
//   -q INT   quality threshold for the %low / %high split [20].
//            "-q0" prints the full per-quality distribution instead.
//
// The "-o/--output FILE" flag wired in cmd/ is the project-wide Go-port
// convenience and does not affect parity.

package seqtk

import (
	"bufio"
	"fmt"
	"io"
	"math"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fastq"
)

// FqchkOptions configures Fqchk. Only -q is exposed (matches upstream).
type FqchkOptions struct {
	// QThres is upstream's -q: bases with quality < QThres count toward
	// "%low". Set to 0 to print the full per-quality distribution.
	QThres int
}

// DefaultFqchkQThres mirrors upstream's default at seqtk.c:1874.
const DefaultFqchkQThres = 20

// posStat is the per-position accumulator. Quality histogram is sized
// to 94 entries (q ∈ [0, 93]) and base counts to 5 (A / C / G / T / N).
type posStat struct {
	q [94]int64
	b [5]int64
}

// Fqchk reads FASTQ records from r and writes the per-position report
// to w. The report follows upstream "seqtk fqchk" byte-for-byte.
//
// The caller is responsible for opening / closing r (use the same
// OpenInput helper the other subcommands use; it transparently handles
// gzip and "-" for stdin).
func Fqchk(r io.Reader, w io.Writer, opts FqchkOptions) error {
	const offset = 33 // PHRED+33; upstream is not configurable
	reader := fastq.NewReader(r, fastq.Phred33)

	// perr[k] = 10 ** (-k/10), with perr[0..3] = 0.5 to match upstream
	// seqtk.c:1895. The 0.5 clamp keeps the "errQ" formula well-defined
	// for very low quality scores that would otherwise blow up.
	var perr [94]float64
	for k := 0; k <= 93; k++ {
		perr[k] = math.Pow(10, -0.1*float64(k))
	}
	perr[0], perr[1], perr[2], perr[3] = 0.5, 0.5, 0.5, 0.5

	var (
		pos    []posStat
		all    posStat
		nReads int64
		totLen int64
		minLen int64 = math.MaxInt64
		maxLen int64
	)

	for {
		rec, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return fmt.Errorf("fqchk: read: %w", err)
		}
		if len(rec.Quality) == 0 {
			continue
		}
		nReads++
		l := int64(len(rec.Sequence))
		totLen += l
		if l < minLen {
			minLen = l
		}
		if l > maxLen {
			maxLen = l
		}
		if int(maxLen) > len(pos) {
			// Grow to at least maxLen entries. We don't bother
			// imitating upstream's kroundup64 growth — Go's
			// append handles amortised growth.
			grown := make([]posStat, maxLen)
			copy(grown, pos)
			pos = grown
		}
		for i := 0; i < len(rec.Quality); i++ {
			q := int(rec.Quality[i]) - offset
			if q < 0 {
				q = 0
			}
			if q > 93 {
				q = 93
			}
			b := int(seqNT6Table[rec.Sequence[i]])
			// upstream: b = b? b - 1 : 4 (1..4 -> 0..3, else -> 4)
			if b == 0 || b == 5 {
				b = 4
			} else {
				b--
			}
			pos[i].q[q]++
			pos[i].b[b]++
		}
	}

	// If there were no qualified records at all, fall back to zero
	// values; upstream would divide by n==0 in the avg_len line and
	// produce "nan" — match that for full parity.
	bw := bufio.NewWriter(w)

	// Aggregate the per-position arrays into `all`.
	for i := 0; i < int(maxLen); i++ {
		for k := 0; k <= 93; k++ {
			all.q[k] += pos[i].q[k]
		}
		for k := 0; k <= 4; k++ {
			all.b[k] += pos[i].b[k]
		}
	}
	nDiffQ := 0
	for k := 0; k <= 93; k++ {
		if all.q[k] > 0 {
			nDiffQ++
		}
	}
	if minLen == math.MaxInt64 {
		minLen = 0
	}
	avg := float64(0)
	if nReads > 0 {
		avg = float64(totLen) / float64(nReads)
	}
	if _, err := fmt.Fprintf(bw, "min_len: %d; max_len: %d; avg_len: %.2f; %d distinct quality values\n",
		minLen, maxLen, avg, nDiffQ); err != nil {
		return err
	}

	// Header row.
	if _, err := bw.WriteString("POS\t#bases\t%A\t%C\t%G\t%T\t%N\tavgQ\terrQ"); err != nil {
		return err
	}
	if opts.QThres <= 0 {
		for k := 0; k <= 93; k++ {
			if all.q[k] > 0 {
				if _, err := fmt.Fprintf(bw, "\t%%Q%d", k); err != nil {
					return err
				}
			}
		}
	} else {
		if _, err := bw.WriteString("\t%low\t%high"); err != nil {
			return err
		}
	}
	if err := bw.WriteByte('\n'); err != nil {
		return err
	}

	// ALL row, then per-position rows.
	if err := fqchkAux(bw, &all, 0, &all, &perr, opts.QThres); err != nil {
		return err
	}
	for i := 0; i < int(maxLen); i++ {
		if err := fqchkAux(bw, &pos[i], int64(i+1), &all, &perr, opts.QThres); err != nil {
			return err
		}
	}
	return bw.Flush()
}

// fqchkAux writes one row to bw — either the "ALL" aggregate row (when
// pos <= 0) or a per-position row. Mirrors upstream's static fqc_aux.
// allRow.q is used to decide which %Qk columns to print when qthres <=
// 0; perr supplies the per-quality error probabilities.
func fqchkAux(bw *bufio.Writer, p *posStat, pos int64, allRow *posStat, perr *[94]float64, qthres int) error {
	var (
		sum    int64
		qsum   int64
		sumLow int64
		psum   float64
	)
	if pos <= 0 {
		if _, err := bw.WriteString("ALL"); err != nil {
			return err
		}
	} else {
		if _, err := fmt.Fprintf(bw, "%d", pos); err != nil {
			return err
		}
	}
	for k := 0; k <= 4; k++ {
		sum += p.b[k]
	}
	if _, err := fmt.Fprintf(bw, "\t%d", sum); err != nil {
		return err
	}
	for k := 0; k <= 4; k++ {
		if _, err := fmt.Fprintf(bw, "\t%.1f", 100.0*float64(p.b[k])/float64(sum)); err != nil {
			return err
		}
	}
	for k := 0; k <= 93; k++ {
		qsum += p.q[k] * int64(k)
		psum += float64(p.q[k]) * perr[k]
		if k < qthres {
			sumLow += p.q[k]
		}
	}
	avgQ := float64(qsum) / float64(sum)
	errQ := -4.343 * math.Log((psum+1e-6)/(float64(sum)+1e-6))
	if _, err := fmt.Fprintf(bw, "\t%.1f\t%.1f", avgQ, errQ); err != nil {
		return err
	}
	if qthres <= 0 {
		for k := 0; k <= 93; k++ {
			if allRow.q[k] > 0 {
				if _, err := fmt.Fprintf(bw, "\t%.2f", 100.0*float64(p.q[k])/float64(sum)); err != nil {
					return err
				}
			}
		}
	} else {
		if _, err := fmt.Fprintf(bw, "\t%.1f\t%.1f",
			100.0*float64(sumLow)/float64(sum),
			100.0*float64(sum-sumLow)/float64(sum)); err != nil {
			return err
		}
	}
	return bw.WriteByte('\n')
}
