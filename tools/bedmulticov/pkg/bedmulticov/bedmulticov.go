// Package bedmulticov implements `bedtools multicov`: for each interval in
// a primary BED file (`-bed`), it reports the count of overlapping records
// from each of N input files (`-files` / `-bams`). The output is the
// original A columns followed by N integer columns — one per input file —
// holding the overlap count.
//
// Upstream supports BED *and* indexed BAM inputs. This port currently
// covers BED inputs only; the CLI surfaces a clear error if a `.bam` path
// is supplied. See tools/bedmulticov/README.md for the `t.Skip` rationale
// on the BAM cases from `reference_code/bedtools/test/multicov`.
//
// Internally each input file is loaded into a per-chromosome interval
// tree (`pkg/bioformats/bed.IntervalTree`), and the A file is streamed
// line-by-line. Optional strand filters (-s same / -S opposite),
// fraction-of-A (-f), fraction-of-B (-F), and reciprocal (-r) thresholds
// mirror upstream's semantics.
package bedmulticov

import (
	"bufio"
	"fmt"
	"io"
	"sort"
	"strconv"
	"strings"

	bedpkg "github.com/yassineS/bio_ai_experiment/pkg/bioformats/bed"
)

// Options configures Run.
type Options struct {
	// FractionA is the minimum fraction of A that must be covered by a
	// single B record for it to count as overlapping (mirrors `-f`).
	// 0 = any positive overlap counts. Range: (0, 1].
	FractionA float64
	// FractionB is the analogous minimum fraction of B covered by A
	// (mirrors `-F`). 0 = unconstrained.
	FractionB float64
	// Reciprocal mirrors `-r`: when set together with FractionA, the
	// same threshold is also applied to B (equivalent to FractionB=FractionA).
	Reciprocal bool
	// SameStrand mirrors `-s`: only count B records on the same strand as A.
	SameStrand bool
	// OppositeStrand mirrors `-S`: only count B records on the opposite
	// strand from A.
	OppositeStrand bool
}

// Run reads A from aR and the N B files from bRs in order, indexes each B
// into per-chromosome interval trees, then streams A and emits one row per
// A record with one count column appended per B file. Returns the number
// of A records processed.
func Run(aR io.Reader, bRs []io.Reader, out io.Writer, opts Options) (int, error) {
	if opts.SameStrand && opts.OppositeStrand {
		return 0, fmt.Errorf("cannot combine -s and -S")
	}
	if opts.FractionA < 0 || opts.FractionA > 1 {
		return 0, fmt.Errorf("-f must be in [0,1], got %g", opts.FractionA)
	}
	if opts.FractionB < 0 || opts.FractionB > 1 {
		return 0, fmt.Errorf("-F must be in [0,1], got %g", opts.FractionB)
	}
	if opts.Reciprocal && opts.FractionA <= 0 {
		return 0, fmt.Errorf("-r requires -f to be specified")
	}
	// Reciprocal: apply FractionA threshold to B as well.
	effFracB := opts.FractionB
	if opts.Reciprocal && effFracB < opts.FractionA {
		effFracB = opts.FractionA
	}

	trees := make([]map[string]*bedpkg.IntervalTree, len(bRs))
	for i, br := range bRs {
		t, err := indexB(br)
		if err != nil {
			return 0, fmt.Errorf("file %d: %w", i+1, err)
		}
		trees[i] = t
	}

	bw := bufio.NewWriter(out)
	defer bw.Flush()

	sc := bufio.NewScanner(aR)
	sc.Buffer(make([]byte, 64*1024), 16*1024*1024)
	count := 0
	lineNo := 0
	for sc.Scan() {
		lineNo++
		raw := sc.Text()
		trimmed := strings.TrimSpace(raw)
		if trimmed == "" || strings.HasPrefix(trimmed, "#") ||
			strings.HasPrefix(trimmed, "track") || strings.HasPrefix(trimmed, "browser") {
			continue
		}
		fields := strings.Split(raw, "\t")
		if len(fields) < 3 {
			return count, fmt.Errorf("line %d: BED record needs >=3 columns: %q", lineNo, raw)
		}
		rec, err := parseRecord(fields)
		if err != nil {
			return count, fmt.Errorf("line %d: %w", lineNo, err)
		}
		// Emit A's original columns verbatim, then one count per B file.
		if _, err := bw.WriteString(strings.Join(fields, "\t")); err != nil {
			return count, err
		}
		for _, t := range trees {
			n := countOverlaps(rec, t[rec.Chrom], opts, effFracB)
			if _, err := fmt.Fprintf(bw, "\t%d", n); err != nil {
				return count, err
			}
		}
		if err := bw.WriteByte('\n'); err != nil {
			return count, err
		}
		count++
	}
	if err := sc.Err(); err != nil {
		return count, err
	}
	return count, nil
}

// countOverlaps returns the number of B records that overlap a after
// applying strand and fraction filters.
func countOverlaps(a *bedpkg.Record, t *bedpkg.IntervalTree, opts Options, effFracB float64) int {
	if t == nil {
		return 0
	}
	cand := t.Query(a)
	if len(cand) == 0 {
		return 0
	}
	lenA := a.ChromEnd - a.ChromStart
	n := 0
	for _, b := range cand {
		if !strandOK(a, b, opts) {
			continue
		}
		overlap := overlapLen(a, b)
		if overlap <= 0 {
			continue
		}
		if opts.FractionA > 0 && lenA > 0 {
			if float64(overlap)/float64(lenA) < opts.FractionA {
				continue
			}
		}
		if effFracB > 0 {
			lenB := b.ChromEnd - b.ChromStart
			if lenB <= 0 {
				continue
			}
			if float64(overlap)/float64(lenB) < effFracB {
				continue
			}
		}
		n++
	}
	return n
}

// overlapLen returns the length of the intersection of a and b's spans.
// 0 if disjoint.
func overlapLen(a, b *bedpkg.Record) int {
	start := a.ChromStart
	if b.ChromStart > start {
		start = b.ChromStart
	}
	end := a.ChromEnd
	if b.ChromEnd < end {
		end = b.ChromEnd
	}
	if end <= start {
		return 0
	}
	return end - start
}

// strandOK applies the -s / -S filters. Missing strand on either side is
// treated as "no match" under a strand filter (matches upstream).
func strandOK(a, b *bedpkg.Record, opts Options) bool {
	if opts.SameStrand {
		if a.Strand == "" || b.Strand == "" {
			return false
		}
		return a.Strand == b.Strand
	}
	if opts.OppositeStrand {
		if a.Strand == "" || b.Strand == "" {
			return false
		}
		return a.Strand != b.Strand
	}
	return true
}

// parseRecord parses the minimum subset of a BED line we need for overlap
// + strand filtering. Extra columns are preserved by the caller as raw
// fields.
func parseRecord(fields []string) (*bedpkg.Record, error) {
	start, err := strconv.Atoi(strings.TrimSpace(fields[1]))
	if err != nil {
		return nil, fmt.Errorf("invalid chromStart %q: %v", fields[1], err)
	}
	end, err := strconv.Atoi(strings.TrimSpace(fields[2]))
	if err != nil {
		return nil, fmt.Errorf("invalid chromEnd %q: %v", fields[2], err)
	}
	r := &bedpkg.Record{
		Chrom:      fields[0],
		ChromStart: start,
		ChromEnd:   end,
	}
	if len(fields) >= 6 {
		r.Strand = fields[5]
	}
	return r, nil
}

// indexB reads a B file fully into memory and returns a per-chrom tree.
func indexB(r io.Reader) (map[string]*bedpkg.IntervalTree, error) {
	rd := bedpkg.NewReader(r)
	all, err := rd.ReadAll()
	if err != nil {
		return nil, err
	}
	byChrom := map[string][]*bedpkg.Record{}
	for _, x := range all {
		byChrom[x.Chrom] = append(byChrom[x.Chrom], x)
	}
	out := make(map[string]*bedpkg.IntervalTree, len(byChrom))
	for chrom, recs := range byChrom {
		sort.SliceStable(recs, func(i, j int) bool {
			if recs[i].ChromStart != recs[j].ChromStart {
				return recs[i].ChromStart < recs[j].ChromStart
			}
			return recs[i].ChromEnd < recs[j].ChromEnd
		})
		out[chrom] = bedpkg.NewIntervalTree(recs)
	}
	return out, nil
}
