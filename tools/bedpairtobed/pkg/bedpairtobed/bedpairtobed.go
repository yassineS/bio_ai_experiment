// Package bedpairtobed implements `bedtools pairtobed`: report overlaps
// between a BEDPE A file (paired-end intervals) and a regular BED B file.
//
// Upstream reference:
//
//	reference_code/bedtools/src/pairToBed/pairToBed.cpp
//	reference_code/bedtools/src/pairToBed/pairToBedMain.cpp
//
// Algorithm (mirrors upstream FindOverlaps / FindOneOrMoreOverlaps):
//
//   - For each A pair, find B records whose coordinates overlap end1 (chrom1,
//     start1..end1) AND records that overlap end2.
//   - An overlap is counted only if the fraction of end-i covered by B is
//     >= MinFractionA (upstream's -f).
//   - The Type field decides which A records are emitted and in what shape:
//     either / both / notboth / neither / xor / notxor (notxor is our extra
//     completeness option; upstream supports the first 5).
//
// We build one per-chromosome interval tree over B once for an O(log n + k)
// query path, as in bedintersect.
//
// Output: each emitted A pair is followed by a tab + the full BED line of the
// hitting B record (for "either" and "both" and "notboth"-with-hits and "xor"
// modes). "neither" and the no-hit branch of "notboth"/"notxor" emit just the
// BEDPE line.

package bedpairtobed

import (
	"fmt"
	"io"
	"sort"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bed"
)

// Type enumerates the -type modes supported by `bedtools pairtobed`.
type Type string

// The five upstream modes plus our completeness option notxor.
const (
	TypeEither  Type = "either"
	TypeBoth    Type = "both"
	TypeNotboth Type = "notboth"
	TypeNeither Type = "neither"
	TypeXor     Type = "xor"
	TypeNotxor  Type = "notxor"
)

// ValidTypes lists the accepted -type values for CLI validation.
func ValidTypes() []string {
	return []string{
		string(TypeEither), string(TypeBoth), string(TypeNotboth),
		string(TypeNeither), string(TypeXor), string(TypeNotxor),
	}
}

// Options controls a Run invocation.
type Options struct {
	// Type selects which A records (and how) are emitted. Default: TypeEither.
	Type Type
	// MinFractionA requires this fraction of an A end's length to be covered
	// by the B record. Default upstream is 1e-9 (effectively 1bp).
	MinFractionA float64
	// SameStrand requires the BEDPE end and the BED hit to share strand.
	SameStrand bool
	// OppositeStrand requires the BEDPE end and the BED hit to have
	// different strands. Mutually exclusive with SameStrand.
	OppositeStrand bool
	// IgnoreStrand: if set, all strand checks are skipped (mirrors -is).
	IgnoreStrand bool
}

// Stats records the number of input pairs and emitted output lines.
type Stats struct {
	APairs       int
	BRecords     int
	OutputLines  int
	EmittedPairs int // number of unique A pairs that produced at least one output line
}

// Run reads BEDPE pairs from a, BED records from b and writes results to w.
func Run(a, b io.Reader, w io.Writer, opts Options) (*Stats, error) {
	if opts.Type == "" {
		opts.Type = TypeEither
	}
	if !isValidType(opts.Type) {
		return nil, fmt.Errorf("invalid -type %q (want one of %s)", opts.Type, strings.Join(ValidTypes(), ", "))
	}
	if opts.SameStrand && opts.OppositeStrand {
		return nil, fmt.Errorf("-s and -S are mutually exclusive")
	}
	if opts.MinFractionA <= 0 {
		opts.MinFractionA = 1e-9
	}

	// Index B by chrom and build per-chrom interval trees.
	bReader := bed.NewReader(b)
	bAll, err := bReader.ReadAll()
	if err != nil {
		return nil, fmt.Errorf("read B: %w", err)
	}
	stats := &Stats{BRecords: len(bAll)}
	byChrom := make(map[string][]*bed.Record)
	for _, r := range bAll {
		byChrom[r.Chrom] = append(byChrom[r.Chrom], r)
	}
	trees := make(map[string]*bed.IntervalTree, len(byChrom))
	for chrom, recs := range byChrom {
		sort.SliceStable(recs, func(i, j int) bool { return recs[i].ChromStart < recs[j].ChromStart })
		trees[chrom] = bed.NewIntervalTree(recs)
	}

	// Buffer output to a string builder per record so we can attribute writes
	// to Stats.EmittedPairs accurately.
	out := &writer{w: w}
	aReader := bed.NewBEDPEReader(a)
	for {
		pair, err := aReader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return stats, err
		}
		stats.APairs++

		hits1 := findEndHits(pair.Chrom1, pair.Start1, pair.End1, pair.Strand1, trees, &opts)
		hits2 := findEndHits(pair.Chrom2, pair.Start2, pair.End2, pair.Strand2, trees, &opts)

		emitted, lines, err := emit(out, pair, hits1, hits2, opts.Type)
		if err != nil {
			return stats, err
		}
		if emitted {
			stats.EmittedPairs++
		}
		stats.OutputLines += lines
	}
	return stats, out.flush()
}

func isValidType(t Type) bool {
	for _, v := range ValidTypes() {
		if string(t) == v {
			return true
		}
	}
	return false
}

// findEndHits returns the B records overlapping (chrom, start, end) by at
// least opts.MinFractionA of the end length. An unaligned end (chrom == "."
// or any of start/end < 0) returns no hits.
func findEndHits(chrom string, start, end int, strand string, trees map[string]*bed.IntervalTree, opts *Options) []*bed.Record {
	if chrom == "" || chrom == "." || start < 0 || end < 0 || end <= start {
		return nil
	}
	tree, ok := trees[chrom]
	if !ok {
		return nil
	}
	q := &bed.Record{Chrom: chrom, ChromStart: start, ChromEnd: end}
	cands := tree.Query(q)
	if len(cands) == 0 {
		return nil
	}
	endLen := end - start
	out := make([]*bed.Record, 0, len(cands))
	for _, c := range cands {
		if !strandOK(strand, c.Strand, opts) {
			continue
		}
		ov := overlap(start, end, c.ChromStart, c.ChromEnd)
		if ov <= 0 {
			continue
		}
		if float64(ov)/float64(endLen) < opts.MinFractionA {
			continue
		}
		out = append(out, c)
	}
	return out
}

func strandOK(aStrand, bStrand string, opts *Options) bool {
	if opts.IgnoreStrand {
		return true
	}
	if opts.SameStrand {
		// Treat empty / "." as wildcards: upstream only enforces when both
		// sides have a strand annotation.
		if aStrand == "" || bStrand == "" || aStrand == "." || bStrand == "." {
			return true
		}
		return aStrand == bStrand
	}
	if opts.OppositeStrand {
		if aStrand == "" || bStrand == "" || aStrand == "." || bStrand == "." {
			return true
		}
		return aStrand != bStrand
	}
	return true
}

func overlap(aS, aE, bS, bE int) int {
	if aS < bS {
		aS = bS
	}
	if aE > bE {
		aE = bE
	}
	if aE <= aS {
		return 0
	}
	return aE - aS
}

// emit dispatches by Type and writes output. Returns (anyEmittedForThisPair,
// numLines, err).
func emit(w *writer, pair *bed.BEDPE, hits1, hits2 []*bed.Record, t Type) (bool, int, error) {
	switch t {
	case TypeEither:
		// One line per hit, on either end.
		lines := 0
		for _, h := range hits1 {
			if err := w.writePairWithBed(pair, h); err != nil {
				return false, lines, err
			}
			lines++
		}
		for _, h := range hits2 {
			if err := w.writePairWithBed(pair, h); err != nil {
				return false, lines, err
			}
			lines++
		}
		return lines > 0, lines, nil

	case TypeBoth:
		// Both ends must hit at least one B record. Upstream emits one line
		// per hit on each end (not the cross product).
		if len(hits1) == 0 || len(hits2) == 0 {
			return false, 0, nil
		}
		lines := 0
		for _, h := range hits1 {
			if err := w.writePairWithBed(pair, h); err != nil {
				return false, lines, err
			}
			lines++
		}
		for _, h := range hits2 {
			if err := w.writePairWithBed(pair, h); err != nil {
				return false, lines, err
			}
			lines++
		}
		return true, lines, nil

	case TypeNotboth:
		// Upstream's notboth: emit the pair with no B when both ends miss;
		// emit pair+hit lines for the side that has hits when only one side
		// hits; emit nothing when both ends hit.
		if len(hits1) == 0 && len(hits2) == 0 {
			if err := w.writePair(pair); err != nil {
				return false, 0, err
			}
			return true, 1, nil
		}
		if len(hits1) > 0 && len(hits2) > 0 {
			return false, 0, nil
		}
		lines := 0
		for _, h := range hits1 {
			if err := w.writePairWithBed(pair, h); err != nil {
				return false, lines, err
			}
			lines++
		}
		for _, h := range hits2 {
			if err := w.writePairWithBed(pair, h); err != nil {
				return false, lines, err
			}
			lines++
		}
		return true, lines, nil

	case TypeNeither:
		if len(hits1) == 0 && len(hits2) == 0 {
			if err := w.writePair(pair); err != nil {
				return false, 0, err
			}
			return true, 1, nil
		}
		return false, 0, nil

	case TypeXor:
		// Exactly one end hits.
		if (len(hits1) > 0) == (len(hits2) > 0) {
			return false, 0, nil
		}
		lines := 0
		for _, h := range hits1 {
			if err := w.writePairWithBed(pair, h); err != nil {
				return false, lines, err
			}
			lines++
		}
		for _, h := range hits2 {
			if err := w.writePairWithBed(pair, h); err != nil {
				return false, lines, err
			}
			lines++
		}
		return lines > 0, lines, nil

	case TypeNotxor:
		// Either both hit, or neither hits. Mirrors upstream behaviour when
		// users want the complement of xor (deviation: upstream doesn't ship
		// this; documented in README).
		switch {
		case len(hits1) == 0 && len(hits2) == 0:
			if err := w.writePair(pair); err != nil {
				return false, 0, err
			}
			return true, 1, nil
		case len(hits1) > 0 && len(hits2) > 0:
			lines := 0
			for _, h := range hits1 {
				if err := w.writePairWithBed(pair, h); err != nil {
					return false, lines, err
				}
				lines++
			}
			for _, h := range hits2 {
				if err := w.writePairWithBed(pair, h); err != nil {
					return false, lines, err
				}
				lines++
			}
			return true, lines, nil
		}
		return false, 0, nil
	}
	return false, 0, fmt.Errorf("unhandled type %q", t)
}

// writer wraps an io.Writer and tracks errors. We intentionally avoid a
// dedicated bed/bedpe formatter because the output mixes BEDPE-shape and
// BED-shape lines.
type writer struct {
	w   io.Writer
	err error
}

func (w *writer) flush() error { return w.err }

func (w *writer) writePair(p *bed.BEDPE) error {
	if w.err != nil {
		return w.err
	}
	_, err := io.WriteString(w.w, p.String()+"\n")
	w.err = err
	return err
}

func (w *writer) writePairWithBed(p *bed.BEDPE, b *bed.Record) error {
	if w.err != nil {
		return w.err
	}
	_, err := io.WriteString(w.w, p.String()+"\t"+formatBED(b)+"\n")
	w.err = err
	return err
}

// formatBED renders a *bed.Record into the same trailing tab-delimited shape
// upstream prints: chrom, start, end, then any populated optional columns.
// We deliberately do not promote zero-valued optional fields to a printed "0"
// — matching `bed.Writer.Write`'s logic.
func formatBED(r *bed.Record) string {
	var sb strings.Builder
	sb.WriteString(r.Chrom)
	sb.WriteByte('\t')
	sb.WriteString(strconv.Itoa(r.ChromStart))
	sb.WriteByte('\t')
	sb.WriteString(strconv.Itoa(r.ChromEnd))
	if r.Name != "" {
		sb.WriteByte('\t')
		sb.WriteString(r.Name)
		if r.Score != 0 || r.Strand != "" {
			sb.WriteByte('\t')
			sb.WriteString(strconv.Itoa(r.Score))
			if r.Strand != "" {
				sb.WriteByte('\t')
				sb.WriteString(r.Strand)
			}
		}
	}
	for _, e := range r.ExtraFields {
		sb.WriteByte('\t')
		sb.WriteString(e)
	}
	return sb.String()
}
