// Package bedpairtopair implements `bedtools pairtopair`: report overlaps
// between two BEDPE files. For each A pair, find B pairs where the chosen
// search type (both/either/neither/notboth) is satisfied with respect to A's
// two ends.
//
// Upstream reference:
//
//	reference_code/bedtools/src/pairToPair/pairToPair.cpp
//	reference_code/bedtools/src/pairToPair/pairToPairMain.cpp
//
// Algorithm (mirrors upstream's FindOverlaps / FindHitsOnBothEnds /
// FindHitsOnEitherEnd):
//
// For each A pair we ask, for each of A's two ends, which B records does it
// overlap — testing both B end1 and B end2 since A and B are unordered.
// Concretely we compute four hit sets (per upstream nomenclature):
//
//	hitsA1B1: A end1 vs B end1
//	hitsA1B2: A end1 vs B end2
//	hitsA2B1: A end2 vs B end1
//	hitsA2B2: A end2 vs B end2
//
// "both" mode means the same B pair must have one end hit by A1 and the
// other end hit by A2 — i.e. either (A1∩B1 ∧ A2∩B2) or (A1∩B2 ∧ A2∩B1) on
// the same B line.
//
// We optionally apply slop (-slop) to every A end before computing the
// overlaps; -ss makes slop strand-aware (extend in the direction of the
// strand only).
//
// We require min(overlap/aLen) >= MinFraction on each end-vs-end test (the
// upstream "_overlapFraction" check).

package bedpairtopair

import (
	"fmt"
	"io"
	"sort"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/bed"
)

// Type enumerates the -type modes upstream supports.
type Type string

const (
	TypeBoth    Type = "both"
	TypeNotboth Type = "notboth"
	TypeEither  Type = "either"
	TypeNeither Type = "neither"
)

// ValidTypes returns the accepted -type values.
func ValidTypes() []string {
	return []string{string(TypeBoth), string(TypeNotboth), string(TypeEither), string(TypeNeither)}
}

// Options configures Run.
type Options struct {
	// Type selects the overlap-search semantics. Default: TypeBoth.
	Type Type
	// MinFraction is the minimum fraction of A-end length that must be
	// covered by the B end to count as a hit (upstream's -f, default 1e-9).
	MinFraction float64
	// IgnoreStrand skips strand checks (upstream's -is).
	IgnoreStrand bool
	// SameStrand requires matching strands. Mutually exclusive with
	// OppositeStrand. Upstream pairtopair has neither flag explicitly, but
	// the more general bedtools convention is to expose them — we accept
	// them and they short-circuit the default-strand-enforced behaviour.
	SameStrand     bool
	OppositeStrand bool
	// Slop bp added to each A end (upstream -slop). Default 0.
	Slop int
	// StrandedSlop applies slop in the direction of the strand only
	// (upstream -ss).
	StrandedSlop bool
	// RequireDifferentNames mirrors upstream -rdn: a B pair whose name
	// matches the A pair is filtered out (used to avoid self-hits when
	// the two BEDPE files are derived from the same dataset).
	RequireDifferentNames bool
}

// Stats counts processed pairs and emitted output lines.
type Stats struct {
	APairs       int
	BPairs       int
	OutputLines  int
	EmittedPairs int
}

// Run walks A and emits matches against B per the selected Type.
func Run(a, b io.Reader, w io.Writer, opts Options) (*Stats, error) {
	if opts.Type == "" {
		opts.Type = TypeBoth
	}
	if !isValid(opts.Type) {
		return nil, fmt.Errorf("invalid -type %q (want one of %v)", opts.Type, ValidTypes())
	}
	if opts.SameStrand && opts.OppositeStrand {
		return nil, fmt.Errorf("-s and -S are mutually exclusive")
	}
	if opts.MinFraction <= 0 {
		opts.MinFraction = 1e-9
	}
	if opts.StrandedSlop && opts.Slop == 0 {
		return nil, fmt.Errorf("-ss requires -slop > 0")
	}

	// Read B fully and index it twice: once by chrom of end1, once by chrom
	// of end2. Each indexed Record carries the B pair's row index so we can
	// look up the originating BEDPE.
	bReader := bed.NewBEDPEReader(b)
	bPairs, err := bReader.ReadAll()
	if err != nil {
		return nil, fmt.Errorf("read B: %w", err)
	}
	idx1 := buildIndex(bPairs, true)
	idx2 := buildIndex(bPairs, false)

	stats := &Stats{BPairs: len(bPairs)}

	aReader := bed.NewBEDPEReader(a)
	out := &writer{w: w}

	for {
		pair, err := aReader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return stats, err
		}
		stats.APairs++

		// Apply slop to A ends.
		s1, e1 := applySlop(pair.Start1, pair.End1, pair.Strand1, opts)
		s2, e2 := applySlop(pair.Start2, pair.End2, pair.Strand2, opts)

		// Compute the four hit sets.
		hitsA1B1 := idx1.findHits(pair.Chrom1, s1, e1, pair.Strand1, pair.Name, &opts)
		hitsA1B2 := idx2.findHits(pair.Chrom1, s1, e1, pair.Strand1, pair.Name, &opts)
		hitsA2B1 := idx1.findHits(pair.Chrom2, s2, e2, pair.Strand2, pair.Name, &opts)
		hitsA2B2 := idx2.findHits(pair.Chrom2, s2, e2, pair.Strand2, pair.Name, &opts)

		emitted, lines, err := emit(out, pair, bPairs, hitsA1B1, hitsA1B2, hitsA2B1, hitsA2B2, opts.Type, opts.RequireDifferentNames)
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

func isValid(t Type) bool {
	for _, v := range ValidTypes() {
		if string(t) == v {
			return true
		}
	}
	return false
}

// applySlop returns the (start, end) for an A end after applying the slop
// configuration. Mirrors upstream behaviour:
//   - -ss with strand "+": end += slop only
//   - -ss with strand "-": start -= slop only (clamped to 0)
//   - non-stranded: start -= slop (clamped) and end += slop
func applySlop(start, end int, strand string, opts Options) (int, int) {
	if opts.Slop == 0 {
		return start, end
	}
	if opts.StrandedSlop {
		switch strand {
		case "+":
			return start, end + opts.Slop
		case "-":
			s := start - opts.Slop
			if s < 0 {
				s = 0
			}
			return s, end
		}
	}
	s := start - opts.Slop
	if s < 0 {
		s = 0
	}
	return s, end + opts.Slop
}

// pairIndex maps chrom -> sorted Records, plus an IntervalTree, with each
// Record's Score field overloaded to carry the originating B-pair line index.
type pairIndex struct {
	trees   map[string]*bed.IntervalTree
	records map[string][]*bed.Record
}

// buildIndex builds an index over the requested end of every B pair.
// When useEnd1 is true we index B's end1; otherwise B's end2.
func buildIndex(bPairs []*bed.BEDPE, useEnd1 bool) *pairIndex {
	byChrom := make(map[string][]*bed.Record)
	for i, p := range bPairs {
		chrom, start, end, strand := p.Chrom2, p.Start2, p.End2, p.Strand2
		if useEnd1 {
			chrom, start, end, strand = p.Chrom1, p.Start1, p.End1, p.Strand1
		}
		if chrom == "" || chrom == "." || start < 0 || end <= start {
			continue
		}
		r := &bed.Record{
			Chrom: chrom, ChromStart: start, ChromEnd: end,
			Strand: strand,
			Score:  i, // hijack Score to store the B-pair row index.
		}
		byChrom[chrom] = append(byChrom[chrom], r)
	}
	trees := make(map[string]*bed.IntervalTree, len(byChrom))
	for chrom, recs := range byChrom {
		sort.SliceStable(recs, func(i, j int) bool { return recs[i].ChromStart < recs[j].ChromStart })
		trees[chrom] = bed.NewIntervalTree(recs)
	}
	return &pairIndex{trees: trees, records: byChrom}
}

// findHits returns the B-pair row indices whose indexed end overlaps the
// query window, with strand checks and the min-fraction filter applied.
// aName / opts.RequireDifferentNames implements the -rdn filter.
func (px *pairIndex) findHits(chrom string, start, end int, aStrand, aName string, opts *Options) []int {
	if chrom == "" || chrom == "." || end <= start {
		return nil
	}
	tree, ok := px.trees[chrom]
	if !ok {
		return nil
	}
	q := &bed.Record{Chrom: chrom, ChromStart: start, ChromEnd: end}
	cands := tree.Query(q)
	if len(cands) == 0 {
		return nil
	}
	aLen := end - start
	var out []int
	for _, c := range cands {
		if !strandMatch(aStrand, c.Strand, opts) {
			continue
		}
		ov := minInt(end, c.ChromEnd) - maxInt(start, c.ChromStart)
		if ov <= 0 {
			continue
		}
		if float64(ov)/float64(aLen) < opts.MinFraction {
			continue
		}
		out = append(out, c.Score)
	}
	return out
}

func strandMatch(aStrand, bStrand string, opts *Options) bool {
	if opts.IgnoreStrand {
		return true
	}
	if opts.SameStrand {
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
	// Default upstream behaviour: enforce same-strand on each end-vs-end
	// match unless -is is provided. Empty/dot strands are treated as
	// wildcards.
	if aStrand == "" || bStrand == "" || aStrand == "." || bStrand == "." {
		return true
	}
	return aStrand == bStrand
}

func maxInt(a, b int) int {
	if a > b {
		return a
	}
	return b
}
func minInt(a, b int) int {
	if a < b {
		return a
	}
	return b
}

// emit dispatches by Type and writes output. Returns (emittedAny, numLines, err).
func emit(w *writer, a *bed.BEDPE, bPairs []*bed.BEDPE, h11, h12, h21, h22 []int, t Type, requireDifferentNames bool) (bool, int, error) {
	// Apply -rdn (require different names) by filtering each hit list early.
	// Upstream wires this through FindOverlapsPerBin so a B that shares its
	// name with A never even enters the hit set; filtering on entry gives
	// the same observable behaviour.
	if requireDifferentNames {
		filter := func(in []int) []int {
			out := in[:0]
			for _, i := range in {
				if applyRDN(a.Name, bPairs[i]) {
					out = append(out, i)
				}
			}
			return out
		}
		h11 = filter(h11)
		h12 = filter(h12)
		h21 = filter(h21)
		h22 = filter(h22)
	}

	// "Same B pair hit on both ends" requires that B-row appear in
	// (h11 AND h22) OR (h12 AND h21) — the two valid orientations.
	bothMatches := intersect(h11, h22)
	bothMatches = appendUnique(bothMatches, intersect(h12, h21)...)
	hasBoth := len(bothMatches) > 0

	switch t {
	case TypeBoth:
		lines := 0
		for _, bIdx := range bothMatches {
			if err := w.writePairPair(a, bPairs[bIdx]); err != nil {
				return false, lines, err
			}
			lines++
		}
		return lines > 0, lines, nil

	case TypeNotboth:
		// Emit the A pair if no B pair has both ends matched.
		if !hasBoth {
			if err := w.writePair(a); err != nil {
				return false, 0, err
			}
			return true, 1, nil
		}
		return false, 0, nil

	case TypeNeither:
		// Emit when NO end of A overlaps any end of any B.
		if len(h11)+len(h12)+len(h21)+len(h22) == 0 {
			if err := w.writePair(a); err != nil {
				return false, 0, err
			}
			return true, 1, nil
		}
		return false, 0, nil

	case TypeEither:
		// Upstream's FindHitsOnEitherEnd emits one line per matched B row:
		// both-end matches are emitted once; single-end matches are also
		// emitted once.
		emitted := make(map[int]bool)
		lines := 0
		// Pass 1: any B with both-end match.
		for _, bIdx := range bothMatches {
			if emitted[bIdx] {
				continue
			}
			emitted[bIdx] = true
			if err := w.writePairPair(a, bPairs[bIdx]); err != nil {
				return false, lines, err
			}
			lines++
		}
		// Pass 2: single-end B hits, in line-number order for determinism.
		for _, bIdx := range mergeSorted(h11, h12, h21, h22) {
			if emitted[bIdx] {
				continue
			}
			emitted[bIdx] = true
			if err := w.writePairPair(a, bPairs[bIdx]); err != nil {
				return false, lines, err
			}
			lines++
		}
		return lines > 0, lines, nil
	}
	return false, 0, fmt.Errorf("unhandled type %q", t)
}

// intersect returns sorted-unique row indices that appear in both inputs.
// Inputs may be unsorted and may contain duplicates.
func intersect(a, b []int) []int {
	if len(a) == 0 || len(b) == 0 {
		return nil
	}
	set := make(map[int]struct{}, len(a))
	for _, x := range a {
		set[x] = struct{}{}
	}
	seen := make(map[int]struct{})
	var out []int
	for _, x := range b {
		if _, ok := set[x]; ok {
			if _, dup := seen[x]; !dup {
				seen[x] = struct{}{}
				out = append(out, x)
			}
		}
	}
	sort.Ints(out)
	return out
}

// appendUnique appends bs into a, skipping any element already in a.
func appendUnique(a []int, bs ...int) []int {
	seen := make(map[int]struct{}, len(a))
	for _, x := range a {
		seen[x] = struct{}{}
	}
	for _, x := range bs {
		if _, ok := seen[x]; !ok {
			seen[x] = struct{}{}
			a = append(a, x)
		}
	}
	sort.Ints(a)
	return a
}

// mergeSorted merges any number of int slices into a single sorted-unique
// slice (using a set).
func mergeSorted(slices ...[]int) []int {
	seen := make(map[int]struct{})
	var all []int
	for _, s := range slices {
		for _, x := range s {
			if _, ok := seen[x]; !ok {
				seen[x] = struct{}{}
				all = append(all, x)
			}
		}
	}
	sort.Ints(all)
	return all
}

// writer wraps an io.Writer and tracks the first error.
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

func (w *writer) writePairPair(a, b *bed.BEDPE) error {
	if w.err != nil {
		return w.err
	}
	_, err := io.WriteString(w.w, a.String()+"\t"+b.String()+"\n")
	w.err = err
	return err
}

// applyRDN filters out B-pairs whose Name matches the A-pair's Name.
// Called from emit() when opts.RequireDifferentNames is true. Returns true if
// the B-pair should be kept (names differ).
func applyRDN(aName string, b *bed.BEDPE) bool {
	return aName != b.Name
}
