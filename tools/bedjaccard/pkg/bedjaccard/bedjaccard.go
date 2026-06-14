// Package bedjaccard computes the Jaccard similarity between two sorted BED
// files, mirroring `bedtools jaccard`.
//
// Given two BED files A and B that are pre-sorted by (chrom, start), it
// computes:
//
//   - intersection: total bases shared between A and B
//   - union:        |A| + |B| - intersection
//   - jaccard:      intersection / union (0 if union is zero)
//   - n:            number of (A, B) interval pairs that overlap
//
// The algorithm performs a single linear sweep, holding only a small "active"
// window of B intervals at any time, so memory use is independent of file
// size beyond the sort prerequisite.
package bedjaccard

import (
	"errors"
	"fmt"
	"io"
	"sort"
	"strconv"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bed"
)

// Options configures the Jaccard computation.
type Options struct {
	// SameStrand ("-s") considers only same-strand A/B pairs; records with an
	// unknown ("." or empty) strand are dropped, matching upstream's
	// SAME_STRAND_EITHER merge mode.
	SameStrand bool
	// StrandFilter ("-S <+|->") restricts BOTH inputs to records on the given
	// strand before the merge/intersection, matching upstream's
	// SAME_STRAND_FORWARD / SAME_STRAND_REVERSE modes. It is mutually
	// exclusive with SameStrand. (Upstream jaccard has no opposite-strand
	// mode; -S takes a strand argument.)
	StrandFilter string
	// FractionA: require at least this fraction of A to overlap B for the
	// pair to count (0..1). Zero disables the check.
	FractionA float64
	// FractionB: require at least this fraction of B to overlap A.
	FractionB float64
	// Split ("-split") treats BED12 records as their individual blocks
	// (exon-aware): each record is expanded into one interval per block
	// before merging/intersection, matching upstream bedtools jaccard -split.
	Split bool
}

// Result is the one-line summary written by Run.
type Result struct {
	Intersection int
	Union        int
	Jaccard      float64
	N            int
}

// Run reads sorted BED records from a and b, computes the Jaccard summary, and
// writes a two-line tab-separated table (header then values) to w. The Result
// is also returned for programmatic use.
func Run(a, b io.Reader, w io.Writer, opts Options) (*Result, error) {
	if opts.SameStrand && opts.StrandFilter != "" {
		return nil, errors.New("-s and -S are mutually exclusive")
	}
	if opts.StrandFilter != "" && opts.StrandFilter != "+" && opts.StrandFilter != "-" {
		return nil, fmt.Errorf("-S must be + or -, got %q", opts.StrandFilter)
	}
	if opts.FractionA < 0 || opts.FractionA > 1 {
		return nil, fmt.Errorf("-f must be in [0,1], got %v", opts.FractionA)
	}
	if opts.FractionB < 0 || opts.FractionB > 1 {
		return nil, fmt.Errorf("-F must be in [0,1], got %v", opts.FractionB)
	}

	res, err := jaccard(a, b, opts)
	if err != nil {
		return nil, err
	}

	if _, err := fmt.Fprintln(w, "intersection\tunion\tjaccard\tn_intersections"); err != nil {
		return nil, err
	}
	if _, err := fmt.Fprintf(w, "%d\t%d\t%s\t%d\n",
		res.Intersection, res.Union, formatJaccard(res.Jaccard), res.N); err != nil {
		return nil, err
	}
	return res, nil
}

// formatJaccard renders the ratio with C++ ostream's default precision
// (6 significant digits with %g-style trimming), which is what upstream
// `bedtools jaccard` uses when it prints the ratio via `cout`.
func formatJaccard(j float64) string {
	return strconv.FormatFloat(j, 'g', 6, 64)
}

// jaccard does the streaming sweep.
//
// Upstream `bedtools jaccard` calls `setUseMergedIntervals(true)` on its
// context (see reference_code/bedtools/src/utils/Contexts/ContextJaccard.cpp),
// which means both A and B are streamed through a per-chromosome merge
// before intersection/union are computed. Without that step our counts
// (and `n_intersections`) diverge from upstream whenever either side
// contains overlapping or adjacent records. We replicate the pre-merge
// by wrapping each input reader in a `mergingReader` that emits merged
// records on the fly. When `-s` (`SameStrand`) or `-S` (`OppositeStrand`)
// is in play, the merge runs per-strand so cross-strand records don't
// collapse into one.
func jaccard(aReader, bReader io.Reader, opts Options) (*Result, error) {
	// -S restricts both inputs to one strand before merging; -s merges
	// per-strand and keeps only same-strand pairs.
	perStrand := opts.SameStrand || opts.StrandFilter != ""
	var aSrc, bSrc recordReader = bed.NewReader(aReader), bed.NewReader(bReader)
	if opts.Split {
		// -split expands each BED12 record into its blocks before the sweep.
		aSrc = &blockSplitReader{in: aSrc}
		bSrc = &blockSplitReader{in: bSrc}
	}
	ra := newMergingReader(aSrc, perStrand, opts.StrandFilter)
	rb := newMergingReader(bSrc, perStrand, opts.StrandFilter)

	var active []*bed.Record
	var (
		lastA, lastB *bed.Record
		bExhausted   bool

		totalA, totalB, totalIntersect, totalPairs int
	)

	for {
		recA, err := ra.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, fmt.Errorf("reading A: %w", err)
		}
		if lastA != nil && !sortedAfter(lastA, recA) {
			return nil, fmt.Errorf("input A is not sorted: %s:%d..%d before %s:%d..%d",
				lastA.Chrom, lastA.ChromStart, lastA.ChromEnd,
				recA.Chrom, recA.ChromStart, recA.ChromEnd)
		}
		lastA = recA
		totalA += recA.ChromEnd - recA.ChromStart

		// Pull more B until the active window covers recA's territory.
		for !bExhausted && needMoreB(active, recA) {
			recB, err := rb.Read()
			if err == io.EOF {
				bExhausted = true
				break
			}
			if err != nil {
				return nil, fmt.Errorf("reading B: %w", err)
			}
			if lastB != nil && !sortedAfter(lastB, recB) {
				return nil, fmt.Errorf("input B is not sorted: %s:%d..%d before %s:%d..%d",
					lastB.Chrom, lastB.ChromStart, lastB.ChromEnd,
					recB.Chrom, recB.ChromStart, recB.ChromEnd)
			}
			lastB = recB
			totalB += recB.ChromEnd - recB.ChromStart
			active = append(active, recB)
		}

		// Drop B intervals that can no longer match recA or any later A.
		active = pruneActive(active, recA)

		// Score recA against the current active window.
		for _, b := range active {
			if b.Chrom != recA.Chrom {
				continue
			}
			if !strandOK(recA, b, opts) {
				continue
			}
			start := recA.ChromStart
			if b.ChromStart > start {
				start = b.ChromStart
			}
			end := recA.ChromEnd
			if b.ChromEnd < end {
				end = b.ChromEnd
			}
			if end <= start {
				continue
			}
			ov := end - start
			if !fractionOK(recA, b, ov, opts) {
				continue
			}
			totalIntersect += ov
			totalPairs++
		}
	}

	// Drain any remaining B records so we account for |B|.
	for !bExhausted {
		recB, err := rb.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, fmt.Errorf("reading B: %w", err)
		}
		if lastB != nil && !sortedAfter(lastB, recB) {
			return nil, fmt.Errorf("input B is not sorted: %s:%d..%d before %s:%d..%d",
				lastB.Chrom, lastB.ChromStart, lastB.ChromEnd,
				recB.Chrom, recB.ChromStart, recB.ChromEnd)
		}
		lastB = recB
		totalB += recB.ChromEnd - recB.ChromStart
	}

	union := totalA + totalB - totalIntersect
	jacc := 0.0
	if union > 0 {
		jacc = float64(totalIntersect) / float64(union)
	}
	return &Result{
		Intersection: totalIntersect,
		Union:        union,
		Jaccard:      jacc,
		N:            totalPairs,
	}, nil
}

// needMoreB returns true if more B records should be pulled to fully cover
// candidates for recA: when active is empty, the last B is on an earlier chrom,
// or the last B starts before recA ends on the same chrom.
func needMoreB(active []*bed.Record, recA *bed.Record) bool {
	if len(active) == 0 {
		return true
	}
	last := active[len(active)-1]
	if last.Chrom != recA.Chrom {
		return last.Chrom < recA.Chrom
	}
	return last.ChromStart < recA.ChromEnd
}

// pruneActive removes B intervals that can no longer overlap recA or any
// subsequent A (A is sorted).
func pruneActive(active []*bed.Record, recA *bed.Record) []*bed.Record {
	out := active[:0]
	for _, b := range active {
		if b.Chrom < recA.Chrom {
			continue
		}
		if b.Chrom == recA.Chrom && b.ChromEnd <= recA.ChromStart {
			continue
		}
		out = append(out, b)
	}
	return out
}

// sortedAfter checks whether next comes at or after prev in (chrom, start) order.
func sortedAfter(prev, next *bed.Record) bool {
	if prev.Chrom != next.Chrom {
		return prev.Chrom < next.Chrom
	}
	return next.ChromStart >= prev.ChromStart
}

// strandOK applies the per-pair -s filter. With -s only same-strand pairs
// count (an unknown "." or empty strand on either side never matches). The
// -S single-strand filter is applied earlier (records on the other strand
// are dropped before merging), so no per-pair check is needed for it.
func strandOK(a, b *bed.Record, opts Options) bool {
	if opts.SameStrand {
		if a.Strand == "" || a.Strand == "." || b.Strand == "" || b.Strand == "." {
			return false
		}
		return a.Strand == b.Strand
	}
	return true
}

// fractionOK applies -f / -F overlap-fraction filters.
func fractionOK(a, b *bed.Record, overlap int, opts Options) bool {
	if opts.FractionA > 0 {
		lenA := a.ChromEnd - a.ChromStart
		if lenA == 0 || float64(overlap)/float64(lenA) < opts.FractionA {
			return false
		}
	}
	if opts.FractionB > 0 {
		lenB := b.ChromEnd - b.ChromStart
		if lenB == 0 || float64(overlap)/float64(lenB) < opts.FractionB {
			return false
		}
	}
	return true
}

// mergingReader wraps a *bed.Reader and emits merged records on the fly:
// adjacent or overlapping records on the same chromosome are coalesced
// into a single output record before being returned by Read. This
// mirrors `bedtools jaccard`'s upstream pre-merge step (which is enabled
// in `ContextJaccard` via `setUseMergedIntervals(true)`).
//
// When perStrand is true, the merge runs per-strand: records with
// different strand values do NOT merge, even if they overlap. This
// matches upstream's behaviour under `-s`. With perStrand=false, strand
// is ignored for merging.
//
// The reader buffers any pending records (at most one per strand bucket)
// and replays them in input order; it also remembers a single
// "lookahead" record per stream so that the next Read can complete the
// in-progress merge.
// recordReader is the minimal record source the merging reader consumes;
// both *bed.Reader and *blockSplitReader satisfy it.
type recordReader interface {
	Read() (*bed.Record, error)
}

// blockSplitReader wraps a recordReader and, when a BED12 record carries
// block information, emits one sub-record per block (the exon-aware -split
// behaviour). Records without blocks pass through unchanged. Blocks are
// emitted in their input order, which is ascending by start within a record
// (BlockStarts is non-decreasing).
type blockSplitReader struct {
	in     recordReader
	queued []*bed.Record
}

// Read returns the next (possibly block-split) record.
func (s *blockSplitReader) Read() (*bed.Record, error) {
	for {
		if len(s.queued) > 0 {
			out := s.queued[0]
			s.queued = s.queued[1:]
			return out, nil
		}
		rec, err := s.in.Read()
		if err != nil {
			return nil, err
		}
		if rec.BlockCount <= 0 || len(rec.BlockSizes) == 0 {
			return rec, nil
		}
		for i := 0; i < len(rec.BlockSizes); i++ {
			start := rec.ChromStart
			if i < len(rec.BlockStarts) {
				start += rec.BlockStarts[i]
			}
			end := start + rec.BlockSizes[i]
			block := *rec
			block.ChromStart = start
			block.ChromEnd = end
			block.BlockCount = 0
			block.BlockSizes = nil
			block.BlockStarts = nil
			s.queued = append(s.queued, &block)
		}
	}
}

type mergingReader struct {
	in        recordReader
	perStrand bool
	// strandFilter, when non-empty ("+" or "-"), drops every raw record on a
	// different strand before merging (upstream -S single-strand filter).
	strandFilter string

	// Pending merges, keyed by strand bucket. When perStrand is false the
	// only key in use is "" (single bucket). For perStrand we use the
	// raw strand string ("+", "-", "." or "").
	pending map[string]*bed.Record

	// queued holds completed merged records waiting to be emitted in
	// input order (FIFO). Read drains this first.
	queued []*bed.Record

	// done signals the underlying reader hit EOF and pending should be
	// flushed before returning io.EOF to the caller.
	done bool

	// lastIn tracks the last RAW record pulled from the underlying
	// reader so that we can enforce sort order on the original input
	// (not on our merged output, which is sorted by construction).
	lastIn *bed.Record
}

func newMergingReader(r recordReader, perStrand bool, strandFilter string) *mergingReader {
	return &mergingReader{
		in:           r,
		perStrand:    perStrand,
		strandFilter: strandFilter,
		pending:      make(map[string]*bed.Record),
	}
}

// Read returns the next merged record. It blocks on the underlying
// reader as needed to complete in-progress merges, and returns io.EOF
// only after every buffered/pending merged record has been delivered.
func (m *mergingReader) Read() (*bed.Record, error) {
	for {
		// Drain anything already queued first so output order matches
		// the input stream's chrom/start order.
		if len(m.queued) > 0 {
			out := m.queued[0]
			m.queued = m.queued[1:]
			return out, nil
		}
		if m.done {
			// Flush any final pending merges, in (chrom, start) order so
			// downstream sweep sortedness checks still hold.
			if len(m.pending) > 0 {
				out := m.flushPending()
				if len(out) == 0 {
					return nil, io.EOF
				}
				m.queued = out
				continue
			}
			return nil, io.EOF
		}

		rec, err := m.in.Read()
		if err == io.EOF {
			m.done = true
			continue
		}
		if err != nil {
			return nil, err
		}
		// Enforce sort order on the underlying raw stream. Our merged
		// output is sorted by construction, but we still want to surface
		// upstream-style "input is not sorted" diagnostics for the
		// caller's benefit.
		if m.lastIn != nil && !sortedAfter(m.lastIn, rec) {
			return nil, fmt.Errorf("input is not sorted: %s:%d..%d before %s:%d..%d",
				m.lastIn.Chrom, m.lastIn.ChromStart, m.lastIn.ChromEnd,
				rec.Chrom, rec.ChromStart, rec.ChromEnd)
		}
		m.lastIn = rec

		// -S single-strand filter: drop every record not on the wanted
		// strand before it can be merged or counted.
		if m.strandFilter != "" && rec.Strand != m.strandFilter {
			continue
		}

		// Under a stranded merge, upstream's FileRecordMergeMgr drops
		// any record with an UNKNOWN (".") or missing strand (see
		// reference_code/bedtools/src/utils/FileRecordTools/FileRecordMergeMgr.cpp
		// lines 47-58 + 96-129). bedtools jaccard is one of the tools
		// that drives that mgr with `-s` (SAME_STRAND_EITHER). Replicate
		// that behaviour here so |A|/|B| match upstream byte-for-byte.
		if m.perStrand {
			if rec.Strand == "" || rec.Strand == "." {
				continue
			}
		}

		key := m.strandKey(rec)
		cur, ok := m.pending[key]
		if !ok {
			m.pending[key] = rec
			continue
		}
		// New record on a different chrom from the pending one? Flush
		// EVERY pending bucket because the chrom changed in input — by
		// sortedness no later record can extend any old bucket.
		if rec.Chrom != cur.Chrom {
			m.queued = append(m.queued, m.flushPending()...)
			m.pending[key] = rec
			continue
		}
		// Same chrom & same strand bucket: extend if overlapping/
		// adjacent, otherwise emit cur and start a new run.
		if rec.ChromStart <= cur.ChromEnd {
			if rec.ChromEnd > cur.ChromEnd {
				cur.ChromEnd = rec.ChromEnd
			}
			continue
		}
		m.queued = append(m.queued, cur)
		m.pending[key] = rec
	}
}

// strandKey returns the bucket key for rec under the current mode.
func (m *mergingReader) strandKey(rec *bed.Record) string {
	if !m.perStrand {
		return ""
	}
	return rec.Strand
}

// flushPending returns every pending record sorted by (chrom, start, end)
// so the caller's sweep continues to see a sorted stream. The pending
// map is cleared.
func (m *mergingReader) flushPending() []*bed.Record {
	out := make([]*bed.Record, 0, len(m.pending))
	for _, r := range m.pending {
		out = append(out, r)
	}
	m.pending = make(map[string]*bed.Record)
	sort.Slice(out, func(i, j int) bool {
		if out[i].Chrom != out[j].Chrom {
			return out[i].Chrom < out[j].Chrom
		}
		if out[i].ChromStart != out[j].ChromStart {
			return out[i].ChromStart < out[j].ChromStart
		}
		return out[i].ChromEnd < out[j].ChromEnd
	})
	return out
}
