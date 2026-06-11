// Raw, column-preserving implementation of the upstream-parity intersect
// output modes (default intersection, -wa, -wb, -c, -v). Unlike the legacy
// typed bed.Record path, this echoes A's (and B's, for -wb) original input
// columns verbatim and re-encodes clipped coordinates per the source format,
// matching `bedtools intersect` byte-for-byte. It also supports BED/VCF/GFF/BAM
// inputs via readInRecords.
package bedintersect

import (
	"bufio"
	"fmt"
	"io"
	"strconv"
)

// intersectRaw runs the default / -wa / -wb / -c / -v output modes over the
// raw inRecord model. Returns the number of output lines written.
func intersectRaw(readerA, readerB io.Reader, writer io.Writer, opts IntersectOptions) (int, error) {
	bRecords, err := readInRecords(readerB)
	if err != nil {
		return 0, fmt.Errorf("error reading B intervals: %w", err)
	}
	// B is indexed by chromosome but kept in original file order within each
	// chromosome: upstream echoes overlapping B records in the order they appear
	// in the B file (not sorted), so the output matches byte-for-byte.
	byChrom := make(map[string][]*inRecord)
	for _, b := range bRecords {
		byChrom[b.chrom] = append(byChrom[b.chrom], b)
	}

	aRecords, err := readInRecords(readerA)
	if err != nil {
		return 0, fmt.Errorf("error reading A intervals: %w", err)
	}

	bw := bufio.NewWriter(writer)
	count := 0
	emit := func(s string) error {
		if _, err := bw.WriteString(s); err != nil {
			return fmt.Errorf("error writing result: %w", err)
		}
		if err := bw.WriteByte('\n'); err != nil {
			return fmt.Errorf("error writing result: %w", err)
		}
		return nil
	}

	for _, a := range aRecords {
		hits := rawOverlaps(a, byChrom[a.chrom], opts)

		switch {
		case opts.NoOverlap: // -v: emit A (verbatim) only when there are no hits
			if len(hits) == 0 {
				if err := emit(a.line); err != nil {
					return count, err
				}
				count++
			}
		case opts.Count: // -c: A (verbatim) + the overlap count as a trailing column
			if err := emit(a.line + "\t" + strconv.Itoa(len(hits))); err != nil {
				return count, err
			}
			count++
		default:
			for _, h := range hits {
				var out string
				switch {
				case opts.WriteA && opts.WriteB:
					out = a.line + "\t" + h.b.line
				case opts.WriteB:
					// -wb alone: A clipped to the overlap, then full original B.
					out = a.clippedLine(h.start, h.end) + "\t" + h.b.line
				case opts.WriteA:
					out = a.line
				default:
					// Default: A clipped to the overlap span, columns verbatim.
					out = a.clippedLine(h.start, h.end)
				}
				if err := emit(out); err != nil {
					return count, err
				}
				count++
			}
		}
	}
	if err := bw.Flush(); err != nil {
		return count, fmt.Errorf("error flushing output: %w", err)
	}
	return count, nil
}

// intersectRawWithStats is intersectRaw plus per-A hit accounting for the -S
// stats summary. It supports the same output modes (default / -wa / -wb / -c /
// -v) over the raw inRecord model.
func intersectRawWithStats(readerA, readerB io.Reader, writer io.Writer, opts IntersectOptions) (*Stats, error) {
	bRecords, err := readInRecords(readerB)
	if err != nil {
		return nil, fmt.Errorf("error reading B intervals: %w", err)
	}
	byChrom := make(map[string][]*inRecord)
	for _, b := range bRecords {
		byChrom[b.chrom] = append(byChrom[b.chrom], b)
	}
	aRecords, err := readInRecords(readerA)
	if err != nil {
		return nil, fmt.Errorf("error reading A intervals: %w", err)
	}

	stats := &Stats{IntervalsB: len(bRecords)}
	bw := bufio.NewWriter(writer)
	emit := func(s string) error {
		if _, err := bw.WriteString(s); err != nil {
			return fmt.Errorf("error writing result: %w", err)
		}
		if err := bw.WriteByte('\n'); err != nil {
			return fmt.Errorf("error writing result: %w", err)
		}
		return nil
	}

	for _, a := range aRecords {
		stats.IntervalsA++
		hits := rawOverlaps(a, byChrom[a.chrom], opts)
		if len(hits) > 0 {
			stats.IntervalsAHit++
			stats.Overlaps += len(hits)
		} else {
			stats.IntervalsAMiss++
		}

		switch {
		case opts.NoOverlap:
			if len(hits) == 0 {
				if err := emit(a.line); err != nil {
					return stats, err
				}
			}
		case opts.Count:
			if err := emit(a.line + "\t" + strconv.Itoa(len(hits))); err != nil {
				return stats, err
			}
		default:
			for _, h := range hits {
				var out string
				switch {
				case opts.WriteA && opts.WriteB:
					out = a.line + "\t" + h.b.line
				case opts.WriteB:
					out = a.clippedLine(h.start, h.end) + "\t" + h.b.line
				case opts.WriteA:
					out = a.line
				default:
					out = a.clippedLine(h.start, h.end)
				}
				if err := emit(out); err != nil {
					return stats, err
				}
			}
		}
	}
	if err := bw.Flush(); err != nil {
		return stats, fmt.Errorf("error flushing output: %w", err)
	}
	return stats, nil
}

// rawHit records a B record overlapping A, with the 0-based overlap span used
// for the default-mode clip.
type rawHit struct {
	b     *inRecord
	start int
	end   int
}

// rawOverlaps returns the B records overlapping A (in B position order),
// applying the -s strand filter and the -f/-F/-r fraction tests with
// split-aware block math when opts.Split is set. The overlap span carried back
// is the whole-record intersection (used to clip A in default mode); under
// -split the fraction tests use the non-redundant block overlap but the clip
// span is still the whole-span intersection, matching upstream.
func rawOverlaps(a *inRecord, bRecords []*inRecord, opts IntersectOptions) []rawHit {
	var hits []rawHit
	aLen := a.end - a.start
	for _, b := range bRecords {
		if a.chrom != b.chrom {
			continue
		}
		if opts.StrandSpec && !sameStrandMatch(a.strand, b.strand) {
			continue
		}
		overlapStart := max(a.start, b.start)
		overlapEnd := min(a.end, b.end)
		overlapLen := overlapEnd - overlapStart
		if overlapLen <= 0 {
			continue
		}
		if overlapLen < opts.MinOverlap {
			continue
		}

		fracOverlap := overlapLen
		lenA := aLen
		lenB := b.end - b.start
		if opts.Split {
			blockOverlap, aSum, bSum := splitBlockOverlapIn(a, b)
			if blockOverlap <= 0 {
				continue
			}
			fracOverlap = blockOverlap
			lenA = aSum
			lenB = bSum
		}
		if opts.FractionA > 0 && fraction(fracOverlap, lenA) < opts.FractionA {
			continue
		}
		if opts.FractionB > 0 && fraction(fracOverlap, lenB) < opts.FractionB {
			continue
		}
		hits = append(hits, rawHit{b: b, start: overlapStart, end: overlapEnd})
	}
	return hits
}

// splitBlockOverlapIn mirrors splitBlockOverlap (the typed-record version) but
// operates on inRecords, expanding BED12 / BAM block columns into absolute
// intervals and returning the non-redundant block overlap and the summed block
// lengths used by the fraction tests.
func splitBlockOverlapIn(a, b *inRecord) (overlap, aSum, bSum int) {
	aBlocks := inBlocks(a)
	bBlocks := inBlocks(b)
	var ovs []block
	for _, hb := range bBlocks {
		for _, kb := range aBlocks {
			s := max(kb.start, hb.start)
			e := min(kb.end, hb.end)
			if e > s {
				ovs = append(ovs, block{s, e})
			}
		}
	}
	return nonRedundantOverlap(ovs), blockSum(aBlocks), blockSum(bBlocks)
}

// inBlocks expands an inRecord into absolute block intervals. BED12 and BAM
// records carry blockCount / blockSizes / blockStarts columns (fields 9..11);
// any other record (or a malformed BED12) is treated as a single whole-span
// block.
func inBlocks(r *inRecord) []block {
	switch r.format {
	case fmtBED, fmtBAM:
		if blks, ok := bed12BlocksFromFields(r.fields, r.start); ok {
			return blks
		}
	}
	// Single whole-span block, with the zero-length [p,p] -> [p-1,p+1] expansion
	// so a zero-length record can still intersect under -split (matching
	// blocksOf in the join path and upstream's adjustZeroLength).
	s, e, _ := effectiveBounds(r.start, r.end)
	return []block{{s, e}}
}

// bed12BlocksFromFields expands a BED12 record's blocks into absolute
// [start,end) ranges from its raw fields, returning ok=false when the record is
// not parseable as BED12 (fewer than 12 columns or malformed block columns).
func bed12BlocksFromFields(fields []string, recStart int) ([]block, bool) {
	if len(fields) < 12 {
		return nil, false
	}
	blockCount, err := strconv.Atoi(fields[9])
	if err != nil || blockCount <= 0 {
		return nil, false
	}
	sizes := splitCSV(fields[10])
	starts := splitCSV(fields[11])
	if len(sizes) != blockCount || len(starts) != blockCount {
		return nil, false
	}
	blks := make([]block, 0, blockCount)
	for i := 0; i < blockCount; i++ {
		off, err := strconv.Atoi(starts[i])
		if err != nil {
			return nil, false
		}
		size, err := strconv.Atoi(sizes[i])
		if err != nil {
			return nil, false
		}
		s := recStart + off
		blks = append(blks, block{s, s + size})
	}
	return blks, true
}
