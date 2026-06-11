// Package bedintersect provides functionality to find intersecting intervals between two BED files.
package bedintersect

import (
	"fmt"
	"io"
	"sort"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bed"
)

// IntersectOptions contains options for finding intersections.
type IntersectOptions struct {
	MinOverlap int     // Minimum overlap required (default: 1bp)
	FractionA  float64 // Minimum fraction of A that must overlap (0.0-1.0)
	FractionB  float64 // Minimum fraction of B that must overlap (0.0-1.0)
	StrandSpec bool    // Only report hits on same strand
	NoOverlap  bool    // Report A entries with no overlap with B
	WriteA     bool    // Write original A entry (default: write intersection)
	WriteB     bool    // Write B entry instead of A
	Count      bool    // For each A, report count of B overlaps
	Reciprocal bool    // Require reciprocal overlap (both -f and -F must be satisfied)
	Distance   bool    // Report distance to nearest B feature
	Closest    bool    // Report closest B feature for each A
	UseTree    bool    // Use interval tree for large B files

	// LeftJoin enables -loj: report every A record, appending the overlapping
	// B record (or a null B placeholder when there are no overlaps).
	LeftJoin bool
	// WriteOverlap enables -wo: report A and B for each overlap, followed by
	// the number of overlapping bases.
	WriteOverlap bool
	// WriteAllOverlap enables -wao: like WriteOverlap, but also report A with a
	// null B and an overlap count of 0 when A has no overlaps.
	WriteAllOverlap bool
	// Split enables -split: treat BED12 blocks as separate intervals when
	// determining overlaps and counting overlapping bases.
	Split bool
}

// Intersect finds intervals in A that overlap with intervals in B.
func Intersect(readerA, readerB io.Reader, writer io.Writer, opts IntersectOptions) (int, error) {
	// The join/overlap output modes (-loj, -wo, -wao, -wa -wb) echo A and B
	// columns verbatim and in B-file order, so they use a separate raw,
	// line-preserving code path.
	if opts.usesJoinMode() {
		return intersectJoin(readerA, readerB, writer, opts)
	}
	// Every upstream-parity output mode (default intersection, -wa, -wb, -c, -v)
	// uses the raw column-preserving path so input columns echo verbatim and
	// BAM/VCF/GFF inputs are supported; opts.UseTree selects the interval-tree
	// index there. Only the bedintersect-only distance/closest extensions
	// (-d/-k) keep the legacy typed bed.Record path below.
	if !opts.Distance && !opts.Closest {
		return intersectRaw(readerA, readerB, writer, opts)
	}
	return intersectClosest(readerA, readerB, writer, opts)
}

// intersectClosest implements the bedintersect-only -d/-k distance/closest
// extensions over the typed bed.Record model. It is not part of the
// upstream-parity surface (upstream `bedtools intersect` has no distance mode);
// the dedicated `bedclosest` tool covers `bedtools closest`.
func intersectClosest(readerA, readerB io.Reader, writer io.Writer, opts IntersectOptions) (int, error) {
	chromIndex, err := readChromIndex(readerB)
	if err != nil {
		return 0, err
	}

	bedReaderA := bed.NewReader(readerA)
	bedWriter := bed.NewWriter(writer)
	count := 0

	for {
		recordA, err := bedReaderA.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return 0, fmt.Errorf("error reading A intervals: %w", err)
		}

		closest, dist := findClosest(recordA, chromIndex[recordA.Chrom], opts)
		result, write := closestResult(recordA, closest, dist, opts)
		if write {
			if err := bedWriter.Write(result); err != nil {
				return 0, fmt.Errorf("error writing result: %w", err)
			}
			count++
		}
	}

	if err := bedWriter.Flush(); err != nil {
		return 0, fmt.Errorf("error flushing output: %w", err)
	}
	return count, nil
}

// readChromIndex reads every B record and buckets it by chromosome, sorted by
// start, for the typed distance/closest path.
func readChromIndex(readerB io.Reader) (map[string][]*bed.Record, error) {
	bedReaderB := bed.NewReader(readerB)
	var intervalsB []*bed.Record
	for {
		record, err := bedReaderB.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, fmt.Errorf("error reading B intervals: %w", err)
		}
		intervalsB = append(intervalsB, record)
	}
	sort.Slice(intervalsB, func(i, j int) bool {
		if intervalsB[i].Chrom != intervalsB[j].Chrom {
			return intervalsB[i].Chrom < intervalsB[j].Chrom
		}
		return intervalsB[i].ChromStart < intervalsB[j].ChromStart
	})
	chromIndex := make(map[string][]*bed.Record)
	for _, interval := range intervalsB {
		chromIndex[interval.Chrom] = append(chromIndex[interval.Chrom], interval)
	}
	return chromIndex, nil
}

// closestResult builds the output record for the -d/-k path and reports whether
// it should be written. -d emits A's coordinates with the distance in the name
// field (or -1 when no B is on the chromosome); -k emits the closest B record
// (and nothing when there is none).
func closestResult(recordA, closest *bed.Record, dist int, opts IntersectOptions) (*bed.Record, bool) {
	if closest == nil {
		if opts.Distance {
			return &bed.Record{
				Chrom:      recordA.Chrom,
				ChromStart: recordA.ChromStart,
				ChromEnd:   recordA.ChromEnd,
				Name:       "-1",
			}, true
		}
		return nil, false
	}
	if opts.Distance {
		return &bed.Record{
			Chrom:      recordA.Chrom,
			ChromStart: recordA.ChromStart,
			ChromEnd:   recordA.ChromEnd,
			Name:       fmt.Sprintf("%d", dist),
		}, true
	}
	// opts.Closest
	return closest, true
}

// findClosest finds the closest interval in B to A and returns the distance.
// Returns the closest interval and distance (0 if overlapping, positive if upstream/downstream).
func findClosest(a *bed.Record, bIntervals []*bed.Record, opts IntersectOptions) (*bed.Record, int) {
	var closest *bed.Record
	minDist := -1

	for _, b := range bIntervals {
		// Skip if chromosomes don't match
		if a.Chrom != b.Chrom {
			continue
		}

		// Check strand if required
		if opts.StrandSpec && !sameStrandMatch(a.Strand, b.Strand) {
			continue
		}

		// Calculate distance
		var dist int
		if a.ChromEnd <= b.ChromStart {
			// A is upstream of B
			dist = b.ChromStart - a.ChromEnd
		} else if b.ChromEnd <= a.ChromStart {
			// B is upstream of A
			dist = a.ChromStart - b.ChromEnd
		} else {
			// Overlapping
			dist = 0
		}

		// Update closest if this is closer
		if minDist == -1 || dist < minDist {
			minDist = dist
			closest = b
		}
	}

	return closest, minDist
}

// Stats contains statistics about the intersect operation.
type Stats struct {
	IntervalsA     int
	IntervalsB     int
	Overlaps       int
	IntervalsAHit  int
	IntervalsAMiss int
}

// IntersectWithStats performs intersection and returns detailed statistics.
func IntersectWithStats(readerA, readerB io.Reader, writer io.Writer, opts IntersectOptions) (*Stats, error) {
	// Every upstream-parity output mode uses the raw, column-preserving path so
	// the output matches `bedtools intersect` byte-for-byte (and supports
	// BAM/VCF/GFF inputs). Only the bedintersect-only -d/-k extensions fall
	// through to the legacy typed path below.
	if !opts.Distance && !opts.Closest {
		return intersectRawWithStats(readerA, readerB, writer, opts)
	}
	return intersectClosestWithStats(readerA, readerB, writer, opts)
}

// intersectClosestWithStats is intersectClosest plus per-A hit accounting for
// the -S stats summary, over the typed bed.Record -d/-k path.
func intersectClosestWithStats(readerA, readerB io.Reader, writer io.Writer, opts IntersectOptions) (*Stats, error) {
	chromIndex, err := readChromIndex(readerB)
	if err != nil {
		return nil, err
	}
	intervalsB := 0
	for _, recs := range chromIndex {
		intervalsB += len(recs)
	}

	stats := &Stats{IntervalsB: intervalsB}
	bedReaderA := bed.NewReader(readerA)
	bedWriter := bed.NewWriter(writer)

	for {
		recordA, err := bedReaderA.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return nil, fmt.Errorf("error reading A intervals: %w", err)
		}

		stats.IntervalsA++
		closest, dist := findClosest(recordA, chromIndex[recordA.Chrom], opts)
		if closest != nil {
			stats.IntervalsAHit++
		} else {
			stats.IntervalsAMiss++
		}
		result, write := closestResult(recordA, closest, dist, opts)
		if write {
			if err := bedWriter.Write(result); err != nil {
				return nil, fmt.Errorf("error writing result: %w", err)
			}
		}
	}

	if err := bedWriter.Flush(); err != nil {
		return nil, fmt.Errorf("error flushing output: %w", err)
	}
	return stats, nil
}

// sameStrandMatch reports whether two strand columns count as the same strand
// under `-s`, matching upstream Record::sameChromIntersects. Only "+" and "-"
// are real strands; ".", "*", and a missing column are all UNKNOWN and can
// never satisfy a same-strand requirement, so a hit is reported only when both
// strands are a known value and equal.
func sameStrandMatch(a, b string) bool {
	if a != "+" && a != "-" {
		return false
	}
	return a == b
}

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}
