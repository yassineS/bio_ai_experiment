// Package bedintersect provides functionality to find intersecting intervals between two BED files.
package bedintersect

import (
	"fmt"
	"io"
	"sort"
	"strconv"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bed"
)

// IntersectOptions contains options for finding intersections.
type IntersectOptions struct {
	MinOverlap int     // Minimum overlap required (default: 1bp)
	FractionA  float64 // Minimum fraction of A that must overlap (0.0-1.0)
	FractionB  float64 // Minimum fraction of B that must overlap (0.0-1.0)
	StrandSpec bool    // -s: only report hits on the same strand
	// ForceOpposite enables -S: only report hits on the opposite strand.
	ForceOpposite bool
	NoOverlap     bool // -v: report A entries with no overlap with B
	WriteA        bool // -wa: write original A entry (default: write intersection)
	WriteB        bool // -wb: write B entry alongside A
	Count         bool // -c: for each A, report total count of B overlaps
	// CountEach enables -C: for each A, report the overlap count with each B
	// file on a separate line.
	CountEach bool
	// Unique enables -u: report each A record once if it has any overlap.
	Unique     bool
	Reciprocal bool // -r: require reciprocal overlap (both -f and -F)
	// EitherFraction enables -e: satisfying the -f OR the -F fraction test
	// suffices (rather than requiring both).
	EitherFraction bool
	Distance       bool // Report distance to nearest B feature (bedintersect extension)
	Closest        bool // Report closest B feature for each A (bedintersect extension)
	UseTree        bool // Use interval tree for large B files

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

	// Names holds the per-B-file aliases for -names. When set (and more than one
	// B file is present), the alias is printed as the DB-id column instead of a
	// numeric file id in -wb/-loj/-wao/-C output.
	Names []string
	// FileNames enables -filenames: print each B file's name as the DB-id column
	// instead of a numeric file id.
	FileNames bool
	// FilePaths holds the B file paths in the same order as the B readers, used
	// to render the DB-id column under -filenames.
	FilePaths []string
	// SortOut enables -sortout: sort each A record's DB hits by position across
	// all B files before printing (rather than grouping them by B file).
	SortOut bool
	// Header enables -header: echo the A file's header lines (lines before the
	// first data record) verbatim ahead of the results.
	Header bool

	// Sorted enables -sorted: validate that the inputs are coordinate-sorted and
	// error otherwise, matching upstream's chromsweep precondition. It does not
	// change the output order (which already matches upstream's bin order).
	Sorted bool
	// GenomeOrder is the ordered list of chromosome names from a -g genome file.
	// When set (with -sorted), each input file's chromosomes must appear in this
	// order, else upstream errors. Empty means no genome-order constraint.
	GenomeOrder []string
	// GenomeFile is the -g genome file path, echoed verbatim in the -sorted
	// genome-order error message.
	GenomeFile string
	// NoNameCheck enables -nonamecheck: suppress the chromosome naming-convention
	// warning for sorted data. It is accepted for compatibility; this port does
	// not emit that warning regardless, so the flag is a no-op.
	NoNameCheck bool
	// NameA and NameB hold the A and (first) B file names, used only to render
	// upstream's -sorted out-of-order error messages.
	NameA string
	NameB string
	// Warnings, when non-nil, receives upstream-style stderr warnings (currently
	// the chromosome naming-convention warning). It is separate from the output
	// writer so warnings never pollute the data stream.
	Warnings io.Writer
}

// multiDB reports whether the DB-id column must be emitted, i.e. whether there
// is more than one B file, matching upstream's behaviour of prefixing each B
// record with a file id / name only when multiple databases are supplied.
func (opts IntersectOptions) multiDB() bool {
	return len(opts.FilePaths) > 1
}

// dbLabel returns the DB-id column value for B file i (0-based), following
// upstream precedence: an explicit -names alias, else the file name under
// -filenames, else the 1-based numeric file id.
func (opts IntersectOptions) dbLabel(i int) string {
	if i < len(opts.Names) {
		return opts.Names[i]
	}
	if opts.FileNames && i < len(opts.FilePaths) {
		return opts.FilePaths[i]
	}
	return strconv.Itoa(i + 1)
}

// Intersect finds intervals in A that overlap with intervals in B. It is the
// single-B-file entry point; IntersectMulti generalises it to multiple B files.
func Intersect(readerA, readerB io.Reader, writer io.Writer, opts IntersectOptions) (int, error) {
	return IntersectMulti(readerA, []io.Reader{readerB}, writer, opts)
}

// IntersectMulti finds intervals in A that overlap intervals in any of the B
// readers. Multiple B files reproduce upstream `bedtools intersect -b f1 f2 ...`
// semantics: each B record is tagged with its file id (0-based) so the DB-id
// column (-names/-filenames/numeric) and per-file -C counts can be emitted.
func IntersectMulti(readerA io.Reader, readersB []io.Reader, writer io.Writer, opts IntersectOptions) (int, error) {
	// The bedintersect-only distance/closest extensions (-d/-k) keep the legacy
	// typed bed.Record path and only support a single B file.
	if opts.Distance || opts.Closest {
		return intersectClosest(readerA, readersB[0], writer, opts)
	}
	// All upstream-parity output modes use the raw column-preserving path so
	// input columns echo verbatim and BAM/VCF/GFF inputs are supported. The
	// per-A-record output shape (default / -wa / -wb / -c / -C / -u / -v / join
	// modes) is selected inside intersectRaw.
	return intersectRaw(readerA, readersB, writer, opts)
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
		return intersectRawWithStats(readerA, []io.Reader{readerB}, writer, opts)
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

// strandOK applies the strand filter for one A/B pair. Without -s or -S any
// pair passes. With -s only same-strand pairs pass (sameStrandMatch); with -S
// only opposite-strand pairs pass. For -S, as for -s, ".", "*" and a missing
// strand are UNKNOWN and can never match (both records must carry a real "+"
// or "-" strand), mirroring upstream Record::sameChromIntersects when
// useDiffStrands is set.
func strandOK(opts IntersectOptions, a, b string) bool {
	switch {
	case opts.StrandSpec:
		return sameStrandMatch(a, b)
	case opts.ForceOpposite:
		if (a != "+" && a != "-") || (b != "+" && b != "-") {
			return false
		}
		return a != b
	default:
		return true
	}
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
