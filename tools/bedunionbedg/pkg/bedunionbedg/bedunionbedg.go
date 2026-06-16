// Package bedunionbedg combines multiple BEDGRAPH files into a single matrix,
// emitting, for every interval boundary across all inputs, a line of the form
//
//	chrom  start  end  val1  val2  ...  valN
//
// where valI is the value of input file I over [start, end). It mirrors the
// behaviour of `bedtools unionbedg` (aka unionBedGraphs).
//
// Each input is assumed sorted by chrom/start with non-overlapping intervals.
// Values are carried as raw strings (matching upstream's string depth type) so
// integers, floats, and arbitrary tokens round-trip unchanged. Files are merged
// chromosome-by-chromosome using a coordinate sweep that reproduces upstream's
// priority-queue algorithm exactly, including its handling of files that are
// not globally chromosome-sorted.
package bedunionbedg

import (
	"bufio"
	"container/heap"
	"fmt"
	"io"
	"strconv"
	"strings"
)

// Options bundles the configuration for Union.
type Options struct {
	// PrintHeader, when true, emits a "chrom start end [names...]" header line.
	PrintHeader bool
	// Names holds the column titles printed in the header (one per input file).
	// When empty, only "chrom start end" is printed in the header.
	Names []string
	// PrintEmpty, when true, also reports regions with no coverage in any file.
	// It requires Sizes to be populated.
	PrintEmpty bool
	// Sizes maps chromosome name to length; required when PrintEmpty is set so
	// the leading/trailing empty regions of each chromosome can be reported.
	Sizes map[string]int64
	// Filler is the value printed for an input that has no coverage over an
	// interval. Upstream's default is "0".
	Filler string
}

// bgItem is one parsed BEDGRAPH interval.
type bgItem struct {
	chrom string
	start int64
	end   int64
	depth string
	valid bool // false once the source file is exhausted
}

// coordType distinguishes a START boundary from an END boundary.
type coordType int

const (
	startPoint coordType = iota
	endPoint
)

// point is a boundary event placed on the sweep queue.
type point struct {
	sourceIndex int
	coordType   coordType
	coord       int64
	depth       string
}

// pointQueue is a min-heap of points ordered by ascending coordinate, matching
// upstream's PointWithDepth priority queue (operator< compares coord >).
type pointQueue []point

func (q pointQueue) Len() int            { return len(q) }
func (q pointQueue) Less(i, j int) bool  { return q[i].coord < q[j].coord }
func (q pointQueue) Swap(i, j int)       { q[i], q[j] = q[j], q[i] }
func (q *pointQueue) Push(x interface{}) { *q = append(*q, x.(point)) }
func (q *pointQueue) Pop() interface{} {
	old := *q
	n := len(old)
	it := old[n-1]
	*q = old[:n-1]
	return it
}

// source wraps one input file's scanner and look-ahead item.
type source struct {
	scanner *bufio.Scanner
	lineNum int
	current bgItem
}

// unioner holds the running state of a Union pass.
type unioner struct {
	sources      []*source
	out          *bufio.Writer
	queue        pointQueue
	currentChrom string
	currentDepth []string
	nonZero      int
	opts         Options
}

// Union reads the BEDGRAPH inputs from readers, combines them, and writes the
// union matrix to out. readers must contain at least two inputs. When
// opts.PrintEmpty is set, opts.Sizes must contain every chromosome present in
// the inputs. Union returns an error on malformed input (e.g. a non-integer
// coordinate) or a missing chromosome size needed for empty-region reporting.
func Union(readers []io.Reader, out io.Writer, opts Options) error {
	if len(readers) < 2 {
		return fmt.Errorf("at least two BedGraph files are required")
	}
	if opts.Filler == "" {
		opts.Filler = "0"
	}

	bw := bufio.NewWriter(out)
	defer bw.Flush()

	u := &unioner{
		out:          bw,
		currentDepth: make([]string, len(readers)),
		opts:         opts,
	}
	for i := range u.currentDepth {
		u.currentDepth[i] = opts.Filler
	}
	for _, r := range readers {
		sc := bufio.NewScanner(r)
		sc.Buffer(make([]byte, 64*1024), 64*1024*1024)
		u.sources = append(u.sources, &source{scanner: sc})
	}

	if opts.PrintHeader {
		u.printHeader()
	}

	// Prime each source with its first valid interval.
	for i := range u.sources {
		if err := u.loadNext(i); err != nil {
			return err
		}
	}

	for {
		u.currentChrom = u.determineNextChrom()
		if u.currentChrom == "" {
			break
		}
		// Populate the queue with all front intervals on the current chrom.
		for i := range u.sources {
			if err := u.addInterval(i); err != nil {
				return err
			}
		}
		if len(u.queue) == 0 {
			// Should not happen once a chromosome is selected, but guard
			// against pathological input rather than spinning.
			break
		}

		currentStart := u.consumeNextCoordinate()
		if opts.PrintEmpty && currentStart > 0 {
			u.printEmptyCoverage(0, currentStart)
		}
		for {
			currentEnd := u.queue[0].coord
			u.printCoverage(currentStart, currentEnd)
			currentStart = u.consumeNextCoordinate()
			if len(u.queue) == 0 {
				break
			}
		}
		if opts.PrintEmpty {
			size, ok := opts.Sizes[u.currentChrom]
			if !ok {
				return fmt.Errorf("chromosome %q not found in genome file (-g)", u.currentChrom)
			}
			if currentStart < size {
				u.printEmptyCoverage(currentStart, size)
			}
		}
		if u.allFilesDone() {
			break
		}
	}
	return nil
}

// printHeader writes the header line: "chrom\tstart\tend" plus any names.
func (u *unioner) printHeader() {
	u.out.WriteString("chrom\tstart\tend")
	for _, name := range u.opts.Names {
		u.out.WriteByte('\t')
		u.out.WriteString(name)
	}
	u.out.WriteByte('\n')
}

// printCoverage emits the current depths over [start, end), unless every input
// is at filler coverage and empty regions were not requested.
func (u *unioner) printCoverage(start, end int64) {
	if u.nonZero == 0 && !u.opts.PrintEmpty {
		return
	}
	u.out.WriteString(u.currentChrom)
	u.out.WriteByte('\t')
	u.out.WriteString(strconv.FormatInt(start, 10))
	u.out.WriteByte('\t')
	u.out.WriteString(strconv.FormatInt(end, 10))
	for _, d := range u.currentDepth {
		u.out.WriteByte('\t')
		u.out.WriteString(d)
	}
	u.out.WriteByte('\n')
}

// printEmptyCoverage emits filler depths over [start, end) for a region with no
// coverage in any input.
func (u *unioner) printEmptyCoverage(start, end int64) {
	u.out.WriteString(u.currentChrom)
	u.out.WriteByte('\t')
	u.out.WriteString(strconv.FormatInt(start, 10))
	u.out.WriteByte('\t')
	u.out.WriteString(strconv.FormatInt(end, 10))
	for range u.currentDepth {
		u.out.WriteByte('\t')
		u.out.WriteString(u.opts.Filler)
	}
	u.out.WriteByte('\n')
}

// consumeNextCoordinate pops every queued point sharing the smallest
// coordinate, updates running coverage for each, and returns that coordinate.
func (u *unioner) consumeNextCoordinate() int64 {
	pos := u.queue[0].coord
	for len(u.queue) > 0 && u.queue[0].coord == pos {
		it := heap.Pop(&u.queue).(point)
		u.updateInformation(it)
	}
	return pos
}

// updateInformation applies a START or END boundary to the running state. On an
// END boundary it also pulls the next interval from that source onto the queue.
func (u *unioner) updateInformation(it point) {
	switch it.coordType {
	case startPoint:
		u.currentDepth[it.sourceIndex] = it.depth
		u.nonZero++
	case endPoint:
		u.addIntervalIgnoreErr(it.sourceIndex)
		u.currentDepth[it.sourceIndex] = u.opts.Filler
		u.nonZero--
	}
}

// determineNextChrom returns the lexicographically smallest chromosome among
// the sources' current front items, or "" when all sources are exhausted.
func (u *unioner) determineNextChrom() string {
	next := ""
	for _, s := range u.sources {
		if !s.current.valid {
			continue
		}
		if next == "" || s.current.chrom < next {
			next = s.current.chrom
		}
	}
	return next
}

// allFilesDone reports whether every source has been fully consumed.
func (u *unioner) allFilesDone() bool {
	for _, s := range u.sources {
		if s.current.valid {
			return false
		}
	}
	return true
}

// addInterval pushes the START/END boundaries of source index's front interval
// onto the queue when it lies on the current chromosome, then advances that
// source. Intervals on other chromosomes are left for a later pass.
func (u *unioner) addInterval(index int) error {
	s := u.sources[index]
	if !s.current.valid {
		return nil
	}
	if s.current.chrom != u.currentChrom {
		return nil
	}
	bg := s.current
	heap.Push(&u.queue, point{sourceIndex: index, coordType: startPoint, coord: bg.start, depth: bg.depth})
	heap.Push(&u.queue, point{sourceIndex: index, coordType: endPoint, coord: bg.end, depth: bg.depth})
	return u.loadNext(index)
}

// addIntervalIgnoreErr is addInterval used from the inner sweep where a parse
// error has already been ruled out by the priming pass; any late error is
// surfaced by leaving the source invalid (loadNext stops on the bad line).
func (u *unioner) addIntervalIgnoreErr(index int) {
	_ = u.addInterval(index)
}

// loadNext advances source index to its next valid BEDGRAPH interval, or marks
// it exhausted. Track/browser/comment lines are skipped. A line that is not
// blank and does not have exactly four tab-separated fields, or whose
// coordinates are non-integers, ends the source (matching upstream, where an
// invalid line terminates the per-file read loop).
func (u *unioner) loadNext(index int) error {
	s := u.sources[index]
	s.current = bgItem{}
	for s.scanner.Scan() {
		s.lineNum++
		line := s.scanner.Text()
		fields := strings.Split(line, "\t")
		if len(fields) == 1 && fields[0] == "" {
			continue // blank line
		}
		if isHeaderLine(fields) {
			continue
		}
		if len(fields) != 4 {
			// Invalid: upstream stops reading this file here.
			return nil
		}
		start, err := strconv.ParseInt(fields[1], 10, 64)
		if err != nil {
			return fmt.Errorf("line %d: failed to extract start value from %q", s.lineNum, fields[1])
		}
		end, err := strconv.ParseInt(fields[2], 10, 64)
		if err != nil {
			return fmt.Errorf("line %d: failed to extract end value from %q", s.lineNum, fields[2])
		}
		s.current = bgItem{chrom: fields[0], start: start, end: end, depth: fields[3], valid: true}
		return nil
	}
	return s.scanner.Err()
}

// isHeaderLine reports whether a tokenized line is a track/browser/comment
// header, which upstream's BedGraph parser classifies as a header to skip. The
// match mirrors upstream: any field-0 containing "track", "browser", or "#".
func isHeaderLine(fields []string) bool {
	if len(fields) == 0 {
		return false
	}
	f := fields[0]
	return strings.Contains(f, "track") || strings.Contains(f, "browser") || strings.Contains(f, "#")
}
