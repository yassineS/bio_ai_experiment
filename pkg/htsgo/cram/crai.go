package cram

import (
	"bufio"
	"compress/gzip"
	"io"
	"os"
	"sort"
	"strconv"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/region"
)

// CRAIEntry is one line of a CRAM index (.crai) file. A .crai locates a
// single CRAM slice within the file: which reference span the slice
// covers and where to seek to read it. A region query returns the
// entries whose alignment span overlaps the wanted region.
type CRAIEntry struct {
	// RefID is the reference sequence id the slice's records align to:
	// >= 0 indexes the SAM @SQ lines, -1 marks an unmapped-reads slice
	// and -2 marks a slice spanning multiple references.
	RefID int32
	// AlignmentStart is the 1-based reference coordinate of the slice's
	// first aligned base. It is 0 for an unmapped-reads slice.
	AlignmentStart int32
	// AlignmentSpan is the number of reference bases the slice covers.
	AlignmentSpan int32
	// ContainerOffset is the absolute byte offset, from the start of the
	// CRAM file, of the container that holds the slice.
	ContainerOffset int64
	// SliceOffset is the byte offset of the slice within its container.
	SliceOffset int64
	// SliceSize is the size in bytes of the slice's blocks.
	SliceSize int64
}

// alignmentEnd returns the 1-based inclusive coordinate of the slice's
// last covered reference base. A zero-span slice covers AlignmentStart.
func (e CRAIEntry) alignmentEnd() int32 {
	if e.AlignmentSpan <= 0 {
		return e.AlignmentStart
	}
	return e.AlignmentStart + e.AlignmentSpan - 1
}

// CRAIIndex is the in-memory form of a parsed .crai index: the ordered
// list of slice entries the file contained.
type CRAIIndex struct {
	// Entries holds every index entry in the order the .crai listed them.
	Entries []CRAIEntry
}

// ReadCRAI parses a CRAM index (.crai) from r. A .crai is a
// gzip-compressed, tab-separated text file; each line carries six
// integers: reference id, alignment start, alignment span, the
// container's absolute byte offset, the slice's byte offset within the
// container and the slice's size in bytes. A malformed line is reported
// as an error; ReadCRAI never panics on malformed input.
func ReadCRAI(r io.Reader) (*CRAIIndex, error) {
	gz, err := gzip.NewReader(r)
	if err != nil {
		return nil, wrapf(err, "reading the .crai gzip stream")
	}
	defer gz.Close()
	idx := &CRAIIndex{}
	sc := bufio.NewScanner(gz)
	sc.Buffer(make([]byte, 0, 64*1024), 1<<20)
	line := 0
	for sc.Scan() {
		line++
		text := sc.Text()
		if text == "" {
			continue // tolerate a blank trailing line.
		}
		entry, err := parseCRAILine(text)
		if err != nil {
			return nil, errFormat(".crai line %d: %v", line, err)
		}
		idx.Entries = append(idx.Entries, entry)
	}
	if err := sc.Err(); err != nil {
		return nil, wrapf(err, "reading the .crai stream")
	}
	return idx, nil
}

// OpenCRAI opens and parses the named .crai index file.
func OpenCRAI(path string) (*CRAIIndex, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	return ReadCRAI(f)
}

// parseCRAILine parses one tab-separated .crai record into a CRAIEntry.
// It requires exactly six integer fields and rejects anything else,
// returning an error rather than panicking so the fuzz target's malformed
// inputs surface as clean failures.
func parseCRAILine(text string) (CRAIEntry, error) {
	fields := splitTabs(text)
	if len(fields) != 6 {
		return CRAIEntry{}, errFormat("want 6 tab-separated fields, got %d", len(fields))
	}
	var nums [6]int64
	for i, f := range fields {
		v, err := strconv.ParseInt(f, 10, 64)
		if err != nil {
			return CRAIEntry{}, errFormat("field %d %q is not an integer", i+1, f)
		}
		nums[i] = v
	}
	// The container/slice offsets and the size are file positions: a
	// negative value is never valid and would mislead a seeking caller.
	if nums[3] < 0 || nums[4] < 0 || nums[5] < 0 {
		return CRAIEntry{}, errFormat("negative container offset, slice offset or slice size")
	}
	// Reference id, alignment start and span must fit a 32-bit field —
	// they index the @SQ lines and SAM coordinates, which are 32-bit.
	if nums[0] < minInt32 || nums[0] > maxInt32 ||
		nums[1] < minInt32 || nums[1] > maxInt32 ||
		nums[2] < minInt32 || nums[2] > maxInt32 {
		return CRAIEntry{}, errFormat("reference id, start or span out of 32-bit range")
	}
	if nums[2] < 0 {
		return CRAIEntry{}, errFormat("negative alignment span %d", nums[2])
	}
	return CRAIEntry{
		RefID:           int32(nums[0]),
		AlignmentStart:  int32(nums[1]),
		AlignmentSpan:   int32(nums[2]),
		ContainerOffset: nums[3],
		SliceOffset:     nums[4],
		SliceSize:       nums[5],
	}, nil
}

// minInt32 and maxInt32 bound the .crai fields that index 32-bit SAM
// quantities (reference id, alignment start, alignment span).
const (
	minInt32 = -1 << 31
	maxInt32 = 1<<31 - 1
)

// splitTabs splits a line on tab characters; it is the .crai field
// splitter. The fields are subslices of s — no copy.
func splitTabs(s string) []string {
	var out []string
	start := 0
	for i := 0; i < len(s); i++ {
		if s[i] == '\t' {
			out = append(out, s[start:i])
			start = i + 1
		}
	}
	out = append(out, s[start:])
	return out
}

// Query returns the index entries whose alignment range overlaps the
// half-open 0-based interval [beg0, end0) on reference refID. The result
// is sorted by container offset then slice offset so a caller can seek
// through the file in order. An end0 <= 0 means "to the end of the
// reference".
func (idx *CRAIIndex) Query(refID int32, beg0, end0 int64) []CRAIEntry {
	if end0 <= 0 {
		end0 = 1 << 62 // effectively unbounded for SAM coordinates.
	}
	var hits []CRAIEntry
	for _, e := range idx.Entries {
		if e.RefID != refID {
			continue
		}
		// The .crai start is 1-based inclusive; convert the slice span to
		// a half-open 0-based interval [sliceBeg0, sliceEnd0) to compare
		// with the query's [beg0, end0).
		sliceBeg0 := int64(e.AlignmentStart) - 1
		if sliceBeg0 < 0 {
			sliceBeg0 = 0
		}
		// alignmentEnd is the 1-based inclusive last base; the half-open
		// 0-based exclusive end is that value (1-based inclusive == 0-based
		// exclusive). A zero-span slice still covers its single start base.
		sliceEnd0 := int64(e.alignmentEnd())
		if sliceEnd0 < sliceBeg0+1 {
			sliceEnd0 = sliceBeg0 + 1
		}
		if sliceBeg0 < end0 && sliceEnd0 > beg0 {
			hits = append(hits, e)
		}
	}
	sort.Slice(hits, func(a, b int) bool {
		if hits[a].ContainerOffset != hits[b].ContainerOffset {
			return hits[a].ContainerOffset < hits[b].ContainerOffset
		}
		return hits[a].SliceOffset < hits[b].SliceOffset
	})
	return hits
}

// QueryRegion is the region-package-typed form of Query: it returns the
// index entries overlapping a ResolvedRegion. It is the entry point a
// region-query CLI path uses after resolving a chrom name to a reference
// id against the CRAM's SAM header.
func (idx *CRAIIndex) QueryRegion(r region.ResolvedRegion) []CRAIEntry {
	return idx.Query(int32(r.RefID), int64(r.Beg0), int64(r.End0))
}
