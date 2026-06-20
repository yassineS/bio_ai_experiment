package giab

import (
	"bufio"
	"compress/gzip"
	"io"
	"os"
	"sort"
	"strconv"
	"strings"
)

// interval is a half-open [start, end) BED interval (0-based start).
type interval struct{ start, end int }

// RegionSet is a set of genomic intervals keyed by chromosome, used to restrict
// the concordance comparison to the high-confidence BED. Intervals are stored
// sorted and merged per chromosome so membership tests are a binary search.
type RegionSet struct {
	byChrom map[string][]interval
	n       int
}

// Empty reports whether the region set contains no intervals.
func (rs *RegionSet) Empty() bool { return rs == nil || rs.n == 0 }

// Count returns the number of intervals in the set.
func (rs *RegionSet) Count() int {
	if rs == nil {
		return 0
	}
	return rs.n
}

// Contains reports whether a 1-based VCF position on chrom falls inside any
// interval. VCF POS is 1-based; BED is 0-based half-open, so a VCF position p
// is inside [start, end) when start <= p-1 < end, i.e. start < p <= end.
func (rs *RegionSet) Contains(chrom string, pos int) bool {
	if rs.Empty() {
		return false
	}
	ivs := rs.byChrom[chrom]
	if len(ivs) == 0 {
		return false
	}
	p0 := pos - 1 // to 0-based
	// Find the first interval whose end > p0.
	i := sort.Search(len(ivs), func(i int) bool { return ivs[i].end > p0 })
	return i < len(ivs) && ivs[i].start <= p0
}

// ParseBED reads a BED stream (at least three columns: chrom, start, end) into
// a RegionSet. Track/browser/comment lines and blank lines are skipped. Extra
// columns are ignored. Intervals are merged per chromosome.
func ParseBED(r io.Reader) (*RegionSet, error) {
	tmp := map[string][]interval{}
	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 0, 1<<20), 1<<26)
	for sc.Scan() {
		line := strings.TrimRight(sc.Text(), "\r")
		if line == "" || strings.HasPrefix(line, "#") ||
			strings.HasPrefix(line, "track") || strings.HasPrefix(line, "browser") {
			continue
		}
		f := strings.Split(line, "\t")
		if len(f) < 3 {
			f = strings.Fields(line)
		}
		if len(f) < 3 {
			continue
		}
		start, err := strconv.Atoi(strings.TrimSpace(f[1]))
		if err != nil {
			continue
		}
		end, err := strconv.Atoi(strings.TrimSpace(f[2]))
		if err != nil {
			continue
		}
		if end < start {
			start, end = end, start
		}
		tmp[f[0]] = append(tmp[f[0]], interval{start, end})
	}
	if err := sc.Err(); err != nil {
		return nil, err
	}
	rs := &RegionSet{byChrom: map[string][]interval{}}
	for chrom, ivs := range tmp {
		merged := mergeIntervals(ivs)
		rs.byChrom[chrom] = merged
		rs.n += len(merged)
	}
	return rs, nil
}

// ParseBEDFile parses a BED file by path, transparently decompressing a
// gzip/bgzip-framed .gz (GIAB stratification BEDs are commonly bgzipped, which
// is a valid gzip stream the stdlib reader decodes).
func ParseBEDFile(path string) (*RegionSet, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer f.Close()
	var r io.Reader = f
	if strings.HasSuffix(path, ".gz") {
		gz, err := gzip.NewReader(f)
		if err != nil {
			return nil, err
		}
		gz.Multistream(true) // bgzip is a concatenation of gzip members
		defer gz.Close()
		r = gz
	}
	return ParseBED(r)
}

// mergeIntervals sorts and coalesces overlapping/adjacent intervals.
func mergeIntervals(ivs []interval) []interval {
	if len(ivs) == 0 {
		return nil
	}
	sort.Slice(ivs, func(i, j int) bool {
		if ivs[i].start != ivs[j].start {
			return ivs[i].start < ivs[j].start
		}
		return ivs[i].end < ivs[j].end
	})
	out := ivs[:1]
	for _, iv := range ivs[1:] {
		last := &out[len(out)-1]
		if iv.start <= last.end {
			if iv.end > last.end {
				last.end = iv.end
			}
			continue
		}
		out = append(out, iv)
	}
	return out
}
