package tabix

import (
	"bufio"
	"sort"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
)

// targetInterval is a single half-open [Beg, End) interval on a chromosome,
// expressed in 0-based coordinates. It mirrors the internal representation
// htslib's regidx uses, except that End is exclusive here for symmetry with
// the rest of this package.
type targetInterval struct {
	beg, end int // 0-based half-open
}

// Targets holds a set of intervals parsed from a `-T/--targets` file, grouped
// by chromosome and sorted by start coordinate. It is used as a strict
// post-filter: a record is kept only when it overlaps at least one interval.
//
// Coordinate interpretation matches htslib's regidx: a file whose name ends in
// `.bed`, `.bed.gz`, or `.bed.bgz` is read as BED (0-based, half-open); any
// other file is read as the default tab format (1-based, inclusive). A line
// holding only a chromosome name selects that entire chromosome.
type Targets struct {
	byChrom map[string][]targetInterval
}

// isBEDName reports whether name should be parsed using BED coordinate
// conventions, mirroring htslib's regidx_init filename sniffing.
func isBEDName(name string) bool {
	lower := strings.ToLower(name)
	return strings.HasSuffix(lower, ".bed") ||
		strings.HasSuffix(lower, ".bed.gz") ||
		strings.HasSuffix(lower, ".bed.bgz")
}

// LoadTargets reads an interval file at path into a Targets set. The file may
// be plain text or bgzip/gzip compressed. The coordinate convention (BED vs
// tab) is chosen from the filename exactly as htslib does.
func LoadTargets(path string) (*Targets, error) {
	r, err := iohelper.OpenReader(path)
	if err != nil {
		return nil, err
	}
	defer r.Close()

	bed := isBEDName(path)
	t := &Targets{byChrom: make(map[string][]targetInterval)}
	scan := bufio.NewScanner(r)
	scan.Buffer(make([]byte, 0, 1<<16), 1<<24)
	for scan.Scan() {
		line := strings.TrimSpace(scan.Text())
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Fields(line)
		chrom := fields[0]
		iv, ok := parseTargetInterval(fields, bed)
		if !ok {
			continue
		}
		t.byChrom[chrom] = append(t.byChrom[chrom], iv)
	}
	if err := scan.Err(); err != nil {
		return nil, err
	}
	for c := range t.byChrom {
		ivs := t.byChrom[c]
		sort.Slice(ivs, func(a, b int) bool { return ivs[a].beg < ivs[b].beg })
		t.byChrom[c] = ivs
	}
	return t, nil
}

// parseTargetInterval converts the whitespace-split fields of one interval
// line into a 0-based half-open interval. bed selects BED conventions; ok is
// false when the line cannot be parsed (and should be skipped).
func parseTargetInterval(fields []string, bed bool) (targetInterval, bool) {
	if len(fields) == 1 {
		// Chromosome-only line: select the whole chromosome.
		return targetInterval{beg: 0, end: maxTargetCoor}, true
	}
	begRaw, ok := parseInt(fields[1])
	if !ok {
		return targetInterval{}, false
	}
	if bed {
		// BED: col2 is 0-based start, col3 is exclusive end.
		beg := begRaw
		end := beg + 1
		if len(fields) >= 3 {
			if e, ok2 := parseInt(fields[2]); ok2 {
				end = e
			}
		}
		if end <= beg {
			end = beg + 1
		}
		return targetInterval{beg: beg, end: end}, true
	}
	// Tab (default): col2 is 1-based start, col3 is 1-based inclusive end.
	if begRaw <= 0 {
		return targetInterval{}, false
	}
	beg := begRaw - 1
	end := beg + 1 // single-position record when no usable end column
	if len(fields) >= 3 {
		if e, ok2 := parseInt(fields[2]); ok2 && e > 0 {
			end = e // 1-based inclusive end == 0-based exclusive end
		}
	}
	if end <= beg {
		end = beg + 1
	}
	return targetInterval{beg: beg, end: end}, true
}

// maxTargetCoor is the open upper bound used for chromosome-only target lines.
const maxTargetCoor = 1 << 62

// parseInt parses a non-negative decimal integer, returning ok=false on any
// parse error so the caller can skip malformed lines.
func parseInt(s string) (int, bool) {
	n := 0
	if s == "" {
		return 0, false
	}
	for i := 0; i < len(s); i++ {
		c := s[i]
		if c < '0' || c > '9' {
			return 0, false
		}
		n = n*10 + int(c-'0')
	}
	return n, true
}

// Overlaps reports whether the half-open record interval [beg, end) on chrom
// overlaps any target interval. This mirrors htslib's regidx_overlap test
// applied by `tabix -T` (overlap on the 0-based inclusive record span
// [beg, end-1]).
func (t *Targets) Overlaps(chrom string, beg, end int) bool {
	ivs := t.byChrom[chrom]
	if len(ivs) == 0 {
		return false
	}
	// Intervals are sorted by beg; binary-search for the first interval that
	// starts at or after end (those cannot overlap), then scan the earlier
	// candidates. Because intervals are not merged, a linear scan over the
	// candidate prefix is the simplest correct approach and is fine for
	// typical target-file sizes.
	lo := sort.Search(len(ivs), func(i int) bool { return ivs[i].beg >= end })
	for i := 0; i < lo; i++ {
		if ivs[i].end > beg {
			return true
		}
	}
	return false
}
