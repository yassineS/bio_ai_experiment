// Join output modes for bedintersect (-loj, -wo, -wao, -wa+-wb), plus the
// -split block-aware variants. These modes echo the original input columns of
// A and B verbatim, in the original B-file order, exactly like upstream
// `bedtools intersect`. They therefore operate on raw, line-preserving records
// rather than the typed bed.Record used by the default code path.
package bedintersect

import (
	"bufio"
	"fmt"
	"io"
	"strconv"
	"strings"
)

// block is a half-open sub-interval [start,end) of a record. For non-split
// records (or records without BED12 blocks) the single block equals the whole
// record span.
type block struct {
	start int
	end   int
}

// blocksOf returns the sub-intervals of a record. With split=false (or a record
// that is not valid BED12) the whole span is the single block. With split=true
// and a valid BED12 record, the blockStarts/blockSizes columns are expanded
// into absolute coordinates, mirroring BlockMgr::getBlocksFromBed12. A
// zero-length whole-span block is expanded to [p-1,p+1] so it can still
// intersect, matching upstream's adjustZeroLength (the split overlap counter
// reports the expanded width, with no zero-length correction).
func blocksOf(rec *inRecord, split bool) []block {
	wholeSpan := func() []block {
		s, e, _ := effectiveBounds(rec.start, rec.end)
		return []block{{s, e}}
	}
	if !split {
		return wholeSpan()
	}
	// Only BED12 and BAM records carry block columns; VCF/GFF are always a
	// single span (matching the raw path's inBlocks so -split agrees across the
	// join and default output modes for the same input).
	if rec.format != fmtBED && rec.format != fmtBAM {
		return wholeSpan()
	}
	blks, ok := bed12BlocksFromFields(rec.fields, rec.start)
	if !ok {
		return wholeSpan()
	}
	return blks
}

// splitCSV splits a comma-separated list, dropping a single trailing comma
// (BED12 block columns conventionally end with one).
func splitCSV(s string) []string {
	s = strings.TrimSuffix(s, ",")
	if s == "" {
		return nil
	}
	return strings.Split(s, ",")
}

// blockSum totals the lengths of a set of blocks.
func blockSum(blks []block) int {
	total := 0
	for _, b := range blks {
		total += b.end - b.start
	}
	return total
}

// joinHit records a B record that passed all overlap filters for a given A,
// together with the total overlapping bases used by -wo/-wao.
type joinHit struct {
	b            *inRecord
	overlapBases int
}

// findJoinHits returns the B records overlapping A (in B-file order) and the
// per-hit overlap-base counts, applying the -s strand filter and the -f/-F/-r
// fraction filters with split-aware block math when opts.Split is set. The
// candidate scan uses full-record spans (as upstream's bin index does); block
// refinement only affects hit determination and overlap counts.
func findJoinHits(a *inRecord, bRecords []*inRecord, opts IntersectOptions) []joinHit {
	if opts.Split {
		return findJoinHitsSplit(a, bRecords, opts)
	}
	var hits []joinHit
	aLen := a.end - a.start
	aStart, aEnd, aZero := effectiveBounds(a.start, a.end)
	for _, b := range bRecords {
		if a.chrom != b.chrom {
			continue
		}
		if opts.StrandSpec && !sameStrandMatch(a.strand, b.strand) {
			continue
		}
		bStart, bEnd, bZero := effectiveBounds(b.start, b.end)
		overlapStart := max(aStart, bStart)
		overlapEnd := min(aEnd, bEnd)
		overlapBases := overlapEnd - overlapStart
		if overlapBases <= 0 {
			continue
		}
		// MinOverlap is a bedintersect-specific knob (upstream intersect has no
		// -m); honour it here so the join modes match the default rawOverlaps
		// path, which also filters on the detected overlap length. This is
		// tested before the zero-length correction so a zero-length hit (which
		// reports 0 overlapping bases but is still a valid intersection under
		// the default -m 1) is not dropped.
		if overlapBases < opts.MinOverlap {
			continue
		}
		// Zero-length records were expanded by 1bp on each side for detection;
		// undo that when reporting the overlap base count (upstream
		// RecordOutputMgr::reportOverlapDetail does maxStart++ / minEnd--).
		if aZero || bZero {
			overlapBases = (overlapEnd - 1) - (overlapStart + 1)
			if overlapBases < 0 {
				overlapBases = 0
			}
		}
		// Per-record fraction tests over whole spans, mirroring
		// Record::sameChromIntersects (default !eitherFraction: both must hold).
		if opts.FractionA > 0 {
			if fraction(overlapBases, aLen) < opts.FractionA {
				continue
			}
		}
		if opts.FractionB > 0 || opts.Reciprocal {
			bLen := b.end - b.start
			fracB := fraction(overlapBases, bLen)
			if opts.FractionB > 0 && fracB < opts.FractionB {
				continue
			}
			if opts.Reciprocal && fracB < opts.FractionA {
				continue
			}
		}
		hits = append(hits, joinHit{b: b, overlapBases: overlapBases})
	}
	return hits
}

// findJoinHitsSplit implements the -split hit selection, mirroring
// BlockMgr::findBlockedOverlaps: a B record is a candidate if any of its blocks
// overlaps any A block; the per-hit overlap count is the total block overlap for
// that B; and the -f/-F/-r fraction tests are applied ONCE across all hits
// combined (non-redundant overlap over the A block-sum and the summed B block
// lengths). If a combined test fails, every hit for this A is dropped.
func findJoinHitsSplit(a *inRecord, bRecords []*inRecord, opts IntersectOptions) []joinHit {
	aBlocks := blocksOf(a, true)
	aBlockSum := blockSum(aBlocks)
	var hits []joinHit
	var allOverlaps []block
	hitBlockSum := 0
	for _, b := range bRecords {
		if a.chrom != b.chrom {
			continue
		}
		if opts.StrandSpec && !sameStrandMatch(a.strand, b.strand) {
			continue
		}
		// Hit determination is purely block-vs-block below (blocksOf already
		// applies the zero-length expansion), so no whole-span pre-filter here:
		// a raw-span check would wrongly reject zero-length records.
		bBlocks := blocksOf(b, true)
		overlapBases := 0
		hadOverlap := false
		for _, hb := range bBlocks {
			for _, kb := range aBlocks {
				s := max(kb.start, hb.start)
				e := min(kb.end, hb.end)
				if e > s {
					overlapBases += e - s
					allOverlaps = append(allOverlaps, block{s, e})
					hadOverlap = true
				}
			}
		}
		if !hadOverlap {
			continue
		}
		hitBlockSum += blockSum(bBlocks)
		hits = append(hits, joinHit{b: b, overlapBases: overlapBases})
	}

	if len(hits) > 0 && (opts.FractionA > 0 || opts.FractionB > 0 || opts.Reciprocal) {
		uniq := nonRedundantOverlap(allOverlaps)
		if opts.FractionA > 0 && fraction(uniq, aBlockSum) < opts.FractionA {
			return nil
		}
		if opts.FractionB > 0 && fraction(uniq, hitBlockSum) < opts.FractionB {
			return nil
		}
		if opts.Reciprocal && fraction(uniq, hitBlockSum) < opts.FractionA {
			return nil
		}
	}
	return hits
}

// effectiveBounds returns the [start,end) used for overlap detection, expanding
// a zero-length record [p,p] to [p-1,p+1] exactly as upstream's
// Record::adjustZeroLength does so it can still intersect. The third return
// value reports whether the record was zero-length.
func effectiveBounds(start, end int) (int, int, bool) {
	if start == end {
		return start - 1, end + 1, true
	}
	return start, end, false
}

// fraction returns part/whole as a float64, guarding division by zero.
func fraction(part, whole int) float64 {
	if whole <= 0 {
		return 0
	}
	return float64(part) / float64(whole)
}

// nonRedundantOverlap merges a set of overlap intervals and returns the total
// covered length, so overlapping block pairs are not double counted in the
// fraction tests (mirrors BlockMgr::getNonRedundantOverlap).
func nonRedundantOverlap(intervals []block) int {
	if len(intervals) == 0 {
		return 0
	}
	// Insertion sort by start (interval counts are tiny: blocks-of-A * blocks-of-B).
	for i := 1; i < len(intervals); i++ {
		for j := i; j > 0 && intervals[j].start < intervals[j-1].start; j-- {
			intervals[j], intervals[j-1] = intervals[j-1], intervals[j]
		}
	}
	total := 0
	curStart := intervals[0].start
	curEnd := intervals[0].end
	for _, iv := range intervals[1:] {
		if iv.start <= curEnd {
			if iv.end > curEnd {
				curEnd = iv.end
			}
		} else {
			total += curEnd - curStart
			curStart = iv.start
			curEnd = iv.end
		}
	}
	total += curEnd - curStart
	return total
}

// nullDBString returns the placeholder columns upstream emits for a database
// record on an A-with-no-hits line under -loj/-wao. The shape depends on the
// detected record type of the B file, exactly as RecordOutputMgr::null does.
func nullDBString(dbType dbRecordType, numFields int) string {
	switch dbType {
	case dbBed3:
		return ".\t-1\t-1"
	case dbBed4, dbBedGraph:
		return ".\t-1\t-1\t."
	case dbBed5:
		return ".\t-1\t-1\t.\t-1"
	case dbBed6:
		return ".\t-1\t-1\t.\t-1\t."
	case dbBed12:
		return ".\t-1\t-1\t.\t-1\t.\t.\t.\t.\t.\t.\t."
	case dbVCF:
		// VcfRecord::printNull: chrom ".", POS "-1", then "." for every other
		// column (the affected reference span is not echoed for a null VCF).
		var sb strings.Builder
		sb.WriteString(".\t-1")
		for i := 2; i < numFields; i++ {
			sb.WriteString("\t.")
		}
		return sb.String()
	case dbGFF:
		// GffRecord::printNull: seqid/source/type ".", start/end "-1", then "."
		// for every remaining column (score, strand, frame, attributes).
		var sb strings.Builder
		sb.WriteString(".\t.\t.\t-1\t-1")
		for i := 5; i < numFields; i++ {
			sb.WriteString("\t.")
		}
		return sb.String()
	default: // dbBedPlus / dbBed6Plus: "." "-1" "-1" then "." for each extra col.
		var sb strings.Builder
		sb.WriteString(".\t-1\t-1")
		for i := 3; i < numFields; i++ {
			sb.WriteString("\t.")
		}
		return sb.String()
	}
}

// dbRecordType enumerates the BED record classifications relevant to the null
// placeholder shape (a subset of upstream's FileRecordTypeChecker types).
type dbRecordType int

const (
	dbBed3 dbRecordType = iota
	dbBed4
	dbBedGraph
	dbBed5
	dbBed6
	dbBed12
	dbBed6Plus
	dbBedPlus
	dbVCF
	dbGFF
)

// classifyDB determines the database record type from the first data record of
// B, mirroring FileRecordTypeChecker::isBedFormat's column-count + content
// rules. It also returns the field count, which the BED+ null shape depends on.
func classifyDB(bRecords []*inRecord) (dbRecordType, int) {
	if len(bRecords) == 0 {
		return dbBed3, 3
	}
	f := bRecords[0].fields
	n := len(f)
	// VCF and GFF B files have format-specific null placeholder shapes
	// (mirroring VcfRecord::printNull / GffRecord::printNull).
	switch bRecords[0].format {
	case fmtVCF:
		return dbVCF, n
	case fmtGFF:
		return dbGFF, n
	}
	switch {
	case n == 3:
		return dbBed3, 3
	case n == 4:
		if isNumericField(f[3]) {
			return dbBedGraph, 4
		}
		return dbBed4, 4
	case n == 5 && isNumericField(f[4]):
		return dbBed5, 5
	case n == 6 && isStrandField(f[5]):
		return dbBed6, 6
	case n == 12 && isStrandField(f[5]) && isNumericField(f[6]) && isNumericField(f[7]) && isNumericField(f[9]):
		return dbBed12, 12
	default:
		if n >= 6 && isStrandField(f[5]) {
			return dbBed6Plus, n
		}
		return dbBedPlus, n
	}
}

// isStrandField reports whether s is a valid strand token (+, -, ., *).
func isStrandField(s string) bool {
	return s == "+" || s == "-" || s == "." || s == "*"
}

// isNumericField mirrors upstream ParseTools::isNumeric: the string may contain
// only digits, sign, decimal point and exponent characters, and must contain at
// least one digit.
func isNumericField(s string) bool {
	hasDigit := false
	for _, c := range s {
		switch {
		case c >= '0' && c <= '9':
			hasDigit = true
		case c == '-' || c == '+' || c == '.' || c == 'e' || c == 'E':
		default:
			return false
		}
	}
	return hasDigit
}

// intersectJoin implements the -loj / -wo / -wao / (-wa -wb) output modes,
// echoing A and B columns verbatim in B-file order. It returns the number of
// output lines written.
func intersectJoin(readerA, readerB io.Reader, w io.Writer, opts IntersectOptions) (int, error) {
	bRecords, err := readInRecords(readerB)
	if err != nil {
		return 0, fmt.Errorf("error reading B intervals: %w", err)
	}
	dbType, dbFields := classifyDB(bRecords)
	nullB := nullDBString(dbType, dbFields)

	// Index B by chromosome, preserving file order within each chromosome.
	byChrom := make(map[string][]*inRecord)
	for _, b := range bRecords {
		byChrom[b.chrom] = append(byChrom[b.chrom], b)
	}

	out := bufio.NewWriter(w)
	bw := &bufWriter{w: out}
	count := 0

	aRecords, err := readInRecords(readerA)
	if err != nil {
		return 0, fmt.Errorf("error reading A intervals: %w", err)
	}
	for _, a := range aRecords {
		hits := findJoinHits(a, byChrom[a.chrom], opts)

		if len(hits) == 0 {
			// -wao and -loj still emit the A record (with a null B) when there
			// are no hits; -wo and (-wa -wb) emit nothing.
			switch {
			case opts.WriteAllOverlap:
				bw.writeString(a.line)
				bw.writeString("\t")
				bw.writeString(nullB)
				bw.writeString("\t0\n")
				count++
			case opts.LeftJoin:
				bw.writeString(a.line)
				bw.writeString("\t")
				bw.writeString(nullB)
				bw.writeString("\n")
				count++
			}
			continue
		}

		for _, h := range hits {
			bw.writeString(a.line)
			bw.writeString("\t")
			bw.writeString(h.b.line)
			if opts.WriteOverlap || opts.WriteAllOverlap {
				bw.writeString("\t")
				bw.writeString(strconv.Itoa(h.overlapBases))
			}
			bw.writeString("\n")
			count++
		}
	}
	if bw.err != nil {
		return 0, bw.err
	}
	if err := out.Flush(); err != nil {
		return 0, fmt.Errorf("error flushing output: %w", err)
	}
	return count, nil
}

// bufWriter accumulates the first write error so the hot loop stays branch-light.
type bufWriter struct {
	w   *bufio.Writer
	err error
}

func (b *bufWriter) writeString(s string) {
	if b.err != nil {
		return
	}
	_, b.err = b.w.WriteString(s)
}

// usesJoinMode reports whether any option requires the raw, B-ordered join
// output path (rather than the default typed path).
func (opts IntersectOptions) usesJoinMode() bool {
	return opts.LeftJoin || opts.WriteOverlap || opts.WriteAllOverlap || (opts.WriteA && opts.WriteB)
}
