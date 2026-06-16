// Shared block, fraction and null-placeholder helpers for the upstream-parity
// intersect output modes. The -split block math, the B-file record-type
// classification used to shape -loj/-wao null placeholders, and the
// usesJoinMode dispatch all live here; the per-A emission lives in
// rawintersect.go (the emitter type).
package bedintersect

import (
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
	// single span.
	if rec.format != fmtBED && rec.format != fmtBAM {
		return wholeSpan()
	}
	blks, ok := bed12BlocksFromFields(rec.fields, rec.start)
	if !ok {
		return wholeSpan()
	}
	return blks
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

// usesJoinMode reports whether any option requires the join output path (which
// echoes both A and B columns), namely -loj, -wo, -wao, or -wa together with
// -wb.
func (opts IntersectOptions) usesJoinMode() bool {
	return opts.LeftJoin || opts.WriteOverlap || opts.WriteAllOverlap || (opts.WriteA && opts.WriteB)
}
