// Package bedclosest finds, for each interval in A, the closest interval(s) in
// B (mirrors `bedtools closest`).
//
// Both inputs MUST be sorted on (chrom, start). For each A interval the closest
// B on the same chromosome is reported. The output line is A's columns + the
// chosen B's columns, optionally followed by a distance column (only when `-d`
// or `-D` is requested). The selection, ordering, tie handling, k-closest, and
// directional behaviours mirror upstream bedtools' CloseSweep precisely, and
// the no-hit placeholder ("null") matches the per-record-type shape that
// upstream's RecordOutputMgr::printNull emits.
package bedclosest

import (
	"bufio"
	"fmt"
	"io"
	"sort"
	"strconv"
	"strings"
)

// DistanceMode selects how the distance between A and B is computed and
// reported. It only takes effect when ReportDistance is set.
type DistanceMode int

const (
	// DistanceRef reports an unsigned absolute distance (upstream `-d`). The
	// sign is never applied in this mode.
	DistanceRef DistanceMode = iota
	// DistanceSignedRef computes the sign on the reference: downstream of A
	// (to the right) is positive, upstream (to the left) is negative
	// (upstream `-D ref`).
	DistanceSignedRef
	// DistanceA computes the sign relative to A's strand (BED6 col 6). On a
	// '-' strand A, what was "downstream on the reference" becomes "upstream"
	// and the sign is flipped (upstream `-D a`).
	DistanceA
	// DistanceB computes the sign relative to B's strand (upstream `-D b`).
	DistanceB
)

// signDistance reports whether the active distance mode applies a sign to the
// distance (i.e. one of the `-D` modes).
func (d DistanceMode) signDistance() bool {
	return d == DistanceSignedRef || d == DistanceA || d == DistanceB
}

// TieBreak controls how ties among multiple equally-close B intervals are
// handled.
type TieBreak int

const (
	// TieAll emits one output row per tied B.
	TieAll TieBreak = iota
	// TieFirst emits only the first tied B (lowest position in B's input).
	TieFirst
	// TieLast emits only the last tied B.
	TieLast
)

// MultiDBMode selects how multiple -b databases are combined.
type MultiDBMode int

const (
	// MultiDBEach reports the closest B from each database independently
	// (upstream's default `-mdb each`).
	MultiDBEach MultiDBMode = iota
	// MultiDBAll treats every database as one combined set and reports the
	// single overall closest (upstream `-mdb all`).
	MultiDBAll
)

// Options configures Closest.
type Options struct {
	// ReportDistance controls whether a trailing distance column is emitted.
	// It is set when the user passes `-d` or `-D`. Upstream omits the column
	// otherwise.
	ReportDistance bool
	// DistanceMode controls the distance column's sign convention. Ignored
	// unless ReportDistance is set.
	DistanceMode DistanceMode

	// KClosest is the number of closest hits to report per A interval
	// (upstream `-k`). It defaults to 1 when zero.
	KClosest int

	// TieBreak controls how ties (multiple equally-close B intervals) are
	// resolved. Default TieAll.
	TieBreak TieBreak

	// IgnoreUpstream (`-iu`) drops B features upstream of A; IgnoreDownstream
	// (`-id`) is the symmetric flag; IgnoreOverlaps (`-io`) drops overlapping
	// B features. Upstream/downstream are determined relative to the active
	// `-D` mode.
	IgnoreUpstream   bool
	IgnoreDownstream bool
	IgnoreOverlaps   bool

	// ForceUpstream (`-fu`) reports the closest upstream B before any
	// equally-or-closer downstream/overlapping B; ForceDownstream (`-fd`) is
	// the symmetric flag.
	ForceUpstream   bool
	ForceDownstream bool

	// SameStrand (`-s`) restricts candidate B intervals to those on the same
	// strand as A; OppositeStrand (`-S`) restricts to the opposite strand.
	SameStrand     bool
	OppositeStrand bool

	// DifferentNames (`-N`) requires the reported closest B to have a
	// different name (BED column 4) than A.
	DifferentNames bool

	// MultiDBMode controls how multiple -b databases are combined.
	MultiDBMode MultiDBMode

	// DBLabels, when non-nil, supplies the label printed in the database
	// column for each database, in -b order (from `-names`/`-filenames`).
	// When nil, ClosestMulti prints the 1-based database index instead.
	DBLabels []string

	// PrintHeader (`-header`) echoes A's leading header/comment lines to the
	// output before the result rows.
	PrintHeader bool

	// WarnWriter, when non-nil, enables upstream's cross-file chromosome
	// sort-order and naming-convention validation (the closest sortAndNaming
	// checks). WARNING and ERROR text is written here (upstream's stderr). When
	// nil the validation is skipped entirely and only the legacy per-file
	// (chrom,start) sort check applies. The CLI always sets this to os.Stderr.
	WarnWriter io.Writer

	// QueryName and DBNames provide the file-name labels printed in the
	// validation WARNING/ERROR messages (the query file and each -b database,
	// in -b order). They are only consulted when WarnWriter is non-nil. When a
	// name is empty a generic placeholder is used.
	QueryName string
	DBNames   []string
}

// Validate reports configuration errors that cannot be captured by the type
// system, mirroring upstream's mutually-exclusive flag checks.
func (o Options) Validate() error {
	if o.SameStrand && o.OppositeStrand {
		return fmt.Errorf("-s and -S are mutually exclusive")
	}
	return nil
}

// kVal returns the effective number of closest hits to report.
func (o Options) kVal() int {
	if o.KClosest <= 0 {
		return 1
	}
	return o.KClosest
}

// Strand value classification mirroring upstream Record::strandType.
const (
	strandForward = iota
	strandReverse
	strandUnknown
)

// Row is the parsed representation of a BED line. The original fields are
// preserved so the full input column count round-trips to the output.
type Row struct {
	Fields  []string
	Chrom   string
	Start   int
	End     int
	StrandV int // strandForward / strandReverse / strandUnknown
}

// nameOf returns a row's BED name (column 4) or "" when absent.
func nameOf(r *Row) string {
	if len(r.Fields) >= 4 {
		return r.Fields[3]
	}
	return ""
}

// strandValOf maps a BED strand field to upstream's strand classification.
// Records with fewer than 6 columns have no strand and are treated as forward,
// matching upstream's default for unstranded intervals.
func strandValOf(fields []string) int {
	if len(fields) >= 6 {
		switch fields[5] {
		case "+":
			return strandForward
		case "-":
			return strandReverse
		default:
			return strandUnknown
		}
	}
	return strandForward
}

// ReadAll parses every BED record from r into Row values. Lines that begin
// with '#', "track", or "browser" are skipped (collected separately as header
// lines), as are blank lines. The returned header slice holds the leading
// header lines (before the first record) for `-header` support.
func ReadAll(r io.Reader) (rows []*Row, header []string, err error) {
	sc := bufio.NewScanner(r)
	sc.Buffer(make([]byte, 64*1024), 64*1024*1024)
	lineNum := 0
	seenRecord := false
	for sc.Scan() {
		lineNum++
		raw := sc.Text()
		trimmed := strings.TrimSpace(raw)
		if trimmed == "" {
			continue
		}
		if strings.HasPrefix(trimmed, "#") ||
			strings.HasPrefix(trimmed, "track") || strings.HasPrefix(trimmed, "browser") {
			if !seenRecord {
				header = append(header, raw)
			}
			continue
		}
		fields := strings.Split(raw, "\t")
		if len(fields) < 3 {
			return nil, nil, fmt.Errorf("line %d: BED record must have at least 3 fields", lineNum)
		}
		start, e := strconv.Atoi(strings.TrimSpace(fields[1]))
		if e != nil {
			return nil, nil, fmt.Errorf("line %d: invalid start %q: %v", lineNum, fields[1], e)
		}
		end, e := strconv.Atoi(strings.TrimSpace(fields[2]))
		if e != nil {
			return nil, nil, fmt.Errorf("line %d: invalid end %q: %v", lineNum, fields[2], e)
		}
		if end < start {
			return nil, nil, fmt.Errorf("line %d: end < start (%d < %d)", lineNum, end, start)
		}
		rows = append(rows, &Row{Fields: fields, Chrom: fields[0], Start: start, End: end, StrandV: strandValOf(fields)})
		seenRecord = true
	}
	if e := sc.Err(); e != nil {
		return nil, nil, e
	}
	return rows, header, nil
}

// CheckSorted returns an error naming the offending line if rows are not sorted
// on (chrom, start). Equal chromosomes must have non-decreasing starts.
func CheckSorted(rows []*Row, label string) error {
	for i := 1; i < len(rows); i++ {
		prev, cur := rows[i-1], rows[i]
		if cur.Chrom == prev.Chrom {
			if cur.Start < prev.Start {
				return fmt.Errorf("%s is not sorted on (chrom, start) at record %d: %s\t%d came after %s\t%d",
					label, i+1, cur.Chrom, cur.Start, prev.Chrom, prev.Start)
			}
		}
	}
	return nil
}

// RecordType classifies a database's BED record type, mirroring upstream's
// FileRecordTypeChecker. It determines the null placeholder shape.
type RecordType int

const (
	recBed3 RecordType = iota
	recBed4
	recBed5
	recBed6
	recBed12
	recBedGraph
	recBedPlus
)

// isNumericField mirrors upstream isNumeric: digits with optional sign,
// decimal point, and exponent markers; at least one digit required.
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

// isStrandField mirrors upstream isStrandField.
func isStrandField(s string) bool {
	return s == "+" || s == "-" || s == "." || s == "*"
}

// detectRecordType determines the BED record type from the first record,
// mirroring FileRecordTypeChecker::handleBed. nFields is the column count.
func detectRecordType(first *Row) (RecordType, int) {
	if first == nil {
		return recBed3, 3
	}
	f := first.Fields
	n := len(f)
	switch {
	case n == 3:
		return recBed3, 3
	case n == 4:
		if isNumericField(f[3]) {
			return recBedGraph, 4
		}
		return recBed4, 4
	case n == 5 && isNumericField(f[4]):
		return recBed5, 5
	case n == 6 && isStrandField(f[5]):
		return recBed6, 6
	case n == 12 && isStrandField(f[5]) && isNumericField(f[6]) && isNumericField(f[7]) && isNumericField(f[9]):
		return recBed12, 12
	case n > 3:
		// BED6_PLUS / BED_PLUS: Bed3 null + dots for every extra field.
		return recBedPlus, n
	default:
		return recBed3, 3
	}
}

// nullFields returns the per-record-type null placeholder fields, matching
// upstream's printNull for each BED record type.
func nullFields(rt RecordType, n int) []string {
	switch rt {
	case recBed3:
		return []string{".", "-1", "-1"}
	case recBed4:
		return []string{".", "-1", "-1", "."}
	case recBed5:
		return []string{".", "-1", "-1", ".", "-1"}
	case recBed6:
		return []string{".", "-1", "-1", ".", "-1", "."}
	case recBed12:
		return []string{".", "-1", "-1", ".", "-1", ".", ".", ".", ".", ".", ".", "."}
	case recBedGraph:
		return []string{".", "-1", "-1", "."}
	default: // recBedPlus: ". -1 -1" then a dot per extra column.
		out := []string{".", "-1", "-1"}
		for i := 3; i < n; i++ {
			out = append(out, ".")
		}
		return out
	}
}

// streamDir classifies a non-overlapping B hit relative to A.
type streamDir int

const (
	streamOverlap streamDir = iota
	streamUpstream
	streamDownstream
)

// chromB holds one database's rows for a single chromosome, in input order.
type chromB struct {
	rows []*Row
}

// db is a single database: rows bucketed by chromosome plus its null shape.
type db struct {
	byChrom  map[string]*chromB
	nullFlds []string
}

// indexDB buckets a database's rows by chromosome (preserving input order) and
// records the null placeholder shape from the first record.
func indexDB(bRows []*Row) *db {
	byChrom := make(map[string]*chromB, 16)
	for _, b := range bRows {
		c := byChrom[b.Chrom]
		if c == nil {
			c = &chromB{}
			byChrom[b.Chrom] = c
		}
		c.rows = append(c.rows, b)
	}
	var first *Row
	if len(bRows) > 0 {
		first = bRows[0]
	}
	rt, n := detectRecordType(first)
	return &db{byChrom: byChrom, nullFlds: nullFields(rt, n)}
}

// hit holds a chosen B together with its (possibly signed) reported distance.
type hit struct {
	b    *Row
	dist int64
}

// classify returns the stream direction of B relative to A and the absolute
// distance, mirroring CloseSweep::considerRecord. Overlap returns dist 0.
func classify(a, b *Row, opts Options) (streamDir, int64) {
	if intersects(a, b) {
		return streamOverlap, 0
	}
	if b.Start >= a.End {
		// HIT IS TO THE RIGHT OF THE QUERY.
		dist := int64(b.Start-a.End) + 1
		if opts.DistanceMode.signDistance() {
			if (opts.DistanceMode == DistanceA && a.StrandV == strandReverse) ||
				(opts.DistanceMode == DistanceB && b.StrandV == strandForward) {
				return streamUpstream, dist
			}
		}
		return streamDownstream, dist
	}
	// HIT IS TO THE LEFT OF THE QUERY.
	dist := int64(a.Start-b.End) + 1
	if opts.DistanceMode.signDistance() {
		if (opts.DistanceMode == DistanceA && a.StrandV == strandReverse) ||
			(opts.DistanceMode == DistanceB && b.StrandV == strandForward) {
			return streamDownstream, dist
		}
	}
	return streamUpstream, dist
}

// intersects mirrors Record::sameChromIntersects for the no-fraction case,
// including the zero-length and touching-interval edge cases.
func intersects(a, b *Row) bool {
	maxStart := a.Start
	if b.Start > maxStart {
		maxStart = b.Start
	}
	minEnd := a.End
	if b.End < minEnd {
		minEnd = b.End
	}
	if minEnd < maxStart {
		return false
	}
	otherZeroLen := b.Start == b.End
	localZeroLen := a.Start == a.End
	if minEnd == maxStart && !otherZeroLen && !localZeroLen {
		return false
	}
	return true
}

// strandBad reports whether B should be excluded for A under -s/-S, mirroring
// the badStrand test in CloseSweep::tryToAddRecord.
func strandBad(a, b *Row, opts Options) bool {
	if !opts.SameStrand && !opts.OppositeStrand {
		return false
	}
	hasUnknown := a.StrandV == strandUnknown || b.StrandV == strandUnknown
	if opts.SameStrand {
		return hasUnknown || a.StrandV != b.StrandV
	}
	// OppositeStrand
	return hasUnknown || a.StrandV == b.StrandV
}

// distGroup holds all candidate B's tied at one absolute distance, in B input
// order (push_back order in upstream's RecDistList).
type distGroup struct {
	dist  int64
	elems []*Row
}

// recDistList accumulates candidate B's bucketed and sorted by absolute
// distance, capped at the k closest distinct distances, mirroring upstream's
// RecDistList.
type recDistList struct {
	k      int
	groups []distGroup // sorted ascending by dist
}

func newRecDistList(k int) *recDistList { return &recDistList{k: k} }

// addRec inserts a B at the given absolute distance, keeping at most k distinct
// distances, mirroring RecDistList::addRec. A new distance larger than all kept
// is dropped once k distinct distances are already held.
func (l *recDistList) addRec(dist int64, rec *Row) {
	// Find existing group or insertion point (groups sorted ascending).
	idx := sort.Search(len(l.groups), func(i int) bool { return l.groups[i].dist >= dist })
	if idx < len(l.groups) && l.groups[idx].dist == dist {
		l.groups[idx].elems = append(l.groups[idx].elems, rec)
		return
	}
	if len(l.groups) >= l.k {
		// Already have k distinct distances. Only insert if smaller than max.
		if idx >= len(l.groups) {
			return
		}
		// Insert and drop the largest.
		l.groups = append(l.groups, distGroup{})
		copy(l.groups[idx+1:], l.groups[idx:])
		l.groups[idx] = distGroup{dist: dist, elems: []*Row{rec}}
		l.groups = l.groups[:l.k]
		return
	}
	// Room to add a new distinct distance.
	l.groups = append(l.groups, distGroup{})
	copy(l.groups[idx+1:], l.groups[idx:])
	l.groups[idx] = distGroup{dist: dist, elems: []*Row{rec}}
}

// hitsForA returns the closest hits for A in one database, applying the full
// selection/ordering/k/tie/directional logic of CloseSweep.
func hitsForA(a *Row, d *db, opts Options) []hit {
	k := opts.kVal()
	up := newRecDistList(k)
	down := newRecDistList(k)
	over := newRecDistList(k)

	if c := d.byChrom[a.Chrom]; c != nil {
		for _, b := range c.rows {
			stream, dist := classify(a, b, opts)
			// Exclusion (badStrand / badNames / badStream).
			bad := strandBad(a, b, opts) ||
				(opts.DifferentNames && nameOf(a) == nameOf(b))
			switch stream {
			case streamOverlap:
				if opts.IgnoreOverlaps || bad {
					continue
				}
				over.addRec(0, b)
			case streamUpstream:
				if opts.IgnoreUpstream || bad {
					continue
				}
				up.addRec(dist, b)
			case streamDownstream:
				if opts.IgnoreDownstream || bad {
					continue
				}
				down.addRec(dist, b)
			}
		}
	}

	return finalize(up, down, over, opts)
}

// finalize merges the upstream/downstream/overlap candidate lists into the
// ordered output hit list, mirroring CloseSweep::finalizeSelections.
func finalize(up, down, over *recDistList, opts Options) []hit {
	k := opts.kVal()
	var out []hit
	used := 0

	emitGroup := func(g distGroup, signed int64) {
		elems := g.elems
		switch opts.TieBreak {
		case TieFirst:
			out = append(out, hit{b: elems[0], dist: signed})
			used++
		case TieLast:
			out = append(out, hit{b: elems[len(elems)-1], dist: signed})
			used++
		default:
			for _, e := range elems {
				out = append(out, hit{b: e, dist: signed})
				used++
			}
		}
	}

	upI, downI := 0, 0

	if opts.ForceUpstream {
		for upI < len(up.groups) && used < k {
			emitGroup(up.groups[upI], 0-up.groups[upI].dist)
			upI++
		}
	}
	if opts.ForceDownstream {
		for downI < len(down.groups) && used < k {
			emitGroup(down.groups[downI], down.groups[downI].dist)
			downI++
		}
	}

	// Overlaps (distance 0) next.
	if used < k && len(over.groups) > 0 {
		emitGroup(over.groups[0], 0)
	}

	// Merge upstream/downstream by increasing distance until k reached.
	const maxDist = int64(1) << 62
	for used < k {
		upDist := maxDist
		downDist := maxDist
		if upI < len(up.groups) {
			upDist = up.groups[upI].dist
		}
		if downI < len(down.groups) {
			downDist = down.groups[downI].dist
		}
		if upDist == maxDist && downDist == maxDist {
			break
		}
		tie := upDist == downDist
		usedUp, usedDown := false, false
		if upDist < downDist || (tie && opts.TieBreak != TieLast) {
			emitGroup(up.groups[upI], 0-upDist)
			upI++
			usedUp = true
		}
		if downDist < upDist || (tie && opts.TieBreak != TieFirst) {
			emitGroup(down.groups[downI], downDist)
			downI++
			usedDown = true
		}
		if tie {
			if usedUp && !usedDown {
				downI++
			} else if usedDown && !usedUp {
				upI++
			}
		}
	}

	return out
}

// Closest reads BED records from readerA and readerB, finds the closest B for
// each A, and writes the result to writer. Both inputs MUST be sorted on
// (chrom, start). Returns the number of output rows. Closest is the
// single-database entry point with no database-label column.
func Closest(readerA, readerB io.Reader, writer io.Writer, opts Options) (int, error) {
	return ClosestMulti(readerA, []io.Reader{readerB}, writer, opts)
}

// ClosestMulti reads BED records from readerA and one or more B databases,
// finds the closest B(s) for each A, and writes the result to writer. Every
// input MUST be sorted on (chrom, start).
//
// With a single database and no labels, the output omits the database column
// (matching upstream). With multiple databases or with -names/-filenames, a
// database-label column is inserted between A's and B's columns.
func ClosestMulti(readerA io.Reader, readersB []io.Reader, writer io.Writer, opts Options) (int, error) {
	if err := opts.Validate(); err != nil {
		return 0, err
	}
	if len(readersB) == 0 {
		return 0, fmt.Errorf("at least one -b database is required")
	}
	if opts.DBLabels != nil && len(opts.DBLabels) != len(readersB) {
		return 0, fmt.Errorf("number of labels (%d) must match the number of -b files (%d)", len(opts.DBLabels), len(readersB))
	}

	aRows, header, err := ReadAll(readerA)
	if err != nil {
		return 0, fmt.Errorf("error reading A: %w", err)
	}

	dbs := make([]*db, len(readersB))
	bRowsList := make([][]*Row, len(readersB))
	for i, rb := range readersB {
		bRows, _, err := ReadAll(rb)
		if err != nil {
			return 0, fmt.Errorf("error reading B[%d]: %w", i, err)
		}
		bRowsList[i] = bRows
		dbs[i] = indexDB(bRows)
	}

	// Validation strategy. With WarnWriter set (the CLI path) we run upstream's
	// cross-file chromosome sort-order / naming-convention checks, which decide
	// how many query records get emitted and whether the run aborts (matching
	// upstream's WARNING/ERROR text and exit code). Without it we fall back to
	// the legacy per-file (chrom,start) sort check.
	emitLimit := len(aRows)
	var validateErr error
	// dbReachable gates a database's hits to the chromosomes the upstream sweep
	// could actually reach. When validation is off it is nil (no gating).
	var dbReachable func(dbIdx int, chrom string) bool
	if opts.WarnWriter != nil {
		v := newSweepValidator(aRows, bRowsList, validationFileNames(opts, len(readersB)), opts.WarnWriter)
		validateErr = v.runSweep()
		if validateErr != nil {
			emitLimit = v.emitted
		}
		dbReachable = v.dbReachable
	} else {
		if err := CheckSorted(aRows, "A"); err != nil {
			return 0, err
		}
		for i, bRows := range bRowsList {
			if err := CheckSorted(bRows, fmt.Sprintf("B[%d]", i)); err != nil {
				return 0, err
			}
		}
	}

	// Whether to print the database-label column. Upstream prints it whenever
	// there is more than one database or a -names/-filenames label is set.
	printDBCol := len(readersB) > 1 || opts.DBLabels != nil

	label := func(i int) string {
		if opts.DBLabels != nil {
			return opts.DBLabels[i]
		}
		return strconv.Itoa(i + 1)
	}

	// out collects output through a flush-boundary model mirroring upstream's
	// RecordOutputMgr 16K buffer: a mid-stream ERROR loses the unflushed
	// remainder, so only the committed prefix reaches stdout, whereas a clean
	// run or a destructor-time ERROR flushes everything.
	out := newFlushBuffer()

	if opts.PrintHeader {
		for _, h := range header {
			out.writeString(h)
			out.writeByte('\n')
		}
	}

	// signedForOutput renders the reported distance: signed for `-D` modes,
	// absolute for `-d`, mirroring RecordOutputMgr's
	// `dist = signDistance() ? dist : abs(dist)`.
	signedForOutput := func(dist int64) int64 {
		if opts.DistanceMode.signDistance() {
			return dist
		}
		if dist < 0 {
			return -dist
		}
		return dist
	}

	count := 0
	emit := func(a *Row, lbl string, bFields []string, dist int64, hasDist bool) error {
		writeRowBuf(out, a, printDBCol, lbl, bFields, dist, hasDist)
		count++
		return nil
	}

	// emitNull writes the single no-hit placeholder row for A. With multiple
	// databases (or labels) the database column is a literal ".", and the null
	// shape comes from database 0, matching upstream's null(false, true).
	emitNull := func(a *Row) error {
		return emit(a, ".", dbs[0].nullFlds, -1, opts.ReportDistance)
	}

	// emitLimit caps how many query records are written. It is len(aRows) for a
	// clean run or a destructor-time abort (all output produced, then exit 1),
	// and the pre-abort processed count for a mid-stream abort.
	// reachable reports whether database dbIdx may serve hits for A's chrom. It
	// is true when validation is off (no gating), and otherwise reflects the
	// sweep's chromosome reachability so a database stuck on a later chromosome
	// yields null exactly as upstream does.
	reachable := func(dbIdx int, chrom string) bool {
		return dbReachable == nil || dbReachable(dbIdx, chrom)
	}

	for _, a := range aRows[:emitLimit] {
		if opts.MultiDBMode == MultiDBAll && len(dbs) > 1 {
			if err := emitAllMode(a, dbs, opts, reachable, signedForOutput, emit, emitNull); err != nil {
				return count, err
			}
			continue
		}
		// MultiDBEach (or single db): accumulate hits across all databases.
		// Only when no database yields any hit do we emit a single null row.
		any := false
		for dbIdx, d := range dbs {
			if !reachable(dbIdx, a.Chrom) {
				continue
			}
			for _, h := range hitsForA(a, d, opts) {
				any = true
				if err := emit(a, label(dbIdx), h.b.Fields, signedForOutput(h.dist), opts.ReportDistance); err != nil {
					return count, err
				}
			}
		}
		if !any {
			if err := emitNull(a); err != nil {
				return count, err
			}
		}
	}

	// midStreamAbort means upstream's exit(1) fired inside the sweep before the
	// output buffer's destructor flush, so only the already-committed prefix is
	// visible on stdout.
	midStreamAbort := false
	if validateErr != nil {
		if ve, ok := validateErr.(*validationError); ok {
			midStreamAbort = ve.midStream
		}
	}
	var payload []byte
	if midStreamAbort {
		payload = out.committed()
	} else {
		payload = out.all()
	}
	if _, err := writer.Write(payload); err != nil {
		return count, fmt.Errorf("error writing output: %w", err)
	}
	if validateErr != nil {
		return count, validateErr
	}
	return count, nil
}

// validationFileNames assembles the file-name labels for the validator: the
// query name at index 0, then each database name in -b order. Missing names
// fall back to generic placeholders so the validator never panics.
func validationFileNames(opts Options, numDBs int) []string {
	names := make([]string, 0, numDBs+1)
	if opts.QueryName != "" {
		names = append(names, opts.QueryName)
	} else {
		names = append(names, "A")
	}
	for i := 0; i < numDBs; i++ {
		if i < len(opts.DBNames) && opts.DBNames[i] != "" {
			names = append(names, opts.DBNames[i])
		} else {
			names = append(names, fmt.Sprintf("B[%d]", i))
		}
	}
	return names
}

// allCand carries an -mdb all candidate hit with its source database index.
type allCand struct {
	dbIdx int
	b     *Row
	abs   int64
	neg   bool
}

// emitAllMode implements `-mdb all`: gather every database's hits, sort by
// absolute distance (ties broken by B's lessThan order), then take the k
// closest, honouring the tie mode, mirroring CloseSweep::checkMultiDbs.
func emitAllMode(
	a *Row, dbs []*db, opts Options,
	reachable func(dbIdx int, chrom string) bool,
	signedForOutput func(int64) int64,
	emit func(a *Row, lbl string, bFields []string, dist int64, hasDist bool) error,
	emitNull func(a *Row) error,
) error {
	var cands []allCand
	for dbIdx, d := range dbs {
		if !reachable(dbIdx, a.Chrom) {
			continue
		}
		for _, h := range hitsForA(a, d, opts) {
			abs := h.dist
			neg := abs < 0
			if abs < 0 {
				abs = -abs
			}
			cands = append(cands, allCand{dbIdx: dbIdx, b: h.b, abs: abs, neg: neg})
		}
	}
	if len(cands) == 0 {
		// No databases had any chromosome match: emit a single null row.
		return emitNull(a)
	}

	// Sort by absolute distance, then by B's (chrom, start, end) order
	// (Record::lessThan). Upstream uses std::sort (not stable); the lessThan
	// tiebreak makes the order deterministic.
	sort.SliceStable(cands, func(i, j int) bool {
		if cands[i].abs != cands[j].abs {
			return cands[i].abs < cands[j].abs
		}
		return lessThan(cands[i].b, cands[j].b)
	})

	signed := func(c allCand) int64 {
		if c.neg {
			return signedForOutput(-c.abs)
		}
		return signedForOutput(c.abs)
	}
	label := func(dbIdx int) string {
		if opts.DBLabels != nil {
			return opts.DBLabels[dbIdx]
		}
		return strconv.Itoa(dbIdx + 1)
	}

	// Walk distinct absolute distances, applying the tie mode.
	k := opts.kVal()
	used := 0
	for i := 0; i < len(cands) && used < k; {
		dist := cands[i].abs
		numTies := 1
		for i+numTies < len(cands) && cands[i+numTies].abs == dist {
			numTies++
		}
		switch opts.TieBreak {
		case TieFirst:
			c := cands[i]
			if err := emit(a, label(c.dbIdx), c.b.Fields, signed(c), opts.ReportDistance); err != nil {
				return err
			}
			used++
		case TieLast:
			c := cands[i+numTies-1]
			if err := emit(a, label(c.dbIdx), c.b.Fields, signed(c), opts.ReportDistance); err != nil {
				return err
			}
			used++
		default:
			for j := i; j < i+numTies; j++ {
				c := cands[j]
				if err := emit(a, label(c.dbIdx), c.b.Fields, signed(c), opts.ReportDistance); err != nil {
					return err
				}
				used++
			}
		}
		i += numTies
	}
	return nil
}

// lessThan mirrors Record::lessThan for same-chromosome ordering. For -mdb all
// all candidates share A's chromosome, so the chromosome comparison is omitted.
func lessThan(a, b *Row) bool {
	if a.Chrom != b.Chrom {
		return a.Chrom < b.Chrom
	}
	if a.Start != b.Start {
		return a.Start < b.Start
	}
	return a.End < b.End
}

// flushBufThreshold mirrors RecordOutputMgr's flush trigger: with buffered
// output (the default) the 16 KiB buffer is flushed once it reaches 90% full
// (16384 * 0.9 = 14745.6, i.e. >= 14746 bytes). On a mid-stream exit(1) the
// unflushed remainder is lost, so only the committed prefix reaches stdout.
const flushBufThreshold = 14746

// flushBuffer models upstream's RecordOutputMgr output buffer: bytes are
// appended to a pending buffer and committed in whole-buffer chunks once the
// pending size crosses flushBufThreshold (checked after each emitted record,
// matching upstream's `if (needsFlush()) flush()` at every newline). committed
// returns only the flushed prefix (what survives a mid-stream abort); all
// returns the full output (a clean run or a destructor-time abort flushes the
// remainder).
type flushBuffer struct {
	out     []byte // committed (flushed) bytes.
	pending []byte // buffered, not yet flushed.
}

func newFlushBuffer() *flushBuffer { return &flushBuffer{} }

func (b *flushBuffer) writeString(s string) { b.pending = append(b.pending, s...) }
func (b *flushBuffer) writeByte(c byte)     { b.pending = append(b.pending, c) }

// endRecord mirrors the per-record flush check: commit the pending buffer when
// it has reached the flush threshold.
func (b *flushBuffer) endRecord() {
	if len(b.pending) >= flushBufThreshold {
		b.out = append(b.out, b.pending...)
		b.pending = b.pending[:0]
	}
}

// committed returns the flushed prefix only.
func (b *flushBuffer) committed() []byte { return b.out }

// all returns every emitted byte (committed plus pending).
func (b *flushBuffer) all() []byte {
	res := make([]byte, 0, len(b.out)+len(b.pending))
	res = append(res, b.out...)
	res = append(res, b.pending...)
	return res
}

// writeRowBuf writes one output row into the flush buffer: A's columns, an
// optional database-label column, the chosen B's columns, and (when hasDist)
// the distance column, then runs the per-record flush check.
func writeRowBuf(b *flushBuffer, a *Row, printDBCol bool, label string, bFields []string, dist int64, hasDist bool) {
	b.writeString(strings.Join(a.Fields, "\t"))
	if printDBCol {
		b.writeByte('\t')
		b.writeString(label)
	}
	b.writeByte('\t')
	b.writeString(strings.Join(bFields, "\t"))
	if hasDist {
		b.writeByte('\t')
		b.writeString(strconv.FormatInt(dist, 10))
	}
	b.writeByte('\n')
	b.endRecord()
}
