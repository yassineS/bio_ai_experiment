package bedclosest

// This file ports bedtools' cross-file chromosome sort-order and
// naming-convention validation, which the `closest` sortAndNaming sub-suite
// (upstream test ids closest.t01-t06, t10/t11, t14-t17, t19) exercises.
//
// Upstream performs the checks while it streams records through the chromosome
// sweep (NewChromSweep). Two cooperating pieces of state drive it:
//
//   - ContextBase::testNameConventions (src/utils/Contexts/ContextBase.cpp)
//     emits the "inconsistent naming convention" (chr-prefix) and
//     "(leading zero)" WARNINGs the first time any file disagrees with the
//     convention established by the first record seen across all files.
//
//   - NewChromSweep::testChromOrder / testThatAllDbChromsExistInQuery /
//     testLexicoQueryAfterDb (src/utils/NewChromsweep/NewChromsweep.cpp)
//     run a streaming "_lexicoAssumed / _lexicoDisproven" state machine that
//     emits the out-of-order and lexicographic ERRORs and, after all output
//     is produced, the "Database file ... contains chromosome X, but the
//     query file does not" ERROR.
//
// The port loads every file up front and buckets it by chromosome, so this
// file faithfully *re-simulates* the upstream sweep's record advancement
// (init -> per-query next()/masterScan()/chromChange() -> closeOut() ->
// destructor) purely to reproduce the exact sequence of testChromOrder and
// testNameConventions calls. The simulation never computes hits; it only
// determines which WARNING/ERROR fires, in what order, and how many query
// records were already emitted to stdout when a mid-stream ERROR aborts the
// run (so the caller can reproduce upstream's partial-stdout-then-exit-1
// behaviour byte-for-byte).

import (
	"fmt"
	"io"
	"sort"
)

// validationError is returned by ClosestMulti when the upstream sort-order /
// naming-convention validation aborts the run. It carries the exact upstream
// stderr text (already written to the warn writer).
type validationError struct {
	msg string
	// midStream is true when upstream aborted inside the chromosome sweep
	// (exit(1) before the output buffer's destructor flush), so only the
	// already-committed stdout prefix is visible. It is false for the
	// destructor-time "Database file ... contains chromosome" abort, which
	// fires after all output has been flushed.
	midStream bool
}

// Error implements the error interface.
func (e *validationError) Error() string { return e.msg }

// IsValidationError reports whether err is a cross-file sort-order /
// naming-convention validation ERROR from ClosestMulti. Such errors have
// already written their exact upstream stderr text to the configured
// WarnWriter, so the CLI must not re-print them with a generic prefix; it only
// needs to propagate the non-zero exit code.
func IsValidationError(err error) bool {
	_, ok := err.(*validationError)
	return ok
}

// convState classifies whether a file uses a given naming convention,
// mirroring ContextBase's UNTESTED/YES/NO tri-state.
type convState int

const (
	convUntested convState = iota
	convYes
	convNo
)

// hasChrInChromName mirrors Record::hasChrInChromName: the first three
// characters are c, h, r (case-insensitive).
func hasChrInChromName(chrom string) bool {
	if len(chrom) < 3 {
		return false
	}
	c0, c1, c2 := chrom[0], chrom[1], chrom[2]
	return (c0 == 'c' || c0 == 'C') &&
		(c1 == 'h' || c1 == 'H') &&
		(c2 == 'r' || c2 == 'R')
}

// hasLeadingZeroInChromName mirrors Record::hasLeadingZeroInChromName: the
// fourth character (index 3) is '0', provided the name carries the "chr"
// convention (either already established by the caller or detected here).
func hasLeadingZeroInChromName(chrom string, chrKnown bool) bool {
	return len(chrom) >= 4 && chrom[3] == '0' && (chrKnown || hasChrInChromName(chrom))
}

// orderTrack records the first-seen order of chromosomes within one file,
// mirroring NewChromSweep::_fileTracks (a map<string,int> keyed by chrom name
// with the insertion index as value).
type orderTrack struct {
	order map[string]int
}

func newOrderTrack() *orderTrack { return &orderTrack{order: make(map[string]int)} }

// findOrInsert mirrors findChromOrder: returns the chrom's order index,
// assigning the next index the first time the chrom is seen.
func (t *orderTrack) findOrInsert(chrom string) int {
	if v, ok := t.order[chrom]; ok {
		return v
	}
	v := len(t.order)
	t.order[chrom] = v
	return v
}

// chromsLexico returns the file's chromosome names in lexicographic order,
// matching the iteration order of upstream's std::map<string,int> in
// testThatAllDbChromsExistInQuery (which reports the first such chrom).
func (t *orderTrack) chromsLexico() []string {
	out := make([]string, 0, len(t.order))
	for c := range t.order {
		out = append(out, c)
	}
	sort.Strings(out)
	return out
}

// sweepValidator re-simulates the upstream NewChromSweep validation state
// across the query file and every database file. It mirrors the C++ members
// one-for-one so the testChromOrder / testNameConventions call sequence (and
// therefore which message fires first) matches upstream exactly.
type sweepValidator struct {
	// Input records, in file order: query plus one slice per database.
	query []*Row
	dbs   [][]*Row

	// File names as they should appear in messages: index 0 is the query, the
	// rest are databases in -b order. Mirrors getInputFileName(fileIdx).
	fileNames []string
	numFiles  int // len(fileNames) == numDBs + 1

	warn io.Writer // upstream's stderr.

	// --- naming-convention state (ContextBase) ---
	nameWarningTripped      bool
	fileHasChr              []convState
	allFilesHaveChr         convState
	fileHasLeadingZero      []convState
	allFilesHaveLeadingZero convState

	// --- sort-order state (NewChromSweep) ---
	fileTracks      []*orderTrack
	filePrevChrom   []string // valid only when filePrevSet[i] is true.
	filePrevSet     []bool   // mirrors _filePrevChrom[i] != NULL.
	lexicoDisproven bool
	lexicoAssumed   bool

	// --- abort bookkeeping ---
	emitted   int    // query records processed before an abort.
	aborted   bool   // an ERROR fired.
	midStream bool   // the ERROR fired mid-stream (vs in the destructor).
	abortMsg  string // its stderr text.

	// reachable[dbIdx][chrom] records, per database, the query chromosomes the
	// sweep actually advanced the database onto. When lexicographic ordering
	// has been disproven a database can get "stuck" on a chromosome that sorts
	// after the query's, so later query chromosomes never see the database's
	// matching records and must report null even though a naive per-chromosome
	// bucket would match them. For consistently sorted inputs every matching
	// chromosome is reachable, so this gate is a no-op.
	reachable []map[string]bool
}

// newSweepValidator builds a validator for the given query rows, database row
// slices, and the file-name labels (query first, then databases).
func newSweepValidator(query []*Row, dbs [][]*Row, fileNames []string, warn io.Writer) *sweepValidator {
	n := len(dbs) + 1
	v := &sweepValidator{
		query:                   query,
		dbs:                     dbs,
		fileNames:               fileNames,
		numFiles:                n,
		warn:                    warn,
		fileHasChr:              make([]convState, n),
		allFilesHaveChr:         convUntested,
		fileHasLeadingZero:      make([]convState, n),
		allFilesHaveLeadingZero: convUntested,
		fileTracks:              make([]*orderTrack, n),
		filePrevChrom:           make([]string, n),
		filePrevSet:             make([]bool, n),
	}
	for i := 0; i < n; i++ {
		v.fileTracks[i] = newOrderTrack()
	}
	v.reachable = make([]map[string]bool, len(dbs))
	for i := range dbs {
		v.reachable[i] = make(map[string]bool)
	}
	return v
}

// markReachable records that database dbIdx advanced onto chromosome chrom
// during the sweep, so its records there are eligible to produce hits.
func (v *sweepValidator) markReachable(dbIdx int, chrom string) {
	v.reachable[dbIdx][chrom] = true
}

// dbReachable reports whether database dbIdx ever advanced onto chromosome
// chrom. It is consulted by the hit engine to null out matches the upstream
// sweep could not have produced because the database was stuck earlier.
func (v *sweepValidator) dbReachable(dbIdx int, chrom string) bool {
	return v.reachable[dbIdx][chrom]
}

// dbFileIdx maps a 0-based database index to its file index (databases follow
// the query at index 0).
func dbFileIdx(dbIdx int) int { return dbIdx + 1 }

// testNameConventions mirrors ContextBase::testNameConventions: it records the
// first chr/leading-zero convention seen and emits the WARNING when a later
// record disagrees. The warning is emitted at most once for the whole run.
func (v *sweepValidator) testNameConventions(rec *Row, fileIdx int) {
	if v.nameWarningTripped || rec == nil {
		return
	}
	chrom := rec.Chrom

	// chr-prefix convention.
	hasChr := hasChrInChromName(chrom)
	if v.fileHasChr[fileIdx] == convUntested {
		v.fileHasChr[fileIdx] = boolToConv(hasChr)
	}
	if (v.allFilesHaveChr == convYes && !hasChr) || (v.allFilesHaveChr == convNo && hasChr) {
		v.nameConventionWarning(rec, fileIdx, " has inconsistent naming convention for record:\n")
	}
	if v.allFilesHaveChr == convUntested {
		v.allFilesHaveChr = boolToConv(hasChr)
	}

	// leading-zero convention.
	zero := hasLeadingZeroInChromName(chrom, hasChr)
	if v.fileHasLeadingZero[fileIdx] == convUntested {
		v.fileHasLeadingZero[fileIdx] = boolToConv(zero)
	}
	if (v.allFilesHaveLeadingZero == convYes && !zero) || (v.allFilesHaveLeadingZero == convNo && zero) {
		v.nameConventionWarning(rec, fileIdx,
			" has a record where naming convention (leading zero) is inconsistent with other files:\n")
	}
	if v.allFilesHaveLeadingZero == convUntested {
		v.allFilesHaveLeadingZero = boolToConv(zero)
	}
}

// boolToConv maps a detected-convention boolean to the YES/NO tri-state.
func boolToConv(b bool) convState {
	if b {
		return convYes
	}
	return convNo
}

// nameConventionWarning mirrors ContextBase::nameConventionWarning: it builds
// the WARNING message, prints it to stderr (with the trailing endl upstream
// adds), and trips the once-only flag.
func (v *sweepValidator) nameConventionWarning(rec *Row, fileIdx int, message string) {
	msg := "***** WARNING: File " + v.fileNames[fileIdx] + message + printRecord(rec) + "\n"
	v.nameWarningTripped = true
	// Upstream prints the message followed by std::endl (an extra newline).
	fmt.Fprintln(v.warn, msg)
}

// printRecord mirrors Record::print for the Bed3 form used in error/warning
// messages: "chrom\tstart\tend".
func printRecord(rec *Row) string {
	return fmt.Sprintf("%s\t%d\t%d", rec.Chrom, rec.Start, rec.End)
}

// testChromOrder mirrors NewChromSweep::testChromOrder. It returns false when
// it triggers an abort (so callers stop advancing the simulation).
func (v *sweepValidator) testChromOrder(rec *Row, fileIdx int) bool {
	if v.aborted {
		return false
	}
	if rec == nil {
		return true
	}
	chrom := rec.Chrom
	v.fileTracks[fileIdx].findOrInsert(chrom)

	if !v.filePrevSet[fileIdx] {
		v.filePrevChrom[fileIdx] = chrom
		v.filePrevSet[fileIdx] = true
		return true
	}
	if v.filePrevChrom[fileIdx] == chrom {
		return true
	}
	prevChrom := v.filePrevChrom[fileIdx]
	v.filePrevChrom[fileIdx] = chrom

	if v.verifyChromOrderMismatch(chrom, prevChrom, fileIdx) {
		v.abort(fmt.Sprintf(
			"ERROR: chromosome sort ordering for file %s is inconsistent with other files. Record was:\n%s\n",
			v.fileNames[fileIdx], printRecord(rec)), true)
		return false
	}

	if !v.lexicoDisproven && chrom < prevChrom {
		if v.lexicoAssumed {
			v.abort(fmt.Sprintf(
				"ERROR: Sort order was unspecified, and file %s is not sorted lexicographically.\n"+
					"       Please rerun with the -g option for a genome file.\n"+
					"       See documentation for details.\n",
				v.fileNames[fileIdx]), true)
			return false
		}
		v.lexicoDisproven = true
	}
	return true
}

// verifyChromOrderMismatch mirrors NewChromSweep::verifyChromOrderMismatch:
// for every other file that has seen both chrom and prevChrom, report a
// mismatch if it ordered them the opposite way (curr before prev).
func (v *sweepValidator) verifyChromOrderMismatch(chrom, prevChrom string, skipFile int) bool {
	for i := 0; i < v.numFiles; i++ {
		if i == skipFile {
			continue
		}
		track := v.fileTracks[i]
		currOrder, ok := track.order[chrom]
		if !ok {
			continue
		}
		prevOrder, ok := track.order[prevChrom]
		if !ok {
			continue
		}
		if currOrder < prevOrder {
			return true
		}
	}
	return false
}

// testLexicoQueryAfterDb mirrors NewChromSweep::testLexicoQueryAfterDb: when
// the query file lacks the database chrom, fall back to a lexicographic
// comparison and, the first time the query is lexicographically greater, set
// the _lexicoAssumed flag.
func (v *sweepValidator) testLexicoQueryAfterDb(queryChrom, dbChrom string) bool {
	if v.lexicoDisproven {
		return false
	}
	queryGreater := queryChrom > dbChrom
	if !v.lexicoAssumed && queryGreater {
		v.lexicoAssumed = true
	}
	return queryGreater
}

// abort records an ERROR: it writes the message to stderr and captures it so
// the caller can return a validationError. midStream distinguishes a sweep-time
// abort (output buffer not flushed) from a destructor-time abort (all output
// already flushed).
func (v *sweepValidator) abort(msg string, midStream bool) {
	if v.aborted {
		return
	}
	v.aborted = true
	v.midStream = midStream
	v.abortMsg = msg
	fmt.Fprint(v.warn, msg)
}

// dbCursor mirrors the per-database streaming pointer state of the sweep: the
// index of the next unread record, the currently fetched record (mirroring
// _currDbRecs[i], nil when none), and the cache of records held for the
// current query chrom (mirroring _caches[i], whose contents do not themselves
// trigger sort-order checks).
type dbCursor struct {
	rows  []*Row
	pos   int    // index of next record to fetch.
	cur   *Row   // current record (nil == NULL).
	cache []*Row // held records (chrom-order-neutral).
}

// nextDB mirrors nextRecord(false, i): fetch the next database record, or set
// cur to nil at EOF. Returns whether a record was fetched.
func (c *dbCursor) nextDB() bool {
	if c.pos < len(c.rows) {
		c.cur = c.rows[c.pos]
		c.pos++
		return true
	}
	c.cur = nil
	return false
}

// dbEOF reports whether the database file is exhausted (no more records to
// fetch), mirroring FileRecordMgr::eof for the streamed db.
func (c *dbCursor) dbEOF() bool { return c.pos >= len(c.rows) }

// runSweep replays the upstream NewChromSweep lifecycle (init, the per-query
// next()/masterScan()/chromChange() loop, closeOut, and the destructor's
// testThatAllDbChromsExistInQuery) over the pre-loaded records, purely to
// reproduce the WARNING/ERROR sequence and the count of query records emitted
// before any mid-stream abort. It returns a *validationError when the run
// terminates with a non-zero exit (after writing the exact stderr text), and
// sets v.emitted to the number of query records whose output upstream produced
// before that abort. For runs that abort only in the destructor (db chrom not
// in query), v.emitted equals the full query count, matching upstream's
// "produce all output, then exit 1".
func (v *sweepValidator) runSweep() error {
	cursors := make([]*dbCursor, len(v.dbs))
	for i := range v.dbs {
		cursors[i] = &dbCursor{rows: v.dbs[i]}
	}

	// init(): fetch the first record of each db and test its chrom order.
	for i := range cursors {
		cursors[i].nextDB()
		if !v.testChromOrder(cursors[i].cur, dbFileIdx(i)) {
			return v.finishAbort()
		}
	}

	// The per-query loop. We track currQueryRec via qi; -1 means before the
	// first record. needTestSortOrder mirrors the first-query special case.
	qi := -1
	prevQueryChrom := ""
	havePrevQueryChrom := false

	for {
		// next(): advance to the next query record.
		needTestSortOrder := qi < 0
		qi++
		if qi >= len(v.query) {
			break // query EOF.
		}
		curQuery := v.query[qi]

		if needTestSortOrder {
			if !v.testChromOrder(curQuery, 0) {
				return v.finishAbort()
			}
		}

		// allCurrDBrecsNull && allCachesEmpty && !runToQueryEnd short-circuit:
		// upstream returns false (stops the sweep) without printing this query
		// record, deferring a final testChromOrder to closeOut. closest does
		// not set runToQueryEnd, so honour the short-circuit.
		if v.allDBNull(cursors) && v.allCachesEmpty(cursors) {
			// _testLastQueryRec = true; closeOut will test this record.
			v.closeOut(curQuery, qi, cursors, true)
			return v.finishAbort()
		}

		curQueryChrom := curQuery.Chrom

		// masterScan(): per-db chromChange + advance.
		for i := range cursors {
			if v.dbFinished(cursors[i]) {
				continue
			}
			cont := v.chromChange(i, cursors, curQuery, curQueryChrom, prevQueryChrom, havePrevQueryChrom)
			if v.aborted {
				return v.finishAbort()
			}
			if !cont {
				continue
			}
			v.advanceDB(i, cursors, curQuery)
		}

		// This query record's output is produced by the caller now.
		v.emitted = qi + 1

		prevQueryChrom = curQueryChrom
		havePrevQueryChrom = true
	}

	// closeOut(true): test any remaining records for sort order.
	v.closeOut(nil, len(v.query), cursors, false)
	if v.aborted {
		return v.finishAbort()
	}

	// Destructor: testThatAllDbChromsExistInQuery.
	v.testThatAllDbChromsExistInQuery()
	if v.aborted {
		// Upstream produced all output before this destructor error.
		v.emitted = len(v.query)
		return v.finishAbort()
	}
	return nil
}

// finishAbort converts a recorded abort into a *validationError (or nil).
func (v *sweepValidator) finishAbort() error {
	if v.aborted {
		return &validationError{msg: v.abortMsg, midStream: v.midStream}
	}
	return nil
}

// allDBNull mirrors allCurrDBrecsNull.
func (v *sweepValidator) allDBNull(cursors []*dbCursor) bool {
	for _, c := range cursors {
		if c.cur != nil {
			return false
		}
	}
	return true
}

// allCachesEmpty mirrors allCachesEmpty.
func (v *sweepValidator) allCachesEmpty(cursors []*dbCursor) bool {
	for _, c := range cursors {
		if len(c.cache) != 0 {
			return false
		}
	}
	return true
}

// dbFinished mirrors dbFinished: the db is done when its current record is nil
// and its cache is empty.
func (v *sweepValidator) dbFinished(c *dbCursor) bool {
	return c.cur == nil && len(c.cache) == 0
}

// chromChange mirrors NewChromSweep::chromChange's validation-relevant
// behaviour. It runs the name-convention/chrom-order tests for the query (on a
// chrom change) and the current db record, then fast-forwards the database past
// the query when the query is ahead, testing each new db chrom. It returns
// false when the database is ahead of the query (the caller then skips the
// inner advance loop), true otherwise.
func (v *sweepValidator) chromChange(
	dbIdx int, cursors []*dbCursor, curQuery *Row,
	curQueryChrom, prevQueryChrom string, havePrevQueryChrom bool,
) bool {
	c := cursors[dbIdx]
	dbRec := c.cur

	if curQuery != nil && (!havePrevQueryChrom || curQueryChrom != prevQueryChrom) {
		v.testNameConventions(curQuery, 0)
		if !v.testChromOrder(curQuery, 0) {
			return false
		}
	}
	if dbRec != nil {
		v.testNameConventions(dbRec, dbFileIdx(dbIdx))
		if !v.testChromOrder(dbRec, dbFileIdx(dbIdx)) {
			return false
		}
	}

	// If query and db are on the same chrom, keep scanning (return false in
	// upstream means "don't skip"; the C++ returns false here too -> caller
	// proceeds to scan). Upstream returns false to mean "do not continue/skip".
	if dbRec != nil && curQuery != nil && curQuery.Chrom == dbRec.Chrom {
		v.markReachable(dbIdx, curQuery.Chrom)
		return true
	}
	if dbRec == nil || curQuery == nil {
		return true
	}

	if v.queryChromAfterDbRec(dbIdx, cursors, curQuery, dbRec) {
		// Query is ahead: fast-forward the db, testing each new chrom.
		oldDbChrom := dbRec.Chrom
		for c.cur != nil && v.queryChromAfterDbRec(dbIdx, cursors, curQuery, c.cur) {
			if !c.nextDB() {
				break
			}
			if c.cur == nil {
				break
			}
			newChrom := c.cur.Chrom
			if newChrom != oldDbChrom {
				if !v.testChromOrder(c.cur, dbFileIdx(dbIdx)) {
					return false
				}
				oldDbChrom = newChrom
			}
		}
		c.cache = nil // clearCache
		// After fast-forwarding, the database may now sit on the query's chrom.
		if c.cur != nil && c.cur.Chrom == curQuery.Chrom {
			v.markReachable(dbIdx, curQuery.Chrom)
		}
		return true
	}
	// Database is ahead of the query: skip this db for the current query.
	return false
}

// advanceDB mirrors masterScan's inner advance loop: it consumes db records up
// to and past the query on the same chrom, moving them to the cache or
// dropping them. None of this triggers sort-order checks.
func (v *sweepValidator) advanceDB(dbIdx int, cursors []*dbCursor, curQuery *Row) {
	c := cursors[dbIdx]
	// scanCache: drop cached records that no longer match the query chrom.
	kept := c.cache[:0]
	for _, rec := range c.cache {
		if curQuery.Chrom == rec.Chrom && !after(curQuery, rec) {
			kept = append(kept, rec)
		}
	}
	c.cache = kept

	for c.cur != nil && curQuery.Chrom == c.cur.Chrom && !after(c.cur, curQuery) {
		v.markReachable(dbIdx, curQuery.Chrom)
		if after(curQuery, c.cur) {
			// query is past this db rec: drop it.
		} else {
			c.cache = append(c.cache, c.cur)
		}
		c.nextDB()
	}
	// Records held in the cache for the current query chrom are also reachable.
	for _, rec := range c.cache {
		if rec.Chrom == curQuery.Chrom {
			v.markReachable(dbIdx, curQuery.Chrom)
			break
		}
	}
}

// after mirrors Record::after: a is after b when same chrom and a.start >=
// b.end.
func after(a, b *Row) bool {
	return a.Chrom == b.Chrom && a.Start >= b.End
}

// queryChromAfterDbRec mirrors NewChromSweep::queryChromAfterDbRec for the
// no-genome-file case: compare the global first-seen order of the query and db
// chroms within the query file's track, falling back to the lexicographic
// assumption when the query has not seen the db chrom.
func (v *sweepValidator) queryChromAfterDbRec(dbIdx int, cursors []*dbCursor, curQuery, dbRec *Row) bool {
	qChrom := curQuery.Chrom
	dbChrom := dbRec.Chrom
	qTrack := v.fileTracks[0]
	qOrder, ok := qTrack.order[qChrom]
	if !ok {
		// The query record's chrom was always inserted by testChromOrder
		// before masterScan runs, so this should not happen; insert to match
		// upstream's find semantics defensively.
		qOrder = qTrack.findOrInsert(qChrom)
	}
	dbOrder, ok := qTrack.order[dbChrom]
	if !ok {
		return v.testLexicoQueryAfterDb(qChrom, dbChrom)
	}
	return qOrder > dbOrder
}

// closeOut mirrors NewChromSweep::closeOut(true): after the per-query loop it
// re-tests the current/last query record and then streams every remaining db
// record through testChromOrder. testLastQueryRec controls whether the last
// query record still needs a sort-order test.
func (v *sweepValidator) closeOut(lastQuery *Row, _ int, cursors []*dbCursor, testLastQueryRec bool) {
	if v.aborted {
		return
	}
	if testLastQueryRec && lastQuery != nil {
		if !v.testChromOrder(lastQuery, 0) {
			return
		}
	}
	// Upstream then drains the rest of the query file; our query is fully
	// consumed by the loop, so there is nothing left to drain here.

	for i := range cursors {
		c := cursors[i]
		for !c.dbEOF() {
			if !v.testChromOrder(c.cur, dbFileIdx(i)) {
				return
			}
			c.nextDB()
		}
		if !v.testChromOrder(c.cur, dbFileIdx(i)) {
			return
		}
	}
}

// testThatAllDbChromsExistInQuery mirrors
// NewChromSweep::testThatAllDbChromsExistInQuery: once lexicographic order has
// been disproven, every database chrom must also appear in the query file, or
// upstream errors out (after all output is produced).
func (v *sweepValidator) testThatAllDbChromsExistInQuery() {
	if !v.lexicoDisproven {
		return
	}
	qTrack := v.fileTracks[0]
	for i := 1; i < v.numFiles; i++ {
		dbTrack := v.fileTracks[i]
		// Iterate db chroms in first-seen (insertion) order to match upstream's
		// map iteration determinism for the *which chrom* reported. Upstream's
		// std::map iterates in lexicographic key order, so replicate that.
		for _, chrom := range dbTrack.chromsLexico() {
			if chrom == "" {
				continue
			}
			if _, ok := qTrack.order[chrom]; !ok {
				v.abort(fmt.Sprintf(
					"ERROR: Database file %s contains chromosome %s, but the query file does not.\n"+
						"       Please rerun with the -g option for a genome file.\n"+
						"       See documentation for details.\n",
					v.fileNames[i], chrom), false)
				return
			}
		}
	}
}
