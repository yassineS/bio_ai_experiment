// -sorted input-order validation, mirroring upstream
// FileRecordMgr::testInputSortOrder. With -sorted, every input file must be
// coordinate-sorted; with -g a genome file additionally fixes the chromosome
// order. A violation aborts with the same message upstream prints (and exit 1).
package bedintersect

import (
	"errors"
	"fmt"
	"math"
)

// sortError carries upstream's exact out-of-order message plus the offending
// record line, so the CLI can print it verbatim to stderr and exit 1.
type sortError struct {
	msg string
}

func (e *sortError) Error() string { return e.msg }

// IsSortError reports whether err is an upstream-style -sorted input-order
// violation, whose message must be printed verbatim (it already carries the
// "Error: ..." prefix and the offending record line).
func IsSortError(err error) bool {
	_, ok := err.(*sortError)
	return ok
}

// fieldCountError mirrors upstream's exact two-line message for a data line with
// the wrong number of fields (e.g. an extra trailing tab or a stray DOS
// newline). It is surfaced verbatim by the CLI.
type fieldCountError struct{}

func (e *fieldCountError) Error() string {
	return "Error: Type checker found wrong number of fields while tokenizing data line.\n" +
		"Perhaps you have extra TAB at the end of your line? Check with \"cat -t\""
}

// IsVerbatimError reports whether err carries an upstream-format message that the
// CLI must print verbatim (with no added "Error finding intersections:" prefix).
// It unwraps the fmt.Errorf("...: %w") chains the readers build around the
// underlying field-count error.
func IsVerbatimError(err error) bool {
	if IsSortError(err) {
		return true
	}
	var fce *fieldCountError
	return errors.As(err, &fce)
}

// VerbatimMessage returns the verbatim upstream-format message for an error the
// CLI must print as-is. For a wrapped field-count error it unwraps to the bare
// upstream two-line text; otherwise it returns the error's own message.
func VerbatimMessage(err error) string {
	var fce *fieldCountError
	if errors.As(err, &fce) {
		return fce.Error()
	}
	return err.Error()
}

// checkSorted runs the -sorted input-order validation against the A records and
// then the B records, when opts.Sorted is set. It is a no-op otherwise. B
// records are validated per-file (using the dbID tag) so a multi-database -b
// reports the offending file by name.
func checkSorted(opts IntersectOptions, aRecords, bRecords []*inRecord) error {
	if !opts.Sorted {
		return nil
	}
	if err := validateSortOrder(aRecords, opts.NameA, opts.GenomeOrder, opts.GenomeFile); err != nil {
		return err
	}
	// Validate each B file's records in isolation, in file order.
	maxID := 0
	for _, b := range bRecords {
		if b.dbID > maxID {
			maxID = b.dbID
		}
	}
	for id := 0; id <= maxID; id++ {
		var fileRecs []*inRecord
		for _, b := range bRecords {
			if b.dbID == id {
				fileRecs = append(fileRecs, b)
			}
		}
		name := opts.NameB
		if id < len(opts.FilePaths) {
			name = opts.FilePaths[id]
		}
		if err := validateSortOrder(fileRecs, name, opts.GenomeOrder, opts.GenomeFile); err != nil {
			return err
		}
	}
	return nil
}

// validateSortOrder checks that recs are coordinate-sorted, mirroring
// FileRecordMgr::testInputSortOrder. filename names the file for the error
// message; genomeOrder, when non-empty, gives the required chromosome order from
// the -g genome file named genomeFile. It returns a *sortError on the first
// violation.
func validateSortOrder(recs []*inRecord, filename string, genomeOrder []string, genomeFile string) error {
	genomeIdx := map[string]int{}
	for i, c := range genomeOrder {
		genomeIdx[c] = i
	}
	hasGenome := len(genomeOrder) > 0

	found := map[string]bool{}
	prevChrom := ""
	prevStart := int64(math.MaxInt64)
	prevChromID := -1

	for _, r := range recs {
		currStart := int64(r.start)
		if r.start == r.end { // zero-length record: upstream bumps the start by 1
			currStart++
		}
		if r.chrom != prevChrom {
			if found[r.chrom] {
				// Chromosome reappears after a different one: out of order.
				return outOfOrderErr(filename, r)
			}
			if hasGenome {
				currChromID, ok := genomeIdx[r.chrom]
				if !ok {
					currChromID = len(genomeOrder) // unknown chroms sort last
				}
				if currChromID < prevChromID {
					return genomeOrderErr(filename, genomeFile, r)
				}
				prevChromID = currChromID
			}
			found[r.chrom] = true
			prevChrom = r.chrom
			prevStart = int64(math.MaxInt64)
		} else if currStart < prevStart {
			return outOfOrderErr(filename, r)
		}
		prevStart = currStart
	}
	return nil
}

// outOfOrderErr builds the non-genome "out of order record" message.
func outOfOrderErr(filename string, r *inRecord) error {
	return &sortError{msg: fmt.Sprintf(
		"Error: Sorted input specified, but the file %s has the following out of order record\n%s",
		filename, r.line)}
}

// genomeOrderErr builds the genome-file "different sort order" message.
func genomeOrderErr(filename, genomeFile string, r *inRecord) error {
	return &sortError{msg: fmt.Sprintf(
		"Error: Sorted input specified, but the file %s has the following record with a different sort order than the genomeFile %s\n%s",
		filename, genomeFile, r.line)}
}
