// Chromosome naming-convention warning, mirroring upstream
// ContextBase::testNameConventions. When the input files mix "chr"-prefixed and
// bare chromosome names, or mix a leading zero after "chr" (chr1 vs chr01),
// upstream prints a one-time stderr WARNING naming the first offending record.
// It is suppressed by -nonamecheck and never affects the data output.
package bedintersect

import (
	"io"
	"strings"
)

// emitNameWarning writes the chromosome naming-convention warning (if any) to
// opts.Warnings, mirroring upstream's one-time stderr WARNING. It is a no-op
// when -nonamecheck is set, no inconsistency exists, or no warnings sink was
// provided. Upstream prints the message twice; that is reproduced so full
// stderr output matches.
func emitNameWarning(opts IntersectOptions, aRecords, bRecords []*inRecord) {
	if opts.Warnings == nil {
		return
	}
	msg := nameConventionWarning(opts, aRecords, bRecords)
	if msg == "" {
		return
	}
	// "<< msg << endl" appends a trailing newline; upstream emits the block
	// twice (trip plus reprint), so write it twice for byte parity.
	io.WriteString(opts.Warnings, msg+"\n")
	io.WriteString(opts.Warnings, msg+"\n")
}

// tristate models upstream's UNTESTED/YES/NO convention state.
type tristate int

const (
	untested tristate = iota
	yes
	no
)

// hasChrPrefix reports whether a chromosome name begins with "chr"
// (case-insensitive), mirroring Record::hasChrInChromName.
func hasChrPrefix(chrom string) bool {
	if len(chrom) < 3 {
		return false
	}
	c := strings.ToLower(chrom[:3])
	return c == "chr"
}

// hasLeadingZero reports whether the chromosome name has a '0' immediately after
// a "chr" prefix (e.g. chr01), mirroring Record::hasLeadingZeroInChromName.
func hasLeadingZero(chrom string) bool {
	return len(chrom) >= 4 && chrom[3] == '0' && hasChrPrefix(chrom)
}

// nameChecker reproduces upstream's running naming-convention state across all
// processed records, tripping (once) the first time a record's convention
// disagrees with the established global convention.
type nameChecker struct {
	allChr     tristate
	allZero    tristate
	tripped    bool
	warningMsg string
	disabled   bool
}

// check processes one record from file named filename. On the first
// inconsistency it records the upstream warning message and trips, ignoring all
// later records (and unmapped BAM records, which upstream skips).
func (n *nameChecker) check(rec *inRecord, filename string) {
	if n.disabled || n.tripped || rec.unmapped {
		return
	}
	hasChr := hasChrPrefix(rec.chrom)
	if (n.allChr == yes && !hasChr) || (n.allChr == no && hasChr) {
		n.trip(rec, filename, " has inconsistent naming convention for record:\n")
		return
	}
	if n.allChr == untested {
		n.allChr = boolToTri(hasChr)
	}

	zero := hasLeadingZero(rec.chrom)
	if (n.allZero == yes && !zero) || (n.allZero == no && zero) {
		n.trip(rec, filename, " has a record where naming convention (leading zero) is inconsistent with other files:\n")
		return
	}
	if n.allZero == untested {
		n.allZero = boolToTri(zero)
	}
}

// trip builds upstream's exact warning text (message + the offending record
// line) and marks the checker tripped.
func (n *nameChecker) trip(rec *inRecord, filename, message string) {
	n.warningMsg = "***** WARNING: File " + filename + message + rec.line + "\n"
	n.tripped = true
}

func boolToTri(b bool) tristate {
	if b {
		return yes
	}
	return no
}

// nameConventionWarning returns the upstream naming-convention warning message
// (empty when none), processing records in the order upstream does: with
// -sorted the chromsweep interleaves A and B, which for these inputs reduces to
// A then B; without -sorted the database (B) is loaded before A is swept, so B
// is processed first. -nonamecheck suppresses the check entirely.
func nameConventionWarning(opts IntersectOptions, aRecords, bRecords []*inRecord) string {
	if opts.NoNameCheck {
		return ""
	}
	nc := &nameChecker{}
	aName := opts.NameA
	bNameFor := func(id int) string {
		if id < len(opts.FilePaths) {
			return opts.FilePaths[id]
		}
		return opts.NameB
	}
	process := func(recs []*inRecord, name func(*inRecord) string) {
		for _, r := range recs {
			nc.check(r, name(r))
			if nc.tripped {
				return
			}
		}
	}
	bName := func(r *inRecord) string { return bNameFor(r.dbID) }
	aNameFn := func(*inRecord) string { return aName }
	if opts.Sorted {
		process(aRecords, aNameFn)
		if !nc.tripped {
			process(bRecords, bName)
		}
	} else {
		process(bRecords, bName)
		if !nc.tripped {
			process(aRecords, aNameFn)
		}
	}
	return nc.warningMsg
}
