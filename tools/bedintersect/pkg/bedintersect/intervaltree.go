// Package bedintersect provides functionality to find intersecting intervals between two BED files.
//
// The IntervalTree type used here was lifted into `pkg/bioformats/bed` so it
// can be reused by bedcoverage, bedmap, and any other tool that needs
// efficient overlap queries. The names below are kept as thin type aliases so
// existing imports of bedintersect.IntervalTree / NewIntervalTree continue to
// compile.
package bedintersect

import (
	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/bed"
)

// IntervalNode is an alias for bed.IntervalNode (see pkg/bioformats/bed/intervaltree.go).
type IntervalNode = bed.IntervalNode

// IntervalTree is an alias for bed.IntervalTree (see pkg/bioformats/bed/intervaltree.go).
type IntervalTree = bed.IntervalTree

// NewIntervalTree builds a balanced interval tree from the given records.
// Thin wrapper over bed.NewIntervalTree.
func NewIntervalTree(intervals []*bed.Record) *IntervalTree {
	return bed.NewIntervalTree(intervals)
}
