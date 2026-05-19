// Interval tree for efficient overlap queries against a set of BED records.
//
// This is the same balanced augmented-binary-search-tree used by
// `tools/bedintersect`, lifted into the shared `bed` package so that other
// tools (bedcoverage, bedmap, ...) can share one implementation rather than
// each copying the same structure or importing bedintersect's package.
//
// The tree is built once from a slice of *Record (one per chromosome typically
// — callers index by Chrom first) and answers Query in amortised O(log n + k)
// where k is the number of overlapping records.

package bed

// IntervalNode is a single node in the augmented BST. Max is the maximum
// ChromEnd across all records in this subtree, used to prune queries.
type IntervalNode struct {
	Interval *Record
	Max      int
	Left     *IntervalNode
	Right    *IntervalNode
}

// IntervalTree provides efficient interval-overlap queries over a fixed set
// of BED records. Build it once with NewIntervalTree, then call Query.
//
// The tree does not track which chromosome each record belongs to. Callers
// that need cross-chromosome data should build one tree per Chrom.
type IntervalTree struct {
	Root *IntervalNode
}

// NewIntervalTree builds a balanced interval tree from the given records.
// The records slice should be sorted by ChromStart before calling for the
// resulting tree to be balanced.
func NewIntervalTree(intervals []*Record) *IntervalTree {
	if len(intervals) == 0 {
		return &IntervalTree{}
	}
	return &IntervalTree{Root: buildIntervalTree(intervals, 0, len(intervals)-1)}
}

// buildIntervalTree recursively builds a balanced subtree from intervals
// between indices start and end (inclusive).
func buildIntervalTree(intervals []*Record, start, end int) *IntervalNode {
	if start > end {
		return nil
	}
	mid := (start + end) / 2
	node := &IntervalNode{
		Interval: intervals[mid],
		Max:      intervals[mid].ChromEnd,
	}
	node.Left = buildIntervalTree(intervals, start, mid-1)
	node.Right = buildIntervalTree(intervals, mid+1, end)
	if node.Left != nil && node.Left.Max > node.Max {
		node.Max = node.Left.Max
	}
	if node.Right != nil && node.Right.Max > node.Max {
		node.Max = node.Right.Max
	}
	return node
}

// Query returns all records in the tree that overlap with query.
// Overlap uses half-open semantics: a.start < b.end && a.end > b.start.
// The Chrom field of the query is ignored — callers must select the right
// per-chromosome tree themselves.
func (t *IntervalTree) Query(query *Record) []*Record {
	if t.Root == nil {
		return nil
	}
	var results []*Record
	t.queryNode(t.Root, query, &results)
	return results
}

// queryNode recursively searches for overlapping intervals under node.
func (t *IntervalTree) queryNode(node *IntervalNode, query *Record, results *[]*Record) {
	if node == nil {
		return
	}
	if query.ChromStart >= node.Max {
		return
	}
	if node.Left != nil {
		t.queryNode(node.Left, query, results)
	}
	if intervalsOverlap(query, node.Interval) {
		*results = append(*results, node.Interval)
	}
	if node.Right != nil && query.ChromEnd > node.Interval.ChromStart {
		t.queryNode(node.Right, query, results)
	}
}

// intervalsOverlap is the half-open-interval overlap predicate used by the
// interval tree. It does NOT compare Chrom: the tree is per-chrom.
func intervalsOverlap(a, b *Record) bool {
	return a.ChromStart < b.ChromEnd && a.ChromEnd > b.ChromStart
}
