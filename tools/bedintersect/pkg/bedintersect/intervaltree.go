// Package bedintersect provides an interval tree for efficient interval queries.
package bedintersect

import (
	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/bed"
)

// IntervalNode represents a node in the interval tree.
type IntervalNode struct {
	Interval *bed.Record
	Max      int // Maximum end coordinate in subtree
	Left     *IntervalNode
	Right    *IntervalNode
}

// IntervalTree provides efficient interval overlap queries.
type IntervalTree struct {
	Root *IntervalNode
}

// NewIntervalTree creates a new interval tree from a slice of intervals.
func NewIntervalTree(intervals []*bed.Record) *IntervalTree {
	if len(intervals) == 0 {
		return &IntervalTree{}
	}

	// Build balanced tree recursively
	root := buildTree(intervals, 0, len(intervals)-1)
	return &IntervalTree{Root: root}
}

// buildTree recursively builds a balanced interval tree.
func buildTree(intervals []*bed.Record, start, end int) *IntervalNode {
	if start > end {
		return nil
	}

	mid := (start + end) / 2
	node := &IntervalNode{
		Interval: intervals[mid],
		Max:      intervals[mid].ChromEnd,
	}

	// Build left and right subtrees
	node.Left = buildTree(intervals, start, mid-1)
	node.Right = buildTree(intervals, mid+1, end)

	// Update max to be maximum of this node and children
	if node.Left != nil && node.Left.Max > node.Max {
		node.Max = node.Left.Max
	}
	if node.Right != nil && node.Right.Max > node.Max {
		node.Max = node.Right.Max
	}

	return node
}

// Query finds all intervals that overlap with the query interval.
func (tree *IntervalTree) Query(query *bed.Record) []*bed.Record {
	if tree.Root == nil {
		return nil
	}

	var results []*bed.Record
	tree.queryNode(tree.Root, query, &results)
	return results
}

// queryNode recursively searches the tree for overlapping intervals.
func (tree *IntervalTree) queryNode(node *IntervalNode, query *bed.Record, results *[]*bed.Record) {
	if node == nil {
		return
	}

	// If query starts after max end in this subtree, no overlap possible
	if query.ChromStart >= node.Max {
		return
	}

	// Check left subtree
	if node.Left != nil {
		tree.queryNode(node.Left, query, results)
	}

	// Check current node for overlap
	if intervalsOverlap(query, node.Interval) {
		*results = append(*results, node.Interval)
	}

	// Check right subtree only if query might overlap
	if node.Right != nil && query.ChromEnd > node.Interval.ChromStart {
		tree.queryNode(node.Right, query, results)
	}
}

// intervalsOverlap checks if two intervals overlap.
func intervalsOverlap(a, b *bed.Record) bool {
	return a.ChromStart < b.ChromEnd && a.ChromEnd > b.ChromStart
}
