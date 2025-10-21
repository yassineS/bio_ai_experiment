package bedintersect

import (
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/bed"
)

func TestIntervalTreeBasic(t *testing.T) {
	intervals := []*bed.Record{
		{Chrom: "chr1", ChromStart: 100, ChromEnd: 200},
		{Chrom: "chr1", ChromStart: 300, ChromEnd: 400},
		{Chrom: "chr1", ChromStart: 500, ChromEnd: 600},
	}
	
	tree := NewIntervalTree(intervals)
	
	// Query overlapping with first interval
	query := &bed.Record{Chrom: "chr1", ChromStart: 150, ChromEnd: 250}
	results := tree.Query(query)
	
	if len(results) != 1 {
		t.Errorf("Expected 1 result, got %d", len(results))
	}
	if len(results) > 0 && results[0].ChromStart != 100 {
		t.Errorf("Expected interval starting at 100, got %d", results[0].ChromStart)
	}
}

func TestIntervalTreeMultipleOverlaps(t *testing.T) {
	intervals := []*bed.Record{
		{Chrom: "chr1", ChromStart: 100, ChromEnd: 300},
		{Chrom: "chr1", ChromStart: 200, ChromEnd: 400},
		{Chrom: "chr1", ChromStart: 350, ChromEnd: 500},
	}
	
	tree := NewIntervalTree(intervals)
	
	// Query overlapping with first two intervals
	query := &bed.Record{Chrom: "chr1", ChromStart: 250, ChromEnd: 350}
	results := tree.Query(query)
	
	if len(results) != 2 {
		t.Errorf("Expected 2 results, got %d", len(results))
	}
}

func TestIntervalTreeNoOverlap(t *testing.T) {
	intervals := []*bed.Record{
		{Chrom: "chr1", ChromStart: 100, ChromEnd: 200},
		{Chrom: "chr1", ChromStart: 300, ChromEnd: 400},
	}
	
	tree := NewIntervalTree(intervals)
	
	// Query with no overlaps
	query := &bed.Record{Chrom: "chr1", ChromStart: 210, ChromEnd: 290}
	results := tree.Query(query)
	
	if len(results) != 0 {
		t.Errorf("Expected 0 results, got %d", len(results))
	}
}

func TestIntervalTreeEmpty(t *testing.T) {
	tree := NewIntervalTree(nil)
	
	query := &bed.Record{Chrom: "chr1", ChromStart: 100, ChromEnd: 200}
	results := tree.Query(query)
	
	if len(results) != 0 {
		t.Errorf("Expected 0 results from empty tree, got %d", len(results))
	}
}

func TestIntervalTreeAllOverlap(t *testing.T) {
	intervals := []*bed.Record{
		{Chrom: "chr1", ChromStart: 100, ChromEnd: 200},
		{Chrom: "chr1", ChromStart: 150, ChromEnd: 250},
		{Chrom: "chr1", ChromStart: 180, ChromEnd: 280},
		{Chrom: "chr1", ChromStart: 220, ChromEnd: 320},
	}
	
	tree := NewIntervalTree(intervals)
	
	// Query overlapping with all intervals (use 221 to overlap with last interval)
	query := &bed.Record{Chrom: "chr1", ChromStart: 180, ChromEnd: 221}
	results := tree.Query(query)
	
	if len(results) != 4 {
		t.Errorf("Expected 4 results, got %d", len(results))
	}
}
