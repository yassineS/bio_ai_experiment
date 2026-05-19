// Tests for the BAI-specific UnionChunks helper. These were split out
// of pkg/htsgo/region/region_test.go in PR-F because they exercise
// BAI machinery (BAIBuilder, BAIIndex) that lives in this package —
// keeping them here avoids a region→bam cyclic dependency in the
// region package's tests.

package bam

import (
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/region"
)

func TestUnionChunks(t *testing.T) {
	b := NewBAIBuilder(2)
	// chr1 has two records each in its own bin.
	b.AddRecord(0, 100, 200, 0x100, 0x200, true)
	b.AddRecord(0, 300, 400, 0x300, 0x400, true)
	// chr2 has one record.
	b.AddRecord(1, 50, 100, 0x500, 0x600, true)
	idx := b.Finish()

	resolved := []region.ResolvedRegion{
		{RefID: 0, Beg0: 50, End0: 250}, // covers the chr1 first record
		{RefID: 1, Beg0: 0, End0: 200},  // covers the chr2 record
	}
	chunks := UnionChunks(idx, resolved)
	if len(chunks) == 0 {
		t.Fatal("UnionChunks: got 0 chunks, want ≥1")
	}
	// Both source chunks should be represented.
	var sawChr1, sawChr2 bool
	for _, c := range chunks {
		if c.Beg == 0x100 {
			sawChr1 = true
		}
		if c.Beg == 0x500 {
			sawChr2 = true
		}
	}
	if !sawChr1 {
		t.Error("missing chr1 chunk in union")
	}
	if !sawChr2 {
		t.Error("missing chr2 chunk in union")
	}
}

func TestUnionChunksEmpty(t *testing.T) {
	idx := &BAIIndex{}
	if got := UnionChunks(idx, nil); got != nil {
		t.Errorf("expected nil for empty input, got %v", got)
	}
}
