package bedmerge

import (
	"bytes"
	"testing"
)

// TestUnitStrandQueueOrdering pins the per-strand min-heap behaviour the merge
// state machine relies on: top()/pop() return the global (chrom,start,end)
// minimum across strands, while top(strand)/pop(strand) restrict to one strand.
func TestUnitStrandQueueOrdering(t *testing.T) {
	q := &strandQueue{}
	// Push out of order; strands interleaved.
	q.push(record{chrom: "chr1", start: 10, end: 75, strand: "-"})
	q.push(record{chrom: "chr1", start: 10, end: 60, strand: "-"})
	q.push(record{chrom: "chr1", start: 5, end: 9, strand: "+"})
	q.push(record{chrom: "chr1", start: 10, end: 50, strand: "+"})

	// Global minimum is the +strand record at start=5.
	top, ok := q.top()
	if !ok || top.start != 5 || top.strand != "+" {
		t.Fatalf("global top = %+v, ok=%v; want start 5 '+'", top, ok)
	}
	q.pop()

	// Restricted to '-': the two minus records come out end-ascending (60, 75).
	m1, ok1 := q.topStrand("-")
	if !ok1 || m1.end != 60 {
		t.Fatalf("minus top = %+v; want end 60", m1)
	}
	q.popStrand("-")
	m2, _ := q.topStrand("-")
	if m2.end != 75 {
		t.Fatalf("minus second = %+v; want end 75", m2)
	}
}

// TestUnitSortRecordsTieBreakPreservesInput pins that sortRecords keeps input
// order on equal (chrom, start) keys (no chromEnd tie-break), matching the
// stream order upstream bedtools merge consumes.
func TestUnitSortRecordsTieBreakPreservesInput(t *testing.T) {
	recs := []record{
		{chrom: "chr1", start: 10, end: 100, fields: []string{"chr1", "10", "100", "a"}},
		{chrom: "chr1", start: 10, end: 50, fields: []string{"chr1", "10", "50", "b"}},
		{chrom: "chr1", start: 10, end: 75, fields: []string{"chr1", "10", "75", "c"}},
		{chrom: "chr1", start: 5, end: 9, fields: []string{"chr1", "5", "9", "z"}},
	}
	sortRecords(recs)
	// start=5 first, then the start=10 group in INPUT order a, b, c.
	want := []string{"z", "a", "b", "c"}
	for i, w := range want {
		got := recs[i].fields[3]
		if got != w {
			t.Fatalf("position %d = %q, want %q", i, got, w)
		}
	}
}

// TestUnitMergeCollapseInputOrder is a binary-free check that the default merge
// path collapses an equal-start group's values in input order.
func TestUnitMergeCollapseInputOrder(t *testing.T) {
	in := "chr1\t10\t100\ta\nchr1\t10\t50\tb\nchr1\t10\t75\tc\nchr1\t10\t60\td\n"
	var out bytes.Buffer
	co, err := ParseColumnOps("4", "collapse")
	if err != nil {
		t.Fatalf("ParseColumnOps: %v", err)
	}
	if _, err := Merge(bytes.NewReader([]byte(in)), &out, MergeOptions{ColumnOps: co}); err != nil {
		t.Fatalf("Merge: %v", err)
	}
	want := "chr1\t10\t100\ta,b,c,d\n"
	if out.String() != want {
		t.Fatalf("collapse order.\nwant: %q\ngot:  %q", want, out.String())
	}
}
