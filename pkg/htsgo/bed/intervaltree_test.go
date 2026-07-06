package bed

import "testing"

// makeRecs is a tiny helper for building per-chrom record slices in tests.
func makeRecs(chrom string, ranges [][2]int) []*Record {
	out := make([]*Record, len(ranges))
	for i, r := range ranges {
		out[i] = &Record{Chrom: chrom, ChromStart: r[0], ChromEnd: r[1]}
	}
	return out
}

func TestIntervalTree_EmptyTree(t *testing.T) {
	tree := NewIntervalTree(nil)
	if tree.Root != nil {
		t.Fatalf("expected nil Root for empty tree")
	}
	got := tree.Query(&Record{ChromStart: 0, ChromEnd: 100})
	if got != nil {
		t.Errorf("Query on empty tree: expected nil, got %v", got)
	}
}

func TestIntervalTree_SingleRecord(t *testing.T) {
	recs := makeRecs("chr1", [][2]int{{100, 200}})
	tree := NewIntervalTree(recs)

	// Overlapping query.
	got := tree.Query(&Record{ChromStart: 150, ChromEnd: 250})
	if len(got) != 1 {
		t.Fatalf("expected 1 hit, got %d", len(got))
	}

	// Non-overlapping query (before).
	got = tree.Query(&Record{ChromStart: 0, ChromEnd: 50})
	if len(got) != 0 {
		t.Errorf("expected 0 hits before, got %d", len(got))
	}

	// Non-overlapping query (after).
	got = tree.Query(&Record{ChromStart: 300, ChromEnd: 400})
	if len(got) != 0 {
		t.Errorf("expected 0 hits after, got %d", len(got))
	}

	// Touching (half-open: start==end means no overlap).
	got = tree.Query(&Record{ChromStart: 200, ChromEnd: 300})
	if len(got) != 0 {
		t.Errorf("touching half-open intervals should not overlap, got %d", len(got))
	}
}

func TestIntervalTree_MultipleOverlaps(t *testing.T) {
	recs := makeRecs("chr1", [][2]int{
		{0, 100},
		{50, 150},
		{120, 180},
		{200, 300},
		{250, 280},
	})
	tree := NewIntervalTree(recs)

	got := tree.Query(&Record{ChromStart: 60, ChromEnd: 130})
	// Should match [0,100), [50,150), [120,180).
	if len(got) != 3 {
		t.Fatalf("expected 3 hits, got %d", len(got))
	}

	got = tree.Query(&Record{ChromStart: 270, ChromEnd: 290})
	// Should match [200,300), [250,280).
	if len(got) != 2 {
		t.Fatalf("expected 2 hits, got %d", len(got))
	}

	got = tree.Query(&Record{ChromStart: 180, ChromEnd: 200})
	if len(got) != 0 {
		t.Errorf("gap query should hit nothing, got %d", len(got))
	}
}

func TestIntervalTree_QueryBeforeMaxPruning(t *testing.T) {
	// Build a tree where the leftmost record has the largest end. Queries
	// before that record should still prune correctly because the subtree's
	// Max is propagated up. This exercises the `query.ChromStart >= node.Max`
	// pruning branch.
	recs := makeRecs("chr1", [][2]int{
		{0, 1000},
		{10, 20},
		{30, 40},
	})
	tree := NewIntervalTree(recs)
	got := tree.Query(&Record{ChromStart: 5000, ChromEnd: 6000})
	if len(got) != 0 {
		t.Errorf("far-away query should hit nothing, got %d", len(got))
	}
}

func TestIntervalTree_DeepTree(t *testing.T) {
	// 64 sequential 10bp records ⇒ tree of depth ~7. Each tenth record
	// overlaps the query.
	var ranges [][2]int
	for i := 0; i < 64; i++ {
		ranges = append(ranges, [2]int{i * 10, i*10 + 10})
	}
	tree := NewIntervalTree(makeRecs("chr1", ranges))
	// Query [105, 145) overlaps records starting at 100, 110, 120, 130, 140.
	got := tree.Query(&Record{ChromStart: 105, ChromEnd: 145})
	if len(got) != 4 {
		// [100,110) overlaps 105 ⇒ hit
		// [110,120) overlaps ⇒ hit
		// [120,130) overlaps ⇒ hit
		// [130,140) overlaps ⇒ hit
		// [140,150) starts at 140 < 145 ⇒ hit (5 total)
		t.Logf("got: %v", got)
		// Recount expected.
		expect := 0
		for _, r := range ranges {
			if r[0] < 145 && r[1] > 105 {
				expect++
			}
		}
		if len(got) != expect {
			t.Fatalf("expected %d hits, got %d", expect, len(got))
		}
	}
}

// TestIntervalTree_OverlapsMatchesQuery cross-checks the allocation-free
// Overlaps predicate against (len(Query(...)) > 0) across a spread of queries,
// including the boundary (half-open touching) cases, an empty tree and a
// deep tree. Overlaps must return the same yes/no answer Query would.
func TestIntervalTree_OverlapsMatchesQuery(t *testing.T) {
	empty := NewIntervalTree(nil)
	if empty.Overlaps(0, 100) {
		t.Errorf("Overlaps on empty tree should be false")
	}

	recs := makeRecs("chr1", [][2]int{
		{0, 100}, {50, 150}, {120, 180}, {200, 300}, {250, 280}, {1000, 1000},
	})
	tree := NewIntervalTree(recs)

	// A grid of query ranges, deliberately including touching boundaries.
	queries := [][2]int{
		{-10, 0}, {0, 1}, {60, 130}, {100, 120}, {180, 200}, {200, 200},
		{270, 290}, {299, 301}, {300, 400}, {999, 1001}, {5000, 6000},
	}
	for _, q := range queries {
		want := len(tree.Query(&Record{ChromStart: q[0], ChromEnd: q[1]})) > 0
		got := tree.Overlaps(q[0], q[1])
		if got != want {
			t.Errorf("Overlaps(%d,%d)=%v; Query says %v", q[0], q[1], got, want)
		}
	}

	// Deep tree: 64 sequential 10bp records, sweep every 5bp query.
	var ranges [][2]int
	for i := 0; i < 64; i++ {
		ranges = append(ranges, [2]int{i * 10, i*10 + 10})
	}
	deep := NewIntervalTree(makeRecs("chr1", ranges))
	for s := -20; s < 700; s += 5 {
		q := [2]int{s, s + 7}
		want := len(deep.Query(&Record{ChromStart: q[0], ChromEnd: q[1]})) > 0
		if got := deep.Overlaps(q[0], q[1]); got != want {
			t.Errorf("deep Overlaps(%d,%d)=%v; Query says %v", q[0], q[1], got, want)
		}
	}
}

func TestIntervalsOverlap(t *testing.T) {
	cases := []struct {
		a, b [2]int
		want bool
	}{
		{[2]int{0, 100}, [2]int{50, 150}, true},
		{[2]int{0, 100}, [2]int{100, 200}, false}, // half-open touching
		{[2]int{0, 100}, [2]int{99, 101}, true},
		{[2]int{50, 60}, [2]int{0, 100}, true}, // contained
		{[2]int{0, 100}, [2]int{200, 300}, false},
	}
	for _, c := range cases {
		a := &Record{ChromStart: c.a[0], ChromEnd: c.a[1]}
		b := &Record{ChromStart: c.b[0], ChromEnd: c.b[1]}
		if got := intervalsOverlap(a, b); got != c.want {
			t.Errorf("intervalsOverlap(%v,%v) = %v; want %v", c.a, c.b, got, c.want)
		}
	}
}
