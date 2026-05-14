package samtools

import (
	"testing"
)

func TestParseRegion(t *testing.T) {
	cases := []struct {
		in   string
		want Region
		bad  bool
	}{
		{"chr1", Region{Chrom: "chr1", Beg: 1, End: 0}, false},
		{"chr1:100", Region{Chrom: "chr1", Beg: 100, End: 0}, false},
		{"chr1:100-200", Region{Chrom: "chr1", Beg: 100, End: 200}, false},
		{"chr1:100-", Region{Chrom: "chr1", Beg: 100, End: 0}, false},
		{"chr1:", Region{Chrom: "chr1", Beg: 1, End: 0}, false},
		{"chrX:1,000-2,000", Region{Chrom: "chrX", Beg: 1000, End: 2000}, false},
		{"", Region{}, true},
		{":100-200", Region{}, true},
		{"chr1:200-100", Region{}, true},
		{"chr1:abc-200", Region{}, true},
	}
	for _, c := range cases {
		got, err := ParseRegion(c.in)
		if c.bad {
			if err == nil {
				t.Errorf("ParseRegion(%q): expected error, got %+v", c.in, got)
			}
			continue
		}
		if err != nil {
			t.Errorf("ParseRegion(%q): %v", c.in, err)
			continue
		}
		if got != c.want {
			t.Errorf("ParseRegion(%q): got %+v, want %+v", c.in, got, c.want)
		}
	}
}

func TestResolveRegions(t *testing.T) {
	lookup := func(name string) int {
		switch name {
		case "chr1":
			return 0
		case "chr2":
			return 1
		}
		return -1
	}
	resolved, unknown, err := ResolveRegions([]string{"chr1:100-200", "chrZ:1-2", "chr2:50"}, lookup)
	if err != nil {
		t.Fatalf("ResolveRegions: %v", err)
	}
	if len(resolved) != 2 {
		t.Errorf("resolved count: got %d, want 2", len(resolved))
	}
	if len(unknown) != 1 || unknown[0] != "chrZ" {
		t.Errorf("unknown: got %v, want [chrZ]", unknown)
	}
	// First resolved: chr1:100-200 → 0-based [99, 200).
	if resolved[0].RefID != 0 || resolved[0].Beg0 != 99 || resolved[0].End0 != 200 {
		t.Errorf("resolved[0]: got %+v", resolved[0])
	}
	// Second resolved: chr2:50 → open-ended, End0 becomes the "infinity" sentinel.
	if resolved[1].RefID != 1 || resolved[1].Beg0 != 49 {
		t.Errorf("resolved[1]: got %+v", resolved[1])
	}
	if resolved[1].End0 < 1<<29 {
		t.Errorf("open-ended End0 should be very large, got %d", resolved[1].End0)
	}
}

func TestResolveRegionsBadSpec(t *testing.T) {
	if _, _, err := ResolveRegions([]string{""}, func(string) int { return 0 }); err == nil {
		t.Error("expected error on bad region spec")
	}
}

func TestUnionChunks(t *testing.T) {
	b := NewBAIBuilder(2)
	// chr1 has two records each in its own bin.
	b.AddRecord(0, 100, 200, 0x100, 0x200, true)
	b.AddRecord(0, 300, 400, 0x300, 0x400, true)
	// chr2 has one record.
	b.AddRecord(1, 50, 100, 0x500, 0x600, true)
	idx := b.Finish()

	resolved := []ResolvedRegion{
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

// TestRegionOverlapsRef exercises the OverlapsRef helper directly.
func TestRegionOverlapsRef(t *testing.T) {
	r := Region{Chrom: "chr1", Beg: 100, End: 200}
	// Refs in test: 0 = chr1, 1 = chr2.
	if !r.OverlapsRef(0, 0, 50, 75) {
		t.Error("record 50..125 should overlap region 100..200")
	}
	if r.OverlapsRef(0, 0, 250, 10) {
		t.Error("record 250..260 should NOT overlap region 100..200")
	}
	if r.OverlapsRef(0, 1, 100, 50) {
		t.Error("record on chr2 should NOT overlap a chr1 region")
	}
	// Open-ended region (End == 0).
	rOpen := Region{Chrom: "chr1", Beg: 100, End: 0}
	if !rOpen.OverlapsRef(0, 0, 1_000_000, 50) {
		t.Error("open-ended region should overlap far-downstream record")
	}
}
