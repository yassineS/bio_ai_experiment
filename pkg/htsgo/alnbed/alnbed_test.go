package alnbed

import (
	"strings"
	"testing"
)

const samFixture = `@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
@SQ	SN:chr2	LN:500
r1	0	chr1	11	60	10M	*	0	0	ACGTACGTAC	IIIIIIIIII
r2	16	chr1	21	60	5M10N5M	*	0	0	ACGTACGTAC	IIIIIIIIII
unmapped	4	*	0	0	*	*	0	0	ACGT	IIII
r3	0	chr2	1	60	3M2D4M	*	0	0	ACGTACG	IIIIIII
`

func TestReader_References(t *testing.T) {
	r, err := NewReader(strings.NewReader(samFixture))
	if err != nil {
		t.Fatalf("NewReader: %v", err)
	}
	refs := r.References()
	if len(refs) != 2 {
		t.Fatalf("want 2 refs, got %d", len(refs))
	}
	if refs[0].Name != "chr1" || refs[0].Length != 1000 {
		t.Errorf("ref0 = %+v", refs[0])
	}
	if refs[1].Name != "chr2" || refs[1].Length != 500 {
		t.Errorf("ref1 = %+v", refs[1])
	}
}

func TestReader_ConvertsAndSkipsUnmapped(t *testing.T) {
	r, err := NewReader(strings.NewReader(samFixture))
	if err != nil {
		t.Fatalf("NewReader: %v", err)
	}
	var recs []string
	for {
		rec, err := r.Read()
		if err != nil {
			break
		}
		recs = append(recs, rec.Chrom)

		switch rec.Name {
		case "r1":
			// 10M at 1-based pos 11 -> 0-based [10,20), single block, '+'.
			if rec.ChromStart != 10 || rec.ChromEnd != 20 {
				t.Errorf("r1 span = [%d,%d), want [10,20)", rec.ChromStart, rec.ChromEnd)
			}
			if rec.Strand != "+" || rec.BlockCount != 1 {
				t.Errorf("r1 strand/blocks = %s/%d", rec.Strand, rec.BlockCount)
			}
		case "r2":
			// 5M10N5M at pos 21 -> [20,25) and [35,40); whole span [20,40),
			// reverse strand.
			if rec.ChromStart != 20 || rec.ChromEnd != 40 {
				t.Errorf("r2 span = [%d,%d), want [20,40)", rec.ChromStart, rec.ChromEnd)
			}
			if rec.Strand != "-" {
				t.Errorf("r2 strand = %s, want -", rec.Strand)
			}
			if rec.BlockCount != 2 || len(rec.BlockSizes) != 2 {
				t.Fatalf("r2 blocks = %d %v", rec.BlockCount, rec.BlockSizes)
			}
			if rec.BlockSizes[0] != 5 || rec.BlockSizes[1] != 5 {
				t.Errorf("r2 block sizes = %v, want [5 5]", rec.BlockSizes)
			}
			if rec.BlockStarts[0] != 0 || rec.BlockStarts[1] != 15 {
				t.Errorf("r2 block starts = %v, want [0 15]", rec.BlockStarts)
			}
		case "r3":
			// 3M2D4M at pos 1 -> deletion does NOT split: single block [0,9).
			if rec.ChromStart != 0 || rec.ChromEnd != 9 {
				t.Errorf("r3 span = [%d,%d), want [0,9)", rec.ChromStart, rec.ChromEnd)
			}
			if rec.BlockCount != 1 {
				t.Errorf("r3 should be one block (D does not split), got %d", rec.BlockCount)
			}
		}
	}
	// r1, r2, r3 emitted; unmapped skipped.
	if len(recs) != 3 {
		t.Fatalf("want 3 mapped records, got %d (%v)", len(recs), recs)
	}
}
