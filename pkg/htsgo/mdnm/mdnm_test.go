package mdnm

import (
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// mustCigar parses a CIGAR string into a sam.Cigar, failing the test on a
// malformed string. It keeps the table cases readable.
func mustCigar(t *testing.T, s string) sam.Cigar {
	t.Helper()
	c, err := sam.ParseCigar(s)
	if err != nil {
		t.Fatalf("ParseCigar(%q): %v", s, err)
	}
	return c
}

// TestCompute exercises the MD/NM walk over a fixed reference for a spread of
// CIGAR shapes: perfect matches, substitutions, insertions, deletions,
// reference skips, soft/hard clips, N bases, adjacent mismatches, and a
// slice-local reference span (refOffset > 0).
func TestCompute(t *testing.T) {
	// 0-based:    0         1         2
	//             0123456789012345678901234567
	const ref = "ACGTACGTACGTACGTACGTACGTACGT"

	cases := []struct {
		name      string
		pos       int // 1-based POS
		cigar     string
		seq       string
		refOffset int
		wantMD    string
		wantNM    int
	}{
		{
			name:   "perfect match",
			pos:    1,
			cigar:  "8M",
			seq:    "ACGTACGT",
			wantMD: "8",
			wantNM: 0,
		},
		{
			name:   "single substitution",
			pos:    1,
			cigar:  "8M",
			seq:    "ACGTACGA", // last T->A, ref base T at index 7
			wantMD: "7T0",
			wantNM: 1,
		},
		{
			name:   "single internal mismatch",
			pos:    1,
			cigar:  "4M",
			seq:    "AGGT", // idx1 ref C vs read G mismatch; others match
			wantMD: "1C2",
			wantNM: 1,
		},
		{
			name:   "leading and adjacent mismatches",
			pos:    1,
			cigar:  "4M",
			seq:    "TTGT", // idx0 A->T, idx1 C->T, idx2 G==G, idx3 T==T
			wantMD: "0A0C2",
			wantNM: 2,
		},
		{
			name:   "insertion adds to NM only",
			pos:    1,
			cigar:  "4M2I4M",
			seq:    "ACGTNNACGT", // 4M ref ACGT, 2I (NN, not counted vs ref), 4M ref ACGT (pos5..8)
			wantMD: "8",
			wantNM: 2,
		},
		{
			name:   "deletion emits caret bases",
			pos:    1,
			cigar:  "4M2D4M",
			seq:    "ACGTGTAC", // 4M=ACGT, delete ref idx4..5 (AC), 4M ref idx6..9 (GTAC)
			wantMD: "4^AC4",
			wantNM: 2,
		},
		{
			name:   "reference skip does not affect MD/NM",
			pos:    1,
			cigar:  "4M4N4M",
			seq:    "ACGTACGT", // 4M ref idx0..3, skip idx4..7, 4M ref idx8..11 (ACGT)
			wantMD: "8",
			wantNM: 0,
		},
		{
			name:   "soft clip consumes query only",
			pos:    1,
			cigar:  "2S6M",
			seq:    "XXACGTAC", // 2S skipped, 6M ref idx0..5 (ACGTAC)
			wantMD: "6",
			wantNM: 0,
		},
		{
			name:   "hard clip consumes nothing",
			pos:    1,
			cigar:  "2H6M",
			seq:    "ACGTAC",
			wantMD: "6",
			wantNM: 0,
		},
		{
			name:   "N in read is a mismatch",
			pos:    1,
			cigar:  "4M",
			seq:    "ANGT", // idx1 ref C vs read N -> mismatch
			wantMD: "1C2",
			wantNM: 1,
		},
		{
			name:      "slice-local reference span via refOffset",
			pos:       9, // ref idx 8 = 'A'
			cigar:     "4M",
			seq:       "ACGA", // idx8..11 ref = ACGT; last T->A
			refOffset: 8,      // ref[] passed below starts at 0-based 8
			wantMD:    "3T0",
			wantNM:    1,
		},
		{
			name:   "deletion then mismatch emits 0 after caret",
			pos:    1,
			cigar:  "4M2D4M",
			seq:    "ACGTATAC", // after del, 4M ref idx6..9 = GTAC; read ATAC: idx6 G->A mismatch
			wantMD: "4^AC0G3",
			wantNM: 3,
		},
	}

	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			rec := &sam.Record{
				Pos:   int64(tc.pos),
				Cigar: mustCigar(t, tc.cigar),
				Seq:   tc.seq,
			}
			r := []byte(ref)
			if tc.refOffset > 0 {
				r = r[tc.refOffset:]
			}
			md, nm := Compute(rec, r, tc.refOffset)
			if md != tc.wantMD {
				t.Errorf("MD = %q, want %q", md, tc.wantMD)
			}
			if nm != tc.wantNM {
				t.Errorf("NM = %d, want %d", nm, tc.wantNM)
			}
		})
	}
}
