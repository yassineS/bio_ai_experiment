package baq

import (
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// makeRec builds a minimal mapped Record for SamCapMapq tests: a single-CIGAR
// alignment at the given 1-based pos with uniform base quality q.
func makeRec(pos int32, cigar, seq string, q byte) *sam.Record {
	c, err := sam.ParseCigar(cigar)
	if err != nil {
		panic(err)
	}
	qual := make([]byte, len(seq))
	for i := range qual {
		qual[i] = q
	}
	return &sam.Record{
		QName: "r", RName: "c", Pos: pos, MapQ: 60,
		Cigar: c, Seq: seq, Qual: qual,
	}
}

// TestSamCapMapq checks SamCapMapq against values computed from htslib's
// sam_cap_mapq formula. The mismatch / quality / length statistics that drive
// the formula are reproduced exactly, so these expected caps are byte-faithful.
func TestSamCapMapq(t *testing.T) {
	// 10M read, every base quality 13. The reference matches the read
	// except where noted; each differing base is one mismatch contributing
	// q += 13. Note that upstream's `len` accumulator counts each matched
	// base inside the inner loop AND adds the whole op length again after
	// it, so for a 10M op the effective len in the formula is 20. The
	// expected caps below are taken directly from htslib's sam_cap_mapq.
	cases := []struct {
		name    string
		seq     string
		ref     string
		thres   int
		wantCap int
	}{
		{"one_mismatch_thres40", "AAAAAAAAAA", "CAAAAAAAAA", 40, 40},
		{"two_mismatch_thres40", "AAAAAAAAAA", "CCAAAAAAAA", 40, 38},
		{"three_mismatch_thres40", "AAAAAAAAAA", "CCCAAAAAAA", 40, 36},
		{"one_mismatch_thres20", "AAAAAAAAAA", "CAAAAAAAAA", 20, 20},
		// A perfectly-matching read has no mismatches: t = 0, which is
		// <= thres, so cap = round(sqrt(1)*40) = 40.
		{"perfect_match", "AAAAAAAAAA", "AAAAAAAAAA", 40, 40},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			rec := makeRec(1, "10M", tc.seq, 13)
			got := SamCapMapq(rec, []byte(tc.ref), tc.thres)
			if got != tc.wantCap {
				t.Errorf("cap = %d, want %d", got, tc.wantCap)
			}
		})
	}
}

// TestSamCapMapqLowQualIgnored checks that bases with quality below 13 are
// excluded from the mismatch counting, matching upstream's `qual[z] >= 13`
// guard: a read of all-low-quality bases has no countable mismatches.
func TestSamCapMapqLowQualIgnored(t *testing.T) {
	rec := makeRec(1, "10M", "AAAAAAAAAA", 5) // qual 5 < 13
	// No bases are counted, so mm = len = q = 0 and t = 0: cap == thres.
	if got := SamCapMapq(rec, []byte("CCCCCCCCCC"), 40); got != 40 {
		t.Errorf("all-low-quality cap = %d, want 40", got)
	}
}
