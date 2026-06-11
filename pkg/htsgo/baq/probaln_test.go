package baq

import "testing"

// residues translates an ASCII A/C/G/T string into the 0/1/2/3/4 encoding
// ProbalnGlocal expects.
func residues(s string) []byte {
	out := make([]byte, len(s))
	for i := 0; i < len(s); i++ {
		out[i] = asciiToResidue[s[i]]
	}
	return out
}

// TestProbalnGlocalVectors checks ProbalnGlocal against scores, posterior
// quality arrays and decoded states captured from htslib's own probaln.c
// (compiled standalone). These pin the numerical core byte-for-byte.
func TestProbalnGlocalVectors(t *testing.T) {
	cases := []struct {
		ref, query string
		qv         byte
		bw         int
		wantScore  int
		wantQ      []byte
		wantState  []int
	}{
		{
			ref: "acttc", query: "attc", qv: 30, bw: 10,
			wantScore: 38,
			wantQ:     []byte{4, 17, 19, 20},
			wantState: []int{4, 8, 12, 16},
		},
		{
			ref: "acttc", query: "attc", qv: 10, bw: 10,
			wantScore: 21,
			wantQ:     []byte{13, 14, 14, 14},
			wantState: []int{4, 8, 12, 16},
		},
		{
			ref: "GGGCATTAGCC", query: "GGGCATAGCC", qv: 25, bw: 7,
			wantScore: 33,
			wantQ:     []byte{32, 33, 36, 58, 33, 3, 33, 54, 36, 33},
			wantState: []int{0, 4, 8, 12, 16, 20, 28, 32, 36, 40},
		},
		{
			ref: "ACGTACGTACGT", query: "ACGTACGTACGT", qv: 40, bw: 10,
			wantScore: 5,
			wantQ:     []byte{36, 52, 67, 74, 74, 74, 74, 74, 74, 67, 52, 36},
			wantState: []int{0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44},
		},
		{
			ref: "TTTTTTTTTT", query: "TTTTTATTTT", qv: 20, bw: 5,
			wantScore: 30,
			wantQ:     []byte{12, 12, 12, 12, 12, 9, 12, 12, 12, 12},
			wantState: []int{0, 4, 8, 12, 16, 20, 24, 28, 32, 36},
		},
	}

	for i, tc := range cases {
		ref := residues(tc.ref)
		query := residues(tc.query)
		iqual := make([]byte, len(query))
		for j := range iqual {
			iqual[j] = tc.qv
		}
		state := make([]int, len(query))
		q := make([]byte, len(query))
		score, err := ProbalnGlocal(ref, query, iqual, Par{D: 0.001, E: 0.1, BW: tc.bw}, state, q)
		if err != nil {
			t.Fatalf("case %d: ProbalnGlocal error: %v", i, err)
		}
		if score != tc.wantScore {
			t.Errorf("case %d: score = %d, want %d", i, score, tc.wantScore)
		}
		for j := range tc.wantQ {
			if q[j] != tc.wantQ[j] {
				t.Errorf("case %d: q[%d] = %d, want %d", i, j, q[j], tc.wantQ[j])
			}
		}
		for j := range tc.wantState {
			if state[j] != tc.wantState[j] {
				t.Errorf("case %d: state[%d] = %d, want %d", i, j, state[j], tc.wantState[j])
			}
		}
	}
}

// TestProbalnGlocalHomopolymerFloatWidth pins the single-ULP rounding case
// that motivated matching htslib's float-width g_qual2prob table. htslib
// stores the per-base error probabilities as C `float` (32-bit) and only
// promotes them to double inside the forward/backward recurrences; keeping
// the table at full float64 precision instead perturbs the scaled DP enough
// to flip the integer Phred score by one unit at long homopolymer columns.
//
// This ref/query pair is one such case: it straddles many 'A' homopolymer
// runs with a single-base deletion, at the low base quality (q=12) where the
// emission probabilities matter most. The upstream probaln.c self-test
// (compiled standalone) returns 140 for these inputs. The pre-fix float64
// table returned 141; matching the float-width table returns 140 like
// htslib. See docs/PARITY_ROADMAP.md (mpileup probaln residual, now closed).
func TestProbalnGlocalHomopolymerFloatWidth(t *testing.T) {
	ref := residues("TCTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATACGAAAAAAAAAAAAAAAAAAAAAAAAAATCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATAAAAAAAAAAAAGCAAAAAAAAAAAAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAAAAAAAAAAAAACAAAAAGAAAAAAAAAACCAAAAAAAAAAAAAAAAAAAAAAGCTAAAAAAAAAAAAAGAAAAAAAAAAAAAAAAGTTCAAAAAAAATAAAAAAAAAAAAAAAAAAAAAAAAAAC")
	query := residues("TCTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATACGAAAAAAAAAAAAAAAAAAAAAAAAATCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATAAAAAAAAAAAAGCAAAAAAAAAAAAGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGAAAAAAAAAAAAACAAAAAGAAAAAAAAAACCAAAAAAAAAAAAAAAAAAAAAAGCTAAAAAAAAAAAAAGAAAAAAAAAAAAAAAAGTTCAAAAAAAATAAAAAAAAAAAAAAAAAAAAAAAAAAC")
	iqual := make([]byte, len(query))
	for i := range iqual {
		iqual[i] = 12
	}
	score, err := ProbalnGlocal(ref, query, iqual, Par{D: 0.001, E: 0.1, BW: 13}, nil, nil)
	if err != nil {
		t.Fatalf("ProbalnGlocal error: %v", err)
	}
	if score != 140 {
		t.Errorf("homopolymer float-width score = %d, want 140 (htslib float-table value; 141 is the pre-fix float64 drift)", score)
	}
}

// TestProbalnGlocalDegenerate exercises the early-return paths: empty ref or
// query, and the likelihood-only mode (nil state/q).
func TestProbalnGlocalDegenerate(t *testing.T) {
	if s, err := ProbalnGlocal(nil, residues("ACGT"), nil, Par{D: 0.001, E: 0.1, BW: 10}, nil, nil); err != nil || s != 0 {
		t.Errorf("empty ref: score=%d err=%v, want 0/nil", s, err)
	}
	if s, err := ProbalnGlocal(residues("ACGT"), nil, nil, Par{D: 0.001, E: 0.1, BW: 10}, nil, nil); err != nil || s != 0 {
		t.Errorf("empty query: score=%d err=%v, want 0/nil", s, err)
	}
	// Likelihood-only mode must match the score from the full-decode run.
	ref, query := residues("acttc"), residues("attc")
	iqual := []byte{30, 30, 30, 30}
	full, _ := ProbalnGlocal(ref, query, iqual, Par{D: 0.001, E: 0.1, BW: 10}, make([]int, 4), make([]byte, 4))
	lk, err := ProbalnGlocal(ref, query, iqual, Par{D: 0.001, E: 0.1, BW: 10}, nil, nil)
	if err != nil || lk != full {
		t.Errorf("likelihood-only score=%d err=%v, want %d/nil", lk, err, full)
	}
}
