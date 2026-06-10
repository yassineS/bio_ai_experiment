package edlib

import (
	"strings"
	"testing"
)

func TestAlignNWBasic(t *testing.T) {
	cases := []struct {
		name     string
		query    string
		target   string
		wantDist int
	}{
		{"identical", "ACGT", "ACGT", 0},
		{"single substitution", "ACGT", "ACAT", 1},
		{"single insertion in target", "ACGT", "ACGGT", 1},
		{"single deletion in target", "ACGT", "ACT", 1},
		{"two substitutions", "AAAA", "AGGA", 2},
		{"completely different", "AAAA", "CCCC", 4},
		{"long match", strings.Repeat("AC", 50), strings.Repeat("AC", 50), 0},
		{"long with one mismatch", strings.Repeat("AC", 50), strings.Repeat("AC", 25) + "GC" + strings.Repeat("AC", 24), 1},
		{"empty query empty target", "", "", 0},
		{"empty query nonempty target", "", "ACGT", 4},
		{"nonempty query empty target", "ACGT", "", 4},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			r, err := Align([]byte(tc.query), []byte(tc.target), Config{K: -1, Mode: ModeNW, Task: TaskDistance})
			if err != nil {
				t.Fatalf("Align: %v", err)
			}
			if r.EditDistance != tc.wantDist {
				t.Fatalf("dist = %d want %d", r.EditDistance, tc.wantDist)
			}
		})
	}
}

func TestAlignHWBasic(t *testing.T) {
	// HW = infix: query should embed in target with leading/trailing gaps free.
	cases := []struct {
		name     string
		query    string
		target   string
		wantDist int
		wantEnd  int // expected last end position (0-based, inclusive)
	}{
		{"exact substring", "ACGT", "TTACGTAA", 0, 5},
		{"substring middle", "GGG", "AAAGGGTTT", 0, 5},
		{"one mismatch inside", "GGG", "AAAGCGTTT", 1, 5},
		{"prefix only", "ACGT", "ACGTAAAA", 0, 3},
		{"suffix only", "ACGT", "AAAAACGT", 0, 7},
		{"single base", "C", "AAAACAAAA", 0, 4},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			r, err := Align([]byte(tc.query), []byte(tc.target), Config{K: -1, Mode: ModeHW, Task: TaskLoc})
			if err != nil {
				t.Fatalf("Align: %v", err)
			}
			if r.EditDistance != tc.wantDist {
				t.Fatalf("dist = %d want %d", r.EditDistance, tc.wantDist)
			}
			if len(r.EndLocations) == 0 {
				t.Fatalf("no end locations")
			}
			found := false
			for _, e := range r.EndLocations {
				if e == tc.wantEnd {
					found = true
					break
				}
			}
			if !found {
				t.Fatalf("end locations %v, want one to be %d", r.EndLocations, tc.wantEnd)
			}
		})
	}
}

func TestAlignSHWBasic(t *testing.T) {
	// SHW = prefix: trailing target is free.
	cases := []struct {
		name     string
		query    string
		target   string
		wantDist int
	}{
		{"exact prefix", "AACT", "AACTGGC", 0},
		{"prefix with one mismatch", "AAGT", "AACTGGC", 1},
		{"longer target same start", "ACGT", "ACGTNNNNNNNNNNN", 0},
		{"query == target", "ACGT", "ACGT", 0},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			r, err := Align([]byte(tc.query), []byte(tc.target), Config{K: -1, Mode: ModeSHW, Task: TaskDistance})
			if err != nil {
				t.Fatalf("Align: %v", err)
			}
			if r.EditDistance != tc.wantDist {
				t.Fatalf("dist = %d want %d", r.EditDistance, tc.wantDist)
			}
		})
	}
}

func TestAlignKCutoff(t *testing.T) {
	// With k=0, anything other than an exact match should yield -1.
	r, err := Align([]byte("ACGT"), []byte("AGGT"), Config{K: 0, Mode: ModeNW, Task: TaskDistance})
	if err != nil {
		t.Fatalf("Align: %v", err)
	}
	if r.EditDistance != -1 {
		t.Fatalf("expected -1 (over k), got %d", r.EditDistance)
	}

	// k=2 should accept distance 1.
	r, err = Align([]byte("ACGT"), []byte("AGGT"), Config{K: 2, Mode: ModeNW, Task: TaskDistance})
	if err != nil {
		t.Fatalf("Align: %v", err)
	}
	if r.EditDistance != 1 {
		t.Fatalf("expected 1, got %d", r.EditDistance)
	}
}

func TestAlignHWStartLoc(t *testing.T) {
	// "GCT" embedded in "AAAGCTTTT": should start at 3, end at 5.
	r, err := Align([]byte("GCT"), []byte("AAAGCTTTT"), Config{K: -1, Mode: ModeHW, Task: TaskLoc})
	if err != nil {
		t.Fatalf("Align: %v", err)
	}
	if r.EditDistance != 0 {
		t.Fatalf("dist = %d want 0", r.EditDistance)
	}
	if len(r.StartLocations) == 0 || len(r.EndLocations) == 0 {
		t.Fatalf("missing locations")
	}
	if r.StartLocations[0] != 3 || r.EndLocations[0] != 5 {
		t.Fatalf("start=%v end=%v, want start=[3] end=[5]", r.StartLocations, r.EndLocations)
	}
}

func TestAlignNWPath(t *testing.T) {
	// q="ACGT" vs t="AGT" => expected path: match A, mismatch (or delete) ...
	// distance is 1 — best CIGAR is 1=1D2= (match, delete C from query, match GT).
	r, err := Align([]byte("ACGT"), []byte("AGT"), Config{K: -1, Mode: ModeNW, Task: TaskPath})
	if err != nil {
		t.Fatalf("Align: %v", err)
	}
	if r.EditDistance != 1 {
		t.Fatalf("dist = %d want 1", r.EditDistance)
	}
	if r.Alignment == nil {
		t.Fatalf("alignment nil")
	}
	// Verify the alignment opcodes consume q and t correctly.
	pq, pt := 0, 0
	for _, op := range r.Alignment {
		switch op {
		case OpMatch, OpMismatch:
			pq++
			pt++
		case OpInsert: // consume target only
			pt++
		case OpDelete: // consume query only
			pq++
		}
	}
	if pq != 4 || pt != 3 {
		t.Fatalf("alignment did not consume q,t fully: pq=%d pt=%d", pq, pt)
	}
}

func TestAlignHWPath(t *testing.T) {
	// Embed "GAT" inside "AAAGCTAAA"; closest match is at position 3-5 with
	// edit distance 1 (substitute A->C).
	r, err := Align([]byte("GAT"), []byte("AAAGCTAAA"), Config{K: -1, Mode: ModeHW, Task: TaskPath})
	if err != nil {
		t.Fatalf("Align: %v", err)
	}
	if r.EditDistance != 1 {
		t.Fatalf("dist = %d want 1", r.EditDistance)
	}
	if r.Alignment == nil {
		t.Fatalf("alignment nil")
	}
	// Sum of consumed query positions must equal queryLen.
	consumedQ := 0
	for _, op := range r.Alignment {
		if op != OpInsert {
			consumedQ++
		}
	}
	if consumedQ != 3 {
		t.Fatalf("alignment consumed %d query bases, want 3", consumedQ)
	}
}

func TestAlignNWPathExact(t *testing.T) {
	r, err := Align([]byte("ACGT"), []byte("ACGT"), Config{K: -1, Mode: ModeNW, Task: TaskPath})
	if err != nil {
		t.Fatalf("Align: %v", err)
	}
	if r.EditDistance != 0 {
		t.Fatalf("dist = %d want 0", r.EditDistance)
	}
	for _, op := range r.Alignment {
		if op != OpMatch {
			t.Fatalf("expected all matches, got op %d in %v", op, r.Alignment)
		}
	}
	if len(r.Alignment) != 4 {
		t.Fatalf("alignment length %d, want 4", len(r.Alignment))
	}
}

func TestAlignLong(t *testing.T) {
	// Stress the multi-block path (qLen > 64).
	q := strings.Repeat("ACGT", 30) // 120 bases
	tgt := "NNN" + q + "NNN"
	r, err := Align([]byte(q), []byte(tgt), Config{K: -1, Mode: ModeHW, Task: TaskLoc})
	if err != nil {
		t.Fatalf("Align: %v", err)
	}
	if r.EditDistance != 0 {
		t.Fatalf("dist = %d want 0", r.EditDistance)
	}
	if r.StartLocations[0] != 3 || r.EndLocations[0] != 122 {
		t.Fatalf("start=%v end=%v, want 3 / 122", r.StartLocations, r.EndLocations)
	}
}

func TestAlignLongMismatches(t *testing.T) {
	// Substitute a single base in a long query: distance should be 1.
	q := []byte(strings.Repeat("A", 200))
	tgt := []byte(strings.Repeat("A", 100) + "G" + strings.Repeat("A", 99))
	r, err := Align(q, tgt, Config{K: -1, Mode: ModeNW, Task: TaskDistance})
	if err != nil {
		t.Fatalf("Align: %v", err)
	}
	if r.EditDistance != 1 {
		t.Fatalf("dist = %d want 1", r.EditDistance)
	}
}

func TestAlignHWEmptyQuery(t *testing.T) {
	r, err := Align(nil, []byte("ACGT"), Config{K: -1, Mode: ModeHW, Task: TaskLoc})
	if err != nil {
		t.Fatalf("Align: %v", err)
	}
	if r.EditDistance != 0 {
		t.Fatalf("dist = %d want 0 for empty query", r.EditDistance)
	}
}

func TestAlignNWEmptyInputs(t *testing.T) {
	r, err := Align(nil, nil, Config{K: -1, Mode: ModeNW, Task: TaskDistance})
	if err != nil {
		t.Fatalf("Align: %v", err)
	}
	if r.EditDistance != 0 {
		t.Fatalf("dist = %d want 0", r.EditDistance)
	}
}
