package fastp

import "testing"

func TestOverrepSteps(t *testing.T) {
	got := overrepSteps(100)
	want := [5]int{10, 20, 40, 100, 98}
	if got != want {
		t.Fatalf("overrepSteps(100) = %v, want %v", got, want)
	}
	// seqLen-2 caps the final step below 150.
	if s := overrepSteps(30); s[4] != 28 {
		t.Fatalf("overrepSteps(30) last step = %d, want 28", s[4])
	}
	// seqLen large: final step capped at 150.
	if s := overrepSteps(300); s[4] != 150 {
		t.Fatalf("overrepSteps(300) last step = %d, want 150", s[4])
	}
}

func TestOverRepPassed(t *testing.T) {
	// s * count must exceed the length-specific bar.
	cases := []struct {
		seqLen   int
		count    int64
		sampling int
		want     bool
	}{
		{10, 51, 10, true},  // 510 > 500
		{10, 50, 10, false}, // 500 not > 500
		{20, 21, 10, true},  // 210 > 200
		{40, 11, 10, true},  // 110 > 100
		{100, 6, 10, true},  // 60 > 50
		{50, 3, 10, true},   // default bar 20: 30 > 20
		{50, 2, 10, false},  // 20 not > 20
	}
	for _, c := range cases {
		seq := rep('A', c.seqLen)
		if got := overRepPassed(seq, c.count, c.sampling); got != c.want {
			t.Errorf("overRepPassed(len=%d,count=%d,s=%d) = %v, want %v",
				c.seqLen, c.count, c.sampling, got, c.want)
		}
	}
}

func TestBuildHotSeqsAndSampling(t *testing.T) {
	// Build an input where a 40bp motif recurs frequently enough to clear
	// the length-40 threshold (>= 20 occurrences in the build pass).
	motif := rep('A', 20) + rep('C', 20) // 40bp, high "complexity" for the map
	var seqs []string
	for i := 0; i < 40; i++ {
		seqs = append(seqs, motif+rep('G', 60)) // 100bp reads
	}
	a := newOverrepAnalyzer(seqs, 1, 100)
	if len(a.counts) == 0 {
		t.Fatalf("expected at least one candidate hot sequence")
	}
	// Sample the same reads; the motif (and its windows) should accumulate.
	for _, s := range seqs {
		a.sampleRead(s)
	}
	passed := a.passedSequences()
	if len(passed) == 0 {
		t.Fatalf("expected at least one passing overrepresented sequence")
	}
}

func TestSortedSeqs(t *testing.T) {
	m := map[string]int64{"CCCC": 1, "AAAA": 2, "BBBB": 3}
	got := sortedSeqs(m)
	want := []string{"AAAA", "BBBB", "CCCC"}
	for i := range want {
		if got[i] != want[i] {
			t.Fatalf("sortedSeqs = %v, want %v", got, want)
		}
	}
}
