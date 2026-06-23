package fastp

import (
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fastq"
)

// TestTrimByOverlapAnalysis is a regression test for the PE read-2 adapter bug
// the GIAB real-data test surfaced: in paired-end overlap mode our pipeline
// trimmed only read-1's read-through adapter (both mates used the read-1 adapter
// sequence), leaving the read-2 adapter in place. The fix ports upstream's
// AdapterTrimmer::trimByOverlapAnalysis, which trims BOTH mates to the insert
// length when the overlap shows read-through (Offset < 0).
func TestTrimByOverlapAnalysis(t *testing.T) {
	r1 := &fastq.Record{
		Sequence: []byte(strings.Repeat("A", 40) + strings.Repeat("C", 20)),
		Quality:  []byte(strings.Repeat("I", 60)),
	}
	r2 := &fastq.Record{
		Sequence: []byte(strings.Repeat("G", 40) + strings.Repeat("T", 20)),
		Quality:  []byte(strings.Repeat("I", 60)),
	}
	ov := OverlapAnalysisResult{Overlapped: true, Offset: -20, OverlapLen: 40}

	a1, a2, did := trimByOverlapAnalysis(r1, r2, ov, 0, 0)
	if !did {
		t.Fatal("expected a trim for read-through (Offset < 0)")
	}
	if len(r1.Sequence) != 40 || len(r2.Sequence) != 40 {
		t.Fatalf("trimmed lengths R1=%d R2=%d, want 40/40", len(r1.Sequence), len(r2.Sequence))
	}
	if len(r1.Quality) != 40 || len(r2.Quality) != 40 {
		t.Fatalf("quality not resized to 40: R1=%d R2=%d", len(r1.Quality), len(r2.Quality))
	}
	if a1 != strings.Repeat("C", 20) || a2 != strings.Repeat("T", 20) {
		t.Fatalf("trimmed adapters a1=%q a2=%q", a1, a2)
	}
}

// TestTrimByOverlapAnalysisNoReadThrough verifies no trim happens when the insert
// is at least the read length (Offset >= 0) or there is no overlap.
func TestTrimByOverlapAnalysisNoReadThrough(t *testing.T) {
	mk := func() (*fastq.Record, *fastq.Record) {
		return &fastq.Record{Sequence: []byte(strings.Repeat("A", 60)), Quality: []byte(strings.Repeat("I", 60))},
			&fastq.Record{Sequence: []byte(strings.Repeat("G", 60)), Quality: []byte(strings.Repeat("I", 60))}
	}
	r1, r2 := mk()
	if _, _, did := trimByOverlapAnalysis(r1, r2, OverlapAnalysisResult{Overlapped: true, Offset: 5, OverlapLen: 55}, 0, 0); did || len(r1.Sequence) != 60 {
		t.Fatalf("must not trim when Offset >= 0 (did=%v len=%d)", did, len(r1.Sequence))
	}
	r1, r2 = mk()
	if _, _, did := trimByOverlapAnalysis(r1, r2, OverlapAnalysisResult{Overlapped: false}, 0, 0); did || len(r2.Sequence) != 60 {
		t.Fatalf("must not trim when not overlapped (did=%v len=%d)", did, len(r2.Sequence))
	}
}
