package fastp

import (
	"bytes"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fastq"
)

func TestDupTrackerEmpty(t *testing.T) {
	tr := NewDupTracker(3)
	if tr.Rate() != 0 {
		t.Fatalf("empty tracker should have zero rate, got %v", tr.Rate())
	}
	if tr.Observe(nil) {
		t.Fatalf("nil sequence should not be flagged duplicate")
	}
	if tr.Total() != 0 {
		t.Fatalf("nil sequence should not bump total, got %d", tr.Total())
	}
	if hist := tr.Histogram(); len(hist) != 0 {
		t.Fatalf("empty tracker should have empty histogram, got %v", hist)
	}
}

func TestDupTrackerAccuracyClamping(t *testing.T) {
	low := NewDupTracker(-5)
	if low.Accuracy() != dupAccuracyMin {
		t.Fatalf("accuracy floor: want %d, got %d", dupAccuracyMin, low.Accuracy())
	}
	high := NewDupTracker(100)
	if high.Accuracy() != dupAccuracyMax {
		t.Fatalf("accuracy ceiling: want %d, got %d", dupAccuracyMax, high.Accuracy())
	}
}

func TestDupTrackerCountsDuplicates(t *testing.T) {
	tr := NewDupTracker(3)
	// Three reads, only one unique sequence -> two duplicates.
	seq := []byte("ACGTACGTACGTACGT")
	if dup := tr.Observe(seq); dup {
		t.Fatalf("first observation should not be flagged duplicate")
	}
	if dup := tr.Observe(seq); !dup {
		t.Fatalf("second observation of same key should be duplicate")
	}
	if dup := tr.Observe(seq); !dup {
		t.Fatalf("third observation of same key should be duplicate")
	}
	// One distinct read; should not be flagged.
	if dup := tr.Observe([]byte("TTTTGGGGCCCCAAAA")); dup {
		t.Fatalf("distinct sequence should not collide for small N")
	}
	if got, want := tr.Total(), int64(4); got != want {
		t.Fatalf("total: want %d, got %d", want, got)
	}
	// 2 out of 4 observations were duplicates.
	if rate := tr.Rate(); rate < 0.49 || rate > 0.51 {
		t.Fatalf("dup rate: want ~0.5, got %v", rate)
	}
}

func TestDupTrackerHistogramShape(t *testing.T) {
	tr := NewDupTracker(3)
	for i := 0; i < 5; i++ {
		tr.Observe([]byte("AAAAAAAAAAAAAAAA"))
	}
	tr.Observe([]byte("CCCCCCCCCCCCCCCC"))
	tr.Observe([]byte("GGGGGGGGGGGGGGGG"))
	hist := tr.Histogram()
	// One cell has count 5 (so contributes 5 reads to hist[5]); two
	// cells have count 1 each (contributing 1 read each to hist[1]).
	if got := hist[5]; got != 5 {
		t.Fatalf("hist[5]: want 5 reads, got %d", got)
	}
	if got := hist[1]; got != 2 {
		t.Fatalf("hist[1]: want 2 reads, got %d", got)
	}
}

func TestDupTrackerShortSequence(t *testing.T) {
	tr := NewDupTracker(2)
	tr.Observe([]byte("ACG"))
	if tr.Observe([]byte("ACG")) != true {
		t.Fatalf("short identical sequences should still be detected as duplicates")
	}
}

func TestProcessSingleEndDedup(t *testing.T) {
	// Three identical reads + one distinct -> after dedup we keep the
	// first occurrence and the distinct read; the two repeats are dropped.
	input := strings.Join([]string{
		"@r1", "ACGTACGTACGTACGTACGT", "+", "IIIIIIIIIIIIIIIIIIII",
		"@r2", "ACGTACGTACGTACGTACGT", "+", "IIIIIIIIIIIIIIIIIIII",
		"@r3", "ACGTACGTACGTACGTACGT", "+", "IIIIIIIIIIIIIIIIIIII",
		"@r4", "TTTTGGGGCCCCAAAATTTT", "+", "IIIIIIIIIIIIIIIIIIII",
		"",
	}, "\n")
	var out bytes.Buffer
	opts := DefaultProcessOptions()
	opts.MinLength = 1
	opts.LengthRequired = 1
	opts.QualThreshold = 0
	opts.DupCalcAccuracy = 3
	opts.Dedup = true

	stats, err := ProcessSingleEnd(strings.NewReader(input), &out, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("ProcessSingleEnd: %v", err)
	}
	if stats.TotalReads != 4 {
		t.Fatalf("total reads: want 4, got %d", stats.TotalReads)
	}
	if stats.DedupDropped != 2 {
		t.Fatalf("dedup dropped: want 2, got %d", stats.DedupDropped)
	}
	if stats.CleanReads != 2 {
		t.Fatalf("clean reads: want 2, got %d", stats.CleanReads)
	}
	if stats.DupRate <= 0 {
		t.Fatalf("dup rate should be > 0, got %v", stats.DupRate)
	}
	if stats.DupHist == nil {
		t.Fatalf("dup hist should be populated when accuracy > 0")
	}
}

func TestProcessSingleEndDupCalcOnly(t *testing.T) {
	// Same setup but without --dedup: stats should still report rate
	// and histogram, but no reads are dropped.
	input := strings.Join([]string{
		"@r1", "ACGTACGTACGTACGTACGT", "+", "IIIIIIIIIIIIIIIIIIII",
		"@r2", "ACGTACGTACGTACGTACGT", "+", "IIIIIIIIIIIIIIIIIIII",
		"@r3", "TTTTGGGGCCCCAAAATTTT", "+", "IIIIIIIIIIIIIIIIIIII",
		"",
	}, "\n")
	var out bytes.Buffer
	opts := DefaultProcessOptions()
	opts.MinLength = 1
	opts.LengthRequired = 1
	opts.QualThreshold = 0
	opts.DupCalcAccuracy = 2

	stats, err := ProcessSingleEnd(strings.NewReader(input), &out, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("ProcessSingleEnd: %v", err)
	}
	if stats.CleanReads != 3 {
		t.Fatalf("clean reads (no dedup): want 3, got %d", stats.CleanReads)
	}
	if stats.DedupDropped != 0 {
		t.Fatalf("dedup dropped (no dedup): want 0, got %d", stats.DedupDropped)
	}
	if stats.DupTotal != 3 {
		t.Fatalf("dup total: want 3, got %d", stats.DupTotal)
	}
}
