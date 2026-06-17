package fastp

// Binary-free unit tests for (1) the cut_tail sliding-window boundary that the
// parity pipeline flagged at scale, and (2) the plumbing of the new
// --disable_quality_filtering / --disable_length_filtering / --disable_trim_poly_g
// flags. These run without the upstream binary.

import (
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fastq"
)

// TestUnitCutTailBoundary pins the exact cut position cut_tail computes for
// hand-checked length+quality windows. The cut_tail walk (filter.cpp:180-209)
// scans a window of size w from the 3' end toward the 5' end; on the first
// window (ending at index t) whose mean quality is >= threshold, IF that
// window is not already the 3'-most one (t < l-1) it rewinds t = t-w+1 to the
// window's START — dropping w-1 extra bases — then skips any trailing N's and
// keeps [front, t+1). The strict loop bound `t-w >= front` makes the 5'-most
// window the floor. The expected values below are exactly what upstream
// produces (verified byte-for-byte by TestParity_Fastp_Case18); these cases
// lock the rolling sum, the t-w+1 rewind, and the trailing-N skip so a
// regression is caught without the upstream binary.
func TestUnitCutTailBoundary(t *testing.T) {
	enc := fastq.Phred33
	// Phred33 ascii: 'I'=40, '5'=20, '#'=2; 'N' is an N base.
	cases := []struct {
		name   string
		seq    string
		qual   string
		window int
		mq     int
		wantLo int
		wantHi int
	}{
		{
			// All high quality: the 3'-most window already passes, nothing is
			// cut, keep the whole read.
			name: "all_high_keep_all", window: 4, mq: 20,
			seq:    strings.Repeat("A", 20),
			qual:   strings.Repeat("I", 20),
			wantLo: 0, wantHi: 20,
		},
		{
			// High body (16 'I') then a 4-base low-q tail ('#'=2). Scanning
			// 3'->5', the first window with mean >= 20 ends at t=15 (window
			// [12..15] = I,I,I,# -> (40*3+2)/4 = 30.5). Upstream then rewinds
			// t = t-w+1 = 12 and keeps [0, t+1) = [0,15). That rewind drops
			// w-1 = 3 extra bases of the passing window — the load-bearing
			// quirk that makes the cut 5bp, not the naive 4bp.
			name: "low_tail_rewind", window: 4, mq: 20,
			seq:    strings.Repeat("A", 20),
			qual:   strings.Repeat("I", 16) + strings.Repeat("#", 4),
			wantLo: 0, wantHi: 15,
		},
		{
			// Same qualities; trailing bases are 'N'. The kept region after the
			// quality rewind ([0,15)) is all 'A', so the trailing-N skip is a
			// no-op and the boundary is identical to low_tail_rewind.
			name: "low_tail_rewind_with_N_tail", window: 4, mq: 20,
			seq:    strings.Repeat("A", 16) + strings.Repeat("N", 4),
			qual:   strings.Repeat("I", 16) + strings.Repeat("#", 4),
			wantLo: 0, wantHi: 15,
		},
		{
			// Kept region ends in N's: the trailing-N skip pulls the boundary
			// left. Body is 11 'A' + 5 'N' (still hi-Q) + 4 'A' + 4 '#' tail.
			// The quality rewind lands at index 14 (kept [0,15)); the N-skip
			// then walks left over the N's, leaving [0,11).
			name: "trailing_N_skip", window: 4, mq: 20,
			seq:    strings.Repeat("A", 11) + strings.Repeat("N", 5) + strings.Repeat("A", 4),
			qual:   strings.Repeat("I", 16) + strings.Repeat("#", 4),
			wantLo: 0, wantHi: 11,
		},
		{
			// Window mean exactly on the threshold (mean == mq) counts as PASS
			// (>=). The 3'-most window is quality '5'(=20) with mq=20, so it
			// passes immediately at the natural end with no rewind: keep all.
			name: "mean_exactly_threshold_keeps", window: 4, mq: 20,
			seq:    strings.Repeat("A", 12),
			qual:   strings.Repeat("I", 8) + strings.Repeat("5", 4),
			wantLo: 0, wantHi: 12,
		},
		{
			// Entirely low quality: no window ever reaches the threshold, the
			// loop runs to its floor and t lands at front; rlen = t-front+1 = 1,
			// so upstream keeps a single 5' base ([0,1)). The MinLength filter
			// downstream then rejects such a 1bp read.
			name: "all_low_keeps_one", window: 4, mq: 20,
			seq:    strings.Repeat("A", 12),
			qual:   strings.Repeat("#", 12),
			wantLo: 0, wantHi: 1,
		},
	}
	for _, c := range cases {
		c := c
		t.Run(c.name, func(t *testing.T) {
			opts := DefaultProcessOptions()
			opts.CutTail = true
			opts.CutWindowSize = c.window
			opts.CutMeanQuality = c.mq
			lo, hi := slidingWindowCut([]byte(c.seq), []byte(c.qual), enc, opts)
			if lo != c.wantLo || hi != c.wantHi {
				t.Fatalf("cut_tail [%d,%d), want [%d,%d) (kept %q)",
					lo, hi, c.wantLo, c.wantHi, safeSlice(c.seq, lo, hi))
			}
		})
	}
}

// safeSlice returns seq[lo:hi] guarding against the degenerate drop range.
func safeSlice(seq string, lo, hi int) string {
	if lo < 0 || hi > len(seq) || lo >= hi {
		return ""
	}
	return seq[lo:hi]
}

// TestUnitDisableQualityFiltering checks that DisableQualityFiltering bypasses
// BOTH the N-base/N-percent check and the quality-percentage check (upstream
// gates all of these on qualfilter.enabled, filter.cpp:43-50), while a read
// that would otherwise be dropped survives unchanged.
func TestUnitDisableQualityFiltering(t *testing.T) {
	enc := fastq.Phred33
	// 20bp read with 10 N's — well over the default n_base_limit of 5 — and a
	// quality string that also fails the default quality-percentage filter.
	seq := strings.Repeat("N", 10) + strings.Repeat("A", 10)
	qual := strings.Repeat("#", 20) // phred 2, all below threshold
	rec := &fastq.Record{ID: "r", Sequence: []byte(seq), Quality: []byte(qual)}

	t.Run("filter_on_drops", func(t *testing.T) {
		opts := DefaultProcessOptions()
		opts.MinLength = 1
		opts.LengthRequired = 1
		_, pass := filterRecord(clone(rec), opts, &ProcessStats{}, enc)
		if pass {
			t.Fatalf("expected the N-heavy/low-q read to be dropped by the quality filter")
		}
	})
	t.Run("disable_quality_keeps", func(t *testing.T) {
		opts := DefaultProcessOptions()
		opts.MinLength = 1
		opts.LengthRequired = 1
		opts.DisableQualityFiltering = true
		out, pass := filterRecord(clone(rec), opts, &ProcessStats{}, enc)
		if !pass {
			t.Fatalf("expected the read to survive with --disable_quality_filtering")
		}
		if string(out.Sequence) != seq {
			t.Fatalf("survivor bytes changed: got %q want %q", out.Sequence, seq)
		}
	})
}

// TestUnitDisableLengthFiltering checks that DisableLengthFiltering bypasses
// the too-short (length_required) and too-long (max_length) checks
// (filter.cpp:52-57) while leaving the quality filter intact.
func TestUnitDisableLengthFiltering(t *testing.T) {
	enc := fastq.Phred33
	// 8bp high-quality read, shorter than the default length_required of 15.
	rec := &fastq.Record{ID: "r", Sequence: []byte("ACGTACGT"), Quality: []byte("IIIIIIII")}

	t.Run("filter_on_drops_short", func(t *testing.T) {
		opts := DefaultProcessOptions() // LengthRequired/MinLength = 15
		_, pass := filterRecord(clone(rec), opts, &ProcessStats{}, enc)
		if pass {
			t.Fatalf("expected the 8bp read to be dropped by the length filter")
		}
	})
	t.Run("disable_length_keeps_short", func(t *testing.T) {
		opts := DefaultProcessOptions()
		opts.DisableLengthFiltering = true
		out, pass := filterRecord(clone(rec), opts, &ProcessStats{}, enc)
		if !pass {
			t.Fatalf("expected the short read to survive with --disable_length_filtering")
		}
		if string(out.Sequence) != "ACGTACGT" {
			t.Fatalf("survivor bytes changed: got %q", out.Sequence)
		}
	})
	t.Run("too_long_bypassed", func(t *testing.T) {
		long := &fastq.Record{ID: "r", Sequence: []byte(strings.Repeat("A", 50)), Quality: []byte(strings.Repeat("I", 50))}
		opts := DefaultProcessOptions()
		opts.MaxLength = 20 // would drop a 50bp read
		if _, pass := filterRecord(clone(long), opts, &ProcessStats{}, enc); pass {
			t.Fatalf("expected the 50bp read to be dropped by max_length")
		}
		opts.DisableLengthFiltering = true
		if _, pass := filterRecord(clone(long), opts, &ProcessStats{}, enc); !pass {
			t.Fatalf("expected the 50bp read to survive with --disable_length_filtering")
		}
	})
}

// TestUnitDisableTrimPolyG checks that DisableTrimPolyG suppresses poly-G
// trimming even when TrimPolyG is requested (the -G/--disable_trim_poly_g
// override).
func TestUnitDisableTrimPolyG(t *testing.T) {
	enc := fastq.Phred33
	seq := strings.Repeat("A", 20) + strings.Repeat("G", 15)
	qual := strings.Repeat("I", 35)
	mk := func() *fastq.Record {
		return &fastq.Record{ID: "r", Sequence: []byte(seq), Quality: []byte(qual)}
	}

	opts := DefaultProcessOptions()
	opts.TrimPolyG = true
	opts.PolyGMinLen = 10
	trimmed := trimRecord(mk(), opts, &ProcessStats{}, enc)
	if len(trimmed.Sequence) != 20 {
		t.Fatalf("expected poly-G trim to 20bp, got %d", len(trimmed.Sequence))
	}

	opts.DisableTrimPolyG = true
	untrimmed := trimRecord(mk(), opts, &ProcessStats{}, enc)
	if len(untrimmed.Sequence) != 35 {
		t.Fatalf("expected --disable_trim_poly_g to keep all 35bp, got %d", len(untrimmed.Sequence))
	}
}

// clone returns a deep copy of rec so a filterRecord/trimRecord call that
// re-slices does not mutate the shared fixture across sub-tests.
func clone(rec *fastq.Record) *fastq.Record {
	s := make([]byte, len(rec.Sequence))
	q := make([]byte, len(rec.Quality))
	copy(s, rec.Sequence)
	copy(q, rec.Quality)
	return &fastq.Record{ID: rec.ID, Description: rec.Description, Sequence: s, Quality: q}
}
