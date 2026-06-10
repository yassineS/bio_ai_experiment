// Overrepresented-sequence analysis (upstream fastp -p / -P).
//
// This is a verbatim Go port of upstream fastp's two-phase ORA:
//
//   1. Evaluator::computeOverRepSeq (evaluator.cpp:65-156) scans the first
//      ~1.51M bases of the input, collecting fixed-length substrings at
//      step sizes {10,20,40,100,min(150,seqLen-2)} and keeping those that
//      exceed length-specific count thresholds, then dropping substrings
//      that are dominated by a longer hot sequence. The survivors form the
//      candidate "hot sequence" map.
//
//   2. Stats::statRead (stats.cpp:311-329) samples 1-in-N reads during the
//      main pass; for each sampled read it scans the same step sizes and
//      increments the per-sequence hit count (and a per-cycle distribution)
//      whenever a candidate hot sequence is matched.
//
// The JSON output applies Stats::overRepPassed (stats.cpp:551-565) to keep
// only sequences whose sampled count, scaled by the sampling rate, clears a
// length-specific bar.

package fastp

import (
	"sort"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fastq"
)

// overrepAnalyzer accumulates overrepresented-sequence statistics for a
// single read stream (R1 or R2). It is seeded with a candidate hot-sequence
// map from buildHotSeqs and updated per sampled read via sampleRead.
type overrepAnalyzer struct {
	sampling     int
	evaluatedLen int
	// counts maps a candidate hot sequence to its sampled hit count.
	counts map[string]int64
	// dist maps a candidate hot sequence to its per-cycle hit distribution
	// (length evaluatedLen), populated when JSON/HTML reporting needs it.
	dist map[string][]int64
	// reads counts reads observed (sampled and unsampled) so the 1-in-N
	// sampling cadence matches upstream's mReads % sampling.
	reads int64
}

// overrepSteps returns upstream's five substring step sizes for a read of
// evaluated length seqLen: {10, 20, 40, 100, min(150, seqLen-2)}.
func overrepSteps(seqLen int) [5]int {
	last := 150
	if seqLen-2 < last {
		last = seqLen - 2
	}
	return [5]int{10, 20, 40, 100, last}
}

// buildHotSeqs is a verbatim port of Evaluator::computeOverRepSeq
// (evaluator.cpp:65-156). It scans records (already-read sequences) up to
// the upstream 1.51M-base limit, tallies fixed-length substrings, keeps
// those above the length-specific count thresholds, and removes substrings
// dominated by a longer hot sequence (count/count2 < 10). seqLen is the
// evaluated read length for this stream.
func buildHotSeqs(seqs []string, seqLen int) map[string]int64 {
	const baseLimit = 151 * 10000

	seqCounts := make(map[string]int64)
	var bases int64
	steps := overrepSteps(seqLen)

	for _, seq := range seqs {
		if bases >= baseLimit {
			break
		}
		rlen := len(seq)
		bases += int64(rlen)
		for _, step := range steps {
			if step <= 0 {
				continue
			}
			for i := 0; i < rlen-step; i++ {
				seqCounts[seq[i:i+step]]++
			}
		}
	}

	hot := make(map[string]int64)
	for seq, count := range seqCounts {
		switch {
		case len(seq) >= seqLen-1:
			if count >= 3 {
				hot[seq] = count
			}
		case len(seq) >= 100:
			if count >= 5 {
				hot[seq] = count
			}
		case len(seq) >= 40:
			if count >= 20 {
				hot[seq] = count
			}
		case len(seq) >= 20:
			if count >= 100 {
				hot[seq] = count
			}
		case len(seq) >= 10:
			if count >= 500 {
				hot[seq] = count
			}
		}
	}

	// Remove substrings dominated by a longer hot sequence: drop seq if some
	// other hot seq2 contains it and count/count2 < 10.
	for seq, count := range hot {
		for seq2, count2 := range hot {
			if seq == seq2 {
				continue
			}
			if count2 != 0 && contains(seq2, seq) && count/count2 < 10 {
				delete(hot, seq)
				break
			}
		}
	}
	return hot
}

// contains reports whether substr occurs within s. It is a small stdlib-free
// helper used by buildHotSeqs (strings.Contains would also work; kept inline
// to mirror the C++ find != npos check exactly).
func contains(s, substr string) bool {
	if len(substr) == 0 {
		return true
	}
	for i := 0; i+len(substr) <= len(s); i++ {
		if s[i:i+len(substr)] == substr {
			return true
		}
	}
	return false
}

// newOverrepAnalyzer builds an analyzer seeded with the candidate hot
// sequences from seqs. evaluatedLen is the per-stream read length used for
// step sizing and the per-cycle distribution width.
func newOverrepAnalyzer(seqs []string, sampling, evaluatedLen int) *overrepAnalyzer {
	hot := buildHotSeqs(seqs, evaluatedLen)
	a := &overrepAnalyzer{
		sampling:     sampling,
		evaluatedLen: evaluatedLen,
		counts:       make(map[string]int64, len(hot)),
		dist:         make(map[string][]int64, len(hot)),
	}
	for seq := range hot {
		a.counts[seq] = 0
		a.dist[seq] = make([]int64, evaluatedLen)
	}
	return a
}

// sampleRead is a verbatim port of the per-read ORA update in
// Stats::statRead (stats.cpp:311-329). It increments the read counter and,
// on every 1-in-sampling read, scans the upstream step sizes and bumps the
// candidate hit count and per-cycle distribution for any matched hot
// sequence.
func (a *overrepAnalyzer) sampleRead(seq string) {
	if a == nil {
		return
	}
	if a.reads%int64(a.sampling) == 0 {
		steps := overrepSteps(a.evaluatedLen)
		rlen := len(seq)
		for _, step := range steps {
			if step <= 0 {
				continue
			}
			for i := 0; i < rlen-step; i++ {
				sub := seq[i : i+step]
				if _, ok := a.counts[sub]; ok {
					a.counts[sub]++
					for p := i; p < i+step && p < a.evaluatedLen; p++ {
						a.dist[sub][p]++
					}
					i += step
				}
			}
		}
	}
	a.reads++
}

// overRepPassed is a verbatim port of Stats::overRepPassed
// (stats.cpp:551-565). It reports whether a sampled sequence's count, scaled
// by the sampling rate, clears the length-specific reporting bar.
func overRepPassed(seq string, count int64, sampling int) bool {
	s := int64(sampling)
	switch len(seq) {
	case 10:
		return s*count > 500
	case 20:
		return s*count > 200
	case 40:
		return s*count > 100
	case 100:
		return s*count > 50
	default:
		return s*count > 20
	}
}

// passedSequences returns the overrepresented sequences (and their sampled
// counts) that clear overRepPassed, sorted by sequence for deterministic
// output. Returns nil when none pass.
func (a *overrepAnalyzer) passedSequences() map[string]int64 {
	if a == nil {
		return nil
	}
	out := make(map[string]int64)
	for seq, count := range a.counts {
		if overRepPassed(seq, count, a.sampling) {
			out[seq] = count
		}
	}
	if len(out) == 0 {
		return nil
	}
	return out
}

// sortedSeqs returns the keys of the passed-sequence map in deterministic
// (sorted) order, used by the report writers.
func sortedSeqs(m map[string]int64) []string {
	keys := make([]string, 0, len(m))
	for k := range m {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	return keys
}

// recordSeqs collects the sequences of records as strings, bounded by the
// upstream 1.51M-base evaluation limit, for seeding buildHotSeqs.
func recordSeqs(records []*fastq.Record) []string {
	const baseLimit = 151 * 10000
	out := make([]string, 0, len(records))
	var bases int64
	for _, r := range records {
		if r == nil || bases >= baseLimit {
			break
		}
		out = append(out, string(r.Sequence))
		bases += int64(len(r.Sequence))
	}
	return out
}
