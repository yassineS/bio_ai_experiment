package difffuzz

import "time"

// stillDiverges is the predicate the minimizer shrinks against: it reports
// whether candidate input still triggers a divergence of the SAME class as the
// original (so we don't accidentally minimize toward a different bug).
type stillDiverges func(candidate []byte) bool

// Minimize delta-debugs input down to a smaller reproducer that still triggers
// the same divergence, per the predicate. It is a deterministic, dependency-
// free shrinker that interleaves two strategies:
//
//  1. line-granularity ddmin: remove chunks of whole lines (records), halving
//     the chunk size when no removal preserves the divergence — this collapses
//     a large fixture down to the few records that matter;
//  2. byte-granularity ddmin on the result: the same chunk-removal at byte
//     granularity, plus a trailing-truncation pass, to squeeze the final
//     reproducer.
//
// The predicate is assumed to be reasonably cheap (it re-runs both binaries);
// maxSteps bounds the number of predicate evaluations so a pathological input
// cannot make minimization run unbounded. Minimize never returns something the
// predicate rejects: if no shrink helps, it returns the original input.
func Minimize(input []byte, pred stillDiverges, maxSteps int) []byte {
	if maxSteps <= 0 {
		maxSteps = 2000
	}
	steps := 0
	tick := func() bool {
		steps++
		return steps <= maxSteps
	}

	cur := input
	if !pred(cur) {
		// The predicate doesn't even fire on the original; nothing to do.
		return input
	}

	// Pass 1: line-level ddmin.
	cur = ddmin(cur, splitKeepEmpty, joinLines, pred, tick)
	// Pass 2: byte-level ddmin.
	cur = ddminBytes(cur, pred, tick)
	// Pass 3: trailing truncation (cheap last squeeze).
	cur = truncShrink(cur, pred, tick)
	return cur
}

// ddmin removes chunks of units (lines) from cur while the predicate holds,
// classic delta-debugging: try removing each chunk at the current granularity;
// on any success restart at coarse granularity; when nothing can be removed at
// the finest granularity, stop.
func ddmin(cur []byte, split func([]byte) [][]byte, join func([][]byte) []byte,
	pred stillDiverges, tick func() bool) []byte {
	units := split(cur)
	n := 2
	for len(units) >= 2 {
		chunk := len(units) / n
		if chunk == 0 {
			chunk = 1
		}
		removedAny := false
		for start := 0; start < len(units); start += chunk {
			if !tick() {
				return join(units)
			}
			end := start + chunk
			if end > len(units) {
				end = len(units)
			}
			cand := make([][]byte, 0, len(units)-(end-start))
			cand = append(cand, units[:start]...)
			cand = append(cand, units[end:]...)
			if len(cand) == 0 {
				continue
			}
			if pred(join(cand)) {
				units = cand
				removedAny = true
				n = 2 // restart coarse after a successful removal
				break
			}
		}
		if removedAny {
			continue
		}
		if n >= len(units) {
			break
		}
		if chunk == 1 {
			break
		}
		n *= 2
	}
	return join(units)
}

// ddminBytes runs ddmin at byte granularity (each unit is a single byte).
func ddminBytes(cur []byte, pred stillDiverges, tick func() bool) []byte {
	split := func(b []byte) [][]byte {
		out := make([][]byte, len(b))
		for i := range b {
			out[i] = b[i : i+1]
		}
		return out
	}
	join := func(units [][]byte) []byte {
		out := make([]byte, 0, len(units))
		for _, u := range units {
			out = append(out, u...)
		}
		return out
	}
	return ddmin(cur, split, join, pred, tick)
}

// truncShrink repeatedly tries to drop the buffer's tail (halving the cut)
// while the divergence persists. Cheap and effective on truncation-triggered
// parser divergences.
func truncShrink(cur []byte, pred stillDiverges, tick func() bool) []byte {
	for cut := len(cur) / 2; cut >= 1; cut /= 2 {
		for len(cur) > cut {
			if !tick() {
				return cur
			}
			cand := cur[:len(cur)-cut]
			if len(cand) == 0 || !pred(cand) {
				break
			}
			cur = cand
		}
	}
	return cur
}

// timeoutDeadline is a small helper so callers can express a per-minimize wall
// budget instead of (or alongside) a step budget. It returns a tick function
// that also stops once the deadline passes.
func timeoutDeadline(maxSteps int, budget time.Duration) func() bool {
	deadline := time.Now().Add(budget)
	steps := 0
	return func() bool {
		steps++
		if maxSteps > 0 && steps > maxSteps {
			return false
		}
		return time.Now().Before(deadline)
	}
}
