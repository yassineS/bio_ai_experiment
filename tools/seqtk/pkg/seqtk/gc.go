// gc.go: implementation of `seqtk gc` — find GC-rich (or AT-rich) regions in
// a FASTA file using upstream seqtk's X-dropoff scoring scheme.
//
// Upstream source: reference_code/seqtk/seqtk.c::stk_gc.
//
// The algorithm is NOT a sliding window; it is a one-pass greedy scan where
// every "hit" position contributes +q = (1-frac)/frac to the score and every
// non-hit contributes -1. A region runs from the first hit position (`start`)
// to the position that achieved the running maximum (`max_i`). The scan ends
// the region when the score either drops below zero or falls X (the X-dropoff)
// below the running maximum — at which point, if the region is long enough,
// it is emitted and the scan resumes from `max_i + 1`.
//
// Hit definition uses upstream's seq_nt16_table (see comp.go's seqNT16Table):
//   - GC mode (default): a hit is C(2), G(4) or S(6).
//   - AT mode (`-w`):    a hit is A(1), T(8) or W(9).

package seqtk

import (
	"bufio"
	"fmt"
	"io"
)

// GCOptions configures GC. The defaults match upstream `seqtk gc` (v1.5).
type GCOptions struct {
	// MinLength is the minimum length of a reported region (`-l`, default 20).
	MinLength int
	// MinFrac is the target hit fraction for the X-dropoff score: every hit
	// contributes +(1-MinFrac)/MinFrac, every non-hit contributes -1
	// (`-f`, default 0.6).
	MinFrac float64
	// XDropoff is the maximum running-score drop below the current maximum
	// before a region is ended (`-x`, default 10.0).
	XDropoff float64
	// IsAT toggles the predicate from GC-rich to AT-rich (`-w`).
	IsAT bool
}

// Defaults for the GC subcommand, matching upstream v1.5.
const (
	DefaultGCMinLength = 20
	DefaultGCMinFrac   = 0.6
	DefaultGCXDropoff  = 10.0
)

// GC writes one BED4 row per high-GC (or high-AT, with opts.IsAT) region
// detected in r, using upstream `seqtk gc`'s X-dropoff scoring algorithm.
// Output format is chrom\tstart\tend\thits\n, 0-based half-open, matching
// upstream byte-for-byte (the trailing column is the number of hit positions
// in [start, end), i.e. the C `max_hits - start_hits + 1` expression).
//
// Returns an error if opts.MinFrac is not in (0, 1), opts.MinLength < 1, or
// on I/O errors. The X-dropoff is unconstrained (upstream accepts negative
// values too, where the loop will essentially never end a region).
func GC(r io.Reader, w io.Writer, opts GCOptions) error {
	if opts.MinLength < 1 {
		return fmt.Errorf("gc: -l/--min-length must be >= 1 (got %d)", opts.MinLength)
	}
	if !(opts.MinFrac > 0 && opts.MinFrac < 1) {
		return fmt.Errorf("gc: -f/--min-frac must be in (0, 1) (got %g)", opts.MinFrac)
	}
	q := (1.0 - opts.MinFrac) / opts.MinFrac
	bw := bufio.NewWriter(w)

	err := scanFASTA(r, func(name string, seq []byte) error {
		var (
			start, maxI      int64
			nHits, startHits int64
			maxHits          int64
			sc, maxScore     float64
		)
		for i := 0; i < len(seq); i++ {
			c := seqNT16Table[seq[i]]
			var hit bool
			if opts.IsAT {
				hit = c == 1 || c == 8 || c == 9
			} else {
				hit = c == 2 || c == 4 || c == 6
			}
			if hit {
				nHits++
				if sc == 0 {
					start = int64(i)
					startHits = nHits
				}
				sc += q
				if sc > maxScore {
					maxScore = sc
					maxI = int64(i)
					maxHits = nHits
				}
			} else if sc > 0 {
				sc += -1.0
				if sc < 0 || maxScore-sc > opts.XDropoff {
					if maxI+1-start >= int64(opts.MinLength) {
						if err := writeBED4Int(bw, name, int(start), int(maxI+1), int(maxHits-startHits+1)); err != nil {
							return err
						}
					}
					sc = 0
					maxScore = 0
					// Upstream resets the scan index to max_i then ++i
					// fires at the next loop iteration. Replicate that by
					// rewinding i.
					i = int(maxI)
				}
			}
		}
		if maxScore > 0 && maxI+1-start >= int64(opts.MinLength) {
			if err := writeBED4Int(bw, name, int(start), int(maxI+1), int(maxHits-startHits+1)); err != nil {
				return err
			}
		}
		return nil
	})
	if err != nil {
		return err
	}
	return bw.Flush()
}
