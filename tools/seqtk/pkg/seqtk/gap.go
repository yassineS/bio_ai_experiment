// gap.go: implementation of `seqtk gap`.
//
// Upstream source: reference_code/seqtk/seqtk.c::stk_gap.
// Behaviour: walks every FASTA record and emits a BED3 row
// (chrom\tstart\tend, 0-based half-open) for every maximal run of "gap" bytes
// (non-ACGT, case-insensitive — see IsGapByte) whose length is at least
// opts.MinSize. The end-of-sequence is treated as a terminator that flushes
// any open run, matching upstream's `i == len` sentinel branch.

package seqtk

import (
	"bufio"
	"io"
)

// GapOptions configures Gap. MinSize is the minimum length of a gap run to
// report (upstream's `-l`, default 50).
type GapOptions struct {
	MinSize int
}

// DefaultGapMinSize is the upstream default for `seqtk gap -l`.
const DefaultGapMinSize = 50

// Gap writes one BED3 row per qualifying gap (run of non-ACGT bytes of length
// >= opts.MinSize) found in any record of the FASTA stream r. Output goes to
// w in upstream "seqtk gap" format: chrom\tstart\tend\n, 0-based half-open.
//
// opts.MinSize <= 0 is silently treated like upstream `stk_gap`'s
// `l > 0 && l >= min_size` guard: every non-zero gap run is reported.
func Gap(r io.Reader, w io.Writer, opts GapOptions) error {
	bw := bufio.NewWriter(w)
	err := scanFASTA(r, func(name string, seq []byte) error {
		l := 0
		// Iterate i == len once as a virtual non-gap terminator, so the
		// last run is flushed exactly like upstream's `i == len` branch.
		for i := 0; i <= len(seq); i++ {
			isGap := i < len(seq) && IsGapByte(seq[i])
			if !isGap {
				// Mirror upstream's `if (l > 0 && l >= min_size)` guard:
				// never emit a zero-length BED row even when MinSize <= 0.
				if l > 0 && l >= opts.MinSize {
					if err := writeBED3(bw, name, i-l, i); err != nil {
						return err
					}
				}
				l = 0
			} else {
				l++
			}
		}
		return nil
	})
	if err != nil {
		return err
	}
	return bw.Flush()
}
