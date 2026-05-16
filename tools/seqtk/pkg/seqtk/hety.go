// hety.go: implementation of `seqtk hety`.
//
// Upstream source: reference_code/seqtk/seqtk.c::stk_hety (v1.5-r133,
// lines 567-628). Behaviour: walks each FASTA record in non-overlapping
// windows of -w bases stepped at win_size/n_start, counts heterozygous
// IUPAC codes (2-base ambiguity codes R, Y, S, W, K, M) per window, and
// emits one TSV line per window that contains at least one ACGT or
// 2-base IUPAC code.
//
// Per-line columns (tab-separated, matches upstream stk_hety byte-for-byte):
//
//	name\tstart\tend\t<scaled-het>\t<n_hom+n_het>\t<n_het>
//
// where the third column is "(cnt[2] / (cnt[1] + cnt[2])) * win_size"
// printed with %.2f, cnt[1] is the homozygous count and cnt[2] is the
// heterozygous count over the window. start / end are 0-based; end is
// the index of the byte AFTER the window's last base (so the very last
// emit per record uses end == record-length, even when that length is
// not a multiple of win_size).
//
// Upstream getopt surface (seqtk.c:584):  "w:t:m"
//
//   -w INT   window size [50000]
//   -t INT   number of start positions in a window [5]
//   -m       treat lowercase bases as masked (i.e. as N)
//
// The "-o/--output FILE" flag added in cmd/ is the project-wide Go-port
// convenience and does not affect parity.

package seqtk

import (
	"bufio"
	"fmt"
	"io"
)

// HetyOptions configures Hety. Defaults mirror upstream's stk_hety
// (win_size=50000, n_start=5, is_lower_mask=false).
type HetyOptions struct {
	WinSize     int  // -w: window size; must be > 0
	NStart      int  // -t: number of start positions inside a window; must be > 0
	IsLowerMask bool // -m: treat lowercase bases as masked (count as N)
}

// Default values for HetyOptions, matching upstream defaults at seqtk.c:572.
const (
	DefaultHetyWinSize = 50000
	DefaultHetyNStart  = 5
)

// hetyClass classifies a single byte into upstream's per-position "z"
// value:
//
//	0  -> not counted (N, X, 3-/4-base IUPAC codes, or lowercase when
//	      IsLowerMask=true)
//	1  -> homozygous unambiguous base (A, C, G, T — bitcnt == 1)
//	2  -> heterozygous 2-base IUPAC code (R, Y, S, W, K, M — bitcnt == 2)
//
// The classification mirrors upstream lines 614-619: it looks the byte
// up in seq_nt16_table, takes its bit-count, and maps {>2 -> 0, 2 -> 2,
// 1 -> 1}.
func hetyClass(c byte, isLowerMask bool) byte {
	if isLowerMask && c >= 'a' && c <= 'z' {
		c = 'N'
	}
	x := bitCntTable[seqNT16Table[c]]
	switch {
	case x > 2:
		return 0
	case x == 2:
		return 2
	default:
		return 1
	}
}

// Hety scans every FASTA record in r and writes one TSV line per
// non-overlapping window to w using upstream "seqtk hety" formatting
// (see the package-level comment). Empty windows (zero ACGT + zero
// 2-base IUPAC codes) are silently dropped.
//
// opts.WinSize and opts.NStart must both be > 0; otherwise an error is
// returned (matching upstream's implicit guarantee that win_step ==
// win_size / n_start does not divide by zero or yield zero).
func Hety(r io.Reader, w io.Writer, opts HetyOptions) error {
	if opts.WinSize <= 0 {
		return fmt.Errorf("hety: -w window size must be > 0 (got %d)", opts.WinSize)
	}
	if opts.NStart <= 0 {
		return fmt.Errorf("hety: -t n_start must be > 0 (got %d)", opts.NStart)
	}
	winStep := opts.WinSize / opts.NStart
	if winStep <= 0 {
		return fmt.Errorf("hety: window step (win_size/n_start) must be > 0; got win=%d n_start=%d", opts.WinSize, opts.NStart)
	}

	bw := bufio.NewWriter(w)
	buf := make([]byte, opts.WinSize)

	err := scanFASTA(r, func(name string, seq []byte) error {
		var cnt [3]int64
		next := 0
		l := len(seq)
		// Iterate i in [0, l]; the i == l iteration is the upstream
		// "virtual" terminator that flushes any open window and is
		// the only path through the special l >= win_size tail
		// adjustment.
		for i := 0; i <= l; i++ {
			// Emit-window check: at every win_step boundary past
			// the first full window, AND at i == l.
			if (i >= opts.WinSize && i%winStep == 0) || i == l {
				// Tail adjustment at end-of-sequence: when the
				// sequence is at least one full window long,
				// upstream "rolls back" the count for the bytes
				// that were already evicted from buf but are
				// still inside the final window. The loop below
				// walks the indices [l - win_size, next) and
				// decrements cnt[buf[y % win_size]] for each.
				if i == l && l >= opts.WinSize {
					for y := l - opts.WinSize; y < next; y++ {
						cnt[buf[y%opts.WinSize]]--
					}
				}
				if cnt[1]+cnt[2] > 0 {
					scaled := float64(cnt[2]) / float64(cnt[1]+cnt[2]) * float64(opts.WinSize)
					if _, err := fmt.Fprintf(bw, "%s\t%d\t%d\t%.2f\t%d\t%d\n",
						name, next, i, scaled, cnt[1]+cnt[2], cnt[2]); err != nil {
						return err
					}
				}
				next = i
			}
			if i < l {
				y := i % opts.WinSize
				z := hetyClass(seq[i], opts.IsLowerMask)
				if i >= opts.WinSize {
					cnt[buf[y]]--
				}
				buf[y] = z
				cnt[z]++
			}
		}
		return nil
	})
	if err != nil {
		return err
	}
	return bw.Flush()
}
