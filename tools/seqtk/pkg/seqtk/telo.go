// telo.go: implementation of `seqtk telo`.
//
// Upstream source: reference_code/seqtk/seqtk.c::stk_telo (v1.5-r133,
// lines 1969-2062).
//
// Behaviour: for every FASTA record, run a banded-extension scan from
// both ends looking for the telomeric motif (default "CCCTAA"). The 5'
// scan walks the sequence left-to-right, encoding each base as a 2-bit
// rolling integer and querying a hash set built from every cyclic
// rotation of the motif. The 3' scan walks right-to-left and instead
// queries rotations of the motif's reverse complement (the same hash
// set, fed with the (4-c)/complement-encoded byte). Scoring: +1 for a
// motif rotation match, -penalty otherwise, started only after the
// rolling window has been filled (i.e. after i >= motif length). The
// scan halts when the running score drops more than `-d max_drop`
// below its running maximum; if the maximum reaches `-s min_score`,
// the corresponding interval is emitted.
//
// Output formats (TAB-separated):
//
//	"<name>\t0\t<5' end pos>\t<seq len>\n"          (5' hit)
//	"<name>\t<3' start pos>\t<seq len>\t<seq len>\n" (3' hit)
//
// With -P, instead of BED intervals each iteration emits a per-step
// "profile" row to stdout (one row per i, two rows per record per
// end): "P\t<name>\t<i>\t<score>\t<max>\n" for the 5' scan and
// "Q\t<name>\t<seq.l - i>\t<score>\t<max>\n" for the 3' scan.
//
// After all records, a summary "<sum_telo>\t<sum_input>\n" line is
// written to stderr (NOT stdout) — exactly like upstream.
//
// Upstream getopt surface (seqtk.c:1978): "m:p:d:s:P"
//
//	-m STR   motif [CCCTAA]
//	-p INT   penalty per non-hit position [1]; negative values are
//	         flipped to positive (upstream `if (penalty < 0) penalty = -penalty`).
//	-d INT   max score drop before the scan aborts [2000]
//	-s INT   min running max for an interval to be emitted [300]
//	-P       emit per-position scoring profile instead of BED intervals
//
// The "-o/--output FILE" flag wired in cmd/ is the project-wide
// Go-port convenience and does not affect parity.

package seqtk

import (
	"bufio"
	"fmt"
	"io"
)

// TeloOptions configures Telo. The defaults mirror upstream
// (motif=CCCTAA, penalty=1, max-drop=2000, min-score=300, P=false).
type TeloOptions struct {
	Motif       string // -m
	Penalty     int    // -p
	MaxDrop     int    // -d
	MinScore    int    // -s
	ShowProfile bool   // -P
}

// Default values for TeloOptions, matching upstream defaults at
// seqtk.c:1974.
const (
	DefaultTeloMotif    = "CCCTAA"
	DefaultTeloPenalty  = 1
	DefaultTeloMaxDrop  = 2000
	DefaultTeloMinScore = 300
)

// Telo writes BED-style intervals (or, with -P, a scoring profile) to
// stdoutW and the upstream "<sum_telo>\t<sum_input>\n" summary line to
// stderrW. The split between stdout and stderr is deliberate and
// matches upstream byte-for-byte.
func Telo(r io.Reader, stdoutW, stderrW io.Writer, opts TeloOptions) error {
	if opts.Motif == "" {
		return fmt.Errorf("telo: motif must be non-empty")
	}
	mlen := len(opts.Motif)
	// Upstream stores the rolling window in a uint64, so the motif
	// must encode into ≤ 64 bits (i.e. mlen ≤ 32). Practically motifs
	// are 6-8 bp; reject anything pathological to keep behaviour
	// well-defined.
	if mlen > 32 {
		return fmt.Errorf("telo: motif length must be ≤ 32 (got %d)", mlen)
	}
	penalty := opts.Penalty
	if penalty < 0 {
		penalty = -penalty
	}

	var mask uint64
	if mlen == 32 {
		mask = ^uint64(0)
	} else {
		mask = (uint64(1) << uint(2*mlen)) - 1
	}

	// Build a set containing every cyclic rotation of the motif (the
	// 2-bit-per-base encoding). Upstream uses a khash; a Go map[uint64]struct{}
	// is functionally identical.
	rot := make(map[uint64]struct{}, mlen*2)
	for i := 0; i < mlen; i++ {
		var x uint64
		for j := 0; j < mlen; j++ {
			c := seqNT6Table[opts.Motif[(i+j)%mlen]]
			if c < 1 || c > 4 {
				return fmt.Errorf("telo: motif contains non-ACGT byte %q at position %d", opts.Motif[(i+j)%mlen], (i+j)%mlen)
			}
			x = x<<2 | uint64(c-1)
		}
		rot[x] = struct{}{}
	}

	bw := bufio.NewWriter(stdoutW)
	var sumInput, sumTelo int64

	err := scanFASTA(r, func(name string, seq []byte) error {
		l := len(seq)
		sumInput += int64(l)

		// --- 5' scan (forward, left-to-right) ---
		var score, max int64
		var maxI int = -1
		var st int = 0
		{
			var x uint64
			var k int
			for i := 0; i < l; i++ {
				hit := 0
				c := seqNT6Table[seq[i]]
				if c >= 1 && c <= 4 {
					x = (x<<2 | uint64(c-1)) & mask
					k++
					if k >= mlen {
						if _, ok := rot[x]; ok {
							hit = 1
						}
					}
				} else {
					k = 0
					x = 0
				}
				if opts.ShowProfile {
					if _, err := fmt.Fprintf(bw, "P\t%s\t%d\t%d\t%d\n", name, i, score, max); err != nil {
						return err
					}
				}
				if i >= mlen {
					if hit != 0 {
						score++
					} else {
						score -= int64(penalty)
					}
				}
				if score > max {
					max = score
					maxI = i
				} else if max-score > int64(opts.MaxDrop) {
					break
				}
			}
			if max >= int64(opts.MinScore) {
				if !opts.ShowProfile {
					if _, err := fmt.Fprintf(bw, "%s\t0\t%d\t%d\n", name, maxI+1, l); err != nil {
						return err
					}
				}
				sumTelo += int64(maxI + 1)
				st = maxI + 1
			}
		}

		// --- 3' scan (reverse, right-to-left) ---
		score, max = 0, 0
		maxI = -1
		{
			var x uint64
			var k int
			for i := l - 1; i >= st; i-- {
				hit := 0
				c := seqNT6Table[seq[i]]
				if c >= 1 && c <= 4 {
					// Reverse complement: shift in (4-c) at the low end.
					x = (x<<2 | uint64(4-c)) & mask
					k++
					if k >= mlen {
						if _, ok := rot[x]; ok {
							hit = 1
						}
					}
				} else {
					k = 0
					x = 0
				}
				if opts.ShowProfile {
					if _, err := fmt.Fprintf(bw, "Q\t%s\t%d\t%d\t%d\n", name, l-i, score, max); err != nil {
						return err
					}
				}
				if l-i >= mlen {
					if hit != 0 {
						score++
					} else {
						score -= int64(penalty)
					}
				}
				if score > max {
					max = score
					maxI = i
				} else if max-score > int64(opts.MaxDrop) {
					break
				}
			}
			if max >= int64(opts.MinScore) {
				if !opts.ShowProfile {
					if _, err := fmt.Fprintf(bw, "%s\t%d\t%d\t%d\n", name, maxI, l, l); err != nil {
						return err
					}
				}
				sumTelo += int64(l - maxI)
			}
		}
		return nil
	})
	if err != nil {
		return err
	}
	if err := bw.Flush(); err != nil {
		return err
	}
	// Upstream prints this with `%ld\t%ld\n` to stderr (telomere bases first).
	_, err = fmt.Fprintf(stderrW, "%d\t%d\n", sumTelo, sumInput)
	return err
}
