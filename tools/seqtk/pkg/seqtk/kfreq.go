// kfreq.go: implementation of `seqtk kfreq`.
//
// Upstream source: reference_code/seqtk/seqtk.c::stk_kfreq (v1.5-r133,
// lines 1777-1839).
//
// Behaviour: given a single ACGT k-mer on the command line, walk every
// FASTA record and count two things per record:
//
//   - the number of exact occurrences of the k-mer (and of its
//     reverse complement) — upstream's cnt[0] / cnt[1].
//   - the number of "neighbour" k-mer occurrences (and of its
//     reverse-complement neighbours) — upstream's cnt_nei[0] /
//     cnt_nei[1]. A neighbour is any k-mer at Hamming distance ≤ 1
//     from the target (i.e. the target itself plus all 3·k single-base
//     substitutions; upstream's neighbour set inclusion is implicit
//     because the target is also reachable by substituting any base
//     with itself).
//
// One TSV row per record is emitted:
//
//	name\tlen\t<strand>\t<cnt_nei[which]>\t<cnt[which]>
//
// where `which` is 0 (forward, '+') if cnt_nei[0] > cnt_nei[1], else
// 1 (reverse, '-'). Empty records (`seq.l == 0`) still emit a row with
// counts 0 and `which = 1` ('-'); this matches upstream where both
// neighbour counts start at 0 and the tie-break falls through to 1.
//
// Upstream getopt surface: none — the program takes two positional
// arguments, `<kmer>` and `<in.fa>` (use `-` for stdin). Upstream
// `assert(c >= 1 && c <= 4)` aborts on a kmer with a non-ACGT byte;
// we return a typed error so the Go CLI can exit cleanly.
//
// The "-o/--output FILE" flag wired in cmd/ is the project-wide
// Go-port convenience and does not affect parity.

package seqtk

import (
	"bufio"
	"fmt"
	"io"
)

// KfreqOptions carries the parsed CLI inputs for Kfreq. There are no
// real upstream flags; the k-mer is a required positional argument and
// is exposed here as a typed string.
type KfreqOptions struct {
	// Kmer is the target k-mer; every byte must be A, C, G or T
	// (case-insensitive). Empty or non-ACGT k-mers are rejected.
	Kmer string
}

// Kfreq scans every FASTA record in r and writes one TSV row per
// record to w, matching upstream `seqtk kfreq` byte-for-byte.
func Kfreq(r io.Reader, w io.Writer, opts KfreqOptions) error {
	l := len(opts.Kmer)
	if l == 0 {
		return fmt.Errorf("kfreq: kmer must be non-empty")
	}
	// 2*l bits of state — upstream uses a signed int for the rolling
	// window and the neighbour bitmap, which caps l at 15 in practice.
	// Reject larger kmers explicitly rather than silently overflowing.
	if l > 15 {
		return fmt.Errorf("kfreq: kmer length must be ≤ 15 (got %d)", l)
	}

	// Encode the k-mer into a 2-bit-per-base integer using the same
	// (c-1) shift upstream uses (A=0, C=1, G=2, T=3). Any non-ACGT
	// byte aborts upstream; we return an error.
	var kmer int
	for i := 0; i < l; i++ {
		c := seqNT6Table[opts.Kmer[i]]
		if c < 1 || c > 4 {
			return fmt.Errorf("kfreq: kmer contains non-ACGT byte %q at position %d", opts.Kmer[i], i)
		}
		kmer = kmer<<2 | int(c-1)
	}
	mask := (1 << (2 * l)) - 1

	// Build the neighbour set: for each base position i in [0, l), set
	// nei[kmer with base i replaced by b] = 1 for each b in [0, 4).
	// This makes the target k-mer itself a neighbour too.
	nei := make([]byte, 1<<(2*l))
	for i := 0; i < l; i++ {
		x := kmer &^ (3 << (2 * i))
		for j := 0; j < 4; j++ {
			nei[x|j<<(2*i)] = 1
		}
	}

	bw := bufio.NewWriter(w)
	err := scanFASTA(r, func(name string, seq []byte) error {
		var x0, x1 int      // forward / reverse rolling encodings
		var k int           // number of consecutive ACGT bases seen
		var cnt [2]int64    // exact hits, [forward, reverse]
		var cntNei [2]int64 // neighbour hits, [forward, reverse]
		for i := 0; i < len(seq); i++ {
			c := seqNT6Table[seq[i]]
			if c >= 1 && c <= 4 {
				x0 = (x0<<2 | int(c-1)) & mask
				// Reverse complement rolling encoding: shift in the
				// complement (4-c) at the high end. Upstream writes
				// `x[1] = x[1] >> 2 | (4-c) << 2*(l-1)`; we mask the
				// result to stay inside 2*l bits.
				x1 = (x1>>2 | int(4-c)<<(2*(l-1))) & mask
				if k < l {
					k++
				}
				if k == l {
					if x0 == kmer {
						cnt[0]++
					} else if x1 == kmer {
						cnt[1]++
					}
					if nei[x0] != 0 {
						cntNei[0]++
					} else if nei[x1] != 0 {
						cntNei[1]++
					}
				}
			} else {
				k = 0
			}
		}
		// Tie-break matches upstream: `cnt_nei[0] > cnt_nei[1] ? 0 : 1`.
		// Equal counts (including the all-zero case) pick reverse ('-').
		which := 1
		if cntNei[0] > cntNei[1] {
			which = 0
		}
		strand := byte('-')
		if which == 0 {
			strand = '+'
		}
		_, err := fmt.Fprintf(bw, "%s\t%d\t%c\t%d\t%d\n",
			name, len(seq), strand, cntNei[which], cnt[which])
		return err
	})
	if err != nil {
		return err
	}
	return bw.Flush()
}
