// listhet.go: implementation of `seqtk listhet`.
//
// Upstream source: reference_code/seqtk/seqtk.c::stk_listhet (v1.5-r133,
// lines 1004-1029). Despite the name's VCF-ish ring, this subcommand reads
// FASTA only and uses IUPAC ambiguity codes to identify heterozygous sites
// — exactly the 2-base codes R, Y, S, W, K, M (uppercase or lowercase).
// Per-position output is one TSV row per heterozygous base:
//
//	name\t<1-based pos>\t<byte as-is>\n
//
// Upstream getopt surface: none (positional `<in.fa>` only —
// `stk_listhet` does not call `getopt`).
//
// The Go cmd/ layer adds an `-o/--output FILE` flag as the project-wide
// "write to file instead of stdout" convenience. That's the only addition
// beyond the upstream surface; the algorithm itself is a verbatim port.

package seqtk

import (
	"bufio"
	"fmt"
	"io"
)

// ListHet walks every FASTA record in r and writes one TSV row per byte
// whose IUPAC popcount is exactly 2 (i.e. R, Y, S, W, K, M and their
// lowercase counterparts). Output goes to w in the upstream format
// `name\tpos1based\tbyte\n`. The byte is emitted in its original case;
// the position is 1-based, matching upstream `printf("%s\t%ld\t%c\n",
// name, i+1, b)`.
//
// Empty sequences emit no rows. Records with no heterozygous bytes emit
// no rows. The 3-/4-base ambiguity codes (B, D, H, V, N, X) and the
// unambiguous bases (A, C, G, T) are silently skipped — they have
// popcounts of 3, 4, 1 and 0/4 respectively, none of which equal 2.
func ListHet(r io.Reader, w io.Writer) error {
	bw := bufio.NewWriter(w)
	err := scanFASTA(r, func(name string, seq []byte) error {
		for i := 0; i < len(seq); i++ {
			b := seq[i]
			if bitCntTable[seqNT16Table[b]] == 2 {
				if _, err := fmt.Fprintf(bw, "%s\t%d\t%c\n", name, i+1, b); err != nil {
					return err
				}
			}
		}
		return nil
	})
	if err != nil {
		return err
	}
	return bw.Flush()
}
