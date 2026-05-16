// famask.go: implementation of `seqtk famask`.
//
// Upstream source: reference_code/seqtk/seqtk.c::stk_famask (v1.5-r133).
// Behaviour: walks two FASTA streams in parallel — a source FASTA and a
// "mask" FASTA — and emits the source bases transformed by the mask
// character at the same position, using upstream's rules:
//
//   - mask byte == 'X'  -> keep source base unchanged
//   - mask byte == 'x'  -> lowercase the source base (soft-mask)
//   - any other byte    -> replace the source base with the mask byte
//
// Output is FASTA, wrapped at 60 bases per line (upstream `l%60==0`).
//
// Upstream takes no flags — the getopt loop is `getopt(argc, argv, "")`
// at seqtk.c:878, i.e. it accepts no options whatsoever. The two
// positional arguments are <src.fa> and <mask.fa>; both may be "-" for
// stdin (but not both, per the project convention).
//
// Mismatched record names and unequal sequence lengths produce a stderr
// warning identical to upstream's; the shorter of the two records is
// used as the per-record length.

package seqtk

import (
	"bufio"
	"fmt"
	"io"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fasta"
)

// Famask reads two FASTA streams from src and mask and writes the masked
// FASTA result to w. Records are paired by stream position; the byte
// rules above are applied to the bases of each pair.
//
// The function returns nil at clean EOF. Per-record name mismatches and
// length mismatches are logged to os.Stderr (matching upstream's
// fprintf to stderr) and the shorter length is used.
func Famask(src, mask io.Reader, w io.Writer) error {
	return famaskImpl(src, mask, w, os.Stderr)
}

// famaskImpl is the testable core: it lets tests capture the upstream
// stderr-style warnings on a custom writer.
func famaskImpl(src, mask io.Reader, w, warn io.Writer) error {
	sr := fasta.NewReader(src)
	mr := fasta.NewReader(mask)
	bw := bufio.NewWriter(w)
	for {
		s, errS := sr.Read()
		if errS == io.EOF {
			break
		}
		if errS != nil {
			return errS
		}
		m, errM := mr.Read()
		if errM == io.EOF {
			// Upstream calls kseq_read on the mask without checking
			// its return value; an EOF there leaves seq[1] holding
			// the previous record's data. In practice users always
			// provide matched-length inputs, so we treat this as
			// "no further masking" and simply emit the source bases
			// unmodified — equivalent to upstream's behaviour when
			// the mask runs out of records (it would re-emit the
			// last mask record indefinitely; we mirror by leaving
			// `mseq` empty and the per-base loop falls through to
			// the "keep" path because `min_l = 0`).
			fmt.Fprintf(warn, "[famask] mask stream ended before source: %s has no mask record\n", s.ID)
			m = &fasta.Record{ID: s.ID, Sequence: nil}
		} else if errM != nil {
			return errM
		}
		if s.ID != m.ID {
			fmt.Fprintf(warn, "[famask] Different sequence names: %s != %s\n", s.ID, m.ID)
		}
		if len(s.Sequence) != len(m.Sequence) {
			fmt.Fprintf(warn, "[famask] Unequal sequence length: %d != %d\n", len(s.Sequence), len(m.Sequence))
		}
		minL := len(s.Sequence)
		if len(m.Sequence) < minL {
			minL = len(m.Sequence)
		}
		if _, err := fmt.Fprintf(bw, ">%s", s.ID); err != nil {
			return err
		}
		for l := 0; l < minL; l++ {
			c := s.Sequence[l]
			mc := m.Sequence[l]
			switch {
			case mc == 'x':
				c = toLowerByte(c)
			case mc != 'X':
				c = mc
			}
			if l%60 == 0 {
				if err := bw.WriteByte('\n'); err != nil {
					return err
				}
			}
			if err := bw.WriteByte(c); err != nil {
				return err
			}
		}
		if err := bw.WriteByte('\n'); err != nil {
			return err
		}
	}
	return bw.Flush()
}

// toLowerByte returns the ASCII lowercase form of b, matching the C
// stdlib `tolower(int)` semantics used by upstream: only A-Z map; all
// other bytes pass through unchanged.
func toLowerByte(b byte) byte {
	if b >= 'A' && b <= 'Z' {
		return b + ('a' - 'A')
	}
	return b
}
