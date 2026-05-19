// scanner.go: shared sequence-scanning helpers used by the FASTA "interval"
// subcommands (`gap`, `gc`). Each subcommand walks the sequence byte stream
// and emits BED-style intervals where some predicate holds; the helpers in
// this file factor out the common iteration and BED-writing plumbing so the
// individual subcommands focus on their predicate / scoring logic.
//
// All helpers assume FASTA input — that matches upstream `seqtk gap` and
// `seqtk gc`, both of which call `kseq_read` and only ever look at the
// sequence bytes (no quality data).

package seqtk

import (
	"bufio"
	"fmt"
	"io"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
)

// seqNT6Table folds an ASCII byte to upstream seqtk's 6-state nucleotide
// alphabet (the `seq_nt6_table` array in reference_code/seqtk/seqtk.c).
// The mapping is: A/a -> 1, C/c -> 2, G/g -> 3, T/t -> 4, 0-byte -> 0,
// everything else (including N, U/u (RNA uracil), IUPAC codes and any
// other non-DNA byte) -> 5. `seqtk gap` uses this table to detect "gap"
// runs: a gap is any maximal run of bytes that map to 5, regardless of
// whether they are literal N's.
var seqNT6Table = func() [256]byte {
	var t [256]byte
	for i := range t {
		t[i] = 5
	}
	t[0] = 0
	t['A'], t['a'] = 1, 1
	t['C'], t['c'] = 2, 2
	t['G'], t['g'] = 3, 3
	t['T'], t['t'] = 4, 4
	// U/u (uracil) stays at 5 to match upstream seq_nt6_table[85] == 5
	// (seqtk.c:208). Previously the port mapped U/u -> 4, which silently
	// hid uracil from `seqtk gap`; reviewer-caught regression on PR #112.
	return t
}()

// IsGapByte reports whether b is part of a "gap" under upstream seqtk's
// definition (anything that is not A, C, G or T, case-insensitive — note
// U/u is also a gap, matching upstream).
func IsGapByte(b byte) bool { return seqNT6Table[b] == 5 }

// scanFASTA reads every record from r and calls emit for each one. It is a
// thin wrapper used by the gap / gc subcommands to keep their bodies focused
// on the per-sequence predicate. Errors from the FASTA reader are surfaced;
// emit may return an error to short-circuit the scan.
func scanFASTA(r io.Reader, emit func(name string, seq []byte) error) error {
	fr := fasta.NewReader(r)
	for {
		rec, err := fr.Read()
		if err == io.EOF {
			return nil
		}
		if err != nil {
			return err
		}
		if err := emit(rec.ID, rec.Sequence); err != nil {
			return err
		}
	}
}

// writeBED3 writes one chrom\tstart\tend\n BED record to w. Coordinates are
// 0-based half-open, matching upstream `seqtk gap` output.
func writeBED3(w *bufio.Writer, chrom string, start, end int) error {
	_, err := fmt.Fprintf(w, "%s\t%d\t%d\n", chrom, start, end)
	return err
}

// writeBED4Int writes one chrom\tstart\tend\textra\n BED-like record to w
// where the extra column is an integer. Used by `seqtk gc` to emit its
// hit-count column.
func writeBED4Int(w *bufio.Writer, chrom string, start, end, extra int) error {
	_, err := fmt.Fprintf(w, "%s\t%d\t%d\t%d\n", chrom, start, end, extra)
	return err
}
