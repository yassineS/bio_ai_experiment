// comp.go: per-record nucleotide composition statistics that mirror upstream
// "seqtk comp" output.
//
// Upstream emits one tab-separated row per input record with the columns:
//
//	chr, length, #A, #C, #G, #T, #2, #3, #4, #CpG, #tv, #ts, #CpG-ts
//
// where #2 / #3 / #4 are counts of two-/three-/four-base IUPAC ambiguity
// codes, #CpG counts CpG-overlapping positions, and #tv / #ts / #CpG-ts
// count transversions / transitions / CpG transitions versus the preceding
// base in the sequence.
//
// This file re-uses upstream's seq_nt16 + bitcnt tables verbatim so the
// numbers match upstream byte-for-byte.

package seqtk

import (
	"bufio"
	"fmt"
	"io"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fastq"
)

// seqNT16Table maps an ASCII byte to its IUPAC 4-bit index used by upstream
// seqtk. Indexing is by byte value (0-255). Unknown / non-IUPAC bytes map to
// 15 (N), matching upstream's seqtk.c.
var seqNT16Table = func() [256]byte {
	var t [256]byte
	for i := range t {
		t[i] = 15
	}
	// Uppercase letters.
	t['A'] = 1
	t['B'] = 14
	t['C'] = 2
	t['D'] = 13
	t['G'] = 4
	t['H'] = 11
	t['K'] = 12
	t['M'] = 3
	t['N'] = 15
	t['R'] = 5
	t['S'] = 6
	t['T'] = 8
	t['U'] = 15 // upstream seq_nt16_table[85] = 15 (U is treated as unknown)
	t['V'] = 7
	t['W'] = 9
	t['X'] = 0
	t['Y'] = 10
	// Lowercase letters.
	t['a'] = 1
	t['b'] = 14
	t['c'] = 2
	t['d'] = 13
	t['g'] = 4
	t['h'] = 11
	t['k'] = 12
	t['m'] = 3
	t['n'] = 15
	t['r'] = 5
	t['s'] = 6
	t['t'] = 8
	t['u'] = 15 // upstream seq_nt16_table[117] = 15 (U is treated as unknown)
	t['v'] = 7
	t['w'] = 9
	t['x'] = 0
	t['y'] = 10
	return t
}()

// seqNT16To4Table folds the 16-state IUPAC code down to a 4-state index used
// by upstream when counting unambiguous bases. The 4 returned for ambiguous
// codes is unused at the call site (guarded by bitcnt == 1).
var seqNT16To4Table = [16]byte{4, 0, 1, 4, 2, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4}

// bitCntTable maps a 4-bit IUPAC code to the number of bases it represents.
// 1-bit codes are unambiguous (A/C/G/T), 2-bit codes are R/Y/S/W/K/M, 3-bit
// codes are B/D/H/V, and 4 maps to N (X collides with N for the bitcnt of 4).
var bitCntTable = [16]byte{4, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4}

// Comp writes per-record nucleotide-composition rows for every record in the
// input stream (auto-detected as FASTA or FASTQ), one row per record, in
// upstream "seqtk comp" format: tab-separated
//
//	name\tlen\t#A\t#C\t#G\t#T\t#2\t#3\t#4\t#CpG\t#tv\t#ts\t#CpG-ts
//
// Counts are case-insensitive. Unknown bytes are treated as N.
//
// This is a 1:1 port of upstream's stk_comp inner loop (in
// reference_code/seqtk/seqtk.c) — see the package-level comment above for
// the field definitions.
func Comp(in io.Reader, w io.Writer) error {
	return CompWithRegions(in, w, nil)
}

// CompWithRegions is Comp with optional region restriction (upstream `comp -r
// FILE`). When regions is non-nil, only records whose name appears in the
// region hash are reported, and for each listed [beg,end) interval one output
// row is emitted with the columns "name\tbeg\tend" followed by the 11 counts
// computed over that interval (matching upstream's stk_comp). When regions is
// nil the whole-record path is used and the leading columns are "name\tlen".
func CompWithRegions(in io.Reader, w io.Writer, regions *regHash) error {
	br, isFastq := peekIsFastq(in)
	bw := bufio.NewWriter(w)

	// countInterval computes the 11 composition counts over seq[beg:end],
	// reading neighbour bytes outside the interval for CpG context exactly as
	// upstream does (seq[beg-1] lookback, seq[i+1] lookahead).
	countInterval := func(seq []byte, beg, end int) [11]int64 {
		var cnt [11]int64
		var lb int = -1
		if beg > 0 {
			lb = int(seqNT16Table[seq[beg-1]])
		}
		for i := beg; i < end; i++ {
			b := int(seqNT16Table[seq[i]])
			c := int(bitCntTable[b])
			nb := -1
			if i+1 < len(seq) {
				nb = int(seqNT16Table[seq[i+1]])
			}
			isCpG := false
			if b == 2 || b == 10 { // C or Y
				if nb == 4 || nb == 5 { // G or R
					isCpG = true
				}
			} else if b == 4 || b == 5 { // G or R
				if lb == 2 || lb == 10 { // previous C or Y
					isCpG = true
				}
			}
			if c > 1 {
				cnt[c+2]++ // #2 (idx 4), #3 (idx 5), #4 (idx 6)
			}
			if c == 1 {
				cnt[seqNT16To4Table[b]]++ // #A/#C/#G/#T
			}
			if b == 10 || b == 5 {
				cnt[9]++ // transition (Y, R)
			} else if c == 2 {
				cnt[8]++ // transversion
			}
			if isCpG {
				cnt[7]++ // CpG
				if b == 10 || b == 5 {
					cnt[10]++ // CpG transition
				}
			}
			lb = b
		}
		return cnt
	}

	writeCounts := func(cnt [11]int64) error {
		for _, v := range cnt {
			if _, err := fmt.Fprintf(bw, "\t%d", v); err != nil {
				return err
			}
		}
		return bw.WriteByte('\n')
	}

	emit := func(name string, seq []byte) error {
		if regions != nil {
			intervals, ok := regions.m[name]
			if !ok {
				return nil
			}
			for _, iv := range intervals {
				beg, end := int(iv[0]), int(iv[1])
				if beg < 0 {
					beg = 0
				}
				if end > len(seq) {
					end = len(seq)
				}
				cnt := countInterval(seq, beg, end)
				if _, err := fmt.Fprintf(bw, "%s\t%d\t%d", name, iv[0], iv[1]); err != nil {
					return err
				}
				if err := writeCounts(cnt); err != nil {
					return err
				}
			}
			return nil
		}
		cnt := countInterval(seq, 0, len(seq))
		if _, err := fmt.Fprintf(bw, "%s\t%d", name, len(seq)); err != nil {
			return err
		}
		return writeCounts(cnt)
	}

	if isFastq {
		r := fastq.NewReader(br, fastq.Phred33)
		for {
			rec, err := r.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}
			if err := emit(rec.ID, rec.Sequence); err != nil {
				return err
			}
		}
	} else {
		r := fasta.NewReader(br)
		for {
			rec, err := r.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}
			if err := emit(rec.ID, rec.Sequence); err != nil {
				return err
			}
		}
	}
	return bw.Flush()
}
