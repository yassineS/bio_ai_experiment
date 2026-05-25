// mutations.go: subcommands that mutate sequences in place.
//
// This file contains the Mutfa and Randbase functions implementing the seqtk
// mutfa and randbase subcommands, plus the IUPAC ambiguity-code lookup table
// used by Randbase.

package seqtk

import (
	"bufio"
	"fmt"
	"io"
	"math/rand"
	"os"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
)

// Mutation describes a single point mutation on the forward strand.
// Pos is the 0-based position on the reference sequence; Base is the
// substitute base (a single ASCII byte such as 'A'/'C'/'G'/'T', case-preserved
// from the mutation file).
type Mutation struct {
	Pos  int
	Base byte
}

// MutationSet maps a sequence/chromosome name to the list of mutations to
// apply on that sequence, in input order.
type MutationSet map[string][]Mutation

// ParseMutfile parses a TSV mutation file with at least three whitespace- or
// tab-separated columns per non-empty, non-comment line: chrom, 1-based pos,
// substitute base. Additional trailing columns are ignored to remain
// compatible with the four-column upstream seqtk format (chrom, pos, ref, alt)
// — in that case the third upstream column ("ref") is read as the new base
// only when there are exactly three columns; with four or more columns the
// fourth column is treated as the new base.
//
// Lines whose new-base field is not a single alphabetic ASCII byte are
// skipped with a warning to stderr. Positions <= 0 are skipped with a warning.
// Header / comment lines starting with '#' are ignored.
func ParseMutfile(r io.Reader) (MutationSet, error) {
	out := make(MutationSet)
	scanner := bufio.NewScanner(r)
	buf := make([]byte, 0, 64*1024)
	scanner.Buffer(buf, 10*1024*1024)

	lineNo := 0
	for scanner.Scan() {
		lineNo++
		line := strings.TrimRight(scanner.Text(), "\r\n")
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) < 3 {
			fmt.Fprintf(os.Stderr, "[seqtk mutfa] warning: line %d: expected at least 3 fields, got %d; skipped\n", lineNo, len(fields))
			continue
		}
		chrom := fields[0]
		pos, err := strconv.Atoi(fields[1])
		if err != nil {
			fmt.Fprintf(os.Stderr, "[seqtk mutfa] warning: line %d: invalid position %q; skipped\n", lineNo, fields[1])
			continue
		}
		if pos <= 0 {
			fmt.Fprintf(os.Stderr, "[seqtk mutfa] warning: line %d: position %d is not >= 1; skipped\n", lineNo, pos)
			continue
		}
		// Pick the new-base field: column 3 for a 3-column file, column 4 for >=4 columns
		// (matching upstream seqtk's "chrom  pos  ref  alt" format).
		var baseField string
		if len(fields) >= 4 {
			baseField = fields[3]
		} else {
			baseField = fields[2]
		}
		if len(baseField) != 1 || !isAlpha(baseField[0]) {
			fmt.Fprintf(os.Stderr, "[seqtk mutfa] warning: line %d: invalid base %q; skipped\n", lineNo, baseField)
			continue
		}
		out[chrom] = append(out[chrom], Mutation{Pos: pos - 1, Base: baseField[0]})
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}
	return out, nil
}

// isAlpha reports whether b is an ASCII letter.
func isAlpha(b byte) bool {
	return (b >= 'A' && b <= 'Z') || (b >= 'a' && b <= 'z')
}

// Mutfa applies point mutations to a FASTA stream. mutR provides the mutation
// file (TSV: chrom, 1-based pos, [optional ref,] new base). Output is FASTA
// written to w with the same line-width layout as the input — i.e. the
// physical line breaks of each input sequence are preserved on output.
//
// For each input record, all mutations targeting that record's name are
// applied to its sequence on the forward strand. Out-of-range positions and
// names not present in the input are skipped with a warning on stderr; they
// are not fatal.
func Mutfa(in io.Reader, mutR io.Reader, w io.Writer) error {
	muts, err := ParseMutfile(mutR)
	if err != nil {
		return err
	}

	// Track which chromosomes from the mutation file were actually seen so
	// we can warn about ones that were never found.
	seen := make(map[string]bool, len(muts))

	br, _ := peekIsFastq(in) // upstream seqtk mutfa assumes FASTA; the peek just normalises whitespace
	scanner := bufio.NewScanner(br)
	buf := make([]byte, 0, 64*1024)
	scanner.Buffer(buf, 10*1024*1024)
	bw := bufio.NewWriter(w)

	// Per-record state: when we see a '>' header line, we flush any pending
	// state, write the header verbatim, and start tracking position in the
	// new sequence. Each sequence line is then written back out as-is, with
	// individual bases substituted where mutations target their offsets.
	var (
		curName string
		curMuts []Mutation
		curPos  int  // 0-based offset of the next base to be written
		inSeq   bool // true if we've seen the header for the current record
	)

	applyAndWriteLine := func(line []byte) error {
		if !inSeq {
			return nil
		}
		// We may have to mutate one or more bytes in this line. Make a
		// modifiable copy only if needed.
		var mutated []byte
		lineEnd := curPos + len(line)
		for _, m := range curMuts {
			if m.Pos < curPos || m.Pos >= lineEnd {
				continue
			}
			if mutated == nil {
				mutated = make([]byte, len(line))
				copy(mutated, line)
			}
			mutated[m.Pos-curPos] = m.Base
		}
		out := line
		if mutated != nil {
			out = mutated
		}
		if _, err := bw.Write(out); err != nil {
			return err
		}
		if err := bw.WriteByte('\n'); err != nil {
			return err
		}
		curPos = lineEnd
		return nil
	}

	finishRecord := func() {
		if !inSeq {
			return
		}
		// Warn about out-of-range positions for the just-finished record.
		for _, m := range curMuts {
			if m.Pos >= curPos {
				fmt.Fprintf(os.Stderr, "[seqtk mutfa] warning: position %d on %q is past sequence end (%d); skipped\n", m.Pos+1, curName, curPos)
			}
		}
	}

	for scanner.Scan() {
		raw := scanner.Bytes()
		// scanner.Bytes() is only valid until the next Scan; copy when storing.
		if len(raw) > 0 && raw[0] == '>' {
			// New header: emit any trailing warnings for the previous record.
			finishRecord()

			header := string(raw[1:])
			fields := strings.Fields(header)
			name := ""
			if len(fields) > 0 {
				name = fields[0]
			}
			curName = name
			curMuts = muts[name]
			if curMuts != nil {
				seen[name] = true
			}
			curPos = 0
			inSeq = true

			// Write header verbatim.
			if _, err := bw.Write(raw); err != nil {
				return err
			}
			if err := bw.WriteByte('\n'); err != nil {
				return err
			}
			continue
		}
		// Sequence line (possibly empty — keep blank lines verbatim if not in a record).
		if !inSeq {
			// Stray content before any header — pass through to preserve content.
			if _, err := bw.Write(raw); err != nil {
				return err
			}
			if err := bw.WriteByte('\n'); err != nil {
				return err
			}
			continue
		}
		if err := applyAndWriteLine(raw); err != nil {
			return err
		}
	}
	if err := scanner.Err(); err != nil {
		return err
	}

	finishRecord()

	// Warn about chromosomes from the mutation file that were never seen.
	for name := range muts {
		if !seen[name] {
			fmt.Fprintf(os.Stderr, "[seqtk mutfa] warning: chromosome %q not found in input FASTA; %d mutation(s) skipped\n", name, len(muts[name]))
		}
	}

	return bw.Flush()
}

// iupacExpansions maps each IUPAC ambiguity code to the set of unambiguous
// bases it represents. Keys are uppercase; case is handled at lookup time.
//
// Note that only the two-base codes (R/Y/S/W/K/M) are listed here on purpose:
// upstream `seqtk randbase` only randomises positions whose IUPAC bit-count
// is exactly 2, leaving the three-base codes (B/D/H/V) and the four-base
// code (N) untouched. Randomising those would silently diverge from
// upstream — see `int stk_randbase` in reference_code/seqtk/seqtk.c
// (the `if (a == 2)` guard).
var iupacExpansions = map[byte][]byte{
	'R': {'A', 'G'},
	'Y': {'C', 'T'},
	'S': {'G', 'C'},
	'W': {'A', 'T'},
	'K': {'G', 'T'},
	'M': {'A', 'C'},
}

// pickIUPAC returns one of the bases that the given IUPAC code expands to,
// chosen uniformly at random via rng. The case of the returned base matches
// the case of the input byte: lowercase code => lowercase result. If the byte
// is not an IUPAC ambiguity code (or is an unambiguous base), it is returned
// unchanged.
func pickIUPAC(b byte, rng *rand.Rand) byte {
	upper := b
	lower := false
	if upper >= 'a' && upper <= 'z' {
		upper -= 'a' - 'A'
		lower = true
	}
	choices, ok := iupacExpansions[upper]
	if !ok {
		return b
	}
	out := choices[rng.Intn(len(choices))]
	if lower {
		out += 'a' - 'A'
	}
	return out
}

// Randbase replaces every two-base IUPAC ambiguity code (R/Y/S/W/K/M)
// in a FASTA stream with one of its two underlying bases, preserving
// case. Three- and four-base codes (B/D/H/V/N) pass through unchanged
// — this matches upstream `stk_randbase`'s `if (a == 2)` gate (only
// IUPAC codes whose bit-count is exactly 2 are randomised).
//
// Output layout mirrors upstream byte-for-byte: each record's sequence
// is wrapped to 60 columns regardless of the input's physical line
// breaks (upstream uses `if (i%60 == 0) putchar('\n')` inside
// `stk_randbase`).
//
// The `seed` parameter is accepted for API stability with earlier
// versions of the port; upstream's randbase has no -s flag and always
// uses glibc's default drand48 state (X0 = 0), so for byte parity we
// ignore the seed. Use a separate dedicated PRNG path if seedable
// behaviour is required.
func Randbase(in io.Reader, w io.Writer, seed int64) error {
	_ = seed // intentionally unused; see doc comment.
	d := &drand48State{}

	br, _ := peekIsFastq(in)
	reader := fasta.NewReader(br)
	bw := bufio.NewWriter(w)

	for {
		rec, err := reader.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}
		if _, err := bw.WriteString(">"); err != nil {
			return err
		}
		if _, err := bw.WriteString(rec.Description); err != nil {
			return err
		}
		// Match upstream's output formatting: the header is written
		// without a trailing newline; the per-base loop then runs
		// `if (i%60 == 0) putchar('\n')` which (at i=0) emits the
		// header's newline AND starts the first wrapped sequence line.
		// After the loop a final `putchar('\n')` closes the last
		// (possibly short) sequence line.
		seq := rec.Sequence
		for i, b := range seq {
			if i%60 == 0 {
				if err := bw.WriteByte('\n'); err != nil {
					return err
				}
			}
			if err := bw.WriteByte(pickIUPACDrand48(b, d)); err != nil {
				return err
			}
		}
		if err := bw.WriteByte('\n'); err != nil {
			return err
		}
	}
	return bw.Flush()
}

// pickIUPACDrand48 replicates upstream `stk_randbase`'s per-byte
// substitution rule (seqtk.c:540-559): for 2-base IUPAC codes draw one
// glibc-drand48 sample, pick "ACGT"[j] (or "acgt"[j]) where j is the
// index of the m-th set bit in the IUPAC bitmask. m = (drand48() <
// 0.5) ? 1 : 0. Returns the input byte unchanged for non-2-base codes.
func pickIUPACDrand48(b byte, d *drand48State) byte {
	// seq_nt16_table compresses to a 4-bit mask: bit 0=A, 1=C, 2=G,
	// 3=T. We compute the mask for the upper-case letter then encode
	// case at the very end.
	upper := b
	lower := false
	if upper >= 'a' && upper <= 'z' {
		upper -= 'a' - 'A'
		lower = true
	}
	var mask byte
	switch upper {
	case 'A':
		mask = 1
	case 'C':
		mask = 2
	case 'G':
		mask = 4
	case 'T':
		mask = 8
	case 'R':
		mask = 1 | 4 // A|G
	case 'Y':
		mask = 2 | 8 // C|T
	case 'S':
		mask = 2 | 4 // C|G
	case 'W':
		mask = 1 | 8 // A|T
	case 'K':
		mask = 4 | 8 // G|T
	case 'M':
		mask = 1 | 2 // A|C
	default:
		return b
	}
	// Population count: only randomise when exactly 2 bits set
	// (upstream's `if (a == 2)` gate).
	if bitCount4(mask) != 2 {
		return b
	}
	m := 0
	if d.next() < 0.5 {
		m = 1
	}
	k := 0
	j := 0
	for j = 0; j < 4; j++ {
		if mask&(1<<uint(j)) == 0 {
			continue
		}
		if k == m {
			break
		}
		k++
	}
	out := "ACGT"[j]
	if lower {
		out += 'a' - 'A'
	}
	return out
}

// bitCount4 returns the number of set bits in the low nibble of x.
// Used to mirror upstream's `bitcnt_table` lookup in stk_randbase.
func bitCount4(x byte) int {
	n := 0
	for i := 0; i < 4; i++ {
		if x&(1<<uint(i)) != 0 {
			n++
		}
	}
	return n
}
