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

// drand48 is a faithful port of the POSIX/glibc drand48 random number
// generator used by upstream `seqtk randbase`. It is the classic 48-bit
// linear congruential generator
//
//	X_{n+1} = (a * X_n + c) mod 2^48,  a = 0x5DEECE66D, c = 0xB
//
// with the glibc default initial state X_0 = 0 (i.e. no srand48 call). The
// floating-point value returned by next is X / 2^48 in [0, 1), bit-for-bit
// identical to glibc's drand48 for the default seed, which is what
// stk_randbase relies on.
type drand48 struct {
	x uint64
}

const (
	drand48A    = 0x5DEECE66D
	drand48C    = 0xB
	drand48Mask = (uint64(1) << 48) - 1
)

// next advances the generator and returns the next value in [0, 1), matching
// glibc's drand48().
func (d *drand48) next() float64 {
	d.x = (drand48A*d.x + drand48C) & drand48Mask
	return float64(d.x) / 281474976710656.0 // == X / 2^48
}

// iupacBits maps each IUPAC code to its 4-bit set (A=1, C=2, G=4, T=8),
// mirroring upstream's seq_nt16_table. Codes not present here map to 15 (N),
// which has a bit count of 4 and is therefore never randomised.
var iupacBits = func() [256]uint8 {
	var t [256]uint8
	for i := range t {
		t[i] = 15
	}
	set := func(b byte, v uint8) { t[b] = v; t[b|0x20] = v } // upper + lower
	set('A', 1)
	set('C', 2)
	set('G', 4)
	set('T', 8)
	set('M', 1|2)
	set('R', 1|4)
	set('W', 1|8)
	set('S', 2|4)
	set('Y', 2|8)
	set('K', 4|8)
	set('V', 1|2|4)
	set('H', 1|2|8)
	set('D', 1|4|8)
	set('B', 2|4|8)
	set('N', 15)
	return t
}()

// bitCount4 returns the number of set bits in the low nibble of c, mirroring
// upstream's bitcnt_table.
func bitCount4(c uint8) int {
	n := 0
	for j := uint(0); j < 4; j++ {
		if c&(1<<j) != 0 {
			n++
		}
	}
	return n
}

// RandbaseUpstream is a byte-for-byte port of upstream's stk_randbase. For
// every base whose IUPAC code has exactly two possible nucleotides
// (R/Y/S/W/K/M), it draws one of the two using the glibc drand48 generator
// (m = (drand48() < 0.5)) and replaces the base, preserving case; all other
// bases — including the three-base codes (B/D/H/V), the four-base code (N) and
// the unambiguous bases — pass through unchanged.
//
// Output exactly reproduces stk_randbase's layout: the header is written as
// ">name" (the comment is dropped), the sequence is wrapped at 60 columns
// (with a newline emitted before column 0, i.e. right after the header), and a
// trailing newline closes each record. The drand48 stream starts from glibc's
// default seed (state 0), so the output is deterministic and matches upstream
// regardless of any CLI seed.
func RandbaseUpstream(in io.Reader, w io.Writer) error {
	br, _ := peekIsFastq(in)
	r := fasta.NewReader(br)
	bw := bufio.NewWriter(w)
	d := &drand48{}

	for {
		rec, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}
		if _, err := bw.WriteString(">"); err != nil {
			return err
		}
		if _, err := bw.WriteString(rec.ID); err != nil {
			return err
		}
		seq := rec.Sequence
		for i := 0; i < len(seq); i++ {
			b := seq[i]
			c := iupacBits[b]
			if bitCount4(c) == 2 {
				m := 0
				if d.next() < 0.5 {
					m = 1
				}
				// Walk the set bits in order; pick the m-th one.
				k, j := 0, 0
				for j = 0; j < 4; j++ {
					if c&(1<<uint(j)) == 0 {
						continue
					}
					if k == m {
						break
					}
					k++
				}
				if b >= 'a' && b <= 'z' {
					b = "acgt"[j]
				} else {
					b = "ACGT"[j]
				}
			}
			if i%60 == 0 {
				if err := bw.WriteByte('\n'); err != nil {
					return err
				}
			}
			if err := bw.WriteByte(b); err != nil {
				return err
			}
		}
		if err := bw.WriteByte('\n'); err != nil {
			return err
		}
	}
	return bw.Flush()
}

// Randbase replaces every IUPAC ambiguity base (R/Y/S/W/K/M/B/D/H/V/N) in a
// FASTA stream with a uniform random sample from its expansion, preserving
// case. The output is written to w as FASTA, preserving the input's line
// widths (physical line breaks). Non-ambiguity bases are passed through
// unchanged.
//
// If seed != 0 the random source is seeded deterministically with seed; if
// seed == 0 a time-based seed is used (caller's responsibility).
//
// This is the seeded extension used when the CLI is given -s. The default
// (seedless) randbase path uses RandbaseUpstream, which matches upstream
// byte-for-byte.
func Randbase(in io.Reader, w io.Writer, seed int64) error {
	rng := rand.New(rand.NewSource(seed))

	br, _ := peekIsFastq(in)
	scanner := bufio.NewScanner(br)
	buf := make([]byte, 0, 64*1024)
	scanner.Buffer(buf, 10*1024*1024)
	bw := bufio.NewWriter(w)

	for scanner.Scan() {
		raw := scanner.Bytes()
		if len(raw) > 0 && raw[0] == '>' {
			if _, err := bw.Write(raw); err != nil {
				return err
			}
			if err := bw.WriteByte('\n'); err != nil {
				return err
			}
			continue
		}
		// Sequence line: walk byte-by-byte, substituting IUPAC codes.
		out := make([]byte, len(raw))
		for i, b := range raw {
			out[i] = pickIUPAC(b, rng)
		}
		if _, err := bw.Write(out); err != nil {
			return err
		}
		if err := bw.WriteByte('\n'); err != nil {
			return err
		}
	}
	if err := scanner.Err(); err != nil {
		return err
	}
	return bw.Flush()
}
