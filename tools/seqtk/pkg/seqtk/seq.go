// seq.go: a faithful, byte-for-byte port of upstream `seqtk seq`
// (the `stk_seq` function in reference_code/seqtk/seqtk.c, version
// 1.5-r133). The transformation order below mirrors stk_seq exactly:
// length filter -> random sampling -> odd/even selection -> -S strip
// whitespace -> quality masking -> -U/-x case conversion -> -A/-F
// quality handling -> -C drop comments -> -M region masking ->
// -r/-R reverse complement -> -V quality shift -> -N drop ambiguous
// -> print. Getting that order right is what makes the output match
// the genuine binary.
//
// Rather than reuse the shared pkg/htsgo/fasta and fastq readers (which
// trim whitespace inside sequence lines and split name/comment with
// strings.Fields), this file ships its own kseq-compatible reader so
// that -S (strip whitespace) and -C (drop comments) behave exactly like
// upstream's kseq.h.

package seqtk

import (
	"bufio"
	"io"
)

// SeqOptions mirrors the option bundle parsed by stk_seq. Zero values
// match upstream defaults (see the stk_seq initialisers in seqtk.c:1385).
type SeqOptions struct {
	ForceFasta  bool    // -A / -a: force FASTA output (discard quality)
	DropComment bool    // -C: drop comments on header lines
	RevComp     bool    // -r: reverse complement
	BothStrands bool    // -R: output both forward and reverse complement
	MaskComp    bool    // -c: mask the complement of -M regions
	Odd         bool    // -1: output the 2n-1 (odd) reads only
	Even        bool    // -2: output the 2n (even) reads only
	ShiftQual   bool    // -V: shift quality by (-Q) - 33
	DropAmbig   bool    // -N: drop sequences containing ambiguous bases
	Uppercase   bool    // -U: convert all bases to uppercase
	StripSpace  bool    // -S: strip white space in sequences
	LowerToMask bool    // -x: convert lowercase bases to the -n char
	QualThres   int     // -q: mask bases with quality lower than INT [0]
	MaxQual     int     // -X: mask bases with quality higher than INT [255]
	MaskChar    byte    // -n: masked bases converted to CHAR; 0 => lowercase
	LineLen     int     // -l: residues per line; 0 => no wrap (2^32-1)
	QualShift   int     // -Q: quality shift; ASCII-INT gives base quality [33]
	MinLen      int     // -L: drop sequences shorter than INT [0]
	FakeQual    int     // -F: fake FASTQ quality char; <0 means unset
	Seed        int64   // -s: random seed (effective with -f) [11]
	Frac        float64 // -f: sample FLOAT fraction of sequences [1]
	MaskFile    string  // -M: mask regions in BED or name-list FILE
	maskRegions regHash
}

// NewSeqOptions returns a SeqOptions populated with the same defaults as
// upstream stk_seq.
func NewSeqOptions() SeqOptions {
	return SeqOptions{
		MaxQual:   255,
		QualShift: 33,
		FakeQual:  -1,
		Seed:      11,
		Frac:      1.0,
	}
}

// compTab is upstream's comp_tab (seqtk.c:226): a 128-entry table mapping
// each ASCII byte to its IUPAC complement, preserving case. Bytes >= 128
// map to themselves (matching the C cast (int)seq->seq.s[i] on signed
// char, but kseq only stores bytes; we treat the high range as identity).
var compTab = func() [256]byte {
	var t [256]byte
	for i := range t {
		t[i] = byte(i)
	}
	// Lower 128 from upstream comp_tab.
	src := []byte{
		0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
		16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
		32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47,
		48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63,
		64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
		'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z', 91, 92, 93, 94, 95,
		64, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
		'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127,
	}
	copy(t[:], src)
	return t
}()

// Note: seqNT16Table (comp.go) and seqNT16To4Table (comp.go) are the
// upstream seq_nt16_table / seq_nt16to4_table; -N reuses them.

func isSpaceByte(b byte) bool {
	switch b {
	case ' ', '\t', '\n', '\v', '\f', '\r':
		return true
	}
	return false
}

// toLowerByte is defined in famask.go.

func toUpperByte(b byte) byte {
	if b >= 'a' && b <= 'z' {
		return b - 32
	}
	return b
}

func isLowerByte(b byte) bool { return b >= 'a' && b <= 'z' }

// kseqRecord is the data kseq_read produces for one record.
type kseqRecord struct {
	name    []byte
	comment []byte
	seq     []byte
	qual    []byte // empty for FASTA
}

// kseqReader is a faithful port of kseq.h's kseq_read for the subset of
// behaviour seqtk relies on: it splits the header into name (up to the
// first whitespace) and comment (the rest of the line), concatenates
// sequence lines dropping only newlines, and reads quality for FASTQ.
type kseqReader struct {
	br       *bufio.Reader
	lastChar byte // 0 means "seek next header"
	started  bool
}

func newKseqReader(r io.Reader) *kseqReader {
	return &kseqReader{br: bufio.NewReaderSize(r, 64*1024)}
}

// readByte returns the next byte and a bool that is false at EOF.
func (k *kseqReader) readByte() (byte, bool) {
	b, err := k.br.ReadByte()
	if err != nil {
		return 0, false
	}
	return b, true
}

// read fetches the next record. ok is false at end of input.
func (k *kseqReader) read() (rec kseqRecord, ok bool, err error) {
	if !k.started || k.lastChar == 0 {
		// Jump to the next header line.
		for {
			c, more := k.readByte()
			if !more {
				return rec, false, nil
			}
			if c == '>' || c == '@' {
				k.lastChar = c
				break
			}
		}
	}
	k.started = true

	// Read the name: bytes up to (but not including) the first whitespace,
	// matching ks_getuntil with KS_SEP_SPACE (delimiter 0). The terminating
	// delimiter is consumed; record whether it was a newline.
	var name []byte
	endedAtNewline := false
	for {
		c, more := k.readByte()
		if !more {
			endedAtNewline = true
			break
		}
		if isSpaceByte(c) {
			endedAtNewline = c == '\n'
			break
		}
		name = append(name, c)
	}

	// Read the comment (rest of the line) if the name did not end at a
	// newline. KS_SEP_LINE strips a trailing '\r'.
	var comment []byte
	if !endedAtNewline {
		for {
			c, more := k.readByte()
			if !more || c == '\n' {
				break
			}
			comment = append(comment, c)
		}
		if n := len(comment); n > 0 && comment[n-1] == '\r' {
			comment = comment[:n-1]
		}
	}

	// Read sequence lines until the next header / '+' separator / EOF.
	var seq []byte
	stopped := byte(0) // the header / '+' char that ended the loop, 0 at EOF
	for {
		c, more := k.readByte()
		if !more {
			break
		}
		if c == '>' || c == '+' || c == '@' {
			stopped = c
			break
		}
		if c == '\n' {
			continue
		}
		// First char of a sequence line, then the rest until newline.
		seq = append(seq, c)
		for {
			c2, more2 := k.readByte()
			if !more2 || c2 == '\n' {
				break
			}
			seq = append(seq, c2)
		}
	}

	rec.name = name
	rec.comment = comment
	rec.seq = seq

	if stopped != '+' {
		// FASTA record (or EOF). stopped is the next header char, or 0 at
		// EOF; either way it becomes lastChar so the next read() resumes
		// correctly (0 makes it seek the next header).
		k.lastChar = stopped
		return rec, true, nil
	}

	// FASTQ: skip the rest of the '+' line.
	for {
		c, more := k.readByte()
		if !more {
			// No quality string at EOF: upstream returns -2 (error). We
			// surface it but still report the record is unavailable.
			return rec, false, errNoQuality
		}
		if c == '\n' {
			break
		}
	}

	// Read quality until we have len(seq) bytes (dropping newlines).
	var qual []byte
	for len(qual) < len(seq) {
		c, more := k.readByte()
		if !more {
			break
		}
		if c == '\n' {
			continue
		}
		qual = append(qual, c)
		for len(qual) < len(seq) {
			c2, more2 := k.readByte()
			if !more2 || c2 == '\n' {
				break
			}
			qual = append(qual, c2)
		}
	}
	rec.qual = qual
	k.lastChar = 0 // not yet at the next header line
	return rec, true, nil
}

// errNoQuality mirrors upstream kseq_read returning -2 when a FASTQ
// record lacks a quality string.
var errNoQuality = errSeq("malformed FASTQ: missing quality string")

type errSeq string

func (e errSeq) Error() string { return string(e) }

// Seq processes input according to opts and writes the transformed
// records to output, matching `seqtk seq`. The transformation order is
// identical to stk_seq.
func Seq(input io.Reader, output io.Writer, opts SeqOptions) error {
	w := bufio.NewWriterSize(output, 64*1024)

	// Resolve line length: 0 means no wrap (upstream sets it to UINT_MAX).
	lineLen := opts.LineLen

	var kr *krand
	if opts.Frac < 1.0 {
		kr = newKrand(uint64(opts.Seed))
	}

	r := newKseqReader(input)
	var nSeqs int64
	for {
		rec, ok, err := r.read()
		if err != nil {
			return err
		}
		if !ok {
			break
		}
		nSeqs++

		// 1. Length filter (before random sampling, matching upstream).
		if len(rec.seq) < opts.MinLen {
			continue
		}
		// 2. Random sampling with -f.
		if opts.Frac < 1.0 && kr.drand() >= opts.Frac {
			continue
		}
		// 3. Odd/even (-1/-2) selection.
		if opts.Odd || opts.Even {
			if opts.Odd && nSeqs&1 == 0 {
				continue
			}
			if opts.Even && nSeqs&1 == 1 {
				continue
			}
		}
		// 4. -S: squeeze out white space (qual reindexed by seq positions).
		if opts.StripSpace {
			if len(rec.qual) > 0 {
				k := 0
				for i := 0; i < len(rec.seq); i++ {
					if !isSpaceByte(rec.seq[i]) {
						rec.qual[k] = rec.qual[i]
						k++
					}
				}
				rec.qual = rec.qual[:k]
			}
			k := 0
			for i := 0; i < len(rec.seq); i++ {
				if !isSpaceByte(rec.seq[i]) {
					rec.seq[k] = rec.seq[i]
					k++
				}
			}
			rec.seq = rec.seq[:k]
		}
		// 5. Quality masking (-q/-X). qual_thres in stk_seq is
		// opts.QualThres + qual_shift; the guard is qual_thres > qual_shift.
		qualThres := opts.QualThres + opts.QualShift
		if len(rec.qual) > 0 && qualThres > opts.QualShift {
			for i := 0; i < len(rec.seq); i++ {
				q := int(rec.qual[i])
				if q < qualThres || q > opts.MaxQual {
					if opts.MaskChar != 0 {
						rec.seq[i] = opts.MaskChar
					} else {
						rec.seq[i] = toLowerByte(rec.seq[i])
					}
				}
			}
		}
		// 6. -U uppercase, else -x lowercase->mask char.
		if opts.Uppercase {
			for i := 0; i < len(rec.seq); i++ {
				rec.seq[i] = toUpperByte(rec.seq[i])
			}
		} else if opts.LowerToMask && opts.MaskChar > 0 {
			for i := 0; i < len(rec.seq); i++ {
				if isLowerByte(rec.seq[i]) {
					rec.seq[i] = opts.MaskChar
				}
			}
		}
		// 7. -A force FASTA, else -F fake quality.
		if opts.ForceFasta {
			rec.qual = nil
		} else if opts.FakeQual >= 33 && opts.FakeQual <= 127 {
			q := make([]byte, len(rec.seq))
			for i := range q {
				q[i] = byte(opts.FakeQual)
			}
			rec.qual = q
		}
		// 8. -C drop comments.
		if opts.DropComment {
			rec.comment = nil
		}
		// 9. -M region masking.
		if opts.maskRegions != nil {
			maskRegion(&rec, opts.maskRegions, opts.MaskComp, opts.MaskChar)
		}
		// 10. -r / -R reverse complement.
		if opts.RevComp || opts.BothStrands {
			if opts.BothStrands {
				if err := printSeqSuffix(w, &rec, lineLen, "+"); err != nil {
					return err
				}
			}
			reverseComplementInPlace(rec.seq)
			if len(rec.qual) > 0 {
				reverseInPlace(rec.qual)
			}
		}
		// 11. -V quality shift.
		if opts.ShiftQual && len(rec.qual) > 0 && opts.QualShift != 33 {
			delta := byte(opts.QualShift - 33)
			for i := range rec.qual {
				rec.qual[i] -= delta
			}
		}
		// 12. -N drop ambiguous (last step before printing).
		if opts.DropAmbig {
			ambiguous := false
			for i := 0; i < len(rec.seq); i++ {
				if seqNT16To4Table[seqNT16Table[rec.seq[i]]] > 3 {
					ambiguous = true
					break
				}
			}
			if ambiguous {
				continue
			}
		}
		// 13. Print.
		if opts.BothStrands {
			if err := printSeqSuffix(w, &rec, lineLen, "-"); err != nil {
				return err
			}
		} else {
			if err := printSeq(w, &rec, lineLen); err != nil {
				return err
			}
		}
	}
	return w.Flush()
}

func reverseComplementInPlace(s []byte) {
	n := len(s)
	for i := 0; i < n>>1; i++ {
		c0 := compTab[s[i]]
		c1 := compTab[s[n-1-i]]
		s[i] = c1
		s[n-1-i] = c0
	}
	if n&1 == 1 {
		s[n>>1] = compTab[s[n>>1]]
	}
}

func reverseInPlace(s []byte) {
	for i, j := 0, len(s)-1; i < j; i, j = i+1, j-1 {
		s[i], s[j] = s[j], s[i]
	}
}

// printStr writes the sequence/quality body with optional line wrapping,
// matching stk_printstr (seqtk.c:237).
func printStr(w *bufio.Writer, s []byte, lineLen int) error {
	if lineLen != 0 {
		rest := len(s)
		for i := 0; i < len(s); i += lineLen {
			if err := w.WriteByte('\n'); err != nil {
				return err
			}
			if rest > lineLen {
				if _, err := w.Write(s[i : i+lineLen]); err != nil {
					return err
				}
			} else {
				if _, err := w.Write(s[i:]); err != nil {
					return err
				}
			}
			rest -= lineLen
		}
		return w.WriteByte('\n')
	}
	if err := w.WriteByte('\n'); err != nil {
		return err
	}
	if _, err := w.Write(s); err != nil {
		return err
	}
	return w.WriteByte('\n')
}

func writeHeader(w *bufio.Writer, rec *kseqRecord, nameSuffix string) error {
	if len(rec.qual) > 0 {
		if err := w.WriteByte('@'); err != nil {
			return err
		}
	} else if err := w.WriteByte('>'); err != nil {
		return err
	}
	if _, err := w.Write(rec.name); err != nil {
		return err
	}
	if nameSuffix != "" {
		if _, err := w.WriteString(nameSuffix); err != nil {
			return err
		}
	}
	if len(rec.comment) > 0 {
		if err := w.WriteByte(' '); err != nil {
			return err
		}
		if _, err := w.Write(rec.comment); err != nil {
			return err
		}
	}
	return nil
}

// printSeq mirrors stk_printseq.
func printSeq(w *bufio.Writer, rec *kseqRecord, lineLen int) error {
	if err := writeHeader(w, rec, ""); err != nil {
		return err
	}
	if err := printStr(w, rec.seq, lineLen); err != nil {
		return err
	}
	if len(rec.qual) > 0 {
		if err := w.WriteByte('+'); err != nil {
			return err
		}
		if err := printStr(w, rec.qual, lineLen); err != nil {
			return err
		}
	}
	return nil
}

// printSeqSuffix mirrors stk_printseq_suffix (used by -R).
func printSeqSuffix(w *bufio.Writer, rec *kseqRecord, lineLen int, suffix string) error {
	if err := writeHeader(w, rec, suffix); err != nil {
		return err
	}
	if err := printStr(w, rec.seq, lineLen); err != nil {
		return err
	}
	if len(rec.qual) > 0 {
		if err := w.WriteByte('+'); err != nil {
			return err
		}
		if err := printStr(w, rec.qual, lineLen); err != nil {
			return err
		}
	}
	return nil
}
