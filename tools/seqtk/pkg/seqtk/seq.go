// seq.go: an upstream-faithful port of `seqtk seq`, the common FASTA/Q
// transformation subcommand (stk_seq in reference_code/seqtk/seqtk.c).
//
// The earlier high-level helpers in seqtk.go (ReverseComplement, Filter,
// ConvertFastqToFasta) cover a few isolated transformations but do not match
// upstream's flag set or its byte-exact output layout. SeqRun reproduces
// stk_seq's behaviour for the option tail that matters for parity:
//
//	-A  force FASTA output (discard quality)
//	-C  drop the comment at header lines
//	-r  reverse complement
//	-l INT  number of residues per line (0 => no wrapping)
//	-q INT  mask bases with quality below INT
//	-X INT  mask bases with quality above INT
//	-n CHAR mask masked bases to CHAR (0 => lowercase)
//	-L INT  drop sequences shorter than INT
//	-U  convert all bases to uppercase
//	-V  shift quality by (-Q) - 33
//	-M FILE mask regions named in a BED or name-list file
//	-c  mask the complement of the -M regions
//	-N  drop sequences containing ambiguous bases
//
// The implementation mirrors upstream's ordering of operations precisely so
// the output is byte-for-byte identical.

package seqtk

import (
	"bufio"
	"io"
	"strconv"
	"strings"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fastq"
)

// compTab is upstream seqtk's complement table (comp_tab): it maps an ASCII
// byte to its complement, preserving case and leaving non-nucleotide bytes
// unchanged. It is copied verbatim from reference_code/seqtk/seqtk.c so the
// reverse-complement output matches upstream byte-for-byte.
var compTab = func() [256]byte {
	var t [256]byte
	for i := range t {
		t[i] = byte(i)
	}
	pairs := map[byte]byte{
		'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'U': 'A',
		'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K',
		'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D', 'N': 'N',
	}
	for from, to := range pairs {
		t[from] = to
		t[from|0x20] = to | 0x20 // lowercase
	}
	return t
}()

// SeqOptions holds the option-tail flags for SeqRun, mirroring stk_seq.
type SeqOptions struct {
	ForceFASTA   bool   // -A/-a: discard quality, emit FASTA
	DropComment  bool   // -C: drop the header comment
	RevComp      bool   // -r: reverse complement
	LineLen      int    // -l: residues per line (0 => no wrapping)
	QualThres    int    // -q: mask bases with quality below this (0 => disabled)
	MaxQual      int    // -X: mask bases with quality above this (default 255)
	QualShift    int    // -Q: ASCII offset for quality (default 33)
	MaskChar     byte   // -n: replacement for masked bases (0 => lowercase)
	MinLen       int    // -L: drop sequences shorter than this
	Uppercase    bool   // -U: convert all bases to uppercase
	ShiftQual    bool   // -V: shift quality by (-Q) - 33
	DropAmbig    bool   // -N: drop sequences containing ambiguous bases
	MaskComplent bool   // -c: mask the complement of -M regions
	MaskFile     string // -M: path to a BED or name-list mask file ("" => none)
}

// DefaultSeqOptions returns SeqOptions with upstream's defaults.
func DefaultSeqOptions() SeqOptions {
	return SeqOptions{MaxQual: 255, QualShift: 33}
}

// seqRecord is a minimal kseq-equivalent: a record with its name, comment,
// sequence and (optional) quality kept separate so output layout matches
// upstream exactly.
type seqRecord struct {
	name    string
	comment string
	seq     []byte
	qual    []byte // empty for FASTA
}

// splitNameComment splits a full FASTA/Q header line into the upstream "name"
// (up to the first run of whitespace) and "comment" (everything after that
// whitespace). It matches kseq.h's name/comment split.
func splitNameComment(header string) (name, comment string) {
	i := strings.IndexAny(header, " \t")
	if i < 0 {
		return header, ""
	}
	name = header[:i]
	rest := header[i:]
	rest = strings.TrimLeft(rest, " \t")
	return name, rest
}

// SeqRun reads a FASTA/FASTQ stream from in (auto-detected via the first
// non-whitespace byte) and writes the transformed records to w, applying the
// options in opts in the same order as upstream's stk_seq. It is a faithful
// port of `seqtk seq`'s option tail.
func SeqRun(in io.Reader, w io.Writer, opts SeqOptions) error {
	if opts.MaxQual == 0 {
		opts.MaxQual = 255
	}
	if opts.QualShift == 0 {
		opts.QualShift = 33
	}

	var maskRegions *regHash
	if opts.MaskFile != "" {
		mr, err := readRegHashFile(opts.MaskFile)
		if err != nil {
			return err
		}
		maskRegions = mr
	}

	br, isFastq := peekIsFastq(in)
	bw := bufio.NewWriter(w)

	// qualThres is the absolute quality cutoff in ASCII space.
	qualThres := opts.QualThres + opts.QualShift

	emit := func(rec *seqRecord) error {
		if len(rec.seq) < opts.MinLen {
			return nil
		}

		// -q / -X: quality masking (FASTQ only).
		if len(rec.qual) > 0 && qualThres > opts.QualShift {
			for i := 0; i < len(rec.seq); i++ {
				if int(rec.qual[i]) < qualThres || int(rec.qual[i]) > opts.MaxQual {
					if opts.MaskChar != 0 {
						rec.seq[i] = opts.MaskChar
					} else {
						rec.seq[i] = toLowerByte(rec.seq[i])
					}
				}
			}
		}

		// -U: uppercase everything.
		if opts.Uppercase {
			for i := 0; i < len(rec.seq); i++ {
				rec.seq[i] = toUpperByte(rec.seq[i])
			}
		}

		// -A: drop quality (fastq -> fasta).
		if opts.ForceFASTA {
			rec.qual = nil
		}

		// -C: drop comment.
		if opts.DropComment {
			rec.comment = ""
		}

		// -M: region masking.
		if maskRegions != nil {
			maskRecord(rec, maskRegions, opts.MaskComplent, opts.MaskChar)
		}

		// -r: reverse complement.
		if opts.RevComp {
			reverseComplementInPlace(rec.seq)
			if len(rec.qual) > 0 {
				reverseBytes(rec.qual)
			}
		}

		// -V: shift quality by (-Q) - 33.
		if opts.ShiftQual && len(rec.qual) > 0 && opts.QualShift != 33 {
			delta := byte(opts.QualShift - 33)
			for i := 0; i < len(rec.qual); i++ {
				rec.qual[i] -= delta
			}
		}

		// -N: drop sequences containing ambiguous bases (last step).
		if opts.DropAmbig {
			for i := 0; i < len(rec.seq); i++ {
				if seqNT16To4Table[seqNT16Table[rec.seq[i]]] > 3 {
					return nil
				}
			}
		}

		return writeSeqRecord(bw, rec, opts.LineLen)
	}

	if isFastq {
		r := fastq.NewReader(br, fastq.Phred33)
		fr := &fastq.Record{}
		rec := &seqRecord{}
		for {
			err := r.ReadInto(fr)
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}
			rec.name, rec.comment = splitNameComment(fr.Description)
			rec.seq, rec.qual = fr.Sequence, fr.Quality
			if err := emit(rec); err != nil {
				return err
			}
		}
	} else {
		r := fasta.NewReader(br)
		for {
			fr, err := r.Read()
			if err == io.EOF {
				break
			}
			if err != nil {
				return err
			}
			name, comment := splitNameComment(fr.Description)
			rec := &seqRecord{name: name, comment: comment, seq: fr.Sequence}
			if err := emit(rec); err != nil {
				return err
			}
		}
	}
	return bw.Flush()
}

// writeSeqRecord writes a single record in upstream layout: '>' or '@' marker,
// name, optional " comment", wrapped sequence, and (for FASTQ) the '+' line and
// wrapped quality.
func writeSeqRecord(bw *bufio.Writer, rec *seqRecord, lineLen int) error {
	marker := byte('>')
	if len(rec.qual) > 0 {
		marker = '@'
	}
	if err := bw.WriteByte(marker); err != nil {
		return err
	}
	if _, err := bw.WriteString(rec.name); err != nil {
		return err
	}
	if rec.comment != "" {
		if err := bw.WriteByte(' '); err != nil {
			return err
		}
		if _, err := bw.WriteString(rec.comment); err != nil {
			return err
		}
	}
	if err := writeWrapped(bw, rec.seq, lineLen); err != nil {
		return err
	}
	if len(rec.qual) > 0 {
		if _, err := bw.WriteString("+"); err != nil {
			return err
		}
		if err := writeWrapped(bw, rec.qual, lineLen); err != nil {
			return err
		}
	}
	return nil
}

// writeWrapped writes a leading newline, then s wrapped at lineLen columns
// (lineLen <= 0 => no wrapping), then a trailing newline — matching upstream's
// stk_printstr.
func writeWrapped(bw *bufio.Writer, s []byte, lineLen int) error {
	if lineLen <= 0 {
		if err := bw.WriteByte('\n'); err != nil {
			return err
		}
		if _, err := bw.Write(s); err != nil {
			return err
		}
		return bw.WriteByte('\n')
	}
	for i := 0; i < len(s); i += lineLen {
		if err := bw.WriteByte('\n'); err != nil {
			return err
		}
		end := i + lineLen
		if end > len(s) {
			end = len(s)
		}
		if _, err := bw.Write(s[i:end]); err != nil {
			return err
		}
	}
	return bw.WriteByte('\n')
}

// reverseComplementInPlace reverse-complements seq in place using compTab,
// matching upstream's stk_seq inner loop (including the middle base for
// odd-length sequences).
func reverseComplementInPlace(seq []byte) {
	n := len(seq)
	for i := 0; i < n/2; i++ {
		c0 := compTab[seq[i]]
		c1 := compTab[seq[n-1-i]]
		seq[i] = c1
		seq[n-1-i] = c0
	}
	if n&1 == 1 {
		seq[n/2] = compTab[seq[n/2]]
	}
}

// reverseBytes reverses b in place.
func reverseBytes(b []byte) {
	for i, j := 0, len(b)-1; i < j; i, j = i+1, j-1 {
		b[i], b[j] = b[j], b[i]
	}
}

func toUpperByte(b byte) byte {
	if b >= 'a' && b <= 'z' {
		return b - 32
	}
	return b
}

// regHash maps a sequence name to a list of [beg,end) intervals (0-based,
// half-open), reproducing upstream's khash(reg) built by stk_reg_read.
type regHash struct {
	m map[string][][2]int64
}

const regEndMax = int64(^uint64(0) >> 1) // INT64_MAX, upstream's "to end" sentinel

// RegionSet is an opaque set of named [begin, end) intervals parsed from a BED
// or name-list file, used by `comp -r` and `seq -M`.
type RegionSet = regHash

// ReadRegionFile parses a BED or name-list file (path, "-" for stdin, .gz
// supported) into a RegionSet, mirroring upstream seqtk's stk_reg_read. It is
// exported so the CLI can build the set for `comp -r`.
func ReadRegionFile(path string) (*RegionSet, error) {
	return readRegHashFile(path)
}

// readRegHashFile parses a BED or name-list file the way upstream's
// stk_reg_read does.
func readRegHashFile(path string) (*regHash, error) {
	rc, err := OpenInput(path)
	if err != nil {
		return nil, err
	}
	defer rc.Close()
	return parseRegHash(rc)
}

// parseRegHash reads a region/name file from r into a regHash, mirroring
// upstream's stk_reg_read: each line's first whitespace-delimited token is the
// name; if a second integer token is present it is the (0-based) begin and a
// third is the end. A single integer column is treated as a 1-based position
// (beg-1, beg). A name with no coordinates covers the whole sequence
// ([0, INT64_MAX)).
func parseRegHash(r io.Reader) (*regHash, error) {
	h := &regHash{m: make(map[string][][2]int64)}
	scanner := bufio.NewScanner(r)
	buf := make([]byte, 0, 64*1024)
	scanner.Buffer(buf, 10*1024*1024)
	for scanner.Scan() {
		line := strings.TrimRight(scanner.Text(), "\r")
		if line == "" {
			continue
		}
		fields := strings.Fields(line)
		if len(fields) == 0 {
			continue
		}
		name := fields[0]
		beg, end := int64(-1), int64(-1)
		if len(fields) >= 2 {
			if v, err := strconv.ParseInt(fields[1], 10, 64); err == nil {
				beg = v
				if len(fields) >= 3 {
					if w, err := strconv.ParseInt(fields[2], 10, 64); err == nil {
						end = w
						if end < 0 {
							end = -1
						}
					}
				}
			}
		}
		if end < 0 && beg > 0 {
			end = beg
			beg = beg - 1
		}
		if beg < 0 {
			beg = 0
			end = regEndMax
		}
		h.m[name] = append(h.m[name], [2]int64{beg, end})
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}
	return h, nil
}

// maskRecord applies region masking to rec.seq, mirroring upstream's stk_mask:
// without -c, the listed intervals are masked; with -c, everything OUTSIDE the
// listed intervals is masked (and a record absent from the hash is fully
// masked). maskChr == 0 means soft-mask (lowercase); otherwise bytes are set to
// maskChr.
func maskRecord(rec *seqRecord, h *regHash, complement bool, maskChr byte) {
	regions, ok := h.m[rec.name]
	n := int64(len(rec.seq))
	if !ok {
		if complement {
			for j := int64(0); j < n; j++ {
				if maskChr != 0 {
					rec.seq[j] = maskChr
				} else {
					rec.seq[j] = toLowerByte(rec.seq[j])
				}
			}
		}
		return
	}
	if !complement {
		for _, reg := range regions {
			beg, end := reg[0], reg[1]
			if beg >= n {
				continue
			}
			if end > n {
				end = n
			}
			for j := beg; j < end; j++ {
				if maskChr != 0 {
					rec.seq[j] = maskChr
				} else {
					rec.seq[j] = toLowerByte(rec.seq[j])
				}
			}
		}
		return
	}
	// Complement masking: mask every base not covered by a region.
	covered := make([]bool, n)
	for _, reg := range regions {
		beg, end := reg[0], reg[1]
		if end >= n {
			end = n
		}
		for j := beg; j < end; j++ {
			if j >= 0 && j < n {
				covered[j] = true
			}
		}
	}
	for j := int64(0); j < n; j++ {
		if !covered[j] {
			if maskChr != 0 {
				rec.seq[j] = maskChr
			} else {
				rec.seq[j] = toLowerByte(rec.seq[j])
			}
		}
	}
}
