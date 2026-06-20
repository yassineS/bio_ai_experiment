package sam

import (
	"bufio"
	"errors"
	"fmt"
	"io"
	"strconv"
	"strings"
)

// ErrNoSAMHeader reports that a text stream could not be read as SAM because
// it has no recognisable header and no leading alignment record — i.e. it is
// empty, blank, or arbitrary non-SAM text.
//
// This mirrors htslib's behaviour: hts_detect_format only classifies a
// header-less text stream as SAM when its first line is a complete (11-column)
// SAM record whose column types match (QNAME, integer FLAG, …). Anything else
// is "empty" or "unknown text", for which sam_hdr_read returns NULL and
// `samtools view` exits non-zero with "fail to read the header". A SAM with a
// valid header but zero alignment records is NOT covered by this error.
var ErrNoSAMHeader = errors.New("sam: no SAM header or leading alignment record (empty or non-SAM input)")

// Reader is the interface implemented by both SAM and BAM readers so that
// callers can iterate records without caring about the underlying format.
type Reader interface {
	// Header returns the parsed header. The result is valid for the lifetime
	// of the Reader.
	Header() *Header
	// Read returns the next alignment record, or io.EOF when the stream is
	// exhausted. The returned record is freshly allocated and is owned by
	// the caller.
	Read() (*Record, error)
}

// SAMReader reads alignment records from a text SAM stream.
type SAMReader struct {
	br     *bufio.Reader
	hdr    *Header
	err    error
	lineNo int
}

// NewSAMReader returns a SAMReader that consumes text SAM from r.
// The header is parsed eagerly on construction.
func NewSAMReader(r io.Reader) (*SAMReader, error) {
	br, ok := r.(*bufio.Reader)
	if !ok {
		br = bufio.NewReader(r)
	}
	hdr, _, err := ParseHeader(br)
	if err != nil {
		return nil, err
	}
	// Match htslib: a header-less text stream is only valid SAM when its very
	// first line is a complete alignment record. ParseHeader leaves that line
	// unconsumed in br, so peek (without consuming) and validate it. An empty
	// stream, a blank first line, or arbitrary non-SAM text is rejected here
	// the same way upstream's hts_detect_format declines to call it SAM.
	if len(hdr.Lines) == 0 {
		if err := validateLeadingSAMRecord(br); err != nil {
			return nil, err
		}
	}
	return &SAMReader{br: br, hdr: hdr}, nil
}

// validateLeadingSAMRecord peeks the first line still buffered in br (the line
// ParseHeader stopped at) and reports an error if it is not a parseable SAM
// alignment record. It does not consume any input, so the subsequent Read()
// still observes the record. It is only called when no header lines were seen.
func validateLeadingSAMRecord(br *bufio.Reader) error {
	// Grow the peek window until we have a full first line or hit EOF. Most
	// records are well under 4 KiB; long-CIGAR lines may exceed it, so keep
	// expanding rather than truncating the line under inspection.
	// Peek the bytes already buffered (ParseHeader has filled br by reading the
	// stream up to this line). Asking for the buffer's full capacity returns
	// everything available; a short read yields io.EOF, a long first line
	// yields ErrBufferFull.
	buf, err := br.Peek(br.Size())
	if nl := indexNewline(buf); nl >= 0 {
		buf = buf[:nl]
	} else if err == bufio.ErrBufferFull {
		// The first line is longer than br's buffer (e.g. a very long CIGAR or
		// tag list). It clearly contains real record data, so defer to Read()
		// rather than risk a false rejection on a truncated line.
		return nil
	}
	// buf is now the first line (newline found) or, at EOF with no trailing
	// newline, the entire remaining input.
	line := strings.TrimRight(string(buf), "\r\n")
	if line == "" {
		return ErrNoSAMHeader
	}
	if _, perr := parseSAMRecord(line); perr != nil {
		return ErrNoSAMHeader
	}
	return nil
}

// indexNewline returns the index of the first '\n' in b, or -1 if absent.
func indexNewline(b []byte) int {
	for i, c := range b {
		if c == '\n' {
			return i
		}
	}
	return -1
}

// Header returns the parsed SAM header.
func (sr *SAMReader) Header() *Header { return sr.hdr }

// Read returns the next alignment record from the stream, or io.EOF when no
// more records are available. Blank lines are skipped.
func (sr *SAMReader) Read() (*Record, error) {
	if sr.err != nil {
		return nil, sr.err
	}
	for {
		line, err := sr.br.ReadString('\n')
		if err != nil && err != io.EOF {
			sr.err = err
			return nil, err
		}
		line = strings.TrimRight(line, "\r\n")
		if line == "" {
			if err == io.EOF {
				sr.err = io.EOF
				return nil, io.EOF
			}
			continue
		}
		sr.lineNo++
		rec, perr := parseSAMRecord(line)
		if perr != nil {
			sr.err = perr
			return nil, perr
		}
		if err == io.EOF {
			// Defer the EOF so the caller observes this record on this call
			// and gets EOF on the next.
			sr.err = io.EOF
		}
		return rec, nil
	}
}

// parseSAMRecord parses one tab-delimited text SAM record line.
func parseSAMRecord(line string) (*Record, error) {
	fields := strings.Split(line, "\t")
	if len(fields) < 11 {
		return nil, fmt.Errorf("%w: only %d fields", errBadRecord, len(fields))
	}
	rec := &Record{}
	rec.QName = fields[0]
	flag, err := strconv.ParseUint(fields[1], 10, 16)
	if err != nil {
		return nil, fmt.Errorf("sam: bad FLAG %q: %w", fields[1], err)
	}
	rec.Flag = uint16(flag)
	if fields[2] == "*" {
		rec.RName = ""
	} else {
		rec.RName = fields[2]
	}
	pos, err := strconv.ParseInt(fields[3], 10, 64)
	if err != nil {
		return nil, fmt.Errorf("sam: bad POS %q: %w", fields[3], err)
	}
	rec.Pos = pos
	mapq, err := strconv.ParseUint(fields[4], 10, 8)
	if err != nil {
		return nil, fmt.Errorf("sam: bad MAPQ %q: %w", fields[4], err)
	}
	rec.MapQ = uint8(mapq)
	if fields[5] == "*" {
		rec.Cigar = nil
	} else {
		cig, err := ParseCigar(fields[5])
		if err != nil {
			return nil, err
		}
		rec.Cigar = cig
	}
	if fields[6] == "*" {
		rec.RNext = ""
	} else if fields[6] == "=" {
		rec.RNext = "="
	} else {
		rec.RNext = fields[6]
	}
	pnext, err := strconv.ParseInt(fields[7], 10, 64)
	if err != nil {
		return nil, fmt.Errorf("sam: bad PNEXT %q: %w", fields[7], err)
	}
	rec.PNext = pnext
	tlen, err := strconv.ParseInt(fields[8], 10, 64)
	if err != nil {
		return nil, fmt.Errorf("sam: bad TLEN %q: %w", fields[8], err)
	}
	rec.TLen = tlen
	if fields[9] != "*" {
		rec.Seq = fields[9]
	}
	if fields[10] != "*" {
		rec.Qual = phredASCIIToBytes(fields[10])
	}
	for _, f := range fields[11:] {
		a, err := ParseAux(f)
		if err != nil {
			return nil, err
		}
		rec.Aux = append(rec.Aux, a)
	}
	return rec, nil
}

// phredASCIIToBytes turns the SAM quality string (ASCII-33) into raw Phred
// scores.
func phredASCIIToBytes(s string) []byte {
	out := make([]byte, len(s))
	for i := 0; i < len(s); i++ {
		out[i] = s[i] - 33
	}
	return out
}
