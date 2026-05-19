package sam

import (
	"bufio"
	"fmt"
	"io"
	"strconv"
	"strings"
)

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
	return &SAMReader{br: br, hdr: hdr}, nil
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
	pos, err := strconv.ParseInt(fields[3], 10, 32)
	if err != nil {
		return nil, fmt.Errorf("sam: bad POS %q: %w", fields[3], err)
	}
	rec.Pos = int32(pos)
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
	pnext, err := strconv.ParseInt(fields[7], 10, 32)
	if err != nil {
		return nil, fmt.Errorf("sam: bad PNEXT %q: %w", fields[7], err)
	}
	rec.PNext = int32(pnext)
	tlen, err := strconv.ParseInt(fields[8], 10, 32)
	if err != nil {
		return nil, fmt.Errorf("sam: bad TLEN %q: %w", fields[8], err)
	}
	rec.TLen = int32(tlen)
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
