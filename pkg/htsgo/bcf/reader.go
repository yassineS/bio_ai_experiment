package bcf

import (
	"bufio"
	"encoding/binary"
	"io"
)

// Reader sequentially decodes BCF records from a stream. Callers obtain a
// Reader from NewReader (which auto-reads the header) or by hand using
// ReadHeader + NewReaderWithHeader.
//
// The reader does not own the underlying io.Reader; callers are responsible
// for closing the source.
type Reader struct {
	br     *bufio.Reader
	header *Header
}

// NewReader reads the magic + text header from r and returns a Reader
// positioned at the first record. r should already be BGZF-decompressed —
// pkg/htsgo/iohelper.OpenReader is the easy way to get that.
func NewReader(r io.Reader) (*Reader, error) {
	hdr, err := ReadHeader(r)
	if err != nil {
		return nil, err
	}
	return &Reader{br: hdr.tailReader(), header: hdr}, nil
}

// NewReaderWithHeader returns a Reader that uses an already-parsed Header.
// Use this when you need to inspect the header before constructing the
// Reader (e.g. to set up filters that depend on the dictionary contents).
func NewReaderWithHeader(hdr *Header) *Reader {
	return &Reader{br: hdr.tailReader(), header: hdr}
}

// Header returns the parsed header. It is safe to call concurrently with
// Read because the header is immutable after parse.
func (r *Reader) Header() *Header { return r.header }

// Read returns the next record. It returns io.EOF when the stream is
// exhausted.
func (r *Reader) Read() (*Record, error) {
	var lShared, lIndiv uint32
	if err := binary.Read(r.br, binary.LittleEndian, &lShared); err != nil {
		if err == io.EOF {
			return nil, io.EOF
		}
		return nil, wrapEOF(err)
	}
	if err := binary.Read(r.br, binary.LittleEndian, &lIndiv); err != nil {
		return nil, wrapEOF(err)
	}

	sharedBuf := make([]byte, lShared)
	if _, err := io.ReadFull(r.br, sharedBuf); err != nil {
		return nil, wrapEOF(err)
	}
	indivBuf := make([]byte, lIndiv)
	if _, err := io.ReadFull(r.br, indivBuf); err != nil {
		return nil, wrapEOF(err)
	}

	rec := &Record{}
	if err := decodeShared(sharedBuf, r.header, rec); err != nil {
		return nil, err
	}
	if err := decodeIndiv(indivBuf, r.header, rec); err != nil {
		return nil, err
	}
	return rec, nil
}

// ReadAll consumes the rest of the stream and returns the decoded records.
func (r *Reader) ReadAll() ([]*Record, error) {
	var out []*Record
	for {
		rec, err := r.Read()
		if err == io.EOF {
			return out, nil
		}
		if err != nil {
			return out, err
		}
		out = append(out, rec)
	}
}
