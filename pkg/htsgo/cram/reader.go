package cram

import (
	"bufio"
	"fmt"
	"io"
	"os"
)

// Container is one fully-parsed CRAM container: its header and the
// blocks that follow it. The first container of a file holds the SAM
// text header; every later container holds a compression-header block
// followed by one or more slices (a slice-header block plus its data
// blocks).
type Container struct {
	// Header is the container's parsed header.
	Header ContainerHeader
	// Blocks holds the container's blocks in file order. It is nil for
	// the EOF marker container, whose Header.IsEOF is set.
	Blocks []Block
	// Index is the container's zero-based position in the file.
	Index int
	// Major is the CRAM major version of the file the container came
	// from. It is threaded into the slice-header parse because the
	// record-counter field's width (ITF-8 for v2, LTF-8 for v3+) depends
	// on it. Containers produced by Reader.Next carry it; a zero value
	// (only possible for a hand-built Container) is treated as v3+.
	Major uint8
}

// Reader walks the structural tree of a CRAM stream: its file
// definition, then each container and the blocks within. It validates
// the per-structure CRC32 checksums that CRAM v3+ embeds. It does not
// decode the alignment data series; call Block.Decompress to obtain a
// block's uncompressed payload.
//
// A Reader is not safe for concurrent use; Next advances shared state.
type Reader struct {
	def     FileDefinition
	r       *bufio.Reader
	closer  io.Closer
	count   int
	done    bool
	sawEOF  bool
	lastErr error
}

// NewReader reads and validates the CRAM file definition from r and
// returns a Reader positioned at the first container. It returns an
// error if r does not begin with a recognised CRAM file definition.
func NewReader(r io.Reader) (*Reader, error) {
	br := bufio.NewReader(r)
	def, err := readFileDefinition(br)
	if err != nil {
		return nil, err
	}
	return &Reader{def: def, r: br}, nil
}

// Open opens the named CRAM file and returns a Reader over it. The
// caller must call Close to release the underlying file handle.
func Open(path string) (*Reader, error) {
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	rd, err := NewReader(f)
	if err != nil {
		f.Close()
		return nil, err
	}
	rd.closer = f
	return rd, nil
}

// Close releases the underlying file handle if the Reader was created
// by Open. It is a no-op for a Reader created by NewReader.
func (rd *Reader) Close() error {
	if rd.closer != nil {
		return rd.closer.Close()
	}
	return nil
}

// FileDefinition returns the parsed CRAM file definition.
func (rd *Reader) FileDefinition() FileDefinition { return rd.def }

// Next reads and returns the next container from the stream, including
// all of its blocks. It returns io.EOF once the file is exhausted: the
// CRAM EOF marker container is consumed and reported as io.EOF rather
// than as a data container. A file that runs out before the EOF marker
// is reported as a truncation error, not io.EOF, because a well-formed
// CRAM file always ends with the marker.
func (rd *Reader) Next() (*Container, error) {
	if rd.done {
		if rd.lastErr != nil {
			return nil, rd.lastErr
		}
		return nil, io.EOF
	}
	hdr, err := readContainerHeader(rd.r, rd.def)
	if err == io.EOF {
		// A well-formed CRAM v3 file ends with the EOF marker
		// container; reaching the end of the stream without having
		// seen it means the file is truncated or incomplete.
		rd.done = true
		if !rd.sawEOF {
			rd.lastErr = fmt.Errorf("cram: stream ended after %d container(s) without the CRAM EOF marker (file truncated or incomplete)",
				rd.count)
			return nil, rd.lastErr
		}
		return nil, io.EOF
	}
	if err != nil {
		rd.done = true
		rd.lastErr = err
		return nil, err
	}
	if hdr.IsEOF {
		rd.done = true
		rd.sawEOF = true
		return nil, io.EOF
	}
	c := &Container{Header: hdr, Index: rd.count, Major: rd.def.Major}
	rd.count++
	// The container header's Length field is the byte size of all the
	// blocks that follow. Walking by NumBlocks is the simplest correct
	// approach: each block is self-delimiting (its header carries its
	// own compressed size), and slices contribute their slice-header
	// block plus data blocks to the same flat count.
	limited := &io.LimitedReader{R: rd.r, N: int64(hdr.Length)}
	// Each block occupies at least one byte of the container body, so a
	// declared block count larger than the body is corrupt.
	if int64(hdr.NumBlocks) > int64(hdr.Length) {
		rd.done = true
		rd.lastErr = fmt.Errorf("cram: container %d declares %d blocks but only %d body bytes",
			c.Index, hdr.NumBlocks, hdr.Length)
		return nil, rd.lastErr
	}
	// c.Blocks grows by append rather than being pre-sized to NumBlocks:
	// NumBlocks is untrusted (htslib even documents that it can be
	// stale) and may dwarf the bytes actually present, so a pre-sized
	// make would risk an enormous allocation on malformed input.
	for i := int32(0); i < hdr.NumBlocks; i++ {
		blk, err := readBlock(limited, rd.def, limited.N)
		if err != nil {
			rd.done = true
			rd.lastErr = fmt.Errorf("cram: container %d, block %d: %w", c.Index, i, err)
			return nil, rd.lastErr
		}
		c.Blocks = append(c.Blocks, blk)
	}
	if limited.N != 0 {
		rd.done = true
		rd.lastErr = fmt.Errorf("cram: container %d declared length %d but %d bytes unconsumed after %d blocks",
			c.Index, hdr.Length, limited.N, hdr.NumBlocks)
		return nil, rd.lastErr
	}
	return c, nil
}

// Containers reads every remaining container from the stream and
// returns them as a slice. It is a convenience wrapper over repeated
// Next calls; it returns whatever containers it parsed alongside any
// error encountered before the EOF.
func (rd *Reader) Containers() ([]*Container, error) {
	var out []*Container
	for {
		c, err := rd.Next()
		if err == io.EOF {
			return out, nil
		}
		if err != nil {
			return out, err
		}
		out = append(out, c)
	}
}
