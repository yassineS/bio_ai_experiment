package hfile

import (
	"errors"
	"io"
	"os"
	"strings"
	"sync"
)

// SeekHandle is a read-seekable, randomly-addressable, closable file handle.
// A local *os.File satisfies it directly; a remote object is served by a
// read-ahead-buffered adapter over an hfile.Handle.
type SeekHandle interface {
	io.ReadSeeker
	io.ReaderAt
	io.Closer
}

// defaultWindow is the read-ahead window size for remote seekable handles. A
// BGZF/BAM/CRAM decoder issues many small sequential reads per block; fetching
// a large window per ranged GET coalesces them into a handful of HTTP requests
// instead of one request per field. 4 MiB comfortably spans several maximum
// (64 KiB) BGZF blocks while keeping memory modest.
const defaultWindow = 4 << 20

// OpenSeekable opens name for buffered random-access reading and returns a
// SeekHandle. A remote URL (http(s)://, s3://, gs://) is opened through hfile
// and wrapped in a read-ahead buffer so the many small sequential reads a
// BGZF/BAM/CRAM decoder performs coalesce into a few large ranged GETs; a
// local path (or file:// URL) is opened as an *os.File, whose ReadAt is
// already cheap and needs no buffering. The caller must Close the result.
func OpenSeekable(name string) (SeekHandle, error) {
	if !IsRemote(name) {
		return os.Open(strings.TrimPrefix(name, "file://"))
	}
	h, err := Open(name)
	if err != nil {
		return nil, err
	}
	size, err := h.Size()
	if err != nil {
		h.Close()
		return nil, err
	}
	return &bufSeekReader{h: h, size: size, window: defaultWindow, bufOff: -1}, nil
}

// bufSeekReader adapts a remote hfile.Handle into a SeekHandle with a single
// sliding read-ahead window. Because a decoder marches forward through a
// region, nearly every small read after the first in a window is served from
// memory, so the underlying handle sees only a few large ranged GETs.
//
// It is safe for concurrent use: the buffer is guarded by a mutex, so ReadAt
// may be called from multiple goroutines (each call still returns correct
// bytes, refilling the shared window as needed).
type bufSeekReader struct {
	h      Handle
	size   int64
	window int

	mu     sync.Mutex
	pos    int64  // logical offset for sequential Read/Seek
	buf    []byte // cached window contents
	bufOff int64  // file offset of buf[0]; -1 when empty
}

// fill ensures the window covers off and returns the slice of the window
// starting at off. The caller must hold r.mu. A returned slice may be shorter
// than the window near end-of-file.
func (r *bufSeekReader) fill(off int64) ([]byte, error) {
	if r.bufOff >= 0 && off >= r.bufOff && off < r.bufOff+int64(len(r.buf)) {
		return r.buf[off-r.bufOff:], nil
	}
	if off >= r.size {
		return nil, io.EOF
	}
	n := int64(r.window)
	if off+n > r.size {
		n = r.size - off
	}
	if cap(r.buf) < int(n) {
		r.buf = make([]byte, n)
	} else {
		r.buf = r.buf[:n]
	}
	got, err := r.h.ReadAt(r.buf, off)
	r.buf = r.buf[:got]
	r.bufOff = off
	if err != nil && err != io.EOF {
		return nil, err
	}
	if got == 0 {
		return nil, io.EOF
	}
	return r.buf, nil
}

// ReadAt implements io.ReaderAt by serving from (and refilling) the window. It
// follows the os.File.ReadAt contract: a short read returns a non-nil error,
// and reads at or past end-of-file return io.EOF.
func (r *bufSeekReader) ReadAt(p []byte, off int64) (int, error) {
	if off < 0 {
		return 0, errors.New("hfile: ReadAt: negative offset")
	}
	r.mu.Lock()
	defer r.mu.Unlock()
	n := 0
	for n < len(p) {
		win, err := r.fill(off + int64(n))
		if err != nil {
			return n, err
		}
		c := copy(p[n:], win)
		n += c
	}
	return n, nil
}

// Read implements io.Reader using the logical offset advanced by Seek.
func (r *bufSeekReader) Read(p []byte) (int, error) {
	r.mu.Lock()
	defer r.mu.Unlock()
	if r.pos >= r.size {
		return 0, io.EOF
	}
	win, err := r.fill(r.pos)
	if err != nil {
		return 0, err
	}
	n := copy(p, win)
	r.pos += int64(n)
	return n, nil
}

// Seek implements io.Seeker over the logical read offset.
func (r *bufSeekReader) Seek(offset int64, whence int) (int64, error) {
	r.mu.Lock()
	defer r.mu.Unlock()
	var abs int64
	switch whence {
	case io.SeekStart:
		abs = offset
	case io.SeekCurrent:
		abs = r.pos + offset
	case io.SeekEnd:
		abs = r.size + offset
	default:
		return 0, errors.New("hfile: Seek: invalid whence")
	}
	if abs < 0 {
		return 0, errors.New("hfile: Seek: negative position")
	}
	r.pos = abs
	return abs, nil
}

// Close releases the underlying handle.
func (r *bufSeekReader) Close() error { return r.h.Close() }

var _ SeekHandle = (*bufSeekReader)(nil)
