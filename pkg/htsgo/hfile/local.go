package hfile

import (
	"io"
	"os"
)

// localHandle wraps an *os.File and adapts it to the Handle interface.
// *os.File already provides Read, ReadAt and Close with the required
// semantics; Size is derived from os.File.Stat.
type localHandle struct {
	f *os.File
}

// openLocal opens the local file at path for reading.
func openLocal(path string) (Handle, error) {
	if path == "-" {
		// Reading from stdin has no ReaderAt/Size semantics; callers that
		// need streaming stdin should not route through hfile.Open.
		return &localHandle{f: os.Stdin}, nil
	}
	f, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	return &localHandle{f: f}, nil
}

// Read reads sequentially from the underlying file.
func (h *localHandle) Read(p []byte) (int, error) {
	return h.f.Read(p)
}

// ReadAt reads len(p) bytes starting at byte offset off. It obeys the
// os.File.ReadAt contract, returning io.EOF at end of file.
func (h *localHandle) ReadAt(p []byte, off int64) (int, error) {
	return h.f.ReadAt(p, off)
}

// Size returns the size of the underlying file in bytes.
func (h *localHandle) Size() (int64, error) {
	fi, err := h.f.Stat()
	if err != nil {
		return 0, err
	}
	return fi.Size(), nil
}

// Close closes the underlying file. Closing the stdin handle is a no-op.
func (h *localHandle) Close() error {
	if h.f == os.Stdin {
		return nil
	}
	return h.f.Close()
}

var _ Handle = (*localHandle)(nil)
var _ io.ReaderAt = (*localHandle)(nil)
