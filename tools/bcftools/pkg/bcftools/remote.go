package bcftools

import (
	"bytes"
	"io"
	"os"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/hfile"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/tabix"
)

// readTabixIndex loads the sibling .tbi index for path, transparently handling
// remote URLs by downloading the (small) index bytes through hfile.
func readTabixIndex(path string) (*tabix.Index, error) {
	tbiPath := path + ".tbi"
	if hfile.IsRemote(path) {
		raw, err := hfile.ReadFile(tbiPath)
		if err != nil {
			return nil, err
		}
		return tabix.ReadGz(bytes.NewReader(raw))
	}
	return tabix.ReadFile(tbiPath)
}

// seekableCloser is a read-seekable handle that must be closed when done. A
// local *os.File satisfies it directly; a remote hfile handle is adapted to it
// via an io.SectionReader over the handle's ranged ReadAt. It mirrors the
// helper of the same name in tools/samtools and tools/tabix (other packages).
type seekableCloser interface {
	io.ReadSeeker
	io.Closer
}

// remoteSeekable adapts an hfile.Handle to a seekable, closable reader by
// layering an io.SectionReader (which provides Seek over ReadAt) on top of it
// while delegating Close to the underlying handle.
type remoteSeekable struct {
	*io.SectionReader
	h hfile.Handle
}

func (r *remoteSeekable) Close() error { return r.h.Close() }

// openSeekable opens path for indexed random access. A remote URL (http(s)://,
// s3://, gs://) is opened through hfile and wrapped so it presents the same
// io.ReadSeeker the local *os.File path provides; any other path is opened from
// disk. The caller must Close the result.
func openSeekable(path string) (seekableCloser, error) {
	if hfile.IsRemote(path) {
		h, err := hfile.Open(path)
		if err != nil {
			return nil, err
		}
		size, err := h.Size()
		if err != nil {
			h.Close()
			return nil, err
		}
		return &remoteSeekable{SectionReader: io.NewSectionReader(h, 0, size), h: h}, nil
	}
	return os.Open(path)
}

// siblingExists reports whether the sibling index at path exists and is
// readable. For local paths it is a cheap os.Stat; for remote URLs it probes
// the object with a single ranged read through hfile (the same probing htslib
// performs), treating any error as "absent". The HTTP request that 404s when
// the index is missing is the accepted cost of remote index discovery.
func siblingExists(path string) bool {
	if hfile.IsRemote(path) {
		h, err := hfile.Open(path)
		if err != nil {
			return false
		}
		defer h.Close()
		// A successful Size() means the object exists and is reachable.
		if _, err := h.Size(); err != nil {
			return false
		}
		return true
	}
	_, err := os.Stat(path)
	return err == nil
}
