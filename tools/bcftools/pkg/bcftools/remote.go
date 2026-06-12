package bcftools

import (
	"bytes"
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

// openSeekable opens path for indexed random access. A remote URL (http(s)://,
// s3://, gs://) is opened through hfile with read-ahead buffering so the many
// small sequential reads the BGZF decoder performs coalesce into a few large
// ranged GETs; a local path is opened from disk. The caller must Close the
// result.
func openSeekable(path string) (hfile.SeekHandle, error) {
	return hfile.OpenSeekable(path)
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
