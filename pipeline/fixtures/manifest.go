package fixtures

import (
	"crypto/sha256"
	"encoding/hex"
	"encoding/json"
	"fmt"
	"io"
	"os"
	"path/filepath"
)

// manifestVersion is bumped whenever the generation logic changes in a way that
// invalidates cached fixtures.
const manifestVersion = 1

// Manifest records what was generated for one scale tier so the runner can
// resolve fixture paths and decide whether a cached set is still valid.
type Manifest struct {
	Version int               `json:"version"`
	Scale   string            `json:"scale"`
	Seed    int64             `json:"seed"`
	Params  Params            `json:"params"`
	Files   map[string]string `json:"files"`   // logical name -> absolute path
	Sizes   map[string]int64  `json:"sizes"`   // logical name -> bytes
	Digests map[string]string `json:"digests"` // logical name -> sha256 (raw text inputs only)
}

// Path returns the absolute path for a logical fixture name (e.g. "bam",
// "vcf", "fasta"), or "" if absent.
func (m *Manifest) Path(name string) string { return m.Files[name] }

// manifestPath returns the manifest file location for a fixture dir.
func manifestPath(dir string) string { return filepath.Join(dir, "manifest.json") }

// loadManifest reads a manifest, returning (nil, nil) when absent.
func loadManifest(dir string) (*Manifest, error) {
	b, err := os.ReadFile(manifestPath(dir))
	if os.IsNotExist(err) {
		return nil, nil
	}
	if err != nil {
		return nil, err
	}
	var m Manifest
	if err := json.Unmarshal(b, &m); err != nil {
		return nil, err
	}
	return &m, nil
}

// valid reports whether a cached manifest matches the requested generation and
// all referenced files still exist with the recorded sizes.
func (m *Manifest) valid(scale Scale, seed int64) bool {
	if m == nil || m.Version != manifestVersion || m.Scale != string(scale) || m.Seed != seed {
		return false
	}
	for name, p := range m.Files {
		st, err := os.Stat(p)
		if err != nil {
			return false
		}
		if want, ok := m.Sizes[name]; ok && st.Size() != want {
			return false
		}
	}
	return true
}

// save writes the manifest as indented JSON.
func (m *Manifest) save(dir string) error {
	b, err := json.MarshalIndent(m, "", "  ")
	if err != nil {
		return err
	}
	return os.WriteFile(manifestPath(dir), b, 0o644)
}

// fileSize returns the size of p in bytes.
func fileSize(p string) (int64, error) {
	st, err := os.Stat(p)
	if err != nil {
		return 0, err
	}
	return st.Size(), nil
}

// sha256File returns the hex sha256 of a file's contents.
func sha256File(p string) (string, error) {
	f, err := os.Open(p)
	if err != nil {
		return "", err
	}
	defer f.Close()
	h := sha256.New()
	if _, err := io.Copy(h, f); err != nil {
		return "", err
	}
	return hex.EncodeToString(h.Sum(nil)), nil
}

// recordFile stamps a generated file's size (and digest for text inputs) into
// the manifest.
func (m *Manifest) recordFile(name, path string, withDigest bool) error {
	sz, err := fileSize(path)
	if err != nil {
		return fmt.Errorf("stat %s: %w", path, err)
	}
	m.Files[name] = path
	m.Sizes[name] = sz
	if withDigest {
		d, err := sha256File(path)
		if err != nil {
			return err
		}
		m.Digests[name] = d
	}
	return nil
}
