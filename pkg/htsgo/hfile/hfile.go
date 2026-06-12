// Package hfile is a pure-Go port of htslib's hfile virtual-filesystem
// remote backends. It provides read access to local files and to remote
// objects addressed by HTTP(S), Amazon S3 (s3://) and Google Cloud Storage
// (gs://) URLs, mirroring the behaviour of htslib's hfile_libcurl.c,
// hfile_s3.c and hfile_gcs.c.
//
// The package depends only on the Go standard library. AWS Signature
// Version 4 is implemented by hand using crypto/hmac and crypto/sha256, so
// no third-party SDKs are required.
//
// All handles expose io.ReaderAt for ranged random access, which is the
// access pattern used by the BGZF/BAM/CRAM index path. ReadAt is safe for
// concurrent use on the remote backends: each call issues an independent
// ranged GET and shares no mutable offset.
package hfile

import (
	"fmt"
	"io"
	"os"
	"strings"
)

// Handle is the common interface implemented by every hfile backend. It
// combines sequential reading (io.Reader), ranged random access
// (io.ReaderAt), closing (io.Closer) and a Size accessor reporting the
// total size of the underlying resource in bytes.
//
// ReadAt follows the os.File.ReadAt contract: it always returns a non-nil
// error when it returns fewer bytes than requested, and it returns io.EOF
// once the end of the resource is reached. This is the contract the BGZF
// seek reader relies upon.
type Handle interface {
	io.Reader
	io.ReaderAt
	io.Closer

	// Size returns the total size in bytes of the underlying resource, or
	// an error if the size cannot be determined.
	Size() (int64, error)
}

// Open opens name for reading and returns a Handle. The backend is chosen
// by URL scheme:
//
//	http://  https://   -> HTTP(S) backend
//	s3://bucket/key      -> Amazon S3 backend
//	gs://bucket/object   -> Google Cloud Storage backend
//	file:// or bare path -> local file backed by *os.File
//
// Write access is not implemented; remote handles are read-only.
func Open(name string) (Handle, error) {
	switch SchemeOf(name) {
	case "http", "https":
		return openHTTP(name)
	case "s3":
		return openS3(name)
	case "gs":
		return openGCS(name)
	case "file":
		return openLocal(strings.TrimPrefix(name, "file://"))
	case "":
		return openLocal(name)
	default:
		return nil, fmt.Errorf("hfile: unsupported URL scheme in %q", name)
	}
}

// ReadFile reads the entire contents of name and returns them. It is the
// hfile analogue of os.ReadFile: a remote URL is downloaded in full through
// the appropriate backend, and a local path is read from disk. It is intended
// for small sibling files such as BAI/CSI/TBI/CRAI indexes that accompany a
// remote alignment object.
func ReadFile(name string) ([]byte, error) {
	if !IsRemote(name) {
		return os.ReadFile(strings.TrimPrefix(name, "file://"))
	}
	h, err := Open(name)
	if err != nil {
		return nil, err
	}
	defer h.Close()
	// Prefer a single un-ranged GET: it is one request and stays correct even
	// against a server that ignores Range (returns 200 with the whole body),
	// which the ranged-ReadAt path would mis-stitch. Falls back to ReadAll for
	// any handle that does not support a whole-resource fetch.
	if wr, ok := h.(wholeReader); ok {
		return wr.readWhole()
	}
	return io.ReadAll(h)
}

// IsRemote reports whether name refers to a remote resource handled by one
// of the network backends (http, https, s3 or gs). Local paths and file://
// URLs return false.
func IsRemote(name string) bool {
	switch SchemeOf(name) {
	case "http", "https", "s3", "gs":
		return true
	default:
		return false
	}
}

// SchemeOf returns the lower-cased URL scheme of name, or the empty string
// if name carries no scheme (i.e. it is a bare filesystem path). A scheme
// is the run of characters before "://"; a bare relative or absolute path
// is therefore reported as having no scheme.
func SchemeOf(name string) string {
	i := strings.Index(name, "://")
	if i <= 0 {
		return ""
	}
	scheme := name[:i]
	// A valid scheme contains only ASCII letters, digits, '+', '-' or '.'.
	for _, r := range scheme {
		switch {
		case r >= 'a' && r <= 'z', r >= 'A' && r <= 'Z',
			r >= '0' && r <= '9', r == '+', r == '-', r == '.':
		default:
			return ""
		}
	}
	return strings.ToLower(scheme)
}
