package cram

import (
	"crypto/md5"
	"encoding/hex"
	"os"
	"path/filepath"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
)

// ReferenceSource supplies the reference bases a reference-backed CRAM
// needs to reconstruct its mapped reads. A CRAM mapped read stores its
// sequence as a copy of a reference span plus the read features; the
// decoder fetches that span from a ReferenceSource.
//
// Fetch returns the half-open, 0-based base range [start, end) of the
// named reference sequence as upper-cased ASCII. The name is the SAM
// @SQ contig name. Implementations must return an error — never a short
// or padded result — when the range cannot be satisfied.
type ReferenceSource interface {
	// Fetch returns reference bases [start, end) of contig name, 0-based
	// and half-open, as upper-cased ASCII.
	Fetch(name string, start, end int64) ([]byte, error)
}

// FASTAReference is a ReferenceSource backed by an indexed FASTA file
// (the --reference equivalent). It uses faidx random access so a large
// reference is not held wholly in memory: each Fetch seeks and reads
// only the span it needs.
type FASTAReference struct {
	ra *fasta.RandomAccess
}

// OpenFASTAReference opens the named FASTA file as a ReferenceSource,
// loading (or building) its sibling .fai index for random access. The
// caller must Close it to release the file handle.
func OpenFASTAReference(path string) (*FASTAReference, error) {
	ra, err := fasta.OpenRandomAccess(path)
	if err != nil {
		return nil, wrapf(err, "opening reference FASTA %q", path)
	}
	return &FASTAReference{ra: ra}, nil
}

// Fetch returns the 0-based half-open base range [start, end) of the
// named contig, upper-cased.
func (f *FASTAReference) Fetch(name string, start, end int64) ([]byte, error) {
	return f.ra.Fetch(name, start, end)
}

// Close releases the underlying FASTA file handle.
func (f *FASTAReference) Close() error { return f.ra.Close() }

// RefCacheReference is a ReferenceSource backed by the htslib local
// reference cache: the directory named by the REF_CACHE environment
// variable, in which each reference sequence is stored in a single file
// named by the hex of its MD5 digest, laid out as %2s/%2s/%s — the
// first two hex characters, then the next two, then the remaining 28.
//
// REF_CACHE resolution is keyed by MD5, so a contig is resolved through
// ResolveByMD5 rather than by name. The network REF_PATH URL-fetch
// mechanism is deliberately not implemented; a cache miss is a clear
// error naming the missing digest.
type RefCacheReference struct {
	dir string
}

// OpenRefCache returns a RefCacheReference rooted at dir (typically the
// value of the REF_CACHE environment variable). It does not verify the
// directory exists; a missing entry surfaces on lookup.
func OpenRefCache(dir string) *RefCacheReference {
	return &RefCacheReference{dir: dir}
}

// RefCacheFromEnv returns a RefCacheReference for the REF_CACHE
// environment variable, or ok=false when REF_CACHE is unset or empty.
func RefCacheFromEnv() (*RefCacheReference, bool) {
	dir := os.Getenv("REF_CACHE")
	if dir == "" {
		return nil, false
	}
	return OpenRefCache(dir), true
}

// refCachePath returns the on-disk path of the cache file for the
// reference sequence with the given MD5, applying the htslib %2s/%2s/%s
// directory layout: the first two hex digits name the first directory,
// the next two name the second, and the remaining 28 name the file.
func (c *RefCacheReference) refCachePath(md5hex string) string {
	if len(md5hex) < 4 {
		// A short digest cannot be split; join it whole so the resulting
		// lookup fails with a clear not-found rather than a panic.
		return filepath.Join(c.dir, md5hex)
	}
	return filepath.Join(c.dir, md5hex[:2], md5hex[2:4], md5hex[4:])
}

// ResolveByMD5 returns the full reference sequence whose upper-cased,
// base-only MD5 digest is the given 16-byte value. It looks the digest
// up as a file under the REF_CACHE directory and reads it whole. The
// returned bytes are upper-cased so they can be hashed and sliced
// directly. A cache miss is reported as an error naming the digest.
func (c *RefCacheReference) ResolveByMD5(digest [16]byte) ([]byte, error) {
	hexStr := hex.EncodeToString(digest[:])
	path := c.refCachePath(hexStr)
	data, err := os.ReadFile(path)
	if err != nil {
		if os.IsNotExist(err) {
			return nil, errFormat("reference with MD5 %s not found in REF_CACHE %q", hexStr, c.dir)
		}
		return nil, wrapf(err, "reading REF_CACHE entry for MD5 %s", hexStr)
	}
	// A cache file holds the raw sequence (no FASTA header); upper-case
	// it and strip any stray whitespace so the bytes are decode-ready.
	return normalizeReferenceBases(data), nil
}

// normalizeReferenceBases returns a copy of data with every ASCII
// lower-case letter upper-cased and every whitespace byte dropped. This
// is the CRAM/SAM reference-MD5 normalisation: the digest in a slice
// header is md5 of the upper-cased, base-only sequence.
func normalizeReferenceBases(data []byte) []byte {
	out := make([]byte, 0, len(data))
	for _, b := range data {
		switch {
		case b == '\n' || b == '\r' || b == ' ' || b == '\t':
			continue
		case b >= 'a' && b <= 'z':
			out = append(out, b-('a'-'A'))
		default:
			out = append(out, b)
		}
	}
	return out
}

// referenceMD5 computes the CRAM reference MD5 of a base span: the md5
// digest of the upper-cased, base-only bytes, exactly the value a slice
// header carries. The input is assumed already upper-cased and
// whitespace-free (the FASTA faidx Fetch and the REF_CACHE normaliser
// both guarantee that); it is hashed verbatim.
func referenceMD5(bases []byte) [16]byte {
	return md5.Sum(bases)
}

// verifyReferenceMD5 checks that the md5 of the upper-cased base span
// matches the expected digest a CRAM slice header recorded. A mismatch
// is always a hard error: decoding against the wrong reference would
// silently produce wrong sequence, so the decoder refuses to proceed.
// An all-zero expected digest means the writer recorded no MD5 and the
// check is skipped.
func verifyReferenceMD5(bases []byte, expected [16]byte, contig string, start, span int32) error {
	if expected == ([16]byte{}) {
		return nil // no MD5 was recorded; nothing to verify against.
	}
	got := referenceMD5(bases)
	if got != expected {
		return errFormat("reference MD5 mismatch for %s:%d-%d: slice header expects %s but the supplied reference hashes to %s",
			contig, start, start+span-1,
			hex.EncodeToString(expected[:]), hex.EncodeToString(got[:]))
	}
	return nil
}

// formatMD5 renders a 16-byte digest as the conventional 32-character
// lower-case hex string.
func formatMD5(d [16]byte) string { return hex.EncodeToString(d[:]) }
