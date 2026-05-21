package cram

import (
	"crypto/md5"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// writeFASTA writes a one-contig FASTA (and its .fai) to a temp file and
// returns the FASTA path.
func writeFASTA(t *testing.T, name, seq string) string {
	t.Helper()
	dir := t.TempDir()
	path := filepath.Join(dir, "ref.fa")
	// 60 bases per line, the samtools faidx convention.
	var b strings.Builder
	b.WriteString(">")
	b.WriteString(name)
	b.WriteString("\n")
	for i := 0; i < len(seq); i += 60 {
		end := i + 60
		if end > len(seq) {
			end = len(seq)
		}
		b.WriteString(seq[i:end])
		b.WriteString("\n")
	}
	if err := os.WriteFile(path, []byte(b.String()), 0o644); err != nil {
		t.Fatalf("write FASTA: %v", err)
	}
	return path
}

// TestNormalizeReferenceBases checks the CRAM reference normalisation:
// lower-case is upper-cased and whitespace is stripped.
func TestNormalizeReferenceBases(t *testing.T) {
	got := normalizeReferenceBases([]byte("ac g\tt\nAC\rGT"))
	if string(got) != "ACGTACGT" {
		t.Errorf("normalizeReferenceBases = %q, want ACGTACGT", got)
	}
}

// TestReferenceMD5 checks the digest matches the CRAM/SAM convention:
// md5 of the upper-cased, base-only sequence.
func TestReferenceMD5(t *testing.T) {
	seq := []byte("ACGTACGTAC")
	got := referenceMD5(seq)
	want := md5.Sum([]byte("ACGTACGTAC"))
	if got != want {
		t.Errorf("referenceMD5 = %x, want %x", got, want)
	}
}

// TestVerifyReferenceMD5 covers the verify path: a matching digest
// passes, a mismatch is a hard error, and an all-zero expected digest
// (no MD5 recorded) skips the check.
func TestVerifyReferenceMD5(t *testing.T) {
	seq := []byte("ACGTACGT")
	good := md5.Sum(seq)
	if err := verifyReferenceMD5(seq, good, "chr", 1, 8); err != nil {
		t.Errorf("a matching MD5 must verify: %v", err)
	}
	var bad [16]byte
	bad[0] = 0xff
	if err := verifyReferenceMD5(seq, bad, "chr", 1, 8); err == nil {
		t.Error("a mismatched MD5 must be a hard error")
	}
	if err := verifyReferenceMD5(seq, [16]byte{}, "chr", 1, 8); err != nil {
		t.Errorf("an all-zero expected MD5 must skip the check: %v", err)
	}
}

// TestFASTAReference exercises the FASTA-backed reference source.
func TestFASTAReference(t *testing.T) {
	path := writeFASTA(t, "chr1", "ACGTACGTACGTACGTACGT")
	ref, err := OpenFASTAReference(path)
	if err != nil {
		t.Fatalf("OpenFASTAReference: %v", err)
	}
	defer ref.Close()
	bases, err := ref.Fetch("chr1", 4, 8)
	if err != nil {
		t.Fatalf("Fetch: %v", err)
	}
	if string(bases) != "ACGT" {
		t.Errorf("Fetch[4,8) = %q, want ACGT", bases)
	}
	if _, err := ref.Fetch("chrX", 0, 1); err == nil {
		t.Error("Fetch of an unknown contig should error")
	}
}

// TestOpenFASTAReferenceMissing checks a missing FASTA errors cleanly.
func TestOpenFASTAReferenceMissing(t *testing.T) {
	if _, err := OpenFASTAReference(filepath.Join(t.TempDir(), "nope.fa")); err == nil {
		t.Error("OpenFASTAReference on a missing file should error")
	}
}

// TestRefCachePathLayout checks the htslib %2s/%2s/%s cache path layout.
func TestRefCachePathLayout(t *testing.T) {
	c := OpenRefCache("/refs")
	// A 32-hex-char MD5 splits into 2 + 2 + 28.
	got := c.refCachePath("0123456789abcdef0123456789abcdef")
	want := filepath.Join("/refs", "01", "23", "456789abcdef0123456789abcdef")
	if got != want {
		t.Errorf("refCachePath = %q, want %q", got, want)
	}
	// A short digest cannot be split; it is joined whole so the lookup
	// fails cleanly rather than panicking.
	if c.refCachePath("ab") != filepath.Join("/refs", "ab") {
		t.Errorf("short-digest path = %q", c.refCachePath("ab"))
	}
}

// TestRefCacheResolveByMD5 writes a reference into a REF_CACHE-layout
// directory and resolves it by digest.
func TestRefCacheResolveByMD5(t *testing.T) {
	dir := t.TempDir()
	seq := []byte("ACGTACGTACGTACGT")
	digest := md5.Sum(seq)
	c := OpenRefCache(dir)
	path := c.refCachePath(formatMD5(digest))
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		t.Fatalf("mkdir: %v", err)
	}
	if err := os.WriteFile(path, seq, 0o644); err != nil {
		t.Fatalf("write cache file: %v", err)
	}
	got, err := c.ResolveByMD5(digest)
	if err != nil {
		t.Fatalf("ResolveByMD5: %v", err)
	}
	if string(got) != string(seq) {
		t.Errorf("ResolveByMD5 = %q, want %q", got, seq)
	}
	// A digest with no cache file is a clear not-found error.
	var missing [16]byte
	missing[0] = 0x99
	if _, err := c.ResolveByMD5(missing); err == nil {
		t.Error("ResolveByMD5 of a missing digest should error")
	}
}

// TestRefCacheFromEnv checks the REF_CACHE environment-variable factory.
func TestRefCacheFromEnv(t *testing.T) {
	t.Setenv("REF_CACHE", "")
	if _, ok := RefCacheFromEnv(); ok {
		t.Error("an empty REF_CACHE must report ok=false")
	}
	t.Setenv("REF_CACHE", "/some/cache")
	c, ok := RefCacheFromEnv()
	if !ok || c == nil {
		t.Fatal("a set REF_CACHE must yield a cache")
	}
	if c.dir != "/some/cache" {
		t.Errorf("cache dir = %q, want /some/cache", c.dir)
	}
}

// TestRefCacheLowercaseNormalised checks a cache file with lower-case or
// whitespace content is normalised on read.
func TestRefCacheLowercaseNormalised(t *testing.T) {
	dir := t.TempDir()
	upper := []byte("ACGTACGT")
	digest := md5.Sum(upper)
	c := OpenRefCache(dir)
	path := c.refCachePath(formatMD5(digest))
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		t.Fatalf("mkdir: %v", err)
	}
	// Store the lower-case, whitespace-padded form.
	if err := os.WriteFile(path, []byte("acgt\nacgt\n"), 0o644); err != nil {
		t.Fatalf("write: %v", err)
	}
	got, err := c.ResolveByMD5(digest)
	if err != nil {
		t.Fatalf("ResolveByMD5: %v", err)
	}
	if string(got) != "ACGTACGT" {
		t.Errorf("ResolveByMD5 normalised = %q, want ACGTACGT", got)
	}
}
