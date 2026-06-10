package cram

import (
	"crypto/md5"
	"os"
	"path/filepath"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// stubReference is a hand-built ReferenceSource: a single named contig
// served from an in-memory string. It exercises the custom-source path
// of SetReference / referenceResolver without an on-disk FASTA.
type stubReference struct {
	name string
	seq  string
	err  error
}

func (s *stubReference) Fetch(name string, start, end int64) ([]byte, error) {
	if s.err != nil {
		return nil, s.err
	}
	if name != s.name || start < 0 || end > int64(len(s.seq)) || end < start {
		return nil, errFormat("stub reference: range %s:%d-%d out of bounds", name, start, end)
	}
	return []byte(s.seq[start:end]), nil
}

// TestResolverNoSource checks that a resolver with no source returns a
// nil span (the C4b 'N'-fill fallback) rather than an error.
func TestResolverNoSource(t *testing.T) {
	var rr *referenceResolver
	if rr.hasSource() {
		t.Error("a nil resolver must report no source")
	}
	bases, err := (&referenceResolver{}).sliceReference(&SliceHeader{}, "chr", "")
	if err != nil || bases != nil {
		t.Errorf("a sourceless resolver should return nil,nil; got %v,%v", bases, err)
	}
}

// TestResolverCustomSource checks the custom (non-FASTA, non-cache)
// ReferenceSource path: the span is fetched by name and MD5-verified.
func TestResolverCustomSource(t *testing.T) {
	seq := "ACGTACGTACGT"
	stub := &stubReference{name: "chr1", seq: seq}
	// Slice covering 1-based 3..8 of the contig (span GTACGT).
	span := seq[2:8]
	var digest [16]byte = md5.Sum([]byte(span))
	sh := &SliceHeader{RefSeqID: 0, AlignmentStart: 3, AlignmentSpan: 6, ReferenceMD5: digest}

	r := &referenceResolver{custom: stub}
	bases, err := r.sliceReference(sh, "chr1", "")
	if err != nil {
		t.Fatalf("sliceReference: %v", err)
	}
	if string(bases) != span {
		t.Errorf("custom-source span = %q, want %q", bases, span)
	}

	// A wrong MD5 must be a hard error.
	sh.ReferenceMD5[0] ^= 0xff
	if _, err := r.sliceReference(sh, "chr1", ""); err == nil {
		t.Error("a custom-source MD5 mismatch must error")
	}
}

// TestResolverInvalidSpan checks the span-validation guards.
func TestResolverInvalidSpan(t *testing.T) {
	r := &referenceResolver{custom: &stubReference{name: "c", seq: "ACGT"}}
	if _, err := r.sliceReference(&SliceHeader{AlignmentStart: 0, AlignmentSpan: 4}, "c", ""); err == nil {
		t.Error("an alignment start < 1 must error")
	}
	if _, err := r.sliceReference(&SliceHeader{AlignmentStart: 1, AlignmentSpan: -1}, "c", ""); err == nil {
		t.Error("a negative span must error")
	}
	// A zero-span slice yields an empty span with no error.
	bases, err := r.sliceReference(&SliceHeader{AlignmentStart: 1, AlignmentSpan: 0}, "c", "")
	if err != nil || len(bases) != 0 {
		t.Errorf("a zero-span slice should yield an empty span; got %v,%v", bases, err)
	}
}

// TestResolverCacheKey checks cacheKey's precedence: a well-formed @SQ
// M5 wins, a slice MD5 is the fallback, and neither is an error.
func TestResolverCacheKey(t *testing.T) {
	contigHex := "0123456789abcdef0123456789abcdef"
	var sliceMD5 [16]byte
	sliceMD5[0] = 0x42

	// The @SQ M5 takes precedence.
	k, err := cacheKey(contigHex, sliceMD5, "c", 1, 10)
	if err != nil {
		t.Fatalf("cacheKey: %v", err)
	}
	if formatMD5(k) != contigHex {
		t.Errorf("cacheKey chose %s, want the @SQ M5 %s", formatMD5(k), contigHex)
	}

	// A malformed @SQ M5 falls back to the slice MD5.
	k, err = cacheKey("not-hex", sliceMD5, "c", 1, 10)
	if err != nil || k != sliceMD5 {
		t.Errorf("cacheKey should fall back to the slice MD5; got %v,%v", k, err)
	}

	// Neither available is an error.
	if _, err := cacheKey("", [16]byte{}, "c", 1, 10); err == nil {
		t.Error("cacheKey with no @SQ M5 and no slice MD5 must error")
	}
}

// TestResolverCacheTooShort checks that a REF_CACHE file shorter than
// the slice span is a hard error.
func TestResolverCacheTooShort(t *testing.T) {
	dir := t.TempDir()
	short := []byte("ACGT")
	digest := md5.Sum(short)
	c := OpenRefCache(dir)
	path := c.refCachePath(formatMD5(digest))
	if err := os.MkdirAll(filepath.Dir(path), 0o755); err != nil {
		t.Fatalf("mkdir: %v", err)
	}
	if err := os.WriteFile(path, short, 0o644); err != nil {
		t.Fatalf("write: %v", err)
	}
	r := &referenceResolver{cache: c}
	// Ask for a span of 10 bases when the cache file is only 4.
	sh := &SliceHeader{AlignmentStart: 1, AlignmentSpan: 10, ReferenceMD5: digest}
	if _, err := r.sliceReference(sh, "c", formatMD5(digest)); err == nil {
		t.Error("a too-short REF_CACHE file must error")
	}
}

// TestResolverFASTAFallsBackToCache checks that when an explicit FASTA
// cannot serve a contig, a configured REF_CACHE is consulted next.
func TestResolverFASTAFallsBackToCache(t *testing.T) {
	// FASTA contains contig chrA only; the slice wants chrB.
	faPath := writeFASTA(t, "chrA", "ACGTACGT")
	fa, err := OpenFASTAReference(faPath)
	if err != nil {
		t.Fatalf("OpenFASTAReference: %v", err)
	}
	defer fa.Close()

	dir := t.TempDir()
	whole := []byte("TTTTGGGGCCCCAAAA")
	digest := md5.Sum(whole)
	c := OpenRefCache(dir)
	cachePath := c.refCachePath(formatMD5(digest))
	if err := os.MkdirAll(filepath.Dir(cachePath), 0o755); err != nil {
		t.Fatalf("mkdir: %v", err)
	}
	if err := os.WriteFile(cachePath, whole, 0o644); err != nil {
		t.Fatalf("write: %v", err)
	}

	r := &referenceResolver{fasta: fa, cache: c}
	span := whole[1:5] // 1-based 2..5
	spanMD5 := md5.Sum(span)
	sh := &SliceHeader{RefSeqID: 0, AlignmentStart: 2, AlignmentSpan: 4, ReferenceMD5: spanMD5}
	bases, err := r.sliceReference(sh, "chrB", formatMD5(digest))
	if err != nil {
		t.Fatalf("sliceReference (FASTA-miss, cache-hit): %v", err)
	}
	if string(bases) != string(span) {
		t.Errorf("cache-fallback span = %q, want %q", bases, span)
	}
}

// TestUseRefCacheFromEnv exercises the REF_CACHE env-var attach path.
func TestUseRefCacheFromEnv(t *testing.T) {
	rr := &RecordReader{}
	t.Setenv("REF_CACHE", "")
	if rr.UseRefCacheFromEnv() {
		t.Error("UseRefCacheFromEnv must report false when REF_CACHE is unset")
	}
	t.Setenv("REF_CACHE", t.TempDir())
	if !rr.UseRefCacheFromEnv() {
		t.Error("UseRefCacheFromEnv must report true when REF_CACHE is set")
	}
	if rr.refResolver == nil || rr.refResolver.cache == nil {
		t.Error("UseRefCacheFromEnv must attach a cache to the resolver")
	}
}

// TestSetReferenceCustom checks SetReference with a custom source and
// that refNameByID rejects an out-of-range id.
func TestSetReferenceCustom(t *testing.T) {
	rr := &RecordReader{refNames: []string{"chr1"}}
	rr.SetReference(&stubReference{name: "chr1", seq: "ACGT"})
	if rr.refResolver == nil || rr.refResolver.custom == nil {
		t.Error("SetReference with a custom source must attach it")
	}
	if _, err := rr.refNameByID(5); err == nil {
		t.Error("refNameByID must reject an out-of-range id")
	}
	if name, err := rr.refNameByID(0); err != nil || name != "chr1" {
		t.Errorf("refNameByID(0) = %q,%v; want chr1,nil", name, err)
	}
}

// TestSetReferenceFASTAMissing checks SetReferenceFASTA surfaces an
// error for a missing file rather than attaching a broken source.
func TestSetReferenceFASTAMissing(t *testing.T) {
	rr := &RecordReader{}
	if err := rr.SetReferenceFASTA(filepath.Join(t.TempDir(), "nope.fa")); err == nil {
		t.Error("SetReferenceFASTA on a missing file must error")
	}
}

// TestContigMD5 checks the @SQ M5-tag lookup: present, absent and
// out-of-range cases.
func TestContigMD5(t *testing.T) {
	rr := &RecordReader{}
	if rr.contigMD5(0) != "" {
		t.Error("contigMD5 with no header must return empty")
	}
	rr.header = &sam.Header{Refs: []sam.Reference{
		{Name: "chr1", Extra: []sam.HeaderField{{Tag: "M5", Value: "abc123"}}},
		{Name: "chr2"},
	}}
	if got := rr.contigMD5(0); got != "abc123" {
		t.Errorf("contigMD5(0) = %q, want abc123", got)
	}
	if got := rr.contigMD5(1); got != "" {
		t.Errorf("contigMD5(1) (no M5) = %q, want empty", got)
	}
	if got := rr.contigMD5(99); got != "" {
		t.Errorf("contigMD5(99) (out of range) = %q, want empty", got)
	}
}

// TestFASTASpanMemoised checks fastaSpan reuses its one-contig memo for
// a run of slices on the same contig.
func TestFASTASpanMemoised(t *testing.T) {
	faPath := writeFASTA(t, "chr1", "ACGTACGTACGTACGTACGT")
	fa, err := OpenFASTAReference(faPath)
	if err != nil {
		t.Fatalf("OpenFASTAReference: %v", err)
	}
	defer fa.Close()
	r := &referenceResolver{fasta: fa}
	first, err := r.fastaSpan("chr1", 1, 4)
	if err != nil {
		t.Fatalf("fastaSpan first: %v", err)
	}
	if string(first) != "ACGT" {
		t.Errorf("first span = %q, want ACGT", first)
	}
	if r.lastContig != "chr1" {
		t.Error("fastaSpan must memoise the contig")
	}
	// A second span on the same contig must come from the memo.
	second, err := r.fastaSpan("chr1", 5, 4)
	if err != nil {
		t.Fatalf("fastaSpan second: %v", err)
	}
	if string(second) != "ACGT" {
		t.Errorf("second span = %q, want ACGT", second)
	}
	// A span past the contig end is a hard error.
	if _, err := r.fastaSpan("chr1", 1, 999); err == nil {
		t.Error("a span past the contig end must error")
	}
}

// TestResolveSliceReferenceUnmapped checks that an unmapped (-1) or
// multi-reference (-2) slice resolves to a nil span even with a source
// attached — those slices fall back to the per-record 'N' path.
func TestResolveSliceReferenceUnmapped(t *testing.T) {
	rr := &RecordReader{
		refNames:    []string{"chr1"},
		refResolver: &referenceResolver{custom: &stubReference{name: "chr1", seq: "ACGT"}},
	}
	for _, id := range []int32{-1, -2} {
		sl := &Slice{
			Header:   &SliceHeader{RefSeqID: id, AlignmentStart: 1, AlignmentSpan: 4, EmbeddedRefID: -1},
			external: map[int32]*Block{},
		}
		bases, _, err := rr.resolveSliceReference(sl)
		if err != nil || bases != nil {
			t.Errorf("RefSeqID %d should resolve to nil span; got %v,%v", id, bases, err)
		}
	}
}
