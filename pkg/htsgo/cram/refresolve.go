package cram

import (
	"encoding/hex"
	"errors"
)

// referenceResolver resolves the reference bases a slice's mapped
// records need, threading together the two reference sources C5
// supports: an explicit FASTA file and the htslib REF_CACHE directory.
// It also memoises the last whole sequence it loaded so that a run of
// slices on the same contig does not re-read the reference.
//
// A nil *referenceResolver means no reference was supplied: the decoder
// then falls back to the C4b behaviour of filling reference-derived
// bases with 'N' and setting NeedsReference.
type referenceResolver struct {
	// fasta is the explicit --reference FASTA, or nil.
	fasta *FASTAReference
	// cache is the REF_CACHE local cache, or nil.
	cache *RefCacheReference
	// custom is a caller-supplied ReferenceSource that is neither a
	// FASTAReference nor a RefCacheReference; it is name-addressed like
	// the FASTA path. It is nil unless SetReference was given one.
	custom ReferenceSource

	// lastContig and lastBases memoise the most recently fetched whole
	// reference sequence (FASTA path) so consecutive slices on the same
	// contig reuse it instead of re-reading.
	lastContig string
	lastBases  []byte
}

// hasSource reports whether the resolver can supply any reference at
// all. A resolver with no FASTA, REF_CACHE or custom source behaves as
// if no reference were supplied.
func (rr *referenceResolver) hasSource() bool {
	return rr != nil && (rr.fasta != nil || rr.cache != nil || rr.custom != nil)
}

// sliceReference resolves and MD5-verifies the reference span a slice
// covers. It returns the upper-cased bases for the half-open 0-based
// range [start-1, start-1+span) of the slice's contig, indexed so that
// span[0] is reference position AlignmentStart.
//
// Resolution order: an explicit FASTA is consulted by contig name; a
// REF_CACHE is consulted by the contig's whole-sequence MD5 — the @SQ
// header's M5 tag (contigMD5), the key htslib's reference cache uses —
// falling back to the slice header's own MD5 when no M5 tag is present.
// Whichever source supplies the bases, the slice header's MD5 (when
// non-zero) is verified against the span — a mismatch is always a hard
// error.
func (rr *referenceResolver) sliceReference(sh *SliceHeader, contig, contigMD5 string) ([]byte, error) {
	if !rr.hasSource() {
		return nil, nil
	}
	start := sh.AlignmentStart
	span := sh.AlignmentSpan
	if start < 1 || span < 0 {
		return nil, errFormat("slice declares an invalid reference span (start %d, span %d)", start, span)
	}
	if span == 0 {
		return []byte{}, nil
	}

	// priorErr accumulates the errors from sources tried before the
	// cache, so a final cache failure reports every path that was tried
	// rather than just the last one.
	var priorErr error

	// Prefer the explicit FASTA: it is name-addressed and can serve any
	// span without the whole-sequence MD5 the cache needs.
	if rr.fasta != nil {
		bases, err := rr.fastaSpan(contig, start, span)
		if err == nil {
			if verr := verifyReferenceMD5(bases, sh.ReferenceMD5, contig, start, span); verr != nil {
				return nil, verr
			}
			return bases, nil
		}
		// Fall through to the cache only when a cache exists; otherwise
		// the FASTA error is the most informative one to return.
		if rr.cache == nil {
			return nil, err
		}
		priorErr = err
	}

	// A custom name-addressed source: fetch the span directly and verify.
	if rr.custom != nil {
		bases, err := rr.custom.Fetch(contig, int64(start-1), int64(start-1)+int64(span))
		if err == nil {
			bases = normalizeReferenceBases(bases)
			if verr := verifyReferenceMD5(bases, sh.ReferenceMD5, contig, start, span); verr != nil {
				return nil, verr
			}
			return bases, nil
		}
		if rr.cache == nil {
			return nil, err
		}
		priorErr = errors.Join(priorErr, err)
	}

	// REF_CACHE: the whole reference sequence is stored in a file named
	// by its MD5. The key is the @SQ M5 tag (the whole-sequence digest)
	// when the header carries one; failing that, the slice header's MD5.
	if rr.cache != nil {
		key, err := cacheKey(contigMD5, sh.ReferenceMD5, contig, start, span)
		if err != nil {
			return nil, errors.Join(priorErr, err)
		}
		whole, err := rr.cache.ResolveByMD5(key)
		if err != nil {
			return nil, errors.Join(priorErr, err)
		}
		end := int(start-1) + int(span)
		if end > len(whole) || start < 1 {
			return nil, errFormat("REF_CACHE reference for MD5 %s is %d bases, too short for slice span %d-%d",
				formatMD5(key), len(whole), start, start+span-1)
		}
		spanBases := whole[start-1 : end]
		// The cache key proved the whole sequence; the slice-span MD5
		// still must verify that the slice header agrees on the span.
		if verr := verifyReferenceMD5(spanBases, sh.ReferenceMD5, contig, start, sh.AlignmentSpan); verr != nil {
			return nil, verr
		}
		return spanBases, nil
	}
	return nil, errFormat("no reference source could resolve %s:%d-%d", contig, start, start+span-1)
}

// cacheKey selects the MD5 a REF_CACHE lookup keys on for a slice: the
// @SQ M5 tag (the whole-sequence digest, contigMD5Hex) when present and
// well-formed, otherwise the slice header's own MD5. It errors only when
// neither is available — a slice that cannot be addressed at all.
func cacheKey(contigMD5Hex string, sliceMD5 [16]byte, contig string, start, span int32) ([16]byte, error) {
	if contigMD5Hex != "" {
		raw, err := hex.DecodeString(contigMD5Hex)
		if err == nil && len(raw) == 16 {
			var k [16]byte
			copy(k[:], raw)
			return k, nil
		}
	}
	if sliceMD5 != ([16]byte{}) {
		return sliceMD5, nil
	}
	return [16]byte{}, errFormat("cannot resolve %s:%d-%d from REF_CACHE: no @SQ M5 tag and no slice reference MD5",
		contig, start, start+span-1)
}

// fastaSpan fetches a reference span from the FASTA source. It memoises
// the whole contig: a CRAM file's slices are contig-ordered, so a run of
// slices on one contig reuses a single Fetch of the whole sequence
// rather than seeking per slice. The memo is bounded to one contig.
func (rr *referenceResolver) fastaSpan(contig string, start, span int32) ([]byte, error) {
	end := int64(start-1) + int64(span)
	if rr.lastContig == contig && rr.lastBases != nil {
		if end <= int64(len(rr.lastBases)) {
			return rr.lastBases[start-1 : end], nil
		}
	}
	// Fetch the whole contig once and memoise it. faidx random access
	// keeps this a single seek+read; holding one contig is acceptable
	// for v1 and avoids a per-slice seek.
	whole, err := rr.fasta.Fetch(contig, 0, contigLength(rr.fasta, contig))
	if err != nil {
		// Fall back to fetching just the span when the whole-contig
		// length is unavailable (the index may lack the contig).
		bases, ferr := rr.fasta.Fetch(contig, int64(start-1), end)
		if ferr != nil {
			return nil, ferr
		}
		return bases, nil
	}
	rr.lastContig = contig
	rr.lastBases = whole
	if end > int64(len(whole)) {
		return nil, errFormat("reference contig %q is %d bases, too short for slice span %d-%d",
			contig, len(whole), start, int64(start)+int64(span)-1)
	}
	return whole[start-1 : end], nil
}

// contigLength returns the length of a contig in the FASTA index, or 0
// when the contig is unknown (the caller then falls back to a span
// fetch).
func contigLength(f *FASTAReference, contig string) int64 {
	n := f.ra.Length(contig)
	if n < 0 {
		return 0
	}
	return n
}
