package alnio

import (
	"bufio"
	"io"

	bgzf "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/cram"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// NewReaderThreaded is NewReader with block-parallel BGZF input decode wired to
// a thread count. When threads >= 2 and the input is a BGZF-wrapped BAM (the
// common on-disk case), the BGZF blocks are inflated concurrently across up to
// threads worker goroutines via bgzf.NewMultiReader, and the decoded BAM body is
// handed to the sam reader. The decoded byte stream — and therefore every record
// and every downstream output — is byte-for-byte identical to the
// single-threaded path for any thread count; threading only changes decode
// throughput, never the data. This mirrors the parallel WRITER already wired to
// -@ on the output side.
//
// threads < 2 (or input that is not BGZF) behaves exactly like NewReader: a
// plain SAM, a raw BAM body, or a plain-gzip SAM is read through the existing
// single-threaded path. CRAM carries its own container framing (slices, not BGZF
// blocks) and is decoded single-threaded here; parallel CRAM slice decode is a
// separate, larger piece of work (see the package docs) and is intentionally out
// of scope for the BGZF read-threading path.
//
// The returned sam.Reader does NOT own r. When the parallel reader is engaged
// its worker goroutines are torn down when the surrounding stream is fully read;
// callers that may abandon the stream early should prefer OpenReaderThreaded,
// whose Close releases the parallel reader's goroutines deterministically.
//
// referenceFASTA is honoured for CRAM input only (SAM and BAM carry their
// sequence inline); pass "" when no reference is available. It is the
// thread-aware drop-in for NewReaderWithReference.
func NewReaderThreaded(r io.Reader, referenceFASTA string, threads int) (sam.Reader, error) {
	if threads < 2 {
		return NewReaderWithReference(r, referenceFASTA)
	}
	br := bufio.NewReader(r)
	head, _ := br.Peek(16)
	if looksLikeCRAM(head) {
		rr, err := cram.NewRecordReader(br)
		if err != nil {
			return nil, err
		}
		if referenceFASTA != "" {
			if err := rr.SetReferenceFASTA(referenceFASTA); err != nil {
				rr.Close()
				return nil, err
			}
		}
		rr.UseRefCacheFromEnv()
		rr.UseRefPathFromEnv()
		return rr, nil
	}
	if !looksLikeBGZF(head) {
		// Plain SAM, raw BAM, or plain-gzip SAM: nothing BGZF to parallelise.
		dec, err := decompressStream(br)
		if err != nil {
			return nil, err
		}
		return sam.NewReader(dec)
	}
	mr, err := bgzf.NewMultiReader(br, threads)
	if err != nil {
		return nil, err
	}
	// mr inflates the BGZF stream into the raw BAM body ("BAM\1" magic …);
	// sam.NewReader detects that magic and skips its own BGZF layer. Because
	// sam.NewReader does not own mr, wrap the reader so Close tears down the
	// parallel decode worker goroutines.
	sr, err := sam.NewReader(mr)
	if err != nil {
		mr.Close()
		return nil, err
	}
	return &mrCloseReader{Reader: sr, mr: mr}, nil
}

// mrCloseReader wraps a sam.Reader whose BGZF layer is a parallel
// bgzf.MultiReader that sam does not own. Its Close releases the MultiReader's
// decode worker goroutines. It does not own the underlying source stream — the
// caller of NewReaderThreaded keeps that responsibility.
type mrCloseReader struct {
	sam.Reader
	mr *bgzf.MultiReader
}

// Close tears down the parallel BGZF decode workers.
func (m *mrCloseReader) Close() error { return m.mr.Close() }

// OpenReaderThreaded opens the alignment file at path and returns a Reader for
// it, engaging block-parallel BGZF input decode when threads >= 2 and the file
// is a BGZF-wrapped BAM. It is the thread-aware analogue of OpenReader: SAM, BAM
// and CRAM are auto-detected, referenceFASTA is honoured for CRAM, and a path of
// "-" or "" reads standard input.
//
// The decoded record stream is identical for any thread count; -@ only affects
// throughput. CRAM is decoded single-threaded (see NewReaderThreaded). The
// returned Reader's Close releases the underlying handle and, for the parallel
// BGZF path, tears down the decode worker goroutines.
func OpenReaderThreaded(path, referenceFASTA string, threads int) (Reader, error) {
	if threads < 2 {
		return OpenReader(path, referenceFASTA)
	}
	if isStdin(path) {
		return newThreadedReaderFromStream(io.NopCloser(stdinReader()), referenceFASTA, threads)
	}
	f, err := openAlnSource(path)
	if err != nil {
		return nil, err
	}
	rc, err := newThreadedReaderFromStream(f, referenceFASTA, threads)
	if err != nil {
		f.Close()
		return nil, err
	}
	return rc, nil
}

// newThreadedReaderFromStream builds a Reader over an already-open stream,
// sniffing its format and engaging the parallel BGZF reader for BGZF-wrapped BAM.
// The stream's Close — and the parallel reader's, when present — is chained into
// the returned Reader's Close so every handle and goroutine is released exactly
// once.
func newThreadedReaderFromStream(rc io.ReadCloser, referenceFASTA string, threads int) (Reader, error) {
	br := bufio.NewReader(rc)
	head, _ := br.Peek(16)
	if looksLikeCRAM(head) {
		// CRAM: single-threaded decode through the reference-aware reader, reusing
		// the buffered bytes already peeked.
		rr, err := cram.NewRecordReader(br)
		if err != nil {
			return nil, err
		}
		if referenceFASTA != "" {
			if err := rr.SetReferenceFASTA(referenceFASTA); err != nil {
				rr.Close()
				return nil, err
			}
		}
		rr.UseRefCacheFromEnv()
		rr.UseRefPathFromEnv()
		return &cramReader{rr: rr, src: rc}, nil
	}
	if !looksLikeBGZF(head) {
		dec, err := decompressStream(br)
		if err != nil {
			return nil, err
		}
		sr, err := sam.NewReader(dec)
		if err != nil {
			return nil, err
		}
		return &samReader{Reader: sr, src: rc}, nil
	}
	mr, err := bgzf.NewMultiReader(br, threads)
	if err != nil {
		return nil, err
	}
	sr, err := sam.NewReader(mr)
	if err != nil {
		mr.Close()
		return nil, err
	}
	return &threadedSamReader{Reader: sr, mr: mr, src: rc}, nil
}

// threadedSamReader wraps a sam.Reader whose BGZF layer is a parallel
// bgzf.MultiReader. Its Close tears down the decode workers and then releases the
// underlying source handle.
type threadedSamReader struct {
	sam.Reader
	mr  *bgzf.MultiReader
	src io.Closer
}

// Close releases the parallel BGZF decode workers and the underlying file
// handle. The first non-nil error is returned.
func (s *threadedSamReader) Close() error {
	err := s.mr.Close()
	if cerr := s.src.Close(); cerr != nil && err == nil {
		err = cerr
	}
	return err
}
