// Package alnio is the format-aware entry point for reading SAM-family
// alignment streams. It sits one layer above pkg/htsgo/sam and
// pkg/htsgo/cram and chooses the right decoder for a stream — text SAM,
// BGZF-wrapped BAM, or CRAM — by sniffing its leading bytes.
//
// # Why this package exists
//
// pkg/htsgo/cram produces sam.Record values, so it imports pkg/htsgo/sam.
// That means pkg/htsgo/sam must not import pkg/htsgo/cram, and so the CRAM
// detection cannot live in the sam package. alnio is the layer that is
// allowed to depend on both: it routes CRAM streams to the CRAM record
// reader and everything else to sam.NewReader.
//
// Both cram.RecordReader and the sam package's readers expose the same
// Header()/Read() pair, so alnio.NewReader returns a sam.Reader regardless
// of the underlying format. Callers (notably the samtools subcommands)
// iterate records without caring which format they came from.
package alnio

import (
	"io"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/cram"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// NewReader returns a sam.Reader for r, auto-detecting whether the input is
// text SAM, BAM, or CRAM by sniffing its leading bytes. It is a drop-in,
// CRAM-aware replacement for sam.NewReader: a SAM or BAM stream is handed to
// sam.NewReader unchanged, while a CRAM stream is decoded by
// pkg/htsgo/cram and adapted to the sam.Reader interface.
//
// A CRAM stream is decoded without an external reference, so a
// reference-backed CRAM yields 'N' for reference-derived bases (the
// documented fallback behaviour). Use OpenReader when a reference is
// available, or drive cram.RecordReader directly for finer control.
func NewReader(r io.Reader) (sam.Reader, error) {
	format, sniffed, err := iohelper.DetectFormat(r)
	if err != nil {
		return nil, err
	}
	if format == iohelper.FormatCRAM {
		return cram.NewRecordReader(sniffed)
	}
	// SAM/BAM streams may be plain-gzip-compressed at the file level;
	// decompressStream strips that (a BGZF-wrapped BAM body is left for
	// sam.NewReader's own BAM reader to decode).
	dec, err := decompressStream(sniffed)
	if err != nil {
		return nil, err
	}
	return sam.NewReader(dec)
}

// NewReaderWithReference is NewReader with an explicit CRAM decode
// reference. It behaves identically to NewReader for SAM and BAM input —
// referenceFASTA is ignored, since SAM and BAM carry their sequence inline.
//
// For CRAM input, a non-empty referenceFASTA names a FASTA file (with a
// sibling .fai) used to reconstruct reference-backed mapped reads; the
// REF_CACHE environment variable is honoured as an additional reference
// source. When referenceFASTA is empty a reference-backed CRAM still
// decodes, with reference-derived bases filled with 'N' — the documented
// fallback. The reference must be attached before the first record is
// read, which this function guarantees.
func NewReaderWithReference(r io.Reader, referenceFASTA string) (sam.Reader, error) {
	format, sniffed, err := iohelper.DetectFormat(r)
	if err != nil {
		return nil, err
	}
	if format == iohelper.FormatCRAM {
		rr, err := cram.NewRecordReader(sniffed)
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
		return rr, nil
	}
	dec, err := decompressStream(sniffed)
	if err != nil {
		return nil, err
	}
	return sam.NewReader(dec)
}

// OpenReader opens the alignment file at path and returns a sam.Reader for
// it, auto-detecting SAM, BAM, or CRAM. A path of "-" or "" reads from
// standard input.
//
// referenceFASTA, when non-empty, names a FASTA file used as the decode
// reference for a reference-backed CRAM input; it is ignored for SAM and
// BAM input, which carry their sequence inline. When referenceFASTA is
// empty a reference-backed CRAM still decodes, with reference-derived bases
// filled with 'N' (the documented fallback). The REF_CACHE environment
// variable is honoured as an additional CRAM reference source.
//
// The returned ReadCloser releases the underlying file (and, for CRAM, any
// reference FASTA handle) when closed. SAM and BAM streams are read through
// iohelper, so a gzip- or BGZF-compressed SAM/BAM file is handled
// transparently; CRAM carries its own framing and is read raw.
func OpenReader(path, referenceFASTA string) (Reader, error) {
	if isStdin(path) {
		return newReaderFromStream(io.NopCloser(stdinReader()), referenceFASTA)
	}
	f, err := openAlnSource(path)
	if err != nil {
		return nil, err
	}
	rc, err := newReaderFromStream(f, referenceFASTA)
	if err != nil {
		f.Close()
		return nil, err
	}
	return rc, nil
}

// Reader is a sam.Reader that owns one or more underlying handles and must
// be closed by the caller. It is what OpenReader returns.
type Reader interface {
	sam.Reader
	io.Closer
}

// newReaderFromStream builds a Reader over an already-open stream, sniffing
// its format. The stream's Close is chained into the returned Reader's
// Close so the file handle is released exactly once.
func newReaderFromStream(rc io.ReadCloser, referenceFASTA string) (Reader, error) {
	format, sniffed, err := iohelper.DetectFormat(rc)
	if err != nil {
		return nil, err
	}
	switch format {
	case iohelper.FormatCRAM:
		rr, err := cram.NewRecordReader(sniffed)
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
		return &cramReader{rr: rr, src: rc}, nil
	default:
		// SAM/BAM streams may themselves be gzip/BGZF-compressed at the
		// file level; route through iohelper so that is handled. For a
		// raw (uncompressed) SAM/BAM stream iohelper is a transparent
		// pass-through.
		dec, err := decompressStream(sniffed)
		if err != nil {
			return nil, err
		}
		sr, err := sam.NewReader(dec)
		if err != nil {
			return nil, err
		}
		return &samReader{Reader: sr, src: rc}, nil
	}
}

// cramReader adapts a cram.RecordReader, which already satisfies
// sam.Reader, into a Reader whose Close also releases the source handle.
type cramReader struct {
	rr  *cram.RecordReader
	src io.Closer
}

// Header returns the SAM header parsed from the CRAM file.
func (c *cramReader) Header() *sam.Header { return c.rr.Header() }

// Read returns the next reconstructed alignment record, or io.EOF at end of
// stream.
func (c *cramReader) Read() (*sam.Record, error) { return c.rr.Read() }

// Close releases the CRAM reader (including any reference FASTA handle) and
// the underlying file.
func (c *cramReader) Close() error {
	err := c.rr.Close()
	if cerr := c.src.Close(); cerr != nil && err == nil {
		err = cerr
	}
	return err
}

// samReader wraps a sam.Reader so its Close releases the source handle.
type samReader struct {
	sam.Reader
	src io.Closer
}

// Close releases the underlying file handle.
func (s *samReader) Close() error { return s.src.Close() }
