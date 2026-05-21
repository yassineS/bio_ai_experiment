package alnio

import (
	"fmt"
	"io"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/cram"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// CRAMWriteOptions configures the CRAM writer adapter. Its zero value is
// the default writer — CRAM v3.0 with no lossy quality binning — exactly
// what NewCRAMWriter produces.
//
// The adapter always writes CRAM v3.0; the output version is not
// exposed here because the samtools CLI has no v3.1-output flag. A
// caller needing v3.1 output uses cram.NewRecordWriterOpts directly.
type CRAMWriteOptions struct {
	// QualityBinning selects the lossy quality-binning scheme applied to
	// each record's QUAL before encoding. The zero value, cram.BinningNone,
	// disables binning and keeps the writer losslessly exact.
	QualityBinning cram.QualityBinning
}

// ParseQualityBinning maps a CRAM quality-binning option string to a
// cram.QualityBinning scheme. It is the parser behind samtools view's
// `--output-fmt-option qbin=...`. Accepted values are "0"/"none"/""
// (no binning), "8"/"illumina-8", "4"/"illumina-4" and "2"/"illumina-2".
// An unrecognised value is an error.
func ParseQualityBinning(s string) (cram.QualityBinning, error) {
	switch s {
	case "", "0", "none", "off":
		return cram.BinningNone, nil
	case "8", "illumina-8", "illumina8":
		return cram.BinningIllumina8, nil
	case "4", "illumina-4", "illumina4":
		return cram.BinningIllumina4, nil
	case "2", "illumina-2", "illumina2":
		return cram.BinningIllumina2, nil
	default:
		return cram.BinningNone, fmt.Errorf("alnio: unknown CRAM quality-binning scheme %q (want none, 8, 4 or 2)", s)
	}
}

// cramWriter adapts a cram.RecordWriter to the sam.Writer interface so
// CRAM output can be wired interchangeably with the SAM and BAM writers.
//
// # Why this adapter exists
//
// sam.Writer has a three-step shape — WriteHeader, then Write per record,
// then Close — whereas cram.NewRecordWriter takes the SAM header at
// construction time (it must, because the header is the CRAM file's first
// container). cramWriter bridges the two: WriteHeader constructs and holds
// the underlying RecordWriter, Write delegates record by record, and Close
// finalises the CRAM file. It is the write-side mirror of alnio's
// CRAM-aware reader.
type cramWriter struct {
	w    io.Writer
	rw   *cram.RecordWriter
	opts CRAMWriteOptions
	err  error
}

// NewCRAMWriter returns a sam.Writer that encodes records to w as a CRAM
// v3.0 file. The header is not written until WriteHeader is called, since
// cram.RecordWriter needs the SAM header to open the file; WriteHeader
// must therefore be called once, before any Write, with a non-nil header.
//
// The returned writer produces reference-free CRAM (a self-contained file
// that decodes without an external reference FASTA), matching
// cram.RecordWriter's behaviour.
func NewCRAMWriter(w io.Writer) sam.Writer {
	return &cramWriter{w: w}
}

// NewCRAMWriterOpts is NewCRAMWriter with explicit configuration — it
// returns a sam.Writer that encodes CRAM to w using opts. The zero-value
// opts behaves exactly like NewCRAMWriter. When opts.QualityBinning
// selects a real scheme the writer applies lossy quality-score binning and
// records a provenance line in the embedded SAM header.
func NewCRAMWriterOpts(w io.Writer, opts CRAMWriteOptions) sam.Writer {
	return &cramWriter{w: w, opts: opts}
}

// WriteHeader constructs the underlying cram.RecordWriter from h and emits
// the CRAM file definition and SAM-header container. CRAM cannot be
// written without a header, so a nil h is an error. WriteHeader must be
// called exactly once and before the first Write.
func (cw *cramWriter) WriteHeader(h *sam.Header) error {
	if cw.err != nil {
		return cw.err
	}
	if cw.rw != nil {
		return fmt.Errorf("alnio: CRAM WriteHeader called twice")
	}
	if h == nil {
		cw.err = fmt.Errorf("alnio: CRAM output requires a SAM header")
		return cw.err
	}
	rw, err := cram.NewRecordWriterOpts(cw.w, h, cram.WriterOptions{
		Binning: cw.opts.QualityBinning,
	})
	if err != nil {
		cw.err = err
		return err
	}
	cw.rw = rw
	return nil
}

// Write encodes one record into the CRAM stream. It fails if WriteHeader
// has not been called yet.
func (cw *cramWriter) Write(rec *sam.Record) error {
	if cw.err != nil {
		return cw.err
	}
	if cw.rw == nil {
		return fmt.Errorf("alnio: CRAM Write before WriteHeader")
	}
	return cw.rw.Write(rec)
}

// Close flushes the final container and appends the CRAM EOF marker. It is
// safe to call Close when WriteHeader was never reached — for example for
// a header-only run that produced an error — in which case it is a no-op.
func (cw *cramWriter) Close() error {
	if cw.rw == nil {
		return cw.err
	}
	return cw.rw.Close()
}
