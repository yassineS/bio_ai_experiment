package sam

import (
	"bufio"
	"io"
	"strconv"
)

// Writer is the common interface for SAM and BAM writers.
type Writer interface {
	WriteHeader(h *Header) error
	Write(rec *Record) error
	Close() error
}

// SAMWriter emits records as text SAM.
type SAMWriter struct {
	bw         *bufio.Writer
	headerDone bool
}

// NewSAMWriter wraps w in a SAMWriter. Writers are responsible for calling
// WriteHeader before Write; passing a nil header is allowed for headerless
// output but is not generally compatible with downstream tools.
func NewSAMWriter(w io.Writer) *SAMWriter {
	return &SAMWriter{bw: bufio.NewWriter(w)}
}

// WriteHeader emits the SAM header to the underlying writer. It is safe to
// call with a nil h (no header emitted, header marked as done).
func (sw *SAMWriter) WriteHeader(h *Header) error {
	sw.headerDone = true
	if h == nil {
		return nil
	}
	_, err := h.WriteTo(sw.bw)
	return err
}

// Write emits one record as a tab-delimited SAM line.
func (sw *SAMWriter) Write(rec *Record) error {
	if !sw.headerDone {
		sw.headerDone = true // allow header-less output
	}
	if _, err := sw.bw.WriteString(rec.QName); err != nil {
		return err
	}
	if err := sw.bw.WriteByte('\t'); err != nil {
		return err
	}
	sw.bw.WriteString(strconv.FormatUint(uint64(rec.Flag), 10))
	sw.bw.WriteByte('\t')
	if rec.RName == "" {
		sw.bw.WriteByte('*')
	} else {
		sw.bw.WriteString(rec.RName)
	}
	sw.bw.WriteByte('\t')
	sw.bw.WriteString(strconv.FormatInt(int64(rec.Pos), 10))
	sw.bw.WriteByte('\t')
	sw.bw.WriteString(strconv.FormatUint(uint64(rec.MapQ), 10))
	sw.bw.WriteByte('\t')
	sw.bw.WriteString(rec.Cigar.String())
	sw.bw.WriteByte('\t')
	if rec.RNext == "" {
		sw.bw.WriteByte('*')
	} else {
		sw.bw.WriteString(rec.RNext)
	}
	sw.bw.WriteByte('\t')
	sw.bw.WriteString(strconv.FormatInt(int64(rec.PNext), 10))
	sw.bw.WriteByte('\t')
	sw.bw.WriteString(strconv.FormatInt(int64(rec.TLen), 10))
	sw.bw.WriteByte('\t')
	if rec.Seq == "" {
		sw.bw.WriteByte('*')
	} else {
		sw.bw.WriteString(rec.Seq)
	}
	sw.bw.WriteByte('\t')
	if len(rec.Qual) == 0 || allQualUnknown(rec.Qual) {
		sw.bw.WriteByte('*')
	} else {
		writePhredASCII(sw.bw, rec.Qual)
	}
	for _, a := range rec.Aux {
		sw.bw.WriteByte('\t')
		sw.bw.WriteString(a.FormatSAM())
	}
	return sw.bw.WriteByte('\n')
}

// Close flushes the writer. It does not close the underlying io.Writer.
func (sw *SAMWriter) Close() error { return sw.bw.Flush() }

// writePhredASCII writes Phred quality bytes back to SAM ASCII-33.
func writePhredASCII(w *bufio.Writer, q []byte) {
	for _, b := range q {
		w.WriteByte(b + 33)
	}
}

// allQualUnknown reports whether every quality byte is 0xff, which by SAM
// convention means "no quality information" and is encoded as "*" in text.
func allQualUnknown(q []byte) bool {
	for _, b := range q {
		if b != 0xff {
			return false
		}
	}
	return true
}
