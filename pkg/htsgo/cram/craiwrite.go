package cram

import (
	"bufio"
	"compress/gzip"
	"io"
	"strconv"
)

// WriteCRAI serialises index entries as a CRAM index (.crai) to w. The
// .crai is a gzip-compressed, tab-separated text file; each line carries
// the six integers ReadCRAI parses: reference id, alignment start,
// alignment span, the container's absolute byte offset, the slice's byte
// offset within the container and the slice's size in bytes.
//
// WriteCRAI is the inverse of ReadCRAI: ReadCRAI(WriteCRAI(entries))
// recovers the same entries. The entries are written in the order given;
// callers that want a seek-ordered index should sort by ContainerOffset
// then SliceOffset first (BuildCRAI already emits entries in file order).
func WriteCRAI(w io.Writer, entries []CRAIEntry) error {
	gz := gzip.NewWriter(w)
	bw := bufio.NewWriter(gz)
	for _, e := range entries {
		if err := writeCRAILine(bw, e); err != nil {
			return wrapf(err, "writing the .crai stream")
		}
	}
	if err := bw.Flush(); err != nil {
		return wrapf(err, "flushing the .crai stream")
	}
	if err := gz.Close(); err != nil {
		return wrapf(err, "finishing the .crai gzip stream")
	}
	return nil
}

// writeCRAILine writes one tab-separated, newline-terminated .crai record.
func writeCRAILine(bw *bufio.Writer, e CRAIEntry) error {
	fields := [6]int64{
		int64(e.RefID),
		int64(e.AlignmentStart),
		int64(e.AlignmentSpan),
		e.ContainerOffset,
		e.SliceOffset,
		e.SliceSize,
	}
	for i, v := range fields {
		if i > 0 {
			if err := bw.WriteByte('\t'); err != nil {
				return err
			}
		}
		if _, err := bw.WriteString(strconv.FormatInt(v, 10)); err != nil {
			return err
		}
	}
	return bw.WriteByte('\n')
}
