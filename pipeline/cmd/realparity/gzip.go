package main

import (
	"compress/gzip"
	"io"
)

// newGzipReader wraps r in a multistream gzip reader, so a bgzipped (BGZF) VCF —
// which is a concatenation of independent gzip members — decompresses fully.
func newGzipReader(r io.Reader) (*gzip.Reader, error) {
	gr, err := gzip.NewReader(r)
	if err != nil {
		return nil, err
	}
	gr.Multistream(true)
	return gr, nil
}
