package cram

import (
	"bytes"
	"fmt"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// TestEncodeThreadsByteIdentical verifies the async encode pipeline is
// deterministic: a file encoded with any number of worker threads is
// byte-for-byte identical to the synchronous (single-thread) encoding. Records
// span several containers (recordsPerSlice is lowered) so the ordered-writer
// path is exercised, not just a single container.
func TestEncodeThreadsByteIdentical(t *testing.T) {
	h := writerTestHeader()
	const n = 95
	records := make([]*sam.Record, 0, n)
	for i := 0; i < n; i++ {
		records = append(records, mkRec(
			fmt.Sprintf("r%04d", i), "chr1", int32(100+i), "20M",
			"ACGTACGTACGTACGTACGT"))
	}

	encode := func(threads int) []byte {
		var buf bytes.Buffer
		rw, err := NewRecordWriterOpts(&buf, h, WriterOptions{EncodeThreads: threads})
		if err != nil {
			t.Fatalf("threads=%d: NewRecordWriterOpts: %v", threads, err)
		}
		rw.recordsPerSlice = 10 // force ~10 containers
		for _, r := range records {
			if err := rw.Write(r); err != nil {
				t.Fatalf("threads=%d: Write: %v", threads, err)
			}
		}
		if err := rw.Close(); err != nil {
			t.Fatalf("threads=%d: Close: %v", threads, err)
		}
		return buf.Bytes()
	}

	want := encode(1)
	for _, threads := range []int{2, 4, 8} {
		got := encode(threads)
		if !bytes.Equal(got, want) {
			t.Errorf("EncodeThreads=%d produced %d bytes, differs from the %d-byte serial encoding",
				threads, len(got), len(want))
		}
	}

	// And the parallel encoding must still round-trip.
	rr, err := NewRecordReader(bytes.NewReader(want))
	if err != nil {
		t.Fatalf("NewRecordReader: %v", err)
	}
	out, err := rr.ReadAll()
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if len(out) != n {
		t.Fatalf("round-trip decoded %d records, want %d", len(out), n)
	}
}
