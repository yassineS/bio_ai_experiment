package cram

import (
	"bytes"
	"strings"
	"testing"
)

// TestRecordReaderHeaderError checks that a CRAM stream whose first
// container is not a SAM-header container is rejected with an error.
func TestRecordReaderHeaderError(t *testing.T) {
	data, ok := loadFixture(t, "dat/test_input_1_a.cram")
	if !ok {
		t.Skipf("samtools submodule not initialised; run `git submodule update --init reference_code/samtools`")
	}
	// Truncating the stream right after the file definition leaves no
	// first container; NewRecordReader must surface that as an error.
	if _, err := NewRecordReader(bytes.NewReader(data[:6])); err == nil {
		t.Error("a stream with no header container should error")
	}
}

// TestRecordReaderMidStreamError corrupts a byte inside a data
// container of the real fixture and checks that the resulting decode
// failure is reported as an error rather than a panic.
func TestRecordReaderMidStreamError(t *testing.T) {
	data, ok := loadFixture(t, "dat/test_input_1_a.cram")
	if !ok {
		t.Skipf("samtools submodule not initialised; run `git submodule update --init reference_code/samtools`")
	}
	// Flip several bytes across the body; for each variant the decoder
	// must either error cleanly or finish without panicking.
	for off := 200; off < len(data)-40; off += 137 {
		corrupt := append([]byte(nil), data...)
		corrupt[off] ^= 0xff
		rr, err := NewRecordReader(bytes.NewReader(corrupt))
		if err != nil {
			continue
		}
		var buf bytes.Buffer
		_ = rr.WriteSAM(&buf) // must not panic.
	}
}

// TestReadAfterEOF checks that Read keeps returning io.EOF once the
// stream is exhausted.
func TestReadAfterEOF(t *testing.T) {
	data, ok := loadFixture(t, "dat/test_input_1_a.cram")
	if !ok {
		t.Skipf("samtools submodule not initialised; run `git submodule update --init reference_code/samtools`")
	}
	rr, err := NewRecordReader(bytes.NewReader(data))
	if err != nil {
		t.Fatalf("NewRecordReader: %v", err)
	}
	if _, err := rr.ReadAll(); err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	for i := 0; i < 3; i++ {
		if _, err := rr.Read(); err == nil {
			t.Error("Read after EOF should keep returning io.EOF")
		}
	}
}

// TestWriteSAMHeaderFirst checks the WriteSAM output begins with the
// embedded SAM header and carries all 15 records.
func TestWriteSAMHeaderFirst(t *testing.T) {
	data, ok := loadFixture(t, "dat/test_input_1_a.cram")
	if !ok {
		t.Skipf("samtools submodule not initialised; run `git submodule update --init reference_code/samtools`")
	}
	rr, err := NewRecordReader(bytes.NewReader(data))
	if err != nil {
		t.Fatalf("NewRecordReader: %v", err)
	}
	var buf bytes.Buffer
	if err := rr.WriteSAM(&buf); err != nil {
		t.Fatalf("WriteSAM: %v", err)
	}
	out := buf.String()
	if !strings.HasPrefix(out, "@HD\tVN:1.4\n") {
		t.Errorf("WriteSAM output should start with the @HD line, got %.20q", out)
	}
	if n := strings.Count(out, "\n"); n != 13+15 {
		t.Errorf("WriteSAM emitted %d lines, want %d (13 header + 15 records)", n, 13+15)
	}
}
