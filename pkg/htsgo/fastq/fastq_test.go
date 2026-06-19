package fastq

import (
	"io"
	"strings"
	"testing"
)

const sampleFASTQ = "@read1 desc one\nACGTACGT\n+\nIIIIIIII\n" +
	"@read2\nTTGGCCAA\n+\n!!!!!!!!\n"

func TestReadParsesRecords(t *testing.T) {
	r := NewReader(strings.NewReader(sampleFASTQ), Phred33)

	rec, err := r.Read()
	if err != nil {
		t.Fatalf("Read: %v", err)
	}
	if rec.ID != "read1" {
		t.Errorf("ID = %q, want read1", rec.ID)
	}
	if rec.Description != "read1 desc one" {
		t.Errorf("Description = %q, want %q", rec.Description, "read1 desc one")
	}
	if string(rec.Sequence) != "ACGTACGT" || string(rec.Quality) != "IIIIIIII" {
		t.Errorf("seq/qual = %q/%q", rec.Sequence, rec.Quality)
	}

	rec2, err := r.Read()
	if err != nil {
		t.Fatalf("Read 2: %v", err)
	}
	if rec2.ID != "read2" || rec2.Description != "read2" {
		t.Errorf("rec2 ID/Description = %q/%q", rec2.ID, rec2.Description)
	}

	if _, err := r.Read(); err != io.EOF {
		t.Errorf("expected EOF, got %v", err)
	}
}

// TestReadIntoReusesBuffersAndMatchesRead checks that ReadInto yields the same
// field values as Read while reusing the record's Sequence/Quality backing
// arrays across calls (no growth when lengths are equal).
func TestReadIntoReusesBuffersAndMatchesRead(t *testing.T) {
	r := NewReader(strings.NewReader(sampleFASTQ), Phred33)
	var rec Record

	if err := r.ReadInto(&rec); err != nil {
		t.Fatalf("ReadInto: %v", err)
	}
	seqArr := &rec.Sequence[0]
	if rec.ID != "read1" || rec.Description != "read1 desc one" {
		t.Errorf("ID/Description = %q/%q", rec.ID, rec.Description)
	}
	if string(rec.Sequence) != "ACGTACGT" || string(rec.Quality) != "IIIIIIII" {
		t.Errorf("seq/qual = %q/%q", rec.Sequence, rec.Quality)
	}

	if err := r.ReadInto(&rec); err != nil {
		t.Fatalf("ReadInto 2: %v", err)
	}
	if string(rec.Sequence) != "TTGGCCAA" {
		t.Errorf("seq2 = %q, want TTGGCCAA", rec.Sequence)
	}
	// Equal-length record must reuse the same backing array.
	if &rec.Sequence[0] != seqArr {
		t.Errorf("ReadInto reallocated the Sequence buffer for an equal-length record")
	}

	if err := r.ReadInto(&rec); err != io.EOF {
		t.Errorf("expected EOF, got %v", err)
	}
}

func TestReadIntoRejectsBadHeader(t *testing.T) {
	r := NewReader(strings.NewReader("not-a-header\nACGT\n+\nIIII\n"), Phred33)
	if err := r.ReadInto(&Record{}); err == nil {
		t.Fatal("expected error for missing '@' header")
	}
}
