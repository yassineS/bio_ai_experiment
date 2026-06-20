package sam

import (
	"bytes"
	"io"
	"testing"
)

// buildSampleBAM round-trips sampleSAM through the BAM writer and returns the
// encoded BAM bytes, for tests that need a small in-memory BAM stream.
func buildSampleBAM(t *testing.T) []byte {
	t.Helper()
	r, err := NewSAMReader(bytes.NewReader([]byte(sampleSAM)))
	if err != nil {
		t.Fatalf("NewSAMReader: %v", err)
	}
	var bam bytes.Buffer
	bw := NewBAMWriter(&bam)
	if err := bw.WriteHeader(r.Header()); err != nil {
		t.Fatalf("WriteHeader: %v", err)
	}
	for {
		rec, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("SAM Read: %v", err)
		}
		if err := bw.Write(rec); err != nil {
			t.Fatalf("BAM Write: %v", err)
		}
	}
	if err := bw.Close(); err != nil {
		t.Fatalf("BAM Close: %v", err)
	}
	return bam.Bytes()
}

// TestReadShallowIntoMatchesFullDecode asserts that the shallow decode path
// yields the same fixed-prefix fields (RName/Pos/MapQ/Flag/RNext/PNext/TLen)
// as a full decode on the same BAM, while leaving the variable-length fields
// (QName/Cigar/Seq/Qual/Aux) empty.
func TestReadShallowIntoMatchesFullDecode(t *testing.T) {
	bam := buildSampleBAM(t)

	full, err := NewBAMReader(bytes.NewReader(bam))
	if err != nil {
		t.Fatalf("NewBAMReader (full): %v", err)
	}
	shallow, err := NewBAMReader(bytes.NewReader(bam))
	if err != nil {
		t.Fatalf("NewBAMReader (shallow): %v", err)
	}

	var sRec Record
	n := 0
	for {
		fRec, ferr := full.Read()
		serr := shallow.ReadShallowInto(&sRec)
		if ferr == io.EOF {
			if serr != io.EOF {
				t.Fatalf("record %d: full hit EOF but shallow err=%v", n, serr)
			}
			break
		}
		if ferr != nil {
			t.Fatalf("record %d: full Read: %v", n, ferr)
		}
		if serr != nil {
			t.Fatalf("record %d: ReadShallowInto: %v", n, serr)
		}

		if sRec.RName != fRec.RName {
			t.Errorf("record %d: RName shallow=%q full=%q", n, sRec.RName, fRec.RName)
		}
		if sRec.Pos != fRec.Pos {
			t.Errorf("record %d: Pos shallow=%d full=%d", n, sRec.Pos, fRec.Pos)
		}
		if sRec.MapQ != fRec.MapQ {
			t.Errorf("record %d: MapQ shallow=%d full=%d", n, sRec.MapQ, fRec.MapQ)
		}
		if sRec.Flag != fRec.Flag {
			t.Errorf("record %d: Flag shallow=%d full=%d", n, sRec.Flag, fRec.Flag)
		}
		if sRec.RNext != fRec.RNext {
			t.Errorf("record %d: RNext shallow=%q full=%q", n, sRec.RNext, fRec.RNext)
		}
		if sRec.PNext != fRec.PNext {
			t.Errorf("record %d: PNext shallow=%d full=%d", n, sRec.PNext, fRec.PNext)
		}
		if sRec.TLen != fRec.TLen {
			t.Errorf("record %d: TLen shallow=%d full=%d", n, sRec.TLen, fRec.TLen)
		}

		// The variable-length region must be left empty by the shallow path.
		if sRec.QName != "" {
			t.Errorf("record %d: shallow QName not empty: %q", n, sRec.QName)
		}
		if len(sRec.Cigar) != 0 {
			t.Errorf("record %d: shallow Cigar not empty: %v", n, sRec.Cigar)
		}
		if sRec.Seq != "" {
			t.Errorf("record %d: shallow Seq not empty: %q", n, sRec.Seq)
		}
		if len(sRec.Qual) != 0 {
			t.Errorf("record %d: shallow Qual not empty: %v", n, sRec.Qual)
		}
		if len(sRec.Aux) != 0 {
			t.Errorf("record %d: shallow Aux not empty: %v", n, sRec.Aux)
		}
		n++
	}
	if n == 0 {
		t.Fatal("decoded no records")
	}
}

// TestReadShallowIntoResetsStaleVariableFields asserts that reusing one record
// across a full decode followed by shallow decodes clears the variable-length
// fields left over from the full decode (no stale data leaks).
func TestReadShallowIntoResetsStaleVariableFields(t *testing.T) {
	bam := buildSampleBAM(t)
	br, err := NewBAMReader(bytes.NewReader(bam))
	if err != nil {
		t.Fatalf("NewBAMReader: %v", err)
	}
	var rec Record
	// First do a full decode to populate the variable-length fields.
	if err := br.ReadInto(&rec); err != nil {
		t.Fatalf("ReadInto: %v", err)
	}
	if rec.QName == "" {
		t.Fatal("expected first full decode to populate QName")
	}
	// Now shallow-decode into the same record; the variable fields must clear.
	if err := br.ReadShallowInto(&rec); err != nil {
		t.Fatalf("ReadShallowInto: %v", err)
	}
	if rec.QName != "" {
		t.Errorf("stale QName after shallow decode: %q", rec.QName)
	}
	if len(rec.Cigar) != 0 {
		t.Errorf("stale Cigar after shallow decode: %v", rec.Cigar)
	}
	if rec.Seq != "" {
		t.Errorf("stale Seq after shallow decode: %q", rec.Seq)
	}
	if len(rec.Qual) != 0 {
		t.Errorf("stale Qual after shallow decode: %v", rec.Qual)
	}
	if len(rec.Aux) != 0 {
		t.Errorf("stale Aux after shallow decode: %v", rec.Aux)
	}
}
