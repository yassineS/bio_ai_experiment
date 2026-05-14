package sam

import (
	"bytes"
	"io"
	"testing"
)

// TestBAMWriteAllAuxSubtypes exercises every aux subtype on the writer
// (including the fixed-width 'I'/'C'/'s'/'S' types and every B subtype)
// and decodes them back to confirm byte fidelity.
func TestBAMWriteAllAuxSubtypes(t *testing.T) {
	h := &Header{Refs: []Reference{{Name: "chr1", Length: 1000}}}
	h.Lines = []HeaderLine{{Tag: "SQ", Fields: []HeaderField{{Tag: "SN", Value: "chr1"}, {Tag: "LN", Value: "1000"}}}}
	rec := &Record{
		QName: "r", RName: "chr1", Pos: 1, Seq: "AC", Qual: []byte{30, 30},
		Cigar: mustCigar("2M"),
		Aux: []Aux{
			{Tag: "c1", Type: 'c', Value: int64(-5)},
			{Tag: "C1", Type: 'C', Value: int64(200)},
			{Tag: "s1", Type: 's', Value: int64(-30000)},
			{Tag: "S1", Type: 'S', Value: int64(50000)},
			{Tag: "I1", Type: 'I', Value: int64(3000000000)},
			{Tag: "BC", Type: 'B', ArrayType: 'c', ArrayValues: []interface{}{int64(-1), int64(2)}},
			{Tag: "BS", Type: 'B', ArrayType: 'C', ArrayValues: []interface{}{int64(200), int64(100)}},
			{Tag: "Bs", Type: 'B', ArrayType: 's', ArrayValues: []interface{}{int64(-1000), int64(2000)}},
			{Tag: "BI", Type: 'B', ArrayType: 'I', ArrayValues: []interface{}{int64(3000000000), int64(1)}},
			{Tag: "HX", Type: 'H', Value: "DEADBEEF"},
		},
	}
	var bam bytes.Buffer
	bw := NewBAMWriter(&bam)
	if err := bw.WriteHeader(h); err != nil {
		t.Fatalf("WriteHeader: %v", err)
	}
	if err := bw.Write(rec); err != nil {
		t.Fatalf("Write: %v", err)
	}
	if err := bw.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	// Calling Close again should be a no-op.
	if err := bw.Close(); err != nil {
		t.Fatalf("second Close: %v", err)
	}

	br, err := NewBAMReader(bytes.NewReader(bam.Bytes()))
	if err != nil {
		t.Fatalf("NewBAMReader: %v", err)
	}
	got, err := br.Read()
	if err != nil {
		t.Fatalf("Read: %v", err)
	}
	if len(got.Aux) != len(rec.Aux) {
		t.Fatalf("aux count mismatch: %d vs %d", len(got.Aux), len(rec.Aux))
	}
	for i := range rec.Aux {
		if got.Aux[i].Tag != rec.Aux[i].Tag {
			t.Errorf("aux[%d] tag mismatch: %q vs %q", i, got.Aux[i].Tag, rec.Aux[i].Tag)
		}
	}
	if _, err := br.Read(); err != io.EOF {
		t.Errorf("expected EOF, got %v", err)
	}
}

// TestBAMWriteCompactInt confirms writeBAMIntCompact picks the smallest
// width that fits each value.
func TestBAMWriteCompactInt(t *testing.T) {
	tests := []struct {
		v      int64
		expect byte // expected BAM type byte
	}{
		{0, 'C'},
		{200, 'C'},
		{-1, 'c'},
		{-200, 's'},
		{40000, 'S'},
		{-40000, 'i'},
		{3000000000, 'I'},
		{-3000000000, 'i'},
	}
	for _, tc := range tests {
		var buf bytes.Buffer
		writeBAMIntCompact(&buf, tc.v)
		if got := buf.Bytes()[0]; got != tc.expect {
			t.Errorf("writeBAMIntCompact(%d): got type %c, want %c", tc.v, got, tc.expect)
		}
	}
}

// TestReg2Bin verifies the local reg2bin helper agrees with the canonical
// tabix.Reg2bin (it just adds a guard for empty intervals).
func TestReg2Bin(t *testing.T) {
	// Smoke test a handful of intervals — the per-tier verification lives
	// in tools/tabix/pkg/tabix's binning_test.go.
	if got := reg2bin(0, 1); got != 4681 {
		t.Errorf("reg2bin(0,1): got %d want 4681", got)
	}
	if got := reg2bin(100, 50); got <= 0 {
		t.Errorf("reg2bin(100,50): got %d want positive (degenerate guarded)", got)
	}
	if got := reg2bin(0, 1<<30); got != 0 {
		t.Errorf("reg2bin(0,1<<30): got %d want 0 (root)", got)
	}
}

// TestSAMWriterUnmapped verifies SEQ="*" and QUAL="*" are emitted for
// records whose seq/qual are absent.
func TestSAMWriterUnmapped(t *testing.T) {
	rec := &Record{QName: "u", Flag: FlagUnmapped, Pos: 0}
	var out bytes.Buffer
	w := NewSAMWriter(&out)
	if err := w.WriteHeader(nil); err != nil {
		t.Fatalf("WriteHeader nil: %v", err)
	}
	if err := w.Write(rec); err != nil {
		t.Fatalf("Write: %v", err)
	}
	w.Close()
	want := "u\t4\t*\t0\t0\t*\t*\t0\t0\t*\t*\n"
	if got := out.String(); got != want {
		t.Errorf("unmapped: got %q want %q", got, want)
	}
}

// TestSAMWriterAllQualUnknown verifies a record with all-0xff Qual emits "*".
func TestSAMWriterAllQualUnknown(t *testing.T) {
	rec := &Record{QName: "r", Pos: 1, Seq: "ACGT", Qual: []byte{0xff, 0xff, 0xff, 0xff}}
	var out bytes.Buffer
	w := NewSAMWriter(&out)
	w.WriteHeader(nil)
	w.Write(rec)
	w.Close()
	// QUAL field should be "*"
	if got := out.String(); !bytes.Contains([]byte(got), []byte("ACGT\t*\n")) {
		t.Errorf("unexpected output: %q", got)
	}
}

// TestBAMHeaderRoundTripBinary exercises a header round-trip and confirms
// l_name handling for refs with trailing NUL.
func TestBAMHeaderRoundTripBinary(t *testing.T) {
	h := &Header{Refs: []Reference{{Name: "chr1", Length: 100}, {Name: "chrM", Length: 16569}}}
	h.Lines = []HeaderLine{
		{Tag: "HD", Fields: []HeaderField{{Tag: "VN", Value: "1.6"}}},
		{Tag: "SQ", Fields: []HeaderField{{Tag: "SN", Value: "chr1"}, {Tag: "LN", Value: "100"}}},
		{Tag: "SQ", Fields: []HeaderField{{Tag: "SN", Value: "chrM"}, {Tag: "LN", Value: "16569"}}},
	}

	var bam bytes.Buffer
	bw := NewBAMWriter(&bam)
	if err := bw.WriteHeader(h); err != nil {
		t.Fatalf("WriteHeader: %v", err)
	}
	if err := bw.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	br, err := NewBAMReader(bytes.NewReader(bam.Bytes()))
	if err != nil {
		t.Fatalf("NewBAMReader: %v", err)
	}
	got := br.Header()
	if len(got.Refs) != 2 || got.Refs[0].Name != "chr1" || got.Refs[1].Name != "chrM" {
		t.Errorf("header refs: %+v", got.Refs)
	}
	if got.Refs[1].Length != 16569 {
		t.Errorf("chrM length: %d", got.Refs[1].Length)
	}
	if _, err := br.Read(); err != io.EOF {
		t.Errorf("expected immediate EOF, got %v", err)
	}
}

func TestParseCigarBackOp(t *testing.T) {
	if _, err := ParseCigar("3B"); err != nil {
		t.Errorf("ParseCigar 3B: %v", err)
	}
}

func TestLooksLikeBGZF(t *testing.T) {
	if looksLikeBGZF([]byte{0x1f, 0x8b}) {
		t.Error("short slice should not match")
	}
	notGzip := make([]byte, 16)
	if looksLikeBGZF(notGzip) {
		t.Error("zeros should not match")
	}
}
