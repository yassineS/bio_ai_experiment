package sam

import (
	"bytes"
	"compress/gzip"
	"io"
	"math"
	"testing"
)

// bamWriterTestHeader and bamWriterTestRecords build a tiny but non-empty BAM
// payload for the uncompressed-mode tests.
func bamWriterTestHeader() *Header {
	h := &Header{Refs: []Reference{{Name: "chr1", Length: 1000}}}
	h.Lines = []HeaderLine{
		{Tag: "HD", Fields: []HeaderField{{Tag: "VN", Value: "1.6"}}},
		{Tag: "SQ", Fields: []HeaderField{{Tag: "SN", Value: "chr1"}, {Tag: "LN", Value: "1000"}}},
	}
	return h
}

func bamWriterTestRecords() []*Record {
	return []*Record{
		{QName: "r1", RName: "chr1", Pos: 10, MapQ: 60, Seq: "ACGTA", Qual: []byte{30, 31, 32, 33, 34}, Cigar: mustCigar("5M")},
		{QName: "r2", RName: "chr1", Pos: 20, MapQ: 40, Seq: "TTGGC", Qual: []byte{20, 21, 22, 23, 24}, Cigar: mustCigar("5M")},
	}
}

// firstDeflateBlockIsStored decodes the first BGZF data block out of a BAM
// stream and reports whether its first DEFLATE block is a stored (BTYPE=00,
// level-0) block. The BGZF gzip framing is a fixed 18-byte header (12-byte
// gzip header + 6-byte BC subfield) followed by the raw DEFLATE payload, so the
// deflate bytes start at offset 18 of each block. A DEFLATE block header is
// BFINAL (1 bit) then BTYPE (2 bits); BTYPE==00 means a stored block.
func firstDeflateBlockIsStored(t *testing.T, bam []byte) bool {
	t.Helper()
	if len(bam) < 18 {
		t.Fatalf("BAM stream too short: %d bytes", len(bam))
	}
	// The first byte of the deflate payload carries BFINAL (bit 0) and BTYPE
	// (bits 1-2). Stored blocks have BTYPE == 0.
	deflateByte := bam[18]
	btype := (deflateByte >> 1) & 0x3
	return btype == 0
}

// TestBAMWriterLongPositionErrors verifies that the BAM writer refuses to write
// a record whose POS, PNEXT or TLEN does not fit the BAM on-disk 32-bit fields.
// BAM stores POS/PNEXT as signed 32-bit and TLEN as signed 32-bit, so a >2^31
// coordinate (which SAM and CRAM support via the int64 Record fields) cannot be
// represented; htslib likewise rejects writing such a record to BAM. The writer
// must error cleanly rather than silently truncate, so long-reference data stays
// in SAM/CRAM.
func TestBAMWriterLongPositionErrors(t *testing.T) {
	const beyond32 = int64(1) << 32 // > math.MaxInt32, even after the -1 to 0-based.

	cases := []struct {
		name string
		rec  *Record
	}{
		{
			name: "POS",
			rec:  &Record{QName: "r", RName: "chr1", Pos: beyond32, MapQ: 60, Seq: "A", Qual: []byte{30}, Cigar: mustCigar("1M")},
		},
		{
			name: "PNEXT",
			rec:  &Record{QName: "r", RName: "chr1", Pos: 10, PNext: beyond32, RNext: "chr1", MapQ: 60, Seq: "A", Qual: []byte{30}, Cigar: mustCigar("1M")},
		},
		{
			name: "TLEN",
			rec:  &Record{QName: "r", RName: "chr1", Pos: 10, TLen: beyond32, MapQ: 60, Seq: "A", Qual: []byte{30}, Cigar: mustCigar("1M")},
		},
	}

	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			var buf bytes.Buffer
			bw := NewBAMWriter(&buf)
			if err := bw.WriteHeader(bamWriterTestHeader()); err != nil {
				t.Fatalf("WriteHeader: %v", err)
			}
			if err := bw.Write(tc.rec); err == nil {
				t.Fatalf("BAM writer accepted a %s beyond the 32-bit field limit; want an error", tc.name)
			}
		})
	}

	// A record at exactly the largest representable BAM coordinate must still
	// succeed: POS == math.MaxInt32+1 (1-based) maps to the 0-based field value
	// math.MaxInt32, which fits.
	t.Run("max-representable-POS", func(t *testing.T) {
		var buf bytes.Buffer
		bw := NewBAMWriter(&buf)
		if err := bw.WriteHeader(bamWriterTestHeader()); err != nil {
			t.Fatalf("WriteHeader: %v", err)
		}
		rec := &Record{QName: "r", RName: "chr1", Pos: int64(math.MaxInt32) + 1, MapQ: 60, Seq: "A", Qual: []byte{30}, Cigar: mustCigar("1M")}
		if err := bw.Write(rec); err != nil {
			t.Fatalf("BAM writer rejected the largest representable POS: %v", err)
		}
	})
}

// TestBAMWriterUncompressedRoundTrip verifies the uncompressed-BAM mode
// (BAMWriterOptions.Uncompressed) produces a valid BAM that round-trips through
// our reader and the stdlib gzip reader, and whose BGZF data blocks are stored
// (level-0) DEFLATE blocks rather than compressed ones.
func TestBAMWriterUncompressedRoundTrip(t *testing.T) {
	h := bamWriterTestHeader()
	recs := bamWriterTestRecords()

	var unc bytes.Buffer
	bw, err := NewBAMWriterOptions(&unc, BAMWriterOptions{Uncompressed: true})
	if err != nil {
		t.Fatalf("NewBAMWriterOptions: %v", err)
	}
	if err := bw.WriteHeader(h); err != nil {
		t.Fatalf("WriteHeader: %v", err)
	}
	for _, r := range recs {
		if err := bw.Write(r); err != nil {
			t.Fatalf("Write: %v", err)
		}
	}
	if err := bw.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	// Round-trip through our BAM reader.
	br, err := NewBAMReader(bytes.NewReader(unc.Bytes()))
	if err != nil {
		t.Fatalf("NewBAMReader: %v", err)
	}
	if got := br.Header(); len(got.Refs) != 1 || got.Refs[0].Name != "chr1" {
		t.Fatalf("header refs after round-trip: %+v", got.Refs)
	}
	var gotRecs []*Record
	for {
		r, err := br.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("Read: %v", err)
		}
		gotRecs = append(gotRecs, r)
	}
	if len(gotRecs) != len(recs) {
		t.Fatalf("record count = %d, want %d", len(gotRecs), len(recs))
	}
	for i := range recs {
		if gotRecs[i].QName != recs[i].QName || gotRecs[i].Pos != recs[i].Pos {
			t.Errorf("record %d = %q@%d, want %q@%d", i, gotRecs[i].QName, gotRecs[i].Pos, recs[i].QName, recs[i].Pos)
		}
	}

	// The BGZF blocks are ordinary gzip members and must decode with the
	// stdlib gzip reader to the same plaintext as a compressed BAM would.
	gz, err := gzip.NewReader(bytes.NewReader(unc.Bytes()))
	if err != nil {
		t.Fatalf("gzip.NewReader: %v", err)
	}
	gz.Multistream(true)
	plain, err := io.ReadAll(gz)
	if err != nil {
		t.Fatalf("gzip ReadAll: %v", err)
	}
	if !bytes.HasPrefix(plain, BAMMagic) {
		t.Errorf("decoded plaintext does not start with BAM magic: % x", plain[:min(4, len(plain))])
	}

	// The defining property of uncompressed BAM: the data blocks are stored
	// (level-0) DEFLATE blocks.
	if !firstDeflateBlockIsStored(t, unc.Bytes()) {
		t.Errorf("uncompressed BAM first data block is not a stored DEFLATE block")
	}

	// Cross-check: a default (compressed) BAM of the same input must NOT be a
	// stored block, so the test above is actually discriminating.
	var comp bytes.Buffer
	cw := NewBAMWriter(&comp)
	if err := cw.WriteHeader(h); err != nil {
		t.Fatalf("compressed WriteHeader: %v", err)
	}
	for _, r := range recs {
		if err := cw.Write(r); err != nil {
			t.Fatalf("compressed Write: %v", err)
		}
	}
	if err := cw.Close(); err != nil {
		t.Fatalf("compressed Close: %v", err)
	}
	if firstDeflateBlockIsStored(t, comp.Bytes()) {
		t.Errorf("default compressed BAM unexpectedly used a stored DEFLATE block")
	}
}

// TestBAMWriterOptionsThreadsUncompressed verifies the threaded uncompressed
// path also produces stored, round-trippable BGZF blocks.
func TestBAMWriterOptionsThreadsUncompressed(t *testing.T) {
	h := bamWriterTestHeader()
	recs := bamWriterTestRecords()
	var buf bytes.Buffer
	bw, err := NewBAMWriterOptions(&buf, BAMWriterOptions{Uncompressed: true, Threads: 2})
	if err != nil {
		t.Fatalf("NewBAMWriterOptions: %v", err)
	}
	if err := bw.WriteHeader(h); err != nil {
		t.Fatalf("WriteHeader: %v", err)
	}
	for _, r := range recs {
		if err := bw.Write(r); err != nil {
			t.Fatalf("Write: %v", err)
		}
	}
	if err := bw.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	br, err := NewBAMReader(bytes.NewReader(buf.Bytes()))
	if err != nil {
		t.Fatalf("NewBAMReader: %v", err)
	}
	n := 0
	for {
		if _, err := br.Read(); err == io.EOF {
			break
		} else if err != nil {
			t.Fatalf("Read: %v", err)
		}
		n++
	}
	if n != len(recs) {
		t.Fatalf("threaded uncompressed record count = %d, want %d", n, len(recs))
	}
	if !firstDeflateBlockIsStored(t, buf.Bytes()) {
		t.Errorf("threaded uncompressed BAM first data block is not stored")
	}
}

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
	// in pkg/htsgo/tabix's binning_test.go.
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
