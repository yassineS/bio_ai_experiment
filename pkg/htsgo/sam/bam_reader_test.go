package sam

import (
	"bytes"
	"encoding/binary"
	"io"
	"strings"
	"testing"

	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
)

// buildHexBAM hand-crafts a minimal BAM stream containing a single
// alignment record. We assemble the BAM-payload bytes by hand, then wrap
// them in BGZF using our bgzip.Writer. This validates the BAM parser
// independently of the BAM writer.
func buildHexBAM(t *testing.T) []byte {
	t.Helper()
	var payload bytes.Buffer

	// Header.
	payload.WriteString("BAM\x01")
	headerText := "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n"
	binary.Write(&payload, binary.LittleEndian, int32(len(headerText)))
	payload.WriteString(headerText)
	binary.Write(&payload, binary.LittleEndian, int32(1)) // n_ref
	name := []byte("chr1\x00")
	binary.Write(&payload, binary.LittleEndian, int32(len(name)))
	payload.Write(name)
	binary.Write(&payload, binary.LittleEndian, int32(1000)) // l_ref

	// Record body.
	var body bytes.Buffer
	binary.Write(&body, binary.LittleEndian, int32(0))  // refID
	binary.Write(&body, binary.LittleEndian, int32(99)) // pos (0-based = 99 → SAM 100)
	body.WriteByte(byte(len("readX") + 1))              // l_read_name (incl NUL)
	body.WriteByte(60)                                  // mapq
	binary.Write(&body, binary.LittleEndian, uint16(0)) // bin (decoded value ignored on read)
	binary.Write(&body, binary.LittleEndian, uint16(1)) // n_cigar_op
	binary.Write(&body, binary.LittleEndian, uint16(0)) // flag = 0
	binary.Write(&body, binary.LittleEndian, int32(4))  // l_seq
	binary.Write(&body, binary.LittleEndian, int32(-1)) // next_refID
	binary.Write(&body, binary.LittleEndian, int32(-1)) // next_pos
	binary.Write(&body, binary.LittleEndian, int32(0))  // tlen
	body.WriteString("readX\x00")
	binary.Write(&body, binary.LittleEndian, uint32(4<<4|CigarMatch)) // 4M
	// SEQ ACGT packed: A=1,C=2,G=4,T=8 → 0x12, 0x48.
	body.WriteByte(0x12)
	body.WriteByte(0x48)
	body.Write([]byte{30, 30, 30, 30}) // qual
	// AUX: NM:i:1 (encoded as 'C' since 1 fits)
	body.WriteString("NMC")
	body.WriteByte(1)

	// Append block_size + body to header.
	binary.Write(&payload, binary.LittleEndian, int32(body.Len()))
	payload.Write(body.Bytes())

	// Wrap in BGZF.
	var bgz bytes.Buffer
	w := bgzip.NewWriter(&bgz)
	w.Write(payload.Bytes())
	if err := w.Close(); err != nil {
		t.Fatalf("bgzip close: %v", err)
	}
	return bgz.Bytes()
}

func TestBAMReaderHexCrafted(t *testing.T) {
	raw := buildHexBAM(t)
	r, err := NewBAMReader(bytes.NewReader(raw))
	if err != nil {
		t.Fatalf("NewBAMReader: %v", err)
	}
	hdr := r.Header()
	if len(hdr.Refs) != 1 || hdr.Refs[0].Name != "chr1" || hdr.Refs[0].Length != 1000 {
		t.Fatalf("bad header refs: %+v", hdr.Refs)
	}
	rec, err := r.Read()
	if err != nil {
		t.Fatalf("Read: %v", err)
	}
	if rec.QName != "readX" {
		t.Errorf("QName: %q", rec.QName)
	}
	if rec.Pos != 100 {
		t.Errorf("Pos: got %d want 100 (1-based)", rec.Pos)
	}
	if rec.MapQ != 60 {
		t.Errorf("MapQ: %d", rec.MapQ)
	}
	if rec.Cigar.String() != "4M" {
		t.Errorf("Cigar: %q", rec.Cigar.String())
	}
	if rec.Seq != "ACGT" {
		t.Errorf("Seq: %q", rec.Seq)
	}
	if len(rec.Qual) != 4 || rec.Qual[0] != 30 {
		t.Errorf("Qual: %v", rec.Qual)
	}
	if len(rec.Aux) != 1 || rec.Aux[0].Tag != "NM" || rec.Aux[0].Type != 'C' {
		t.Errorf("Aux: %+v", rec.Aux)
	}
	if v, _ := rec.Aux[0].Int(); v != 1 {
		t.Errorf("Aux NM value: %d", v)
	}
	if _, err := r.Read(); err != io.EOF {
		t.Errorf("expected EOF after single record, got %v", err)
	}
	if err := r.Close(); err != nil {
		t.Errorf("Close: %v", err)
	}
}

func TestBAMReaderBadMagic(t *testing.T) {
	// Build a valid BGZF stream but with wrong magic.
	var bgz bytes.Buffer
	w := bgzip.NewWriter(&bgz)
	w.Write([]byte("NOPE\x00\x00\x00\x00"))
	w.Close()
	_, err := NewBAMReader(bytes.NewReader(bgz.Bytes()))
	if err == nil {
		t.Error("expected error for bad magic")
	}
}

func TestSamBamRoundTrip(t *testing.T) {
	r, err := NewSAMReader(strings.NewReader(sampleSAM))
	if err != nil {
		t.Fatalf("NewSAMReader: %v", err)
	}

	// Read all SAM records.
	var samRecs []*Record
	for {
		rec, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("SAM Read: %v", err)
		}
		samRecs = append(samRecs, rec)
	}

	// Write them out as BAM.
	var bam bytes.Buffer
	bw := NewBAMWriter(&bam)
	if err := bw.WriteHeader(r.Header()); err != nil {
		t.Fatalf("BAM WriteHeader: %v", err)
	}
	for _, rec := range samRecs {
		if err := bw.Write(rec); err != nil {
			t.Fatalf("BAM Write: %v", err)
		}
	}
	if err := bw.Close(); err != nil {
		t.Fatalf("BAM Close: %v", err)
	}

	// Read the BAM back.
	br, err := NewBAMReader(bytes.NewReader(bam.Bytes()))
	if err != nil {
		t.Fatalf("NewBAMReader: %v", err)
	}
	var bamRecs []*Record
	for {
		rec, err := br.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("BAM Read: %v", err)
		}
		bamRecs = append(bamRecs, rec)
	}
	if len(bamRecs) != len(samRecs) {
		t.Fatalf("expected %d records back, got %d", len(samRecs), len(bamRecs))
	}
	for i := range samRecs {
		s := samRecs[i]
		b := bamRecs[i]
		if s.QName != b.QName || s.Flag != b.Flag || s.RName != b.RName || s.Pos != b.Pos {
			t.Errorf("record %d header mismatch: SAM=%+v BAM=%+v", i, s, b)
		}
		if s.Cigar.String() != b.Cigar.String() {
			t.Errorf("record %d cigar mismatch: %q vs %q", i, s.Cigar.String(), b.Cigar.String())
		}
		if s.Seq != b.Seq {
			t.Errorf("record %d seq: %q vs %q", i, s.Seq, b.Seq)
		}
		if !bytes.Equal(s.Qual, b.Qual) {
			t.Errorf("record %d qual: %v vs %v", i, s.Qual, b.Qual)
		}
		if s.MapQ != b.MapQ || s.PNext != b.PNext || s.TLen != b.TLen {
			t.Errorf("record %d mapq/pnext/tlen mismatch", i)
		}
		// Aux contents should match by tag (BAM might compress 'i' to 'C', that's OK).
		if len(s.Aux) != len(b.Aux) {
			t.Errorf("record %d aux count: %d vs %d", i, len(s.Aux), len(b.Aux))
			continue
		}
		for j := range s.Aux {
			sa := s.Aux[j]
			ba := b.Aux[j]
			if sa.Tag != ba.Tag {
				t.Errorf("record %d aux %d tag: %q vs %q", i, j, sa.Tag, ba.Tag)
			}
		}
	}
}

func TestBAMWriteUnknownRef(t *testing.T) {
	var bam bytes.Buffer
	bw := NewBAMWriter(&bam)
	h := &Header{Refs: []Reference{{Name: "chr1", Length: 100}}}
	h.Lines = []HeaderLine{{Tag: "SQ", Fields: []HeaderField{{Tag: "SN", Value: "chr1"}, {Tag: "LN", Value: "100"}}}}
	if err := bw.WriteHeader(h); err != nil {
		t.Fatalf("WriteHeader: %v", err)
	}
	rec := &Record{QName: "r", RName: "chrUNKNOWN", Pos: 1, Seq: "A"}
	if err := bw.Write(rec); err == nil {
		t.Error("expected error for unknown RName")
	}
}

func TestBAMWriteRoundTripAuxTypes(t *testing.T) {
	// Build a record exercising each aux type, write to BAM, read back, verify.
	h := &Header{Refs: []Reference{{Name: "chr1", Length: 1000}}}
	h.Lines = []HeaderLine{{Tag: "SQ", Fields: []HeaderField{{Tag: "SN", Value: "chr1"}, {Tag: "LN", Value: "1000"}}}}
	rec := &Record{
		QName: "r1",
		Flag:  0,
		RName: "chr1",
		Pos:   1,
		MapQ:  60,
		Cigar: mustCigar("4M"),
		Seq:   "ACGT",
		Qual:  []byte{30, 30, 30, 30},
		Aux: []Aux{
			{Tag: "NM", Type: 'i', Value: int64(0)},
			{Tag: "AS", Type: 'i', Value: int64(-12)},
			{Tag: "BI", Type: 'i', Value: int64(70000)},
			{Tag: "XF", Type: 'f', Value: 1.25},
			{Tag: "RG", Type: 'Z', Value: "rg1"},
			{Tag: "XA", Type: 'A', Value: "U"},
			{Tag: "BB", Type: 'B', ArrayType: 'i', ArrayValues: []interface{}{int64(1), int64(2), int64(3)}},
			{Tag: "BF", Type: 'B', ArrayType: 'f', ArrayValues: []interface{}{0.5, 1.5}},
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

	br, err := NewBAMReader(bytes.NewReader(bam.Bytes()))
	if err != nil {
		t.Fatalf("NewBAMReader: %v", err)
	}
	got, err := br.Read()
	if err != nil {
		t.Fatalf("Read: %v", err)
	}
	if len(got.Aux) != len(rec.Aux) {
		t.Fatalf("aux count: got %d, want %d", len(got.Aux), len(rec.Aux))
	}
	for i, want := range rec.Aux {
		if got.Aux[i].Tag != want.Tag {
			t.Errorf("aux[%d] tag: got %q want %q", i, got.Aux[i].Tag, want.Tag)
		}
		switch want.Type {
		case 'A':
			if got.Aux[i].Value != want.Value {
				t.Errorf("aux[%d] A value: got %v want %v", i, got.Aux[i].Value, want.Value)
			}
		case 'Z':
			if got.Aux[i].Value != want.Value {
				t.Errorf("aux[%d] Z value: got %v want %v", i, got.Aux[i].Value, want.Value)
			}
		case 'f':
			if gv, _ := got.Aux[i].Value.(float64); gv != want.Value {
				t.Errorf("aux[%d] f value: got %v want %v", i, gv, want.Value)
			}
		case 'i':
			wv := want.Value.(int64)
			gv, _ := got.Aux[i].Int()
			if gv != wv {
				t.Errorf("aux[%d] int value: got %d want %d", i, gv, wv)
			}
		case 'B':
			if got.Aux[i].ArrayType != want.ArrayType {
				t.Errorf("aux[%d] B subtype: got %c want %c", i, got.Aux[i].ArrayType, want.ArrayType)
			}
			if len(got.Aux[i].ArrayValues) != len(want.ArrayValues) {
				t.Errorf("aux[%d] B count: got %d want %d", i, len(got.Aux[i].ArrayValues), len(want.ArrayValues))
			}
		}
	}
}

func mustCigar(s string) Cigar {
	c, err := ParseCigar(s)
	if err != nil {
		panic(err)
	}
	return c
}

// TestNewReaderRawBAM exercises the code path where iohelper has already
// decompressed the BGZF layer and NewReader sees a raw "BAM\1" body.
func TestNewReaderRawBAM(t *testing.T) {
	raw := buildHexBAM(t)
	// Decompress the BGZF layer ourselves to simulate iohelper's transparent
	// decompression.
	bgz, err := bgzip.NewReader(bytes.NewReader(raw))
	if err != nil {
		t.Fatalf("bgzip.NewReader: %v", err)
	}
	defer bgz.Close()
	decoded, err := io.ReadAll(bgz)
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	r, err := NewReader(bytes.NewReader(decoded))
	if err != nil {
		t.Fatalf("NewReader on raw BAM body: %v", err)
	}
	if _, ok := r.(*BAMReader); !ok {
		t.Fatalf("expected BAMReader for raw BAM body, got %T", r)
	}
	rec, err := r.(*BAMReader).Read()
	if err != nil {
		t.Fatalf("Read: %v", err)
	}
	if rec.QName != "readX" {
		t.Errorf("QName: %q", rec.QName)
	}
	// Close on a raw BAM reader is a no-op.
	if err := r.(*BAMReader).Close(); err != nil {
		t.Errorf("Close: %v", err)
	}
}

func TestNewReaderAutoDetect(t *testing.T) {
	// SAM path.
	r, err := NewReader(strings.NewReader(sampleSAM))
	if err != nil {
		t.Fatalf("NewReader(SAM): %v", err)
	}
	if _, ok := r.(*SAMReader); !ok {
		t.Errorf("expected SAMReader, got %T", r)
	}

	// BAM path: write a minimal BAM and try to detect it.
	raw := buildHexBAM(t)
	r2, err := NewReader(bytes.NewReader(raw))
	if err != nil {
		t.Fatalf("NewReader(BAM): %v", err)
	}
	if _, ok := r2.(*BAMReader); !ok {
		t.Errorf("expected BAMReader, got %T", r2)
	}
}
