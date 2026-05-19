package bcf

import (
	"bytes"
	"encoding/binary"
	"io"
	"testing"
)

// recordSpec captures the inputs to buildRecord. It mirrors buildShared but
// also lets the caller wire in a per-sample portion.
type recordSpec struct {
	chromID    int32
	pos        int32
	rlen       int32
	qual       float32
	qualMiss   bool
	id         string
	alleles    []string
	filters    []int32
	infoKeys   []int32
	infoVals   [][]byte
	nSample    uint32
	nFmt       uint8
	indivBytes []byte
}

// buildRecord encodes one full record (length prefixes + shared + indiv).
func buildRecord(spec recordSpec) []byte {
	// Custom buildShared that lets us control n_sample / n_fmt.
	var shared bytes.Buffer
	binary.Write(&shared, binary.LittleEndian, uint32(spec.chromID))
	binary.Write(&shared, binary.LittleEndian, uint32(spec.pos))
	binary.Write(&shared, binary.LittleEndian, uint32(spec.rlen))
	if spec.qualMiss {
		binary.Write(&shared, binary.LittleEndian, MissingFloat32)
	} else {
		var b [4]byte
		binary.LittleEndian.PutUint32(b[:], floatBits(spec.qual))
		shared.Write(b[:])
	}
	binary.Write(&shared, binary.LittleEndian, uint16(len(spec.infoKeys)))
	binary.Write(&shared, binary.LittleEndian, uint16(len(spec.alleles)))
	packed := (spec.nSample & 0x00FFFFFF) | (uint32(spec.nFmt) << 24)
	binary.Write(&shared, binary.LittleEndian, packed)

	if spec.id == "" {
		shared.Write(EncodeMissing())
	} else {
		shared.Write(EncodeTypedString(spec.id))
	}
	for _, al := range spec.alleles {
		shared.Write(EncodeTypedString(al))
	}
	if len(spec.filters) == 0 {
		shared.Write(EncodeMissing())
	} else if len(spec.filters) == 1 {
		idx := spec.filters[0]
		shared.Write(EncodeTypedInt8(int8(idx)))
	} else {
		shared.Write(EncodeTypedInt32Vec(spec.filters))
	}
	for i, k := range spec.infoKeys {
		shared.Write(EncodeTypedInt8(int8(k)))
		shared.Write(spec.infoVals[i])
	}

	var rec bytes.Buffer
	binary.Write(&rec, binary.LittleEndian, uint32(shared.Len()))
	binary.Write(&rec, binary.LittleEndian, uint32(len(spec.indivBytes)))
	rec.Write(shared.Bytes())
	rec.Write(spec.indivBytes)
	return rec.Bytes()
}

func TestReaderRoundTrip(t *testing.T) {
	// Two records, no per-sample data.
	rec1 := buildRecord(recordSpec{
		chromID:  0,
		pos:      99,
		rlen:     1,
		qual:     30,
		id:       "rs1",
		alleles:  []string{"A", "T"},
		filters:  []int32{0},
		infoKeys: []int32{2},
		infoVals: [][]byte{EncodeTypedInt8(50)},
	})
	rec2 := buildRecord(recordSpec{
		chromID:  1,
		pos:      199,
		rlen:     1,
		qualMiss: true,
		alleles:  []string{"G", "C"},
		filters:  []int32{1}, // q10
	})

	body := append(rec1, rec2...)
	stream := buildBCFStream(t, body)

	r, err := NewReader(bytes.NewReader(stream))
	if err != nil {
		t.Fatalf("NewReader: %v", err)
	}
	if r.Header() == nil {
		t.Fatal("nil header from Reader")
	}
	recs, err := r.ReadAll()
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if len(recs) != 2 {
		t.Fatalf("got %d records, want 2", len(recs))
	}
	if recs[0].ID != "rs1" || recs[1].ID != "" {
		t.Fatalf("ids: %q %q", recs[0].ID, recs[1].ID)
	}
	if recs[0].Pos+1 != 100 || recs[1].Pos+1 != 200 {
		t.Fatalf("positions: %d %d", recs[0].Pos, recs[1].Pos)
	}

	v0 := recs[0].ToVariant(r.Header())
	if v0.Chrom != "chr1" || v0.Pos != 100 || v0.Info["DP"] != "50" {
		t.Fatalf("variant[0]: %+v", v0)
	}
	v1 := recs[1].ToVariant(r.Header())
	if v1.Qual != -1 {
		t.Fatalf("variant[1] qual should be missing, got %v", v1.Qual)
	}
}

func TestReaderEmptyStream(t *testing.T) {
	stream := buildBCFStream(t, nil)
	r, err := NewReader(bytes.NewReader(stream))
	if err != nil {
		t.Fatal(err)
	}
	if _, err := r.Read(); err != io.EOF {
		t.Fatalf("expected EOF, got %v", err)
	}
}

func TestReaderTruncatedLPrefix(t *testing.T) {
	stream := buildBCFStream(t, []byte{0x01, 0x00}) // 2 bytes, less than 8
	r, err := NewReader(bytes.NewReader(stream))
	if err != nil {
		t.Fatal(err)
	}
	if _, err := r.Read(); err == nil {
		t.Fatal("expected truncated error")
	}
}

func TestNewReaderWithHeader(t *testing.T) {
	stream := buildBCFStream(t, nil)
	hdr, err := ReadHeader(bytes.NewReader(stream))
	if err != nil {
		t.Fatal(err)
	}
	r := NewReaderWithHeader(hdr)
	if r.Header() != hdr {
		t.Fatal("Header() should return the same instance")
	}
}
