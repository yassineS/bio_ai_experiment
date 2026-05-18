package bcf

import (
	"bytes"
	"encoding/binary"
	"math"
	"testing"
)

// buildSharedPortion assembles one BCF "shared" portion for testing. It
// covers the fixed prefix, the typed ID/REF/ALT/FILTER fields, and a few
// INFO entries.
type sharedSpec struct {
	chromID     int32
	pos         int32
	rlen        int32
	qual        float32 // pass NaN with bits MissingFloat32 to indicate missing
	qualMissing bool
	id          string
	alleles     []string // index 0 is REF
	filterIdx   []int32  // dictionary indices into INFO/FILTER dict
	infoKeyIdx  []int32  // dictionary indices into INFO dict
	infoValues  [][]byte // already-encoded typed values
}

func buildShared(spec sharedSpec) []byte {
	var b bytes.Buffer
	binary.Write(&b, binary.LittleEndian, uint32(spec.chromID))
	binary.Write(&b, binary.LittleEndian, uint32(spec.pos))
	binary.Write(&b, binary.LittleEndian, uint32(spec.rlen))

	qualBits := uint32(MissingFloat32)
	if !spec.qualMissing {
		qualBits = math.Float32bits(spec.qual)
	}
	binary.Write(&b, binary.LittleEndian, qualBits)

	binary.Write(&b, binary.LittleEndian, uint16(len(spec.infoKeyIdx)))
	binary.Write(&b, binary.LittleEndian, uint16(len(spec.alleles)))

	packed := uint32(0)
	packed |= 0 & 0x00FFFFFF // n_sample = 0
	packed |= 0 << 24        // n_fmt    = 0
	binary.Write(&b, binary.LittleEndian, packed)

	// id
	if spec.id == "" {
		b.Write(EncodeMissing())
	} else {
		b.Write(EncodeTypedString(spec.id))
	}
	// REF + ALTs
	for _, al := range spec.alleles {
		b.Write(EncodeTypedString(al))
	}
	// FILTER vector
	switch {
	case len(spec.filterIdx) == 0:
		b.Write(EncodeMissing())
	case len(spec.filterIdx) == 1:
		idx := spec.filterIdx[0]
		switch {
		case idx >= -128 && idx <= 127:
			b.Write(EncodeTypedInt8(int8(idx)))
		case idx >= -32768 && idx <= 32767:
			b.Write(EncodeTypedInt16(int16(idx)))
		default:
			b.Write(EncodeTypedInt32(idx))
		}
	default:
		b.Write(EncodeTypedInt32Vec(spec.filterIdx))
	}
	// INFO key/value pairs
	for i, key := range spec.infoKeyIdx {
		b.Write(EncodeTypedInt8(int8(key)))
		b.Write(spec.infoValues[i])
	}
	return b.Bytes()
}

func TestDecodeSharedFull(t *testing.T) {
	hdr := mustHeader(t)

	infoDP := EncodeTypedInt8(42)         // DP=42
	infoAF := EncodeTypedFloat(0.25)      // AF=0.25
	infoTag := EncodeTypedString("hello") // TAG=hello
	infoH2 := EncodeMissing()             // H2 (flag); value is "missing"

	shared := buildShared(sharedSpec{
		chromID:    0,  // chr1
		pos:        99, // -> VCF POS 100
		rlen:       1,
		qual:       29.5,
		id:         "rs42",
		alleles:    []string{"A", "T"},
		filterIdx:  []int32{0}, // PASS
		infoKeyIdx: []int32{2, 3, 4, 5},
		infoValues: [][]byte{infoDP, infoAF, infoTag, infoH2},
	})

	rec := &Record{}
	if err := decodeShared(shared, hdr, rec); err != nil {
		t.Fatalf("decodeShared: %v", err)
	}
	if rec.ChromID != 0 || rec.Pos != 99 || rec.Rlen != 1 {
		t.Fatalf("prefix: %+v", rec)
	}
	if rec.Qual != 29.5 {
		t.Fatalf("qual: %v", rec.Qual)
	}
	if rec.ID != "rs42" {
		t.Fatalf("id: %q", rec.ID)
	}
	if len(rec.Alleles) != 2 || rec.Alleles[0] != "A" || rec.Alleles[1] != "T" {
		t.Fatalf("alleles: %+v", rec.Alleles)
	}
	if len(rec.Filters) != 1 || rec.Filters[0] != 0 {
		t.Fatalf("filters: %+v", rec.Filters)
	}
	if len(rec.InfoKeys) != 4 {
		t.Fatalf("info count: %d", len(rec.InfoKeys))
	}

	v := rec.ToVariant(hdr)
	if v.Chrom != "chr1" || v.Pos != 100 {
		t.Fatalf("variant chrom/pos: %+v", v)
	}
	if v.Ref != "A" || v.Alt[0] != "T" {
		t.Fatalf("variant alleles: %+v", v)
	}
	if v.Filter[0] != "PASS" {
		t.Fatalf("variant filter: %+v", v)
	}
	if v.Info["DP"] != "42" {
		t.Fatalf("info DP: %q", v.Info["DP"])
	}
	if v.Info["AF"] != "0.25" {
		t.Fatalf("info AF: %q", v.Info["AF"])
	}
	if v.Info["TAG"] != "hello" {
		t.Fatalf("info TAG: %q", v.Info["TAG"])
	}
}

func TestDecodeSharedMissingQual(t *testing.T) {
	hdr := mustHeader(t)
	shared := buildShared(sharedSpec{
		chromID:     0,
		pos:         0,
		rlen:        1,
		qualMissing: true,
		alleles:     []string{"G", "C"},
	})
	rec := &Record{}
	if err := decodeShared(shared, hdr, rec); err != nil {
		t.Fatal(err)
	}
	v := rec.ToVariant(hdr)
	if v.Qual != -1 {
		t.Fatalf("missing qual should map to -1, got %v", v.Qual)
	}
	if r := rec.QualString(); r != "." {
		t.Fatalf("QualString: %q", r)
	}
}

func TestDecodeSharedTruncated(t *testing.T) {
	hdr := mustHeader(t)
	rec := &Record{}
	if err := decodeShared([]byte{0, 0, 0}, hdr, rec); err == nil {
		t.Fatal("expected error for truncated shared")
	}
}

func TestDecodeSharedIntArrayInfo(t *testing.T) {
	hdr := mustHeader(t)
	infoDP := EncodeTypedInt32Vec([]int32{1, 2, 3})
	shared := buildShared(sharedSpec{
		chromID:    1,
		pos:        499,
		rlen:       2,
		qual:       0,
		alleles:    []string{"AC", "GT"},
		filterIdx:  []int32{1}, // q10
		infoKeyIdx: []int32{2},
		infoValues: [][]byte{infoDP},
	})
	rec := &Record{}
	if err := decodeShared(shared, hdr, rec); err != nil {
		t.Fatal(err)
	}
	v := rec.ToVariant(hdr)
	if v.Chrom != "chr2" || v.Pos != 500 {
		t.Fatalf("chrom/pos: %+v", v)
	}
	if v.Filter[0] != "q10" {
		t.Fatalf("filter: %+v", v.Filter)
	}
	if v.Info["DP"] != "1,2,3" {
		t.Fatalf("info DP: %q", v.Info["DP"])
	}
}

func TestQualStringWhole(t *testing.T) {
	r := &Record{Qual: 30}
	if s := r.QualString(); s != "30" {
		t.Fatalf("got %q", s)
	}
}

func TestQualOrMissing(t *testing.T) {
	r := &Record{Qual: math.Float32frombits(MissingFloat32)}
	q := r.QualOrMissing()
	if !math.IsNaN(float64(q)) {
		t.Fatalf("expected NaN, got %v", q)
	}
	r.Qual = 1.5
	if r.QualOrMissing() != 1.5 {
		t.Fatalf("expected 1.5")
	}
}

func TestIDOrMissing(t *testing.T) {
	r := &Record{}
	if r.IDOrMissing() != "." {
		t.Fatal("empty ID should become .")
	}
	r.ID = "rs1"
	if r.IDOrMissing() != "rs1" {
		t.Fatal("non-empty ID should pass through")
	}
}

func TestDecodeIndivGT(t *testing.T) {
	hdr := mustHeader(t)
	// FMT key 0 = GT, dict says it's String/Number=1 but on the wire GT
	// uses int8 with the (allele<<1)|phased encoding.
	// Two samples: 0/0 and 0|1.
	// 0 => (-1<<1)|0 = ... actually htslib uses: 0 = missing, (allele+1)<<1 | phase.
	// 0/0 => [2,2] (0+1=1, 1<<1=2, phase=0)
	// 0|1 => [2,5] (allele 0 unphased prefix, allele 1 phased: (1+1)<<1|1 = 5)
	gtVals := []int8{2, 2, 2, 5}
	rawGT := make([]byte, 0, 5)
	// Descriptor for FORMAT/GT: type=int8 (low nibble 1), size=2 (high
	// nibble) which is the per-sample diploid dimension. The wire
	// payload spans nSample (=2) × per-sample-size (=2) = 4 int8 bytes.
	rawGT = append(rawGT, 0x21)
	for _, v := range gtVals {
		rawGT = append(rawGT, byte(v))
	}
	// FMT key 1 = DP, int8 vector with per-sample dim = 1 → descriptor
	// 0x11, payload = nSample (=2) × 1 = 2 bytes.
	rawDP := []byte{0x11, 10, 20}

	// The text-header fixture in buildBCFStream is:
	//   FILTER q10 -> IDX 1, INFO {DP,AF,TAG,H2} -> IDX {2,3,4,5},
	//   FORMAT GT -> IDX 6, FORMAT DP -> IDX 7
	// (PASS is implicit IDX=0). The on-wire FORMAT references must use
	// the unified IDX values, not local FmtTags positions.
	var b bytes.Buffer
	b.Write(EncodeTypedInt8(6)) // FMT key = GT (unified IDX 6)
	b.Write(rawGT)
	b.Write(EncodeTypedInt8(7)) // FMT key = DP (unified IDX 7)
	b.Write(rawDP)

	rec := &Record{NSample: 2, NFmt: 2}
	if err := decodeIndiv(b.Bytes(), hdr, rec); err != nil {
		t.Fatal(err)
	}
	if len(rec.FmtKeys) != 2 {
		t.Fatalf("fmt keys: %+v", rec.FmtKeys)
	}

	// Stage the record through ToVariant by also setting the shared
	// fields the converter expects to read.
	rec.ChromID = 0
	rec.Pos = 0
	rec.Alleles = []string{"A", "G"}
	v := rec.ToVariant(hdr)
	if len(v.Samples) != 2 {
		t.Fatalf("samples: %+v", v.Samples)
	}
	if got := v.Samples[0].Data["GT"]; got != "0/0" {
		t.Fatalf("S1 GT: %q", got)
	}
	if got := v.Samples[1].Data["GT"]; got != "0|1" {
		t.Fatalf("S2 GT: %q", got)
	}
	if v.Samples[0].Data["DP"] != "10" || v.Samples[1].Data["DP"] != "20" {
		t.Fatalf("DP: %+v", v.Samples)
	}
}

func TestFormatTypedFlag(t *testing.T) {
	tv := TypedValue{Descriptor: TypeMissing}
	if got := formatTyped(tv, &DictEntry{Type: "Flag"}); got != "" {
		t.Fatalf("flag should format to empty string, got %q", got)
	}
}

func TestFormatTypedFloatMissing(t *testing.T) {
	tv := TypedValue{
		Descriptor: TypeFloat,
		Floats:     []float32{math.Float32frombits(MissingFloat32), 1.5},
	}
	got := formatTyped(tv, &DictEntry{Type: "Float"})
	if got != ".,1.5" {
		t.Fatalf("got %q", got)
	}
}

func TestFormatTypedIntMissing(t *testing.T) {
	tv := TypedValue{Descriptor: TypeInt32, Ints: []int32{MissingInt32, 7}}
	got := formatTyped(tv, &DictEntry{Type: "Integer"})
	if got != ".,7" {
		t.Fatalf("got %q", got)
	}
}

func TestFormatGTAllele(t *testing.T) {
	if formatGTAllele(0) != "." {
		t.Fatal("0 should be .")
	}
	if formatGTAllele(2) != "0" {
		t.Fatal("2 should be allele 0")
	}
	if formatGTAllele(MissingInt32) != "." {
		t.Fatal("missing should be .")
	}
	if formatGTAllele(EndOfVectorInt32) != "" {
		t.Fatal("end-of-vector should be empty string")
	}
}

func mustHeader(t *testing.T) *Header {
	t.Helper()
	stream := buildBCFStream(t, nil)
	hdr, err := ReadHeader(bytes.NewReader(stream))
	if err != nil {
		t.Fatal(err)
	}
	return hdr
}
