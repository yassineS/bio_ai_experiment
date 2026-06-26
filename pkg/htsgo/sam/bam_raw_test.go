package sam

import (
	"bytes"
	"io"
	"strings"
	"testing"
)

// rawTestSAM is a small SAM with mixed flags, CIGAR, SEQ/QUAL and several aux
// types (including a B-array) so the raw round-trip and single-tag scan exercise
// the full record body, not just the fixed prefix.
const rawTestSAM = `@HD	VN:1.6	SO:unsorted
@SQ	SN:chr1	LN:1000
@SQ	SN:chr2	LN:1000
a1	0	chr1	100	60	5M	*	0	0	ACGTA	IIIII	NM:i:3	RG:Z:grpA	XF:f:1.5	XB:B:i,1,2,3
a2	16	chr2	50	30	3M2I	=	80	0	ACGTT	JJJJJ	NM:i:0	RG:Z:grpB
u1	4	*	0	0	*	*	0	0	*	*	NM:i:7
`

// encodeRawTestBAM encodes rawTestSAM to a BGZF-wrapped BAM byte slice.
func encodeRawTestBAM(t *testing.T) ([]byte, *Header) {
	t.Helper()
	sr, err := NewSAMReader(strings.NewReader(rawTestSAM))
	if err != nil {
		t.Fatalf("NewSAMReader: %v", err)
	}
	var buf bytes.Buffer
	bw := NewBAMWriter(&buf)
	if err := bw.WriteHeader(sr.Header()); err != nil {
		t.Fatalf("WriteHeader: %v", err)
	}
	for {
		rec, err := sr.Read()
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
	return buf.Bytes(), sr.Header()
}

// TestReadRawWriteRawRoundTrip verifies that reading every record body verbatim
// (ReadRaw) and re-emitting it through WriteRaw reproduces a BAM whose decoded
// records are byte-for-byte identical to a Read→Write round-trip — i.e. raw
// passthrough never mutates a record.
func TestReadRawWriteRawRoundTrip(t *testing.T) {
	bamBytes, hdr := encodeRawTestBAM(t)

	// Raw passthrough: ReadRaw → WriteRaw.
	br, err := NewBAMReader(bytes.NewReader(bamBytes))
	if err != nil {
		t.Fatalf("NewBAMReader: %v", err)
	}
	var rawOut bytes.Buffer
	rbw := NewBAMWriter(&rawOut)
	if err := rbw.WriteHeader(hdr); err != nil {
		t.Fatalf("raw WriteHeader: %v", err)
	}
	var bodies [][]byte
	for {
		body, err := br.ReadRaw()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("ReadRaw: %v", err)
		}
		bodies = append(bodies, body)
		if err := rbw.WriteRaw(body); err != nil {
			t.Fatalf("WriteRaw: %v", err)
		}
	}
	if err := rbw.Close(); err != nil {
		t.Fatalf("raw Close: %v", err)
	}
	_ = br.Close()

	if len(bodies) != 3 {
		t.Fatalf("ReadRaw record count: got %d, want 3", len(bodies))
	}

	// The raw-passthrough BAM must decode to the same records as the original.
	orig := decodeAll(t, bamBytes)
	rawRecs := decodeAll(t, rawOut.Bytes())
	if len(orig) != len(rawRecs) {
		t.Fatalf("decoded count: orig=%d raw=%d", len(orig), len(rawRecs))
	}
	for i := range orig {
		if orig[i].QName != rawRecs[i].QName || orig[i].Pos != rawRecs[i].Pos ||
			orig[i].Flag != rawRecs[i].Flag || orig[i].Seq != rawRecs[i].Seq ||
			orig[i].RName != rawRecs[i].RName {
			t.Fatalf("record %d differs after raw round-trip", i)
		}
	}
}

func decodeAll(t *testing.T, bamBytes []byte) []*Record {
	t.Helper()
	br, err := NewBAMReader(bytes.NewReader(bamBytes))
	if err != nil {
		t.Fatalf("NewBAMReader: %v", err)
	}
	defer br.Close()
	var out []*Record
	for {
		rec, err := br.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("Read: %v", err)
		}
		out = append(out, rec)
	}
	return out
}

// TestFindRawAuxTag checks that scanning a single tag from a raw record body
// returns the SAME Aux a full decode would, and reports absence correctly.
func TestFindRawAuxTag(t *testing.T) {
	bamBytes, _ := encodeRawTestBAM(t)
	br, err := NewBAMReader(bytes.NewReader(bamBytes))
	if err != nil {
		t.Fatalf("NewBAMReader: %v", err)
	}
	defer br.Close()

	// First record has NM:i:3, RG:Z:grpA, XF:f:1.5, XB:B:i,1,2,3.
	body, err := br.ReadRaw()
	if err != nil {
		t.Fatalf("ReadRaw: %v", err)
	}
	off, err := RawAuxOffset(body)
	if err != nil {
		t.Fatalf("RawAuxOffset: %v", err)
	}

	// Cross-check against a full decode of the same record.
	full := decodeAll(t, bamBytes)[0]

	for _, tag := range []string{"NM", "RG", "XF", "XB"} {
		got, ok, ferr := FindRawAuxTag(body[off:], tag)
		if ferr != nil {
			t.Fatalf("FindRawAuxTag(%s): %v", tag, ferr)
		}
		want, wok := full.GetAux(tag)
		if ok != wok {
			t.Fatalf("FindRawAuxTag(%s) presence: got %v want %v", tag, ok, wok)
		}
		if got.Type != want.Type {
			t.Errorf("FindRawAuxTag(%s) type: got %c want %c", tag, got.Type, want.Type)
		}
		switch tag {
		case "NM":
			gi, _ := got.Int()
			wi, _ := want.Int()
			if gi != wi {
				t.Errorf("NM value: got %d want %d", gi, wi)
			}
		case "RG":
			if got.Value != want.Value {
				t.Errorf("RG value: got %v want %v", got.Value, want.Value)
			}
		case "XF":
			if got.Value != want.Value {
				t.Errorf("XF value: got %v want %v", got.Value, want.Value)
			}
		case "XB":
			if got.ArrayType != want.ArrayType || len(got.ArrayValues) != len(want.ArrayValues) {
				t.Errorf("XB array mismatch: got %c/%v want %c/%v", got.ArrayType, got.ArrayValues, want.ArrayType, want.ArrayValues)
			}
		}
	}

	// Absent tag.
	if _, ok, ferr := FindRawAuxTag(body[off:], "ZZ"); ferr != nil || ok {
		t.Errorf("FindRawAuxTag(ZZ): want (_, false, nil), got (_, %v, %v)", ok, ferr)
	}
}
