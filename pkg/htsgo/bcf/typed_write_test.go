package bcf

import (
	"bytes"
	"testing"
)

// TestEncodeIntsWidthSelection asserts that encodeInts picks the smallest
// width that fits every value. The on-wire descriptor's low nibble tells us
// which TypeInt* was selected.
func TestEncodeIntsWidthSelection(t *testing.T) {
	cases := []struct {
		vs       []int32
		wantType byte
	}{
		{[]int32{1, 2, 3}, TypeInt8},
		{[]int32{200}, TypeInt16},
		{[]int32{100000}, TypeInt32},
	}
	for _, c := range cases {
		got := encodeInts(c.vs)
		if len(got) == 0 {
			t.Fatalf("empty encoding for %v", c.vs)
		}
		gotType := got[0] & 0x0F
		if gotType != c.wantType {
			t.Errorf("encodeInts(%v) descriptor = %d, want %d", c.vs, gotType, c.wantType)
		}
	}
}

// TestEncodeIntsEmpty returns the missing sentinel for an empty slice.
func TestEncodeIntsEmpty(t *testing.T) {
	if got := encodeInts(nil); !bytes.Equal(got, EncodeMissing()) {
		t.Errorf("empty encodeInts: %v", got)
	}
}

// TestEncodeTypedValueRoundTrip serialises each TypeInt* / TypeFloat / TypeChar
// flavour, parses it back via DecodeTyped, and asserts the result matches.
func TestEncodeTypedValueRoundTrip(t *testing.T) {
	t.Run("int8", func(t *testing.T) {
		tv := TypedValue{Descriptor: TypeInt8, Length: 2, Ints: []int32{1, MissingInt32}}
		buf := encodeTypedValue(tv)
		off := 0
		got, err := DecodeTyped(buf, &off)
		if err != nil {
			t.Fatal(err)
		}
		if got.Descriptor != TypeInt8 || len(got.Ints) != 2 || got.Ints[0] != 1 || got.Ints[1] != MissingInt32 {
			t.Errorf("round-trip failed: %+v", got)
		}
	})
	t.Run("int16", func(t *testing.T) {
		tv := TypedValue{Descriptor: TypeInt16, Length: 1, Ints: []int32{200}}
		buf := encodeTypedValue(tv)
		off := 0
		got, err := DecodeTyped(buf, &off)
		if err != nil {
			t.Fatal(err)
		}
		if got.Descriptor != TypeInt16 || got.Ints[0] != 200 {
			t.Errorf("int16 round-trip failed: %+v", got)
		}
	})
	t.Run("int32", func(t *testing.T) {
		tv := TypedValue{Descriptor: TypeInt32, Length: 1, Ints: []int32{200000}}
		buf := encodeTypedValue(tv)
		off := 0
		got, err := DecodeTyped(buf, &off)
		if err != nil {
			t.Fatal(err)
		}
		if got.Descriptor != TypeInt32 || got.Ints[0] != 200000 {
			t.Errorf("int32 round-trip failed: %+v", got)
		}
	})
	t.Run("float", func(t *testing.T) {
		tv := TypedValue{Descriptor: TypeFloat, Length: 1, Floats: []float32{1.5}}
		buf := encodeTypedValue(tv)
		off := 0
		got, err := DecodeTyped(buf, &off)
		if err != nil {
			t.Fatal(err)
		}
		if got.Descriptor != TypeFloat || got.Floats[0] != 1.5 {
			t.Errorf("float round-trip failed: %+v", got)
		}
	})
	t.Run("char", func(t *testing.T) {
		tv := TypedValue{Descriptor: TypeChar, Length: 3, String: "ABC"}
		buf := encodeTypedValue(tv)
		off := 0
		got, err := DecodeTyped(buf, &off)
		if err != nil {
			t.Fatal(err)
		}
		if got.Descriptor != TypeChar || got.String != "ABC" {
			t.Errorf("char round-trip failed: %+v", got)
		}
	})
	t.Run("missing", func(t *testing.T) {
		tv := TypedValue{Descriptor: TypeMissing, Length: 0}
		buf := encodeTypedValue(tv)
		off := 0
		got, err := DecodeTyped(buf, &off)
		if err != nil {
			t.Fatal(err)
		}
		if got.Descriptor != TypeMissing {
			t.Errorf("missing round-trip failed: %+v", got)
		}
	})
}

// TestIntWidthSentinels checks that the missing/end-of-vector sentinels keep
// their int8 width — the smallest viable representation.
func TestIntWidthSentinels(t *testing.T) {
	if got := intWidth(MissingInt32); got != 1 {
		t.Errorf("missing width: %d", got)
	}
	if got := intWidth(EndOfVectorInt32); got != 1 {
		t.Errorf("end-of-vec width: %d", got)
	}
}
