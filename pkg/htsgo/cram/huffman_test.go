package cram

import "testing"

// TestHuffmanDegenerate checks a single-symbol alphabet decodes to its
// symbol without consuming any bits.
func TestHuffmanDegenerate(t *testing.T) {
	tbl, err := newHuffmanTable([]int32{42}, []int32{0})
	if err != nil {
		t.Fatalf("newHuffmanTable: %v", err)
	}
	if !tbl.degenerate {
		t.Fatalf("single-symbol table should be degenerate")
	}
	br := newBitReader(nil)
	for i := 0; i < 5; i++ {
		sym, err := tbl.decode(br)
		if err != nil {
			t.Fatalf("decode #%d: %v", i, err)
		}
		if sym != 42 {
			t.Errorf("decode #%d: got %d want 42", i, sym)
		}
	}
}

// TestHuffmanCanonical builds a multi-symbol canonical table and decodes
// a known bitstream. With symbols A,B,C,D and bit lengths 1,2,3,3 the
// canonical codes are A=0, B=10, C=110, D=111.
func TestHuffmanCanonical(t *testing.T) {
	tbl, err := newHuffmanTable([]int32{'A', 'B', 'C', 'D'}, []int32{1, 2, 3, 3})
	if err != nil {
		t.Fatalf("newHuffmanTable: %v", err)
	}
	var w bitWriter
	// Encode the sequence D C A B A using the canonical codes.
	w.writeBits(0b111, 3) // D
	w.writeBits(0b110, 3) // C
	w.writeBits(0b0, 1)   // A
	w.writeBits(0b10, 2)  // B
	w.writeBits(0b0, 1)   // A
	br := newBitReader(w.bytes())
	want := []int32{'D', 'C', 'A', 'B', 'A'}
	for i, exp := range want {
		got, err := tbl.decode(br)
		if err != nil {
			t.Fatalf("decode #%d: %v", i, err)
		}
		if got != exp {
			t.Errorf("decode #%d: got %d want %d", i, got, exp)
		}
	}
}

// TestHuffmanMismatch rejects an alphabet/bit-length size mismatch and
// an empty alphabet.
func TestHuffmanMismatch(t *testing.T) {
	if _, err := newHuffmanTable([]int32{1, 2}, []int32{1}); err == nil {
		t.Errorf("size mismatch should error")
	}
	if _, err := newHuffmanTable(nil, nil); err == nil {
		t.Errorf("empty alphabet should error")
	}
}

// TestHuffmanZeroLengthMultiSymbol rejects a zero bit length in a
// multi-symbol alphabet.
func TestHuffmanZeroLengthMultiSymbol(t *testing.T) {
	if _, err := newHuffmanTable([]int32{1, 2}, []int32{0, 1}); err == nil {
		t.Errorf("zero bit length in a multi-symbol alphabet should error")
	}
}

// TestHuffmanBadBitLength rejects an out-of-range bit length.
func TestHuffmanBadBitLength(t *testing.T) {
	if _, err := newHuffmanTable([]int32{1, 2}, []int32{1, 99}); err == nil {
		t.Errorf("bit length 99 should be rejected")
	}
}

// TestHuffmanDecodeTruncated checks decode errors on an exhausted
// bitstream rather than looping or panicking.
func TestHuffmanDecodeTruncated(t *testing.T) {
	// Codes A=0, B=10, C=110, D=1110 — every code word ends in a 0, so
	// an all-ones stream never completes a code and runs out of bits.
	tbl, err := newHuffmanTable([]int32{'A', 'B', 'C', 'D'}, []int32{1, 2, 3, 4})
	if err != nil {
		t.Fatalf("newHuffmanTable: %v", err)
	}
	br := newBitReader([]byte{0xff})
	if _, err := tbl.decode(br); err == nil {
		t.Errorf("decode of an incomplete code should error")
	}
}
