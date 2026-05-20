package cram

import (
	"errors"
	"io"
	"testing"
)

// bitWriter is a test helper that builds an MSB-first bitstream the
// bitReader can read back.
type bitWriter struct {
	out []byte
	cur byte
	nb  uint
}

// writeBit appends a single bit.
func (w *bitWriter) writeBit(b uint32) {
	w.cur = (w.cur << 1) | byte(b&1)
	w.nb++
	if w.nb == 8 {
		w.out = append(w.out, w.cur)
		w.cur, w.nb = 0, 0
	}
}

// writeBits appends the low n bits of v, most-significant first.
func (w *bitWriter) writeBits(v uint32, n uint) {
	for i := int(n) - 1; i >= 0; i-- {
		w.writeBit((v >> uint(i)) & 1)
	}
}

// writeUnary appends n 1-bits then a 0-bit.
func (w *bitWriter) writeUnary(n int) {
	for i := 0; i < n; i++ {
		w.writeBit(1)
	}
	w.writeBit(0)
}

// bytes flushes any partial byte (zero-padded) and returns the stream.
func (w *bitWriter) bytes() []byte {
	if w.nb > 0 {
		w.out = append(w.out, w.cur<<(8-w.nb))
		w.cur, w.nb = 0, 0
	}
	return w.out
}

// TestBitReaderReadBits round-trips fixed-width values through the bit
// reader.
func TestBitReaderReadBits(t *testing.T) {
	var w bitWriter
	want := []struct {
		v uint32
		n uint
	}{
		{0b101, 3}, {0b1111_0000, 8}, {1, 1}, {0, 5}, {0x1234, 16},
	}
	for _, c := range want {
		w.writeBits(c.v, c.n)
	}
	br := newBitReader(w.bytes())
	for i, c := range want {
		got, err := br.readBits(c.n)
		if err != nil {
			t.Fatalf("readBits #%d: %v", i, err)
		}
		if got != c.v {
			t.Errorf("readBits #%d: got %d want %d", i, got, c.v)
		}
	}
}

// TestBitReaderZeroBits checks reading zero bits yields zero and
// consumes nothing.
func TestBitReaderZeroBits(t *testing.T) {
	br := newBitReader([]byte{0xff})
	v, err := br.readBits(0)
	if err != nil || v != 0 {
		t.Fatalf("readBits(0) = %d, %v; want 0, nil", v, err)
	}
	if br.exhausted() {
		t.Errorf("zero-bit read should not have consumed the stream")
	}
}

// TestBitReaderTooManyBits checks an over-wide request is rejected.
func TestBitReaderTooManyBits(t *testing.T) {
	br := newBitReader([]byte{0, 0, 0, 0, 0})
	if _, err := br.readBits(33); err == nil {
		t.Errorf("readBits(33) should error")
	}
	br64 := newBitReader(make([]byte, 9))
	if _, err := br64.readBits64(65); err == nil {
		t.Errorf("readBits64(65) should error")
	}
}

// TestBitReaderExhausted checks the reader reports EOF cleanly when out
// of bits.
func TestBitReaderExhausted(t *testing.T) {
	br := newBitReader([]byte{0xab})
	if _, err := br.readBits(8); err != nil {
		t.Fatalf("readBits(8): %v", err)
	}
	if !br.exhausted() {
		t.Errorf("reader should be exhausted")
	}
	if _, err := br.readBit(); !errors.Is(err, io.ErrUnexpectedEOF) {
		t.Errorf("readBit past end: got %v, want ErrUnexpectedEOF", err)
	}
}

// TestBitReaderReadUnary round-trips unary counts.
func TestBitReaderReadUnary(t *testing.T) {
	var w bitWriter
	counts := []int{0, 1, 5, 13, 0}
	for _, c := range counts {
		w.writeUnary(c)
	}
	br := newBitReader(w.bytes())
	for i, want := range counts {
		got, err := br.readUnary(64)
		if err != nil {
			t.Fatalf("readUnary #%d: %v", i, err)
		}
		if got != want {
			t.Errorf("readUnary #%d: got %d want %d", i, got, want)
		}
	}
}

// TestBitReaderUnaryLimit checks a runaway all-ones stream is bounded.
func TestBitReaderUnaryLimit(t *testing.T) {
	allOnes := make([]byte, 16)
	for i := range allOnes {
		allOnes[i] = 0xff
	}
	br := newBitReader(allOnes)
	if _, err := br.readUnary(8); err == nil {
		t.Errorf("readUnary on an all-ones stream should hit the limit")
	}
}

// TestBitReaderAlign checks align discards a partial byte.
func TestBitReaderAlign(t *testing.T) {
	br := newBitReader([]byte{0b1010_0000, 0xcd})
	if _, err := br.readBits(3); err != nil {
		t.Fatalf("readBits(3): %v", err)
	}
	br.align()
	got, err := br.readBits(8)
	if err != nil {
		t.Fatalf("readBits(8) after align: %v", err)
	}
	if got != 0xcd {
		t.Errorf("after align got %#x, want 0xcd", got)
	}
}

// TestBitReaderReadBits64 round-trips a value wider than 32 bits.
func TestBitReaderReadBits64(t *testing.T) {
	var w bitWriter
	const want uint64 = 0x1234_5678_9abc_def0
	for i := 63; i >= 0; i-- {
		w.writeBit(uint32((want >> uint(i)) & 1))
	}
	br := newBitReader(w.bytes())
	got, err := br.readBits64(64)
	if err != nil {
		t.Fatalf("readBits64: %v", err)
	}
	if got != want {
		t.Errorf("readBits64 = %#x, want %#x", got, want)
	}
}
