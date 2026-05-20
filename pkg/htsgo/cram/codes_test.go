package cram

import (
	"math/bits"
	"testing"
)

// encGamma writes the Elias gamma code of v (v >= 1) to w, matching the
// readGamma decoder: a unary length prefix then the value's low bits.
func encGamma(w *bitWriter, v uint32) {
	if v == 0 {
		panic("gamma code is undefined for 0")
	}
	l := bits.Len32(v) - 1 // floor(log2 v)
	w.writeUnary(l)
	w.writeBits(v&((1<<uint(l))-1), uint(l))
}

// TestGammaRoundTrip encodes a range of values and decodes them back.
func TestGammaRoundTrip(t *testing.T) {
	values := []uint32{1, 2, 3, 4, 5, 7, 8, 100, 1000, 65535, 1 << 20}
	var w bitWriter
	for _, v := range values {
		encGamma(&w, v)
	}
	br := newBitReader(w.bytes())
	for i, want := range values {
		// readGamma applies a zero offset here.
		enc := &Encoding{ID: EncodingGamma, Offset: 0}
		got, err := enc.decodeInt(&seriesSource{core: br})
		if err != nil {
			t.Fatalf("gamma #%d: %v", i, err)
		}
		if uint32(got) != want {
			t.Errorf("gamma #%d: got %d want %d", i, got, want)
		}
	}
}

// TestGammaOffset checks the gamma encoding subtracts its offset.
func TestGammaOffset(t *testing.T) {
	var w bitWriter
	encGamma(&w, 10)
	enc := &Encoding{ID: EncodingGamma, Offset: 3}
	got, err := enc.decodeInt(&seriesSource{core: newBitReader(w.bytes())})
	if err != nil {
		t.Fatalf("gamma: %v", err)
	}
	if got != 7 {
		t.Errorf("gamma with offset 3: got %d want 7", got)
	}
}

// encSubexp writes the sub-exponential code of v with order k, matching
// the readSubexp decoder.
func encSubexp(w *bitWriter, v uint32, k int) {
	// Determine bucket: n=0 covers [0, 2^k); n>=1 covers
	// [2^(k+n-1), 2^(k+n)).
	if v < (1 << uint(k)) {
		w.writeUnary(0)
		w.writeBits(v, uint(k))
		return
	}
	n := 1
	for v >= (1 << uint(k+n)) {
		n++
	}
	nbits := k + n - 1
	w.writeUnary(n)
	w.writeBits(v-(1<<uint(nbits)), uint(nbits))
}

// TestSubexpRoundTrip exercises the sub-exponential code across the
// bucket boundaries for several order parameters.
func TestSubexpRoundTrip(t *testing.T) {
	for _, k := range []int{0, 1, 2, 5} {
		values := []uint32{0, 1, 2, 3, 7, 8, 9, 31, 32, 100, 1000, 100000}
		var w bitWriter
		for _, v := range values {
			encSubexp(&w, v, k)
		}
		br := newBitReader(w.bytes())
		for i, want := range values {
			got, err := readSubexp(br, k)
			if err != nil {
				t.Fatalf("subexp k=%d #%d: %v", k, i, err)
			}
			if uint32(got) != want {
				t.Errorf("subexp k=%d #%d: got %d want %d", k, i, got, want)
			}
		}
	}
}

// encGolomb writes the Golomb code of v with divisor m, matching the
// readGolomb decoder (unary quotient + truncated-binary remainder).
func encGolomb(w *bitWriter, v, m int) {
	q := v / m
	r := v % m
	w.writeUnary(q)
	if m == 1 {
		return
	}
	b := 0
	for (1 << uint(b)) < m {
		b++
	}
	cutoff := (1 << uint(b)) - m
	if r < cutoff {
		w.writeBits(uint32(r), uint(b-1))
	} else {
		w.writeBits(uint32(r+cutoff), uint(b))
	}
}

// TestGolombRoundTrip exercises Golomb coding for several divisors,
// including non-power-of-two ones that use truncated binary remainders.
func TestGolombRoundTrip(t *testing.T) {
	for _, m := range []int{1, 2, 3, 5, 7, 8, 13, 64} {
		values := []int{0, 1, 2, 3, 4, 5, 10, 50, 127, 500}
		var w bitWriter
		for _, v := range values {
			encGolomb(&w, v, m)
		}
		br := newBitReader(w.bytes())
		for i, want := range values {
			got, err := readGolomb(br, m)
			if err != nil {
				t.Fatalf("golomb m=%d #%d: %v", m, i, err)
			}
			if int(got) != want {
				t.Errorf("golomb m=%d #%d: got %d want %d", m, i, got, want)
			}
		}
	}
}

// TestGolombRiceRoundTrip exercises the Golomb-Rice path: a Golomb code
// whose divisor is 1<<K.
func TestGolombRiceRoundTrip(t *testing.T) {
	for _, k := range []int{0, 1, 3, 6} {
		m := 1 << uint(k)
		values := []int{0, 1, 7, 8, 100, 255}
		var w bitWriter
		for _, v := range values {
			encGolomb(&w, v, m)
		}
		enc := &Encoding{ID: EncodingGolombRice, K: int32(k), Offset: 0}
		s := &seriesSource{core: newBitReader(w.bytes())}
		for i, want := range values {
			got, err := enc.decodeInt(s)
			if err != nil {
				t.Fatalf("golomb-rice k=%d #%d: %v", k, i, err)
			}
			if int(got) != want {
				t.Errorf("golomb-rice k=%d #%d: got %d want %d", k, i, got, want)
			}
		}
	}
}

// TestGolombInvalidDivisor checks a non-positive divisor errors.
func TestGolombInvalidDivisor(t *testing.T) {
	br := newBitReader([]byte{0})
	if _, err := readGolomb(br, 0); err == nil {
		t.Errorf("readGolomb with m=0 should error")
	}
}

// TestSubexpBadK rejects an out-of-range order parameter.
func TestSubexpBadK(t *testing.T) {
	if _, err := readSubexp(newBitReader([]byte{0}), -1); err == nil {
		t.Errorf("readSubexp with k=-1 should error")
	}
	if _, err := readSubexp(newBitReader([]byte{0}), 99); err == nil {
		t.Errorf("readSubexp with k=99 should error")
	}
}

// TestGammaOverlongPrefix rejects a gamma code whose length prefix
// exceeds 31 bits.
func TestGammaOverlongPrefix(t *testing.T) {
	var w bitWriter
	w.writeUnary(32) // 32 leading ones — over the 31-bit cap
	w.writeBits(0, 1)
	if _, err := readGamma(newBitReader(w.bytes())); err == nil {
		t.Errorf("a 32-bit gamma length should be rejected")
	}
}

// TestCodesTruncated checks the integer codes error — never panic — on
// a bitstream that ends mid-code.
func TestCodesTruncated(t *testing.T) {
	if _, err := readGamma(newBitReader(nil)); err == nil {
		t.Errorf("readGamma on empty stream should error")
	}
	// 0xf0 = 1111 0000: a unary prefix of 4 ones then 0, leaving four
	// bits — short of the k=8 tail this subexp code then needs.
	if _, err := readSubexp(newBitReader([]byte{0xf0}), 8); err == nil {
		t.Errorf("readSubexp on a truncated stream should error")
	}
	if _, err := readGolomb(newBitReader([]byte{0xff}), 5); err == nil {
		t.Errorf("readGolomb on a truncated stream should error")
	}
}
