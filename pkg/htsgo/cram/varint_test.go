package cram

import (
	"bytes"
	"errors"
	"io"
	"testing"
)

// The expected byte sequences below are the big-endian uint7 / zig-zag
// sint7 encodings htscodecs/varint.h (BIG_END branch) produces and
// cram_io.c's uint7_get_* read. They were cross-checked against the
// upstream encoding (e.g. 0x454f46 -> 82 95 9e 46 is the "EOF" marker the
// v4 EOF container stores in its ref_seq_start field).

func TestUint7At32(t *testing.T) {
	cases := []struct {
		val   int32
		bytes []byte
	}{
		{0, []byte{0x00}},
		{1, []byte{0x01}},
		{127, []byte{0x7f}},
		{128, []byte{0x81, 0x00}},
		{300, []byte{0x82, 0x2c}},
		{16384, []byte{0x81, 0x80, 0x00}},
		{0x454f46, []byte{0x82, 0x95, 0x9e, 0x46}},
		{2097152, []byte{0x81, 0x80, 0x80, 0x00}},
		{2147483647, []byte{0x87, 0xff, 0xff, 0xff, 0x7f}}, // 2^31 - 1
		{-1, []byte{0x8f, 0xff, 0xff, 0xff, 0x7f}},         // 0xffffffff as uint32
	}
	for _, c := range cases {
		got, n, err := uint7At32(c.bytes, 0)
		if err != nil {
			t.Errorf("uint7At32(%x): unexpected error %v", c.bytes, err)
			continue
		}
		if got != c.val {
			t.Errorf("uint7At32(%x) = %d, want %d", c.bytes, got, c.val)
		}
		if n != len(c.bytes) {
			t.Errorf("uint7At32(%x) consumed %d bytes, want %d", c.bytes, n, len(c.bytes))
		}
		// Round-trip through the writer for the non-negative values.
		if c.val >= 0 {
			if enc := appendUint7(nil, uint64(uint32(c.val))); !bytes.Equal(enc, c.bytes) {
				t.Errorf("appendUint7(%d) = %x, want %x", c.val, enc, c.bytes)
			}
		}
	}
}

func TestSint7At32(t *testing.T) {
	cases := []struct {
		val   int32
		bytes []byte
	}{
		{0, []byte{0x00}},
		{-1, []byte{0x01}},
		{1, []byte{0x02}},
		{-2, []byte{0x03}},
		{2, []byte{0x04}},
		{63, []byte{0x7e}},
		{-64, []byte{0x7f}},
		{1000, []byte{0x8f, 0x50}},
		{-1000, []byte{0x8f, 0x4f}},
	}
	for _, c := range cases {
		got, n, err := sint7At32(c.bytes, 0)
		if err != nil {
			t.Errorf("sint7At32(%x): unexpected error %v", c.bytes, err)
			continue
		}
		if got != c.val {
			t.Errorf("sint7At32(%x) = %d, want %d", c.bytes, got, c.val)
		}
		if n != len(c.bytes) {
			t.Errorf("sint7At32(%x) consumed %d bytes, want %d", c.bytes, n, len(c.bytes))
		}
		if enc := appendSint7(nil, int64(c.val)); !bytes.Equal(enc, c.bytes) {
			t.Errorf("appendSint7(%d) = %x, want %x", c.val, enc, c.bytes)
		}
	}
}

func TestUint7At64(t *testing.T) {
	cases := []struct {
		val   int64
		bytes []byte
	}{
		{0, []byte{0x00}},
		{0x454f46, []byte{0x82, 0x95, 0x9e, 0x46}},
		{1 << 35, []byte{0x81, 0x80, 0x80, 0x80, 0x80, 0x00}},
		// 2^53 (beyond int32) round-trips through the 64-bit reader.
		{1 << 53, []byte{0x90, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x00}},
	}
	for _, c := range cases {
		got, n, err := uint7At64(c.bytes, 0)
		if err != nil {
			t.Errorf("uint7At64(%x): unexpected error %v", c.bytes, err)
			continue
		}
		if got != c.val {
			t.Errorf("uint7At64(%x) = %d, want %d", c.bytes, got, c.val)
		}
		if n != len(c.bytes) {
			t.Errorf("uint7At64(%x) consumed %d bytes, want %d", c.bytes, n, len(c.bytes))
		}
		if enc := appendUint7(nil, uint64(c.val)); !bytes.Equal(enc, c.bytes) {
			t.Errorf("appendUint7(%d) = %x, want %x", c.val, enc, c.bytes)
		}
	}
}

func TestSint7At64(t *testing.T) {
	cases := []struct {
		val   int64
		bytes []byte
	}{
		{0, []byte{0x00}},
		{-1, []byte{0x01}},
		{1 << 40, []byte{0xc0, 0x80, 0x80, 0x80, 0x80, 0x00}},
		{-(1 << 40), []byte{0xbf, 0xff, 0xff, 0xff, 0xff, 0x7f}},
	}
	for _, c := range cases {
		got, n, err := sint7At64(c.bytes, 0)
		if err != nil {
			t.Errorf("sint7At64(%x): unexpected error %v", c.bytes, err)
			continue
		}
		if got != c.val {
			t.Errorf("sint7At64(%x) = %d, want %d", c.bytes, got, c.val)
		}
		if n != len(c.bytes) {
			t.Errorf("sint7At64(%x) consumed %d bytes, want %d", c.bytes, n, len(c.bytes))
		}
		if enc := appendSint7(nil, c.val); !bytes.Equal(enc, c.bytes) {
			t.Errorf("appendSint7(%d) = %x, want %x", c.val, enc, c.bytes)
		}
	}
}

// TestReadUint7Streaming exercises the io.Reader-based decoders the v4
// container- and block-header parsers use, confirming they agree with the
// byte-slice decoders and report the same byte count.
func TestReadUint7Streaming(t *testing.T) {
	in := []byte{0x82, 0x95, 0x9e, 0x46} // 0x454f46
	v32, n, err := readUint7_32(bytes.NewReader(in))
	if err != nil || v32 != 0x454f46 || n != 4 {
		t.Fatalf("readUint7_32 = %d, %d, %v; want 4542278, 4, nil", v32, n, err)
	}
	v64, n, err := readUint7_64(bytes.NewReader(in))
	if err != nil || v64 != 0x454f46 || n != 4 {
		t.Fatalf("readUint7_64 = %d, %d, %v; want 4542278, 4, nil", v64, n, err)
	}
	s32, n, err := readSint7_32(bytes.NewReader([]byte{0x01}))
	if err != nil || s32 != -1 || n != 1 {
		t.Fatalf("readSint7_32 = %d, %d, %v; want -1, 1, nil", s32, n, err)
	}
}

// TestUint7Truncated confirms a value cut short mid-continuation is
// reported as an unexpected EOF rather than read as a short value.
func TestUint7Truncated(t *testing.T) {
	// 0x81 has the continuation bit set but no following byte.
	if _, _, err := uint7At32([]byte{0x81}, 0); !errors.Is(err, io.ErrUnexpectedEOF) {
		t.Errorf("uint7At32 truncated: err = %v, want io.ErrUnexpectedEOF", err)
	}
	if _, _, err := readUint7_32(bytes.NewReader([]byte{0x81})); !errors.Is(err, io.ErrUnexpectedEOF) {
		t.Errorf("readUint7_32 truncated: err = %v, want io.ErrUnexpectedEOF", err)
	}
	if _, _, err := uint7At32(nil, 0); !errors.Is(err, io.ErrUnexpectedEOF) {
		t.Errorf("uint7At32 empty: err = %v, want io.ErrUnexpectedEOF", err)
	}
}
