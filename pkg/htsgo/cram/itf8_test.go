package cram

import (
	"bytes"
	"io"
	"testing"
)

// TestReadITF8 exercises the ITF-8 reader across the full range of
// 1-to-5-byte encodings, including the boundary values where the byte
// count steps up and the negative-valued (high-bit-set) case.
func TestReadITF8(t *testing.T) {
	cases := []struct {
		name  string
		bytes []byte
		want  int32
		n     int
	}{
		{"zero", []byte{0x00}, 0, 1},
		{"one", []byte{0x01}, 1, 1},
		{"max-1-byte", []byte{0x7f}, 127, 1},
		{"min-2-byte", []byte{0x80, 0x80}, 128, 2},
		{"mid-2-byte", []byte{0xa0, 0x00}, 0x2000, 2},
		{"min-3-byte", []byte{0xc0, 0x40, 0x00}, 0x4000, 3},
		{"min-4-byte", []byte{0xe0, 0x20, 0x00, 0x00}, 0x200000, 4},
		{"five-byte-max", []byte{0xff, 0xff, 0xff, 0xff, 0x0f}, -1, 5},
		{"five-byte-small", []byte{0xf0, 0x00, 0x00, 0x01, 0x00}, 0x10, 5},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			v, n, err := readITF8(bytes.NewReader(c.bytes))
			if err != nil {
				t.Fatalf("readITF8: %v", err)
			}
			if v != c.want {
				t.Errorf("value = %d (%#x), want %d (%#x)", v, uint32(v), c.want, uint32(c.want))
			}
			if n != c.n {
				t.Errorf("byte count = %d, want %d", n, c.n)
			}
		})
	}
}

// TestReadITF8Truncated checks that a truncated ITF-8 value reports an
// unexpected-EOF error rather than a clean EOF or a panic.
func TestReadITF8Truncated(t *testing.T) {
	for _, in := range [][]byte{nil, {0x80}, {0xc0, 0x00}, {0xff, 0x00, 0x00}} {
		if _, _, err := readITF8(bytes.NewReader(in)); err == nil {
			t.Errorf("expected error for truncated ITF-8 %v", in)
		}
	}
}

// TestReadLTF8 exercises the LTF-8 reader across the 1-to-9-byte
// encodings, including the all-ones first byte that signals a full
// 8-byte trailing value.
func TestReadLTF8(t *testing.T) {
	cases := []struct {
		name  string
		bytes []byte
		want  int64
		n     int
	}{
		{"zero", []byte{0x00}, 0, 1},
		{"max-1-byte", []byte{0x7f}, 127, 1},
		{"min-2-byte", []byte{0x80, 0x80}, 128, 2},
		{"min-3-byte", []byte{0xc0, 0x40, 0x00}, 0x4000, 3},
		{"eight-byte-zero", []byte{0xfe, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00}, 0, 8},
		{"nine-byte-one", []byte{0xff, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x01}, 1, 9},
		{"nine-byte-max", []byte{0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff}, -1, 9},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			v, n, err := readLTF8(bytes.NewReader(c.bytes))
			if err != nil {
				t.Fatalf("readLTF8: %v", err)
			}
			if v != c.want {
				t.Errorf("value = %d, want %d", v, c.want)
			}
			if n != c.n {
				t.Errorf("byte count = %d, want %d", n, c.n)
			}
		})
	}
}

// TestReadLTF8Truncated checks that a truncated LTF-8 value reports an
// error rather than a clean EOF or a panic.
func TestReadLTF8Truncated(t *testing.T) {
	for _, in := range [][]byte{nil, {0x80}, {0xff, 0x00, 0x00}} {
		if _, _, err := readLTF8(bytes.NewReader(in)); err == nil {
			t.Errorf("expected error for truncated LTF-8 %v", in)
		}
	}
}

// TestFileDefinitionRoundTrip checks the file-definition fields parse
// out of a hand-built 26-byte header.
func TestFileDefinitionRoundTrip(t *testing.T) {
	buf := make([]byte, fileDefSize)
	copy(buf, "CRAM")
	buf[4] = 3
	buf[5] = 1
	copy(buf[6:], "myfile.cram")
	def, err := readFileDefinition(bytes.NewReader(buf))
	if err != nil {
		t.Fatalf("readFileDefinition: %v", err)
	}
	if def.Major != 3 || def.Minor != 1 {
		t.Errorf("version = %s, want 3.1", def.VersionString())
	}
	if def.FileIDString() != "myfile.cram" {
		t.Errorf("file id = %q, want %q", def.FileIDString(), "myfile.cram")
	}
	if !def.hasCRC() {
		t.Errorf("v3 should report hasCRC")
	}
}

// TestFileDefinitionV2NoCRC checks a v2 file definition reports no CRC.
func TestFileDefinitionV2NoCRC(t *testing.T) {
	buf := make([]byte, fileDefSize)
	copy(buf, "CRAM")
	buf[4] = 2
	def, err := readFileDefinition(bytes.NewReader(buf))
	if err != nil {
		t.Fatalf("readFileDefinition: %v", err)
	}
	if def.hasCRC() {
		t.Errorf("v2 should not report hasCRC")
	}
}

// TestVersionString pins the VersionString formatting.
func TestVersionString(t *testing.T) {
	if got := (FileDefinition{Major: 3, Minor: 0}).VersionString(); got != "3.0" {
		t.Errorf("VersionString = %q, want 3.0", got)
	}
	if got := (FileDefinition{Major: 2, Minor: 1}).VersionString(); got != "2.1" {
		t.Errorf("VersionString = %q, want 2.1", got)
	}
}

// TestFileIDStringEmpty checks an all-NUL file id stringifies to empty.
func TestFileIDStringEmpty(t *testing.T) {
	if got := (FileDefinition{}).FileIDString(); got != "" {
		t.Errorf("empty file id = %q, want \"\"", got)
	}
}

// TestEOFToUnexpected checks the EOF-translation helper turns a clean
// io.EOF mid-structure into an unexpected-EOF-flavoured error.
func TestEOFToUnexpected(t *testing.T) {
	err := eofToUnexpected(io.EOF, "thing")
	if err == nil {
		t.Fatal("expected non-nil error")
	}
	if got := err.Error(); !bytes.Contains([]byte(got), []byte("truncated")) {
		t.Errorf("error %q should mention truncation", got)
	}
}

// TestAppendITF8RoundTrip checks the ITF-8 encoder against the decoder
// at the size-class boundaries (1- through 5-byte encodings) and the
// int32 extremes. itf8At is the same decoder the reader relies on.
func TestAppendITF8RoundTrip(t *testing.T) {
	values := []int32{
		0, 1, 127, 128, 16383, 16384, 2097151, 2097152,
		268435455, 268435456, 1 << 30, 2147483647, -1, -128, -2147483648,
	}
	for _, v := range values {
		enc := appendITF8(nil, v)
		got, n, err := itf8At(enc, 0)
		if err != nil {
			t.Errorf("itf8At(%d-encoding): %v", v, err)
			continue
		}
		if got != v {
			t.Errorf("ITF-8 round-trip: encoded %d, decoded %d", v, got)
		}
		if n != len(enc) {
			t.Errorf("ITF-8 %d: decoder consumed %d bytes, encoding is %d", v, n, len(enc))
		}
	}
}

// TestAppendLTF8RoundTrip checks the LTF-8 encoder against the decoder
// at the size-class boundaries and the int64 extremes.
func TestAppendLTF8RoundTrip(t *testing.T) {
	values := []int64{
		0, 1, 127, 128, 1<<35 - 1, 1 << 35, 1<<56 - 1, 1 << 56,
		9223372036854775807, -1, -9223372036854775808,
	}
	for _, v := range values {
		enc := appendLTF8(nil, v)
		got, n, err := ltf8At(enc, 0)
		if err != nil {
			t.Errorf("ltf8At(%d-encoding): %v", v, err)
			continue
		}
		if got != v {
			t.Errorf("LTF-8 round-trip: encoded %d, decoded %d", v, got)
		}
		if n != len(enc) {
			t.Errorf("LTF-8 %d: decoder consumed %d bytes, encoding is %d", v, n, len(enc))
		}
	}
}
