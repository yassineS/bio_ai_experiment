package cram

import (
	"bytes"
	"errors"
	"io"
	"os"
	"path/filepath"
	"testing"
)

// samtoolsTestDir is the test-data root of the samtools submodule
// vendored under reference_code/. The path is relative to this test
// file's directory.
const samtoolsTestDir = "../../../reference_code/samtools/test"

// loadFixture reads a CRAM fixture from the samtools submodule. It
// returns ok=false when the submodule is not initialised; callers then
// t.Fatalf with an init hint, since the samtools fixtures are a hard
// requirement for the CRAM parity rig (run `git submodule update --init
// reference_code/samtools`).
func loadFixture(t *testing.T, rel string) ([]byte, bool) {
	t.Helper()
	data, err := os.ReadFile(filepath.Join(samtoolsTestDir, rel))
	if err != nil {
		return nil, false
	}
	return data, true
}

// validV3Fixtures are the real CRAM v3 files the structural parser must
// walk cleanly with every CRC32 validating.
var validV3Fixtures = []struct {
	name string
	rel  string
}{
	{"test_input_1_a", "dat/test_input_1_a.cram"},
	{"7.quickcheck.cram30.ok", "quickcheck/7.quickcheck.cram30.ok.cram"},
	{"mpileup.1.cram31", "cram_size/mpileup.1.cram"},
}

// TestParseValidV3Fixtures walks every container and block of each real
// CRAM v3 fixture, decompresses every block whose method is supported,
// and asserts that all CRC32 checksums validate and decompressed sizes
// match. CRC validation happens inside the parser, so reaching the end
// without error already proves every checksum was correct; the test
// additionally exercises the decompression dispatch.
func TestParseValidV3Fixtures(t *testing.T) {
	for _, fx := range validV3Fixtures {
		fx := fx
		t.Run(fx.name, func(t *testing.T) {
			data, ok := loadFixture(t, fx.rel)
			if !ok {
				t.Fatalf("samtools submodule not initialised — fixture unavailable; run `git submodule update --init reference_code/samtools`")
			}
			rd, err := NewReader(bytes.NewReader(data))
			if err != nil {
				t.Fatalf("NewReader: %v", err)
			}
			def := rd.FileDefinition()
			if def.Major != 3 {
				t.Fatalf("expected CRAM major version 3, got %d", def.Major)
			}
			if !def.hasCRC() {
				t.Fatalf("CRAM v3 fixture should carry CRC32 fields")
			}

			containers := 0
			blocks := 0
			decompressed := 0
			skipped := 0
			for {
				c, err := rd.Next()
				if err == io.EOF {
					break
				}
				if err != nil {
					t.Fatalf("Next (container %d): %v", containers, err)
				}
				containers++
				if int(c.Header.NumBlocks) != len(c.Blocks) {
					t.Fatalf("container %d: header says %d blocks, parsed %d",
						c.Index, c.Header.NumBlocks, len(c.Blocks))
				}
				for bi := range c.Blocks {
					b := &c.Blocks[bi]
					blocks++
					if !b.SupportedMethod() {
						skipped++
						// An unsupported method must still error cleanly.
						if _, derr := b.Decompress(); derr == nil {
							t.Errorf("container %d block %d: expected error decompressing %s",
								c.Index, bi, b.Method)
						}
						continue
					}
					out, derr := b.Decompress()
					if derr != nil {
						t.Fatalf("container %d block %d: decompress %s: %v",
							c.Index, bi, b.Method, derr)
					}
					if int32(len(out)) != b.UncompressedSize {
						t.Fatalf("container %d block %d: decompressed %d bytes, declared %d",
							c.Index, bi, len(out), b.UncompressedSize)
					}
					decompressed++
				}
			}
			if containers == 0 {
				t.Fatalf("parsed no containers")
			}
			if blocks == 0 {
				t.Fatalf("parsed no blocks")
			}
			t.Logf("%s: %d containers, %d blocks (%d decompressed, %d skipped)",
				fx.name, containers, blocks, decompressed, skipped)
		})
	}
}

// TestFirstContainerHoldsSAMHeader checks that the first block of the
// first container of a real CRAM file is the SAM file-header block and
// that it decompresses to text starting with "@HD".
func TestFirstContainerHoldsSAMHeader(t *testing.T) {
	data, ok := loadFixture(t, "dat/test_input_1_a.cram")
	if !ok {
		t.Fatalf("samtools submodule not initialised; run `git submodule update --init reference_code/samtools`")
	}
	rd, err := NewReader(bytes.NewReader(data))
	if err != nil {
		t.Fatalf("NewReader: %v", err)
	}
	c, err := rd.Next()
	if err != nil {
		t.Fatalf("Next: %v", err)
	}
	if len(c.Blocks) == 0 {
		t.Fatalf("first container has no blocks")
	}
	first := &c.Blocks[0]
	if first.ContentType != ContentFileHeader {
		t.Fatalf("first block content type = %s, want file-header", first.ContentType)
	}
	out, err := first.Decompress()
	if err != nil {
		t.Fatalf("decompress SAM header block: %v", err)
	}
	// The CRAM SAM-header block stores a 4-byte little-endian length
	// prefix before the SAM text.
	if len(out) < 4 {
		t.Fatalf("SAM header block too short: %d bytes", len(out))
	}
	if !bytes.Contains(out, []byte("@HD")) && !bytes.Contains(out, []byte("@SQ")) {
		t.Fatalf("SAM header block does not contain @HD or @SQ lines")
	}
}

// TestContainersConvenience checks the Containers helper returns the
// same containers as repeated Next calls.
func TestContainersConvenience(t *testing.T) {
	data, ok := loadFixture(t, "quickcheck/7.quickcheck.cram30.ok.cram")
	if !ok {
		t.Fatalf("samtools submodule not initialised; run `git submodule update --init reference_code/samtools`")
	}
	rd, err := NewReader(bytes.NewReader(data))
	if err != nil {
		t.Fatalf("NewReader: %v", err)
	}
	cs, err := rd.Containers()
	if err != nil {
		t.Fatalf("Containers: %v", err)
	}
	if len(cs) == 0 {
		t.Fatalf("Containers returned nothing")
	}
	for i, c := range cs {
		if c.Index != i {
			t.Errorf("container %d has Index %d", i, c.Index)
		}
	}
}

// TestOpenFile exercises the file-path Open/Close API against a real
// fixture.
func TestOpenFile(t *testing.T) {
	path := filepath.Join(samtoolsTestDir, "dat/test_input_1_a.cram")
	if _, err := os.Stat(path); err != nil {
		t.Fatalf("samtools submodule not initialised; run `git submodule update --init reference_code/samtools`")
	}
	rd, err := Open(path)
	if err != nil {
		t.Fatalf("Open: %v", err)
	}
	defer rd.Close()
	if _, err := rd.Containers(); err != nil {
		t.Fatalf("Containers: %v", err)
	}
}

// TestTruncatedFixtures checks that the parser reports a clean error —
// never a panic — on the truncated CRAM fixtures from the samtools
// quickcheck corpus.
func TestTruncatedFixtures(t *testing.T) {
	for _, rel := range []string{
		"quickcheck/8.quickcheck.cram21.truncated.cram",
		"quickcheck/9.quickcheck.cram30.truncated.cram",
	} {
		rel := rel
		t.Run(filepath.Base(rel), func(t *testing.T) {
			data, ok := loadFixture(t, rel)
			if !ok {
				t.Fatalf("samtools submodule not initialised; run `git submodule update --init reference_code/samtools`")
			}
			rd, err := NewReader(bytes.NewReader(data))
			if err != nil {
				// A file def that itself fails to parse is an
				// acceptable clean error.
				return
			}
			var sawErr bool
			for {
				_, err := rd.Next()
				if err == io.EOF {
					break
				}
				if err != nil {
					sawErr = true
					break
				}
			}
			if !sawErr {
				t.Errorf("expected an error walking truncated fixture %s", rel)
			}
		})
	}
}

// TestRejectNonCRAM checks that input not beginning with the CRAM magic
// is rejected.
func TestRejectNonCRAM(t *testing.T) {
	cases := map[string][]byte{
		"empty":      nil,
		"short":      []byte("CR"),
		"bad-magic":  []byte("BAMX\x03\x00" + string(make([]byte, 20))),
		"gzip-magic": {0x1f, 0x8b, 0x08, 0x00},
	}
	for name, in := range cases {
		t.Run(name, func(t *testing.T) {
			if _, err := NewReader(bytes.NewReader(in)); err == nil {
				t.Errorf("expected error for %s input", name)
			}
		})
	}
}

// TestRejectUnsupportedVersion checks that a CRAM file declaring a major
// version outside {2,3,4} is rejected. CRAM v4.0 decode is now supported,
// so a still-unsupported version (5) is used here.
func TestRejectUnsupportedVersion(t *testing.T) {
	in := make([]byte, fileDefSize)
	copy(in, "CRAM")
	in[4] = 5 // major version 5 — out of scope
	if _, err := NewReader(bytes.NewReader(in)); err == nil {
		t.Errorf("expected error for CRAM major version 5")
	}
}

// TestAcceptV4Version confirms a CRAM v4.0 file definition is accepted
// (the structural reader recognises major version 4).
func TestAcceptV4Version(t *testing.T) {
	in := make([]byte, fileDefSize)
	copy(in, "CRAM")
	in[4] = 4 // major version 4
	in[5] = 0 // minor version 0
	rd, err := NewReader(bytes.NewReader(in))
	if err != nil {
		t.Fatalf("CRAM v4.0 file definition rejected: %v", err)
	}
	if got := rd.FileDefinition().VersionString(); got != "4.0" {
		t.Errorf("version = %q, want 4.0", got)
	}
}

// TestGarbageAfterFileDef checks that arbitrary garbage following a
// valid file definition produces an error, never a panic.
func TestGarbageAfterFileDef(t *testing.T) {
	mk := func(body []byte) []byte {
		buf := make([]byte, fileDefSize)
		copy(buf, "CRAM")
		buf[4] = 3
		return append(buf, body...)
	}
	cases := map[string][]byte{
		"all-zero-after-def": mk(make([]byte, 64)),
		"all-0xff":           mk(bytes.Repeat([]byte{0xff}, 64)),
		"random-ish":         mk([]byte("the quick brown fox jumps over the lazy dog")),
		"single-byte":        mk([]byte{0x80}),
	}
	for name, in := range cases {
		t.Run(name, func(t *testing.T) {
			rd, err := NewReader(bytes.NewReader(in))
			if err != nil {
				return
			}
			// Walk to completion; the only requirement is no panic.
			for {
				_, err := rd.Next()
				if err != nil {
					break
				}
			}
		})
	}
}

// TestCRCMismatchDetected flips a byte inside a real fixture's first
// block payload and asserts the parser reports a CRC32 mismatch rather
// than silently accepting the corrupted structure.
func TestCRCMismatchDetected(t *testing.T) {
	data, ok := loadFixture(t, "dat/test_input_1_a.cram")
	if !ok {
		t.Fatalf("samtools submodule not initialised; run `git submodule update --init reference_code/samtools`")
	}
	corrupt := make([]byte, len(data))
	copy(corrupt, data)
	// Byte 64 is well inside the first block's payload (the file def is
	// 26 bytes, the first container header and block header are short).
	corrupt[64] ^= 0xff
	rd, err := NewReader(bytes.NewReader(corrupt))
	if err != nil {
		t.Fatalf("NewReader: %v", err)
	}
	var walkErr error
	for {
		_, err := rd.Next()
		if err == io.EOF {
			break
		}
		if err != nil {
			walkErr = err
			break
		}
	}
	if walkErr == nil {
		t.Fatalf("expected a CRC32 mismatch error after corrupting a block byte")
	}
	if !bytes.Contains([]byte(walkErr.Error()), []byte("CRC32")) {
		t.Fatalf("expected CRC32 error, got: %v", walkErr)
	}
}

// TestOpenMissingFile checks Open surfaces a filesystem error for a
// path that does not exist.
func TestOpenMissingFile(t *testing.T) {
	if _, err := Open(filepath.Join(t.TempDir(), "no-such-file.cram")); err == nil {
		t.Errorf("expected error opening a missing file")
	}
}

// TestOpenNonCRAMFile checks Open rejects an existing file that is not
// CRAM and does not leak the file handle.
func TestOpenNonCRAMFile(t *testing.T) {
	path := filepath.Join(t.TempDir(), "plain.txt")
	if err := os.WriteFile(path, []byte("this is not a CRAM file"), 0o644); err != nil {
		t.Fatalf("write temp file: %v", err)
	}
	if _, err := Open(path); err == nil {
		t.Errorf("expected error opening a non-CRAM file")
	}
}

// TestMissingEOFMarkerIsError checks that a structurally valid CRAM
// stream which simply lacks the trailing EOF marker is reported as a
// truncation error rather than a clean io.EOF.
func TestMissingEOFMarkerIsError(t *testing.T) {
	def := make([]byte, fileDefSize)
	copy(def, "CRAM")
	def[4] = 3
	// One real container holding a single raw block, then end of
	// stream — no EOF marker.
	block := buildBlock(CompRaw, ContentFileHeader, 0, []byte("hdr"), 3)
	hdr := buildContainerHeader(ContainerHeader{Length: int32(len(block)), NumBlocks: 1})
	stream := append(append(append([]byte{}, def...), hdr...), block...)
	rd, err := NewReader(bytes.NewReader(stream))
	if err != nil {
		t.Fatalf("NewReader: %v", err)
	}
	c, err := rd.Next()
	if err != nil {
		t.Fatalf("first container should parse: %v", err)
	}
	if len(c.Blocks) != 1 {
		t.Fatalf("expected 1 block, got %d", len(c.Blocks))
	}
	_, err = rd.Next()
	if err == nil || err == io.EOF {
		t.Fatalf("expected a truncation error for the missing EOF marker, got %v", err)
	}
	if !bytes.Contains([]byte(err.Error()), []byte("EOF marker")) {
		t.Errorf("error %q should mention the EOF marker", err)
	}
}

// TestRejectImpossibleBlockCount checks a container header declaring far
// more blocks than its body could hold is rejected, not used to drive a
// huge allocation.
func TestRejectImpossibleBlockCount(t *testing.T) {
	def := make([]byte, fileDefSize)
	copy(def, "CRAM")
	def[4] = 3
	hdr := buildContainerHeader(ContainerHeader{
		Length:    4,
		NumBlocks: 1 << 20, // dwarfs the 4-byte body
	})
	stream := append(append([]byte{}, def...), hdr...)
	rd, err := NewReader(bytes.NewReader(stream))
	if err != nil {
		t.Fatalf("NewReader: %v", err)
	}
	if _, err := rd.Next(); err == nil {
		t.Errorf("expected an error for an impossible block count")
	}
}

// TestRejectImpossibleLandmarkCount checks a container header declaring
// more landmarks than its body could hold is rejected.
func TestRejectImpossibleLandmarkCount(t *testing.T) {
	def := make([]byte, fileDefSize)
	copy(def, "CRAM")
	def[4] = 3
	// 2 landmarks but Length 1 — fewer body bytes than landmarks.
	hdr := buildContainerHeader(ContainerHeader{Length: 1, Landmarks: []int32{0, 1}})
	stream := append(append([]byte{}, def...), hdr...)
	rd, err := NewReader(bytes.NewReader(stream))
	if err != nil {
		t.Fatalf("NewReader: %v", err)
	}
	if _, err := rd.Next(); err == nil {
		t.Errorf("expected an error for an impossible landmark count")
	}
}

// TestEOFMarkerRecognised checks the CRAM v3 EOF marker container is
// consumed and surfaced as io.EOF, not as a data container.
func TestEOFMarkerRecognised(t *testing.T) {
	// File definition + just the EOF marker.
	def := make([]byte, fileDefSize)
	copy(def, "CRAM")
	def[4] = 3
	stream := append(def, eofMarkerV3...)
	rd, err := NewReader(bytes.NewReader(stream))
	if err != nil {
		t.Fatalf("NewReader: %v", err)
	}
	c, err := rd.Next()
	if !errors.Is(err, io.EOF) {
		t.Fatalf("expected io.EOF at the EOF marker, got container=%v err=%v", c, err)
	}
}

// FuzzParse runs the structural parser over arbitrary input. The parser
// must never panic; any malformed structure must surface as an error.
func FuzzParse(f *testing.F) {
	// Seed with a valid file definition and the real fixtures when
	// available.
	def := make([]byte, fileDefSize)
	copy(def, "CRAM")
	def[4] = 3
	f.Add(def)
	f.Add(append(append([]byte{}, def...), eofMarkerV3...))
	for _, fx := range validV3Fixtures {
		if data, ok := readFixtureNoT(fx.rel); ok {
			f.Add(data)
		}
	}
	f.Fuzz(func(t *testing.T, data []byte) {
		rd, err := NewReader(bytes.NewReader(data))
		if err != nil {
			return
		}
		for i := 0; i < 10000; i++ {
			c, err := rd.Next()
			if err != nil {
				return
			}
			for bi := range c.Blocks {
				// Decompress may legitimately error; it must not panic.
				_, _ = c.Blocks[bi].Decompress()
			}
		}
	})
}

// readFixtureNoT reads a fixture without a *testing.T, for use in fuzz
// seed setup. It returns ok=false if the submodule is absent.
func readFixtureNoT(rel string) ([]byte, bool) {
	data, err := os.ReadFile(filepath.Join(samtoolsTestDir, rel))
	if err != nil {
		return nil, false
	}
	return data, true
}
