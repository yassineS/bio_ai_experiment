package cram

import (
	"errors"
	"io"
	"os"
	"path/filepath"
	"testing"
)

// v31Fixture is a real CRAM v3.1 file from the samtools test corpus. It
// references the full human genome, which is not vendored, so the
// reference-backed bases decode as 'N'; the structural and record-level
// decode paths are exercised regardless.
const v31Fixture = "cram_size/mpileup.1.cram"

// TestV31FileDefinition confirms the v3.1 fixture is recognised as CRAM
// major version 3, minor version 1.
func TestV31FileDefinition(t *testing.T) {
	path := filepath.Join(samtoolsTestDir, v31Fixture)
	if _, err := os.Stat(path); err != nil {
		t.Skipf("samtools submodule not initialised — fixture unavailable; run `git submodule update --init reference_code/samtools`")
	}
	rd, err := Open(path)
	if err != nil {
		t.Fatalf("Open: %v", err)
	}
	defer rd.Close()
	fd := rd.FileDefinition()
	if fd.Major != 3 || fd.Minor != 1 {
		t.Errorf("file definition = CRAM %d.%d, want 3.1", fd.Major, fd.Minor)
	}
}

// TestV31RecordDecode decodes every record of the v3.1 fixture. The v3.1
// record format is shared with v3.0, so the C4b/C5 decoder handles it;
// this locks that in. The reference is unavailable, so a reference-backed
// base decodes as 'N' and NeedsReference reports true — that is the
// documented fallback, not a failure.
func TestV31RecordDecode(t *testing.T) {
	path := filepath.Join(samtoolsTestDir, v31Fixture)
	if _, err := os.Stat(path); err != nil {
		t.Skipf("samtools submodule not initialised — fixture unavailable; run `git submodule update --init reference_code/samtools`")
	}
	rr, err := OpenRecords(path)
	if err != nil {
		t.Fatalf("OpenRecords: %v", err)
	}
	defer rr.Close()
	recs, err := rr.ReadAll()
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if len(recs) == 0 {
		t.Fatal("decoded no records from the v3.1 fixture")
	}
	for i, rec := range recs {
		if rec.QName == "" {
			t.Fatalf("record %d has an empty read name", i)
		}
	}
}

// TestV31UsesRANS4x16 confirms the v3.1 fixture actually exercises the
// rANS 4x16 block codec and that every block decompresses to its
// declared size — proving the v3.1 codec wiring, not just the v3.0
// subset, is on the decode path.
func TestV31UsesRANS4x16(t *testing.T) {
	path := filepath.Join(samtoolsTestDir, v31Fixture)
	if _, err := os.Stat(path); err != nil {
		t.Skipf("samtools submodule not initialised — fixture unavailable; run `git submodule update --init reference_code/samtools`")
	}
	rd, err := Open(path)
	if err != nil {
		t.Fatalf("Open: %v", err)
	}
	defer rd.Close()
	sawRANS4x16 := false
	for {
		c, err := rd.Next()
		if errors.Is(err, io.EOF) {
			break
		}
		if err != nil {
			t.Fatalf("Next: %v", err)
		}
		for i := range c.Blocks {
			b := &c.Blocks[i]
			if b.Method == CompRANS4x16 {
				sawRANS4x16 = true
			}
			out, err := b.Decompress()
			if err != nil {
				t.Fatalf("block (method %d, content id %d): %v", b.Method, b.ContentID, err)
			}
			if int32(len(out)) != b.UncompressedSize {
				t.Fatalf("block decompressed to %d bytes, declared %d", len(out), b.UncompressedSize)
			}
		}
	}
	if !sawRANS4x16 {
		t.Error("the v3.1 fixture used no rANS 4x16 block — the 4x16 decode path was not exercised")
	}
}

// TestNewV31CodecsRejected confirms that the v3.1 block codecs all
// surface a clear decode error on garbage input rather than mis-decoding
// or panicking. The arith_dynamic (method 6), fqzcomp (method 7) and
// name-tokeniser (method 8) codecs are all implemented now; this test
// guards their error paths, not their absence.
func TestNewV31CodecsRejected(t *testing.T) {
	for _, m := range []CompressionMethod{CompArith, CompFQZComp, CompNameTok} {
		b := Block{Method: m, UncompressedSize: 4, Data: []byte{1, 2, 3, 4}}
		out, err := b.Decompress()
		if err == nil {
			t.Errorf("method %d: expected a decode error on garbage data, decoded %d bytes", m, len(out))
		}
	}
}
