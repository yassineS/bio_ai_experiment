package bam

import (
	"bytes"
	"encoding/binary"
	"io"
	"testing"
)

// TestBAIWriteReadRoundTrip exercises a tiny hand-built index through the
// writer and back through the reader.
func TestBAIWriteReadRoundTrip(t *testing.T) {
	idx := &BAIIndex{
		NoCoor: 7,
		Refs: []BAIRef{
			{
				Bins: []BAIBin{
					{BinID: 4681, Chunks: []BAIChunk{{Beg: 100, End: 200}, {Beg: 300, End: 400}}},
					{BinID: MetaBin, Chunks: []BAIChunk{{Beg: 100, End: 400}, {Beg: 12, End: 3}}},
				},
				Linear: []uint64{0, 100, 100, 300},
			},
			// Second reference is empty.
			{Bins: nil, Linear: nil},
		},
	}
	var buf bytes.Buffer
	if err := WriteBAI(&buf, idx); err != nil {
		t.Fatalf("WriteBAI: %v", err)
	}

	// The raw bytes should start with "BAI\1" and contain the n_no_coor
	// trailer at the end.
	raw := buf.Bytes()
	if !bytes.Equal(raw[:4], BAIMagic[:]) {
		t.Fatalf("magic: got % x, want % x", raw[:4], BAIMagic)
	}
	if got := binary.LittleEndian.Uint64(raw[len(raw)-8:]); got != 7 {
		t.Errorf("n_no_coor trailer: got %d, want 7", got)
	}

	got, err := ReadBAI(bytes.NewReader(raw))
	if err != nil {
		t.Fatalf("ReadBAI: %v", err)
	}
	if got.NoCoor != 7 {
		t.Errorf("NoCoor round-trip: got %d, want 7", got.NoCoor)
	}
	if len(got.Refs) != 2 {
		t.Fatalf("Refs len: got %d, want 2", len(got.Refs))
	}
	if len(got.Refs[0].Bins) != 2 {
		t.Errorf("bins[0] len: got %d, want 2", len(got.Refs[0].Bins))
	}
	if got.Refs[0].Linear[3] != 300 {
		t.Errorf("linear[3]: got %d, want 300", got.Refs[0].Linear[3])
	}
}

// TestBAIBadMagic confirms a non-BAI input is rejected.
func TestBAIBadMagic(t *testing.T) {
	bad := bytes.NewReader([]byte{'N', 'O', 'T', 'B', 0, 0, 0, 0})
	if _, err := ReadBAI(bad); err != ErrBadBAIMagic {
		t.Errorf("expected ErrBadBAIMagic, got %v", err)
	}
}

// TestBAIBuilderMetaPseudoBin verifies the builder injects the htslib
// pseudo-bin (37450) with exactly two chunks: the first holding (firstVOff,
// lastVOff), the second holding (mapped, unmapped).
func TestBAIBuilderMetaPseudoBin(t *testing.T) {
	b := NewBAIBuilder(1)
	// Three mapped records on ref 0 + one unmapped-but-placed on ref 0.
	if err := b.AddRecord(0, 100, 200, 0x100, 0x200, true); err != nil {
		t.Fatal(err)
	}
	if err := b.AddRecord(0, 300, 400, 0x250, 0x300, true); err != nil {
		t.Fatal(err)
	}
	if err := b.AddRecord(0, 500, 600, 0x320, 0x400, true); err != nil {
		t.Fatal(err)
	}
	if err := b.AddRecord(0, 700, 800, 0x420, 0x500, false); err != nil {
		t.Fatal(err)
	}
	idx := b.Finish()
	first, last, ok := idx.MetaBounds(0)
	if !ok {
		t.Fatal("MetaBounds: missing pseudo-bin")
	}
	if first != 0x100 || last != 0x500 {
		t.Errorf("MetaBounds: got (%#x, %#x), want (0x100, 0x500)", first, last)
	}
	mapped, unmapped, ok := idx.MetaCounts(0)
	if !ok {
		t.Fatal("MetaCounts: missing pseudo-bin")
	}
	if mapped != 3 || unmapped != 1 {
		t.Errorf("MetaCounts: got mapped=%d unmapped=%d, want 3,1", mapped, unmapped)
	}
}

// TestBAIBuilderNoCoor verifies refID == -1 records bump NoCoor and do not
// land in any reference's bin/linear index.
func TestBAIBuilderNoCoor(t *testing.T) {
	b := NewBAIBuilder(1)
	b.AddRecord(-1, 0, 0, 0, 0, false)
	b.AddRecord(-1, 0, 0, 0, 0, false)
	b.AddRecord(-1, 0, 0, 0, 0, false)
	idx := b.Finish()
	if idx.NoCoor != 3 {
		t.Errorf("NoCoor: got %d, want 3", idx.NoCoor)
	}
	if len(idx.Refs[0].Bins) != 0 {
		t.Errorf("expected empty bins on ref 0, got %d", len(idx.Refs[0].Bins))
	}
}

// TestBAIBuilderLinearMultiTile builds a record spanning multiple 16-Kbp
// tiles and confirms every tile records the same virtual offset.
func TestBAIBuilderLinearMultiTile(t *testing.T) {
	b := NewBAIBuilder(1)
	// A record from pos 0 spanning 50,000 bp covers tiles 0, 1, 2.
	if err := b.AddRecord(0, 0, 50_000, 0xABCDEF, 0xABCDFF, true); err != nil {
		t.Fatal(err)
	}
	idx := b.Finish()
	lin := idx.Refs[0].Linear
	if len(lin) < 3 {
		t.Fatalf("expected ≥3 linear entries, got %d", len(lin))
	}
	for i := 0; i < 3; i++ {
		if lin[i] != 0xABCDEF {
			t.Errorf("linear[%d]: got %#x, want 0xABCDEF", i, lin[i])
		}
	}
}

// TestBAIBuilderChunkMerge confirms consecutive records in the same bin
// have their chunks coalesced.
func TestBAIBuilderChunkMerge(t *testing.T) {
	b := NewBAIBuilder(1)
	// Both records map to the same level-5 bin (tile 0). Their virtual
	// offsets are contiguous, so the chunk list should merge to a single
	// (Beg, End) span.
	b.AddRecord(0, 100, 110, 0x100, 0x200, true)
	b.AddRecord(0, 120, 130, 0x200, 0x300, true)
	idx := b.Finish()
	// Find the level-5 bin (4681 covers [0, 16K)).
	var found *BAIBin
	for i, bin := range idx.Refs[0].Bins {
		if bin.BinID == 4681 {
			found = &idx.Refs[0].Bins[i]
			break
		}
	}
	if found == nil {
		t.Fatal("missing bin 4681")
	}
	if len(found.Chunks) != 1 {
		t.Errorf("expected merged single chunk, got %d", len(found.Chunks))
	}
	if found.Chunks[0].Beg != 0x100 || found.Chunks[0].End != 0x300 {
		t.Errorf("merged chunk: got (%#x, %#x), want (0x100, 0x300)", found.Chunks[0].Beg, found.Chunks[0].End)
	}
}

// TestBAIBuilderRegionChunksLinearTrim confirms that the linear-index lower
// bound trims chunks below the relevant tile.
func TestBAIBuilderRegionChunksLinearTrim(t *testing.T) {
	b := NewBAIBuilder(1)
	// Two records — one at tile 0, one at tile 5 (≈81,000 bp).
	b.AddRecord(0, 100, 200, 0x100, 0x200, true)
	b.AddRecord(0, 81_000, 82_000, 0x500, 0x600, true)
	idx := b.Finish()
	// Region covering tile 5 only: minOff should be 0x500.
	chunks := idx.RegionChunks(0, 81_000, 82_000)
	if len(chunks) != 1 {
		t.Fatalf("expected 1 chunk, got %d", len(chunks))
	}
	if chunks[0].Beg != 0x500 {
		t.Errorf("expected Beg=0x500 after linear trim, got %#x", chunks[0].Beg)
	}
}

// TestBAIRegionChunksUnknownRef ensures querying an out-of-range refID
// returns nil cleanly.
func TestBAIRegionChunksUnknownRef(t *testing.T) {
	b := NewBAIBuilder(1)
	b.AddRecord(0, 0, 100, 1, 2, true)
	idx := b.Finish()
	if got := idx.RegionChunks(7, 0, 1000); got != nil {
		t.Errorf("expected nil for unknown ref, got %v", got)
	}
	if got := idx.RegionChunks(0, 50, 10); got != nil {
		t.Errorf("expected nil for inverted range, got %v", got)
	}
}

// TestBAIBuilderUnknownRef confirms AddRecord rejects refID >= numRefs.
func TestBAIBuilderUnknownRef(t *testing.T) {
	b := NewBAIBuilder(2)
	if err := b.AddRecord(2, 0, 1, 1, 2, true); err == nil {
		t.Error("expected error for refID >= numRefs")
	}
}

// TestBAIReadEOFTooEarly ensures truncated BAI input is detected.
func TestBAIReadEOFTooEarly(t *testing.T) {
	// Just magic + a partial n_ref field.
	var buf bytes.Buffer
	buf.Write(BAIMagic[:])
	buf.Write([]byte{0x00, 0x00}) // only 2 bytes of the 4-byte n_ref
	if _, err := ReadBAI(&buf); err == nil || err == io.EOF {
		// Either io.EOF or wrapped error is fine; we just want a non-nil error.
		if err == nil {
			t.Error("expected non-nil error on truncated BAI")
		}
	}
}

// TestBAIByteLayout decodes a hand-built BAI fixture and verifies every
// field byte-by-byte. This keeps WriteBAI's wire format pinned even when
// nobody is looking.
func TestBAIByteLayout(t *testing.T) {
	idx := &BAIIndex{
		Refs: []BAIRef{
			{
				Bins: []BAIBin{
					{BinID: 4681, Chunks: []BAIChunk{{Beg: 100, End: 200}}},
				},
				Linear: []uint64{50},
			},
		},
	}
	var buf bytes.Buffer
	if err := WriteBAI(&buf, idx); err != nil {
		t.Fatalf("WriteBAI: %v", err)
	}
	raw := buf.Bytes()
	// Expected layout:
	//   4 magic + 4 n_ref + (4 n_bin + 4 bin_id + 4 n_chunk + 16 chunk + 4 n_intv + 8 intv) per ref
	//   + 8 n_no_coor trailer (htslib always emits it, even when zero)
	//   = 4 + 4 + (4 + 4 + 4 + 16 + 4 + 8) + 8 = 56 bytes
	if len(raw) != 56 {
		t.Fatalf("expected 56 raw bytes, got %d (% x)", len(raw), raw)
	}
	if string(raw[0:4]) != "BAI\x01" {
		t.Errorf("magic: got %q", raw[0:4])
	}
	if binary.LittleEndian.Uint32(raw[4:8]) != 1 {
		t.Errorf("n_ref: got %d, want 1", binary.LittleEndian.Uint32(raw[4:8]))
	}
	if binary.LittleEndian.Uint32(raw[8:12]) != 1 {
		t.Errorf("n_bin: got %d, want 1", binary.LittleEndian.Uint32(raw[8:12]))
	}
	if binary.LittleEndian.Uint32(raw[12:16]) != 4681 {
		t.Errorf("bin_id: got %d, want 4681", binary.LittleEndian.Uint32(raw[12:16]))
	}
	if binary.LittleEndian.Uint32(raw[16:20]) != 1 {
		t.Errorf("n_chunk: got %d, want 1", binary.LittleEndian.Uint32(raw[16:20]))
	}
	if binary.LittleEndian.Uint64(raw[20:28]) != 100 {
		t.Errorf("chunk_beg: got %d, want 100", binary.LittleEndian.Uint64(raw[20:28]))
	}
	if binary.LittleEndian.Uint64(raw[28:36]) != 200 {
		t.Errorf("chunk_end: got %d, want 200", binary.LittleEndian.Uint64(raw[28:36]))
	}
	if binary.LittleEndian.Uint32(raw[36:40]) != 1 {
		t.Errorf("n_intv: got %d, want 1", binary.LittleEndian.Uint32(raw[36:40]))
	}
	if binary.LittleEndian.Uint64(raw[40:48]) != 50 {
		t.Errorf("intv[0]: got %d, want 50", binary.LittleEndian.Uint64(raw[40:48]))
	}
	if binary.LittleEndian.Uint64(raw[48:56]) != 0 {
		t.Errorf("n_no_coor trailer: got %d, want 0", binary.LittleEndian.Uint64(raw[48:56]))
	}
}

// TestBAIBuilderRegionChunksFromEmpty checks that calling RegionChunks when
// no records have been added returns nil.
func TestBAIBuilderRegionChunksFromEmpty(t *testing.T) {
	idx := NewBAIBuilder(1).Finish()
	if got := idx.RegionChunks(0, 0, 1000); got != nil {
		t.Errorf("expected nil chunks on empty index, got %v", got)
	}
}

// TestBAIMetaAccessorsMissing exercises the negative paths in MetaCounts/MetaBounds.
func TestBAIMetaAccessorsMissing(t *testing.T) {
	idx := &BAIIndex{Refs: []BAIRef{{}}}
	if _, _, ok := idx.MetaCounts(0); ok {
		t.Error("expected MetaCounts(absent) → false")
	}
	if _, _, ok := idx.MetaBounds(0); ok {
		t.Error("expected MetaBounds(absent) → false")
	}
	if _, _, ok := idx.MetaCounts(-1); ok {
		t.Error("expected MetaCounts(bad refID) → false")
	}
	if _, _, ok := idx.MetaBounds(99); ok {
		t.Error("expected MetaBounds(bad refID) → false")
	}
}
