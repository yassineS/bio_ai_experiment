package sam

import (
	"bytes"
	"encoding/binary"
	"io"
	"os"
	"path/filepath"
	"testing"

	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
)

// TestBAMWriterUpstreamOracle pins the BAMWriter's BGZF block-boundary policy
// against golden BAM files produced by a genuine upstream samtools
// (`samtools view -b --no-PG`).
//
// Upstream htslib forces specific BGZF block boundaries that a naive writer
// does not reproduce:
//
//   - sam_hdr_write calls bgzf_flush after the BAM header (sam.c:2053), so the
//     header always occupies its own block(s) and never shares a block with
//     the first alignment record.
//   - each record write calls bgzf_flush_try(4 + block_len) (sam.c:889) before
//     emitting the record, so a record never straddles a BGZF block boundary —
//     if it would overflow the current block, that block is flushed first and
//     the record starts a fresh one.
//
// The committed goldens were generated with `--no-PG` so the header text
// matches what BAMWriter emits (our `view` does not inject samtools' own @PG
// provenance line). We assert two compression-independent invariants that the
// bug-fix is about:
//
//  1. the decoded (uncompressed) BAM byte stream is byte-for-byte identical,
//     which proves the record encoding (including bin numbers for unmapped
//     reads, see reg2bin) matches upstream; and
//  2. the per-block *uncompressed* sizes match, which proves the block
//     boundaries land exactly where htslib puts them.
//
// The DEFLATE payload bytes themselves are owned and byte-locked by the BGZF
// layer's own oracle (pkg/htsgo/bgzf/oracle_test.go) against libdeflate-linked
// htslib; they are deliberately not re-checked here because the golden may have
// been produced by a zlib-linked samtools, which compresses the identical
// payload to different bytes.
func TestBAMWriterUpstreamOracle(t *testing.T) {
	cases := []struct {
		name string
		// wantBlocks is the expected number of BGZF blocks excluding the EOF
		// block; it documents the boundary structure each fixture exercises.
		minDataBlocks int
	}{
		{name: "small", minDataBlocks: 2}, // header block + one record block
		{name: "multi", minDataBlocks: 3}, // header block + multiple record blocks (FlushTry mid-stream)
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			samPath := filepath.Join("testdata", "oracle", tc.name+".sam")
			bamPath := filepath.Join("testdata", "oracle", tc.name+".bam")

			samFile, err := os.Open(samPath)
			if err != nil {
				t.Fatalf("open sam fixture: %v", err)
			}
			defer samFile.Close()

			golden, err := os.ReadFile(bamPath)
			if err != nil {
				t.Fatalf("read golden bam: %v", err)
			}

			// Round-trip the SAM through our BAMWriter.
			r, err := NewSAMReader(samFile)
			if err != nil {
				t.Fatalf("open sam reader: %v", err)
			}
			hdr := r.Header()
			var out bytes.Buffer
			w := NewBAMWriter(&out)
			if err := w.WriteHeader(hdr); err != nil {
				t.Fatalf("write bam header: %v", err)
			}
			for {
				rec, err := r.Read()
				if err == io.EOF {
					break
				}
				if err != nil {
					t.Fatalf("read sam record: %v", err)
				}
				if err := w.Write(rec); err != nil {
					t.Fatalf("write bam record: %v", err)
				}
			}
			if err := w.Close(); err != nil {
				t.Fatalf("close bam writer: %v", err)
			}
			ours := out.Bytes()

			// (1) decoded BAM stream must be byte-identical.
			goldStream := mustDecodeBGZF(t, golden)
			ourStream := mustDecodeBGZF(t, ours)
			if !bytes.Equal(goldStream, ourStream) {
				diff := firstDiff(goldStream, ourStream)
				t.Fatalf("decoded BAM stream differs from upstream: golden %d bytes, ours %d bytes, first diff at byte %d",
					len(goldStream), len(ourStream), diff)
			}

			// (2) per-block uncompressed sizes (block boundaries) must match.
			goldSizes := blockUncompressedSizes(t, golden)
			ourSizes := blockUncompressedSizes(t, ours)
			if !equalInts(goldSizes, ourSizes) {
				t.Fatalf("BGZF block boundary layout differs from upstream:\n golden uncompressed block sizes = %v\n ours   uncompressed block sizes = %v",
					goldSizes, ourSizes)
			}
			// data blocks = all blocks minus the trailing EOF (size 0) block.
			dataBlocks := len(ourSizes)
			if dataBlocks > 0 && ourSizes[dataBlocks-1] == 0 {
				dataBlocks--
			}
			if dataBlocks < tc.minDataBlocks {
				t.Fatalf("fixture %q produced %d data blocks, expected at least %d (boundary policy under-exercised)",
					tc.name, dataBlocks, tc.minDataBlocks)
			}
		})
	}
}

// mustDecodeBGZF decompresses a complete BGZF stream and returns the
// concatenated uncompressed bytes.
func mustDecodeBGZF(t *testing.T, b []byte) []byte {
	t.Helper()
	br, err := bgzip.NewReader(bytes.NewReader(b))
	if err != nil {
		t.Fatalf("open bgzf reader: %v", err)
	}
	defer br.Close()
	out, err := io.ReadAll(br)
	if err != nil {
		t.Fatalf("decode bgzf: %v", err)
	}
	return out
}

// blockUncompressedSizes walks a BGZF stream block-by-block (without
// decompressing) and returns the ISIZE (uncompressed size) of each block in
// order, including the trailing 0-sized EOF block. The per-block boundaries
// are exactly what the header-flush + per-record flush_try policy controls.
func blockUncompressedSizes(t *testing.T, b []byte) []int {
	t.Helper()
	var sizes []int
	off := 0
	for off < len(b) {
		if off+12 > len(b) || b[off] != 0x1f || b[off+1] != 0x8b {
			t.Fatalf("bad gzip magic at offset %d", off)
		}
		xlen := int(binary.LittleEndian.Uint16(b[off+10 : off+12]))
		extra := b[off+12 : off+12+xlen]
		bsize, ok := findBSIZE(extra)
		if !ok {
			t.Fatalf("BC subfield missing at offset %d", off)
		}
		blockLen := int(bsize) + 1
		isize := int(binary.LittleEndian.Uint32(b[off+blockLen-4 : off+blockLen]))
		sizes = append(sizes, isize)
		off += blockLen
	}
	return sizes
}

// findBSIZE returns the BGZF BC-subfield BSIZE value from a gzip extra field.
func findBSIZE(extra []byte) (uint16, bool) {
	for len(extra) >= 4 {
		slen := int(binary.LittleEndian.Uint16(extra[2:4]))
		if 4+slen > len(extra) {
			return 0, false
		}
		if extra[0] == 'B' && extra[1] == 'C' && slen == 2 {
			return binary.LittleEndian.Uint16(extra[4:6]), true
		}
		extra = extra[4+slen:]
	}
	return 0, false
}

func firstDiff(a, b []byte) int {
	n := len(a)
	if len(b) < n {
		n = len(b)
	}
	for i := 0; i < n; i++ {
		if a[i] != b[i] {
			return i
		}
	}
	return n
}

func equalInts(a, b []int) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}
