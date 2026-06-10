package bgzf

import (
	"bytes"
	"fmt"
	"io"
	"math/rand"
	"os"
	"os/exec"
	"path/filepath"
	"sync"
	"testing"
)

// upstreamHtslibBgzipGzi locates (building if necessary) the htslib `bgzip`
// binary vendored under reference_code/htslib. The build runs at most once per
// test process via sync.Once. These tests assert byte-for-byte parity against
// the live upstream binary, so a missing/unbuildable binary is a hard failure
// (t.Fatalf), never a skip.
var (
	upstreamBgzipGziOnce sync.Once
	upstreamBgzipGziPath string
	upstreamBgzipGziErr  error
)

func upstreamHtslibBgzipGzi(t *testing.T) string {
	t.Helper()
	upstreamBgzipGziOnce.Do(func() {
		htslibDir, err := filepath.Abs("../../../reference_code/htslib")
		if err != nil {
			upstreamBgzipGziErr = err
			return
		}
		bin := filepath.Join(htslibDir, "bgzip")
		if _, statErr := os.Stat(bin); statErr == nil {
			upstreamBgzipGziPath = bin
			return
		}
		// Configure + build bgzip if it is not already present.
		if _, statErr := os.Stat(filepath.Join(htslibDir, "config.mk")); statErr != nil {
			for _, args := range [][]string{
				{"autoreconf", "-i"},
				{"./configure"},
			} {
				cmd := exec.Command(args[0], args[1:]...)
				cmd.Dir = htslibDir
				if out, runErr := cmd.CombinedOutput(); runErr != nil {
					upstreamBgzipGziErr = fmt.Errorf("%v: %v\n%s", args, runErr, out)
					return
				}
			}
		}
		cmd := exec.Command("make", "-j4", "bgzip")
		cmd.Dir = htslibDir
		if out, runErr := cmd.CombinedOutput(); runErr != nil {
			upstreamBgzipGziErr = fmt.Errorf("make bgzip: %v\n%s", runErr, out)
			return
		}
		upstreamBgzipGziPath = bin
	})
	if upstreamBgzipGziErr != nil {
		t.Fatalf("locating/building upstream bgzip: %v", upstreamBgzipGziErr)
	}
	if upstreamBgzipGziPath == "" {
		t.Fatalf("upstream bgzip not available")
	}
	return upstreamBgzipGziPath
}

// makeReproduciblePayload returns a deterministic but varied payload large
// enough to span many BGZF blocks (MaxBlockSize is ~64 KiB).
func makeReproduciblePayload(n int) []byte {
	rng := rand.New(rand.NewSource(0x5eed))
	out := make([]byte, n)
	// Mix of compressible runs and random bytes so blocks differ in size.
	for i := 0; i < n; {
		if rng.Intn(3) == 0 {
			run := rng.Intn(200)
			b := byte('A' + rng.Intn(26))
			for j := 0; j < run && i < n; j++ {
				out[i] = b
				i++
			}
		} else {
			out[i] = byte(rng.Intn(256))
			i++
		}
	}
	return out
}

// bgzipCompress writes payload to a new .gz file at path using our own Writer.
func bgzipCompress(t *testing.T, path string, payload []byte) {
	t.Helper()
	f, err := os.Create(path)
	if err != nil {
		t.Fatalf("create %s: %v", path, err)
	}
	defer f.Close()
	w := NewWriter(f)
	if _, err := w.Write(payload); err != nil {
		t.Fatalf("bgzf write: %v", err)
	}
	if err := w.Close(); err != nil {
		t.Fatalf("bgzf close: %v", err)
	}
}

// TestSeekReader_ParityWithUpstreamBgzipExtraction confirms that our
// SeekReader, driven by an upstream-`bgzip -r`-produced .gzi, extracts the
// exact same region bytes that `bgzip -b N -s M` produces.
func TestSeekReader_ParityWithUpstreamBgzipExtraction(t *testing.T) {
	bgzip := upstreamHtslibBgzipGzi(t)
	dir := t.TempDir()

	payload := makeReproduciblePayload(300 * 1024)
	gzPath := filepath.Join(dir, "data.gz")
	bgzipCompress(t, gzPath, payload)

	// Have upstream build the .gzi via `bgzip -r`.
	reidx := exec.Command(bgzip, "-r", gzPath)
	if out, err := reidx.CombinedOutput(); err != nil {
		t.Fatalf("bgzip -r: %v\n%s", err, out)
	}
	gziPath := gzPath + ".gzi"
	gziBytes, err := os.ReadFile(gziPath)
	if err != nil {
		t.Fatalf("read upstream .gzi: %v", err)
	}
	index, err := ReadGZI(bytes.NewReader(gziBytes))
	if err != nil {
		t.Fatalf("ReadGZI(upstream): %v", err)
	}

	gzFile, err := os.Open(gzPath)
	if err != nil {
		t.Fatalf("open gz: %v", err)
	}
	defer gzFile.Close()

	regions := []struct{ off, size int64 }{
		{0, 100},
		{1, 1},
		{1000, 5000},
		{65000, 2000}, // straddles a block boundary
		{int64(len(payload)) - 10, 10},
		{131072, 131072}, // spans multiple full blocks
		{50, 0},          // size 0 → to end of stream
	}
	for _, rg := range regions {
		rg := rg
		t.Run(fmt.Sprintf("off%d_sz%d", rg.off, rg.size), func(t *testing.T) {
			// Upstream extraction via bgzip -b/-s.
			args := []string{"-c", "-b", fmt.Sprint(rg.off)}
			if rg.size > 0 {
				args = append(args, "-s", fmt.Sprint(rg.size))
			}
			args = append(args, gzPath)
			cmd := exec.Command(bgzip, args...)
			want, err := cmd.Output()
			if err != nil {
				t.Fatalf("upstream bgzip %v: %v", args, err)
			}

			sr := NewSeekReader(gzFile, index)
			got, err := sr.ReadRegion(rg.off, rg.size)
			if err != nil {
				t.Fatalf("ReadRegion(%d,%d): %v", rg.off, rg.size, err)
			}
			if !bytes.Equal(got, want) {
				t.Fatalf("region off=%d size=%d mismatch: got %d bytes, want %d bytes\n got[:32]=%x\nwant[:32]=%x",
					rg.off, rg.size, len(got), len(want), prefix(got), prefix(want))
			}
			// Cross-check against the raw payload too.
			end := rg.off + rg.size
			if rg.size == 0 || end > int64(len(payload)) {
				end = int64(len(payload))
			}
			if !bytes.Equal(got, payload[rg.off:end]) {
				t.Fatalf("region off=%d size=%d does not match source payload slice", rg.off, rg.size)
			}
		})
	}
}

// TestWriteGZI_UsableByUpstreamBgzip confirms a .gzi produced by our WriteGZI
// is accepted by upstream `bgzip -b/-s` for partial extraction — i.e. our
// writer is byte-compatible with htslib's index format.
func TestWriteGZI_UsableByUpstreamBgzip(t *testing.T) {
	bgzip := upstreamHtslibBgzipGzi(t)
	dir := t.TempDir()

	payload := makeReproduciblePayload(250 * 1024)
	gzPath := filepath.Join(dir, "ours.gz")
	bgzipCompress(t, gzPath, payload)

	// Scan our own compressed file and write a .gzi with WriteGZI.
	gzFile, err := os.Open(gzPath)
	if err != nil {
		t.Fatalf("open gz: %v", err)
	}
	offsets, err := Scan(gzFile)
	if err != nil {
		t.Fatalf("Scan: %v", err)
	}
	gzFile.Close()

	gziPath := gzPath + ".gzi"
	gziFile, err := os.Create(gziPath)
	if err != nil {
		t.Fatalf("create gzi: %v", err)
	}
	if err := WriteGZI(gziFile, offsets); err != nil {
		t.Fatalf("WriteGZI: %v", err)
	}
	gziFile.Close()

	// Verify upstream's own `bgzip -r` would produce a byte-identical .gzi.
	cp := filepath.Join(dir, "ref.gz")
	if err := copyFile(gzPath, cp); err != nil {
		t.Fatalf("copy: %v", err)
	}
	if out, err := exec.Command(bgzip, "-r", cp).CombinedOutput(); err != nil {
		t.Fatalf("bgzip -r: %v\n%s", err, out)
	}
	refGzi, err := os.ReadFile(cp + ".gzi")
	if err != nil {
		t.Fatalf("read ref gzi: %v", err)
	}
	ourGzi, err := os.ReadFile(gziPath)
	if err != nil {
		t.Fatalf("read our gzi: %v", err)
	}
	if !bytes.Equal(ourGzi, refGzi) {
		t.Fatalf("WriteGZI bytes differ from upstream bgzip -r:\n our=%x\n ref=%x", ourGzi, refGzi)
	}

	// And confirm upstream can extract using our index.
	for _, rg := range []struct{ off, size int64 }{{0, 64}, {70000, 1234}, {int64(len(payload)) - 5, 5}} {
		want := payload[rg.off : rg.off+rg.size]
		out, err := exec.Command(bgzip, "-c", "-b", fmt.Sprint(rg.off), "-s", fmt.Sprint(rg.size), gzPath).Output()
		if err != nil {
			t.Fatalf("upstream extract with our gzi: %v", err)
		}
		if !bytes.Equal(out, want) {
			t.Fatalf("upstream extract off=%d size=%d mismatch using our gzi", rg.off, rg.size)
		}
	}
}

// TestSeekReader_ParityWithUpstreamOnUpstreamCompressed compresses the payload
// with upstream `bgzip` (not our Writer) and confirms our SeekReader still
// extracts identically — exercising parity on a third-party block layout.
func TestSeekReader_ParityWithUpstreamOnUpstreamCompressed(t *testing.T) {
	bgzip := upstreamHtslibBgzipGzi(t)
	dir := t.TempDir()

	payload := makeReproduciblePayload(180 * 1024)
	rawPath := filepath.Join(dir, "raw.bin")
	if err := os.WriteFile(rawPath, payload, 0o644); err != nil {
		t.Fatalf("write raw: %v", err)
	}
	// `bgzip -i` compresses and writes the .gzi in one shot.
	if out, err := exec.Command(bgzip, "-i", rawPath).CombinedOutput(); err != nil {
		t.Fatalf("bgzip -i: %v\n%s", err, out)
	}
	gzPath := rawPath + ".gz"
	gziBytes, err := os.ReadFile(gzPath + ".gzi")
	if err != nil {
		t.Fatalf("read gzi: %v", err)
	}
	index, err := ReadGZI(bytes.NewReader(gziBytes))
	if err != nil {
		t.Fatalf("ReadGZI: %v", err)
	}
	gzFile, err := os.Open(gzPath)
	if err != nil {
		t.Fatalf("open: %v", err)
	}
	defer gzFile.Close()

	for _, off := range []int64{0, 999, 65535, 65536, int64(len(payload)) - 1} {
		sr := NewSeekReader(gzFile, index)
		got, err := sr.ReadRegion(off, 256)
		if err != nil {
			t.Fatalf("ReadRegion off=%d: %v", off, err)
		}
		end := off + 256
		if end > int64(len(payload)) {
			end = int64(len(payload))
		}
		if !bytes.Equal(got, payload[off:end]) {
			t.Fatalf("off=%d mismatch on upstream-compressed file", off)
		}
	}
}

func prefix(b []byte) []byte {
	if len(b) > 32 {
		return b[:32]
	}
	return b
}

func copyFile(src, dst string) error {
	data, err := os.ReadFile(src)
	if err != nil {
		return err
	}
	return os.WriteFile(dst, data, 0o644)
}

// TestSeekReader_Unit exercises the SeekReader without any upstream binary,
// using our own Writer to produce a multi-block stream and a Scan-derived
// index. It verifies seek-anywhere correctness, sparse-index fallback (where
// the index omits some blocks), and the past-end error.
func TestSeekReader_Unit(t *testing.T) {
	payload := makeReproduciblePayload(200 * 1024)
	var buf bytes.Buffer
	w := NewWriter(&buf)
	if _, err := w.Write(payload); err != nil {
		t.Fatalf("write: %v", err)
	}
	if err := w.Close(); err != nil {
		t.Fatalf("close: %v", err)
	}
	compressed := buf.Bytes()
	ra := bytes.NewReader(compressed)

	full, err := Scan(bytes.NewReader(compressed))
	if err != nil {
		t.Fatalf("Scan: %v", err)
	}
	// Build a .gzi-style index (drop the implicit (0,0) like ReadGZI does).
	index := make([]BlockOffset, len(full))
	copy(index, full)
	if len(index) > 0 && index[0].CompressedOffset == 0 {
		index = index[1:]
	}

	t.Run("seek_anywhere", func(t *testing.T) {
		for _, off := range []int64{0, 1, 100, 65279, 65280, 65281, 130559, int64(len(payload)) - 1} {
			sr := NewSeekReader(ra, index)
			got, err := sr.ReadRegion(off, 300)
			if err != nil {
				t.Fatalf("off=%d: %v", off, err)
			}
			end := off + 300
			if end > int64(len(payload)) {
				end = int64(len(payload))
			}
			if !bytes.Equal(got, payload[off:end]) {
				t.Fatalf("off=%d mismatch", off)
			}
		}
	})

	t.Run("sparse_index_fallback", func(t *testing.T) {
		// Keep only every other index entry so Read must walk the compressed
		// stream to find unindexed blocks.
		var sparse []BlockOffset
		for i, e := range index {
			if i%2 == 0 {
				sparse = append(sparse, e)
			}
		}
		sr := NewSeekReader(ra, sparse)
		got, err := sr.ReadRegion(0, int64(len(payload)))
		if err != nil {
			t.Fatalf("read all with sparse index: %v", err)
		}
		if !bytes.Equal(got, payload) {
			t.Fatalf("sparse-index full read mismatch (%d vs %d bytes)", len(got), len(payload))
		}
	})

	t.Run("empty_index_fallback", func(t *testing.T) {
		// With no index entries at all (only the implicit (0,0) block),
		// every block past the first must be located by walking the
		// compressed stream forward. This pins the regression where
		// advanceBlock re-derived the current block's compressed offset
		// from the index and never advanced past the first block.
		sr := NewSeekReader(ra, nil)
		got, err := sr.ReadRegion(0, int64(len(payload)))
		if err != nil {
			t.Fatalf("read all with empty index: %v", err)
		}
		if !bytes.Equal(got, payload) {
			t.Fatalf("empty-index full read mismatch (%d vs %d bytes)", len(got), len(payload))
		}
		// And a seek into a late block must work too.
		sr2 := NewSeekReader(ra, nil)
		off := int64(len(payload)) - 777
		got2, err := sr2.ReadRegion(off, 777)
		if err != nil {
			t.Fatalf("empty-index late seek: %v", err)
		}
		if !bytes.Equal(got2, payload[off:]) {
			t.Fatalf("empty-index late seek mismatch")
		}
	})

	t.Run("seek_to_end", func(t *testing.T) {
		sr := NewSeekReader(ra, index)
		if err := sr.SeekTo(int64(len(payload))); err != nil {
			t.Fatalf("seek to exact end: %v", err)
		}
		b := make([]byte, 8)
		if _, err := sr.Read(b); err != io.EOF {
			t.Fatalf("read at end: want io.EOF, got %v", err)
		}
	})

	t.Run("seek_past_end", func(t *testing.T) {
		sr := NewSeekReader(ra, index)
		if err := sr.SeekTo(int64(len(payload)) + 100); err != ErrSeekPastEnd {
			t.Fatalf("seek past end: want ErrSeekPastEnd, got %v", err)
		}
	})

	t.Run("offset_tracking", func(t *testing.T) {
		sr := NewSeekReader(ra, index)
		if err := sr.SeekTo(1234); err != nil {
			t.Fatalf("seek: %v", err)
		}
		if sr.Offset() != 1234 {
			t.Fatalf("Offset after seek = %d, want 1234", sr.Offset())
		}
		b := make([]byte, 10)
		if _, err := sr.Read(b); err != nil {
			t.Fatalf("read: %v", err)
		}
		if sr.Offset() != 1244 {
			t.Fatalf("Offset after read = %d, want 1244", sr.Offset())
		}
	})
}
