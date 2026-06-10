package fasta

import (
	"bytes"
	"encoding/binary"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"sync"
	"testing"

	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
)

// makeBGZF returns the BGZF-compressed bytes of payload.
func makeBGZF(t *testing.T, payload []byte) []byte {
	t.Helper()
	var buf bytes.Buffer
	w := bgzip.NewWriter(&buf)
	if _, err := w.Write(payload); err != nil {
		t.Fatalf("bgzip write: %v", err)
	}
	if err := w.Close(); err != nil {
		t.Fatalf("bgzip close: %v", err)
	}
	return buf.Bytes()
}

// writeTempFile writes data to a fresh temp file under t.TempDir() and
// returns the path. The file is cleaned up automatically when the test
// exits.
func writeTempFile(t *testing.T, dir, name string, data []byte) string {
	t.Helper()
	p := filepath.Join(dir, name)
	if err := os.WriteFile(p, data, 0o644); err != nil {
		t.Fatalf("write %s: %v", name, err)
	}
	return p
}

func TestIsBGZFFile(t *testing.T) {
	dir := t.TempDir()

	bgzfData := makeBGZF(t, []byte(">chr1\nACGT\n"))
	bgzfPath := writeTempFile(t, dir, "ref.fa.gz", bgzfData)
	got, err := isBGZFFile(bgzfPath)
	if err != nil {
		t.Fatalf("isBGZFFile(bgzf): %v", err)
	}
	if !got {
		t.Errorf("isBGZFFile on BGZF file: got false, want true")
	}

	plainPath := writeTempFile(t, dir, "plain.fa", []byte(">chr1\nACGT\n"))
	got, err = isBGZFFile(plainPath)
	if err != nil {
		t.Fatalf("isBGZFFile(plain): %v", err)
	}
	if got {
		t.Errorf("isBGZFFile on plain text: got true, want false")
	}

	// Tiny file (<4 bytes) → not BGZF.
	tinyPath := writeTempFile(t, dir, "tiny.bin", []byte{0x1f, 0x8b})
	got, err = isBGZFFile(tinyPath)
	if err != nil {
		t.Fatalf("isBGZFFile(tiny): %v", err)
	}
	if got {
		t.Errorf("isBGZFFile on truncated header: got true, want false")
	}
}

func TestLoadGZIBlockOffsets_RoundTrip(t *testing.T) {
	// Hand-build a .gzi with two entries beyond the implicit (0,0).
	var buf bytes.Buffer
	var hdr [8]byte
	binary.LittleEndian.PutUint64(hdr[:], 2)
	buf.Write(hdr[:])
	var entry [16]byte
	binary.LittleEndian.PutUint64(entry[0:8], 100)
	binary.LittleEndian.PutUint64(entry[8:16], 65280)
	buf.Write(entry[:])
	binary.LittleEndian.PutUint64(entry[0:8], 240)
	binary.LittleEndian.PutUint64(entry[8:16], 130560)
	buf.Write(entry[:])

	entries, err := readGZIBlockOffsets(&buf)
	if err != nil {
		t.Fatalf("readGZIBlockOffsets: %v", err)
	}
	if got, want := len(entries), 2; got != want {
		t.Fatalf("entries: got %d, want %d", got, want)
	}
	if entries[0].CompressedOffset != 100 || entries[0].UncompressedOffset != 65280 {
		t.Errorf("entry[0] = %+v, want {100, 65280}", entries[0])
	}
	if entries[1].CompressedOffset != 240 || entries[1].UncompressedOffset != 130560 {
		t.Errorf("entry[1] = %+v, want {240, 130560}", entries[1])
	}
}

func TestLoadGZIBlockOffsets_ShortFile(t *testing.T) {
	if _, err := readGZIBlockOffsets(bytes.NewReader([]byte{0x01, 0x02})); err == nil {
		t.Fatal("expected error on short .gzi header, got nil")
	}
	// Header says 2 entries, but only one is provided.
	var buf bytes.Buffer
	var hdr [8]byte
	binary.LittleEndian.PutUint64(hdr[:], 2)
	buf.Write(hdr[:])
	var entry [16]byte
	binary.LittleEndian.PutUint64(entry[0:8], 100)
	binary.LittleEndian.PutUint64(entry[8:16], 65280)
	buf.Write(entry[:])
	if _, err := readGZIBlockOffsets(&buf); err == nil {
		t.Fatal("expected error on missing entry, got nil")
	}
}

func TestDecompressBGZF_RoundTrip(t *testing.T) {
	dir := t.TempDir()
	payload := []byte(">chr1\nACGT\nACGT\n>chr2\nNNNN\n")
	p := writeTempFile(t, dir, "ref.fa.gz", makeBGZF(t, payload))
	got, err := decompressBGZF(p)
	if err != nil {
		t.Fatalf("decompressBGZF: %v", err)
	}
	if !bytes.Equal(got, payload) {
		t.Errorf("payload mismatch: got %q want %q", got, payload)
	}
}

func TestOpenRandomAccessBGZF_Fetch(t *testing.T) {
	dir := t.TempDir()
	payload := []byte(">chr1\nACGTACGTAC\nGTACGTACGT\n>chr2\nNNNNNNNNNN\n")
	p := writeTempFile(t, dir, "ref.fa.gz", makeBGZF(t, payload))

	ra, err := OpenRandomAccessBGZF(p)
	if err != nil {
		t.Fatalf("OpenRandomAccessBGZF: %v", err)
	}
	defer ra.Close()

	cases := []struct {
		name       string
		start, end int64
		want       string
	}{
		{"chr1", 0, 4, "ACGT"},
		{"chr1", 8, 12, "ACGT"}, // spans the line break in the file
		{"chr1", 0, 20, "ACGTACGTACGTACGTACGT"},
		{"chr2", 0, 10, "NNNNNNNNNN"},
	}
	for _, c := range cases {
		got, err := ra.Fetch(c.name, c.start, c.end)
		if err != nil {
			t.Errorf("Fetch(%s, %d, %d): %v", c.name, c.start, c.end, err)
			continue
		}
		if string(got) != c.want {
			t.Errorf("Fetch(%s, %d, %d) = %q, want %q", c.name, c.start, c.end, got, c.want)
		}
	}

	// Through-the-roof OpenRandomAccess should auto-route to the BGZF
	// implementation by sniffing the magic.
	ra2, err := OpenRandomAccess(p)
	if err != nil {
		t.Fatalf("OpenRandomAccess (auto BGZF): %v", err)
	}
	defer ra2.Close()
	seq, err := ra2.Fetch("chr1", 0, 4)
	if err != nil {
		t.Fatalf("OpenRandomAccess Fetch: %v", err)
	}
	if string(seq) != "ACGT" {
		t.Errorf("auto-BGZF Fetch = %q, want %q", seq, "ACGT")
	}
}

func TestOpenRandomAccessBGZF_RejectsPlainFASTA(t *testing.T) {
	dir := t.TempDir()
	p := writeTempFile(t, dir, "plain.fa", []byte(">chr1\nACGT\n"))
	if _, err := OpenRandomAccessBGZF(p); err == nil {
		t.Fatal("expected error opening plain FASTA via BGZF API, got nil")
	}
}

func TestOpenRandomAccessBGZFFullHeader(t *testing.T) {
	dir := t.TempDir()
	payload := []byte(">chr1 assembled by consortium X\nACGT\n")
	p := writeTempFile(t, dir, "ref.fa.gz", makeBGZF(t, payload))
	ra, err := OpenRandomAccessBGZFFullHeader(p)
	if err != nil {
		t.Fatalf("OpenRandomAccessBGZFFullHeader: %v", err)
	}
	defer ra.Close()
	got, err := ra.Fetch("chr1 assembled by consortium X", 0, 4)
	if err != nil {
		t.Fatalf("Fetch full-header contig: %v", err)
	}
	if string(got) != "ACGT" {
		t.Errorf("Fetch = %q, want %q", got, "ACGT")
	}
}

func TestOpenRandomAccessBGZF_BadGZI(t *testing.T) {
	dir := t.TempDir()
	payload := []byte(">chr1\nACGT\n")
	p := writeTempFile(t, dir, "ref.fa.gz", makeBGZF(t, payload))
	// Drop a broken .gzi next to it.
	if err := os.WriteFile(p+".gzi", []byte{0x01, 0x02}, 0o644); err != nil {
		t.Fatalf("write bad .gzi: %v", err)
	}
	if _, err := OpenRandomAccessBGZF(p); err == nil {
		t.Fatal("expected error on bad .gzi, got nil")
	}
}

// Sanity-check the BuildIndex{,FullHeader}Bytes exports.
func TestBuildIndexBytes(t *testing.T) {
	payload := []byte(">chr1\nACGT\n>chr2\nNNNN\n")
	idx, err := BuildIndexBytes(payload)
	if err != nil {
		t.Fatalf("BuildIndexBytes: %v", err)
	}
	if got := len(idx.Entries()); got != 2 {
		t.Errorf("entries: got %d, want 2", got)
	}
	idxF, err := BuildIndexFullHeaderBytes([]byte(">chr1 X\nACGT\n"))
	if err != nil {
		t.Fatalf("BuildIndexFullHeaderBytes: %v", err)
	}
	if _, ok := idxF.Get("chr1 X"); !ok {
		t.Errorf("full-header index missing %q", "chr1 X")
	}
}

// upstreamBgzipForFasta locates/builds the htslib bgzip binary used to produce
// real .gz + .gzi sidecars for the partial-decompression Fetch parity test.
var (
	fastaBgzipOnce sync.Once
	fastaBgzipPath string
	fastaBgzipErr  error
)

func upstreamBgzipForFasta(t *testing.T) string {
	t.Helper()
	fastaBgzipOnce.Do(func() {
		htslibDir, err := filepath.Abs("../../../reference_code/htslib")
		if err != nil {
			fastaBgzipErr = err
			return
		}
		bin := filepath.Join(htslibDir, "bgzip")
		if _, statErr := os.Stat(bin); statErr == nil {
			fastaBgzipPath = bin
			return
		}
		if _, statErr := os.Stat(filepath.Join(htslibDir, "config.mk")); statErr != nil {
			for _, args := range [][]string{{"autoreconf", "-i"}, {"./configure"}} {
				cmd := exec.Command(args[0], args[1:]...)
				cmd.Dir = htslibDir
				if out, runErr := cmd.CombinedOutput(); runErr != nil {
					fastaBgzipErr = fmt.Errorf("%v: %v\n%s", args, runErr, out)
					return
				}
			}
		}
		cmd := exec.Command("make", "-j4", "bgzip")
		cmd.Dir = htslibDir
		if out, runErr := cmd.CombinedOutput(); runErr != nil {
			fastaBgzipErr = fmt.Errorf("make bgzip: %v\n%s", runErr, out)
			return
		}
		fastaBgzipPath = bin
	})
	if fastaBgzipErr != nil {
		t.Fatalf("locating/building upstream bgzip: %v", fastaBgzipErr)
	}
	if fastaBgzipPath == "" {
		t.Fatalf("upstream bgzip not available")
	}
	return fastaBgzipPath
}

// TestBGZFFetch_GZIPartialSeekParity confirms that OpenRandomAccessBGZF, when a
// .gzi sidecar (produced by upstream bgzip) is present, serves Fetch via the
// SeekReader partial-decompression path and returns bytes identical to a plain
// uncompressed FASTA Fetch over the same data.
func TestBGZFFetch_GZIPartialSeekParity(t *testing.T) {
	bgzip := upstreamBgzipForFasta(t)
	dir := t.TempDir()

	// Build a multi-contig FASTA whose total size spans several BGZF blocks.
	var fa bytes.Buffer
	contigs := map[string]string{}
	bases := []byte("ACGTN")
	for c := 0; c < 4; c++ {
		name := fmt.Sprintf("chr%d", c+1)
		fmt.Fprintf(&fa, ">%s\n", name)
		var seq bytes.Buffer
		n := 20000 + c*7000
		for i := 0; i < n; i++ {
			seq.WriteByte(bases[(i*7+c)%len(bases)])
		}
		contigs[name] = seq.String()
		// Wrap at 60 columns like samtools faidx expects.
		s := seq.Bytes()
		for i := 0; i < len(s); i += 60 {
			end := i + 60
			if end > len(s) {
				end = len(s)
			}
			fa.Write(s[i:end])
			fa.WriteByte('\n')
		}
	}

	plainPath := filepath.Join(dir, "ref.fa")
	if err := os.WriteFile(plainPath, fa.Bytes(), 0o644); err != nil {
		t.Fatalf("write plain fasta: %v", err)
	}

	// Compress with upstream bgzip -i (writes ref.fa.gz + ref.fa.gz.gzi).
	gzPath := plainPath + ".gz"
	if out, err := exec.Command(bgzip, "-i", "-k", plainPath).CombinedOutput(); err != nil {
		t.Fatalf("bgzip -i: %v\n%s", err, out)
	}
	if _, err := os.Stat(gzPath + ".gzi"); err != nil {
		t.Fatalf("expected .gzi sidecar: %v", err)
	}

	// Build the samtools-style .fai (offsets into the uncompressed stream) by
	// indexing the plain FASTA, and copy it to the .gz.fai location.
	plainIdx, err := BuildIndex(plainPath)
	if err != nil {
		t.Fatalf("BuildIndex: %v", err)
	}
	if err := plainIdx.Save(gzPath + ".fai"); err != nil {
		t.Fatalf("Save index: %v", err)
	}

	// Open both the compressed (.gzi-backed) and plain references.
	bgzfRA, err := OpenRandomAccessBGZF(gzPath)
	if err != nil {
		t.Fatalf("OpenRandomAccessBGZF: %v", err)
	}
	defer bgzfRA.Close()
	// Confirm the partial-seek backend is in use (not the in-memory fallback).
	if _, ok := bgzfRA.r.(*gziReaderAt); !ok {
		t.Fatalf("expected gziReaderAt backend, got %T", bgzfRA.r)
	}

	plainRA, err := OpenRandomAccess(plainPath)
	if err != nil {
		t.Fatalf("OpenRandomAccess plain: %v", err)
	}
	defer plainRA.Close()

	type reg struct {
		name       string
		start, end int64
	}
	regions := []reg{
		{"chr1", 0, 10},
		{"chr1", 100, 200},
		{"chr2", 0, 27000},
		{"chr3", 12345, 12400},
		{"chr4", 40000 - 5, 40000},
		{"chr1", 19990, 20000},
	}
	for _, r := range regions {
		got, err := bgzfRA.Fetch(r.name, r.start, r.end)
		if err != nil {
			t.Fatalf("bgzf Fetch %s:%d-%d: %v", r.name, r.start, r.end, err)
		}
		want, err := plainRA.Fetch(r.name, r.start, r.end)
		if err != nil {
			t.Fatalf("plain Fetch %s:%d-%d: %v", r.name, r.start, r.end, err)
		}
		if !bytes.Equal(got, want) {
			t.Fatalf("%s:%d-%d mismatch:\n got=%s\nwant=%s", r.name, r.start, r.end, got, want)
		}
	}
}
