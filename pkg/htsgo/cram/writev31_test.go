package cram

import (
	"bytes"
	"errors"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// roundTripVersion writes records through a RecordWriter targeting the
// given CRAM version and reads them back, returning the decoded records.
// It is the version-aware sibling of roundTrip.
func roundTripVersion(t *testing.T, h *sam.Header, records []*sam.Record, version Version) []*sam.Record {
	t.Helper()
	var buf bytes.Buffer
	if err := WriteCRAMVersion(&buf, h, records, version); err != nil {
		t.Fatalf("WriteCRAMVersion(%s): %v", version, err)
	}
	rr, err := NewRecordReader(bytes.NewReader(buf.Bytes()))
	if err != nil {
		t.Fatalf("NewRecordReader: %v", err)
	}
	out, err := rr.ReadAll()
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	return out
}

// blockMethodStats walks every container of a CRAM byte stream and
// reports the file definition together with the set of block compression
// methods it observed. It is a writer-side oracle: it confirms a v3.1
// file genuinely reaches the rANS 4x16 path and a v3.0 file never does.
func blockMethodStats(t *testing.T, data []byte) (FileDefinition, map[CompressionMethod]int) {
	t.Helper()
	rd, err := NewReader(bytes.NewReader(data))
	if err != nil {
		t.Fatalf("NewReader: %v", err)
	}
	methods := make(map[CompressionMethod]int)
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
			methods[b.Method]++
			out, err := b.Decompress()
			if err != nil {
				t.Fatalf("block (method %d, content id %d): %v", b.Method, b.ContentID, err)
			}
			if int32(len(out)) != b.UncompressedSize {
				t.Fatalf("block (method %d) decompressed to %d bytes, declared %d",
					b.Method, len(out), b.UncompressedSize)
			}
		}
	}
	return rd.FileDefinition(), methods
}

// v31SampleRecords returns a batch of records broad enough that every
// data series the writer can populate is exercised, and large enough
// (the QS, RN and BF blocks clear the 32-byte codec threshold) that the
// v3.1 writer genuinely picks rANS 4x16 for at least one block. It mixes
// mapped reads with every CIGAR operation, unmapped reads, aux tags
// including B-arrays, and mate pairs across references.
func v31SampleRecords() []*sam.Record {
	var recs []*sam.Record
	// Many plain mapped reads: the repetitive QS and RN series compress
	// well under rANS 4x16, forcing the v3.1 codec onto the write path.
	for i := 0; i < 60; i++ {
		recs = append(recs, mkRec("read", "chr1", int32(100+i), "8M", "ACGTACGT"))
	}
	// Every CIGAR shape.
	recs = append(recs,
		mkRec("ins", "chr1", 200, "4M2I4M", "AAAACCGGGG"),
		mkRec("del", "chr1", 300, "5M3D5M", "ACGTAACGTA"),
		mkRec("soft", "chr1", 400, "3S7M", "TTTACGTACG"),
		mkRec("hard", "chr1", 500, "2H8M", "ACGTACGT"),
		mkRec("skip", "chr1", 600, "4M10N4M", "ACGTACGT"),
		mkRec("pad", "chr1", 700, "4M2P4M", "ACGTACGT"),
		mkRec("mix", "chr1", 800, "2S3M2D3M2H", "TTACGACG"),
	)
	// Unmapped reads.
	for i := 0; i < 5; i++ {
		recs = append(recs, &sam.Record{
			QName: "unmapped", Flag: sam.FlagUnmapped, RName: "*", RNext: "*",
			Seq: "TTTTGGGGCC", Qual: []byte{20, 21, 22, 23, 24, 25, 26, 27, 28, 29},
		})
	}
	// Aux tags of every supported type, including B-arrays.
	aux1 := mkRec("aux1", "chr1", 900, "8M", "ACGTACGT")
	aux1.Aux = []sam.Aux{
		{Tag: "NM", Type: 'i', Value: int64(3)},
		{Tag: "Xc", Type: 'c', Value: int64(-5)},
		{Tag: "XF", Type: 'f', Value: float64(2.5)},
		{Tag: "XA", Type: 'A', Value: "Q"},
		{Tag: "MD", Type: 'Z', Value: "8"},
		{Tag: "XH", Type: 'H', Value: "DEADBEEF"},
	}
	aux2 := mkRec("aux2", "chr2", 1000, "8M", "TTTTGGGG")
	aux2.Aux = []sam.Aux{
		{Tag: "BI", Type: 'B', ArrayType: 'i',
			ArrayValues: []interface{}{int64(1), int64(-2), int64(3)}},
		{Tag: "BC", Type: 'B', ArrayType: 'C',
			ArrayValues: []interface{}{int64(10), int64(20), int64(255)}},
	}
	recs = append(recs, aux1, aux2)
	// A mate pair on one reference and a cross-reference pair.
	p1 := mkRec("pair", "chr1", 1100, "10M", "ACGTACGTAC")
	p1.Flag = sam.FlagPaired | sam.FlagProperPair | sam.FlagRead1 | sam.FlagMateReverse
	p1.RNext, p1.PNext, p1.TLen = "=", 1300, 210
	p2 := mkRec("pair", "chr1", 1300, "10M", "TGCATGCATG")
	p2.Flag = sam.FlagPaired | sam.FlagProperPair | sam.FlagRead2 | sam.FlagReverse
	p2.RNext, p2.PNext, p2.TLen = "=", 1100, -210
	x1 := mkRec("xref", "chr1", 1400, "10M", "ACGTACGTAC")
	x1.Flag = sam.FlagPaired | sam.FlagRead1
	x1.RNext, x1.PNext = "chr2", 5000
	recs = append(recs, p1, p2, x1)
	return recs
}

// TestWriteCRAMVersionV31RoundTrip writes the broad C8 record corpus as
// CRAM v3.1, reads it back through the existing reader, and asserts every
// record survives. The unmapped-read MAPQ/CIGAR caveat applies, exactly
// as for the v3.0 round-trip tests.
func TestWriteCRAMVersionV31RoundTrip(t *testing.T) {
	h := writerTestHeader()
	records := v31SampleRecords()
	out := roundTripVersion(t, h, records, VersionV31)
	if len(out) != len(records) {
		t.Fatalf("round-trip record count = %d, want %d", len(out), len(records))
	}
	for i := range records {
		assertRecordEqual(t, i, out[i], records[i])
	}
}

// TestWriteCRAMV31FileDefinition confirms a v3.1-written file carries
// CRAM major 3 / minor 1 and uses at least one rANS 4x16 (method 5)
// block — proving the v3.1 codec path is genuinely exercised — while
// every block still decompresses to its declared size.
func TestWriteCRAMV31FileDefinition(t *testing.T) {
	h := writerTestHeader()
	var buf bytes.Buffer
	if err := WriteCRAMVersion(&buf, h, v31SampleRecords(), VersionV31); err != nil {
		t.Fatalf("WriteCRAMVersion: %v", err)
	}
	fd, methods := blockMethodStats(t, buf.Bytes())
	if fd.Major != 3 || fd.Minor != 1 {
		t.Errorf("file definition = CRAM %d.%d, want 3.1", fd.Major, fd.Minor)
	}
	if methods[CompRANS4x16] == 0 {
		t.Errorf("v3.1 file used no rANS 4x16 block; method counts = %v", methods)
	}
}

// TestWriteCRAMV30NoRANS confirms a v3.0-written file keeps minor version
// 0 and never emits a rANS 4x16 (method 5) block — that codec is a v3.1
// capability and must not leak into v3.0 output.
func TestWriteCRAMV30NoRANS(t *testing.T) {
	h := writerTestHeader()
	var buf bytes.Buffer
	if err := WriteCRAMVersion(&buf, h, v31SampleRecords(), VersionV30); err != nil {
		t.Fatalf("WriteCRAMVersion: %v", err)
	}
	fd, methods := blockMethodStats(t, buf.Bytes())
	if fd.Major != 3 || fd.Minor != 0 {
		t.Errorf("file definition = CRAM %d.%d, want 3.0", fd.Major, fd.Minor)
	}
	if methods[CompRANS4x16] != 0 {
		t.Errorf("v3.0 file used %d rANS 4x16 block(s); method 5 is a v3.1 codec",
			methods[CompRANS4x16])
	}
}

// TestWriteCRAMVersionDefaultsV30 confirms the convenience wrappers that
// do not take a version still produce v3.0, so existing callers are
// unaffected by the addition of v3.1 support.
func TestWriteCRAMVersionDefaultsV30(t *testing.T) {
	h := writerTestHeader()
	var buf bytes.Buffer
	if err := WriteCRAM(&buf, h, v31SampleRecords()); err != nil {
		t.Fatalf("WriteCRAM: %v", err)
	}
	fd, methods := blockMethodStats(t, buf.Bytes())
	if fd.Minor != 0 {
		t.Errorf("WriteCRAM default minor version = %d, want 0", fd.Minor)
	}
	if methods[CompRANS4x16] != 0 {
		t.Errorf("default writer used a method-5 block, want none")
	}
}

// TestWriterVersionString confirms the Version.String helper reports the
// expected "major.minor" form.
func TestWriterVersionString(t *testing.T) {
	if got := VersionV30.String(); got != "3.0" {
		t.Errorf("VersionV30.String() = %q, want \"3.0\"", got)
	}
	if got := VersionV31.String(); got != "3.1" {
		t.Errorf("VersionV31.String() = %q, want \"3.1\"", got)
	}
}

// TestCreateCRAMVersionV31 exercises the file-path v3.1 constructor and
// confirms the file it produces round-trips and is a v3.1 file.
func TestCreateCRAMVersionV31(t *testing.T) {
	h := writerTestHeader()
	path := t.TempDir() + "/out31.cram"
	rw, err := CreateCRAMVersion(path, h, VersionV31)
	if err != nil {
		t.Fatalf("CreateCRAMVersion: %v", err)
	}
	records := v31SampleRecords()
	for i, rec := range records {
		if err := rw.Write(rec); err != nil {
			t.Fatalf("Write record %d: %v", i, err)
		}
	}
	if err := rw.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	rr, err := OpenRecords(path)
	if err != nil {
		t.Fatalf("OpenRecords: %v", err)
	}
	defer rr.Close()
	out, err := rr.ReadAll()
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if len(out) != len(records) {
		t.Fatalf("got %d records, want %d", len(out), len(records))
	}
	for i := range records {
		assertRecordEqual(t, i, out[i], records[i])
	}
	rd, err := Open(path)
	if err != nil {
		t.Fatalf("Open: %v", err)
	}
	defer rd.Close()
	if fd := rd.FileDefinition(); fd.Major != 3 || fd.Minor != 1 {
		t.Errorf("file definition = CRAM %d.%d, want 3.1", fd.Major, fd.Minor)
	}
}

// TestChooseBlockCompression unit-tests the per-version block codec
// selection directly: a tiny payload always stays raw, v3.0 never picks
// method 5, and v3.1 may pick rANS 4x16 on a payload it shrinks while
// never producing output larger than the input.
func TestChooseBlockCompression(t *testing.T) {
	tiny := []byte{1, 2, 3, 4}
	if m, stored := chooseBlockCompression(VersionV31, tiny); m != CompRaw || !bytes.Equal(stored, tiny) {
		t.Errorf("tiny payload: method = %s, want raw stored verbatim", m)
	}
	// A long, highly repetitive payload compresses well under rANS 4x16.
	repetitive := bytes.Repeat([]byte{'A'}, 4096)
	m30, s30 := chooseBlockCompression(VersionV30, repetitive)
	if m30 == CompRANS4x16 {
		t.Errorf("v3.0 selected rANS 4x16, a v3.1-only codec")
	}
	if len(s30) > len(repetitive) {
		t.Errorf("v3.0 stored %d bytes, larger than the %d-byte input", len(s30), len(repetitive))
	}
	m31, s31 := chooseBlockCompression(VersionV31, repetitive)
	if len(s31) > len(repetitive) {
		t.Errorf("v3.1 stored %d bytes, larger than the %d-byte input", len(s31), len(repetitive))
	}
	// rANS 4x16 should win on this input; if a future tweak makes gzip
	// smaller that is still acceptable, but method 5 must be reachable.
	if m31 != CompRANS4x16 && m31 != CompGzip {
		t.Errorf("v3.1 method = %s, want rans4x16 or gzip", m31)
	}
}

// TestWriteCRAMV31SamtoolsCrossCheck is the v3.1 live oracle. It mirrors
// TestWriteCRAMV30SamtoolsCrossCheck (in writev40_test.go) but targets CRAM
// v3.1, whose distinguishing feature is the rANS 4x16 (method 5) codec. It
// writes the shared cross-check corpus with our Go writer, first proves the
// file is genuinely v3.1 and genuinely reaches the rANS 4x16 encode path,
// then has the vendored upstream `samtools view` decode it and asserts every
// record matches field-by-field. The rANS 4x16 codec itself is already pinned
// byte-exact against the htscodecs vectors at the codec layer; this test
// closes the remaining gap by proving a writer-produced v3.1 file — container
// framing, slice layout, and rANS 4x16 blocks together — round-trips through
// the real upstream decoder. A missing or un-buildable upstream is a hard
// failure, never a skip.
func TestWriteCRAMV31SamtoolsCrossCheck(t *testing.T) {
	samtools := upstreamSamtoolsCram(t)
	h := writerTestHeader()
	records := v40CrossCheckRecords()
	dir := t.TempDir()

	v31Path := filepath.Join(dir, "ours.v31.cram")
	if err := writeRecordsToFile(v31Path, h, records, VersionV31); err != nil {
		t.Fatalf("writing v3.1 CRAM: %v", err)
	}

	// Prove the file is genuinely CRAM 3.1 and genuinely exercises the
	// v3.1-only rANS 4x16 codec, so this is a real v3.1 oracle and not a
	// silently-degraded v3.0 path.
	raw, err := os.ReadFile(v31Path)
	if err != nil {
		t.Fatalf("reading back our v3.1 CRAM: %v", err)
	}
	fd, methods := blockMethodStats(t, raw)
	if fd.Major != 3 || fd.Minor != 1 {
		t.Fatalf("our file reports CRAM %d.%d, want 3.1", fd.Major, fd.Minor)
	}
	if methods[CompRANS4x16] == 0 {
		t.Fatalf("our v3.1 file used no rANS 4x16 block; method counts = %v", methods)
	}

	got := map[string][]string{}
	for _, line := range samtoolsView(t, samtools, v31Path) {
		f := strings.Split(line, "\t")
		if len(f) < 11 {
			t.Fatalf("samtools emitted a short SAM line: %q", line)
		}
		got[f[0]+"/"+f[1]] = f
	}
	if len(got) != len(records) {
		t.Fatalf("samtools decoded %d records from our v3.1 file, want %d", len(got), len(records))
	}
	for i, rec := range records {
		f, ok := got[rec.QName+"/"+flagString(rec.Flag)]
		if !ok {
			t.Fatalf("record %d (%s) absent from samtools' v3.1 decode", i, rec.QName)
		}
		assertSAMFieldsMatch(t, i, f, rec)
	}
}
