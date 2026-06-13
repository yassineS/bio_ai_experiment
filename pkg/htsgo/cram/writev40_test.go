package cram

import (
	"bytes"
	"fmt"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// CRAM v4.0 ENCODE tests. The writer gained a v4.0 target (VersionV40):
// uint7 LEB128 varint framing throughout the container/block/slice headers
// and compression-header maps, 64-bit alignment coordinates, the
// VARINT_UNSIGNED / VARINT_SIGNED integer data-series codecs, and the
// distinct 31-byte v4 EOF marker. These tests prove the v4 output both
// re-reads through our own decoder (the writer<->reader contract) and is
// genuine CRAM 4.0 that upstream `samtools view` reads identically to our
// already-proven v3.0 output.

// TestWriteCRAMV40RoundTrip writes the broad C8 record corpus as CRAM v4.0,
// reads it back through our own RecordReader, and asserts every record
// survives. The unmapped-read MAPQ/CIGAR caveat applies, exactly as for the
// v3.0 / v3.1 round-trip tests. This is the primary writer<->reader
// contract test for v4.
func TestWriteCRAMV40RoundTrip(t *testing.T) {
	h := writerTestHeader()
	records := v31SampleRecords()
	out := roundTripVersion(t, h, records, VersionV40)
	if len(out) != len(records) {
		t.Fatalf("round-trip record count = %d, want %d", len(out), len(records))
	}
	for i := range records {
		assertRecordEqual(t, i, out[i], records[i])
	}
}

// TestWriteCRAMV40CoreShapes round-trips the focused single-record shapes
// (mapped, unmapped, reverse-strand, every CIGAR op, aux tags, mate pairs,
// no-quality, multi-reference) through v4.0, mirroring the v3.0 writer
// tests so a v4-specific framing bug in any single shape is caught in
// isolation.
func TestWriteCRAMV40CoreShapes(t *testing.T) {
	h := writerTestHeader()

	cig, _ := sam.ParseCigar("10M")
	mapped := &sam.Record{
		QName: "read1", Flag: 0, RName: "chr1", Pos: 100, MapQ: 60,
		Cigar: cig, RNext: "*", Seq: "ACGTACGTAC",
		Qual: []byte{30, 31, 32, 33, 34, 35, 36, 37, 38, 39},
	}
	unmapped := &sam.Record{
		QName: "unmapped1", Flag: sam.FlagUnmapped, RName: "*", RNext: "*",
		Seq: "TTTTGGGGCC", Qual: []byte{20, 21, 22, 23, 24, 25, 26, 27, 28, 29},
	}
	rev := mkRec("rev", "chr2", 1000, "12M", "ACGTACGTACGT")
	rev.Flag = sam.FlagReverse

	cigarShapes := []*sam.Record{
		mkRec("ins", "chr1", 200, "4M2I4M", "AAAACCGGGG"),
		mkRec("del", "chr1", 300, "5M3D5M", "ACGTAACGTA"),
		mkRec("soft", "chr1", 400, "3S7M", "TTTACGTACG"),
		mkRec("hard", "chr1", 500, "2H8M", "ACGTACGT"),
		mkRec("skip", "chr1", 600, "4M10N4M", "ACGTACGT"),
		mkRec("pad", "chr1", 700, "4M2P4M", "ACGTACGT"),
		mkRec("mix", "chr1", 800, "2S3M2D3M2H", "TTACGACG"),
	}

	records := append([]*sam.Record{mapped, unmapped, rev}, cigarShapes...)
	out := roundTripVersion(t, h, records, VersionV40)
	if len(out) != len(records) {
		t.Fatalf("got %d records, want %d", len(out), len(records))
	}
	for i := range records {
		assertRecordEqual(t, i, out[i], records[i])
	}
}

// TestWriteCRAMV40FileDefinition confirms a v4.0-written file carries CRAM
// major 4 / minor 0 and that every block decompresses to its declared
// size. (Block compression for v4 is the v3.0 set; the v4 change is the
// integer framing and the data-series codecs, not the block codec.)
func TestWriteCRAMV40FileDefinition(t *testing.T) {
	h := writerTestHeader()
	var buf bytes.Buffer
	if err := WriteCRAMVersion(&buf, h, v31SampleRecords(), VersionV40); err != nil {
		t.Fatalf("WriteCRAMVersion: %v", err)
	}
	fd, _ := blockMethodStats(t, buf.Bytes())
	if fd.Major != 4 || fd.Minor != 0 {
		t.Errorf("file definition = CRAM %d.%d, want 4.0", fd.Major, fd.Minor)
	}
}

// TestWriteCRAMV40UsesV4Codecs confirms the v4 writer genuinely declares
// the v4-only VARINT codecs (VARINT_UNSIGNED and VARINT_SIGNED) in its
// data-series encoding map, so the round-trip above exercises the new
// integer-codec paths rather than passing on EXTERNAL-of-ITF8 framing that
// would be wrong for v4. It also confirms the v3.0 writer never emits them.
func TestWriteCRAMV40UsesV4Codecs(t *testing.T) {
	for _, tc := range []struct {
		version Version
		want    bool
	}{
		{VersionV40, true},
		{VersionV30, false},
	} {
		var buf bytes.Buffer
		if err := WriteCRAMVersion(&buf, writerTestHeader(), v31SampleRecords(), tc.version); err != nil {
			t.Fatalf("WriteCRAMVersion(%s): %v", tc.version, err)
		}
		rd, err := NewReader(bytes.NewReader(buf.Bytes()))
		if err != nil {
			t.Fatalf("NewReader: %v", err)
		}
		seen := map[EncodingID]bool{}
		for {
			c, err := rd.Next()
			if err != nil {
				break
			}
			dc, derr := ParseDataContainer(c)
			if derr != nil {
				continue
			}
			for _, enc := range dc.Compression.DataSeries {
				seen[enc.ID] = true
			}
		}
		hasVarint := seen[EncodingVarintUnsigned] && seen[EncodingVarintSigned]
		if hasVarint != tc.want {
			t.Errorf("%s: VARINT codecs present = %v, want %v (seen=%v)",
				tc.version, hasVarint, tc.want, seen)
		}
	}
}

// TestWriteCRAMV40Empty writes a v4.0 file with no records and confirms it
// reads back as an empty stream with the header intact — exercising the v4
// file-definition, uint7 header-container framing and v4 EOF recognition on
// the empty path.
func TestWriteCRAMV40Empty(t *testing.T) {
	h := writerTestHeader()
	var buf bytes.Buffer
	if err := WriteCRAMVersion(&buf, h, nil, VersionV40); err != nil {
		t.Fatalf("WriteCRAMVersion: %v", err)
	}
	rr, err := NewRecordReader(bytes.NewReader(buf.Bytes()))
	if err != nil {
		t.Fatalf("NewRecordReader: %v", err)
	}
	if len(rr.Header().Refs) != 2 {
		t.Fatalf("header refs = %d, want 2", len(rr.Header().Refs))
	}
	out, err := rr.ReadAll()
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if len(out) != 0 {
		t.Fatalf("got %d records from an empty v4 file, want 0", len(out))
	}
}

// TestWriteCRAMV40EOFMarker confirms the v4 writer terminates the file with
// the distinct 31-byte v4 EOF sentinel, not the 38-byte v3 one, so a
// v4-aware reader recognises the clean end of stream.
func TestWriteCRAMV40EOFMarker(t *testing.T) {
	var buf bytes.Buffer
	if err := WriteCRAMVersion(&buf, writerTestHeader(), nil, VersionV40); err != nil {
		t.Fatalf("WriteCRAMVersion: %v", err)
	}
	data := buf.Bytes()
	if !bytes.HasSuffix(data, eofMarkerV4) {
		t.Fatalf("v4 file does not end with the 31-byte v4 EOF marker")
	}
	if bytes.HasSuffix(data, eofMarkerV3) {
		t.Fatalf("v4 file ended with the v3 EOF marker")
	}
}

// TestCreateCRAMVersionV40 exercises the file-path v4.0 constructor and
// confirms the file it produces round-trips and is a v4.0 file.
func TestCreateCRAMVersionV40(t *testing.T) {
	h := writerTestHeader()
	path := filepath.Join(t.TempDir(), "out40.cram")
	rw, err := CreateCRAMVersion(path, h, VersionV40)
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
	out, err := rr.ReadAll()
	rr.Close()
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
	if fd := rd.FileDefinition(); fd.Major != 4 || fd.Minor != 0 {
		t.Errorf("file definition = CRAM %d.%d, want 4.0", fd.Major, fd.Minor)
	}
}

// TestWriterVersionStringV40 confirms Version.String reports "4.0" for the
// new target.
func TestWriterVersionStringV40(t *testing.T) {
	if got := VersionV40.String(); got != "4.0" {
		t.Errorf("VersionV40.String() = %q, want \"4.0\"", got)
	}
}

// v40CrossCheckRecords is a corpus broad enough that the live samtools
// cross-check exercises every data series the v4 writer emits — mapped
// reads with every CIGAR op, unmapped reads, aux tags of each type, and
// mate pairs — while staying within the simple writer's encodable shapes.
func v40CrossCheckRecords() []*sam.Record {
	// Every record gets a distinct QNAME so the cross-check can key each
	// decoded SAM line back to its original by QNAME+FLAG.
	var recs []*sam.Record
	for i := 0; i < 20; i++ {
		recs = append(recs, mkRec(fmt.Sprintf("read%d", i), "chr1", int32(100+i*5), "8M", "ACGTACGT"))
	}
	recs = append(recs,
		mkRec("ins", "chr1", 200, "4M2I4M", "AAAACCGGGG"),
		mkRec("del", "chr1", 300, "5M3D5M", "ACGTAACGTA"),
		mkRec("soft", "chr1", 400, "3S7M", "TTTACGTACG"),
		mkRec("skip", "chr1", 600, "4M10N4M", "ACGTACGT"),
	)
	for i := 0; i < 4; i++ {
		recs = append(recs, &sam.Record{
			QName: fmt.Sprintf("unmapped%d", i), Flag: sam.FlagUnmapped, RName: "*", RNext: "*",
			Seq: "TTTTGGGGCC", Qual: []byte{20, 21, 22, 23, 24, 25, 26, 27, 28, 29},
		})
	}
	aux1 := mkRec("aux1", "chr1", 900, "8M", "ACGTACGT")
	aux1.Aux = []sam.Aux{
		{Tag: "NM", Type: 'i', Value: int64(3)},
		{Tag: "Xc", Type: 'c', Value: int64(-5)},
		{Tag: "XA", Type: 'A', Value: "Q"},
		{Tag: "MD", Type: 'Z', Value: "8"},
	}
	recs = append(recs, aux1)
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

// TestWriteCRAMV40SamtoolsCrossCheck is the live oracle: it writes a record
// corpus as CRAM v4.0 with our Go writer, then has the vendored upstream
// `samtools view` decode it to SAM text and asserts the decoded records
// match the originals field-by-field (QNAME/FLAG/RNAME/POS/MAPQ/CIGAR/
// RNEXT/PNEXT/TLEN/SEQ/QUAL), normalising the documented CRAM lossiness the
// v3 writer tests already cover (=/X→M — not exercised here, the corpus
// uses M; unmapped reads lose MAPQ and CIGAR; absent quality round-trips to
// absent). A successful samtools decode proves the v4.0 bytes are genuine
// CRAM 4.0 — real uint7 framing, real VARINT/EXTERNAL codecs and the v4
// CORE-block/slice layout samtools requires — not merely self-consistent
// with our own reader. (Upstream prints a "v4.0 is still a draft" warning on
// stderr but reads the file.) A missing or un-buildable upstream is a hard
// failure, never a skip.
//
// The comparison goes against the originals rather than against our own v3
// output because the simple writer's v3 path, while it round-trips through
// our reader, is not itself decodable by upstream samtools (it omits the
// CORE block htslib requires) — so the v4 path is the first writer output
// this tree proves against the live upstream decoder.
func TestWriteCRAMV40SamtoolsCrossCheck(t *testing.T) {
	samtools := upstreamSamtoolsCram(t)
	h := writerTestHeader()
	records := v40CrossCheckRecords()
	dir := t.TempDir()

	v4Path := filepath.Join(dir, "ours.v4.cram")
	if err := writeRecordsToFile(v4Path, h, records, VersionV40); err != nil {
		t.Fatalf("writing v4.0 CRAM: %v", err)
	}

	// Confirm samtools genuinely sees our file as CRAM major version 4.
	if maj := openMajor(t, v4Path); maj != 4 {
		t.Fatalf("our v4 file reports CRAM major %d, want 4", maj)
	}

	// samtools decodes our v4 file to SAM text; key each decoded line by
	// QNAME+FLAG so the comparison is order-independent.
	got := map[string][]string{}
	for _, line := range samtoolsView(t, samtools, v4Path) {
		f := strings.Split(line, "\t")
		if len(f) < 11 {
			t.Fatalf("samtools emitted a short SAM line: %q", line)
		}
		got[f[0]+"/"+f[1]] = f
	}
	if len(got) != len(records) {
		t.Fatalf("samtools decoded %d records from our v4 file, want %d", len(got), len(records))
	}

	for i, rec := range records {
		key := rec.QName + "/" + flagString(rec.Flag)
		f, ok := got[key]
		if !ok {
			t.Fatalf("record %d (%s) absent from samtools' v4 decode", i, key)
		}
		assertSAMFieldsMatch(t, i, f, rec)
	}
}

// assertSAMFieldsMatch compares one samtools-emitted SAM field slice f
// against the original record rec, applying the documented CRAM lossiness:
// an unmapped read loses its MAPQ and CIGAR, and an absent quality decodes
// back to absent ("*").
func assertSAMFieldsMatch(t *testing.T, idx int, f []string, rec *sam.Record) {
	t.Helper()
	// f: QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL ...
	if f[1] != flagString(rec.Flag) {
		t.Errorf("record %d FLAG = %s, want %d", idx, f[1], rec.Flag)
	}
	if normName(f[2]) != normName(rec.RName) {
		t.Errorf("record %d RNAME = %q, want %q", idx, f[2], rec.RName)
	}
	if f[3] != intString(int(rec.Pos)) {
		t.Errorf("record %d POS = %s, want %d", idx, f[3], rec.Pos)
	}
	if normName(f[6]) != normName(rec.RNext) {
		t.Errorf("record %d RNEXT = %q, want %q", idx, f[6], rec.RNext)
	}
	if f[7] != intString(int(rec.PNext)) {
		t.Errorf("record %d PNEXT = %s, want %d", idx, f[7], rec.PNext)
	}
	if f[8] != intString(int(rec.TLen)) {
		t.Errorf("record %d TLEN = %s, want %d", idx, f[8], rec.TLen)
	}
	if normSeq(f[9]) != normSeq(rec.Seq) {
		t.Errorf("record %d SEQ = %q, want %q", idx, f[9], rec.Seq)
	}
	// Quality: samtools prints "*" for absent; our records carry explicit
	// ascending quality, which samtools renders as printable ASCII (Phred+33).
	wantQual := phred33(rec.Qual, len(normSeq(rec.Seq)))
	if f[10] != wantQual {
		t.Errorf("record %d QUAL = %q, want %q", idx, f[10], wantQual)
	}
	if rec.Flag&sam.FlagUnmapped == 0 {
		if f[4] != intString(int(rec.MapQ)) {
			t.Errorf("record %d MAPQ = %s, want %d", idx, f[4], rec.MapQ)
		}
		if f[5] != rec.Cigar.String() {
			t.Errorf("record %d CIGAR = %q, want %q", idx, f[5], rec.Cigar.String())
		}
	}
}

// phred33 renders a quality slice as the SAM Phred+33 ASCII string, or "*"
// when the quality is absent (empty / all-0xff). readLen is the SEQ length.
func phred33(qual []byte, readLen int) string {
	if qualIsAbsent(qual) || readLen == 0 {
		return "*"
	}
	out := make([]byte, len(qual))
	for i, q := range qual {
		out[i] = q + 33
	}
	return string(out)
}

// flagString formats a SAM flag as its decimal string.
func flagString(flag uint16) string { return intString(int(flag)) }

// intString formats an int as its decimal string.
func intString(v int) string { return strconv.Itoa(v) }

// writeRecordsToFile writes h and records to a CRAM file at path with the
// given version.
func writeRecordsToFile(path string, h *sam.Header, records []*sam.Record, version Version) error {
	rw, err := CreateCRAMVersion(path, h, version)
	if err != nil {
		return err
	}
	for _, rec := range records {
		if err := rw.Write(rec); err != nil {
			rw.Close()
			return err
		}
	}
	return rw.Close()
}

// openMajor returns the CRAM major version of the file at path.
func openMajor(t *testing.T, path string) uint8 {
	t.Helper()
	rd, err := Open(path)
	if err != nil {
		t.Fatalf("Open %s: %v", path, err)
	}
	defer rd.Close()
	return rd.FileDefinition().Major
}

// samtoolsView runs `samtools view path` and returns the alignment record
// lines (no header).
func samtoolsView(t *testing.T, samtools, path string) []string {
	t.Helper()
	out, err := exec.Command(samtools, "view", path).Output()
	if err != nil {
		// Surface stderr (the "v4.0 draft" warning, or a real decode error).
		if ee, ok := err.(*exec.ExitError); ok {
			t.Fatalf("samtools view %s: %v\n%s", path, err, ee.Stderr)
		}
		t.Fatalf("samtools view %s: %v", path, err)
	}
	lines := splitNonEmptyLines(string(out))
	cleaned := make([]string, 0, len(lines))
	for _, l := range lines {
		if strings.HasPrefix(l, "@") {
			continue
		}
		cleaned = append(cleaned, l)
	}
	return cleaned
}

// TestWriteCRAMV30SamtoolsCrossCheck proves the v3.0 writer output is now
// decodable by upstream samtools. Before the per-slice CORE block was added,
// `samtools view` failed our v3 output with "Failure to decode slice"; this
// pins the fix using the same field-by-field comparison as the v4 cross-check.
func TestWriteCRAMV30SamtoolsCrossCheck(t *testing.T) {
	samtools := upstreamSamtoolsCram(t)
	h := writerTestHeader()
	records := v40CrossCheckRecords()
	dir := t.TempDir()

	v3Path := filepath.Join(dir, "ours.v3.cram")
	if err := writeRecordsToFile(v3Path, h, records, VersionV30); err != nil {
		t.Fatalf("writing v3.0 CRAM: %v", err)
	}
	if maj := openMajor(t, v3Path); maj != 3 {
		t.Fatalf("our v3 file reports CRAM major %d, want 3", maj)
	}

	got := map[string][]string{}
	for _, line := range samtoolsView(t, samtools, v3Path) {
		f := strings.Split(line, "\t")
		if len(f) < 11 {
			t.Fatalf("samtools emitted a short SAM line: %q", line)
		}
		got[f[0]+"/"+f[1]] = f
	}
	if len(got) != len(records) {
		t.Fatalf("samtools decoded %d records from our v3 file, want %d", len(got), len(records))
	}
	for i, rec := range records {
		f, ok := got[rec.QName+"/"+flagString(rec.Flag)]
		if !ok {
			t.Fatalf("record %d (%s) absent from samtools' v3 decode", i, rec.QName)
		}
		assertSAMFieldsMatch(t, i, f, rec)
	}
}
