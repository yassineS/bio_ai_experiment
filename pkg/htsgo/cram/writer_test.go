package cram

import (
	"bytes"
	"reflect"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// writerTestHeader builds a small SAM header with two references and one
// read group, enough to exercise the writer's reference indexing.
func writerTestHeader() *sam.Header {
	text := "@HD\tVN:1.6\tSO:unsorted\n" +
		"@SQ\tSN:chr1\tLN:100000\n" +
		"@SQ\tSN:chr2\tLN:50000\n" +
		"@RG\tID:rg1\tSM:sample1\n"
	h, err := sam.ParseHeaderText(text)
	if err != nil {
		panic(err)
	}
	return h
}

// roundTrip writes records through a RecordWriter and reads them back
// through a RecordReader, returning the decoded records.
func roundTrip(t *testing.T, h *sam.Header, records []*sam.Record) []*sam.Record {
	t.Helper()
	var buf bytes.Buffer
	if err := WriteCRAM(&buf, h, records); err != nil {
		t.Fatalf("WriteCRAM: %v", err)
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

// assertRecordEqual compares the round-trip-significant fields of two
// records. The unmapped-read MAPQ/CIGAR caveat is applied: CRAM does not
// store those for unmapped reads.
func assertRecordEqual(t *testing.T, idx int, got, want *sam.Record) {
	t.Helper()
	if got.QName != want.QName {
		t.Errorf("record %d QName = %q, want %q", idx, got.QName, want.QName)
	}
	if got.Flag != want.Flag {
		t.Errorf("record %d Flag = %#x, want %#x", idx, got.Flag, want.Flag)
	}
	if normName(got.RName) != normName(want.RName) {
		t.Errorf("record %d RName = %q, want %q", idx, got.RName, want.RName)
	}
	if got.Pos != want.Pos {
		t.Errorf("record %d Pos = %d, want %d", idx, got.Pos, want.Pos)
	}
	if normName(got.RNext) != normName(want.RNext) {
		t.Errorf("record %d RNext = %q, want %q", idx, got.RNext, want.RNext)
	}
	if got.PNext != want.PNext {
		t.Errorf("record %d PNext = %d, want %d", idx, got.PNext, want.PNext)
	}
	if got.TLen != want.TLen {
		t.Errorf("record %d TLen = %d, want %d", idx, got.TLen, want.TLen)
	}
	gotSeq := normSeq(got.Seq)
	wantSeq := normSeq(want.Seq)
	if gotSeq != wantSeq {
		t.Errorf("record %d Seq = %q, want %q", idx, gotSeq, wantSeq)
	}
	if !qualEqual(got.Qual, want.Qual, len(wantSeq)) {
		t.Errorf("record %d Qual = %v, want %v", idx, got.Qual, want.Qual)
	}
	if want.Flag&sam.FlagUnmapped == 0 {
		if got.MapQ != want.MapQ {
			t.Errorf("record %d MapQ = %d, want %d", idx, got.MapQ, want.MapQ)
		}
		if got.Cigar.String() != want.Cigar.String() {
			t.Errorf("record %d Cigar = %q, want %q", idx, got.Cigar.String(), want.Cigar.String())
		}
	}
	if !auxEqual(got.Aux, want.Aux) {
		t.Errorf("record %d Aux = %#v, want %#v", idx, got.Aux, want.Aux)
	}
}

// normSeq normalises the SAM "*" / "" no-sequence forms to "".
func normSeq(s string) string {
	if s == "*" {
		return ""
	}
	return s
}

// normName normalises the SAM "*" no-reference name to the empty string
// the reader produces.
func normName(s string) string {
	if s == "*" {
		return ""
	}
	return s
}

// qualEqual compares two quality slices, treating nil/empty/all-0xff as
// equivalent "no quality".
func qualEqual(got, want []byte, readLen int) bool {
	gn := qualIsAbsent(got)
	wn := qualIsAbsent(want)
	if gn && wn {
		return true
	}
	if gn != wn {
		return false
	}
	if len(got) != len(want) {
		return false
	}
	for i := range got {
		if got[i] != want[i] {
			return false
		}
	}
	return true
}

// qualIsAbsent reports whether a quality slice is the SAM "no quality"
// state: nil, empty, or every byte 0xff.
func qualIsAbsent(q []byte) bool {
	if len(q) == 0 {
		return true
	}
	for _, b := range q {
		if b != 0xff {
			return false
		}
	}
	return true
}

// auxEqual compares two aux-tag slices for round-trip equality, ignoring
// order differences only when the underlying values match by tag.
func auxEqual(got, want []sam.Aux) bool {
	if len(got) != len(want) {
		return false
	}
	for i := range want {
		if !reflect.DeepEqual(got[i], want[i]) {
			return false
		}
	}
	return true
}

// TestWriteCRAMEmptyInput writes a file with no records and confirms it
// reads back as an empty stream with the header intact.
func TestWriteCRAMEmptyInput(t *testing.T) {
	h := writerTestHeader()
	var buf bytes.Buffer
	if err := WriteCRAM(&buf, h, nil); err != nil {
		t.Fatalf("WriteCRAM: %v", err)
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
		t.Fatalf("got %d records from an empty file, want 0", len(out))
	}
}

// TestWriteCRAMSingleMapped round-trips one plain mapped read.
func TestWriteCRAMSingleMapped(t *testing.T) {
	h := writerTestHeader()
	cig, _ := sam.ParseCigar("10M")
	rec := &sam.Record{
		QName: "read1",
		Flag:  0,
		RName: "chr1",
		Pos:   100,
		MapQ:  60,
		Cigar: cig,
		RNext: "*",
		Seq:   "ACGTACGTAC",
		Qual:  []byte{30, 31, 32, 33, 34, 35, 36, 37, 38, 39},
	}
	out := roundTrip(t, h, []*sam.Record{rec})
	if len(out) != 1 {
		t.Fatalf("got %d records, want 1", len(out))
	}
	assertRecordEqual(t, 0, out[0], rec)
}

// TestWriteCRAMSingleUnmapped round-trips one unmapped read.
func TestWriteCRAMSingleUnmapped(t *testing.T) {
	h := writerTestHeader()
	rec := &sam.Record{
		QName: "unmapped1",
		Flag:  sam.FlagUnmapped,
		RName: "*",
		RNext: "*",
		Seq:   "TTTTGGGGCC",
		Qual:  []byte{20, 21, 22, 23, 24, 25, 26, 27, 28, 29},
	}
	out := roundTrip(t, h, []*sam.Record{rec})
	if len(out) != 1 {
		t.Fatalf("got %d records, want 1", len(out))
	}
	assertRecordEqual(t, 0, out[0], rec)
}

// mkRec builds a mapped record from a CIGAR string and a sequence,
// filling quality with a deterministic ascending pattern.
func mkRec(qname, rname string, pos int32, cigarStr, seq string) *sam.Record {
	cig, err := sam.ParseCigar(cigarStr)
	if err != nil {
		panic(err)
	}
	qual := make([]byte, len(seq))
	for i := range qual {
		qual[i] = byte(20 + i%30)
	}
	return &sam.Record{
		QName: qname, Flag: 0, RName: rname, Pos: pos, MapQ: 40,
		Cigar: cig, RNext: "*", Seq: seq, Qual: qual,
	}
}

// TestWriteCRAMCigarShapes round-trips reads whose CIGARs exercise
// insertions, deletions, soft clips, hard clips, reference skips and
// padding.
func TestWriteCRAMCigarShapes(t *testing.T) {
	h := writerTestHeader()
	records := []*sam.Record{
		// 4M2I4M: an internal insertion.
		mkRec("ins", "chr1", 200, "4M2I4M", "AAAACCGGGG"),
		// 5M3D5M: an internal deletion (D consumes no read bases).
		mkRec("del", "chr1", 300, "5M3D5M", "ACGTAACGTA"),
		// 3S7M: a leading soft clip.
		mkRec("soft", "chr1", 400, "3S7M", "TTTACGTACG"),
		// 2H8M: a leading hard clip (H consumes no read bases).
		mkRec("hard", "chr1", 500, "2H8M", "ACGTACGT"),
		// 4M10N4M: a reference skip.
		mkRec("skip", "chr1", 600, "4M10N4M", "ACGTACGT"),
		// 4M2P4M: padding.
		mkRec("pad", "chr1", 700, "4M2P4M", "ACGTACGT"),
		// 2S3M2D3M2H: a mix of clips, a deletion and matches.
		mkRec("mix", "chr1", 800, "2S3M2D3M2H", "TTACGACG"),
	}
	out := roundTrip(t, h, records)
	if len(out) != len(records) {
		t.Fatalf("got %d records, want %d", len(out), len(records))
	}
	for i := range records {
		assertRecordEqual(t, i, out[i], records[i])
	}
}

// TestWriteCRAMReverseStrand round-trips a reverse-strand mapped read.
func TestWriteCRAMReverseStrand(t *testing.T) {
	h := writerTestHeader()
	rec := mkRec("rev", "chr2", 1000, "12M", "ACGTACGTACGT")
	rec.Flag = sam.FlagReverse
	out := roundTrip(t, h, []*sam.Record{rec})
	assertRecordEqual(t, 0, out[0], rec)
}

// TestWriteCRAMAuxTags round-trips reads carrying auxiliary tags of
// every supported value type, including a B-array tag.
func TestWriteCRAMAuxTags(t *testing.T) {
	h := writerTestHeader()
	rec1 := mkRec("aux1", "chr1", 100, "8M", "ACGTACGT")
	rec1.Aux = []sam.Aux{
		{Tag: "NM", Type: 'i', Value: int64(3)},
		{Tag: "AS", Type: 'i', Value: int64(-7)},
		{Tag: "Xc", Type: 'c', Value: int64(-5)},
		{Tag: "XC", Type: 'C', Value: int64(200)},
		{Tag: "Xs", Type: 's', Value: int64(-1000)},
		{Tag: "XS", Type: 'S', Value: int64(60000)},
		{Tag: "XF", Type: 'f', Value: float64(2.5)},
		{Tag: "XA", Type: 'A', Value: "Q"},
		{Tag: "MD", Type: 'Z', Value: "8"},
		{Tag: "XH", Type: 'H', Value: "DEADBEEF"},
	}
	rec2 := mkRec("aux2", "chr1", 200, "8M", "TTTTGGGG")
	rec2.Aux = []sam.Aux{
		{Tag: "BI", Type: 'B', ArrayType: 'i',
			ArrayValues: []interface{}{int64(1), int64(-2), int64(3)}},
		{Tag: "BC", Type: 'B', ArrayType: 'C',
			ArrayValues: []interface{}{int64(10), int64(20), int64(255)}},
		{Tag: "BF", Type: 'B', ArrayType: 'f',
			ArrayValues: []interface{}{float64(1.5), float64(-2.25)}},
	}
	// A record with no tags, interleaved, so the tag dictionary has
	// multiple distinct combinations.
	rec3 := mkRec("aux3", "chr1", 300, "8M", "CCCCAAAA")
	records := []*sam.Record{rec1, rec2, rec3}
	out := roundTrip(t, h, records)
	if len(out) != 3 {
		t.Fatalf("got %d records, want 3", len(out))
	}
	for i := range records {
		assertRecordEqual(t, i, out[i], records[i])
	}
}

// TestWriteCRAMMatePair round-trips a properly-paired read pair on the
// same reference.
func TestWriteCRAMMatePair(t *testing.T) {
	h := writerTestHeader()
	r1 := mkRec("pair", "chr1", 100, "10M", "ACGTACGTAC")
	r1.Flag = sam.FlagPaired | sam.FlagProperPair | sam.FlagRead1 | sam.FlagMateReverse
	r1.RNext = "="
	r1.PNext = 300
	r1.TLen = 210
	r2 := mkRec("pair", "chr1", 300, "10M", "TGCATGCATG")
	r2.Flag = sam.FlagPaired | sam.FlagProperPair | sam.FlagRead2 | sam.FlagReverse
	r2.RNext = "="
	r2.PNext = 100
	r2.TLen = -210
	records := []*sam.Record{r1, r2}
	out := roundTrip(t, h, records)
	if len(out) != 2 {
		t.Fatalf("got %d records, want 2", len(out))
	}
	for i := range records {
		assertRecordEqual(t, i, out[i], records[i])
	}
}

// TestWriteCRAMMateCrossRef round-trips a pair whose mates map to
// different references.
func TestWriteCRAMMateCrossRef(t *testing.T) {
	h := writerTestHeader()
	r1 := mkRec("xref", "chr1", 100, "10M", "ACGTACGTAC")
	r1.Flag = sam.FlagPaired | sam.FlagRead1
	r1.RNext = "chr2"
	r1.PNext = 5000
	r2 := mkRec("xref", "chr2", 5000, "10M", "TGCATGCATG")
	r2.Flag = sam.FlagPaired | sam.FlagRead2
	r2.RNext = "chr1"
	r2.PNext = 100
	records := []*sam.Record{r1, r2}
	out := roundTrip(t, h, records)
	for i := range records {
		assertRecordEqual(t, i, out[i], records[i])
	}
}

// TestWriteCRAMMultiReference round-trips records spanning two
// references in one slice, forcing the multi-reference (RI series) path.
func TestWriteCRAMMultiReference(t *testing.T) {
	h := writerTestHeader()
	records := []*sam.Record{
		mkRec("r1", "chr1", 100, "8M", "ACGTACGT"),
		mkRec("r2", "chr2", 200, "8M", "TTTTGGGG"),
		mkRec("r3", "chr1", 300, "8M", "CCCCAAAA"),
		{QName: "r4", Flag: sam.FlagUnmapped, RName: "*", RNext: "*",
			Seq: "GGGGCCCC", Qual: []byte{20, 21, 22, 23, 24, 25, 26, 27}},
	}
	out := roundTrip(t, h, records)
	if len(out) != len(records) {
		t.Fatalf("got %d records, want %d", len(out), len(records))
	}
	for i := range records {
		assertRecordEqual(t, i, out[i], records[i])
	}
}

// TestWriteCRAMMultipleContainers forces several containers by capping
// the slice size, confirming the record counter and container chaining
// are correct across boundaries.
func TestWriteCRAMMultipleContainers(t *testing.T) {
	h := writerTestHeader()
	var records []*sam.Record
	for i := 0; i < 25; i++ {
		records = append(records, mkRec(
			"read", "chr1", int32(100+i*10), "6M", "ACGTAC"))
	}
	var buf bytes.Buffer
	rw, err := NewRecordWriter(&buf, h)
	if err != nil {
		t.Fatalf("NewRecordWriter: %v", err)
	}
	rw.recordsPerSlice = 10 // force three containers.
	for i, rec := range records {
		if err := rw.Write(rec); err != nil {
			t.Fatalf("Write record %d: %v", i, err)
		}
	}
	if err := rw.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	rr, err := NewRecordReader(bytes.NewReader(buf.Bytes()))
	if err != nil {
		t.Fatalf("NewRecordReader: %v", err)
	}
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
}

// TestWriteCRAMNoQuality round-trips a read whose quality is absent
// (the SAM "*" form), confirming it stays absent.
func TestWriteCRAMNoQuality(t *testing.T) {
	h := writerTestHeader()
	rec := &sam.Record{
		QName: "noqual", Flag: 0, RName: "chr1", Pos: 100, MapQ: 30,
		RNext: "*", Seq: "ACGTACGT",
	}
	rec.Cigar, _ = sam.ParseCigar("8M")
	out := roundTrip(t, h, []*sam.Record{rec})
	if !qualIsAbsent(out[0].Qual) {
		t.Errorf("quality = %v, want absent", out[0].Qual)
	}
	assertRecordEqual(t, 0, out[0], rec)
}

// TestCreateCRAMFile exercises the file-path convenience constructor.
func TestCreateCRAMFile(t *testing.T) {
	h := writerTestHeader()
	path := t.TempDir() + "/out.cram"
	rw, err := CreateCRAM(path, h)
	if err != nil {
		t.Fatalf("CreateCRAM: %v", err)
	}
	rec := mkRec("file1", "chr1", 100, "8M", "ACGTACGT")
	if err := rw.Write(rec); err != nil {
		t.Fatalf("Write: %v", err)
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
	if len(out) != 1 {
		t.Fatalf("got %d records, want 1", len(out))
	}
	assertRecordEqual(t, 0, out[0], rec)
}

// TestWriteCRAMRejectsUnknownReference confirms a record naming a
// reference absent from the header is rejected, not silently mangled.
func TestWriteCRAMRejectsUnknownReference(t *testing.T) {
	h := writerTestHeader()
	var buf bytes.Buffer
	rw, err := NewRecordWriter(&buf, h)
	if err != nil {
		t.Fatalf("NewRecordWriter: %v", err)
	}
	rec := mkRec("bad", "chrX", 100, "8M", "ACGTACGT")
	if err := rw.Write(rec); err == nil {
		t.Fatalf("Write of an unknown-reference record should fail")
	}
}

// TestWriteCRAMErrorPaths exercises the writer's argument-validation and
// lifecycle error paths.
func TestWriteCRAMErrorPaths(t *testing.T) {
	h := writerTestHeader()

	t.Run("nil record", func(t *testing.T) {
		var buf bytes.Buffer
		rw, _ := NewRecordWriter(&buf, h)
		if err := rw.Write(nil); err == nil {
			t.Error("Write(nil) should fail")
		}
	})

	t.Run("write after close", func(t *testing.T) {
		var buf bytes.Buffer
		rw, _ := NewRecordWriter(&buf, h)
		if err := rw.Close(); err != nil {
			t.Fatalf("Close: %v", err)
		}
		if err := rw.Write(mkRec("late", "chr1", 1, "4M", "ACGT")); err == nil {
			t.Error("Write after Close should fail")
		}
	})

	t.Run("double close", func(t *testing.T) {
		var buf bytes.Buffer
		rw, _ := NewRecordWriter(&buf, h)
		if err := rw.Close(); err != nil {
			t.Fatalf("first Close: %v", err)
		}
		if err := rw.Close(); err != nil {
			t.Errorf("second Close: %v", err)
		}
	})

	t.Run("cigar query mismatch", func(t *testing.T) {
		var buf bytes.Buffer
		rw, _ := NewRecordWriter(&buf, h)
		// CIGAR consumes 10 query bases but SEQ is 4 bytes.
		rec := mkRec("bad", "chr1", 1, "4M", "ACGT")
		rec.Cigar, _ = sam.ParseCigar("10M")
		if err := rw.Write(rec); err == nil {
			t.Error("CIGAR/SEQ length mismatch should be rejected")
		}
	})

	t.Run("unknown mate reference", func(t *testing.T) {
		var buf bytes.Buffer
		rw, _ := NewRecordWriter(&buf, h)
		rec := mkRec("m", "chr1", 1, "4M", "ACGT")
		rec.RNext = "chrZ"
		if err := rw.Write(rec); err == nil {
			t.Error("unknown mate reference should be rejected")
		}
	})

	t.Run("nil header tolerated", func(t *testing.T) {
		var buf bytes.Buffer
		rw, err := NewRecordWriter(&buf, nil)
		if err != nil {
			t.Fatalf("NewRecordWriter(nil header): %v", err)
		}
		// An unmapped record needs no @SQ entry.
		rec := &sam.Record{QName: "u", Flag: sam.FlagUnmapped, Seq: "ACGT",
			Qual: []byte{1, 2, 3, 4}}
		if err := rw.Write(rec); err != nil {
			t.Fatalf("Write: %v", err)
		}
		if err := rw.Close(); err != nil {
			t.Fatalf("Close: %v", err)
		}
	})
}

// failingWriter is an io.Writer that returns an error after letting a
// fixed number of bytes through, used to exercise the writer's I/O
// error paths.
type failingWriter struct {
	allow int
}

func (fw *failingWriter) Write(p []byte) (int, error) {
	if fw.allow <= 0 {
		return 0, errSimulatedIO
	}
	if len(p) > fw.allow {
		n := fw.allow
		fw.allow = 0
		return n, errSimulatedIO
	}
	fw.allow -= len(p)
	return len(p), nil
}

// errSimulatedIO is the error a failingWriter returns.
var errSimulatedIO = &simulatedIOError{}

type simulatedIOError struct{}

func (*simulatedIOError) Error() string { return "simulated I/O failure" }

// TestWriteCRAMIOErrors confirms a failure of the underlying writer is
// reported, not swallowed, by Write and Close.
func TestWriteCRAMIOErrors(t *testing.T) {
	h := writerTestHeader()

	t.Run("header write failure", func(t *testing.T) {
		_, err := NewRecordWriter(&failingWriter{allow: 0}, h)
		if err == nil {
			t.Error("NewRecordWriter should report a header write failure")
		}
	})

	t.Run("container flush failure", func(t *testing.T) {
		// Allow the file definition and header container through, then
		// fail when the data container is flushed.
		fw := &failingWriter{allow: 400}
		rw, err := NewRecordWriter(fw, h)
		if err != nil {
			t.Fatalf("NewRecordWriter: %v", err)
		}
		rw.recordsPerSlice = 1
		err = rw.Write(mkRec("r", "chr1", 1, "4M", "ACGT"))
		if err == nil {
			// The flush may instead surface from a later call; force it.
			err = rw.Close()
		}
		if err == nil {
			t.Error("a container flush failure should be reported")
		}
	})

	t.Run("WriteCRAM record error", func(t *testing.T) {
		// A record naming an unknown reference makes WriteCRAM fail with
		// the per-record error wrapper.
		var buf bytes.Buffer
		err := WriteCRAM(&buf, h, []*sam.Record{mkRec("bad", "nope", 1, "4M", "ACGT")})
		if err == nil {
			t.Error("WriteCRAM should report a bad record")
		}
	})
}

// TestWriteCRAMNarrowIntTags round-trips integer tags read at their
// narrow BAM widths ('c', 'C', 's', 'S') to confirm the width is
// preserved through the CRAM tag-encoding map.
func TestWriteCRAMNarrowIntTags(t *testing.T) {
	h := writerTestHeader()
	rec := mkRec("narrow", "chr1", 100, "4M", "ACGT")
	rec.Aux = []sam.Aux{
		{Tag: "T1", Type: 'c', Value: int64(-1)},
		{Tag: "T2", Type: 'C', Value: int64(255)},
		{Tag: "T3", Type: 's', Value: int64(-32768)},
		{Tag: "T4", Type: 'S', Value: int64(65535)},
		{Tag: "T5", Type: 'I', Value: int64(4000000000)},
	}
	out := roundTrip(t, h, []*sam.Record{rec})
	for i, a := range out[0].Aux {
		if a.Type != rec.Aux[i].Type {
			t.Errorf("tag %s type = %c, want %c", a.Tag, a.Type, rec.Aux[i].Type)
		}
	}
	assertRecordEqual(t, 0, out[0], rec)
}

// TestWriteCRAMEmptyTagValue round-trips a record with an empty 'Z'
// string tag, the smallest non-trivial tag value.
func TestWriteCRAMEmptyTagValue(t *testing.T) {
	h := writerTestHeader()
	rec := mkRec("emptytag", "chr1", 100, "4M", "ACGT")
	rec.Aux = []sam.Aux{{Tag: "ZZ", Type: 'Z', Value: ""}}
	out := roundTrip(t, h, []*sam.Record{rec})
	assertRecordEqual(t, 0, out[0], rec)
}

// TestWriteCRAMFixtureRoundTrip reads an existing samtools CRAM fixture,
// re-writes its records through the writer, reads the result back, and
// asserts the records survive the round trip.
func TestWriteCRAMFixtureRoundTrip(t *testing.T) {
	data, ok := loadFixture(t, "dat/test_input_1_a.cram")
	if !ok {
		t.Fatalf("samtools submodule not initialised — fixture unavailable; run `git submodule update --init reference_code/samtools`")
	}
	src, err := NewRecordReader(bytes.NewReader(data))
	if err != nil {
		t.Fatalf("NewRecordReader (fixture): %v", err)
	}
	original, err := src.ReadAll()
	if err != nil {
		t.Fatalf("ReadAll (fixture): %v", err)
	}
	if len(original) == 0 {
		t.Fatalf("fixture decoded to no records; the test_input_1_a.cram fixture must decode to records for the round-trip to be exercised")
	}

	var buf bytes.Buffer
	if err := WriteCRAM(&buf, src.Header(), original); err != nil {
		t.Fatalf("WriteCRAM: %v", err)
	}
	rt, err := NewRecordReader(bytes.NewReader(buf.Bytes()))
	if err != nil {
		t.Fatalf("NewRecordReader (round-trip): %v", err)
	}
	got, err := rt.ReadAll()
	if err != nil {
		t.Fatalf("ReadAll (round-trip): %v", err)
	}
	if len(got) != len(original) {
		t.Fatalf("round-trip record count = %d, want %d", len(got), len(original))
	}
	for i := range original {
		assertRecordEqual(t, i, got[i], original[i])
	}
}
