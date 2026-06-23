package cram

import (
	"bytes"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// knownRefBases is a tiny synthetic reference whose MD5 is fixed and easy to
// reproduce: md5("ACGTACGTAC") == 45aff2fecf7615d56bc0567dffab9fa8. The bases
// are already upper-cased with no whitespace, exactly the canonical form htslib
// hashes for the @SQ M5 tag, so md5.Sum over them reproduces the upstream value.
const (
	knownRefBases = "ACGTACGTAC"
	knownRefMD5   = "45aff2fecf7615d56bc0567dffab9fa8"
	knownRefPath  = "/refs/synthetic.fa"
)

// TestCRAMHeaderM5URInjection asserts the CRAM writer fills the @SQ M5 and UR
// tags from the reference exactly as upstream htslib's cram_write_SAM_hdr does:
//   - M5 is the lower-case hex MD5 of the (upper-cased, whitespace-stripped)
//     reference bases for that contig.
//   - UR is the reference path passed via WriterOptions.ReferencePath.
//
// It also checks the writer never overwrites an M5/UR already present on an
// @SQ, and that the @SQ/@RG/@PG ordering follows upstream's emission order
// (verbatim input order with @HD hoisted first), not a regrouping.
func TestCRAMHeaderM5URInjection(t *testing.T) {
	// Input header order @HD,@SQ,@SQ,@RG,@PG with the second @SQ already
	// carrying its own M5 (which must be preserved untouched).
	const in = "@HD\tVN:1.6\tSO:coordinate\n" +
		"@SQ\tSN:c1\tLN:10\n" +
		"@SQ\tSN:c2\tLN:10\tM5:deadbeefdeadbeefdeadbeefdeadbeef\n" +
		"@RG\tID:rg1\tSM:s1\n" +
		"@PG\tID:p1\tPN:p1\n"
	h, err := sam.ParseHeaderText(in)
	if err != nil {
		t.Fatalf("ParseHeaderText: %v", err)
	}

	// Only c1 has bases supplied; c2 keeps its existing M5. Both get UR.
	ref := map[string][]byte{"c1": []byte(knownRefBases)}

	var buf bytes.Buffer
	rw, err := NewRecordWriterOpts(&buf, h, WriterOptions{
		Reference:     ref,
		ReferencePath: knownRefPath,
	})
	if err != nil {
		t.Fatalf("NewRecordWriterOpts: %v", err)
	}
	if err := rw.Close(); err != nil {
		t.Fatalf("RecordWriter.Close: %v", err)
	}

	rr, err := NewRecordReader(bytes.NewReader(buf.Bytes()))
	if err != nil {
		t.Fatalf("NewRecordReader: %v", err)
	}
	got := rr.Header().Text()

	// The whole embedded header must equal the upstream-correct text: @HD
	// first, @SQ lines keep their input position (not regrouped after @RG), c1
	// gains M5(known)+UR, c2 keeps its existing M5 and gains only UR.
	want := "@HD\tVN:1.6\tSO:coordinate\n" +
		"@SQ\tSN:c1\tLN:10\tM5:" + knownRefMD5 + "\tUR:" + knownRefPath + "\n" +
		"@SQ\tSN:c2\tLN:10\tM5:deadbeefdeadbeefdeadbeefdeadbeef\tUR:" + knownRefPath + "\n" +
		"@RG\tID:rg1\tSM:s1\n" +
		"@PG\tID:p1\tPN:p1\n"
	if got != want {
		t.Fatalf("embedded CRAM header mismatch:\n got=%q\nwant=%q", got, want)
	}
}

// TestCRAMHeaderURWithoutBases checks that ReferencePath alone injects UR onto
// every @SQ even when no bases are supplied (so no M5 can be computed) — the
// two tags are filled independently, mirroring upstream where UR comes from the
// -T path and M5 from the loaded sequence.
func TestCRAMHeaderURWithoutBases(t *testing.T) {
	const in = "@HD\tVN:1.6\n@SQ\tSN:c1\tLN:10\n"
	h, err := sam.ParseHeaderText(in)
	if err != nil {
		t.Fatalf("ParseHeaderText: %v", err)
	}
	var buf bytes.Buffer
	if err := WriteCRAMRefPath(&buf, h, nil, knownRefPath); err != nil {
		t.Fatalf("write CRAM: %v", err)
	}
	rr, err := NewRecordReader(bytes.NewReader(buf.Bytes()))
	if err != nil {
		t.Fatalf("NewRecordReader: %v", err)
	}
	got := rr.Header().Text()
	want := "@HD\tVN:1.6\n@SQ\tSN:c1\tLN:10\tUR:" + knownRefPath + "\n"
	if got != want {
		t.Fatalf("UR-only header mismatch:\n got=%q\nwant=%q", got, want)
	}
}

// TestCRAMReferenceFreeHeaderUnchanged guards the reference-free path: with
// neither bases nor a reference path the embedded header is left byte-identical
// (no M5/UR injection), so existing reference-free callers are unaffected.
func TestCRAMReferenceFreeHeaderUnchanged(t *testing.T) {
	const in = "@HD\tVN:1.6\n@SQ\tSN:c1\tLN:10\n@RG\tID:rg1\n"
	h, err := sam.ParseHeaderText(in)
	if err != nil {
		t.Fatalf("ParseHeaderText: %v", err)
	}
	var buf bytes.Buffer
	if err := WriteCRAM(&buf, h, nil); err != nil {
		t.Fatalf("WriteCRAM: %v", err)
	}
	rr, err := NewRecordReader(bytes.NewReader(buf.Bytes()))
	if err != nil {
		t.Fatalf("NewRecordReader: %v", err)
	}
	if got := rr.Header().Text(); got != in {
		t.Fatalf("reference-free header changed:\n got=%q\nwant=%q", got, in)
	}
}

// TestCRAMHeaderM5RoundTripRecords confirms M5/UR injection does not disturb
// the record stream: a CRAM written with a reference (and thus M5/UR in its
// @SQ lines) decodes to exactly the records that went in. The M5 lives only in
// the header, so the encoded records must be untouched.
func TestCRAMHeaderM5RoundTripRecords(t *testing.T) {
	const in = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c1\tLN:10\n"
	h, err := sam.ParseHeaderText(in)
	if err != nil {
		t.Fatalf("ParseHeaderText: %v", err)
	}
	// A single unmapped read keeps its bases literal, so the round-trip is
	// reference-independent on decode while the writer still injects M5/UR.
	rec := &sam.Record{
		QName: "r1",
		Flag:  sam.FlagUnmapped,
		RName: "*",
		Seq:   "ACGTACGTAC",
		Qual:  []byte{40, 40, 40, 40, 40, 40, 40, 40, 40, 40},
	}
	ref := map[string][]byte{"c1": []byte(knownRefBases)}

	var withRef, refFree bytes.Buffer
	if err := WriteCRAMRefPath(&withRef, h, []*sam.Record{rec}, knownRefPath, ref); err != nil {
		t.Fatalf("write CRAM with ref: %v", err)
	}
	if err := WriteCRAM(&refFree, h, []*sam.Record{rec}); err != nil {
		t.Fatalf("write reference-free CRAM: %v", err)
	}

	// The M5/UR-bearing CRAM must decode to the same record as the plain one.
	got := decodeOne(t, withRef.Bytes())
	wantRec := decodeOne(t, refFree.Bytes())
	if got.QName != wantRec.QName || got.Seq != wantRec.Seq || !bytes.Equal(got.Qual, wantRec.Qual) || got.Flag != wantRec.Flag {
		t.Fatalf("record changed by M5 injection:\n got=%+v\nwant=%+v", got, wantRec)
	}
	if got.Seq != rec.Seq || !bytes.Equal(got.Qual, rec.Qual) {
		t.Fatalf("record not round-tripped: got Seq=%q Qual=%v", got.Seq, got.Qual)
	}

	// And its header must carry the injected M5/UR.
	rr, err := NewRecordReader(bytes.NewReader(withRef.Bytes()))
	if err != nil {
		t.Fatalf("NewRecordReader: %v", err)
	}
	hdr := rr.Header().Text()
	if !strings.Contains(hdr, "M5:"+knownRefMD5) || !strings.Contains(hdr, "UR:"+knownRefPath) {
		t.Fatalf("expected M5/UR in round-tripped header, got %q", hdr)
	}
}

// decodeOne reads exactly one record from a CRAM byte stream.
func decodeOne(t *testing.T, b []byte) *sam.Record {
	t.Helper()
	rr, err := NewRecordReader(bytes.NewReader(b))
	if err != nil {
		t.Fatalf("NewRecordReader: %v", err)
	}
	rec, err := rr.Read()
	if err != nil {
		t.Fatalf("Read: %v", err)
	}
	return rec
}

// WriteCRAMRefPath is a test helper that writes a complete CRAM file with the
// given reference path (for UR) and optional reference base map (for M5). It
// keeps the regression tests terse without widening the public WriteCRAM API.
func WriteCRAMRefPath(w *bytes.Buffer, h *sam.Header, records []*sam.Record, refPath string, refMaps ...map[string][]byte) error {
	var ref map[string][]byte
	if len(refMaps) > 0 {
		ref = refMaps[0]
	}
	rw, err := NewRecordWriterOpts(w, h, WriterOptions{Reference: ref, ReferencePath: refPath})
	if err != nil {
		return err
	}
	for _, rec := range records {
		if err := rw.Write(rec); err != nil {
			return err
		}
	}
	return rw.Close()
}
