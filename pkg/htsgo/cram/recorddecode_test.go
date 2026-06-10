package cram

import (
	"bytes"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// wantInput1ARecords is the expected text-SAM body of decoding
// dat/test_input_1_a.cram. It is byte-for-byte equal to `samtools view`
// of that fixture (the parity oracle): the fields match the sibling
// dat/test_input_1_a.sam except for the RG/PG tag order (RG comes from a
// CRAM data series and is appended after the dictionary tags) and record
// 15 (the unmapped read u1, which CRAM stores without mapping quality or
// CIGAR — decoded as "u1 4 * 0 0 * * 0 0 ..."). The CRAM file is
// reference-free (RR preservation flag false), so every base is
// recovered with no external reference.
// The tag order mirrors samtools view of this fixture: a record's RG tag
// comes from a dedicated CRAM data series (not the tag dictionary), so on
// decode it is appended after the dictionary tags (e.g. PG) rather than
// in the original SAM's position — exactly htslib's behaviour.
var wantInput1ARecords = []string{
	"r000\t99\tinsert\t50\t30\t10M\t=\t80\t30\tATTTAGCTAC\tAAAAAAAAAA\tPG:Z:bull\tRG:Z:cow",
	"r000\t211\tinsert\t80\t30\t10M\t=\t50\t-30\tCCCAATCATT\tAAAAAAAAAA\tPG:Z:bull\tRG:Z:cow",
	"r001\t163\tref1\t7\t30\t8M4I4M1D3M\t=\t37\t39\tTTAGATAAAGAGGATACTG\t*\tXX:B:S,12561,2,20,112\tYY:i:100\tPG:Z:colt\tRG:Z:fish",
	"r002\t0\tref1\t9\t30\t1S2I6M1P1I1P1I4M2I\t*\t0\t0\tAAAAGATAAGGGATAAA\t*\tXA:Z:abc\tXB:i:-10\tPG:Z:colt",
	"r003\t0\tref1\t9\t30\t5H6M\t*\t0\t0\tAGCTAA\t*\tRG:Z:cow",
	"r004\t0\tref1\t16\t30\t6M14N1I5M\t*\t0\t0\tATAGCTCTCAGC\t*\tPG:Z:colt\tRG:Z:colt",
	"r003\t16\tref1\t29\t30\t6H5M\t*\t0\t0\tTAGGC\t*\tPG:Z:colt\tRG:Z:cow",
	"r001\t83\tref1\t37\t30\t9M\t=\t7\t-39\tCAGCGCCAT\t*\tPG:Z:colt\tRG:Z:fish",
	"x1\t0\tref2\t1\t30\t20M\t*\t0\t0\tAGGTTTTATAAAACAAATAA\t*\tPG:Z:bull\tRG:Z:colt",
	"x2\t0\tref2\t2\t30\t21M\t*\t0\t0\tGGTTTTATAAAACAAATAATT\t?????????????????????\tPG:Z:bull\tRG:Z:colt",
	"x3\t0\tref2\t6\t30\t9M4I13M\t*\t0\t0\tTTATAAAACAAATAATTAAGTCTACA\t??????????????????????????\tPG:Z:bull\tRG:Z:fish",
	"x4\t0\tref2\t10\t30\t25M\t*\t0\t0\tCAAATAATTAAGTCTACAGAGCAAC\t?????????????????????????\tPG:Z:bull\tRG:Z:fish",
	"x5\t0\tref2\t12\t30\t24M\t*\t0\t0\tAATAATTAAGTCTACAGAGCAACT\t????????????????????????\tPG:Z:bull\tRG:Z:fish",
	"x6\t0\tref2\t14\t30\t23M\t*\t0\t0\tTAATTAAGTCTACAGAGCAACTA\t???????????????????????\tRG:Z:cow",
	"u1\t4\t*\t0\t0\t*\t*\t0\t0\tTAATTAAGTCTACAGAAAAAAAA\t???????????????????????",
}

// TestDecodeInput1AToSAM is the C4b correctness oracle: it decodes the
// real CRAM v3.0 fixture dat/test_input_1_a.cram record by record and
// asserts the reconstructed SAM matches the expected body. The fixture
// is reference-free, so the decode must reproduce every base, quality,
// CIGAR and tag with no external reference.
func TestDecodeInput1AToSAM(t *testing.T) {
	data, ok := loadFixture(t, "dat/test_input_1_a.cram")
	if !ok {
		t.Skip("samtools submodule not initialised — fixture unavailable")
	}
	rr, err := NewRecordReader(bytes.NewReader(data))
	if err != nil {
		t.Fatalf("NewRecordReader: %v", err)
	}
	if rr.Header() == nil {
		t.Fatal("expected an embedded SAM header")
	}
	// The embedded header carries the four @SQ references the records
	// align against.
	if got := len(rr.Header().Refs); got != 4 {
		t.Fatalf("expected 4 @SQ references, got %d", got)
	}

	recs, err := rr.ReadAll()
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if len(recs) != len(wantInput1ARecords) {
		t.Fatalf("decoded %d records, want %d", len(recs), len(wantInput1ARecords))
	}
	if rr.NeedsReference() {
		t.Error("reference-free CRAM must not need an external reference")
	}
	for i, rec := range recs {
		got := formatRecord(rec)
		if got != wantInput1ARecords[i] {
			t.Errorf("record %d mismatch:\n got %q\nwant %q", i, got, wantInput1ARecords[i])
		}
	}
}

// TestWriteSAMInput1A exercises the WriteSAM convenience path and checks
// that the emitted body equals the per-record oracle.
func TestWriteSAMInput1A(t *testing.T) {
	data, ok := loadFixture(t, "dat/test_input_1_a.cram")
	if !ok {
		t.Skip("samtools submodule not initialised — fixture unavailable")
	}
	rr, err := NewRecordReader(bytes.NewReader(data))
	if err != nil {
		t.Fatalf("NewRecordReader: %v", err)
	}
	var buf bytes.Buffer
	if err := rr.WriteSAM(&buf); err != nil {
		t.Fatalf("WriteSAM: %v", err)
	}
	lines := strings.Split(strings.TrimRight(buf.String(), "\n"), "\n")
	var body []string
	for _, l := range lines {
		if !strings.HasPrefix(l, "@") {
			body = append(body, l)
		}
	}
	if len(body) != len(wantInput1ARecords) {
		t.Fatalf("WriteSAM emitted %d records, want %d", len(body), len(wantInput1ARecords))
	}
	for i, l := range body {
		if l != wantInput1ARecords[i] {
			t.Errorf("WriteSAM record %d:\n got %q\nwant %q", i, l, wantInput1ARecords[i])
		}
	}
	// The first emitted line must be the SAM @HD header.
	if !strings.HasPrefix(buf.String(), "@HD\t") {
		t.Error("WriteSAM output must start with the @HD header line")
	}
}

// TestDecodeQuickcheckStructural structurally decodes the reference-
// backed fixture 7.quickcheck.cram30.ok.cram: every record traverses
// without panic or error, and the decoder reports that an external
// reference is needed (the file's RR preservation flag is true).
func TestDecodeQuickcheckStructural(t *testing.T) {
	data, ok := loadFixture(t, "quickcheck/7.quickcheck.cram30.ok.cram")
	if !ok {
		t.Skip("samtools submodule not initialised — fixture unavailable")
	}
	rr, err := NewRecordReader(bytes.NewReader(data))
	if err != nil {
		t.Fatalf("NewRecordReader: %v", err)
	}
	recs, err := rr.ReadAll()
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if len(recs) == 0 {
		t.Fatal("expected the quickcheck fixture to carry records")
	}
	if !rr.NeedsReference() {
		t.Error("reference-backed CRAM should report NeedsReference")
	}
	// Every record must have a consistent SEQ/CIGAR query length even
	// though the bases an external reference would supply are 'N'.
	for i, rec := range recs {
		if len(rec.Cigar) > 0 && rec.Seq != "" {
			if ql := rec.Cigar.QueryLength(); ql != len(rec.Seq) {
				t.Errorf("record %d: CIGAR query length %d != SEQ length %d", i, ql, len(rec.Seq))
			}
		}
	}
}

// formatRecord renders one sam.Record as a tab-delimited SAM line using
// the sam package's writer, for comparison against the text oracle.
func formatRecord(rec *sam.Record) string {
	var buf bytes.Buffer
	w := sam.NewSAMWriter(&buf)
	if err := w.Write(rec); err != nil {
		return "write-error: " + err.Error()
	}
	w.Close()
	return strings.TrimRight(buf.String(), "\n")
}

// TestParseTagDictionary checks the preservation-map tag-dictionary
// parser against the layouts CRAM uses: a single list, several lists,
// and an empty list.
func TestParseTagDictionary(t *testing.T) {
	td, err := parseTagDictionary([]byte("PGZ\x00XAZXBcPGZ\x00\x00PGZ\x00"))
	if err != nil {
		t.Fatalf("parseTagDictionary: %v", err)
	}
	if len(td) != 4 {
		t.Fatalf("expected 4 dictionary entries, got %d", len(td))
	}
	if len(td[0]) != 1 || td[0][0].String() != "PGZ" {
		t.Errorf("entry 0 = %v, want [PGZ]", td[0])
	}
	if len(td[1]) != 3 {
		t.Errorf("entry 1 has %d tags, want 3", len(td[1]))
	}
	if len(td[2]) != 0 {
		t.Errorf("entry 2 should be empty, got %v", td[2])
	}
	if _, err := parseTagDictionary([]byte("PG")); err == nil {
		t.Error("expected an error for a truncated tag key")
	}
	if td, err := parseTagDictionary(nil); err != nil || td != nil {
		t.Errorf("empty dictionary should yield (nil,nil), got (%v,%v)", td, err)
	}
}

// TestDecodeTagValue checks the per-type tag-value decoder.
func TestDecodeTagValue(t *testing.T) {
	cases := []struct {
		key  tagKey
		raw  []byte
		want string
	}{
		{tagKey{'N', 'M', 'i'}, []byte{3, 0, 0, 0}, "NM:i:3"},
		{tagKey{'X', 'C', 'c'}, []byte{0xfb}, "XC:i:-5"},
		{tagKey{'X', 'A', 'Z'}, []byte("hi\x00"), "XA:Z:hi"},
		{tagKey{'X', 'F', 'f'}, []byte{0x00, 0x00, 0x80, 0x3f}, "XF:f:1"},
		{tagKey{'X', 'M', 'A'}, []byte{'Q'}, "XM:A:Q"},
		{tagKey{'X', 'B', 'B'}, []byte{'S', 2, 0, 0, 0, 0x11, 0x31, 0x02, 0x00}, "XB:B:S,12561,2"},
	}
	for _, c := range cases {
		aux, err := decodeTagValue(c.key, c.raw)
		if err != nil {
			t.Errorf("decodeTagValue(%s): %v", c.key, err)
			continue
		}
		if got := aux.FormatSAM(); got != c.want {
			t.Errorf("decodeTagValue(%s) = %q, want %q", c.key, got, c.want)
		}
	}
	if _, err := decodeTagValue(tagKey{'N', 'M', 'i'}, []byte{1}); err == nil {
		t.Error("expected an error for a truncated 'i' value")
	}
	if _, err := decodeTagValue(tagKey{'X', 'X', '?'}, []byte{0}); err == nil {
		t.Error("expected an error for an unknown tag type")
	}
}

// TestCigarBuilder checks the run-coalescing CIGAR builder.
func TestCigarBuilder(t *testing.T) {
	var b cigarBuilder
	b.add(sam.CigarMatch, 4)
	b.add(sam.CigarMatch, 6) // merges into the previous M.
	b.add(sam.CigarInsertion, 2)
	b.add(sam.CigarDeletion, 0) // dropped: zero length.
	b.add(sam.CigarMatch, 3)
	if got := b.cigar().String(); got != "10M2I3M" {
		t.Errorf("cigar = %q, want %q", got, "10M2I3M")
	}
}

// TestMergeAux checks that the data-series read-group tag is appended
// after the dictionary tags (htslib's order), is suppressed when the
// dictionary already carries an RG tag, and adds nothing when nil.
func TestMergeAux(t *testing.T) {
	rg := &sam.Aux{Tag: "RG", Type: 'Z', Value: "grp"}
	xx := sam.Aux{Tag: "XX", Type: 'i', Value: int64(1)}
	pg := sam.Aux{Tag: "PG", Type: 'Z', Value: "prog"}

	// RG from the data series is appended last, after the dictionary tags.
	got := mergeAux([]sam.Aux{xx, pg}, rg)
	if len(got) != 3 || got[2].Tag != "RG" {
		t.Errorf("RG should be appended after the dictionary tags, got %v", auxTags(got))
	}
	got = mergeAux([]sam.Aux{xx}, rg)
	if len(got) != 2 || got[1].Tag != "RG" {
		t.Errorf("RG should be appended last, got %v", auxTags(got))
	}
	// A dictionary that already carries RG suppresses the data-series RG.
	dictRG := sam.Aux{Tag: "RG", Type: 'Z', Value: "dict"}
	got = mergeAux([]sam.Aux{dictRG, xx}, rg)
	if len(got) != 2 || got[0].Value != "dict" {
		t.Errorf("a dictionary RG should be kept and the data-series RG dropped, got %v", auxTags(got))
	}
	if got := mergeAux([]sam.Aux{xx}, nil); len(got) != 1 {
		t.Errorf("nil RG should add nothing, got %v", auxTags(got))
	}
}

// auxTags lists the tag names of an aux slice for diagnostics.
func auxTags(aux []sam.Aux) []string {
	out := make([]string, len(aux))
	for i, a := range aux {
		out[i] = a.Tag
	}
	return out
}
