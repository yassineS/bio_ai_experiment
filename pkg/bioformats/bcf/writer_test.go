package bcf

import (
	"bytes"
	"math"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/vcf"
)

// makeHeader builds the same in-memory Header used across the test suite,
// constructed by parsing the canonical text header in buildBCFStream.
func makeHeader(t *testing.T) *Header {
	t.Helper()
	const text = "##fileformat=VCFv4.3\n" +
		"##contig=<ID=chr1,length=200>\n" +
		"##contig=<ID=chr2,length=100>\n" +
		"##FILTER=<ID=q10,Description=\"quality below 10\">\n" +
		"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">\n" +
		"##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele freq\">\n" +
		"##INFO=<ID=TAG,Number=1,Type=String,Description=\"A tag\">\n" +
		"##INFO=<ID=H2,Number=0,Type=Flag,Description=\"HapMap2 membership\">\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
		"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
	hdr, err := parseTextHeader(text)
	if err != nil {
		t.Fatalf("parseTextHeader: %v", err)
	}
	hdr.Text = text
	return hdr
}

// TestWriterRoundTripVariant emits two records via the variant path and reads
// them back through Reader; it asserts every observable field is preserved.
func TestWriterRoundTripVariant(t *testing.T) {
	hdr := makeHeader(t)
	var buf bytes.Buffer
	w := NewWriter(&buf, hdr)

	v1 := &vcf.Variant{
		Chrom: "chr1", Pos: 100, ID: "rs1",
		Ref: "A", Alt: []string{"T"},
		Qual: 30, Filter: []string{"PASS"},
		Info: map[string]string{"DP": "50", "AF": "0.25", "H2": ""},
	}
	v2 := &vcf.Variant{
		Chrom: "chr2", Pos: 25, ID: ".",
		Ref: "G", Alt: []string{"C", "T"},
		Qual:   -1, // missing
		Filter: []string{"q10"},
		Info:   map[string]string{"DP": "10"},
	}
	for _, v := range []*vcf.Variant{v1, v2} {
		if err := w.Write(v); err != nil {
			t.Fatalf("Write: %v", err)
		}
	}
	if err := w.Flush(); err != nil {
		t.Fatalf("Flush: %v", err)
	}

	r, err := NewReader(&buf)
	if err != nil {
		t.Fatalf("NewReader: %v", err)
	}
	recs, err := r.ReadAll()
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if len(recs) != 2 {
		t.Fatalf("got %d records, want 2", len(recs))
	}

	gotV1 := recs[0].ToVariant(r.Header())
	if gotV1.Chrom != "chr1" || gotV1.Pos != 100 || gotV1.ID != "rs1" || gotV1.Ref != "A" || gotV1.Alt[0] != "T" {
		t.Errorf("v1 round-trip mismatch: %+v", gotV1)
	}
	if gotV1.Info["DP"] != "50" {
		t.Errorf("v1 INFO/DP: got %q want 50", gotV1.Info["DP"])
	}
	if gotV1.Info["AF"] != "0.25" {
		t.Errorf("v1 INFO/AF: got %q want 0.25", gotV1.Info["AF"])
	}

	gotV2 := recs[1].ToVariant(r.Header())
	if gotV2.Chrom != "chr2" || gotV2.Pos != 25 || gotV2.Ref != "G" {
		t.Errorf("v2 round-trip mismatch: %+v", gotV2)
	}
	if gotV2.Qual != -1 {
		t.Errorf("v2 missing qual got %v", gotV2.Qual)
	}
	if len(gotV2.Alt) != 2 || gotV2.Alt[0] != "C" || gotV2.Alt[1] != "T" {
		t.Errorf("v2 alts: %v", gotV2.Alt)
	}
	if len(gotV2.Filter) == 0 || gotV2.Filter[0] != "q10" {
		t.Errorf("v2 filter: %v", gotV2.Filter)
	}
}

// TestWriterPerSampleRoundTrip exercises GT and FORMAT/DP across two samples.
func TestWriterPerSampleRoundTrip(t *testing.T) {
	hdr := makeHeader(t)
	var buf bytes.Buffer
	w := NewWriter(&buf, hdr)

	v := &vcf.Variant{
		Chrom: "chr1", Pos: 100, ID: ".",
		Ref: "A", Alt: []string{"T"},
		Qual: 99, Filter: []string{"PASS"},
		Format: []string{"GT", "DP"},
		Samples: []vcf.Sample{
			{Name: "S1", Data: map[string]string{"GT": "0/1", "DP": "20"}},
			{Name: "S2", Data: map[string]string{"GT": "1|1", "DP": "30"}},
		},
	}
	if err := w.Write(v); err != nil {
		t.Fatalf("Write: %v", err)
	}
	if err := w.Flush(); err != nil {
		t.Fatalf("Flush: %v", err)
	}

	r, err := NewReader(&buf)
	if err != nil {
		t.Fatal(err)
	}
	recs, err := r.ReadAll()
	if err != nil || len(recs) != 1 {
		t.Fatalf("ReadAll: recs=%d err=%v", len(recs), err)
	}
	got := recs[0].ToVariant(r.Header())
	if len(got.Samples) != 2 {
		t.Fatalf("got %d samples", len(got.Samples))
	}
	if got.Samples[0].Data["GT"] != "0/1" {
		t.Errorf("S1 GT: %q", got.Samples[0].Data["GT"])
	}
	if got.Samples[1].Data["GT"] != "1|1" {
		t.Errorf("S2 GT: %q", got.Samples[1].Data["GT"])
	}
	if got.Samples[0].Data["DP"] != "20" || got.Samples[1].Data["DP"] != "30" {
		t.Errorf("DPs: %q / %q", got.Samples[0].Data["DP"], got.Samples[1].Data["DP"])
	}
}

// TestWriterWriteRecordIdentity checks that WriteRecord re-emits a decoded
// Record without altering observable fields after a parse-write-parse cycle.
func TestWriterWriteRecordIdentity(t *testing.T) {
	// Build a stream via the test fixture builder so we have real Records.
	rec := buildRecord(recordSpec{
		chromID:  0,
		pos:      199,
		rlen:     1,
		qual:     45,
		id:       "rsX",
		alleles:  []string{"C", "G"},
		filters:  []int32{0},
		infoKeys: []int32{2, 3},
		infoVals: [][]byte{EncodeTypedInt8(7), EncodeTypedFloat(0.5)},
	})
	stream := buildBCFStream(t, rec)
	r, err := NewReader(bytes.NewReader(stream))
	if err != nil {
		t.Fatal(err)
	}
	recs, err := r.ReadAll()
	if err != nil || len(recs) != 1 {
		t.Fatalf("read: %v", err)
	}

	// Now write it back out using WriteRecord and re-parse.
	var buf bytes.Buffer
	w := NewWriter(&buf, r.Header())
	if err := w.WriteRecord(recs[0]); err != nil {
		t.Fatalf("WriteRecord: %v", err)
	}
	if err := w.Flush(); err != nil {
		t.Fatal(err)
	}
	r2, err := NewReader(&buf)
	if err != nil {
		t.Fatal(err)
	}
	out, err := r2.ReadAll()
	if err != nil || len(out) != 1 {
		t.Fatalf("re-read: %v", err)
	}
	if out[0].ChromID != recs[0].ChromID || out[0].Pos != recs[0].Pos || out[0].ID != recs[0].ID {
		t.Errorf("identity broken: %+v vs %+v", out[0], recs[0])
	}
	if len(out[0].Alleles) != len(recs[0].Alleles) || out[0].Alleles[0] != "C" {
		t.Errorf("alleles broken: %v", out[0].Alleles)
	}
}

// TestWriterWriteRecord_FormatRoundtrip pins WriteRecord through a
// FORMAT-bearing record so the encoder's per-sample dimension survives
// the parse-write-parse cycle. Regression for the wave-21 review's
// finding that `encodeTypedValue` was using `tv.Length` (per-sample
// dim) as the total payload count and truncating FORMAT payloads.
//
// The fixture is 2 samples × diploid GT + scalar DP: we hand-assemble
// the indiv block matching the header dictionary (PASS=0, q10=1,
// DP/AF/TAG/H2 = 2..5, GT=6, DP_fmt=7).
func TestWriterWriteRecord_FormatRoundtrip(t *testing.T) {
	// 2 samples, diploid GT: s1=0/0 → [2,2]; s2=0|1 → [2,5] (phased on
	// the second allele). FORMAT/DP scalar: [10, 20]. Per-sample dim
	// for GT is 2, for DP is 1; descriptor bytes 0x21 (int8, size=2)
	// and 0x11 (int8, size=1) respectively.
	indiv := []byte{}
	indiv = append(indiv, EncodeTypedInt8(6)...) // FMT key = GT (IDX 6)
	indiv = append(indiv, 0x21, 2, 2, 2, 5)      // descriptor + payload
	indiv = append(indiv, EncodeTypedInt8(7)...) // FMT key = DP (IDX 7)
	indiv = append(indiv, 0x11, 10, 20)          // descriptor + payload

	rec := buildRecord(recordSpec{
		chromID:    0,
		pos:        149,
		rlen:       1,
		qual:       40,
		id:         "rsY",
		alleles:    []string{"A", "T"},
		nSample:    2,
		nFmt:       2,
		indivBytes: indiv,
	})
	stream := buildBCFStream(t, rec)
	r, err := NewReader(bytes.NewReader(stream))
	if err != nil {
		t.Fatal(err)
	}
	recs, err := r.ReadAll()
	if err != nil || len(recs) != 1 {
		t.Fatalf("first read: %v / %d", err, len(recs))
	}
	if got := recs[0].FmtVals; len(got) != 2 || got[0].Length != 2 || len(got[0].Ints) != 4 {
		t.Fatalf("first decode shape: %+v", got)
	}

	// Write the decoded record back out via WriteRecord and re-read.
	// If encodeTypedValue truncates the per-sample payload, the second
	// read will either fail (too few bytes) or yield mangled GT values.
	var buf bytes.Buffer
	w := NewWriter(&buf, r.Header())
	if err := w.WriteRecord(recs[0]); err != nil {
		t.Fatalf("WriteRecord: %v", err)
	}
	if err := w.Flush(); err != nil {
		t.Fatal(err)
	}
	r2, err := NewReader(&buf)
	if err != nil {
		t.Fatal(err)
	}
	out, err := r2.ReadAll()
	if err != nil {
		t.Fatalf("re-read: %v", err)
	}
	if len(out) != 1 {
		t.Fatalf("re-read count: %d", len(out))
	}
	if len(out[0].FmtVals) != 2 {
		t.Fatalf("re-read FmtVals: %d", len(out[0].FmtVals))
	}
	gt := out[0].FmtVals[0]
	if gt.Length != 2 || len(gt.Ints) != 4 {
		t.Errorf("GT shape after roundtrip: Length=%d Ints=%v (want Length=2, len(Ints)=4)", gt.Length, gt.Ints)
	} else {
		want := []int32{2, 2, 2, 5}
		for i, v := range want {
			if gt.Ints[i] != v {
				t.Errorf("GT[%d]: got %d want %d", i, gt.Ints[i], v)
			}
		}
	}
	dp := out[0].FmtVals[1]
	if dp.Length != 1 || len(dp.Ints) != 2 {
		t.Errorf("DP shape after roundtrip: Length=%d Ints=%v (want Length=1, len(Ints)=2)", dp.Length, dp.Ints)
	} else if dp.Ints[0] != 10 || dp.Ints[1] != 20 {
		t.Errorf("DP values: %v want [10 20]", dp.Ints)
	}
}

// TestNewWriterFromVCFHeader exercises the convenience constructor: build a
// vcf.Header, emit a record, parse it back.
func TestNewWriterFromVCFHeader(t *testing.T) {
	vh := &vcf.Header{
		MetaInfo: []string{
			"##fileformat=VCFv4.2",
			"##contig=<ID=chr1>",
			"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">",
		},
		Samples: nil,
	}
	var buf bytes.Buffer
	w, err := NewWriterFromVCFHeader(&buf, vh)
	if err != nil {
		t.Fatalf("NewWriterFromVCFHeader: %v", err)
	}
	if err := w.Write(&vcf.Variant{
		Chrom: "chr1", Pos: 5, ID: ".", Ref: "A", Alt: []string{"T"},
		Qual: 12, Filter: []string{"PASS"}, Info: map[string]string{"DP": "8"},
	}); err != nil {
		t.Fatalf("Write: %v", err)
	}
	if err := w.Flush(); err != nil {
		t.Fatal(err)
	}
	if buf.Len() < 10 {
		t.Fatalf("output too short: %d bytes", buf.Len())
	}
	// Parse it back.
	r, err := NewReader(&buf)
	if err != nil {
		t.Fatal(err)
	}
	recs, err := r.ReadAll()
	if err != nil || len(recs) != 1 {
		t.Fatalf("ReadAll: recs=%d err=%v", len(recs), err)
	}
	got := recs[0].ToVariant(r.Header())
	if got.Pos != 5 {
		t.Errorf("Pos: %d", got.Pos)
	}
}

// TestWriterNilHeader checks the nil-header error path on the convenience
// constructor.
func TestNewWriterFromVCFHeaderNil(t *testing.T) {
	var buf bytes.Buffer
	if _, err := NewWriterFromVCFHeader(&buf, nil); err == nil {
		t.Fatal("expected error for nil vcf.Header")
	}
}

// TestParseGT covers the GT text parser in detail.
func TestParseGT(t *testing.T) {
	cases := []struct {
		in   string
		want []int32
	}{
		{"0/0", []int32{2, 2}},
		{"0/1", []int32{2, 4}},
		{"1|1", []int32{4, 5}},
		{".", []int32{MissingInt32}},
		{"", []int32{MissingInt32}},
		{"./.", []int32{MissingInt32, MissingInt32}},
		{"1", []int32{4}},
	}
	for _, c := range cases {
		got := parseGT(c.in)
		if !int32SlicesEqual(got, c.want) {
			t.Errorf("parseGT(%q) = %v, want %v", c.in, got, c.want)
		}
	}
}

func int32SlicesEqual(a, b []int32) bool {
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

// TestPickIntWidth picks the narrowest width that fits every value.
func TestPickIntWidth(t *testing.T) {
	cases := []struct {
		vs   []int32
		want int
	}{
		{[]int32{0, 1, 2}, 1},
		{[]int32{127, -120}, 1},
		{[]int32{128}, 2},
		{[]int32{32767}, 2},
		{[]int32{32768}, 4},
		{[]int32{MissingInt32}, 1},
		{[]int32{EndOfVectorInt32}, 1},
		{[]int32{-129}, 2},
	}
	for _, c := range cases {
		if got := pickIntWidth(c.vs); got != c.want {
			t.Errorf("pickIntWidth(%v) = %d, want %d", c.vs, got, c.want)
		}
	}
}

// TestEncodeInfoValue covers each INFO type branch.
func TestEncodeInfoValue(t *testing.T) {
	intEntry := DictEntry{ID: "DP", Type: "Integer", Number: "1"}
	if got := encodeInfoValue(intEntry, "42"); len(got) == 0 {
		t.Fatal("integer encoding empty")
	}
	floatEntry := DictEntry{ID: "AF", Type: "Float", Number: "A"}
	if got := encodeInfoValue(floatEntry, "0.5"); len(got) == 0 {
		t.Fatal("float encoding empty")
	}
	flagEntry := DictEntry{ID: "H2", Type: "Flag", Number: "0"}
	got := encodeInfoValue(flagEntry, "")
	if len(got) != 2 || got[0] != 0x11 {
		t.Errorf("flag encoding: %v", got)
	}
	stringEntry := DictEntry{ID: "TAG", Type: "String", Number: "1"}
	if got := encodeInfoValue(stringEntry, "hello"); len(got) == 0 {
		t.Fatal("string encoding empty")
	}
	if got := encodeInfoValue(stringEntry, "."); !bytes.Equal(got, EncodeMissing()) {
		t.Errorf("dot string should be missing: %v", got)
	}
}

// TestEncodeIntsFromText exercises the integer text parser, including missing
// elements and non-numeric strings that fall back to char encoding.
func TestEncodeIntsFromText(t *testing.T) {
	if got := encodeIntsFromText("."); !bytes.Equal(got, EncodeMissing()) {
		t.Errorf("dot int: got %v", got)
	}
	if got := encodeIntsFromText("1,2,3"); len(got) == 0 {
		t.Fatal("multi int empty")
	}
	// Non-numeric falls back to string.
	got := encodeIntsFromText("abc")
	if len(got) == 0 || (got[0]&0x0F) != TypeChar {
		t.Errorf("non-numeric int should fall back to char: %v", got)
	}
	// "." inside the list becomes the missing sentinel.
	if got := encodeIntsFromText("1,.,2"); len(got) == 0 {
		t.Fatal("missing-in-list empty")
	}
}

// TestEncodeFloatsFromText exercises the float text parser similarly.
func TestEncodeFloatsFromText(t *testing.T) {
	if got := encodeFloatsFromText("."); !bytes.Equal(got, EncodeMissing()) {
		t.Errorf("dot float: %v", got)
	}
	if got := encodeFloatsFromText("0.1,0.2"); len(got) == 0 {
		t.Fatal("vec float empty")
	}
	got := encodeFloatsFromText("not-a-number")
	if len(got) == 0 || (got[0]&0x0F) != TypeChar {
		t.Errorf("non-numeric float should fall back to char: %v", got)
	}
}

// TestWriteHeaderTwice asserts WriteHeader is idempotent.
func TestWriteHeaderTwice(t *testing.T) {
	hdr := makeHeader(t)
	var buf bytes.Buffer
	w := NewWriter(&buf, hdr)
	if err := w.WriteHeader(); err != nil {
		t.Fatal(err)
	}
	n1 := buf.Len()
	if err := w.WriteHeader(); err != nil {
		t.Fatal(err)
	}
	if buf.Len() != n1 {
		t.Fatalf("WriteHeader is not idempotent: %d vs %d", n1, buf.Len())
	}
}

// TestBuildBCFTextHeaderSynthesizesFileformat ensures synthesized headers
// always include a fileformat line.
func TestBuildBCFTextHeaderSynthesizesFileformat(t *testing.T) {
	vh := &vcf.Header{MetaInfo: []string{"##contig=<ID=x>"}}
	got := buildBCFTextHeader(vh)
	if !strings.Contains(got, "##fileformat=") {
		t.Errorf("missing fileformat line in: %q", got)
	}
}

// TestEncodeMissingFloatBitPattern guards the exact float missing value the
// reader expects to see on the wire.
func TestEncodeMissingFloatBitPattern(t *testing.T) {
	want := uint32(0x7F800001)
	got := math.Float32bits(math.Float32frombits(MissingFloat32))
	if got != want {
		t.Fatalf("missing float bits got %#x want %#x", got, want)
	}
}

// TestEncodeFormatFieldsEncodesEmpty covers the no-samples path.
func TestEncodeFormatFieldsEncodesEmpty(t *testing.T) {
	entry := DictEntry{ID: "DP", Type: "Integer", Number: "1"}
	got, err := encodeFormatField(entry, "DP", nil, 0)
	if err != nil {
		t.Fatal(err)
	}
	// Zero samples → missing.
	if !bytes.Equal(got, EncodeMissing()) {
		t.Errorf("expected missing for empty sample list, got %v", got)
	}
}
