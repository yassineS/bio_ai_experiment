package bcf

import (
	"bytes"
	"encoding/binary"
	"strings"
	"testing"
)

// buildBCFStream assembles a synthetic BCF byte stream that can be fed
// directly into ReadHeader / NewReader. It is intentionally minimal — one
// fileformat line, two contigs, a couple of INFO/FILTER/FORMAT lines, and
// a #CHROM line with two samples.
func buildBCFStream(t *testing.T, body []byte) []byte {
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
	full := text + "\x00"
	var buf bytes.Buffer
	buf.Write(Magic[:])
	if err := binary.Write(&buf, binary.LittleEndian, uint32(len(full))); err != nil {
		t.Fatal(err)
	}
	buf.WriteString(full)
	if body != nil {
		buf.Write(body)
	}
	return buf.Bytes()
}

func TestReadHeaderBasic(t *testing.T) {
	stream := buildBCFStream(t, nil)
	hdr, err := ReadHeader(bytes.NewReader(stream))
	if err != nil {
		t.Fatalf("ReadHeader: %v", err)
	}
	if !strings.Contains(hdr.Text, "##fileformat=VCFv4.3") {
		t.Fatalf("text missing fileformat line: %q", hdr.Text)
	}
	if len(hdr.Contigs) != 2 || hdr.Contigs[0].ID != "chr1" || hdr.Contigs[1].ID != "chr2" {
		t.Fatalf("contigs: %+v", hdr.Contigs)
	}
	// InfoTags = PASS + q10 + DP + AF + TAG + H2 = 6 entries
	if len(hdr.InfoTags) != 6 {
		t.Fatalf("info tags count: %d, want 6: %+v", len(hdr.InfoTags), hdr.InfoTags)
	}
	if hdr.InfoTags[0].ID != "PASS" || hdr.InfoTags[1].ID != "q10" {
		t.Fatalf("info tag order: %+v", hdr.InfoTags)
	}
	if len(hdr.FmtTags) != 2 || hdr.FmtTags[0].ID != "GT" || hdr.FmtTags[1].ID != "DP" {
		t.Fatalf("fmt tags: %+v", hdr.FmtTags)
	}
	if !reflectStringSlicesEqual(hdr.Samples, []string{"S1", "S2"}) {
		t.Fatalf("samples: %+v", hdr.Samples)
	}
}

func TestReadHeaderBadMagic(t *testing.T) {
	bad := bytes.NewReader([]byte("BCF\x01\x01....."))
	if _, err := ReadHeader(bad); err == nil {
		t.Fatal("expected ErrBadMagic")
	}
}

func TestReadHeaderTruncated(t *testing.T) {
	if _, err := ReadHeader(bytes.NewReader([]byte("BCF"))); err == nil {
		t.Fatal("expected truncated error")
	}
}

func TestReadHeaderEmptyText(t *testing.T) {
	var buf bytes.Buffer
	buf.Write(Magic[:])
	binary.Write(&buf, binary.LittleEndian, uint32(0))
	if _, err := ReadHeader(&buf); err == nil {
		t.Fatal("expected error for empty header text")
	}
}

func TestParseStructured(t *testing.T) {
	e := parseStructured(`<ID=DP,Number=1,Type=Integer,Description="Read depth, total">`)
	if e.ID != "DP" || e.Number != "1" || e.Type != "Integer" {
		t.Fatalf("got %+v", e)
	}
}

func TestSplitStructuredQuoted(t *testing.T) {
	parts := splitStructured(`ID=DP,Description="a, b, c",Type=Integer`)
	want := []string{`ID=DP`, `Description="a, b, c"`, `Type=Integer`}
	if !reflectStringSlicesEqual(parts, want) {
		t.Fatalf("got %v, want %v", parts, want)
	}
}

func TestContigNameOutOfRange(t *testing.T) {
	h := &Header{Contigs: []DictEntry{{ID: "chr1"}}}
	if got := h.ContigName(-1); got != "" {
		t.Fatalf("got %q", got)
	}
	if got := h.ContigName(99); got != "" {
		t.Fatalf("got %q", got)
	}
	if got := h.ContigName(0); got != "chr1" {
		t.Fatalf("got %q", got)
	}
}

func TestInfoTagOutOfRange(t *testing.T) {
	h := &Header{}
	if h.InfoTag(-1) != nil || h.InfoTag(5) != nil {
		t.Fatal("expected nil for out of range index")
	}
}

func TestFmtTagOutOfRange(t *testing.T) {
	h := &Header{}
	if h.FmtTag(-1) != nil || h.FmtTag(5) != nil {
		t.Fatal("expected nil for out of range index")
	}
}

// TestParseTextHeader_NameDedupAcrossINFOAndFORMAT verifies that an
// `##INFO=<ID=DP>` followed by `##FORMAT=<ID=DP>` (with no explicit
// `,IDX=` annotations, mimicking upstream htslib output) ends up with
// both DP entries sharing the same wire IDX. Without this dedup our
// FmtTag(idx) lookup in record.ToVariant silently drops upstream-
// written FORMAT/DP fields whose wire IDX matches the earlier
// INFO/DP IDX.
func TestParseTextHeader_NameDedupAcrossINFOAndFORMAT(t *testing.T) {
	text := "##fileformat=VCFv4.2\n" +
		"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total depth\">\n" +
		"##INFO=<ID=AF,Number=A,Type=Float,Description=\"AF\">\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n" +
		"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Per-sample depth\">\n" +
		"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"GQ\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n"
	hdr, err := parseTextHeader(text)
	if err != nil {
		t.Fatal(err)
	}
	infoDP := findEntry(hdr.InfoTags, "DP")
	fmtDP := findEntry(hdr.FmtTags, "DP")
	if infoDP == nil || fmtDP == nil {
		t.Fatalf("missing DP entries: info=%v fmt=%v", infoDP, fmtDP)
	}
	if infoDP.IDX != fmtDP.IDX {
		t.Errorf("INFO/DP IDX=%d, FORMAT/DP IDX=%d (want equal — htslib name-dedups)", infoDP.IDX, fmtDP.IDX)
	}
	// PASS=0, INFO/DP=1 (deduped), INFO/AF=2, FORMAT/GT=3, FORMAT/DP=1
	// (reuses INFO/DP), FORMAT/GQ=4.
	if want := int32(1); infoDP.IDX != want {
		t.Errorf("INFO/DP IDX=%d, want %d", infoDP.IDX, want)
	}
	fmtGT := findEntry(hdr.FmtTags, "GT")
	if fmtGT == nil || fmtGT.IDX != 3 {
		t.Errorf("FORMAT/GT IDX=%d, want 3", fmtGT.IDX)
	}
	fmtGQ := findEntry(hdr.FmtTags, "GQ")
	if fmtGQ == nil || fmtGQ.IDX != 4 {
		t.Errorf("FORMAT/GQ IDX=%d, want 4", fmtGQ.IDX)
	}
}

func findEntry(slice []DictEntry, id string) *DictEntry {
	for i := range slice {
		if slice[i].ID == id {
			return &slice[i]
		}
	}
	return nil
}

func reflectStringSlicesEqual(a, b []string) bool {
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
