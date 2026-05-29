package bcf

import (
	"bytes"
	"compress/gzip"
	"errors"
	"io"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// buildOracleBCF reads testdata/oracle/basic.vcf and writes a BCF stream
// through the exact path bcftools view -Ou/-Ob uses:
// NewWriterFromVCFHeader -> buildBCFTextHeader -> parseTextHeader.
func buildOracleBCF(t *testing.T) (*Header, []byte) {
	t.Helper()
	f, err := os.Open(filepath.Join("testdata", "oracle", "basic.vcf"))
	if err != nil {
		t.Fatalf("open basic.vcf: %v", err)
	}
	defer f.Close()

	vr := vcf.NewReader(f)
	vh, err := vr.ReadHeader()
	if err != nil {
		t.Fatalf("vcf ReadHeader: %v", err)
	}
	variants, err := vr.ReadAll()
	if err != nil {
		t.Fatalf("vcf ReadAll: %v", err)
	}

	var buf bytes.Buffer
	w, err := NewWriterFromVCFHeader(&buf, vh)
	if err != nil {
		t.Fatalf("NewWriterFromVCFHeader: %v", err)
	}
	for _, v := range variants {
		if err := w.Write(v); err != nil {
			t.Fatalf("Write %s:%d: %v", v.Chrom, v.Pos, err)
		}
	}
	if err := w.Flush(); err != nil {
		t.Fatalf("Flush: %v", err)
	}
	return w.Header(), buf.Bytes()
}

// TestBCFWriterDictParity verifies our BCF writer reproduces htslib's unified
// INFO/FILTER/FORMAT dictionary numbering: a FORMAT tag whose ID already
// appeared as an INFO tag reuses that IDX (FORMAT/DP reuses INFO/DP's IDX),
// contigs use a separate dict numbered from 0, and all records round-trip.
func TestBCFWriterDictParity(t *testing.T) {
	hdr, _ := buildOracleBCF(t)

	// Expected unified INFO/FILTER/FORMAT IDX in first-appearance order.
	wantTag := map[string]int32{
		"PASS":  0,
		"DP":    1, // INFO/DP; FORMAT/DP must reuse this
		"AF":    2,
		"AC":    3,
		"AN":    4,
		"INDEL": 5,
		"q10":   6,
		"lowDP": 7,
		"GT":    8,
		"GQ":    9,
	}

	got := map[string]int32{}
	for _, e := range hdr.InfoTags {
		got[e.ID] = e.IDX
	}
	for _, e := range hdr.FmtTags {
		if prev, ok := got[e.ID]; ok && prev != e.IDX {
			t.Errorf("tag %q has conflicting IDX: INFO/FILTER=%d FORMAT=%d", e.ID, prev, e.IDX)
		}
		got[e.ID] = e.IDX
	}
	for name, idx := range wantTag {
		if got[name] != idx {
			t.Errorf("tag %q IDX = %d, want %d", name, got[name], idx)
		}
	}

	// The shared-name invariant: FORMAT/DP and INFO/DP carry the same IDX.
	var infoDP, fmtDP int32 = -1, -1
	for _, e := range hdr.InfoTags {
		if e.ID == "DP" {
			infoDP = e.IDX
		}
	}
	for _, e := range hdr.FmtTags {
		if e.ID == "DP" {
			fmtDP = e.IDX
		}
	}
	if infoDP != 1 || fmtDP != 1 {
		t.Errorf("INFO/DP IDX=%d FORMAT/DP IDX=%d, both want 1", infoDP, fmtDP)
	}

	// Contigs use their own dict starting at 0.
	if len(hdr.Contigs) != 2 || hdr.Contigs[0].ID != "chr1" || hdr.Contigs[0].IDX != 0 ||
		hdr.Contigs[1].ID != "chr2" || hdr.Contigs[1].IDX != 1 {
		t.Errorf("contigs = %+v, want chr1=0 chr2=1", hdr.Contigs)
	}
}

// TestBCFWriterDictRoundTrip writes the oracle VCF as BCF and reads it back
// with our own reader, asserting all 6 records survive.
func TestBCFWriterDictRoundTrip(t *testing.T) {
	_, stream := buildOracleBCF(t)

	r, err := NewReader(bytes.NewReader(stream))
	if err != nil {
		t.Fatalf("NewReader: %v", err)
	}
	recs, err := r.ReadAll()
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if len(recs) != 6 {
		t.Fatalf("round-trip records = %d, want 6", len(recs))
	}
	// Spot-check the first record's FORMAT/DP decodes (it shares INFO/DP's IDX).
	v := recs[0].ToVariant(r.Header())
	if len(v.Samples) != 3 {
		t.Fatalf("record 0 samples = %d, want 3", len(v.Samples))
	}
	if v.Samples[0].Data["DP"] != "10" {
		t.Errorf("record 0 sample 0 FORMAT/DP = %q, want 10", v.Samples[0].Data["DP"])
	}
}

// TestBCFWriterBodyByteIdentity asserts the BCF record body our writer emits
// for basic.vcf is byte-for-byte identical to the genuine `bcftools view -Ou`
// oracle (testdata/oracle/basic.bcf). We compare only the record body — the
// text header differs because we intentionally omit the ##bcftools provenance
// lines htslib appends.
func TestBCFWriterBodyByteIdentity(t *testing.T) {
	oracle, err := os.ReadFile(filepath.Join("testdata", "oracle", "basic.bcf"))
	if err != nil {
		t.Fatalf("read oracle: %v", err)
	}
	oracleRaw, err := bgzfDecode(oracle)
	if err != nil {
		t.Fatalf("bgzf-decode oracle: %v", err)
	}
	_, oracleBody, err := splitBCF(oracleRaw)
	if err != nil {
		t.Fatalf("split oracle: %v", err)
	}

	_, ourStream := buildOracleBCF(t)
	// buildOracleBCF writes a raw (uncompressed) BCF stream.
	_, ourBody, err := splitBCF(ourStream)
	if err != nil {
		t.Fatalf("split ours: %v", err)
	}

	if !bytes.Equal(oracleBody, ourBody) {
		t.Fatalf("record body differs from oracle: oracle %d bytes, ours %d bytes",
			len(oracleBody), len(ourBody))
	}
}

// splitBCF returns (textHeader, recordBody) from a raw (decompressed) BCF
// stream: 5-byte magic, uint32 l_text, l_text bytes of text, then records.
func splitBCF(raw []byte) (text, body []byte, err error) {
	if len(raw) < 9 || string(raw[:5]) != "BCF\x02\x02" {
		return nil, nil, errors.New("bad BCF magic")
	}
	lText := uint32(raw[5]) | uint32(raw[6])<<8 | uint32(raw[7])<<16 | uint32(raw[8])<<24
	if int(9+lText) > len(raw) {
		return nil, nil, errors.New("truncated BCF text")
	}
	return raw[9 : 9+lText], raw[9+lText:], nil
}

// bgzfDecode inflates a BGZF (gzip-compatible) stream. The genuine bcftools
// oracle is BGZF-wrapped; our writer emits a raw stream, so we only need this
// for the oracle.
func bgzfDecode(data []byte) ([]byte, error) {
	zr, err := gzip.NewReader(bytes.NewReader(data))
	if err != nil {
		return nil, err
	}
	zr.Multistream(true)
	defer zr.Close()
	return io.ReadAll(zr)
}

// TestBCFWriterGenuineBcftoolsReads shells out to the genuine bcftools binary
// (if present) and asserts it reads all 6 records from our BCF output without
// the "Invalid FORMAT id" error that the dictionary bug produced.
func TestBCFWriterGenuineBcftoolsReads(t *testing.T) {
	bin := filepath.Join("..", "..", "..", "reference_code", "bcftools", "bcftools")
	if _, err := os.Stat(bin); err != nil {
		t.Skipf("genuine bcftools not built at %s: %v", bin, err)
	}

	_, stream := buildOracleBCF(t)

	// bcftools wants a BGZF-wrapped BCF on disk to sniff the file type, but it
	// also reads a raw (uncompressed) BCF magic from a regular file. Write the
	// raw stream and let bcftools detect it.
	tmp := filepath.Join(t.TempDir(), "our.bcf")
	if err := os.WriteFile(tmp, stream, 0o644); err != nil {
		t.Fatalf("write tmp bcf: %v", err)
	}

	cmd := exec.Command(bin, "view", "-H", tmp)
	out, err := cmd.CombinedOutput()
	if err != nil {
		t.Fatalf("genuine bcftools view failed: %v\noutput:\n%s", err, out)
	}
	if strings.Contains(string(out), "Invalid FORMAT id") ||
		strings.Contains(string(out), "Bad BCF record") {
		t.Fatalf("genuine bcftools rejected our BCF:\n%s", out)
	}
	lines := 0
	for _, ln := range strings.Split(strings.TrimRight(string(out), "\n"), "\n") {
		if ln != "" && !strings.HasPrefix(ln, "#") {
			lines++
		}
	}
	if lines != 6 {
		t.Fatalf("genuine bcftools read %d records, want 6\noutput:\n%s", lines, out)
	}
}
