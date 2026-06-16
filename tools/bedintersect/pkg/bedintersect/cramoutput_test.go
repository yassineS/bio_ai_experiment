package bedintersect

import (
	"bytes"
	"io"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/cram"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// buildCRAM serialises recs (with hdr) into an in-memory reference-free CRAM
// byte stream, for driving the CRAM-output path in tests.
func buildCRAM(t *testing.T, hdr *sam.Header, recs []*sam.Record) []byte {
	t.Helper()
	var buf bytes.Buffer
	if err := cram.WriteCRAM(&buf, hdr, recs); err != nil {
		t.Fatalf("WriteCRAM: %v", err)
	}
	return buf.Bytes()
}

// decodeCRAMNames decodes a CRAM byte stream and returns its records' QNAMEs in
// order, the simplest surface for asserting which alignments survived.
func decodeCRAMNames(t *testing.T, raw []byte) []string {
	t.Helper()
	rr, err := cram.NewRecordReader(bytes.NewReader(raw))
	if err != nil {
		t.Fatalf("NewRecordReader: %v", err)
	}
	var names []string
	for {
		rec, err := rr.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("Read: %v", err)
		}
		names = append(names, rec.QName)
	}
	return names
}

// TestIntersectBinaryOutput_CRAM drives IntersectBinaryOutput with OutputCRAM
// over the alignment-level output modes and asserts both that the output is
// CRAM-framed and that the right alignments survive: default/-wa/-u keep each A
// with >=1 overlap (once), -v keeps each A with no overlap.
func TestIntersectBinaryOutput_CRAM(t *testing.T) {
	hdr := bamTestHeader(t)
	recs := []*sam.Record{
		bamTestRecord("r1", 0, 110, "20M"),  // chr1:109-129 -> overlaps 100-200
		bamTestRecord("r2", 0, 500, "20M"),  // chr1:499-519 -> no overlap
		bamTestRecord("r3", 16, 150, "20M"), // chr1:149-169 -> overlaps 100-200
	}
	cramBytes := buildCRAM(t, hdr, recs)
	bBED := "chr1\t100\t200\tb1\t0\t+\n"

	tests := []struct {
		name string
		opts IntersectOptions
		want []string
	}{
		{"default", IntersectOptions{MinOverlap: 1}, []string{"r1", "r3"}},
		{"writeA", IntersectOptions{MinOverlap: 1, WriteA: true}, []string{"r1", "r3"}},
		{"unique", IntersectOptions{MinOverlap: 1, Unique: true}, []string{"r1", "r3"}},
		{"invert", IntersectOptions{MinOverlap: 1, NoOverlap: true}, []string{"r2"}},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			var out bytes.Buffer
			n, err := IntersectBinaryOutput(bytes.NewReader(cramBytes),
				[]io.Reader{strings.NewReader(bBED)}, &out, tc.opts,
				AlnOutputOptions{Format: OutputCRAM})
			if err != nil {
				t.Fatalf("IntersectBinaryOutput(CRAM): %v", err)
			}
			raw := out.Bytes()
			if len(raw) < 4 || string(raw[:4]) != "CRAM" {
				t.Fatalf("output is not CRAM-framed (first bytes %q)", firstN(raw, 4))
			}
			got := decodeCRAMNames(t, raw)
			if n != len(tc.want) {
				t.Fatalf("count=%d, want %d", n, len(tc.want))
			}
			if strings.Join(got, ",") != strings.Join(tc.want, ",") {
				t.Fatalf("survivors=%v, want %v", got, tc.want)
			}
		})
	}
}

// TestIntersectBinaryOutput_BAMUnchanged confirms OutputBAM (and the
// IntersectBAMOutput wrapper) still produce BGZF BAM, so the CRAM addition does
// not regress the BAM path.
func TestIntersectBinaryOutput_BAMUnchanged(t *testing.T) {
	hdr := bamTestHeader(t)
	recs := []*sam.Record{bamTestRecord("r1", 0, 110, "20M")}
	cramBytes := buildCRAM(t, hdr, recs)
	bBED := "chr1\t100\t200\tb1\t0\t+\n"

	var out bytes.Buffer
	if _, err := IntersectBinaryOutput(bytes.NewReader(cramBytes),
		[]io.Reader{strings.NewReader(bBED)}, &out, IntersectOptions{MinOverlap: 1},
		AlnOutputOptions{Format: OutputBAM}); err != nil {
		t.Fatalf("IntersectBinaryOutput(BAM): %v", err)
	}
	raw := out.Bytes()
	if len(raw) < 2 || raw[0] != 0x1f || raw[1] != 0x8b {
		t.Fatalf("OutputBAM did not produce a BGZF stream (first bytes %v)", firstN(raw, 2))
	}
	if got := decodeBAMNames(t, raw); strings.Join(got, ",") != "r1" {
		t.Fatalf("BAM survivors=%v, want [r1]", got)
	}
}

// firstN returns up to the first n bytes of b, for diagnostic messages.
func firstN(b []byte, n int) []byte {
	if len(b) < n {
		return b
	}
	return b[:n]
}

// TestClassifyQueryInput distinguishes CRAM from BAM (the distinction that lets
// a CRAM query be re-emitted as CRAM), while text inputs stay QueryText.
func TestClassifyQueryInput(t *testing.T) {
	hdr := bamTestHeader(t)
	bamBytes := buildBAM(t, hdr, []*sam.Record{bamTestRecord("r1", 0, 100, "10M")})
	cramBytes := buildCRAM(t, hdr, []*sam.Record{bamTestRecord("r1", 0, 100, "10M")})

	tests := []struct {
		name string
		data []byte
		want QueryFormat
	}{
		{"bgzf_bam", bamBytes, QueryBAM},
		{"raw_bam_magic", []byte("BAM\x01rest..."), QueryBAM},
		{"cram", cramBytes, QueryCRAM},
		{"cram_magic", []byte("CRAM\x03..."), QueryCRAM},
		{"plain_bed", []byte("chr1\t10\t20\n"), QueryText},
		{"sam_text", []byte("@HD\tVN:1.6\n"), QueryText},
		{"empty", nil, QueryText},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			format, r, err := ClassifyQueryInput(bytes.NewReader(tc.data))
			if err != nil {
				t.Fatalf("ClassifyQueryInput: %v", err)
			}
			if format != tc.want {
				t.Fatalf("format=%d, want %d", format, tc.want)
			}
			replayed, err := io.ReadAll(r)
			if err != nil {
				t.Fatalf("ReadAll(replayed): %v", err)
			}
			if !bytes.Equal(replayed, tc.data) {
				t.Fatalf("replayed %d bytes, want %d", len(replayed), len(tc.data))
			}
		})
	}
}
