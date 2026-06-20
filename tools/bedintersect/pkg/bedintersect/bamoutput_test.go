package bedintersect

import (
	"bytes"
	"compress/gzip"
	"io"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// buildBAM serialises recs (with the given @SQ header) into an in-memory
// BGZF-wrapped BAM byte stream, for driving the BAM-output path in tests.
func buildBAM(t *testing.T, hdr *sam.Header, recs []*sam.Record) []byte {
	t.Helper()
	var buf bytes.Buffer
	w := sam.NewBAMWriter(&buf)
	if err := w.WriteHeader(hdr); err != nil {
		t.Fatalf("WriteHeader: %v", err)
	}
	for _, r := range recs {
		if err := w.Write(r); err != nil {
			t.Fatalf("Write: %v", err)
		}
	}
	if err := w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	return buf.Bytes()
}

// decodeBAMNames decodes a BAM byte stream and returns the QNAMEs of its records
// in order, the simplest surface for asserting which alignments survived.
func decodeBAMNames(t *testing.T, raw []byte) []string {
	t.Helper()
	r, err := sam.NewBAMReader(bytes.NewReader(raw))
	if err != nil {
		t.Fatalf("NewBAMReader: %v", err)
	}
	var names []string
	for {
		rec, err := r.Read()
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

func bamTestHeader(t *testing.T) *sam.Header {
	t.Helper()
	hr, err := sam.NewSAMReader(strings.NewReader("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:10000\n"))
	if err != nil {
		t.Fatalf("NewSAMReader(header): %v", err)
	}
	return hr.Header()
}

func bamTestRecord(qname string, flag uint16, pos int32, cig string) *sam.Record {
	c, _ := sam.ParseCigar(cig)
	seqLen := c.QueryLength()
	return &sam.Record{
		QName: qname, Flag: flag, RName: "chr1", Pos: int64(pos), MapQ: 60,
		Cigar: c, Seq: strings.Repeat("A", seqLen),
		Qual: bytes.Repeat([]byte{30}, seqLen), RNext: "*", PNext: 0,
	}
}

// TestIntersectBAMOutput_Modes drives IntersectBAMOutput over the alignment-level
// output modes and asserts which alignments survive: default/-wa/-u keep each A
// with >=1 overlap (once), -v keeps each A with no overlap. The BAM header is
// preserved and re-emitted.
func TestIntersectBAMOutput_Modes(t *testing.T) {
	hdr := bamTestHeader(t)
	// r1 overlaps b (chr1:100-200); r2 does not; r3 overlaps b.
	recs := []*sam.Record{
		bamTestRecord("r1", 0, 110, "20M"),  // chr1:109-129 -> overlaps 100-200
		bamTestRecord("r2", 0, 500, "20M"),  // chr1:499-519 -> no overlap
		bamTestRecord("r3", 16, 150, "20M"), // chr1:149-169 -> overlaps 100-200
	}
	bamBytes := buildBAM(t, hdr, recs)
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
		// -wb / -loj ignored: fall through to default selection.
		{"writeB_ignored", IntersectOptions{MinOverlap: 1, WriteB: true}, []string{"r1", "r3"}},
		{"leftJoin_ignored", IntersectOptions{MinOverlap: 1, LeftJoin: true}, []string{"r1", "r3"}},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			var out bytes.Buffer
			n, err := IntersectBAMOutput(bytes.NewReader(bamBytes),
				[]io.Reader{strings.NewReader(bBED)}, &out, tc.opts)
			if err != nil {
				t.Fatalf("IntersectBAMOutput: %v", err)
			}
			got := decodeBAMNames(t, out.Bytes())
			if n != len(tc.want) {
				t.Fatalf("count=%d, want %d", n, len(tc.want))
			}
			if strings.Join(got, ",") != strings.Join(tc.want, ",") {
				t.Fatalf("survivors=%v, want %v", got, tc.want)
			}
		})
	}
}

// TestIntersectBAMOutput_Unmapped verifies an unmapped alignment is absent under
// the default mode (it never overlaps) and reported under -v, matching upstream.
func TestIntersectBAMOutput_Unmapped(t *testing.T) {
	hdr := bamTestHeader(t)
	mapped := bamTestRecord("m1", 0, 110, "20M") // overlaps b
	unmapped := &sam.Record{QName: "u1", Flag: sam.FlagUnmapped, RName: "*", Pos: 0,
		Cigar: nil, Seq: "ACGT", Qual: bytes.Repeat([]byte{30}, 4), RNext: "*"}
	bamBytes := buildBAM(t, hdr, []*sam.Record{mapped, unmapped})
	bBED := "chr1\t100\t200\tb1\t0\t+\n"

	t.Run("default", func(t *testing.T) {
		var out bytes.Buffer
		if _, err := IntersectBAMOutput(bytes.NewReader(bamBytes),
			[]io.Reader{strings.NewReader(bBED)}, &out, IntersectOptions{MinOverlap: 1}); err != nil {
			t.Fatalf("IntersectBAMOutput: %v", err)
		}
		if got := decodeBAMNames(t, out.Bytes()); strings.Join(got, ",") != "m1" {
			t.Fatalf("default survivors=%v, want [m1]", got)
		}
	})
	t.Run("invert", func(t *testing.T) {
		var out bytes.Buffer
		if _, err := IntersectBAMOutput(bytes.NewReader(bamBytes),
			[]io.Reader{strings.NewReader(bBED)}, &out, IntersectOptions{MinOverlap: 1, NoOverlap: true}); err != nil {
			t.Fatalf("IntersectBAMOutput: %v", err)
		}
		if got := decodeBAMNames(t, out.Bytes()); strings.Join(got, ",") != "u1" {
			t.Fatalf("invert survivors=%v, want [u1]", got)
		}
	})
}

// TestIsBAMOrCRAMInput classifies BAM, BGZF-wrapped BAM, CRAM, and several text
// streams (including a BGZF-compressed BED) to confirm only true BAM/CRAM
// streams trigger the BAM-output path and the returned reader replays the bytes.
func TestIsBAMOrCRAMInput(t *testing.T) {
	hdr := bamTestHeader(t)
	bamBytes := buildBAM(t, hdr, []*sam.Record{bamTestRecord("r1", 0, 100, "10M")})

	tests := []struct {
		name string
		data []byte
		want bool
	}{
		{"bgzf_bam", bamBytes, true},
		{"raw_bam_magic", []byte("BAM\x01rest..."), true},
		{"cram_magic", []byte("CRAM\x03..."), true},
		{"plain_bed", []byte("chr1\t10\t20\n"), false},
		{"sam_text", []byte("@HD\tVN:1.6\n"), false},
		{"empty", nil, false},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			isBAM, r, err := IsBAMOrCRAMInput(bytes.NewReader(tc.data))
			if err != nil {
				t.Fatalf("IsBAMOrCRAMInput: %v", err)
			}
			if isBAM != tc.want {
				t.Fatalf("isBAM=%v, want %v", isBAM, tc.want)
			}
			// The returned reader must replay every original byte.
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

// TestIsBAMOrCRAMInput_BGZFText confirms a BGZF-compressed BED (the `.bed.gz`
// piped to `-a -` case) is NOT classified as BAM, so it produces text output.
func TestIsBAMOrCRAMInput_BGZFText(t *testing.T) {
	// A BGZF-wrapped text payload: reuse the sam BAMWriter's BGZF backend is
	// overkill; just gzip the text, which shares the 1f 8b magic and exercises
	// the "decompress and check for BAM magic" branch.
	var gzBuf bytes.Buffer
	gw := gzip.NewWriter(&gzBuf)
	if _, err := gw.Write([]byte("chr1\t10\t20\nchr1\t30\t40\n")); err != nil {
		t.Fatalf("gzip write: %v", err)
	}
	if err := gw.Close(); err != nil {
		t.Fatalf("gzip close: %v", err)
	}
	isBAM, r, err := IsBAMOrCRAMInput(bytes.NewReader(gzBuf.Bytes()))
	if err != nil {
		t.Fatalf("IsBAMOrCRAMInput: %v", err)
	}
	if isBAM {
		t.Fatalf("a gzipped BED must not classify as BAM")
	}
	replayed, _ := io.ReadAll(r)
	if !bytes.Equal(replayed, gzBuf.Bytes()) {
		t.Fatalf("replayed bytes differ from input")
	}
}
