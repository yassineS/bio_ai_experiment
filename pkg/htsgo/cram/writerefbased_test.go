package cram

import (
	"bytes"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// TestReferenceBasedEncoding exercises the writer's reference-based path: when a
// reference is supplied, a mapped read is stored as substitution features at its
// mismatches (matched bases reconstructed from the reference on decode) rather
// than carrying the whole sequence. It checks (a) the reads round-trip exactly
// — SEQ/CIGAR/POS/FLAG — through our decoder fed the same reference, and (b) the
// reference-based encoding is materially smaller than the reference-free
// encoding of the same records, proving matched bases are not stored.
func TestReferenceBasedEncoding(t *testing.T) {
	h := writerTestHeader() // chr1 LN:100000

	// A 400 bp reference with base variety so substitutions exercise several
	// rows of the substitution matrix.
	ref := []byte(strings.Repeat("ACGTACGTGGCCATGCTAGC", 20)) // 400 bp
	refMap := map[string][]byte{"chr1": ref}

	// mut copies an aligned window of the reference and applies in-read edits,
	// producing a read that aligns at pos1 (1-based) with the given mismatches.
	mut := func(pos1, n int, edits map[int]byte) string {
		b := append([]byte(nil), ref[pos1-1:pos1-1+n]...)
		for off, base := range edits {
			b[off] = base
		}
		return string(b)
	}

	in := []*sam.Record{
		mkRec("r_exact", "chr1", 50, "30M", mut(50, 30, nil)),                                                         // no mismatch
		mkRec("r_subs", "chr1", 100, "40M", mut(100, 40, map[int]byte{5: 'T', 22: 'A', 31: 'C'})),                     // 3 substitutions
		mkRec("r_indel", "chr1", 200, "20M5D20M", mut(200, 45, map[int]byte{3: 'G', 38: 'T'})[:20]+mut(225, 20, nil)), // deletion + subs
	}

	// Reference-based encode.
	var refBuf bytes.Buffer
	rw, err := NewRecordWriterOpts(&refBuf, h, WriterOptions{Reference: refMap})
	if err != nil {
		t.Fatalf("NewRecordWriterOpts(reference): %v", err)
	}
	for _, rec := range in {
		if err := rw.Write(rec); err != nil {
			t.Fatalf("Write: %v", err)
		}
	}
	if err := rw.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	// Reference-free encode of the same records, for the size comparison.
	var freeBuf bytes.Buffer
	if err := WriteCRAM(&freeBuf, h, in); err != nil {
		t.Fatalf("WriteCRAM(reference-free): %v", err)
	}
	if refBuf.Len() >= freeBuf.Len() {
		t.Errorf("reference-based encoding (%d B) not smaller than reference-free (%d B); matched bases appear still stored",
			refBuf.Len(), freeBuf.Len())
	}

	// Decode the reference-based CRAM, feeding the same reference.
	rr, err := NewRecordReader(bytes.NewReader(refBuf.Bytes()))
	if err != nil {
		t.Fatalf("NewRecordReader: %v", err)
	}
	rr.SetReference(&stubReference{name: "chr1", seq: string(ref)})
	out, err := rr.ReadAll()
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if len(out) != len(in) {
		t.Fatalf("decoded %d records, want %d", len(out), len(in))
	}
	for i := range in {
		if normSeq(out[i].Seq) != normSeq(in[i].Seq) {
			t.Errorf("record %d SEQ mismatch:\n got %s\nwant %s", i, out[i].Seq, in[i].Seq)
		}
		if out[i].Cigar.String() != in[i].Cigar.String() {
			t.Errorf("record %d CIGAR = %s, want %s", i, out[i].Cigar.String(), in[i].Cigar.String())
		}
		if out[i].Pos != in[i].Pos {
			t.Errorf("record %d POS = %d, want %d", i, out[i].Pos, in[i].Pos)
		}
		if out[i].Flag != in[i].Flag {
			t.Errorf("record %d FLAG = %d, want %d", i, out[i].Flag, in[i].Flag)
		}
	}
}
