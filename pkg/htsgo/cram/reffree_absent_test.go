package cram

import (
	"bytes"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// mustCigarRF parses a CIGAR string for the reference-free tests, failing
// the test on a parse error.
func mustCigarRF(t *testing.T, s string) sam.Cigar {
	t.Helper()
	c, err := sam.ParseCigar(s)
	if err != nil {
		t.Fatalf("ParseCigar(%q): %v", s, err)
	}
	return c
}

// absentContigHeader builds a SAM header with a single reference whose
// bases are deliberately NOT in the encode reference map below, modelling
// `samtools view -C -T ref` where the -T reference's contigs do not match
// the BAM (e.g. a chr-prefixed BAM against an unprefixed GRCh37 FASTA).
func absentContigHeader() *sam.Header {
	h, err := sam.ParseHeaderText(
		"@HD\tVN:1.6\tSO:coordinate\n" +
			"@SQ\tSN:chrM\tLN:200\n" +
			"@RG\tID:rg1\tSM:sample1\n")
	if err != nil {
		panic(err)
	}
	return h
}

// TestRefFreeAbsentContigRoundTrips is the regression guard for the
// reference-free-on-absent-contig bug. When a CRAM is encoded with a -T
// reference that does NOT contain the record's contig, the writer encodes
// the slice reference-free (RR=0, bases stored verbatim), exactly as
// upstream htslib does (it warns and proceeds reference-free). The decode
// of such a slice must NOT consult the external reference — the contig is
// absent from it, so a fetch would fail with "contig not in index". This
// test encodes with a reference map that omits the contig, attaches a
// reference source on decode that ERRORS for that contig (modelling the
// fasta "not in index" failure), and asserts the decode succeeds and the
// bases round-trip verbatim.
func TestRefFreeAbsentContigRoundTrips(t *testing.T) {
	h := absentContigHeader()
	records := []*sam.Record{
		{
			QName: "read1", Flag: 0, RName: "chrM", Pos: 1, MapQ: 60,
			Cigar: mustCigarRF(t, "10M"), RNext: "*", PNext: 0, TLen: 0,
			Seq: "ACGTACGTAC", Qual: bytes.Repeat([]byte{'I'}, 10),
		},
		{
			QName: "read2", Flag: 0, RName: "chrM", Pos: 5, MapQ: 60,
			Cigar: mustCigarRF(t, "8M"), RNext: "*", PNext: 0, TLen: 0,
			Seq: "TTTTGGGG", Qual: bytes.Repeat([]byte{'F'}, 8),
		},
	}

	// Encode with a -T reference whose map does NOT contain "chrM" (it has a
	// different, unrelated contig). The writer must fall back to reference-free
	// for chrM and write RR=0; ReferencePath is set so the @SQ UR tag is
	// emitted exactly as `samtools view -C -T` would.
	var buf bytes.Buffer
	rw, err := NewRecordWriterOpts(&buf, h, WriterOptions{
		Reference:     map[string][]byte{"chr1": []byte("ACGTACGTACGTACGT")},
		ReferencePath: "/some/mismatched/ref.fa",
	})
	if err != nil {
		t.Fatalf("NewRecordWriterOpts: %v", err)
	}
	for i, rec := range records {
		if err := rw.Write(rec); err != nil {
			t.Fatalf("Write record %d: %v", i, err)
		}
	}
	if err := rw.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	// Decode with a reference source attached that ERRORS for "chrM" — the
	// decode-side analogue of `fasta: contig "chrM" not in index`. Before the
	// fix the decoder fetched the reference by RefSeqID regardless of RR and
	// surfaced this error; after the fix it honours RR=0 and never asks.
	rr, err := NewRecordReader(bytes.NewReader(buf.Bytes()))
	if err != nil {
		t.Fatalf("NewRecordReader: %v", err)
	}
	rr.SetReference(&stubReference{
		name: "other",
		seq:  "ACGTACGTACGT",
		err:  errFormat(`fasta: contig "chrM" not in index`),
	})
	out, err := rr.ReadAll()
	if err != nil {
		t.Fatalf("ReadAll of a reference-free CRAM whose contig is absent from the -T reference must succeed, got: %v", err)
	}
	if rr.NeedsReference() {
		t.Error("a reference-free slice must not report NeedsReference")
	}
	if len(out) != len(records) {
		t.Fatalf("decoded %d records, want %d", len(out), len(records))
	}
	// The bases were stored verbatim, so SEQ must round-trip exactly without
	// any reference involvement.
	for i := range records {
		if got, want := normSeq(out[i].Seq), records[i].Seq; got != want {
			t.Errorf("record %d SEQ = %q, want %q (verbatim, reference-free)", i, got, want)
		}
		if out[i].Pos != records[i].Pos {
			t.Errorf("record %d Pos = %d, want %d", i, out[i].Pos, records[i].Pos)
		}
		if out[i].Cigar.String() != records[i].Cigar.String() {
			t.Errorf("record %d CIGAR = %q, want %q", i, out[i].Cigar.String(), records[i].Cigar.String())
		}
	}
}

// TestRefFreeAbsentContigSetsRRZero verifies the encode side directly: a
// container whose records are all on a contig absent from the reference map
// must carry RR=0 (reference NOT required) in its compression header, so a
// decoder knows not to fetch the external reference. This is the on-disk
// signal that makes TestRefFreeAbsentContigRoundTrips possible and that
// upstream htslib honours via its no_ref flag.
func TestRefFreeAbsentContigSetsRRZero(t *testing.T) {
	h := absentContigHeader()
	rec := &sam.Record{
		QName: "read1", Flag: 0, RName: "chrM", Pos: 1, MapQ: 60,
		Cigar: mustCigarRF(t, "10M"), RNext: "*", PNext: 0, TLen: 0,
		Seq: "ACGTACGTAC", Qual: bytes.Repeat([]byte{'I'}, 10),
	}
	// Encode a full CRAM stream with a -T reference map that omits "chrM";
	// the writer must fall back to reference-free for that contig and write
	// RR=0 in the data container's compression header.
	var buf bytes.Buffer
	rw, err := NewRecordWriterOpts(&buf, h, WriterOptions{
		Reference:     map[string][]byte{"chr1": []byte("ACGTACGTACGTACGT")},
		ReferencePath: "/some/mismatched/ref.fa",
		EncodeThreads: 1,
	})
	if err != nil {
		t.Fatalf("NewRecordWriterOpts: %v", err)
	}
	if err := rw.Write(rec); err != nil {
		t.Fatalf("Write: %v", err)
	}
	if err := rw.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	rd, err := NewReader(bytes.NewReader(buf.Bytes()))
	if err != nil {
		t.Fatalf("NewReader: %v", err)
	}
	conts, err := rd.Containers()
	if err != nil {
		t.Fatalf("Containers: %v", err)
	}
	// Find the data container (the first one is the file-header container,
	// which carries no compression header).
	sawData := false
	for _, c := range conts {
		dc, derr := ParseDataContainer(c)
		if derr != nil {
			continue // not a data container (e.g. the SAM-header container).
		}
		sawData = true
		if dc.Compression.Preservation.ReferenceRequired {
			t.Error("a container with a contig absent from the reference map must write RR=0 (reference not required)")
		}
	}
	if !sawData {
		t.Fatal("no data container parsed from the CRAM stream")
	}
}
