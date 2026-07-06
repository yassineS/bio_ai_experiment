package cram

import (
	"bytes"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// TestPartialReferenceDisablesM5 reproduces the real-data divergence where a
// BAM carries many @SQ contigs but the -T reference resolves only a subset
// (e.g. a ~2580-contig GIAB BAM against a chr20-only FASTA). Upstream htslib
// (cram_io.c cram_write_SAM_hdr) computes an M5 for every M5-less @SQ, but if
// ANY M5-less @SQ contig is unresolvable in the reference it falls back to
// embed_ref=2 and emits NO computed M5 tags at all — the disable is
// all-or-nothing across the whole header. Our writer must do the same: with
// chr1 resolvable and chr2 absent, NEITHER @SQ may carry a computed M5.
func TestPartialReferenceDisablesM5(t *testing.T) {
	h, err := sam.ParseHeaderText(
		"@HD\tVN:1.6\tSO:coordinate\n" +
			"@SQ\tSN:chr1\tLN:16\n" +
			"@SQ\tSN:chr2\tLN:16\n")
	if err != nil {
		t.Fatalf("ParseHeaderText: %v", err)
	}
	// Reference map resolves ONLY chr1; chr2 is absent (unresolvable).
	var buf bytes.Buffer
	rw, err := NewRecordWriterOpts(&buf, h, WriterOptions{
		Reference:     map[string][]byte{"chr1": []byte("ACGTACGTACGTACGT")},
		ReferencePath: "/some/partial/ref.fa",
		EncodeThreads: 1,
	})
	if err != nil {
		t.Fatalf("NewRecordWriterOpts: %v", err)
	}
	if err := rw.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	rr, err := NewRecordReader(bytes.NewReader(buf.Bytes()))
	if err != nil {
		t.Fatalf("NewRecordReader: %v", err)
	}
	text := rr.Header().Text()
	if strings.Contains(text, "M5:") {
		t.Errorf("partial reference must disable M5 for ALL @SQ, but the embedded header carries an M5 tag:\n%s", text)
	}
	// chr1 must still appear as a bare @SQ (no M5), confirming the resolvable
	// contig did NOT get an M5 it would otherwise have received.
	if !strings.Contains(text, "SN:chr1") || !strings.Contains(text, "SN:chr2") {
		t.Fatalf("both @SQ lines must survive; header:\n%s", text)
	}
}

// TestFullReferenceKeepsM5 is the control: when EVERY M5-less @SQ contig is
// resolvable, the all-or-nothing disable does NOT trigger and each @SQ gains
// its computed M5, exactly as before the partial-reference fix.
func TestFullReferenceKeepsM5(t *testing.T) {
	h, err := sam.ParseHeaderText(
		"@HD\tVN:1.6\tSO:coordinate\n" +
			"@SQ\tSN:chr1\tLN:16\n" +
			"@SQ\tSN:chr2\tLN:16\n")
	if err != nil {
		t.Fatalf("ParseHeaderText: %v", err)
	}
	var buf bytes.Buffer
	rw, err := NewRecordWriterOpts(&buf, h, WriterOptions{
		Reference: map[string][]byte{
			"chr1": []byte("ACGTACGTACGTACGT"),
			"chr2": []byte("TTTTGGGGCCCCAAAA"),
		},
		ReferencePath: "/some/full/ref.fa",
		EncodeThreads: 1,
	})
	if err != nil {
		t.Fatalf("NewRecordWriterOpts: %v", err)
	}
	if err := rw.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	rr, err := NewRecordReader(bytes.NewReader(buf.Bytes()))
	if err != nil {
		t.Fatalf("NewRecordReader: %v", err)
	}
	text := rr.Header().Text()
	if n := strings.Count(text, "M5:"); n != 2 {
		t.Errorf("a fully-resolvable reference must give both @SQ an M5, got %d M5 tags:\n%s", n, text)
	}
}
