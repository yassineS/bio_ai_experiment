package seqtk

import (
	"bufio"
	"bytes"
	"errors"
	"strings"
	"testing"
)

func TestIsGapByte(t *testing.T) {
	// A/C/G/T upper and lower must not be gap bytes. NOTE: U/u
	// (uracil) IS a gap, matching upstream seq_nt6_table['U']==5
	// (seqtk.c:208). The previous test rubber-stamped a bug that
	// mapped U/u to 4; reviewer-caught regression on PR #112.
	for _, b := range []byte("ACGTacgt") {
		if IsGapByte(b) {
			t.Errorf("IsGapByte(%q) = true, want false", b)
		}
	}
	// Every other byte (including U/u) must be a gap byte.
	for _, b := range []byte("UuNnRrSsBb-0?@") {
		if !IsGapByte(b) {
			t.Errorf("IsGapByte(%q) = false, want true", b)
		}
	}
}

func TestScanFASTA_PropagatesEmitError(t *testing.T) {
	const fasta = ">a\nACGT\n>b\nNNNN\n"
	want := errors.New("boom")
	err := scanFASTA(strings.NewReader(fasta), func(name string, seq []byte) error {
		if name == "b" {
			return want
		}
		return nil
	})
	if !errors.Is(err, want) {
		t.Fatalf("scanFASTA error = %v, want %v", err, want)
	}
}

func TestScanFASTA_VisitsEveryRecordInOrder(t *testing.T) {
	const fasta = ">a\nACGT\n>b\nGGCC\n>c\nTTAA\n"
	var seen []string
	err := scanFASTA(strings.NewReader(fasta), func(name string, seq []byte) error {
		seen = append(seen, name+"="+string(seq))
		return nil
	})
	if err != nil {
		t.Fatalf("scanFASTA: %v", err)
	}
	got := strings.Join(seen, ",")
	if got != "a=ACGT,b=GGCC,c=TTAA" {
		t.Fatalf("ordering wrong: %q", got)
	}
}

func TestWriteBED_Helpers(t *testing.T) {
	var buf bytes.Buffer
	bw := bufio.NewWriter(&buf)
	if err := writeBED3(bw, "chr1", 0, 5); err != nil {
		t.Fatalf("writeBED3: %v", err)
	}
	if err := writeBED4Int(bw, "chr2", 10, 20, 7); err != nil {
		t.Fatalf("writeBED4Int: %v", err)
	}
	if err := bw.Flush(); err != nil {
		t.Fatalf("flush: %v", err)
	}
	if got := buf.String(); got != "chr1\t0\t5\nchr2\t10\t20\t7\n" {
		t.Fatalf("BED output mismatch: %q", got)
	}
}
