package sam

import (
	"bufio"
	"bytes"
	"strings"
	"testing"
)

func TestParseHeader(t *testing.T) {
	in := "@HD\tVN:1.6\tSO:coordinate\n" +
		"@SQ\tSN:chr1\tLN:248956422\n" +
		"@SQ\tSN:chr2\tLN:242193529\tM5:abc\n" +
		"@RG\tID:rg1\tSM:sample1\tLB:lib1\n" +
		"@PG\tID:bwa\tPN:bwa\tVN:0.7.17\n" +
		"@CO\tHello, world\n" +
		"read1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGT\tIIII\n"
	br := bufio.NewReader(strings.NewReader(in))
	h, _, err := ParseHeader(br)
	if err != nil {
		t.Fatalf("ParseHeader: %v", err)
	}
	if got := len(h.Lines); got != 6 {
		t.Fatalf("expected 6 header lines, got %d", got)
	}
	if got := len(h.Refs); got != 2 {
		t.Fatalf("expected 2 refs, got %d", got)
	}
	if h.Refs[0].Name != "chr1" || h.Refs[0].Length != 248956422 {
		t.Errorf("bad first ref: %+v", h.Refs[0])
	}
	if h.Refs[1].Name != "chr2" || h.Refs[1].Length != 242193529 {
		t.Errorf("bad second ref: %+v", h.Refs[1])
	}
	if len(h.Refs[1].Extra) != 1 || h.Refs[1].Extra[0].Tag != "M5" {
		t.Errorf("expected extra M5 field on chr2, got %+v", h.Refs[1].Extra)
	}
	if len(h.ReadGroups) != 1 || h.ReadGroups[0].ID != "rg1" {
		t.Errorf("bad RG parse: %+v", h.ReadGroups)
	}
	if len(h.Programs) != 1 || h.Programs[0].ID != "bwa" {
		t.Errorf("bad PG parse: %+v", h.Programs)
	}
	if len(h.Comments) != 1 || h.Comments[0] != "Hello, world" {
		t.Errorf("bad CO parse: %+v", h.Comments)
	}
	// Buffer should now start with the record line.
	peek, _ := br.Peek(5)
	if string(peek) != "read1" {
		t.Errorf("expected to stop at body line, got peek=%q", peek)
	}
}

func TestHeaderRoundTrip(t *testing.T) {
	in := "@HD\tVN:1.6\tSO:coordinate\n" +
		"@SQ\tSN:chr1\tLN:1000\n" +
		"@RG\tID:rg1\tSM:s1\n" +
		"@CO\tcomment text\n"
	br := bufio.NewReader(strings.NewReader(in))
	h, _, err := ParseHeader(br)
	if err != nil {
		t.Fatalf("ParseHeader: %v", err)
	}
	var out bytes.Buffer
	if _, err := h.WriteTo(&out); err != nil {
		t.Fatalf("WriteTo: %v", err)
	}
	if got := out.String(); got != in {
		t.Errorf("header round-trip mismatch:\nin:  %q\nout: %q", in, got)
	}
	if got := h.Text(); got != in {
		t.Errorf("Text() mismatch:\nin:  %q\nout: %q", in, got)
	}
}

func TestParseHeaderErrors(t *testing.T) {
	tests := []struct {
		name string
		in   string
	}{
		{"missing tab", "@HDVN:1.6\n"},
		{"bad field", "@SQ\tSNchr1\n"},
		{"too short", "@H\n"},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			_, _, err := ParseHeader(bufio.NewReader(strings.NewReader(tc.in)))
			if err == nil {
				t.Fatalf("expected error for %q", tc.in)
			}
		})
	}
}

func TestRefIndex(t *testing.T) {
	h := &Header{Refs: []Reference{{Name: "chr1"}, {Name: "chr2"}}}
	if h.RefIndex("chr2") != 1 {
		t.Errorf("RefIndex chr2: got %d, want 1", h.RefIndex("chr2"))
	}
	if h.RefIndex("missing") != -1 {
		t.Errorf("RefIndex missing: got %d, want -1", h.RefIndex("missing"))
	}
}

func TestHeaderLineGet(t *testing.T) {
	h, err := ParseHeaderText("@HD\tVN:1.6\tSO:queryname\n")
	if err != nil {
		t.Fatalf("ParseHeaderText: %v", err)
	}
	hl := h.Lines[0]
	if v, ok := hl.Get("VN"); !ok || v != "1.6" {
		t.Errorf("Get(VN): %q, %v", v, ok)
	}
	if _, ok := hl.Get("XX"); ok {
		t.Error("Get(XX) should have returned false")
	}
}
