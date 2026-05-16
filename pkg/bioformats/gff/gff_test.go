package gff

import (
	"io"
	"strings"
	"testing"
)

func TestParseLine(t *testing.T) {
	const line = "chr1\thavana\tCDS\t100\t199\t.\t+\t0\tID=cds1;Parent=t1;gene_name=GENE"
	f, err := parseLine(line)
	if err != nil {
		t.Fatalf("parseLine: %v", err)
	}
	if f.Seqid != "chr1" || f.Type != "CDS" || f.Start != 100 || f.End != 199 || f.Strand != StrandForward || f.Phase != 0 {
		t.Errorf("unexpected feature: %#v", f)
	}
	if f.ID() != "cds1" || f.Parent() != "t1" {
		t.Errorf("attrs wrong: id=%q parent=%q", f.ID(), f.Parent())
	}
	if got := f.Attributes["gene_name"]; got != "GENE" {
		t.Errorf("gene_name=%q want GENE", got)
	}
}

func TestReaderSkipsCommentsAndBlanks(t *testing.T) {
	const input = `##gff-version 3
# a comment

chr1	src	gene	1	100	.	+	.	ID=g1
chr1	src	mRNA	1	100	.	+	.	ID=t1;Parent=g1
chr1	src	CDS	10	30	.	+	0	ID=c1;Parent=t1
`
	r := NewReader(strings.NewReader(input))
	got, err := r.ReadAll()
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if len(got) != 3 {
		t.Fatalf("expected 3 records, got %d", len(got))
	}
	if got[0].Type != "gene" || got[1].Type != "mRNA" || got[2].Type != "CDS" {
		t.Errorf("types: %s, %s, %s", got[0].Type, got[1].Type, got[2].Type)
	}
}

func TestReaderEOF(t *testing.T) {
	r := NewReader(strings.NewReader(""))
	if _, err := r.Read(); err != io.EOF {
		t.Errorf("expected EOF, got %v", err)
	}
}

func TestParseLineErrors(t *testing.T) {
	cases := []struct {
		name, line string
	}{
		{"too-few-cols", "chr1\tsrc\tCDS\t1"},
		{"bad-start", "chr1\tsrc\tCDS\tone\t10\t.\t+\t.\t."},
		{"bad-end", "chr1\tsrc\tCDS\t1\ttwenty\t.\t+\t.\t."},
		{"bad-phase", "chr1\tsrc\tCDS\t1\t10\t.\t+\tQ\t."},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			if _, err := parseLine(tc.line); err == nil {
				t.Errorf("expected error for %q", tc.line)
			}
		})
	}
}

func TestParseAttributesEmpty(t *testing.T) {
	if m := parseAttributes("."); len(m) != 0 {
		t.Errorf("expected empty map, got %v", m)
	}
	if m := parseAttributes(""); len(m) != 0 {
		t.Errorf("expected empty map, got %v", m)
	}
}

func TestParseAttributesGTF(t *testing.T) {
	m := parseAttributes(`gene_id "ENSG00000001"; transcript_id "ENST00000001"`)
	if m["gene_id"] != "ENSG00000001" {
		t.Errorf("gtf gene_id=%q", m["gene_id"])
	}
	if m["transcript_id"] != "ENST00000001" {
		t.Errorf("gtf transcript_id=%q", m["transcript_id"])
	}
}
