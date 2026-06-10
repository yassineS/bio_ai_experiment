package bcftools

import (
	"bytes"
	"compress/gzip"
	"io"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/gff"
)

// TestDumpCDSExonsTrim3Prime verifies the TRIM_3PRIME coordinate
// adjustment dumpCDSExons applies for incomplete coding lengths, on both
// strands, against hand-computed expectations.
func TestDumpCDSExonsTrim3Prime(t *testing.T) {
	tests := []struct {
		name   string
		strand gff.Strand
		in     []CSQExon
		want   []CSQExon
	}{
		{
			name:   "complete forward unchanged",
			strand: gff.StrandForward,
			in:     []CSQExon{{Start: 100, End: 102}},
			want:   []CSQExon{{Start: 100, End: 102}},
		},
		{
			name:   "forward trims 2 off high end",
			strand: gff.StrandForward,
			in:     []CSQExon{{Start: 100, End: 110}, {Start: 160, End: 171}}, // 11+12=23, 23%3=2
			want:   []CSQExon{{Start: 100, End: 110}, {Start: 160, End: 169}},
		},
		{
			name:   "reverse trims 2 off low end",
			strand: gff.StrandReverse,
			in:     []CSQExon{{Start: 106, End: 110}, {Start: 160, End: 171}}, // 5+12=17, 17%3=2
			want:   []CSQExon{{Start: 108, End: 110}, {Start: 160, End: 171}},
		},
		{
			name:   "reverse spills into adjacent exon",
			strand: gff.StrandReverse,
			in:     []CSQExon{{Start: 100, End: 100}, {Start: 200, End: 210}}, // 1+11=12... make %3!=0
			want:   []CSQExon{{Start: 100, End: 100}, {Start: 200, End: 210}},
		},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			tr := &CSQTranscript{Strand: tc.strand, Coding: true, CDSExons: tc.in}
			got := dumpCDSExons(tr)
			if len(got) != len(tc.want) {
				t.Fatalf("len=%d want %d", len(got), len(tc.want))
			}
			for i := range got {
				if got[i].Start != tc.want[i].Start || got[i].End != tc.want[i].End {
					t.Errorf("exon %d = %d-%d, want %d-%d", i, got[i].Start, got[i].End, tc.want[i].Start, tc.want[i].End)
				}
			}
			// The original slice must be untouched (dumpCDSExons copies).
			if &got[0] == &tr.CDSExons[0] {
				t.Errorf("dumpCDSExons must not alias the transcript's CDSExons")
			}
		})
	}
}

// TestDumpGFFSectionsAndFormat checks DumpGFF emits the five sections in
// upstream's order with the expected per-line format, on a small
// hand-built index (no upstream binary needed).
func TestDumpGFFSectionsAndFormat(t *testing.T) {
	idx := &CSQIndex{
		Transcripts: map[string]*CSQTranscript{},
		Genes: []*CSQGene{
			{ID: "G1", Name: "XYZ", Chrom: "1", Strand: gff.StrandForward, Beg: 90, End: 200, Used: true},
			{ID: "G0", Name: "nc", Chrom: "3", Strand: gff.StrandReverse, Beg: 20, End: 50, Used: false},
		},
	}
	tr := &CSQTranscript{
		ID: "T1", Gene: "XYZ", GeneID: "G1", Biotype: "protein_coding",
		Chrom: "1", Strand: gff.StrandForward, Beg: 90, End: 200, Coding: true, Used: true,
		CDSExons: []CSQExon{{Start: 100, End: 110, Phase: 0}, {Start: 160, End: 171, Phase: 0}},
		Exons:    []CSQExon{{Start: 90, End: 110}, {Start: 160, End: 200}},
		UTRs:     []CSQUTR{{Start: 90, End: 98, Prime5: true}, {Start: 172, End: 200, Prime5: false}},
	}
	nc := &CSQTranscript{
		ID: "T0", Gene: "nc", GeneID: "G0", Biotype: "lincRNA",
		Chrom: "3", Strand: gff.StrandReverse, Beg: 20, End: 50, Coding: false, Used: false,
	}
	idx.Transcripts["T1"] = tr
	idx.Transcripts["T0"] = nc

	var buf bytes.Buffer
	if err := DumpGFF(&buf, idx); err != nil {
		t.Fatalf("DumpGFF: %v", err)
	}
	gr, err := gzip.NewReader(&buf)
	if err != nil {
		t.Fatalf("gzip: %v", err)
	}
	gr.Multistream(true)
	raw, _ := io.ReadAll(gr)
	out := string(raw)

	wantLines := []string{
		"1\t.\tgene\t90\t200\t.\t+\t.\tID=G1;Name=XYZ;used=1",
		"3\t.\tgene\t20\t50\t.\t-\t.\tID=G0;Name=nc;used=0",
		"1\t.\tmRNA\t90\t200\t.\t+\t.\tID=T1;Parent=G1;biotype=protein_coding;used=1",
		"3\t.\tlincRNA\t20\t50\t.\t-\t.\tID=T0;Parent=G0;biotype=lincRNA;used=0",
		"1\t.\tCDS\t100\t110\t.\t+\t0\tParent=T1",
		"1\t.\tCDS\t160\t169\t.\t+\t0\tParent=T1", // 11+12=23, 23%3=2 -> trim 2 off 3'
		"1\t.\tfive_prime_UTR\t90\t98\t.\t+\t.\tParent=T1",
		"1\t.\tthree_prime_UTR\t172\t200\t.\t+\t.\tParent=T1",
		"1\t.\texon\t90\t110\t.\t+\t.\tParent=T1",
		"1\t.\texon\t160\t200\t.\t+\t.\tParent=T1",
	}
	for _, w := range wantLines {
		if !strings.Contains(out, w+"\n") {
			t.Errorf("dump missing line:\n%s\n--- full dump ---\n%s", w, out)
		}
	}
	// Section ordering: genes < transcripts < CDS < UTR < exon.
	order := []string{"\tgene\t", "\tmRNA\t", "\tCDS\t", "_prime_UTR\t", "\texon\t"}
	last := -1
	for _, marker := range order {
		i := strings.Index(out, marker)
		if i < 0 {
			t.Fatalf("marker %q not found", marker)
		}
		if i < last {
			t.Errorf("section %q out of order (at %d, previous %d)", marker, i, last)
		}
		last = i
	}
}

// TestPhaseChar checks the CDS phase column rendering.
func TestPhaseChar(t *testing.T) {
	cases := map[int]string{0: "0", 1: "1", 2: "2", 3: ".", -1: ".", 99: "."}
	for in, want := range cases {
		if got := phaseChar(in); got != want {
			t.Errorf("phaseChar(%d)=%q want %q", in, got, want)
		}
	}
}
