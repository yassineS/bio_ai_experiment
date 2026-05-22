package bcftools

import (
	"bytes"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/gff"
)

// chr1 reference (1-based): positions 1..60
//   1234567890123456789012345678901234567890123456789012345678901234567890
//                  111111111122222222223333333333444444444455555555556666666666
//   "ATGGCGTACAACTAATGAATGAATGAATGAATGAATGAATGAATGAATGAATGAATGAATG"
// CDS exon: 1..30 (forward strand).
//   Codon 1 (pos 1-3) = ATG -> M
//   Codon 2 (pos 4-6) = GCG -> A
//   Codon 3 (pos 7-9) = TAC -> Y
//   Codon 4 (pos 10-12) = AAC -> N
//   Codon 5 (pos 13-15) = TAA -> * (stop)
//
// SNP fixtures:
//   chr1:7  T>A  -> codon 3 TAC->AAC -> Y->N missense
//   chr1:8  A>G  -> codon 3 TAC->TGC -> Y->C missense (random)
//   chr1:9  C>A  -> codon 3 TAC->TAA -> Y->* stop_gained
//   chr1:6  G>C  -> codon 2 GCG->GCC -> A->A synonymous (third codon position)

func buildCSQIndex(t *testing.T) *CSQIndex {
	t.Helper()
	const refSeq = "ATGGCGTACAACTAATGAATGAATGAATGAATGAATGAATGAATGAATGAATGAATGAATG"
	idx := &CSQIndex{
		Refs:        map[string][]byte{"chr1": []byte(refSeq)},
		Transcripts: map[string]*CSQTranscript{},
		ByChrom:     map[string][]*CSQTranscript{},
	}
	tx := &CSQTranscript{
		ID:      "tx1",
		Gene:    "GENE",
		Biotype: "protein_coding",
		Chrom:   "chr1",
		Strand:  gff.StrandForward,
		CDSExons: []CSQExon{
			{Start: 1, End: 30, Phase: 0},
		},
		Exons: []CSQExon{
			{Start: 1, End: 30},
		},
		Beg:    1,
		End:    30,
		Coding: true,
	}
	idx.Transcripts[tx.ID] = tx
	idx.ByChrom["chr1"] = []*CSQTranscript{tx}
	return idx
}

// TestCSQEndToEnd builds a tiny VCF and runs the full CSQ pipeline.
func TestCSQEndToEnd(t *testing.T) {
	idx := buildCSQIndex(t)
	// The haplotype engine requires position-sorted input, like
	// upstream bcftools csq.
	const vcfIn = `##fileformat=VCFv4.2
##contig=<ID=chr1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	6	.	G	C	.	PASS	DP=10
chr1	7	.	T	A	.	PASS	DP=10
chr1	100	.	A	C	.	PASS	DP=10
`
	var out bytes.Buffer
	n, err := CSQ(strings.NewReader(vcfIn), &out, idx, CSQOptions{CustomTag: "BCSQ"})
	if err != nil {
		t.Fatalf("CSQ: %v", err)
	}
	if n != 3 {
		t.Errorf("expected 3 records, got %d", n)
	}
	got := out.String()
	if !strings.Contains(got, "##INFO=<ID=BCSQ,") {
		t.Errorf("missing INFO header for BCSQ:\n%s", got)
	}
	if !strings.Contains(got, "BCSQ=missense|GENE|tx1") {
		t.Errorf("missing missense entry:\n%s", got)
	}
	if !strings.Contains(got, "BCSQ=synonymous|GENE|tx1") {
		t.Errorf("missing synonymous entry:\n%s", got)
	}
	// chr1:100 is outside CDS, should NOT have BCSQ tag.
	lines := strings.Split(got, "\n")
	for _, line := range lines {
		if strings.HasPrefix(line, "chr1\t100\t") && strings.Contains(line, "BCSQ=") {
			t.Errorf("chr1:100 unexpectedly annotated: %s", line)
		}
	}
}

// TestClassifyIndelsInCDS pins that the per-record classifier (this
// slice) now handles indels: a 1bp insertion and a 1bp deletion inside
// the CDS both shift the reading frame and are reported as frameshift.
func TestClassifyIndelsInCDS(t *testing.T) {
	idx := buildCSQIndex(t)
	// Each indel is classified in its own run so the no-GT haplotype
	// engine does not combine them onto a single haplotype: a 1bp
	// insertion and a 1bp deletion inside the CDS are both frameshifts.
	hdr := "##fileformat=VCFv4.2\n##contig=<ID=chr1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	for _, vcfIn := range []string{
		hdr + "chr1\t4\t.\tG\tGA\t.\tPASS\tDP=10\n",
		hdr + "chr1\t8\t.\tAC\tA\t.\tPASS\tDP=10\n",
	} {
		var out bytes.Buffer
		if _, err := CSQ(strings.NewReader(vcfIn), &out, idx, CSQOptions{}); err != nil {
			t.Fatalf("CSQ: %v", err)
		}
		if got := out.String(); !strings.Contains(got, "BCSQ=frameshift|GENE|tx1") {
			t.Errorf("expected indel classified as frameshift:\n%s", got)
		}
	}
}

func TestTranslateCodon(t *testing.T) {
	cases := []struct {
		codon string
		want  byte
	}{
		{"ATG", 'M'},
		{"TAA", '*'},
		{"TAG", '*'},
		{"TGA", '*'},
		{"GCT", 'A'},
		{"NNN", 'X'},
		{"AAA", 'K'},
	}
	for _, tc := range cases {
		var c [3]byte
		copy(c[:], []byte(tc.codon))
		got := translateCodon(c)
		if got != tc.want {
			t.Errorf("translateCodon(%q) = %c want %c", tc.codon, got, tc.want)
		}
	}
}

func TestParseCSQPhase(t *testing.T) {
	cases := []struct {
		in   string
		want byte
		err  bool
	}{
		{"", 'r', false},
		{"a", 'a', false},
		{"m", 'm', false},
		{"r", 'r', false},
		{"R", 'R', false},
		{"s", 's', false},
		{"x", 0, true},
	}
	for _, tc := range cases {
		got, err := ParseCSQPhase(tc.in)
		if tc.err && err == nil {
			t.Errorf("ParseCSQPhase(%q) expected err", tc.in)
		}
		if !tc.err && err != nil {
			t.Errorf("ParseCSQPhase(%q) unexpected err: %v", tc.in, err)
		}
		if !tc.err && got != tc.want {
			t.Errorf("ParseCSQPhase(%q) = %c want %c", tc.in, got, tc.want)
		}
	}
}

// TestLoadCSQIndexFromGFF exercises the GFF parser path. We write a
// tiny GFF + FASTA fixture and call loadCSQIndex.
func TestLoadCSQIndexFromGFF(t *testing.T) {
	dir := t.TempDir()
	fastaPath := dir + "/ref.fa"
	gffPath := dir + "/anno.gff3"
	if err := writeFile(fastaPath, ">chr1\nATGGCGTACAACTAA\n"); err != nil {
		t.Fatal(err)
	}
	if err := writeFile(gffPath, `##gff-version 3
chr1	src	gene	1	15	.	+	.	ID=g1;Name=GENE;biotype=protein_coding
chr1	src	mRNA	1	15	.	+	.	ID=tx1;Parent=g1
chr1	src	CDS	1	15	.	+	0	ID=cds1;Parent=tx1
`); err != nil {
		t.Fatal(err)
	}
	idx, err := loadCSQIndex(fastaPath, gffPath)
	if err != nil {
		t.Fatalf("loadCSQIndex: %v", err)
	}
	if len(idx.Refs["chr1"]) != 15 {
		t.Errorf("expected 15-bp refseq, got %d", len(idx.Refs["chr1"]))
	}
	tx, ok := idx.Transcripts["tx1"]
	if !ok {
		t.Fatalf("tx1 missing")
	}
	if tx.Gene != "GENE" || tx.Biotype != "protein_coding" {
		t.Errorf("unexpected transcript: %#v", tx)
	}
	if len(tx.CDSExons) != 1 || tx.CDSExons[0].Start != 1 || tx.CDSExons[0].End != 15 {
		t.Errorf("unexpected CDS: %#v", tx.CDSExons)
	}
}

// TestCDSOffsetHonoursPhase pins PR #110 review finding #1: cdsOffset
// must consume CDSExon.Phase, the GFF3 5'-leftover-codon-base count.
// With phase=1 on the first CDS exon, position 2 (which would normally
// be CDS offset 1) becomes offset 0 — the start of the reading frame.
func TestCDSOffsetHonoursPhase(t *testing.T) {
	tx := &CSQTranscript{
		Strand: gff.StrandForward,
		CDSExons: []CSQExon{
			{Start: 1, End: 30, Phase: 1},
		},
	}
	// Without phase: pos=2 -> offset 1. With phase=1: pos=2 -> offset 0.
	off, ok := cdsOffset(tx, 2)
	if !ok {
		t.Fatalf("cdsOffset(2): !ok")
	}
	if off != 0 {
		t.Errorf("phase=1: cdsOffset(2) = %d, want 0", off)
	}
	// And pos=1 (the leftover base) -> offset -1, which the caller
	// uses to detect "before frame start" — still returned as ok.
	off, ok = cdsOffset(tx, 1)
	if !ok {
		t.Fatalf("cdsOffset(1): !ok")
	}
	if off != -1 {
		t.Errorf("phase=1: cdsOffset(1) = %d, want -1", off)
	}
}

// TestCSQPhaseRequire pins parity with upstream csq.c (~line 3274):
// an unphased heterozygous genotype under -p require must abort.
func TestCSQPhaseRequire(t *testing.T) {
	idx := buildCSQIndex(t)
	// chr1:7 T>A is a missense inside the CDS; sample S1 carries an
	// unphased het 0/1.
	const vcfIn = `##fileformat=VCFv4.2
##contig=<ID=chr1>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	7	.	T	A	.	PASS	DP=10	GT	0/1
`
	var out bytes.Buffer
	_, err := CSQ(strings.NewReader(vcfIn), &out, idx, CSQOptions{Phase: 'r'})
	if err == nil {
		t.Fatalf("expected error for unphased het under -p require")
	}
	if !strings.Contains(err.Error(), "Unphased heterozygous genotype") {
		t.Errorf("error = %q, want it to mention 'Unphased heterozygous genotype'", err)
	}
	// A phased het (0|1) must NOT error.
	phased := strings.Replace(vcfIn, "GT\t0/1", "GT\t0|1", 1)
	var out2 bytes.Buffer
	if _, err := CSQ(strings.NewReader(phased), &out2, idx, CSQOptions{Phase: 'r'}); err != nil {
		t.Errorf("phased het under -p require unexpectedly errored: %v", err)
	}
}

func writeFile(path, content string) error {
	return writeFileBytes(path, []byte(content))
}
