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
	}
	idx.Transcripts[tx.ID] = tx
	idx.ByChrom["chr1"] = []*CSQTranscript{tx}
	return idx
}

func TestClassifyMissenseSynonymousStop(t *testing.T) {
	idx := buildCSQIndex(t)
	cases := []struct {
		name        string
		pos         int
		refBase     byte
		altBase     byte
		consequence string
		aaChange    string
		dnaChange   string
	}{
		{"missense T>A at 7", 7, 'T', 'A', "missense", "3Y>N", "7T>A"},
		{"stop_gained C>A at 9", 9, 'C', 'A', "stop_gained", "3Y>*", "9C>A"},
		{"synonymous G>C at 6", 6, 'G', 'C', "synonymous", "2A>A", "6G>C"},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			entry, ok := classifyForTranscript(idx.Transcripts["tx1"], idx.Refs["chr1"], tc.pos, tc.refBase, tc.altBase)
			if !ok {
				t.Fatalf("classifyForTranscript returned !ok")
			}
			wantPrefix := tc.consequence + "|GENE|tx1|protein_coding|+|" + tc.aaChange + "|" + tc.dnaChange
			if entry != wantPrefix {
				t.Errorf("entry = %q want %q", entry, wantPrefix)
			}
		})
	}
}

func TestClassifyStartLost(t *testing.T) {
	idx := buildCSQIndex(t)
	// pos 1 is the A in ATG (Met). Mutate it to C, codon becomes CTG -> L.
	entry, ok := classifyForTranscript(idx.Transcripts["tx1"], idx.Refs["chr1"], 1, 'A', 'C')
	if !ok {
		t.Fatalf("classifyForTranscript returned !ok")
	}
	if !strings.HasPrefix(entry, "start_lost|") {
		t.Errorf("expected start_lost, got %q", entry)
	}
}

func TestClassifyReverseStrand(t *testing.T) {
	// Reverse-strand transcript at chr1:31..60.
	// Refseq positions 31..60 (1-based) is "ATGAATGAATGAATGAATGAATGAATGAAT".
	// Reverse complement (the actual coding sequence, 5'->3') is
	// "ATTCATTCATTCATTCATTCATTCATTCAT".
	// First codon ATT -> I (isoleucine).
	// Position 60 is the 1st base of the FIRST codon when reading 3'->5'
	// (genomic) = position 60's complement is the FIRST base of CDS coordinate 0.
	// Genomic base at pos 60 is 'T' -> complement 'A' (the 'A' in ATT).
	refSeq := []byte("ATGGCGTACAACTAATGAATGAATGAATGAATGAATGAATGAATGAATGAATGAATGAATG")
	tx := &CSQTranscript{
		ID:      "tx2",
		Gene:    "GENE2",
		Biotype: "protein_coding",
		Chrom:   "chr1",
		Strand:  gff.StrandReverse,
		CDSExons: []CSQExon{
			{Start: 31, End: 60, Phase: 0},
		},
	}
	// Mutate pos 60 T -> A. Complement(T)=A, complement(A)=T. So
	// first codon ATT becomes TTT -> F. I>F missense.
	entry, ok := classifyForTranscript(tx, refSeq, 60, 'T', 'A')
	if !ok {
		t.Fatalf("classifyForTranscript returned !ok")
	}
	if !strings.HasPrefix(entry, "missense|") {
		t.Errorf("expected missense, got %q", entry)
	}
	if !strings.Contains(entry, "|-|") {
		t.Errorf("expected '-' strand in %q", entry)
	}
}

func TestClassifyOutsideCDS(t *testing.T) {
	idx := buildCSQIndex(t)
	tx := idx.Transcripts["tx1"]
	if cdsCovers(tx, 50) {
		t.Errorf("cdsCovers should be false for pos 50 (CDS is 1..30)")
	}
}

func TestClassifyRefMismatch(t *testing.T) {
	idx := buildCSQIndex(t)
	// Pos 8 is T; passing refBase='G' must fail.
	if _, ok := classifyForTranscript(idx.Transcripts["tx1"], idx.Refs["chr1"], 8, 'G', 'A'); ok {
		t.Errorf("expected !ok when REF does not match FASTA")
	}
}

// TestCSQEndToEnd builds a tiny VCF and runs the full CSQ pipeline.
func TestCSQEndToEnd(t *testing.T) {
	idx := buildCSQIndex(t)
	const vcfIn = `##fileformat=VCFv4.2
##contig=<ID=chr1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	7	.	T	A	.	PASS	DP=10
chr1	6	.	G	C	.	PASS	DP=10
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

func TestClassifyNonSNPSkipped(t *testing.T) {
	idx := buildCSQIndex(t)
	const vcfIn = `##fileformat=VCFv4.2
##contig=<ID=chr1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	7	.	T	TA	.	PASS	DP=10
chr1	7	.	TA	T	.	PASS	DP=10
`
	var out bytes.Buffer
	if _, err := CSQ(strings.NewReader(vcfIn), &out, idx, CSQOptions{}); err != nil {
		t.Fatalf("CSQ: %v", err)
	}
	if strings.Contains(out.String(), "BCSQ=") {
		t.Errorf("indel should not be classified:\n%s", out.String())
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

// TestClassifyStartLost_MidCodon pins PR #110 review finding: start-loss
// must fire for ANY of the 3 ATG positions (codonIdx == 0), not just
// the leading A. Previously the check was `codingOff == 0` which only
// fired on the first base.
func TestClassifyStartLost_MidCodon(t *testing.T) {
	idx := buildCSQIndex(t)
	// pos 2 is the T in ATG. Mutate it to C -> ACG -> T (not Met) -> start_lost.
	entry, ok := classifyForTranscript(idx.Transcripts["tx1"], idx.Refs["chr1"], 2, 'T', 'C')
	if !ok {
		t.Fatalf("classifyForTranscript returned !ok")
	}
	if !strings.HasPrefix(entry, "start_lost|") {
		t.Errorf("pos 2 (T in ATG) expected start_lost, got %q", entry)
	}
	// pos 3 is the G in ATG. Mutate to T -> ATT -> Ile -> start_lost.
	entry, ok = classifyForTranscript(idx.Transcripts["tx1"], idx.Refs["chr1"], 3, 'G', 'T')
	if !ok {
		t.Fatalf("classifyForTranscript returned !ok")
	}
	if !strings.HasPrefix(entry, "start_lost|") {
		t.Errorf("pos 3 (G in ATG) expected start_lost, got %q", entry)
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

func writeFile(path, content string) error {
	return writeFileBytes(path, []byte(content))
}
