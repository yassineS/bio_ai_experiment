package bedgetfasta

// Parity tests against the upstream bedtools getfasta test suite.
//
// Cases mirror reference_code/bedtools/test/getfasta/test-getfasta.sh.
// Inputs live under tools/bedgetfasta/testdata/parity/. We assert
// byte-for-byte equality against the inline expected outputs from the
// upstream test script.

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

func parityPath(name string) string {
	return filepath.Join("..", "..", "testdata", "parity", name)
}

func readParity(t *testing.T, name string) []byte {
	t.Helper()
	data, err := os.ReadFile(parityPath(name))
	if err != nil {
		t.Fatalf("read %s: %v", name, err)
	}
	return data
}

// stagedFasta copies a parity FASTA into a temp dir so we can build the .fai
// next to it without polluting the repo, and so each test gets a clean slate.
func stagedFasta(t *testing.T, name string) string {
	t.Helper()
	dir := t.TempDir()
	body := readParity(t, name)
	target := filepath.Join(dir, name)
	if err := os.WriteFile(target, body, 0o644); err != nil {
		t.Fatalf("stage %s: %v", name, err)
	}
	return target
}

// getfasta.t01 — `chr1\t1\t10` should yield 2 lines (header + seq).
func TestParity_Getfasta_T01_TwoLines(t *testing.T) {
	fa := stagedFasta(t, "t.fa")
	bed := []byte("chr1\t1\t10\n")
	var buf, warn bytes.Buffer
	if _, err := Run(bytes.NewReader(bed), fa, &buf, &warn, Options{}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	got := strings.TrimRight(buf.String(), "\n")
	lines := strings.Count(got, "\n") + 1
	if lines != 2 {
		t.Errorf("expected 2 lines, got %d: %q", lines, buf.String())
	}
	if got != ">chr1:1-10\nggggggggg" {
		t.Errorf("output = %q", got)
	}
}

// getfasta.t02 — `-split` over blocks.bed line 1 yields a 22-base sequence.
func TestParity_Getfasta_T02_SplitLen(t *testing.T) {
	fa := stagedFasta(t, "t.fa")
	bed := readParity(t, "blocks.bed")
	// Use only the first BED12 record (the upstream test does the same:
	// awk(NR==2) → look at line 2 of bedtools output, which is the seq of
	// the first record's split output).
	firstLine := strings.SplitN(string(bed), "\n", 2)[0] + "\n"
	var buf, warn bytes.Buffer
	if _, err := Run(strings.NewReader(firstLine), fa, &buf, &warn, Options{Split: true}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	lines := strings.Split(strings.TrimRight(buf.String(), "\n"), "\n")
	if len(lines) != 2 {
		t.Fatalf("expected 2 lines, got %d: %v", len(lines), lines)
	}
	if len(lines[1]) != 22 {
		t.Errorf("split seq length = %d, want 22", len(lines[1]))
	}
}

// getfasta.t03 — `-split` over both blocks.bed records: line 4 (the seq of
// the second BED12 record) should equal "cta".
func TestParity_Getfasta_T03_SplitSecondRecord(t *testing.T) {
	fa := stagedFasta(t, "t.fa")
	bed := readParity(t, "blocks.bed")
	var buf, warn bytes.Buffer
	if _, err := Run(bytes.NewReader(bed), fa, &buf, &warn, Options{Split: true}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	lines := strings.Split(strings.TrimRight(buf.String(), "\n"), "\n")
	if len(lines) < 4 {
		t.Fatalf("expected >=4 lines, got %d: %v", len(lines), lines)
	}
	if lines[3] != "cta" {
		t.Errorf("line 4 = %q, want \"cta\"", lines[3])
	}
}

// getfasta.t05 — `-split -s` on blocks.bed second record yields "tag".
func TestParity_Getfasta_T05_SplitStrandedRevComp(t *testing.T) {
	fa := stagedFasta(t, "t.fa")
	bed := readParity(t, "blocks.bed")
	var buf, warn bytes.Buffer
	if _, err := Run(bytes.NewReader(bed), fa, &buf, &warn, Options{Split: true, Strand: true}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	lines := strings.Split(strings.TrimRight(buf.String(), "\n"), "\n")
	if len(lines) < 4 {
		t.Fatalf("expected >=4 lines, got %d: %v", len(lines), lines)
	}
	if lines[3] != "tag" {
		t.Errorf("line 4 = %q, want \"tag\"", lines[3])
	}
}

// getfasta.t08 — IUPAC test: default emission preserves case.
func TestParity_Getfasta_T08_IUPACPreservesCase(t *testing.T) {
	fa := stagedFasta(t, "test.iupac.fa")
	bed := readParity(t, "test.iupac.bed")
	var buf, warn bytes.Buffer
	if _, err := Run(bytes.NewReader(bed), fa, &buf, &warn, Options{}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	want := []byte(">1:0-16\nAGCTYRWSKMDVHBXN\n>2:0-16\nagctyrwskmdvhbxn\n")
	if !bytes.Equal(buf.Bytes(), want) {
		t.Errorf("output mismatch.\nwant:\n%s\ngot:\n%s", want, buf.Bytes())
	}
}

// getfasta.t09 — IUPAC test with `-s`: revcomp preserves case.
func TestParity_Getfasta_T09_IUPACRevcompPreservesCase(t *testing.T) {
	fa := stagedFasta(t, "test.iupac.fa")
	bed := readParity(t, "test.iupac.bed")
	var buf, warn bytes.Buffer
	if _, err := Run(bytes.NewReader(bed), fa, &buf, &warn, Options{Strand: true}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	want := []byte(">1:0-16(-)\nNXVDBHKMSWYRAGCT\n>2:0-16(-)\nnxvdbhkmswyragct\n")
	if !bytes.Equal(buf.Bytes(), want) {
		t.Errorf("output mismatch.\nwant:\n%s\ngot:\n%s", want, buf.Bytes())
	}
}

// getfasta.t11 — default header is `chrom:start-end`, both records of
// blocks.bed emit the full 40bp slice.
func TestParity_Getfasta_T11_DefaultHeader(t *testing.T) {
	fa := stagedFasta(t, "t.fa")
	bed := readParity(t, "blocks.bed")
	var buf, warn bytes.Buffer
	if _, err := Run(bytes.NewReader(bed), fa, &buf, &warn, Options{}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	want := []byte(">chr1:0-40\nagggggggggcgggggggggtgggggggggaggggggggg\n>chr1:0-40\nagggggggggcgggggggggtgggggggggaggggggggg\n")
	if !bytes.Equal(buf.Bytes(), want) {
		t.Errorf("output mismatch.\nwant:\n%s\ngot:\n%s", want, buf.Bytes())
	}
}

// getfasta.t12a — `-name` prepends `<name>::` to the header.
func TestParity_Getfasta_T12_NameHeader(t *testing.T) {
	fa := stagedFasta(t, "t.fa")
	bed := readParity(t, "blocks.bed")
	var buf, warn bytes.Buffer
	if _, err := Run(bytes.NewReader(bed), fa, &buf, &warn, Options{Name: true}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	want := []byte(">three_blocks_match::chr1:0-40\nagggggggggcgggggggggtgggggggggaggggggggg\n>three_blocks_match::chr1:0-40\nagggggggggcgggggggggtgggggggggaggggggggg\n")
	if !bytes.Equal(buf.Bytes(), want) {
		t.Errorf("output mismatch.\nwant:\n%s\ngot:\n%s", want, buf.Bytes())
	}
}

// getfasta.t12b — `-nameOnly` emits just the name as header.
func TestParity_Getfasta_T12b_NameOnlyHeader(t *testing.T) {
	fa := stagedFasta(t, "t.fa")
	bed := readParity(t, "blocks.bed")
	var buf, warn bytes.Buffer
	if _, err := Run(bytes.NewReader(bed), fa, &buf, &warn, Options{NameOnly: true}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	want := []byte(">three_blocks_match\nagggggggggcgggggggggtgggggggggaggggggggg\n>three_blocks_match\nagggggggggcgggggggggtgggggggggaggggggggg\n")
	if !bytes.Equal(buf.Bytes(), want) {
		t.Errorf("output mismatch.\nwant:\n%s\ngot:\n%s", want, buf.Bytes())
	}
}

// getfasta.t13a — `-name -s` includes both `<name>::` and `(±)` suffix.
func TestParity_Getfasta_T13_NameStrand(t *testing.T) {
	fa := stagedFasta(t, "t.fa")
	bed := readParity(t, "blocks.bed")
	var buf, warn bytes.Buffer
	if _, err := Run(bytes.NewReader(bed), fa, &buf, &warn, Options{Name: true, Strand: true}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	want := []byte(">three_blocks_match::chr1:0-40(+)\nagggggggggcgggggggggtgggggggggaggggggggg\n>three_blocks_match::chr1:0-40(-)\nccccccccctcccccccccacccccccccgccccccccct\n")
	if !bytes.Equal(buf.Bytes(), want) {
		t.Errorf("output mismatch.\nwant:\n%s\ngot:\n%s", want, buf.Bytes())
	}
}

// getfasta.t16 — `-s -name` on rna.fasta (DNA stored as letters) with two
// records, one minus-strand, one plus-strand.
func TestParity_Getfasta_T16_RNAFixtureStrand(t *testing.T) {
	fa := stagedFasta(t, "rna.fasta")
	bed := readParity(t, "rna.bed")
	var buf, warn bytes.Buffer
	if _, err := Run(bytes.NewReader(bed), fa, &buf, &warn, Options{Strand: true, Name: true}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	want := []byte(">candidate_1::chr1:0-10(-)\ncatcggtcaa\n>candidate_2::chr1:0-10(+)\nuugaccgaug\n")
	if !bytes.Equal(buf.Bytes(), want) {
		t.Errorf("output mismatch.\nwant:\n%s\ngot:\n%s", want, buf.Bytes())
	}
}

// getfasta.t17 — `-rna` replaces T/t with U/u after any revcomp.
func TestParity_Getfasta_T17_RNAOption(t *testing.T) {
	fa := stagedFasta(t, "rna.fasta")
	bed := readParity(t, "rna.bed")
	var buf, warn bytes.Buffer
	if _, err := Run(bytes.NewReader(bed), fa, &buf, &warn, Options{Strand: true, Name: true, RNA: true}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	want := []byte(">candidate_1::chr1:0-10(-)\ncaucggucaa\n>candidate_2::chr1:0-10(+)\nuugaccgaug\n")
	if !bytes.Equal(buf.Bytes(), want) {
		t.Errorf("output mismatch.\nwant:\n%s\ngot:\n%s", want, buf.Bytes())
	}
}

// Upstream getfasta.t07 — warning when chromosome name from BED isn't in the
// FASTA. Our library emits the same WARNING line (with our wording) to the
// supplied warn writer; the wording differs slightly because upstream's
// `t_fH.fa` test relies on the whitespace-aware contig name parsing that
// our `-fullHeader` flag does not yet implement.
func TestParity_Getfasta_T07_FullHeader(t *testing.T) {
	t.Skip("`-fullHeader` not implemented: bedgetfasta uses the same first-token contig naming as samtools faidx, so a BED query for 'chr1 assembled by consortium X' against a FASTA with a multi-word header is treated as a missing chromosome — see PARITY_ROADMAP.md#bedtools")
}

// Upstream getfasta.t18 — input FASTA is bgzipped. Random-access FASTA
// support over BGZF needs a .gzi index; not yet implemented in
// pkg/bioformats/fasta.
func TestParity_Getfasta_T18_BGZF(t *testing.T) {
	t.Skip("BGZF (.fa.gz) FASTA input is not yet supported by pkg/bioformats/fasta random access; needs .gzi index — see PARITY_ROADMAP.md#bedtools")
}
