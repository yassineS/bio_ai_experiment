package fixtures

import (
	"bufio"
	"compress/gzip"
	"math/rand"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// testContigs builds a small contig set with seeded sequences for the writers.
func testContigs(t *testing.T, p Params, seed int64) []contig {
	t.Helper()
	cs := makeContigs(p)
	rng := rand.New(rand.NewSource(seed))
	// writeFasta fills in c.Seq; do it through the FASTA writer for realism.
	dir := t.TempDir()
	if err := writeFasta(filepath.Join(dir, "ref.fa"), cs, rng); err != nil {
		t.Fatalf("writeFasta: %v", err)
	}
	return cs
}

// readLines reads a (possibly gzipped) file into lines, trimming the trailing
// newline.
func readLines(t *testing.T, path string) []string {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("open %s: %v", path, err)
	}
	defer f.Close()
	var r interface {
		ReadString(byte) (string, error)
	}
	if strings.HasSuffix(path, ".gz") {
		gr, err := gzip.NewReader(f)
		if err != nil {
			t.Fatalf("gzip %s: %v", path, err)
		}
		defer gr.Close()
		r = bufio.NewReader(gr)
	} else {
		r = bufio.NewReader(f)
	}
	var out []string
	for {
		line, err := r.ReadString('\n')
		if line != "" {
			out = append(out, strings.TrimRight(line, "\n"))
		}
		if err != nil {
			break
		}
	}
	return out
}

// TestWriteFastqSE checks the single-end FASTQ writer produces a well-formed,
// deterministic file with both plain and gzip variants carrying identical
// content.
func TestWriteFastqSE(t *testing.T) {
	p := Params{NumContigs: 2, ContigLen: 2000, FastqReads: 50, FastqReadLen: 80}
	cs := testContigs(t, p, 1)
	dir := t.TempDir()
	plain := filepath.Join(dir, "r.fastq")
	gz := plain + ".gz"
	rng := rand.New(rand.NewSource(7))
	if err := writeFastqSE(plain, gz, p, rng); err != nil {
		t.Fatalf("writeFastqSE: %v", err)
	}
	_ = cs
	lines := readLines(t, plain)
	if len(lines) != p.FastqReads*4 {
		t.Fatalf("expected %d lines, got %d", p.FastqReads*4, len(lines))
	}
	// Record structure: @name / seq / + / qual; seq and qual same length.
	for i := 0; i < len(lines); i += 4 {
		if !strings.HasPrefix(lines[i], "@read") {
			t.Fatalf("record %d header malformed: %q", i/4, lines[i])
		}
		if lines[i+2] != "+" {
			t.Fatalf("record %d separator = %q, want +", i/4, lines[i+2])
		}
		if len(lines[i+1]) != len(lines[i+3]) {
			t.Fatalf("record %d seq len %d != qual len %d", i/4, len(lines[i+1]), len(lines[i+3]))
		}
	}
	// Plain and gzip carry identical content.
	gzLines := readLines(t, gz)
	if strings.Join(lines, "\n") != strings.Join(gzLines, "\n") {
		t.Fatal("plain and gzip FASTQ content differ")
	}
}

// TestWriteFastqSE_Deterministic checks two runs with the same seed are
// byte-identical and an adapter appears on some reads.
func TestWriteFastqSE_Deterministic(t *testing.T) {
	p := Params{NumContigs: 1, ContigLen: 2000, FastqReads: 200, FastqReadLen: 100}
	gen := func() string {
		dir := t.TempDir()
		plain := filepath.Join(dir, "r.fastq")
		if err := writeFastqSE(plain, "", p, rand.New(rand.NewSource(42))); err != nil {
			t.Fatal(err)
		}
		return strings.Join(readLines(t, plain), "\n")
	}
	a, b := gen(), gen()
	if a != b {
		t.Fatal("same-seed FASTQ generation is not deterministic")
	}
	if !strings.Contains(a, illuminaAdapter[:20]) {
		t.Fatal("expected adapter contamination on some reads")
	}
	if !strings.Contains(a, "N") {
		t.Fatal("expected N bases in low-quality tails")
	}
}

// TestWriteFastqPE checks the paired writer produces matched-name R1/R2 records.
func TestWriteFastqPE(t *testing.T) {
	p := Params{NumContigs: 1, ContigLen: 2000, FastqReads: 30, FastqReadLen: 90}
	dir := t.TempDir()
	r1 := filepath.Join(dir, "r1.fastq")
	r2 := filepath.Join(dir, "r2.fastq")
	if err := writeFastqPE(r1, r2, p, rand.New(rand.NewSource(3))); err != nil {
		t.Fatalf("writeFastqPE: %v", err)
	}
	l1, l2 := readLines(t, r1), readLines(t, r2)
	if len(l1) != p.FastqReads*4 || len(l2) != p.FastqReads*4 {
		t.Fatalf("R1=%d R2=%d lines, want %d each", len(l1), len(l2), p.FastqReads*4)
	}
	for i := 0; i < len(l1); i += 4 {
		n1 := strings.TrimSuffix(l1[i], "/1")
		n2 := strings.TrimSuffix(l2[i], "/2")
		if n1 != n2 {
			t.Fatalf("pair %d name mismatch: %q vs %q", i/4, l1[i], l2[i])
		}
	}
}

// TestWriteGFF checks the GFF3 writer emits a valid header and well-ordered
// gene/mRNA/exon/CDS rows over the contig coordinate space.
func TestWriteGFF(t *testing.T) {
	p := Params{NumContigs: 2, ContigLen: 5000, Genes: 10}
	cs := testContigs(t, p, 1)
	dir := t.TempDir()
	path := filepath.Join(dir, "a.gff3")
	if err := writeGFF(path, cs, p, rand.New(rand.NewSource(9))); err != nil {
		t.Fatalf("writeGFF: %v", err)
	}
	lines := readLines(t, path)
	if lines[0] != "##gff-version 3" {
		t.Fatalf("first line = %q, want ##gff-version 3", lines[0])
	}
	var nGene, nMRNA, nExon, nCDS int
	for _, ln := range lines {
		if strings.HasPrefix(ln, "#") {
			continue
		}
		f := strings.Split(ln, "\t")
		if len(f) != 9 {
			t.Fatalf("GFF row has %d columns, want 9: %q", len(f), ln)
		}
		switch f[2] {
		case "gene":
			nGene++
		case "mRNA":
			nMRNA++
		case "exon":
			nExon++
		case "CDS":
			nCDS++
		}
	}
	if nGene != p.Genes || nMRNA != p.Genes {
		t.Fatalf("genes=%d mRNA=%d, want %d each", nGene, nMRNA, p.Genes)
	}
	if nExon == 0 || nCDS == 0 || nExon != nCDS {
		t.Fatalf("exon=%d CDS=%d, want equal and non-zero", nExon, nCDS)
	}
}

// TestWriteMultiSampleVCF checks the multi-sample VCF carries the requested
// number of sample columns and only valid biallelic genotypes.
func TestWriteMultiSampleVCF(t *testing.T) {
	p := Params{NumContigs: 2, ContigLen: 5000, Variants: 40, MultiSamples: 5}
	cs := testContigs(t, p, 1)
	dir := t.TempDir()
	path := filepath.Join(dir, "m.vcf")
	if err := writeMultiSampleVCF(path, cs, p, rand.New(rand.NewSource(11))); err != nil {
		t.Fatalf("writeMultiSampleVCF: %v", err)
	}
	lines := readLines(t, path)
	var sawHeader bool
	for _, ln := range lines {
		if strings.HasPrefix(ln, "#CHROM") {
			cols := strings.Split(ln, "\t")
			// 9 fixed columns + MultiSamples sample columns.
			if got := len(cols) - 9; got != p.MultiSamples {
				t.Fatalf("header has %d sample columns, want %d", got, p.MultiSamples)
			}
			sawHeader = true
			continue
		}
		if strings.HasPrefix(ln, "#") {
			continue
		}
		cols := strings.Split(ln, "\t")
		for _, s := range cols[9:] {
			gt := strings.SplitN(s, ":", 2)[0]
			switch gt {
			case "0/0", "0/1", "1/1":
			default:
				t.Fatalf("invalid biallelic genotype %q in %q", gt, ln)
			}
		}
	}
	if !sawHeader {
		t.Fatal("no #CHROM header line found")
	}
}
