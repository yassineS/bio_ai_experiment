package runner

import (
	"compress/gzip"
	"os"
	"path/filepath"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pipeline/matrix"
)

// TestStripProvenanceSAM checks @PG/@CO removal without touching data lines.
func TestStripProvenanceSAM(t *testing.T) {
	in := "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n@PG\tID:samtools\tPN:samtools\tVN:1.22\n" +
		"read1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGT\tIIII\n"
	want := "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:100\n" +
		"read1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGT\tIIII\n"
	if got := string(stripProvenance([]byte(in))); got != want {
		t.Errorf("stripProvenance SAM:\n got=%q\nwant=%q", got, want)
	}
}

// TestStripProvenanceVCF checks command/version header removal.
func TestStripProvenanceVCF(t *testing.T) {
	in := "##fileformat=VCFv4.2\n##bcftools_viewCommand=view a.vcf; Date=...\n##contig=<ID=chr1>\nchr1\t1\t.\tA\tG\t60\tPASS\t.\n"
	got := string(stripProvenance([]byte(in)))
	if want := "##fileformat=VCFv4.2\n##contig=<ID=chr1>\nchr1\t1\t.\tA\tG\t60\tPASS\t.\n"; got != want {
		t.Errorf("stripProvenance VCF:\n got=%q\nwant=%q", got, want)
	}
}

// TestStripProvenanceStatsBlock checks the samtools/bcftools stats-style
// provenance comment block ("# This file was produced by ...", the command-line
// echo, the working-directory block, the bare "#" separator, and the gtcheck
// timing line) is removed while the data-describing comment rows survive.
func TestStripProvenanceStatsBlock(t *testing.T) {
	in := "# This file was produced by bcftools stats (1.23+htslib-1.23)\n" +
		"# The command line was:\tbcftools stats a.vcf\n" +
		"#\n" +
		"# ID\t[2]id\t[3]file names\n" +
		"SN\t0\tnumber of records:\t400\n" +
		"INFO\tTime required to process one record .. 0.000003 seconds\n"
	want := "# ID\t[2]id\t[3]file names\n" +
		"SN\t0\tnumber of records:\t400\n"
	if got := string(stripProvenance([]byte(in))); got != want {
		t.Errorf("stripProvenance stats block:\n got=%q\nwant=%q", got, want)
	}
}

// TestStripProvenanceFilterPass checks the auto-inserted ##FILTER=PASS
// boilerplate is dropped (its position differs between ours/upstream) while a
// real ##FILTER definition is preserved.
func TestStripProvenanceFilterPass(t *testing.T) {
	in := "##fileformat=VCFv4.2\n" +
		"##FILTER=<ID=PASS,Description=\"All filters passed\">\n" +
		"##FILTER=<ID=q10,Description=\"Quality below 10\">\n" +
		"chr1\t1\t.\tA\tG\t60\tPASS\t.\n"
	want := "##fileformat=VCFv4.2\n" +
		"##FILTER=<ID=q10,Description=\"Quality below 10\">\n" +
		"chr1\t1\t.\tA\tG\t60\tPASS\t.\n"
	if got := string(stripProvenance([]byte(in))); got != want {
		t.Errorf("stripProvenance FILTER=PASS:\n got=%q\nwant=%q", got, want)
	}
}

// TestCompareByteExact covers match and mismatch.
func TestCompareByteExact(t *testing.T) {
	a := []byte("@PG\tID:x\nchr1\t1\n")
	b := []byte("@PG\tID:y\nchr1\t1\n") // differs only in stripped @PG
	if r := CompareByteExact(a, b); !r.Equal {
		t.Errorf("expected equal after provenance strip, got %+v", r)
	}
	c := []byte("chr1\t2\n")
	if r := CompareByteExact(a, c); r.Equal {
		t.Errorf("expected mismatch, got equal")
	}
}

// TestCompareSimilarity covers numeric tolerance and structural mismatch.
func TestCompareSimilarity(t *testing.T) {
	a := []byte("chr1\t0.1000000\t2\n")
	b := []byte("chr1\t0.1000001\t2\n")
	r := CompareSimilarity(a, b, similarityEpsilon)
	if !r.Equal {
		t.Errorf("expected within-tolerance equal, got %+v", r)
	}
	if r.MaxDeviation == 0 {
		t.Errorf("expected non-zero recorded deviation")
	}

	d := CompareSimilarity([]byte("chr1\t1.0\n"), []byte("chr2\t1.0\n"), similarityEpsilon)
	if d.Equal {
		t.Errorf("expected non-numeric field mismatch to diverge")
	}

	e := CompareSimilarity([]byte("chr1\t1.0\n"), []byte("chr1\t2.0\n"), similarityEpsilon)
	if e.Equal {
		t.Errorf("expected out-of-tolerance numeric to diverge")
	}

	// A per-entry tolerance widens acceptance: a deviation that fails at the
	// default epsilon passes when the entry opts into a looser bound (the
	// bcftools call QUAL libm-last-ULP case).
	f := CompareSimilarity([]byte("chr1\t15.6999\n"), []byte("chr1\t15.6998\n"), resolveEpsilon(2e-5))
	if !f.Equal {
		t.Errorf("expected within widened tolerance, got %+v", f)
	}
	g := CompareSimilarity([]byte("chr1\t15.6999\n"), []byte("chr1\t15.6998\n"), resolveEpsilon(0))
	if g.Equal {
		t.Errorf("expected default tolerance to reject the QUAL last-ULP deviation")
	}
}

// TestCompareOutputFiles_ByteExact covers matching, mismatch, and gzip handling
// of the output-file comparison path used by vcftools and mosdepth.
func TestCompareOutputFiles_ByteExact(t *testing.T) {
	dir := t.TempDir()
	write := func(name, content string) {
		if err := os.WriteFile(filepath.Join(dir, name), []byte(content), 0o644); err != nil {
			t.Fatal(err)
		}
	}
	writeGz := func(name, content string) {
		f, err := os.Create(filepath.Join(dir, name))
		if err != nil {
			t.Fatal(err)
		}
		gw := gzip.NewWriter(f)
		if _, err := gw.Write([]byte(content)); err != nil {
			t.Fatal(err)
		}
		gw.Close()
		f.Close()
	}
	// Plain text, identical.
	write("a.frq", "chr1\t1\tA\n")
	write("b.frq", "chr1\t1\tA\n")
	if r := CompareOutputFiles(filepath.Join(dir, "a"), filepath.Join(dir, "b"), []string{".frq"}, matrix.ByteExact, similarityEpsilon); !r.Equal {
		t.Errorf("identical .frq should match: %+v", r)
	}
	// Gzipped, identical payload (different framing is irrelevant after decode).
	writeGz("a.bed.gz", "chr1\t0\t100\t5\n")
	writeGz("b.bed.gz", "chr1\t0\t100\t5\n")
	if r := CompareOutputFiles(filepath.Join(dir, "a"), filepath.Join(dir, "b"), []string{".bed.gz"}, matrix.ByteExact, similarityEpsilon); !r.Equal {
		t.Errorf("identical gzip payload should match: %+v", r)
	}
	// Mismatch.
	write("a.diff", "x\n")
	write("b.diff", "y\n")
	if r := CompareOutputFiles(filepath.Join(dir, "a"), filepath.Join(dir, "b"), []string{".diff"}, matrix.ByteExact, similarityEpsilon); r.Equal {
		t.Errorf("differing files should diverge")
	}
	// Presence mismatch (one side missing).
	write("a.only", "x\n")
	if r := CompareOutputFiles(filepath.Join(dir, "a"), filepath.Join(dir, "b"), []string{".only"}, matrix.ByteExact, similarityEpsilon); r.Equal {
		t.Errorf("missing-on-one-side should diverge")
	}
}
