package runner

import "testing"

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
	r := CompareSimilarity(a, b)
	if !r.Equal {
		t.Errorf("expected within-tolerance equal, got %+v", r)
	}
	if r.MaxDeviation == 0 {
		t.Errorf("expected non-zero recorded deviation")
	}

	d := CompareSimilarity([]byte("chr1\t1.0\n"), []byte("chr2\t1.0\n"))
	if d.Equal {
		t.Errorf("expected non-numeric field mismatch to diverge")
	}

	e := CompareSimilarity([]byte("chr1\t1.0\n"), []byte("chr1\t2.0\n"))
	if e.Equal {
		t.Errorf("expected out-of-tolerance numeric to diverge")
	}
}
