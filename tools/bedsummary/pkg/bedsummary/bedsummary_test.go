package bedsummary

import (
	"bytes"
	"strings"
	"testing"
)

func TestCompute_Basic(t *testing.T) {
	in := strings.NewReader(`chr1	0	10
chr1	100	200
chr2	0	50
`)
	rows, err := Compute(in, Options{})
	if err != nil {
		t.Fatalf("Compute: %v", err)
	}
	// chr1, chr2, all
	if len(rows) != 3 {
		t.Fatalf("rows = %d, want 3", len(rows))
	}
	if rows[0].Chrom != "chr1" || rows[0].Count != 2 || rows[0].TotalLength != 110 {
		t.Errorf("chr1 row = %+v", rows[0])
	}
	if rows[0].MinLength != 10 || rows[0].MaxLength != 100 {
		t.Errorf("chr1 min/max wrong: %+v", rows[0])
	}
	if rows[0].MeanLength != 55.0 || rows[0].MedianLen != 55.0 {
		t.Errorf("chr1 mean/median wrong: %+v", rows[0])
	}
	if rows[1].Chrom != "chr2" || rows[1].Count != 1 || rows[1].TotalLength != 50 {
		t.Errorf("chr2 row = %+v", rows[1])
	}
	if rows[2].Chrom != "all" || rows[2].Count != 3 || rows[2].TotalLength != 160 {
		t.Errorf("all row = %+v", rows[2])
	}
}

func TestCompute_OddCountMedian(t *testing.T) {
	in := strings.NewReader(`chr1	0	1
chr1	0	5
chr1	0	10
`)
	rows, err := Compute(in, Options{SkipAll: true})
	if err != nil {
		t.Fatalf("Compute: %v", err)
	}
	if len(rows) != 1 {
		t.Fatalf("rows = %d, want 1 (SkipAll)", len(rows))
	}
	if rows[0].MedianLen != 5.0 {
		t.Errorf("median = %v, want 5", rows[0].MedianLen)
	}
}

func TestRun_HeaderAndFormatting(t *testing.T) {
	in := strings.NewReader("chr1\t0\t10\nchr1\t10\t30\n")
	var out bytes.Buffer
	if err := Run(in, &out, Options{}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	got := out.String()
	if !strings.HasPrefix(got, "chrom\tnum_ivls\t") {
		t.Errorf("missing header: %q", got)
	}
	// median of [10,20] is 15 -> integer formatting (no decimal point).
	if !strings.Contains(got, "\t15\n") {
		t.Errorf("expected integer median 15 in output, got: %q", got)
	}
}

func TestRun_NoHeader(t *testing.T) {
	in := strings.NewReader("chr1\t0\t10\n")
	var out bytes.Buffer
	if err := Run(in, &out, Options{NoHeader: true}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	if strings.HasPrefix(out.String(), "chrom\t") {
		t.Errorf("expected no header, got: %q", out.String())
	}
}

func TestRun_FractionalMedian(t *testing.T) {
	// lengths 1, 2, 3, 4 -> median = (2+3)/2 = 2.5
	in := strings.NewReader("chr1\t0\t1\nchr1\t0\t2\nchr1\t0\t3\nchr1\t0\t4\n")
	var out bytes.Buffer
	if err := Run(in, &out, Options{NoHeader: true, SkipAll: true}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	if !strings.Contains(out.String(), "2.500") {
		t.Errorf("expected fractional 2.500 median, got: %q", out.String())
	}
}

func TestCompute_EmptyInput(t *testing.T) {
	rows, err := Compute(strings.NewReader(""), Options{})
	if err != nil {
		t.Fatalf("Compute: %v", err)
	}
	if len(rows) != 0 {
		t.Errorf("rows = %d, want 0", len(rows))
	}
}

func TestCompute_PreservesInputChromOrder(t *testing.T) {
	// Note: chr2 appears before chr1.
	in := strings.NewReader("chr2\t0\t10\nchr1\t0\t5\n")
	rows, err := Compute(in, Options{SkipAll: true})
	if err != nil {
		t.Fatalf("Compute: %v", err)
	}
	if rows[0].Chrom != "chr2" || rows[1].Chrom != "chr1" {
		t.Errorf("chrom order = %v, %v; want chr2, chr1", rows[0].Chrom, rows[1].Chrom)
	}
}

func TestCompute_BadInput(t *testing.T) {
	in := strings.NewReader("chr1\tnotanumber\t100\n")
	if _, err := Compute(in, Options{}); err == nil {
		t.Errorf("expected error on bad start, got nil")
	}
}
