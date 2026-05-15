package bed

import (
	"bytes"
	"io"
	"strings"
	"testing"
)

func TestBEDPEReader_Basic(t *testing.T) {
	in := strings.Join([]string{
		"# header comment",
		"chr1\t10\t20\tchr1\t100\t200\tpairA\t30\t+\t-",
		"chr2\t50\t60\tchr3\t500\t600\tpairB\t.\t-\t+\textra1\textra2",
		"",
		"track name=foo",
		"chr4\t0\t10\tchr4\t20\t30\tpairC\t10\t+\t+",
	}, "\n")
	r := NewBEDPEReader(strings.NewReader(in))
	got, err := r.ReadAll()
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if len(got) != 3 {
		t.Fatalf("expected 3 records, got %d", len(got))
	}
	if got[0].Chrom1 != "chr1" || got[0].Start1 != 10 || got[0].End1 != 20 ||
		got[0].Chrom2 != "chr1" || got[0].Start2 != 100 || got[0].End2 != 200 ||
		got[0].Name != "pairA" || got[0].Score != "30" ||
		got[0].Strand1 != "+" || got[0].Strand2 != "-" {
		t.Errorf("record 0 mismatch: %+v", got[0])
	}
	if len(got[1].Extra) != 2 || got[1].Extra[0] != "extra1" || got[1].Extra[1] != "extra2" {
		t.Errorf("record 1 extras mismatch: %+v", got[1].Extra)
	}
}

func TestBEDPEReader_TooFewFields(t *testing.T) {
	r := NewBEDPEReader(strings.NewReader("chr1\t1\t2\tchr1\t3\t4\tn\t.\t+"))
	if _, err := r.Read(); err == nil {
		t.Fatal("expected error for <10 fields")
	}
}

func TestBEDPEReader_BadInt(t *testing.T) {
	r := NewBEDPEReader(strings.NewReader("chr1\t1\tNaN\tchr1\t3\t4\tn\t.\t+\t-"))
	if _, err := r.Read(); err == nil {
		t.Fatal("expected parse error for non-numeric end1")
	}
}

func TestBEDPEReader_EOFEmpty(t *testing.T) {
	r := NewBEDPEReader(strings.NewReader(""))
	if _, err := r.Read(); err != io.EOF {
		t.Fatalf("expected io.EOF on empty input, got %v", err)
	}
}

func TestBEDPEWriter_Roundtrip(t *testing.T) {
	rec := &BEDPE{
		Chrom1: "chr1", Start1: 10, End1: 20,
		Chrom2: "chr2", Start2: 100, End2: 200,
		Name: "p", Score: "0", Strand1: "+", Strand2: "-",
		Extra: []string{"ext1"},
	}
	var buf bytes.Buffer
	w := NewBEDPEWriter(&buf)
	if err := w.Write(rec); err != nil {
		t.Fatalf("Write: %v", err)
	}
	if err := w.Flush(); err != nil {
		t.Fatalf("Flush: %v", err)
	}
	want := "chr1\t10\t20\tchr2\t100\t200\tp\t0\t+\t-\text1\n"
	if got := buf.String(); got != want {
		t.Errorf("BEDPE roundtrip mismatch:\n got=%q\nwant=%q", got, want)
	}
}

func TestBEDPE_EndRecords(t *testing.T) {
	rec := &BEDPE{
		Chrom1: "chr1", Start1: 10, End1: 20, Strand1: "+",
		Chrom2: "chr2", Start2: 100, End2: 200, Strand2: "-",
		Name: "p",
	}
	e1 := rec.End1Record()
	e2 := rec.End2Record()
	if e1.Chrom != "chr1" || e1.ChromStart != 10 || e1.ChromEnd != 20 || e1.Strand != "+" || e1.Name != "p" {
		t.Errorf("End1Record mismatch: %+v", e1)
	}
	if e2.Chrom != "chr2" || e2.ChromStart != 100 || e2.ChromEnd != 200 || e2.Strand != "-" || e2.Name != "p" {
		t.Errorf("End2Record mismatch: %+v", e2)
	}
}
