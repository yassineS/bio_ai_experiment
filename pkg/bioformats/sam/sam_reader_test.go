package sam

import (
	"io"
	"strings"
	"testing"
)

const sampleSAM = `@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
@SQ	SN:chr2	LN:500
@RG	ID:rg1	SM:s1	LB:lib1
read1	99	chr1	100	60	5M	=	200	105	ACGTA	IIIII	NM:i:0	RG:Z:rg1
read2	147	chr1	200	60	5M	=	100	-105	TGCAT	!!!!!
read3	4	*	0	0	*	*	0	0	*	*
read4	0	chr2	1	30	3M2I	*	0	0	ACGTA	IIIII
`

func TestSAMReader(t *testing.T) {
	r, err := NewSAMReader(strings.NewReader(sampleSAM))
	if err != nil {
		t.Fatalf("NewSAMReader: %v", err)
	}
	hdr := r.Header()
	if len(hdr.Refs) != 2 {
		t.Fatalf("expected 2 refs, got %d", len(hdr.Refs))
	}
	var recs []*Record
	for {
		rec, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("Read: %v", err)
		}
		recs = append(recs, rec)
	}
	if len(recs) != 4 {
		t.Fatalf("expected 4 records, got %d", len(recs))
	}
	if recs[0].QName != "read1" || recs[0].Pos != 100 || recs[0].MapQ != 60 {
		t.Errorf("bad read1: %+v", recs[0])
	}
	if recs[0].Cigar.String() != "5M" {
		t.Errorf("read1 cigar: %q", recs[0].Cigar.String())
	}
	if len(recs[0].Aux) != 2 {
		t.Errorf("read1 aux: %+v", recs[0].Aux)
	}
	if recs[2].RName != "" || !recs[2].IsUnmapped() {
		t.Errorf("read3 should be unmapped with empty rname: %+v", recs[2])
	}
	if recs[2].Seq != "" || recs[2].Qual != nil {
		t.Errorf("read3 should have empty seq/qual: %+v", recs[2])
	}
}

func TestSAMReaderRoundTrip(t *testing.T) {
	r, err := NewSAMReader(strings.NewReader(sampleSAM))
	if err != nil {
		t.Fatalf("NewSAMReader: %v", err)
	}
	var out strings.Builder
	w := NewSAMWriter(&out)
	if err := w.WriteHeader(r.Header()); err != nil {
		t.Fatalf("WriteHeader: %v", err)
	}
	for {
		rec, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("Read: %v", err)
		}
		if err := w.Write(rec); err != nil {
			t.Fatalf("Write: %v", err)
		}
	}
	if err := w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	if got := out.String(); got != sampleSAM {
		t.Errorf("round-trip mismatch:\nin:\n%q\nout:\n%q", sampleSAM, got)
	}
}

func TestSAMReaderErrors(t *testing.T) {
	bad := "@HD\tVN:1.6\nshort\trecord\n"
	r, err := NewSAMReader(strings.NewReader(bad))
	if err != nil {
		t.Fatalf("NewSAMReader: %v", err)
	}
	if _, err := r.Read(); err == nil {
		t.Error("expected error for truncated record")
	}
}

func TestSAMReaderBadFlag(t *testing.T) {
	bad := "@HD\tVN:1.6\nrname\tNOPE\tchr1\t1\t0\t*\t*\t0\t0\t*\t*\n"
	r, _ := NewSAMReader(strings.NewReader(bad))
	if _, err := r.Read(); err == nil {
		t.Error("expected error for bad FLAG")
	}
}
