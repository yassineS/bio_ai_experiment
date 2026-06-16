package bedtobam

import (
	"bytes"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// roundtrip runs Run over bed against the given genome and decodes the BAM back
// into records and a header using the sam package, so unit tests can assert on
// the decoded content without an external decoder.
func roundtrip(t *testing.T, genome, bed string, opts Options) (*sam.Header, []*sam.Record) {
	t.Helper()
	g, err := ReadGenome(strings.NewReader(genome))
	if err != nil {
		t.Fatalf("ReadGenome: %v", err)
	}
	if opts.GenomeFileName == "" {
		opts.GenomeFileName = "genome.txt"
	}
	var buf bytes.Buffer
	if _, err := Run(strings.NewReader(bed), &buf, g, opts); err != nil {
		t.Fatalf("Run: %v", err)
	}
	br, err := sam.NewReader(&buf)
	if err != nil {
		t.Fatalf("open BAM: %v", err)
	}
	var recs []*sam.Record
	for {
		rec, err := br.Read()
		if err != nil {
			break
		}
		recs = append(recs, rec)
	}
	return br.Header(), recs
}

func TestHeader_Format(t *testing.T) {
	h, _ := roundtrip(t, "1\t3000\n2\t2000\n", "1\t10\t20\tr\t0\t+\n", Options{MapQ: 255, GenomeFileName: "/g.txt"})
	text := h.Text()
	want := "@HD\tVN:1.0\tSO:unsorted\n" +
		"@PG\tID:BEDTools_bedToBam\tVN:Vv2.31.1\n" +
		"@SQ\tSN:1\tAS:/g.txt\tLN:3000\n" +
		"@SQ\tSN:2\tAS:/g.txt\tLN:2000\n"
	if text != want {
		t.Errorf("header:\n got %q\nwant %q", text, want)
	}
}

func TestPlainRecord(t *testing.T) {
	_, recs := roundtrip(t, "1\t3000\n", "1\t1000\t2000\tread_name\t255\t+\n", Options{MapQ: 255})
	if len(recs) != 1 {
		t.Fatalf("want 1 record, got %d", len(recs))
	}
	r := recs[0]
	if r.QName != "read_name" || r.RName != "1" || r.Pos != 1001 || r.MapQ != 255 {
		t.Errorf("record = %+v", r)
	}
	if r.Cigar.String() != "1000M" {
		t.Errorf("cigar = %s, want 1000M", r.Cigar.String())
	}
	if r.Flag != 0 {
		t.Errorf("flag = %d, want 0", r.Flag)
	}
}

func TestReverseStrandFlag(t *testing.T) {
	_, recs := roundtrip(t, "1\t3000\n", "1\t10\t50\tr\t0\t-\n", Options{MapQ: 0})
	if recs[0].Flag&sam.FlagReverse == 0 {
		t.Errorf("expected reverse flag, got %d", recs[0].Flag)
	}
}

func TestBED12_Cigar(t *testing.T) {
	bed := "1\t100\t300\tblk\t40\t+\t100\t300\t255,0,0\t3\t10,10,10\t0,20,40\n"
	_, recs := roundtrip(t, "1\t3000\n", bed, Options{MapQ: 255, BED12: true})
	if len(recs) != 1 {
		t.Fatalf("want 1 record, got %d", len(recs))
	}
	// blocks at 0,20,40 sizes 10 -> 10M10N10M10N10M.
	if got := recs[0].Cigar.String(); got != "10M10N10M10N10M" {
		t.Errorf("cigar = %s, want 10M10N10M10N10M", got)
	}
	if recs[0].Pos != 101 {
		t.Errorf("pos = %d, want 101", recs[0].Pos)
	}
}

func TestBED12_LeadingSkip(t *testing.T) {
	// First block starts at offset 5 -> leading 5N.
	bed := "1\t100\t200\tblk\t40\t+\t100\t200\t255,0,0\t1\t10\t5\n"
	_, recs := roundtrip(t, "1\t3000\n", bed, Options{MapQ: 255, BED12: true})
	if got := recs[0].Cigar.String(); got != "5N10M" {
		t.Errorf("cigar = %s, want 5N10M", got)
	}
}

func TestUnknownChrom_MapsToFirstRef(t *testing.T) {
	_, recs := roundtrip(t, "1\t3000\n2\t2000\n", "chrX\t10\t20\tr\t0\t+\n", Options{MapQ: 0})
	if recs[0].RName != "1" {
		t.Errorf("unknown chrom should map to first ref, got %q", recs[0].RName)
	}
}

func TestMissingName_IsError(t *testing.T) {
	g, _ := ReadGenome(strings.NewReader("1\t3000\n"))
	var buf bytes.Buffer
	_, err := Run(strings.NewReader("1\t10\t20\n"), &buf, g, Options{MapQ: 0})
	if err == nil || !strings.Contains(err.Error(), "without name") {
		t.Errorf("expected 'without name' error, got %v", err)
	}
}

func TestMapqRange(t *testing.T) {
	g, _ := ReadGenome(strings.NewReader("1\t3000\n"))
	var buf bytes.Buffer
	if _, err := Run(strings.NewReader("1\t10\t20\tr\t0\t+\n"), &buf, g, Options{MapQ: 300}); err == nil {
		t.Errorf("expected MAPQ range error")
	}
}
