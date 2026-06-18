package bedtag

import (
	"bytes"
	"io"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bed"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// buildBAM writes a minimal single-reference BAM with the given records.
func buildBAM(t *testing.T, recs []*sam.Record) []byte {
	t.Helper()
	hdr := &sam.Header{Refs: []sam.Reference{{Name: "chr1", Length: 1000}}}
	var buf bytes.Buffer
	w := sam.NewBAMWriter(&buf)
	if err := w.WriteHeader(hdr); err != nil {
		t.Fatal(err)
	}
	for _, r := range recs {
		if err := w.Write(r); err != nil {
			t.Fatal(err)
		}
	}
	if err := w.Close(); err != nil {
		t.Fatal(err)
	}
	return buf.Bytes()
}

func mapped(qname string, pos int32, reverse bool) *sam.Record {
	flag := uint16(0)
	if reverse {
		flag = sam.FlagReverse
	}
	return &sam.Record{
		QName: qname, Flag: flag, RName: "chr1", Pos: pos, MapQ: 60,
		Cigar: sam.Cigar{sam.CigarOp(10<<4 | 0)}, // 10M
		Seq:   "ACGTACGTAC", Qual: bytes.Repeat([]byte{30}, 10),
	}
}

// TestTagBAMLabels tags a BAM whose first read overlaps an annotation interval
// and whose second does not; only the overlapping read gets the YB tag.
func TestTagBAMLabels(t *testing.T) {
	in := buildBAM(t, []*sam.Record{
		mapped("hit", 100, false),  // [99,109) overlaps [95,150)
		mapped("miss", 500, false), // [499,509) overlaps nothing
	})
	anno := [][]*bed.Record{{
		{Chrom: "chr1", ChromStart: 95, ChromEnd: 150, Name: "f1", Score: 7, Strand: "+"},
	}}
	var out bytes.Buffer
	if _, err := TagBAM(bytes.NewReader(in), anno, &out, TagBAMOptions{Mode: TagModeLabels, Labels: []string{"iv"}}); err != nil {
		t.Fatal(err)
	}
	tags := readYB(t, out.Bytes())
	if tags["hit"] != "iv" {
		t.Errorf("hit YB = %q, want iv", tags["hit"])
	}
	if _, ok := tags["miss"]; ok {
		t.Errorf("miss should have no YB tag, got %q", tags["miss"])
	}
}

// TestTagBAMNamesAndScores checks the -names and -scores modes join the
// overlapping records' fields, and a custom -tag is honoured.
func TestTagBAMNamesAndScores(t *testing.T) {
	in := buildBAM(t, []*sam.Record{mapped("r", 100, false)})
	anno := [][]*bed.Record{{
		{Chrom: "chr1", ChromStart: 95, ChromEnd: 150, Name: "f1", Score: 7, Strand: "+"},
		{Chrom: "chr1", ChromStart: 98, ChromEnd: 160, Name: "f2", Score: 9, Strand: "+"},
	}}
	var names bytes.Buffer
	if _, err := TagBAM(bytes.NewReader(in), anno, &names, TagBAMOptions{Mode: TagModeNames, Tag: "YK"}); err != nil {
		t.Fatal(err)
	}
	if got := readTag(t, names.Bytes(), "YK")["r"]; got != "f1,f2" {
		t.Errorf("names YK = %q, want f1,f2", got)
	}
	var scores bytes.Buffer
	if _, err := TagBAM(bytes.NewReader(in), anno, &scores, TagBAMOptions{Mode: TagModeScores}); err != nil {
		t.Fatal(err)
	}
	if got := readYB(t, scores.Bytes())["r"]; got != "7,9" {
		t.Errorf("scores YB = %q, want 7,9", got)
	}
}

func readYB(t *testing.T, bamBytes []byte) map[string]string { return readTag(t, bamBytes, "YB") }

func readTag(t *testing.T, bamBytes []byte, tag string) map[string]string {
	t.Helper()
	br, err := sam.NewBAMReader(bytes.NewReader(bamBytes))
	if err != nil {
		t.Fatal(err)
	}
	out := map[string]string{}
	for {
		rec, err := br.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatal(err)
		}
		if a, ok := rec.GetAux(tag); ok {
			if s, ok := a.Value.(string); ok {
				out[rec.QName] = s
			}
		}
	}
	return out
}
