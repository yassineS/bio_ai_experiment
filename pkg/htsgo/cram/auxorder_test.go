package cram

import (
	"bytes"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

func tagsEqual(got, want []string) bool {
	if len(got) != len(want) {
		return false
	}
	for i := range want {
		if got[i] != want[i] {
			return false
		}
	}
	return true
}

// TestAuxForEncoding pins the encoder-side tag migration that makes a CRAM
// v2/v3 file decode with MD, NM and RG at the tail (in that order), matching
// upstream htslib's cram_encode.c routing those tags through dedicated
// handling. CRAM v4 keeps them inline as dictionary placeholders, so reorder
// is false there and the aux slice is returned untouched.
func TestAuxForEncoding(t *testing.T) {
	mapped := func(aux []sam.Aux) *sam.Record {
		return &sam.Record{QName: "r", RName: "chr1", Pos: 1, Cigar: cig([2]int{8, sam.CigarMatch}), Seq: "ACGTACGT", Aux: aux}
	}

	tests := []struct {
		name    string
		rec     *sam.Record
		reorder bool
		want    []string
	}{
		{
			name: "inline MD NM RG migrate to tail in order",
			rec: mapped([]sam.Aux{
				{Tag: "X0", Type: 'i', Value: int64(1)},
				{Tag: "MD", Type: 'Z', Value: "8"},
				{Tag: "PG", Type: 'Z', Value: "p"},
				{Tag: "RG", Type: 'Z', Value: "rg1"},
				{Tag: "AM", Type: 'i', Value: int64(37)},
				{Tag: "NM", Type: 'i', Value: int64(0)},
				{Tag: "XT", Type: 'A', Value: "U"},
			}),
			reorder: true,
			want:    []string{"X0", "PG", "AM", "XT", "MD", "NM", "RG"},
		},
		{
			name: "already-ordered tail is left untouched",
			rec: mapped([]sam.Aux{
				{Tag: "X0", Type: 'i', Value: int64(1)},
				{Tag: "MD", Type: 'Z', Value: "8"},
				{Tag: "NM", Type: 'i', Value: int64(0)},
				{Tag: "RG", Type: 'Z', Value: "rg1"},
			}),
			reorder: true,
			want:    []string{"X0", "MD", "NM", "RG"},
		},
		{
			name: "RG alone migrates to the tail",
			rec: mapped([]sam.Aux{
				{Tag: "RG", Type: 'Z', Value: "rg1"},
				{Tag: "X0", Type: 'i', Value: int64(1)},
			}),
			reorder: true,
			want:    []string{"X0", "RG"},
		},
		{
			name: "unmapped read keeps inline MD and NM, only RG migrates",
			rec: &sam.Record{
				QName: "u", Flag: sam.FlagUnmapped, RName: "*", Seq: "ACGT",
				Aux: []sam.Aux{
					{Tag: "MD", Type: 'Z', Value: "4"},
					{Tag: "RG", Type: 'Z', Value: "rg1"},
					{Tag: "NM", Type: 'i', Value: int64(0)},
					{Tag: "X0", Type: 'i', Value: int64(1)},
				},
			},
			reorder: true,
			want:    []string{"MD", "NM", "X0", "RG"},
		},
		{
			name: "v4 (reorder off) returns the aux unchanged",
			rec: mapped([]sam.Aux{
				{Tag: "X0", Type: 'i', Value: int64(1)},
				{Tag: "MD", Type: 'Z', Value: "8"},
				{Tag: "RG", Type: 'Z', Value: "rg1"},
				{Tag: "NM", Type: 'i', Value: int64(0)},
			}),
			reorder: false,
			want:    []string{"X0", "MD", "RG", "NM"},
		},
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			got := auxForEncoding(tc.rec, tc.reorder)
			names := make([]string, 0, len(got))
			for _, a := range got {
				names = append(names, a.Tag)
			}
			if !tagsEqual(names, tc.want) {
				t.Fatalf("auxForEncoding tags = %v, want %v", names, tc.want)
			}
		})
	}
}

// TestWriteCRAMInlineMDNMReordered is the end-to-end regression for the CRAM
// aux-tag ordering bug: a mapped read whose source aux carries MD and NM
// *inline* (not at the tail) must, after a CRAM v3 encode/decode round-trip,
// emit MD, NM and RG as the final three tags in that order — matching what
// upstream `samtools view` produces. Before the fix the writer kept the inline
// MD/NM (and the inline RG) in their original positions, so the decoded order
// diverged from upstream even though the records were otherwise identical.
func TestWriteCRAMInlineMDNMReordered(t *testing.T) {
	h := writerTestHeader()
	rec := mkRec("read1", "chr1", 100, "8M", "ACGTACGT")
	rec.Aux = []sam.Aux{
		{Tag: "X0", Type: 'i', Value: int64(1)},
		{Tag: "MD", Type: 'Z', Value: "8"},
		{Tag: "PG", Type: 'Z', Value: "MarkDuplicates"},
		{Tag: "RG", Type: 'Z', Value: "rg1"},
		{Tag: "AM", Type: 'i', Value: int64(37)},
		{Tag: "NM", Type: 'i', Value: int64(0)},
		{Tag: "XT", Type: 'A', Value: "U"},
	}

	out := roundTrip(t, h, []*sam.Record{rec})
	if len(out) != 1 {
		t.Fatalf("got %d records, want 1", len(out))
	}
	got := auxTags(out[0].Aux)
	want := []string{"X0", "PG", "AM", "XT", "MD", "NM", "RG"}
	if !tagsEqual(got, want) {
		t.Fatalf("decoded aux order = %v, want %v", got, want)
	}

	// The values must survive the reorder unchanged.
	md, ok := out[0].GetAux("MD")
	if !ok {
		t.Fatal("MD tag lost in round-trip")
	}
	if v, _ := md.String(); v != "8" {
		t.Errorf("MD = %q, want %q", v, "8")
	}
	nm, ok := out[0].GetAux("NM")
	if !ok {
		t.Fatal("NM tag lost in round-trip")
	}
	if v, _ := nm.Int(); v != 0 {
		t.Errorf("NM = %d, want 0", v)
	}
	rg, ok := out[0].GetAux("RG")
	if !ok {
		t.Fatal("RG tag lost in round-trip")
	}
	if v, _ := rg.String(); v != "rg1" {
		t.Errorf("RG = %q, want %q", v, "rg1")
	}
}

// TestWriteCRAMV4KeepsInlineMDNM guards the version gate: a CRAM v4 write must
// NOT migrate inline MD/NM/RG to the tail (htslib v4 keeps them in the tag
// dictionary at their inline position), so the decoded order matches the
// source order verbatim. This pins the boundary between the v2/v3 reorder and
// the v4 in-place behaviour.
func TestWriteCRAMV4KeepsInlineMDNM(t *testing.T) {
	h := writerTestHeader()
	rec := mkRec("read1", "chr1", 100, "8M", "ACGTACGT")
	rec.Aux = []sam.Aux{
		{Tag: "X0", Type: 'i', Value: int64(1)},
		{Tag: "MD", Type: 'Z', Value: "8"},
		{Tag: "RG", Type: 'Z', Value: "rg1"},
		{Tag: "NM", Type: 'i', Value: int64(0)},
	}

	var buf bytes.Buffer
	if err := WriteCRAMVersion(&buf, h, []*sam.Record{rec}, VersionV40); err != nil {
		t.Fatalf("WriteCRAMVersion v4.0: %v", err)
	}
	rr, err := NewRecordReader(&buf)
	if err != nil {
		t.Fatalf("NewRecordReader: %v", err)
	}
	out, err := rr.ReadAll()
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if len(out) != 1 {
		t.Fatalf("got %d records, want 1", len(out))
	}
	got := auxTags(out[0].Aux)
	want := []string{"X0", "MD", "RG", "NM"}
	if !tagsEqual(got, want) {
		t.Fatalf("v4 decoded aux order = %v, want %v (v4 must not reorder)", got, want)
	}
}
