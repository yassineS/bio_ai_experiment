package samtools

import (
	"fmt"
	"reflect"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// TestBamTranslateAuxRepositions is a regression test for the merge aux-tag
// reposition bug the GIAB real-data test surfaced: upstream's bam_translate
// removes the RG (then PG) aux tag and re-appends it at the END of the aux list
// for every record, applying any collision rename. Our merge previously updated
// RG in place and never touched PG, so merged records diverged byte-for-byte.
func TestBamTranslateAuxRepositions(t *testing.T) {
	rec := &sam.Record{Aux: []sam.Aux{
		{Tag: "NM", Type: 'i', Value: int64(0)},
		{Tag: "RG", Type: 'Z', Value: "g1"},
		{Tag: "PG", Type: 'Z', Value: "p1"},
		{Tag: "XS", Type: 'A', Value: "+"},
	}}
	bamTranslateAux(rec, map[string]string{"g1": "g1-AAAA"}, map[string]string{"p1": "p1-BBBB"})

	var got []string
	for _, a := range rec.Aux {
		got = append(got, fmt.Sprintf("%s:%v", a.Tag, a.Value))
	}
	// RG and PG move to the end (RG before PG) with their values translated; the
	// other tags keep their original relative order.
	want := []string{"NM:0", "XS:+", "RG:g1-AAAA", "PG:p1-BBBB"}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("aux after bamTranslateAux = %v, want %v", got, want)
	}
}

// TestBamTranslateAuxIdentity checks that with no rename maps the tags are still
// moved to the end (identity translation), matching upstream which del+appends
// even when the translation is identity.
func TestBamTranslateAuxIdentity(t *testing.T) {
	rec := &sam.Record{Aux: []sam.Aux{
		{Tag: "RG", Type: 'Z', Value: "g1"},
		{Tag: "MD", Type: 'Z', Value: "101"},
		{Tag: "PG", Type: 'Z', Value: "p1"},
		{Tag: "NM", Type: 'i', Value: int64(1)},
	}}
	bamTranslateAux(rec, nil, nil)

	var got []string
	for _, a := range rec.Aux {
		got = append(got, a.Tag)
	}
	want := []string{"MD", "NM", "RG", "PG"}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("tag order after identity translate = %v, want %v", got, want)
	}
}
