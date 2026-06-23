package samtools

import (
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// tagSeq extracts just the per-line @-tags from a header line slice, so a test
// can assert the grouped emission order without enumerating every field.
func tagSeq(lines []sam.HeaderLine) []string {
	out := make([]string, len(lines))
	for i, ln := range lines {
		out[i] = ln.Tag
	}
	return out
}

// TestRegroupMergedHeaderLines checks that regroupMergedHeaderLines reorders a
// merged header into upstream samtools' emission order — @HD, every @SQ, every
// @RG, every @PG, every @CO — when our union path leaves a renamed @RG after
// input 0's @PG block. It also confirms an already-grouped header is left
// untouched and the within-group order is preserved (stable partition).
func TestRegroupMergedHeaderLines(t *testing.T) {
	tests := []struct {
		name string
		in   []sam.HeaderLine
		want []string
	}{
		{
			name: "interleaved RG after PG is regrouped",
			// Mirrors the bug: input 0's header is HD, SQ, RG(1), PG(bwa),
			// PG(MarkDuplicates), PG(samtools); the union then appended a
			// renamed RG and renamed PGs, landing the renamed RG after PG.
			in: []sam.HeaderLine{
				{Tag: "HD"},
				{Tag: "SQ", Fields: []sam.HeaderField{{Tag: "SN", Value: "chr1"}}},
				{Tag: "SQ", Fields: []sam.HeaderField{{Tag: "SN", Value: "chr2"}}},
				{Tag: "RG", Fields: []sam.HeaderField{{Tag: "ID", Value: "1"}}},
				{Tag: "PG", Fields: []sam.HeaderField{{Tag: "ID", Value: "bwa"}}},
				{Tag: "PG", Fields: []sam.HeaderField{{Tag: "ID", Value: "MarkDuplicates"}}},
				// Renamed RG from input 1, appended after input 0's @PG block.
				{Tag: "RG", Fields: []sam.HeaderField{{Tag: "ID", Value: "1-055424A4"}}},
				{Tag: "PG", Fields: []sam.HeaderField{{Tag: "ID", Value: "bwa-3A2CCEF5"}}},
			},
			want: []string{"HD", "SQ", "SQ", "RG", "RG", "PG", "PG", "PG"},
		},
		{
			name: "already grouped is unchanged",
			in: []sam.HeaderLine{
				{Tag: "HD"},
				{Tag: "SQ"},
				{Tag: "RG"},
				{Tag: "PG"},
				{Tag: "CO"},
			},
			want: []string{"HD", "SQ", "RG", "PG", "CO"},
		},
		{
			name: "comment and unknown lines sort after PG",
			in: []sam.HeaderLine{
				{Tag: "HD"},
				{Tag: "CO"},
				{Tag: "SQ"},
				{Tag: "ZZ"}, // user-defined record sorts last
				{Tag: "PG"},
			},
			want: []string{"HD", "SQ", "PG", "CO", "ZZ"},
		},
	}
	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			got := tagSeq(regroupMergedHeaderLines(tc.in))
			if len(got) != len(tc.want) {
				t.Fatalf("length mismatch: got %v want %v", got, tc.want)
			}
			for i := range got {
				if got[i] != tc.want[i] {
					t.Fatalf("order mismatch at %d: got %v want %v", i, got, tc.want)
				}
			}
		})
	}
}

// TestRegroupMergedHeaderLinesPreservesWithinGroupOrder verifies the stable
// partition keeps each group's relative sequence, so input 0's groups precede
// later inputs' renamed groups exactly as upstream emits them.
func TestRegroupMergedHeaderLinesPreservesWithinGroupOrder(t *testing.T) {
	in := []sam.HeaderLine{
		{Tag: "HD"},
		{Tag: "RG", Fields: []sam.HeaderField{{Tag: "ID", Value: "1"}}},
		{Tag: "PG", Fields: []sam.HeaderField{{Tag: "ID", Value: "bwa"}}},
		{Tag: "PG", Fields: []sam.HeaderField{{Tag: "ID", Value: "MarkDuplicates"}}},
		{Tag: "RG", Fields: []sam.HeaderField{{Tag: "ID", Value: "1-AAAA"}}},
		{Tag: "PG", Fields: []sam.HeaderField{{Tag: "ID", Value: "bwa-BBBB"}}},
	}
	got := regroupMergedHeaderLines(in)
	wantRGIDs := []string{"1", "1-AAAA"}
	wantPGIDs := []string{"bwa", "MarkDuplicates", "bwa-BBBB"}
	var rgIDs, pgIDs []string
	for _, ln := range got {
		if len(ln.Fields) == 0 {
			continue
		}
		switch ln.Tag {
		case "RG":
			rgIDs = append(rgIDs, ln.Fields[0].Value)
		case "PG":
			pgIDs = append(pgIDs, ln.Fields[0].Value)
		}
	}
	if len(rgIDs) != len(wantRGIDs) {
		t.Fatalf("RG count: got %v want %v", rgIDs, wantRGIDs)
	}
	for i := range rgIDs {
		if rgIDs[i] != wantRGIDs[i] {
			t.Fatalf("RG order: got %v want %v", rgIDs, wantRGIDs)
		}
	}
	if len(pgIDs) != len(wantPGIDs) {
		t.Fatalf("PG count: got %v want %v", pgIDs, wantPGIDs)
	}
	for i := range pgIDs {
		if pgIDs[i] != wantPGIDs[i] {
			t.Fatalf("PG order: got %v want %v", pgIDs, wantPGIDs)
		}
	}
}
