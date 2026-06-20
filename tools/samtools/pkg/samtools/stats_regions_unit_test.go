package samtools

import (
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// statsRegionTestHeader builds a minimal SAM header with two references, used
// by the binary-free positional-region unit tests. It needs no upstream binary
// and no on-disk file, so it runs with the reference_code submodules
// unpopulated.
func statsRegionTestHeader() *sam.Header {
	return &sam.Header{
		Refs: []sam.Reference{
			{Name: "chr1", Length: 100},
			{Name: "chr2", Length: 80},
		},
	}
}

// TestUnitRegionsFromSpecs covers the positional-region parser
// (regionsFromSpecs): 1-based inclusive interval construction, whole-chromosome
// and open-ended clamping to the header length, the covered-base count that
// drives the "bases inside the target" SN line, interval merging, and the
// unknown-reference error. It is a pure-function test with no upstream binary.
func TestUnitRegionsFromSpecs(t *testing.T) {
	hdr := statsRegionTestHeader()

	tests := []struct {
		name      string
		specs     []string
		wantByRef map[string][]regionInterval
		wantCount int64
	}{
		{
			name:  "single closed interval",
			specs: []string{"chr1:10-20"},
			wantByRef: map[string][]regionInterval{
				"chr1": {{Beg: 10, End: 20}},
			},
			wantCount: 11,
		},
		{
			name:  "whole chromosome clamps to header length",
			specs: []string{"chr2"},
			wantByRef: map[string][]regionInterval{
				"chr2": {{Beg: 1, End: 80}},
			},
			wantCount: 80,
		},
		{
			name:  "open-ended range clamps to header length",
			specs: []string{"chr1:50-"},
			wantByRef: map[string][]regionInterval{
				"chr1": {{Beg: 50, End: 100}},
			},
			wantCount: 51,
		},
		{
			name:  "single-coordinate range is one open-ended interval",
			specs: []string{"chr1:90"},
			wantByRef: map[string][]regionInterval{
				"chr1": {{Beg: 90, End: 100}},
			},
			wantCount: 11,
		},
		{
			name:  "two references combine their counts",
			specs: []string{"chr1:1-60", "chr2"},
			wantByRef: map[string][]regionInterval{
				"chr1": {{Beg: 1, End: 60}},
				"chr2": {{Beg: 1, End: 80}},
			},
			wantCount: 140,
		},
		{
			name:  "overlapping intervals on one ref are merged once",
			specs: []string{"chr1:10-20", "chr1:15-30"},
			wantByRef: map[string][]regionInterval{
				"chr1": {{Beg: 10, End: 30}},
			},
			wantCount: 21,
		},
		{
			name:  "end past the contig length clamps down",
			specs: []string{"chr1:90-200"},
			wantByRef: map[string][]regionInterval{
				"chr1": {{Beg: 90, End: 100}},
			},
			wantCount: 11,
		},
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			tr, err := regionsFromSpecs(tc.specs, hdr)
			if err != nil {
				t.Fatalf("regionsFromSpecs(%v): %v", tc.specs, err)
			}
			if tr.count != tc.wantCount {
				t.Errorf("count = %d, want %d", tr.count, tc.wantCount)
			}
			if len(tr.byRef) != len(tc.wantByRef) {
				t.Fatalf("byRef has %d refs, want %d", len(tr.byRef), len(tc.wantByRef))
			}
			for ref, want := range tc.wantByRef {
				got, ok := tr.byRef[ref]
				if !ok {
					t.Fatalf("missing ref %q", ref)
				}
				if len(got) != len(want) {
					t.Fatalf("ref %q: %d intervals, want %d (%v)", ref, len(got), len(want), got)
				}
				for i := range want {
					if got[i] != want[i] {
						t.Errorf("ref %q interval %d = %+v, want %+v", ref, i, got[i], want[i])
					}
				}
			}
		})
	}
}

// TestUnitRegionsFromSpecsUnknownRef confirms an unknown reference is a hard
// error, mirroring upstream's multi-region iterator refusing to build for a
// reference absent from the header.
func TestUnitRegionsFromSpecsUnknownRef(t *testing.T) {
	hdr := statsRegionTestHeader()
	if _, err := regionsFromSpecs([]string{"chrX:1-10"}, hdr); err == nil {
		t.Fatal("expected error for unknown reference, got nil")
	}
}

// TestUnitStatsIsInRegions exercises the per-read overlap predicate that gates
// every stats counter when positional regions (or -t) are active. The cases
// pin upstream stats.c's is_in_regions overlap rule: even a single-base overlap
// includes a read, a read entirely before/after the interval is excluded, and a
// read on a reference absent from the region set is excluded. No upstream
// binary or file is required.
func TestUnitStatsIsInRegions(t *testing.T) {
	mk := func(rname string, pos1 int32, length int32) *sam.Record {
		return &sam.Record{
			QName: "r",
			RName: rname,
			Pos:   int64(pos1),
			MapQ:  60,
			Cigar: sam.Cigar{sam.CigarOp(uint32(length)<<4 | sam.CigarMatch)},
			RNext: "*",
		}
	}

	tests := []struct {
		name   string
		region regionInterval // on chr1
		recRef string
		recPos int32
		recLen int32
		wantIn bool
	}{
		{"read fully inside", regionInterval{Beg: 10, End: 60}, "chr1", 20, 10, true},
		{"single-base overlap at start", regionInterval{Beg: 19, End: 60}, "chr1", 10, 10, true}, // read 10..19
		{"read ends one before region begins", regionInterval{Beg: 20, End: 60}, "chr1", 10, 10, false},
		{"read begins one after region ends", regionInterval{Beg: 1, End: 9}, "chr1", 10, 10, false},
		{"read on absent reference", regionInterval{Beg: 1, End: 100}, "chr2", 10, 10, false},
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			c := newStatsCounters()
			c.IsSorted = 1
			c.regions = &targetRegions{byRef: map[string][]regionInterval{
				"chr1": {tc.region},
			}}
			c.regCursor = map[string]int{}
			rec := mk(tc.recRef, tc.recPos, tc.recLen)
			in, err := c.isInRegions(rec)
			if err != nil {
				t.Fatalf("isInRegions: %v", err)
			}
			if in != tc.wantIn {
				t.Errorf("isInRegions = %v, want %v", in, tc.wantIn)
			}
		})
	}
}
