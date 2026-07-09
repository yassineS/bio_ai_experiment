package mdnm

import (
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// TestComputeIntoMatchesComputeWithReusedBuffer verifies that ComputeInto
// (the pooled-buffer variant the CRAM decoder uses for per-slice MD/NM
// regeneration) returns byte-identical MD:Z and NM:i to Compute across a batch
// of records while reusing one scratch buffer — exactly the reuse pattern
// regenerateMDNM applies. If the buffer were not truncated before each record,
// or the result string were to alias the reused backing array, a later record
// would corrupt an earlier record's MD; this test would then observe a mismatch
// or a mutated earlier result.
func TestComputeIntoMatchesComputeWithReusedBuffer(t *testing.T) {
	const ref = "ACGTACGTACGTACGTACGTACGTACGT"
	r := []byte(ref)

	type recSpec struct {
		pos   int
		cigar string
		seq   string
	}
	// A spread of shapes with different MD lengths so the reused buffer both
	// grows and later shrinks — the case most likely to leak stale bytes.
	specs := []recSpec{
		{pos: 1, cigar: "8M", seq: "ACGTACGT"},      // perfect: MD "8"
		{pos: 1, cigar: "8M", seq: "AAGTAAGT"},      // two mismatches: MD "1C4C1"
		{pos: 1, cigar: "4M2D4M", seq: "ACGTATAC"},  // deletion: MD "4^AC0G3"
		{pos: 5, cigar: "4M", seq: "ACGT"},          // short: MD "4"
		{pos: 1, cigar: "12M", seq: "ACGTACGTACGT"}, // longer perfect: MD "12"
		{pos: 1, cigar: "3M1I3M", seq: "ACGTACG"},   // insertion
		{pos: 2, cigar: "6M", seq: "CGTNCG"},        // N base mismatch
		{pos: 1, cigar: "8M", seq: "ACGTACGT"},      // back to short after long
	}

	// Compute the reference (unpooled) answers first.
	type want struct {
		md string
		nm int
	}
	wants := make([]want, len(specs))
	for i, s := range specs {
		rec := &sam.Record{Pos: int64(s.pos), Cigar: mustCigar(t, s.cigar), Seq: s.seq}
		md, nm := Compute(rec, r, 0)
		wants[i] = want{md, nm}
	}

	// Now replay through ComputeInto with a single reused buffer, exactly as
	// regenerateMDNM does, and confirm each result matches and prior results
	// stay intact.
	var buf []byte
	got := make([]string, len(specs))
	for i, s := range specs {
		rec := &sam.Record{Pos: int64(s.pos), Cigar: mustCigar(t, s.cigar), Seq: s.seq}
		md, nm, b := ComputeInto(rec, r, 0, buf)
		buf = b
		got[i] = md
		if md != wants[i].md {
			t.Fatalf("record %d: ComputeInto MD = %q, want %q", i, md, wants[i].md)
		}
		if nm != wants[i].nm {
			t.Fatalf("record %d: ComputeInto NM = %d, want %d", i, nm, wants[i].nm)
		}
		// Every earlier result must be unchanged (no aliasing of the reused
		// backing array into a retained string).
		for j := 0; j < i; j++ {
			if got[j] != wants[j].md {
				t.Fatalf("record %d encode corrupted earlier record %d: MD = %q, want %q",
					i, j, got[j], wants[j].md)
			}
		}
	}
}
