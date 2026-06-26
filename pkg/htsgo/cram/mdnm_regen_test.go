package cram

import (
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// cig builds a CIGAR from (length, op) pairs for the regeneration tests.
func cig(pairs ...[2]int) sam.Cigar {
	c := make(sam.Cigar, 0, len(pairs))
	for _, p := range pairs {
		c = append(c, sam.CigarOp(uint32(p[0])<<4|uint32(p[1])))
	}
	return c
}

// auxString returns a record's aux tags rendered as a SAM tab-joined
// string, in order, so a test can assert both the value and the position
// of MD/NM relative to the other tags.
func auxString(rec *sam.Record) string {
	out := ""
	for i, a := range rec.Aux {
		if i > 0 {
			out += "\t"
		}
		out += a.FormatSAM()
	}
	return out
}

// TestRegenerateMDNM is a hermetic regression test for the CRAM-decode
// MD/NM regeneration: a CRAM stores reads relative to the reference and
// drops MD:Z/NM:i, so upstream `samtools view -T ref file.cram`
// recomputes them. regenerateMDNM mirrors that. This test does not depend
// on the upstream binary; it constructs records directly and asserts the
// recomputed MD/NM values, their position in the aux block (spliced
// before a trailing data-series RG, exactly where htslib's cram_to_bam
// writes them), and the gating (unmapped records skipped, already-present
// tags left untouched, no reference suppresses regeneration).
func TestRegenerateMDNM(t *testing.T) {
	// Reference span starting at 1-based coordinate 1: positions 1..20.
	ref := []byte("ACGTACGTACGTACGTACGT")
	const refStart int32 = 1

	t.Run("perfect match before trailing RG", func(t *testing.T) {
		rec := &sam.Record{
			QName: "r1", Flag: 0, RName: "chr1", Pos: 1,
			Cigar: cig([2]int{8, sam.CigarMatch}),
			Seq:   "ACGTACGT",
			Aux: []sam.Aux{
				{Tag: "PG", Type: 'Z', Value: "x"},
				{Tag: "RG", Type: 'Z', Value: "rg1"},
			},
		}
		regenerateMDNM([]*sam.Record{rec}, nil, ref, refStart)
		// MD/NM must land between PG and the data-series RG.
		if got, want := auxString(rec), "PG:Z:x\tMD:Z:8\tNM:i:0\tRG:Z:rg1"; got != want {
			t.Fatalf("aux = %q, want %q", got, want)
		}
	})

	t.Run("single mismatch", func(t *testing.T) {
		// Reference at 1..8 is ACGTACGT; read has a C where ref has G at
		// the 3rd base (1-based 3). MD = "2G5", NM = 1.
		rec := &sam.Record{
			QName: "r2", Flag: 0, RName: "chr1", Pos: 1,
			Cigar: cig([2]int{8, sam.CigarMatch}),
			Seq:   "ACCTACGT",
		}
		regenerateMDNM([]*sam.Record{rec}, nil, ref, refStart)
		md, _ := rec.GetAux("MD")
		nm, _ := rec.GetAux("NM")
		if v, _ := md.String(); v != "2G5" {
			t.Errorf("MD = %q, want %q", v, "2G5")
		}
		if v, _ := nm.Int(); v != 1 {
			t.Errorf("NM = %d, want 1", v)
		}
		// With no trailing RG, MD then NM are appended at the very end.
		if got, want := auxString(rec), "MD:Z:2G5\tNM:i:1"; got != want {
			t.Errorf("aux = %q, want %q", got, want)
		}
	})

	t.Run("deletion", func(t *testing.T) {
		// 4M2D4M against ACGTACGTACGT: ref positions 5..6 (AC) are deleted.
		// MD = "4^AC4", NM = 2 (deletion length).
		rec := &sam.Record{
			QName: "r3", Flag: 0, RName: "chr1", Pos: 1,
			Cigar: cig([2]int{4, sam.CigarMatch}, [2]int{2, sam.CigarDeletion}, [2]int{4, sam.CigarMatch}),
			Seq:   "ACGTGTAC",
		}
		regenerateMDNM([]*sam.Record{rec}, nil, ref, refStart)
		md, _ := rec.GetAux("MD")
		nm, _ := rec.GetAux("NM")
		if v, _ := md.String(); v != "4^AC4" {
			t.Errorf("MD = %q, want %q", v, "4^AC4")
		}
		if v, _ := nm.Int(); v != 2 {
			t.Errorf("NM = %d, want 2", v)
		}
	})

	t.Run("insertion", func(t *testing.T) {
		// 4M2I4M: the 2 inserted bases consume read only and add to NM but
		// not MD. Ref ACGT then ACGT. MD = "8", NM = 2.
		rec := &sam.Record{
			QName: "r4", Flag: 0, RName: "chr1", Pos: 1,
			Cigar: cig([2]int{4, sam.CigarMatch}, [2]int{2, sam.CigarInsertion}, [2]int{4, sam.CigarMatch}),
			Seq:   "ACGTTTACGT",
		}
		regenerateMDNM([]*sam.Record{rec}, nil, ref, refStart)
		md, _ := rec.GetAux("MD")
		nm, _ := rec.GetAux("NM")
		if v, _ := md.String(); v != "8" {
			t.Errorf("MD = %q, want %q", v, "8")
		}
		if v, _ := nm.Int(); v != 2 {
			t.Errorf("NM = %d, want 2", v)
		}
	})

	t.Run("unmapped record left untouched", func(t *testing.T) {
		rec := &sam.Record{
			QName: "u1", Flag: sam.FlagUnmapped, RName: "*", Pos: 0,
			Seq: "ACGTACGT",
		}
		regenerateMDNM([]*sam.Record{rec}, nil, ref, refStart)
		if _, ok := rec.GetAux("MD"); ok {
			t.Error("unmapped record should not gain an MD tag")
		}
		if _, ok := rec.GetAux("NM"); ok {
			t.Error("unmapped record should not gain an NM tag")
		}
	})

	t.Run("existing MD and NM preserved", func(t *testing.T) {
		rec := &sam.Record{
			QName: "r5", Flag: 0, RName: "chr1", Pos: 1,
			Cigar: cig([2]int{8, sam.CigarMatch}),
			Seq:   "ACCTACGT", // a real mismatch, but stored tags win.
			Aux: []sam.Aux{
				{Tag: "MD", Type: 'Z', Value: "stored"},
				{Tag: "NM", Type: 'i', Value: int64(99)},
			},
		}
		regenerateMDNM([]*sam.Record{rec}, nil, ref, refStart)
		md, _ := rec.GetAux("MD")
		nm, _ := rec.GetAux("NM")
		if v, _ := md.String(); v != "stored" {
			t.Errorf("stored MD overwritten: got %q", v)
		}
		if v, _ := nm.Int(); v != 99 {
			t.Errorf("stored NM overwritten: got %d", v)
		}
	})

	t.Run("only NM stored regenerates MD only", func(t *testing.T) {
		// A record carrying NM but not MD must gain MD (recomputed) while
		// its stored NM is left untouched — htslib regenerates each tag
		// independently.
		rec := &sam.Record{
			QName: "r6", Flag: 0, RName: "chr1", Pos: 1,
			Cigar: cig([2]int{8, sam.CigarMatch}),
			Seq:   "ACCTACGT",
			Aux: []sam.Aux{
				{Tag: "NM", Type: 'i', Value: int64(7)},
			},
		}
		regenerateMDNM([]*sam.Record{rec}, nil, ref, refStart)
		md, hasMD := rec.GetAux("MD")
		if !hasMD {
			t.Fatal("MD not regenerated when only NM was stored")
		}
		if v, _ := md.String(); v != "2G5" {
			t.Errorf("MD = %q, want %q", v, "2G5")
		}
		nm, _ := rec.GetAux("NM")
		if v, _ := nm.Int(); v != 7 {
			t.Errorf("stored NM should be preserved, got %d", v)
		}
	})

	t.Run("nil reference suppresses regeneration", func(t *testing.T) {
		rec := &sam.Record{
			QName: "r7", Flag: 0, RName: "chr1", Pos: 1,
			Cigar: cig([2]int{8, sam.CigarMatch}),
			Seq:   "ACGTACGT",
		}
		regenerateMDNM([]*sam.Record{rec}, nil, nil, refStart)
		if _, ok := rec.GetAux("MD"); ok {
			t.Error("no reference must not regenerate MD")
		}
		if _, ok := rec.GetAux("NM"); ok {
			t.Error("no reference must not regenerate NM")
		}
	})

	// The cF "no MD"/"no NM" bits an embed_ref=2 (reduced/consensus) CRAM
	// stores per record mark that the source read carried no MD/NM. htslib
	// then leaves that record bare even though the (embedded) reference is
	// available (cram_decode.c:2117-2122). These cases pin that this decoder
	// reproduces the per-record suppression instead of fabricating tags —
	// the exact regression `samtools view` of an embed_ref=2 file exposes.
	t.Run("cF no-MD bit suppresses MD only", func(t *testing.T) {
		rec := &sam.Record{
			QName: "c1", Flag: 0, RName: "chr1", Pos: 1,
			Cigar: cig([2]int{8, sam.CigarMatch}),
			Seq:   "ACCTACGT", // a real mismatch — regeneration would yield MD:Z:2G5.
		}
		regenerateMDNM([]*sam.Record{rec}, []mdnmSuppress{{md: true}}, ref, refStart)
		if _, ok := rec.GetAux("MD"); ok {
			t.Error("cF no-MD bit must suppress MD regeneration")
		}
		// NM is not suppressed by the no-MD bit, so it is still recomputed.
		nm, ok := rec.GetAux("NM")
		if !ok {
			t.Fatal("NM should still be regenerated when only MD is suppressed")
		}
		if v, _ := nm.Int(); v != 1 {
			t.Errorf("NM = %d, want 1", v)
		}
	})

	t.Run("cF no-NM bit suppresses NM only", func(t *testing.T) {
		rec := &sam.Record{
			QName: "c2", Flag: 0, RName: "chr1", Pos: 1,
			Cigar: cig([2]int{8, sam.CigarMatch}),
			Seq:   "ACCTACGT",
		}
		regenerateMDNM([]*sam.Record{rec}, []mdnmSuppress{{nm: true}}, ref, refStart)
		md, ok := rec.GetAux("MD")
		if !ok {
			t.Fatal("MD should still be regenerated when only NM is suppressed")
		}
		if v, _ := md.String(); v != "2G5" {
			t.Errorf("MD = %q, want %q", v, "2G5")
		}
		if _, ok := rec.GetAux("NM"); ok {
			t.Error("cF no-NM bit must suppress NM regeneration")
		}
	})

	t.Run("cF both bits leave the record bare", func(t *testing.T) {
		// This is the embed_ref=2 fixture case: a perfectly matching read
		// whose source had neither MD nor NM. Before honouring cF the decoder
		// wrongly emitted MD:Z:8 / NM:i:0 here (the bug this test guards).
		rec := &sam.Record{
			QName: "c3", Flag: 0, RName: "chr1", Pos: 1,
			Cigar: cig([2]int{8, sam.CigarMatch}),
			Seq:   "ACGTACGT",
			Aux:   []sam.Aux{{Tag: "RG", Type: 'Z', Value: "rg1"}},
		}
		regenerateMDNM([]*sam.Record{rec}, []mdnmSuppress{{md: true, nm: true}}, ref, refStart)
		if got, want := auxString(rec), "RG:Z:rg1"; got != want {
			t.Fatalf("aux = %q, want %q (cF must leave the record bare)", got, want)
		}
	})

	t.Run("suppress shorter than recs falls back to regenerate", func(t *testing.T) {
		// A short/nil suppress slice must not panic and must regenerate the
		// records it does not cover (the common no-cF path passes nil).
		rec := &sam.Record{
			QName: "c4", Flag: 0, RName: "chr1", Pos: 1,
			Cigar: cig([2]int{8, sam.CigarMatch}),
			Seq:   "ACGTACGT",
		}
		regenerateMDNM([]*sam.Record{rec}, []mdnmSuppress{}, ref, refStart)
		if _, ok := rec.GetAux("MD"); !ok {
			t.Error("a record beyond the suppress slice must still regenerate MD")
		}
	})
}

// TestInsertBeforeTrailingRG asserts the splice position helper directly:
// MD/NM go immediately before a trailing data-series RG, else at the end.
func TestInsertBeforeTrailingRG(t *testing.T) {
	add := []sam.Aux{
		{Tag: "MD", Type: 'Z', Value: "8"},
		{Tag: "NM", Type: 'i', Value: int64(0)},
	}

	t.Run("before trailing RG", func(t *testing.T) {
		aux := []sam.Aux{
			{Tag: "PG", Type: 'Z', Value: "x"},
			{Tag: "RG", Type: 'Z', Value: "rg1"},
		}
		out := insertBeforeTrailingRG(aux, add)
		tags := []string{}
		for _, a := range out {
			tags = append(tags, a.Tag)
		}
		want := []string{"PG", "MD", "NM", "RG"}
		if len(tags) != len(want) {
			t.Fatalf("tags = %v, want %v", tags, want)
		}
		for i := range want {
			if tags[i] != want[i] {
				t.Fatalf("tags = %v, want %v", tags, want)
			}
		}
	})

	t.Run("appended when no trailing RG", func(t *testing.T) {
		aux := []sam.Aux{{Tag: "PG", Type: 'Z', Value: "x"}}
		out := insertBeforeTrailingRG(aux, add)
		tags := []string{}
		for _, a := range out {
			tags = append(tags, a.Tag)
		}
		want := []string{"PG", "MD", "NM"}
		for i := range want {
			if tags[i] != want[i] {
				t.Fatalf("tags = %v, want %v", tags, want)
			}
		}
	})

	t.Run("empty add returns aux unchanged", func(t *testing.T) {
		aux := []sam.Aux{{Tag: "RG", Type: 'Z', Value: "rg1"}}
		out := insertBeforeTrailingRG(aux, nil)
		if len(out) != 1 || out[0].Tag != "RG" {
			t.Fatalf("empty add altered aux: %v", out)
		}
	})
}
