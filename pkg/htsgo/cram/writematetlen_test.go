package cram

import (
	"bytes"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// mkMateRec builds a paired record with explicit mate fields, for the
// slice-span and equal-position TLEN tests. The quality is a fixed ramp so the
// record round-trips losslessly.
func mkMateRec(qname, rname string, pos int32, cigarStr, seq string, flag uint16, pnext int32, tlen int32) *sam.Record {
	var cig sam.Cigar
	if cigarStr != "" {
		c, err := sam.ParseCigar(cigarStr)
		if err != nil {
			panic(err)
		}
		cig = c
	}
	qual := make([]byte, len(seq))
	for i := range qual {
		qual[i] = byte(20 + i%30)
	}
	rnext := "*"
	if flag&sam.FlagPaired != 0 {
		rnext = "="
	}
	return &sam.Record{
		QName: qname, Flag: flag, RName: rname, Pos: int64(pos), MapQ: 40,
		Cigar: cig, RNext: rnext, PNext: int64(pnext), TLen: int64(tlen),
		Seq: seq, Qual: qual,
	}
}

// TestSliceSpanIncludesPlacedUnmapped is the BUG 1 regression guard: a placed
// (ref_id set, POS > 0) but FlagUnmapped read whose POS is below every mapped
// read in the slice must still contribute to the slice's alignment start. When
// it was excluded, an absolute AP stored for it fell before ref_seq_start and
// htslib rejected the slice with "Failure to decode slice".
func TestSliceSpanIncludesPlacedUnmapped(t *testing.T) {
	seq100 := strings.Repeat("ACGT", 25) // 100 bp
	records := []*sam.Record{
		mkMateRec("m1", "chr1", 2000, "100M", seq100, 0, 0, 0),
		// Placed-unmapped: FlagUnmapped set, but a real contig and POS below
		// the mapped reads. No CIGAR, so its reference end is its POS.
		mkMateRec("u1", "chr1", 1500, "", "ACGT", sam.FlagUnmapped, 0, 0),
		mkMateRec("m2", "chr1", 2500, "100M", seq100, 0, 0, 0),
	}

	start, span := sliceSpan(records)
	if start != 1500 {
		t.Errorf("sliceSpan start = %d, want 1500 (placed-unmapped read must set it)", start)
	}
	// Span must reach the rightmost mapped end (2500 + 100 - 1 = 2599).
	if wantSpan := int32(2599 - 1500 + 1); span != wantSpan {
		t.Errorf("sliceSpan span = %d, want %d", span, wantSpan)
	}
}

// TestEqualPosMateTLenOverride is the BUG 2 regression guard: for an equal-
// alignment-start same-reference primary mate pair, upstream htslib keeps the
// pair attached and re-derives each mate's TLEN sign on decode. The writer
// stores every record detached (verbatim TS), so it must pre-apply that
// re-derivation to stay byte-exact with an upstream-written CRAM's decoded
// output. Here the earlier record is the sole rightmost, so it takes the
// negative span and the later record the positive span — the reverse of the
// input BAM's leftmost-positive convention.
func TestEqualPosMateTLenOverride(t *testing.T) {
	// Mirror the real-data pair: two mates at the same POS whose soft clips give
	// them different reference ends. p (READ1) ends one base further right.
	pSeq := strings.Repeat("ACGT", 37) // 148 bp -> matches 13S135M query length.
	crSeq := strings.Repeat("TGCA", 37)
	// p: 13S135M at 1000 -> aend 1134. cr: 134M14S at 1000 -> aend 1133.
	p := mkMateRec("pair", "chr1", 1000, "13S135M", pSeq,
		sam.FlagPaired|sam.FlagProperPair|sam.FlagReverse|sam.FlagRead1, 1000, 135)
	cr := mkMateRec("pair", "chr1", 1000, "134M14S", crSeq,
		sam.FlagPaired|sam.FlagProperPair|sam.FlagMateReverse|sam.FlagRead2, 1000, -135)
	// A non-equal-position pair whose verbatim TLEN must be left untouched.
	nSeq := strings.Repeat("ACGT", 25)
	nA := mkMateRec("nonequal", "chr1", 3000, "100M", nSeq,
		sam.FlagPaired|sam.FlagProperPair|sam.FlagMateReverse|sam.FlagRead1, 3300, 400)
	nB := mkMateRec("nonequal", "chr1", 3300, "100M", nSeq,
		sam.FlagPaired|sam.FlagProperPair|sam.FlagReverse|sam.FlagRead2, 3000, -400)

	enc := &recordEncoder{version: VersionV30}
	enc.computeTLenOverrides([]*sam.Record{p, cr, nA, nB})

	if got, ok := enc.tlenOverride[p]; !ok || got != -135 {
		t.Errorf("earlier equal-pos mate TLEN override = %d (ok=%v), want -135", got, ok)
	}
	if got, ok := enc.tlenOverride[cr]; !ok || got != 135 {
		t.Errorf("later equal-pos mate TLEN override = %d (ok=%v), want 135", got, ok)
	}
	if _, ok := enc.tlenOverride[nA]; ok {
		t.Errorf("non-equal-pos mate A must not be overridden")
	}
	if _, ok := enc.tlenOverride[nB]; ok {
		t.Errorf("non-equal-pos mate B must not be overridden")
	}
}

// TestEqualPosMateTLenOverrideAttachGate verifies the attach gate: when the
// input TLEN signs do NOT match htslib's encode-time convention, htslib would
// store the pair detached and verbatim, so the writer must leave the TLEN
// untouched rather than re-sign it.
func TestEqualPosMateTLenOverrideAttachGate(t *testing.T) {
	pSeq := strings.Repeat("ACGT", 37)
	crSeq := strings.Repeat("TGCA", 37)
	// Same coordinates as the override test, but the input TLENs use the sign
	// the decoder would REDERIVE (p=-135, cr=+135). htslib's encode gate then
	// sees a mismatch against its encode-time sign and detaches -> verbatim.
	p := mkMateRec("pair", "chr1", 1000, "13S135M", pSeq,
		sam.FlagPaired|sam.FlagProperPair|sam.FlagReverse|sam.FlagRead1, 1000, -135)
	cr := mkMateRec("pair", "chr1", 1000, "134M14S", crSeq,
		sam.FlagPaired|sam.FlagProperPair|sam.FlagMateReverse|sam.FlagRead2, 1000, 135)

	enc := &recordEncoder{version: VersionV30}
	enc.computeTLenOverrides([]*sam.Record{p, cr})

	if _, ok := enc.tlenOverride[p]; ok {
		t.Errorf("pair with non-canonical input TLEN must be left verbatim (p overridden)")
	}
	if _, ok := enc.tlenOverride[cr]; ok {
		t.Errorf("pair with non-canonical input TLEN must be left verbatim (cr overridden)")
	}
}

// TestEqualPosMateRoundTrip encodes a slice holding both a placed-unmapped read
// below the mapped reads (BUG 1) and an equal-position mate pair (BUG 2), then
// decodes it and checks every record survives and the pair's TLEN matches
// upstream's re-derived convention (earlier = -span, later = +span).
func TestEqualPosMateRoundTrip(t *testing.T) {
	h := writerTestHeader()
	pSeq := strings.Repeat("ACGT", 37)
	crSeq := strings.Repeat("TGCA", 37)
	records := []*sam.Record{
		mkMateRec("pair", "chr1", 1000, "13S135M", pSeq,
			sam.FlagPaired|sam.FlagProperPair|sam.FlagReverse|sam.FlagRead1, 1000, 135),
		// Placed-unmapped read below the mapped pair.
		mkMateRec("orphan", "chr1", 800, "", "ACGTACGT", sam.FlagUnmapped, 0, 0),
		mkMateRec("pair", "chr1", 1000, "134M14S", crSeq,
			sam.FlagPaired|sam.FlagProperPair|sam.FlagMateReverse|sam.FlagRead2, 1000, -135),
	}

	var buf bytes.Buffer
	if err := WriteCRAM(&buf, h, records); err != nil {
		t.Fatalf("WriteCRAM: %v", err)
	}
	rr, err := NewRecordReader(bytes.NewReader(buf.Bytes()))
	if err != nil {
		t.Fatalf("NewRecordReader: %v", err)
	}
	out, err := rr.ReadAll()
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	if len(out) != len(records) {
		t.Fatalf("decoded %d records, want %d", len(out), len(records))
	}
	// Records preserve input order; index 0 is the earlier mate, 2 the later.
	if out[0].TLen != -135 {
		t.Errorf("earlier mate decoded TLEN = %d, want -135 (upstream re-derived sign)", out[0].TLen)
	}
	if out[2].TLen != 135 {
		t.Errorf("later mate decoded TLEN = %d, want 135 (upstream re-derived sign)", out[2].TLen)
	}
	if out[1].Pos != 800 {
		t.Errorf("placed-unmapped read POS = %d, want 800", out[1].Pos)
	}
}
