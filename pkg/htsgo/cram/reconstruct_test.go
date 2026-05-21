package cram

import (
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// refFreeHeader builds a minimal reference-free compression header that
// names the small set of data-series encodings the synthetic feature
// tests draw from. Every series is a degenerate construct so the tests
// can drive the decoder without a full CRAM file.
func refFreeHeader() *CompressionHeader {
	return &CompressionHeader{
		Preservation: PreservationMap{ReferenceRequired: false},
		DataSeries:   map[dataSeriesKey]*Encoding{},
		Tags:         map[tagKey]*Encoding{},
	}
}

// newReconDecoder builds a recordDecoder over a hand-built series source
// for the reconstruct tests.
func newReconDecoder(refRequired bool) *recordDecoder {
	h := refFreeHeader()
	h.Preservation.ReferenceRequired = refRequired
	return &recordDecoder{
		h:     h,
		src:   &SeriesSource{s: newTestSource(nil, nil)},
		slice: &SliceHeader{},
	}
}

// TestReconstructMappedMatchOnly checks the simplest reconstruction: a
// single base-stretch feature spanning the whole read yields a pure
// match CIGAR and copies the bases verbatim.
func TestReconstructMappedMatchOnly(t *testing.T) {
	rd := newReconDecoder(false)
	feats := []readFeature{{code: featBases, pos: 1, bases: []byte("ACGTAC")}}
	seq, qual, cig, err := rd.reconstructMapped(feats, 6, 1)
	if err != nil {
		t.Fatalf("reconstructMapped: %v", err)
	}
	if string(seq) != "ACGTAC" {
		t.Errorf("seq = %q, want ACGTAC", seq)
	}
	if cig.String() != "6M" {
		t.Errorf("cigar = %q, want 6M", cig.String())
	}
	for i, q := range qual {
		if q != 0xff {
			t.Errorf("qual[%d] = %d, want 0xff (unknown)", i, q)
		}
	}
}

// TestReconstructMappedAllFeatures drives every read-feature code that
// contributes to a CIGAR, checking the reconstructed sequence and CIGAR.
func TestReconstructMappedAllFeatures(t *testing.T) {
	rd := newReconDecoder(false)
	// A read of length 8: 2 stretch bases (M), 1 insert base (I),
	// 1 single base (M), a quality score, a deletion (D), a soft clip
	// of 2 (S), padding (P), and a final single base (M).
	feats := []readFeature{
		{code: featBases, pos: 1, bases: []byte("AC")},
		{code: featInsertBase, pos: 3, base: 'G'},
		{code: featBase, pos: 4, base: 'T', quality: 20},
		{code: featQualityScore, pos: 1, quality: 30},
		{code: featDeletion, pos: 5, length: 3},
		{code: featSoftClip, pos: 5, bases: []byte("NN")},
		{code: featPadding, pos: 7, length: 2},
		{code: featInsertion, pos: 7, bases: []byte("CC")},
	}
	seq, qual, cig, err := rd.reconstructMapped(feats, 8, 1)
	if err != nil {
		t.Fatalf("reconstructMapped: %v", err)
	}
	if string(seq) != "ACGTNNCC" {
		t.Errorf("seq = %q, want ACGTNNCC", seq)
	}
	if cig.String() != "2M1I1M3D2S2P2I" {
		t.Errorf("cigar = %q, want 2M1I1M3D2S2P2I", cig.String())
	}
	if qual[0] != 30 {
		t.Errorf("qual[0] = %d, want 30 (from the Q feature)", qual[0])
	}
	if qual[3] != 20 {
		t.Errorf("qual[3] = %d, want 20 (from the B feature)", qual[3])
	}
}

// TestReconstructQualityStretch checks the "q" feature, a stretch of
// quality scores annotating already-placed bases.
func TestReconstructQualityStretch(t *testing.T) {
	rd := newReconDecoder(false)
	feats := []readFeature{
		{code: featBases, pos: 1, bases: []byte("ACGT")},
		{code: featScores, pos: 1, bases: []byte{10, 11, 12, 13}},
	}
	_, qual, _, err := rd.reconstructMapped(feats, 4, 1)
	if err != nil {
		t.Fatalf("reconstructMapped: %v", err)
	}
	for i, want := range []byte{10, 11, 12, 13} {
		if qual[i] != want {
			t.Errorf("qual[%d] = %d, want %d", i, qual[i], want)
		}
	}
}

// TestReconstructSubstitutionNeedsReference checks that a substitution
// feature fills the base with 'N' and flags the record as needing an
// external reference.
func TestReconstructSubstitutionNeedsReference(t *testing.T) {
	rd := newReconDecoder(true)
	feats := []readFeature{
		{code: featBases, pos: 1, bases: []byte("AC")},
		{code: featSubst, pos: 3, substCode: 1},
		{code: featBases, pos: 4, bases: []byte("GT")},
	}
	seq, _, cig, err := rd.reconstructMapped(feats, 5, 1)
	if err != nil {
		t.Fatalf("reconstructMapped: %v", err)
	}
	if string(seq) != "ACNGT" {
		t.Errorf("seq = %q, want ACNGT", seq)
	}
	if cig.String() != "5M" {
		t.Errorf("cigar = %q, want 5M", cig.String())
	}
	if !rd.needsReference {
		t.Error("a substitution feature must set needsReference")
	}
}

// TestReconstructUncoveredBase checks that a read position no feature
// covers is filled with 'N' and flags the reference requirement.
func TestReconstructUncoveredBase(t *testing.T) {
	rd := newReconDecoder(true)
	// Two base-stretches leave a one-base gap at read position 3, which
	// in a reference-backed file the reference would supply.
	feats := []readFeature{
		{code: featBases, pos: 1, bases: []byte("AC")},
		{code: featBases, pos: 4, bases: []byte("GT")},
	}
	seq, _, cig, err := rd.reconstructMapped(feats, 5, 1)
	if err != nil {
		t.Fatalf("reconstructMapped: %v", err)
	}
	if string(seq) != "ACNGT" {
		t.Errorf("seq = %q, want ACNGT", seq)
	}
	if cig.String() != "5M" {
		t.Errorf("cigar = %q, want 5M", cig.String())
	}
	if !rd.needsReference {
		t.Error("an uncovered base must set needsReference")
	}
}

// TestReconstructErrors checks the malformed-feature error paths.
func TestReconstructErrors(t *testing.T) {
	rd := newReconDecoder(false)
	cases := []struct {
		name  string
		feats []readFeature
		ln    int32
	}{
		{"feature out of range", []readFeature{{code: featBases, pos: 99, bases: []byte("A")}}, 4},
		{"consuming feature behind cursor", []readFeature{
			{code: featBases, pos: 1, bases: []byte("ACGT")},
			{code: featBases, pos: 2, bases: []byte("A")},
		}, 5},
		{"bases overrun the read", []readFeature{{code: featBases, pos: 1, bases: []byte("ACGTACGT")}}, 4},
		{"negative deletion length", []readFeature{
			{code: featBases, pos: 1, bases: []byte("ACGT")},
			{code: featDeletion, pos: 5, length: -1},
		}, 4},
		{"unknown feature code", []readFeature{{code: 'Z', pos: 1}}, 4},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			if _, _, _, err := rd.reconstructMapped(c.feats, c.ln, 1); err == nil {
				t.Error("expected an error")
			}
		})
	}
	if _, _, _, err := rd.reconstructMapped(nil, -1, 1); err == nil {
		t.Error("a negative read length must error")
	}
}

// TestReconstructRefSkipAndHardClip checks the N and H features, which
// add a CIGAR operation without consuming read bases.
func TestReconstructRefSkipAndHardClip(t *testing.T) {
	rd := newReconDecoder(false)
	feats := []readFeature{
		{code: featHardClip, pos: 1, length: 3},
		{code: featBases, pos: 1, bases: []byte("ACG")},
		{code: featRefSkip, pos: 4, length: 10},
		{code: featBases, pos: 4, bases: []byte("TT")},
	}
	_, _, cig, err := rd.reconstructMapped(feats, 5, 1)
	if err != nil {
		t.Fatalf("reconstructMapped: %v", err)
	}
	if cig.String() != "3H3M10N2M" {
		t.Errorf("cigar = %q, want 3H3M10N2M", cig.String())
	}
}

// TestLinkMates checks the mate-resolution field cross-fill and the
// TLEN sign convention.
func TestLinkMates(t *testing.T) {
	up := &sam.Record{QName: "p", RName: "chr1", Pos: 100, Flag: sam.FlagPaired,
		Cigar: sam.Cigar{sam.CigarOp(10<<4 | sam.CigarMatch)}}
	down := &sam.Record{QName: "p", RName: "chr1", Pos: 200, Flag: sam.FlagPaired | sam.FlagReverse,
		Cigar: sam.Cigar{sam.CigarOp(10<<4 | sam.CigarMatch)}}
	linkMates(up, down)
	if up.RNext != "=" || down.RNext != "=" {
		t.Errorf("RNEXT should be '=' for same-reference mates, got %q/%q", up.RNext, down.RNext)
	}
	if up.PNext != 200 || down.PNext != 100 {
		t.Errorf("PNEXT cross-fill wrong: up=%d down=%d", up.PNext, down.PNext)
	}
	if up.Flag&sam.FlagMateReverse == 0 {
		t.Error("up should have the mate-reverse flag set from down")
	}
	// TLEN spans 100..209, i.e. 110, positive for the upstream record.
	if up.TLen != 110 || down.TLen != -110 {
		t.Errorf("TLEN = up:%d down:%d, want 110/-110", up.TLen, down.TLen)
	}
}

// TestResolveMatesOutOfRange checks that a next-fragment distance that
// points outside the slice is reported as an error, not a panic.
func TestResolveMatesOutOfRange(t *testing.T) {
	decoded := []*decodedRecord{
		{rec: &sam.Record{}, mateDownstream: true, nextFragment: 50},
		{rec: &sam.Record{}},
	}
	if err := resolveMates(decoded); err == nil {
		t.Error("expected an error for an out-of-range mate distance")
	}
}
