package cram

import (
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// cigarM builds a single-op M CIGAR of the given length, so a record's
// EndPosition is Pos+n-1.
func cigarM(n int) sam.Cigar {
	return sam.Cigar{sam.CigarOp(uint32(n)<<4 | sam.CigarMatch)}
}

// TestLinkMatesTLen exercises the decode-side TLEN reconstruction against
// htslib's cram_decode.c mate cross-reference tie-break, the convention
// writeencode.go computeTLenOverrides also encodes to. It covers the three
// cases the fix targets: a full-overlap equal-start FR pair (READ1 gets the
// positive span), an unequal-start pair whose lower-POS mate reaches further
// right (the span magnitude must use the max end, not the higher-POS mate's
// end), and a plain non-overlapping pair (leftmost gets the positive span).
func TestLinkMatesTLen(t *testing.T) {
	tests := []struct {
		name                   string
		upPos, downPos         int64
		upLen, downLen         int
		upFlag, downFlag       uint16
		wantUpTLen, wantDnTLen int64
	}{
		{
			// (a) equal-POS full overlap FR pair: identical span [100,199]. On a
			// full overlap htslib resolves the sign by READ1, order-independently;
			// up carries READ1 so it takes the positive span.
			name:  "equal-pos full overlap, READ1 up gets +",
			upPos: 100, downPos: 100,
			upLen: 100, downLen: 100,
			upFlag:     sam.FlagPaired | sam.FlagRead1,
			downFlag:   sam.FlagPaired | sam.FlagRead2 | sam.FlagReverse,
			wantUpTLen: 100, wantDnTLen: -100,
		},
		{
			// (a') same pair but READ1 is the downstream record: the positive span
			// follows READ1 regardless of in-slice order.
			name:  "equal-pos full overlap, READ1 down gets +",
			upPos: 100, downPos: 100,
			upLen: 100, downLen: 100,
			upFlag:     sam.FlagPaired | sam.FlagRead2,
			downFlag:   sam.FlagPaired | sam.FlagRead1 | sam.FlagReverse,
			wantUpTLen: -100, wantDnTLen: 100,
		},
		{
			// (b) unequal-POS pair where the lower-POS mate (up) reaches further
			// right than the higher-POS mate. aright must be max(300,200)=300, so
			// span = 300-100+1 = 201, NOT 200-100+1 = 101 (the old bug).
			name:  "lower-pos mate reaches furthest right",
			upPos: 100, downPos: 150,
			upLen: 201, downLen: 51, // up: [100,300], down: [150,200]
			upFlag:     sam.FlagPaired | sam.FlagRead1,
			downFlag:   sam.FlagPaired | sam.FlagRead2 | sam.FlagReverse,
			wantUpTLen: 201, wantDnTLen: -201,
		},
		{
			// (c) plain non-overlapping pair: up [100,149], down [300,349].
			// span = 349-100+1 = 250; the leftmost (up) gets the positive span.
			name:  "non-overlapping, leftmost gets +",
			upPos: 100, downPos: 300,
			upLen: 50, downLen: 50,
			upFlag:     sam.FlagPaired | sam.FlagRead1,
			downFlag:   sam.FlagPaired | sam.FlagRead2 | sam.FlagReverse,
			wantUpTLen: 250, wantDnTLen: -250,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			up := &sam.Record{
				QName: "r", RName: "chr20", Pos: tt.upPos,
				Flag: tt.upFlag, Cigar: cigarM(tt.upLen),
			}
			down := &sam.Record{
				QName: "r", RName: "chr20", Pos: tt.downPos,
				Flag: tt.downFlag, Cigar: cigarM(tt.downLen),
			}
			linkMates(up, down)
			if up.TLen != tt.wantUpTLen {
				t.Errorf("up.TLen = %d, want %d", up.TLen, tt.wantUpTLen)
			}
			if down.TLen != tt.wantDnTLen {
				t.Errorf("down.TLen = %d, want %d", down.TLen, tt.wantDnTLen)
			}
			// A within-slice pair must sum to zero (equal magnitude, opposite sign).
			if up.TLen+down.TLen != 0 {
				t.Errorf("TLen signs do not cancel: up=%d down=%d", up.TLen, down.TLen)
			}
		})
	}
}

// TestLinkMatesTLenGuard confirms the cross-reference / unmapped / different-RName
// guard keeps TLEN at zero, so detached/cross-slice/mate-unmapped pairs are not
// given a spurious template length.
func TestLinkMatesTLenGuard(t *testing.T) {
	tests := []struct {
		name     string
		up, down *sam.Record
	}{
		{
			name: "mate unmapped",
			up:   &sam.Record{RName: "chr20", Pos: 100, Cigar: cigarM(50)},
			down: &sam.Record{RName: "chr20", Pos: 150, Flag: sam.FlagUnmapped},
		},
		{
			name: "different reference",
			up:   &sam.Record{RName: "chr20", Pos: 100, Cigar: cigarM(50)},
			down: &sam.Record{RName: "chr21", Pos: 150, Cigar: cigarM(50)},
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			linkMates(tt.up, tt.down)
			if tt.up.TLen != 0 || tt.down.TLen != 0 {
				t.Errorf("expected TLen 0/0, got up=%d down=%d", tt.up.TLen, tt.down.TLen)
			}
		})
	}
}
