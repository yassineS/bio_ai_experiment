package strfinder

import (
	"reflect"
	"testing"
)

// encodeSeq turns an uppercase ASCII string into the 2-bit byte slice
// the str_finder algorithm operates on. Lower-case letters keep their
// raw ASCII byte so the lowerOnly filter can spot them, but their low
// two bits are NOT the nucleotide code — tests that care about
// lowerOnly use a separate helper.
func encodeSeq(s string) []byte {
	out := make([]byte, len(s))
	for i, b := range []byte(s) {
		out[i] = EncodeNt(b)
	}
	return out
}

// TestFindSTR_AA: a 2-base homopolymer "AA" yields one rlen=1 hit
// covering both bases. Mirrors str_finder.c:231 — the j>=1 + rlen=1
// branch fires when the new 2-bit base equals the previous one.
func TestFindSTR_AA(t *testing.T) {
	got := FindSTR(encodeSeq("AA"), false)
	want := []RepElement{{Start: 0, End: 1, RepLen: 1}}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("FindSTR(AA) = %+v, want %+v", got, want)
	}
}

// TestFindSTR_AAA: at i=1 the rlen=1 branch fires; its scan-ahead loop
// already extends end past i=2 (str_finder.c:69 `while (cp2 < cp_end)`).
// At i=2 the overlap-dedup guard (str_finder.c:48-50) suppresses a
// second add_rep — net result is a single entry spanning the whole run.
func TestFindSTR_AAA(t *testing.T) {
	got := FindSTR(encodeSeq("AAA"), false)
	want := []RepElement{{Start: 0, End: 2, RepLen: 1}}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("FindSTR(AAA) = %+v, want %+v", got, want)
	}
}

// TestFindSTR_AAAA: an even longer homopolymer. The rlen=1 entry at
// i=1 already covers the whole run because addRep's scan-ahead loop
// (str_finder.c:69) walks cp2 to the end of the run. Subsequent
// addRep calls at i=2, i=3 (both rlen=1 and rlen=2) are all suppressed
// by the "already handled in the previous overlap" guard
// (str_finder.c:48-50) because the existing tail's interval [0..3]
// covers `pos - rlen*2 + 1 .. pos` for each later call. Net result:
// exactly one entry {0,3,1}.
func TestFindSTR_AAAA(t *testing.T) {
	got := FindSTR(encodeSeq("AAAA"), false)
	want := []RepElement{{Start: 0, End: 3, RepLen: 1}}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("FindSTR(AAAA) = %+v, want %+v", got, want)
	}
}

// TestFindSTR_AAAAA hand-traces add_rep for a 5-base homopolymer. The
// first addRep at i=1 / rlen=1 walks cp2 all the way to cp_end=5, so
// end = pos + cp2 - (pos+1) = 4 and start walks back to 0. Every
// subsequent addRep call (rlen=1 at i=2,3,4 and rlen=2 at i=3,4) is
// suppressed by the overlap guard. Net: exactly {0,4,1}.
func TestFindSTR_AAAAA(t *testing.T) {
	got := FindSTR(encodeSeq("AAAAA"), false)
	want := []RepElement{{Start: 0, End: 4, RepLen: 1}}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("FindSTR(AAAAA) = %+v, want %+v", got, want)
	}
}

// TestFindSTR_AC_AC: a dinucleotide repeat "ACAC". Only the rlen=2
// branch at i=3 fires (rlen=1 never matches because adjacent bases
// differ). addRep walks cp1 back two non-pad bases to cp1=1, then to
// cp1=2 after the inner skip-loop runs; cp2=4=cpEnd so the scan-ahead
// loop does not extend it. end = pos+cp2-(pos+1) = 3 and the walk-back
// pulls startPos to 0. Net: exactly one entry {0,3,2}.
func TestFindSTR_ACAC(t *testing.T) {
	got := FindSTR(encodeSeq("ACAC"), false)
	want := []RepElement{{Start: 0, End: 3, RepLen: 2}}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("FindSTR(ACAC) = %+v, want %+v", got, want)
	}
}

// TestFindSTR_AAAACG hand-traces the boundary case where a homopolymer
// run ends before the end of the sequence. addRep at i=1 / rlen=1
// scan-ahead stops when cons[3]=A != cons[4]=C, so end = 3 (not the
// full length). Later addRep calls at i=2,3 (rlen=1) and i=3 (rlen=2)
// are all suppressed by the overlap guard. Net: exactly {0,3,1}.
func TestFindSTR_AAAACG(t *testing.T) {
	got := FindSTR(encodeSeq("AAAACG"), false)
	want := []RepElement{{Start: 0, End: 3, RepLen: 1}}
	if !reflect.DeepEqual(got, want) {
		t.Fatalf("FindSTR(AAAACG) = %+v, want %+v", got, want)
	}
}

// TestFindSTR_Mixed: AGGGGAGG mixes a run of G's and a smaller break.
// The example comment in str_finder.c (lines 312-316) describes this
// shape. We assert that FindSTR returns at least one homopolymer entry
// over the GGGG span.
func TestFindSTR_Mixed(t *testing.T) {
	got := FindSTR(encodeSeq("AGGGGAGG"), false)
	if len(got) == 0 {
		t.Fatal("FindSTR(AGGGGAGG) returned no repeats")
	}
	// Look for a hit that covers the GGGG run starting at index 1.
	found := false
	for _, r := range got {
		if r.Start <= 1 && r.End >= 4 {
			found = true
			break
		}
	}
	if !found {
		t.Fatalf("FindSTR(AGGGGAGG) = %+v, missing GGGG span", got)
	}
}

// TestFindSTR_NoRepeat: a sequence with no adjacent repeats yields an
// empty slice.
func TestFindSTR_NoRepeat(t *testing.T) {
	got := FindSTR(encodeSeq("ACGTA"), false)
	if len(got) != 0 {
		t.Fatalf("FindSTR(ACGTA) = %+v, want empty", got)
	}
}

// TestFindSTR_LowerOnly: when lowerOnly is true the repeat is only
// reported if its span contains a lower-case letter. We feed in a
// pre-encoded byte slice with the lower-case markers preserved by
// stashing them in the same slot as the 2-bit value (b|0x20). The
// addRep lowerOnly check inspects isLowerASCII on the original byte;
// so to exercise it we pass a slice where one base in the run is a raw
// lower-case byte ('a') instead of the 2-bit code. The repeat-detection
// logic uses cons[i]&3, so 'a' (0x61 → low bits 01) is treated as 'C',
// breaking the homopolymer. To keep the run intact we sprinkle a 'g'
// (0x67 → low bits 11, i.e. 'T') — also breaking it. The cleanest test
// is therefore to verify the *negative* case: an all-upper-case repeat
// returns nothing under lowerOnly=true.
func TestFindSTR_LowerOnly_UpperReturnsNothing(t *testing.T) {
	got := FindSTR(encodeSeq("AAAA"), true)
	if len(got) != 0 {
		t.Fatalf("FindSTR(AAAA, lowerOnly=true) = %+v, want empty", got)
	}
}

// TestFindSTR_LowerOnly_Lower verifies the positive lower-case path: a
// custom byte slice carrying the 2-bit code 0 (== 'A') plus the
// lower-case marker bit. We use the raw byte 0x20 (which has 2-bit code
// 0 and is below 'a', so isLowerASCII returns false) for the non-marker
// bases, and the byte 'a' (0x61, low two bits 0b01 = 'C') would
// actually break the run. To inject a "lower" marker without changing
// the 2-bit code we use the byte 0x60 (low two bits 0 = 'A',
// isLowerASCII still false — outside [a..z]). isLowerASCII only matches
// 'a'..'z' so a true positive needs a low-bits-0 letter, of which there
// is exactly one: the letter 'p' (0x70, low bits 0). We test that.
func TestFindSTR_LowerOnly_Lower(t *testing.T) {
	// Build "ApAA" where 'p' is the marker (2-bit code = 0x70 & 3 = 0).
	cons := []byte{0, 'p', 0, 0}
	got := FindSTR(cons, true)
	if len(got) == 0 {
		t.Fatalf("FindSTR with lower-case marker = empty, want at least one hit")
	}
}

// TestFindSTR_PadSkip: the '*' byte is treated as padding and skipped
// while building the rolling word. "A*A" should look like "AA" to
// FindSTR.
func TestFindSTR_PadSkip(t *testing.T) {
	cons := []byte{0, '*', 0}
	got := FindSTR(cons, false)
	if len(got) == 0 {
		t.Fatalf("FindSTR(A*A) = empty, want at least one hit")
	}
}

// TestFindSTR64_LongRepeat: a 28-base tandem repeat (14-mer x 2) lands
// in the rlen=14 branch of find_STR64. The test asserts FindSTR64
// returns a hit covering it.
func TestFindSTR64_LongRepeat(t *testing.T) {
	// 14-mer "ACGTACGTACGTAC" repeated.
	unit := "ACGTACGTACGTAC"
	got := FindSTR64(encodeSeq(unit+unit), false)
	if len(got) == 0 {
		t.Fatalf("FindSTR64(28-mer repeat) = empty")
	}
}

// TestConsMarkSTR_Homopolymer: ConsMarkSTR over a clean homopolymer
// marks every base of the repeat span with a non-zero bit.
func TestConsMarkSTR_Homopolymer(t *testing.T) {
	mask := ConsMarkSTR(encodeSeq("AAAAA"), false)
	for i, b := range mask {
		if b == 0 {
			t.Fatalf("ConsMarkSTR(AAAAA) mask[%d] = 0, want non-zero", i)
		}
	}
}

// TestConsMarkSTR_NoRepeat: ConsMarkSTR returns all zeroes on a
// non-repeating sequence.
func TestConsMarkSTR_NoRepeat(t *testing.T) {
	mask := ConsMarkSTR(encodeSeq("ACGTA"), false)
	for i, b := range mask {
		if b != 0 {
			t.Fatalf("ConsMarkSTR(ACGTA) mask[%d] = %d, want 0", i, b)
		}
	}
}

// TestRepElementSortOrder: FindSTR returns hits in the order they were
// appended by add_rep (DL_APPEND tail), so their End coordinates are
// monotonically non-decreasing.
func TestRepElementSortOrder(t *testing.T) {
	got := FindSTR(encodeSeq("AAAAACCCCC"), false)
	for i := 1; i < len(got); i++ {
		if got[i].End < got[i-1].End {
			t.Fatalf("FindSTR returned out-of-order End: %+v", got)
		}
	}
}
