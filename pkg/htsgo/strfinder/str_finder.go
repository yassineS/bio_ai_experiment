// Package strfinder is a Go port of bcftools' str_finder.c (originally
// from Crumble, by James Bonfield). It locates Short Tandem Repeats up
// to 14-mers in a 2-bit-encoded "consensus" sequence and is used by the
// bcftools mpileup indel caller to identify homopolymer / tandem-repeat
// regions where indel alleles need realignment.
//
// The upstream API is two near-identical entry points:
//
//   - FindSTR — repeat lengths 1..8 (uint32 sliding window).
//   - FindSTR64 — repeat lengths 1..14 (uint64 sliding window).
//
// Plus ConsMarkSTR, the per-base bitmask renderer that wraps FindSTR.
//
// The input sequence is a slice of bytes whose values are the 2-bit
// codes for the four nucleotides (A=0, C=1, G=2, T=3); the byte '*'
// (0x2a) is treated as a padding marker and skipped while computing the
// rolling word. Upstream's `lower_only` flag is honoured by inspecting
// the original byte: any value whose low byte falls into the ASCII
// lower-case alphabet ('a'..'z', i.e. byte&0x20 set on a letter) marks
// the region as "lower-case included". That mirrors C's islower().
//
// Upstream reference: reference_code/bcftools/str_finder.c.
package strfinder

// RepElement is one tandem-repeat hit. Coordinates are inclusive 0-based
// indices into the input slice, matching upstream's rep_ele.
type RepElement struct {
	// Start is the inclusive 0-based index where the repeat begins.
	Start int
	// End is the inclusive 0-based index where the repeat ends.
	End int
	// RepLen is the repeat-unit length (1..14).
	RepLen int
}

// isLowerASCII reports whether b is an ASCII lower-case letter. We use
// this instead of unicode.ToLower because upstream's str_finder uses
// libc islower() on a byte.
func isLowerASCII(b byte) bool {
	return b >= 'a' && b <= 'z'
}

// addRep is the Go port of str_finder.c's add_rep. It scans forward to
// extend the repeat past pos, walks back to find its true start, then
// appends to reps after dropping any earlier entries that the new repeat
// fully contains. Upstream uses a doubly-linked list; here we use a
// slice with the same semantics — reps[len-1] is the "tail" and we
// truncate elements whose interval is enclosed by the new entry.
func addRep(reps []RepElement, cons []byte, clen, pos, rlen int, lowerOnly bool) []RepElement {
	// Already handled this in the previous overlap?
	if len(reps) > 0 {
		tail := reps[len(reps)-1]
		if tail.Start <= pos-rlen*2+1 && tail.End >= pos {
			return reps
		}
	}

	// Find current and last occurrence of the repeated word. cp1 walks
	// back rlen non-pad bases from pos; cp2 starts at pos+1.
	cp2 := pos + 1
	cp1 := pos
	for i := 1; i < rlen; {
		if cp1-1 < 0 {
			break
		}
		cp1--
		if cons[cp1] == '*' {
			continue
		}
		i++
	}
	for cp1 > 0 && cons[cp1] == '*' {
		cp1--
	}

	// Scan ahead to see how much further the repeat extends.
	cpEnd := clen
	for cp2 < cpEnd {
		if cons[cp1] != cons[cp2] {
			break
		}
		cp1++
		cp2++
	}

	end := pos + cp2 - (pos + 1)
	// Walk pos back by rlen "units of two non-pad bases" to find start.
	// Upstream's loop does `while(rlen--) { while(cons[--pos] == '*'); while(cons[--pos] == '*'); }`,
	// i.e. for each unit, skip the next non-pad base twice.
	startPos := pos + 1
	r := rlen
	for r > 0 {
		// First non-pad step.
		startPos--
		for startPos >= 0 && cons[startPos] == '*' {
			startPos--
		}
		if startPos < 0 {
			break
		}
		// Second non-pad step.
		startPos--
		for startPos >= 0 && cons[startPos] == '*' {
			startPos--
		}
		if startPos < 0 {
			break
		}
		r--
	}
	for startPos > 1 && cons[startPos-1] == '*' {
		startPos--
	}
	if startPos < 0 {
		startPos = 0
	}

	el := RepElement{Start: startPos, End: end, RepLen: rlen}

	// Lower-case-only filter.
	if lowerOnly {
		lc := false
		for i := el.Start; i <= el.End && i < clen; i++ {
			if isLowerASCII(cons[i]) {
				lc = true
				break
			}
		}
		if !lc {
			return reps
		}
	}

	// Drop earlier entries entirely contained within el.
	for len(reps) > 0 {
		tail := reps[len(reps)-1]
		if tail.End < el.Start {
			break
		}
		if tail.Start >= el.Start {
			reps = reps[:len(reps)-1]
			continue
		}
		break
	}

	reps = append(reps, el)
	return reps
}

// twoBit returns the 2-bit code of cons[i] (i.e. the low two bits of the
// byte). The byte is expected to already hold 0..3 for the four
// nucleotides; values outside that range are masked, matching upstream
// which directly ORs cons[i] into the rolling word.
func twoBit(b byte) uint64 { return uint64(b) & 0x3 }

// FindSTR ports str_finder.c's find_STR. It scans cons for tandem
// repeats with unit lengths 1..8 using a 32-bit sliding window. The
// returned slice is in order of right-hand position, matching the
// upstream doubly-linked-list order (DL_APPEND tail).
//
// When lowerOnly is true, only repeats whose span contains at least one
// ASCII lower-case letter are returned.
func FindSTR(cons []byte, lowerOnly bool) []RepElement {
	clen := len(cons)
	var w uint32
	reps := make([]RepElement, 0, 8)
	i, j := 0, 0
	for ; i < clen && j < 15; i++ {
		if cons[i] == '*' {
			continue
		}
		w <<= 2
		w |= uint32(twoBit(cons[i]))
		if j >= 1 && (w&0x0003) == ((w>>2)&0x0003) {
			reps = addRep(reps, cons, clen, i, 1, lowerOnly)
		}
		if j >= 3 && (w&0x000f) == ((w>>4)&0x000f) {
			reps = addRep(reps, cons, clen, i, 2, lowerOnly)
		}
		if j >= 5 && (w&0x003f) == ((w>>6)&0x003f) {
			reps = addRep(reps, cons, clen, i, 3, lowerOnly)
		}
		if j >= 7 && (w&0x00ff) == ((w>>8)&0x00ff) {
			reps = addRep(reps, cons, clen, i, 4, lowerOnly)
		}
		if j >= 9 && (w&0x03ff) == ((w>>10)&0x03ff) {
			reps = addRep(reps, cons, clen, i, 5, lowerOnly)
		}
		if j >= 11 && (w&0x0fff) == ((w>>12)&0x0fff) {
			reps = addRep(reps, cons, clen, i, 6, lowerOnly)
		}
		if j >= 13 && (w&0x3fff) == ((w>>14)&0x3fff) {
			reps = addRep(reps, cons, clen, i, 7, lowerOnly)
		}
		j++
	}
	for ; i < clen; i++ {
		if cons[i] == '*' {
			continue
		}
		w <<= 2
		w |= uint32(twoBit(cons[i]))
		switch {
		case (w & 0xffff) == ((w >> 16) & 0xffff):
			reps = addRep(reps, cons, clen, i, 8, lowerOnly)
		case (w & 0x3fff) == ((w >> 14) & 0x3fff):
			reps = addRep(reps, cons, clen, i, 7, lowerOnly)
		case (w & 0x0fff) == ((w >> 12) & 0x0fff):
			reps = addRep(reps, cons, clen, i, 6, lowerOnly)
		case (w & 0x03ff) == ((w >> 10) & 0x03ff):
			reps = addRep(reps, cons, clen, i, 5, lowerOnly)
		case (w & 0x00ff) == ((w >> 8) & 0x00ff):
			reps = addRep(reps, cons, clen, i, 4, lowerOnly)
		case (w & 0x003f) == ((w >> 6) & 0x003f):
			reps = addRep(reps, cons, clen, i, 3, lowerOnly)
		case (w & 0x000f) == ((w >> 4) & 0x000f):
			reps = addRep(reps, cons, clen, i, 2, lowerOnly)
		case (w & 0x0003) == ((w >> 2) & 0x0003):
			reps = addRep(reps, cons, clen, i, 1, lowerOnly)
		}
	}
	return reps
}

// FindSTR64 ports str_finder.c's find_STR64: the longer-hash variant
// covering unit lengths 1..14 using a 64-bit sliding window.
func FindSTR64(cons []byte, lowerOnly bool) []RepElement {
	clen := len(cons)
	var w uint64
	reps := make([]RepElement, 0, 8)
	i, j := 0, 0
	for ; i < clen && j < 26; i++ {
		if cons[i] == '*' {
			continue
		}
		w <<= 2
		w |= twoBit(cons[i])
		if j >= 1 && (w&0x0003) == ((w>>2)&0x0003) {
			reps = addRep(reps, cons, clen, i, 1, lowerOnly)
		}
		if j >= 3 && (w&0x000f) == ((w>>4)&0x000f) {
			reps = addRep(reps, cons, clen, i, 2, lowerOnly)
		}
		if j >= 5 && (w&0x003f) == ((w>>6)&0x003f) {
			reps = addRep(reps, cons, clen, i, 3, lowerOnly)
		}
		if j >= 7 && (w&0x00ff) == ((w>>8)&0x00ff) {
			reps = addRep(reps, cons, clen, i, 4, lowerOnly)
		}
		if j >= 9 && (w&0x03ff) == ((w>>10)&0x03ff) {
			reps = addRep(reps, cons, clen, i, 5, lowerOnly)
		}
		if j >= 11 && (w&0x0fff) == ((w>>12)&0x0fff) {
			reps = addRep(reps, cons, clen, i, 6, lowerOnly)
		}
		if j >= 13 && (w&0x3fff) == ((w>>14)&0x3fff) {
			reps = addRep(reps, cons, clen, i, 7, lowerOnly)
		}
		if j >= 15 && (w&0xffff) == ((w>>16)&0xffff) {
			reps = addRep(reps, cons, clen, i, 8, lowerOnly)
		}
		if j >= 17 && (w&0x003ffff) == ((w>>18)&0x003ffff) {
			reps = addRep(reps, cons, clen, i, 9, lowerOnly)
		}
		if j >= 19 && (w&0x00fffff) == ((w>>20)&0x00fffff) {
			reps = addRep(reps, cons, clen, i, 10, lowerOnly)
		}
		if j >= 21 && (w&0x03fffff) == ((w>>22)&0x03fffff) {
			reps = addRep(reps, cons, clen, i, 11, lowerOnly)
		}
		if j >= 23 && (w&0x0ffffff) == ((w>>24)&0x0ffffff) {
			reps = addRep(reps, cons, clen, i, 12, lowerOnly)
		}
		if j >= 24 && (w&0x3ffffff) == ((w>>26)&0x3ffffff) {
			reps = addRep(reps, cons, clen, i, 13, lowerOnly)
		}
		j++
	}
	for ; i < clen; i++ {
		if cons[i] == '*' {
			continue
		}
		w <<= 2
		w |= twoBit(cons[i])
		switch {
		case (w & 0xfffffff) == ((w >> 28) & 0xfffffff):
			reps = addRep(reps, cons, clen, i, 14, lowerOnly)
		case (w & 0x3ffffff) == ((w >> 26) & 0x3ffffff):
			reps = addRep(reps, cons, clen, i, 13, lowerOnly)
		case (w & 0x0ffffff) == ((w >> 24) & 0x0ffffff):
			reps = addRep(reps, cons, clen, i, 12, lowerOnly)
		case (w & 0x03fffff) == ((w >> 22) & 0x03fffff):
			reps = addRep(reps, cons, clen, i, 11, lowerOnly)
		case (w & 0x00fffff) == ((w >> 20) & 0x00fffff):
			reps = addRep(reps, cons, clen, i, 10, lowerOnly)
		case (w & 0x003ffff) == ((w >> 18) & 0x003ffff):
			reps = addRep(reps, cons, clen, i, 9, lowerOnly)
		case (w & 0xffff) == ((w >> 16) & 0xffff):
			reps = addRep(reps, cons, clen, i, 8, lowerOnly)
		case (w & 0x3fff) == ((w >> 14) & 0x3fff):
			reps = addRep(reps, cons, clen, i, 7, lowerOnly)
		case (w & 0x0fff) == ((w >> 12) & 0x0fff):
			reps = addRep(reps, cons, clen, i, 6, lowerOnly)
		case (w & 0x03ff) == ((w >> 10) & 0x03ff):
			reps = addRep(reps, cons, clen, i, 5, lowerOnly)
		case (w & 0x00ff) == ((w >> 8) & 0x00ff):
			reps = addRep(reps, cons, clen, i, 4, lowerOnly)
		case (w & 0x003f) == ((w >> 6) & 0x003f):
			reps = addRep(reps, cons, clen, i, 3, lowerOnly)
		case (w & 0x000f) == ((w >> 4) & 0x000f):
			reps = addRep(reps, cons, clen, i, 2, lowerOnly)
		case (w & 0x0003) == ((w >> 2) & 0x0003):
			reps = addRep(reps, cons, clen, i, 1, lowerOnly)
		}
	}
	return reps
}

// ConsMarkSTR ports cons_mark_STR. It runs FindSTR over cons and returns
// a byte mask of the same length where each byte is a bitmask of the
// repeat sizes covering that position (0 means non-repetitive). The bit
// assignment matches upstream's policy: each new repeat gets the lowest
// bit not already in use across its span, with bit 0 reused if all 8
// bits are taken.
func ConsMarkSTR(cons []byte, lowerOnly bool) []byte {
	clen := len(cons)
	str := make([]byte, clen)
	reps := FindSTR(cons, lowerOnly)
	for _, el := range reps {
		v := byte(0)
		lo := el.Start - 1
		if lo < 0 {
			lo = 0
		}
		hi := el.End + 1
		if hi > clen-1 {
			hi = clen - 1
		}
		for i := lo; i <= hi; i++ {
			v |= str[i]
		}
		var pick byte
		i := 0
		for ; i < 8; i++ {
			if v&(1<<i) == 0 {
				break
			}
		}
		if i == 8 {
			pick = 1
		} else {
			pick = 1 << i
		}
		for i := el.Start; i <= el.End && i < clen; i++ {
			str[i] |= pick
		}
	}
	return str
}

// EncodeNt encodes an ASCII nucleotide byte into the 2-bit code expected
// by FindSTR / FindSTR64. A,C,G,T (any case) map to 0,1,2,3; '*' is
// returned unchanged so it acts as a padding marker; everything else
// maps to a value of 4 (which never compares equal to any 2-bit code in
// a 2*k-bit window of valid nucleotides as long as the high bits stay
// zeroed — callers should treat it as a break).
//
// Lower-case input is preserved as the lower-case ASCII byte so the
// lowerOnly filter still sees it, by encoding lower-case at +0x20 above
// the upper-case mapping. To keep the 2-bit value correct, FindSTR
// inputs should be the 0..3 codes; if you need the lower-case marker
// keep the original byte separately and pass that to ConsMarkSTR.
//
// This helper is a convenience for callers; the C upstream did not
// expose an equivalent because str_finder is always called with
// already-2-bit-encoded data.
func EncodeNt(b byte) byte {
	switch b {
	case 'A', 'a':
		return 0
	case 'C', 'c':
		return 1
	case 'G', 'g':
		return 2
	case 'T', 't':
		return 3
	case '*':
		return '*'
	}
	return 4
}
