package cram

// substMatrix is the decode-side CRAM base-substitution matrix. A CRAM
// "X" read feature names a substituted base not directly but by a 2-bit
// substitution code relative to the reference base at that position.
// substMatrix.lookup turns a (reference base, code) pair into the read
// base.
//
// The matrix is built from the five-byte SM entry of the preservation
// map. Bases are ordered A, C, G, T, N. For reference base index r, the
// four candidate substituted bases are the bases A, C, G, T, N with
// index r removed, in that order. SM[r] packs four 2-bit fields, MSB
// first: the field at column j (j = 0..3) gives the code that decodes
// to candidate j. Inverting that mapping yields lookup[r][code].
type substMatrix struct {
	// table[refIdx][code] is the substituted read base. refIdx is the
	// index of the reference base in "ACGTN"; code is the 2-bit value the
	// BS data series carried.
	table [5][4]byte
}

// substBases is the canonical CRAM base ordering for the substitution
// matrix: A, C, G, T, N.
var substBases = [5]byte{'A', 'C', 'G', 'T', 'N'}

// substCandidates[r] lists, for reference base index r, the four
// candidate substituted bases — the entries of substBases with index r
// removed, order preserved. It is the column ordering SM[r] indexes.
var substCandidates = func() [5][4]byte {
	var out [5][4]byte
	for r := 0; r < 5; r++ {
		j := 0
		for b := 0; b < 5; b++ {
			if b == r {
				continue
			}
			out[r][j] = substBases[b]
			j++
		}
	}
	return out
}()

// newSubstMatrix builds the decode-side substitution matrix from the
// five raw bytes of the preservation map's SM entry. A nil or short sm
// yields the identity-like default htslib uses when no SM was written:
// code j decodes to candidate j directly.
func newSubstMatrix(sm []byte) substMatrix {
	var m substMatrix
	for r := 0; r < 5; r++ {
		var packed byte
		if r < len(sm) {
			packed = sm[r]
		} else {
			// Default: codes 0..3 map straight to candidates 0..3, i.e.
			// the packed byte 0b00_01_10_11.
			packed = 0x1B
		}
		// SM[r] holds four 2-bit fields MSB-first; field at column j
		// (j=0..3) is the code that decodes to candidate j. Invert it.
		for j := 0; j < 4; j++ {
			code := (packed >> uint((3-j)*2)) & 0x3
			m.table[r][code] = substCandidates[r][j]
		}
	}
	return m
}

// refBaseIndex returns the index of a reference base in "ACGTN". An
// unrecognised base (lower-case is normalised away upstream, so this
// catches IUPAC ambiguity codes and gaps) maps to the N slot — the
// safest fallback, since a substitution off an ambiguous reference base
// cannot be resolved exactly.
func refBaseIndex(b byte) int {
	switch b {
	case 'A':
		return 0
	case 'C':
		return 1
	case 'G':
		return 2
	case 'T':
		return 3
	default:
		return 4
	}
}

// lookup returns the read base a substitution of code yields against the
// given reference base.
func (m substMatrix) lookup(refBase byte, code byte) byte {
	return m.table[refBaseIndex(refBase)][code&0x3]
}

// substCodeFor returns the 2-bit substitution code that names readBase
// relative to refBase under the DEFAULT substitution matrix (SM rows all
// 0x1B), the inverse of lookup with that matrix. The writer emits no SM entry,
// so the decoder's newSubstMatrix(nil) builds exactly this default — code j
// decodes to substCandidates[refIdx][j] — and the encoder must therefore set
// code = the index of readBase among the candidates (the bases ACGTN with the
// reference base removed, order preserved). An unrecognised read base maps to
// the N candidate. Callers only invoke this for a genuine mismatch
// (readBase != refBase), so readBase is always one of the four candidates.
func substCodeFor(refBase, readBase byte) byte {
	r := refBaseIndex(refBase)
	rb := substBases[refBaseIndex(readBase)] // normalise non-ACGTN read base to N
	for j := 0; j < 4; j++ {
		if substCandidates[r][j] == rb {
			return byte(j)
		}
	}
	// readBase == refBase (not a real substitution) — should not happen; fall
	// back to code 0 so the value is at least well-defined.
	return 0
}
