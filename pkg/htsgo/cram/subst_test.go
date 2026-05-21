package cram

import "testing"

// TestSubstMatrixIdentity checks the default / identity substitution
// matrix (the 0x1B packed byte samtools writes): code j against
// reference base r decodes to the j-th candidate base — the bases ACGTN
// with r removed, in order.
func TestSubstMatrixIdentity(t *testing.T) {
	m := newSubstMatrix([]byte{0x1B, 0x1B, 0x1B, 0x1B, 0x1B})
	// Reference base A: candidates C,G,T,N for codes 0,1,2,3.
	cases := []struct {
		ref  byte
		code byte
		want byte
	}{
		{'A', 0, 'C'}, {'A', 1, 'G'}, {'A', 2, 'T'}, {'A', 3, 'N'},
		{'C', 0, 'A'}, {'C', 1, 'G'}, {'C', 2, 'T'}, {'C', 3, 'N'},
		{'G', 0, 'A'}, {'G', 1, 'C'}, {'G', 2, 'T'}, {'G', 3, 'N'},
		{'T', 0, 'A'}, {'T', 1, 'C'}, {'T', 2, 'G'}, {'T', 3, 'N'},
		{'N', 0, 'A'}, {'N', 1, 'C'}, {'N', 2, 'G'}, {'N', 3, 'T'},
	}
	for _, c := range cases {
		if got := m.lookup(c.ref, c.code); got != c.want {
			t.Errorf("lookup(%c, %d) = %c, want %c", c.ref, c.code, got, c.want)
		}
	}
}

// TestSubstMatrixNonIdentity checks a non-identity SM byte. For ref base
// A, candidates are C,G,T,N at columns 0..3. A packed byte 0b11_10_01_00
// (0xE4) maps column j to code 3-j, so code 0 decodes to candidate 3,
// code 1 to candidate 2, and so on.
func TestSubstMatrixNonIdentity(t *testing.T) {
	m := newSubstMatrix([]byte{0xE4, 0xE4, 0xE4, 0xE4, 0xE4})
	// 0xE4 = 11 10 01 00: column 0 (candidate C) gets code 3,
	// column 1 (G) gets code 2, column 2 (T) gets code 1, column 3 (N)
	// gets code 0. Inverting: code 0 -> N, 1 -> T, 2 -> G, 3 -> C.
	want := map[byte]byte{0: 'N', 1: 'T', 2: 'G', 3: 'C'}
	for code, w := range want {
		if got := m.lookup('A', code); got != w {
			t.Errorf("lookup(A, %d) = %c, want %c", code, got, w)
		}
	}
}

// TestSubstMatrixDefault checks that a nil / short SM falls back to the
// identity mapping.
func TestSubstMatrixDefault(t *testing.T) {
	m := newSubstMatrix(nil)
	if got := m.lookup('A', 0); got != 'C' {
		t.Errorf("default matrix lookup(A,0) = %c, want C", got)
	}
	if got := m.lookup('T', 2); got != 'G' {
		t.Errorf("default matrix lookup(T,2) = %c, want G", got)
	}
}

// TestSubstMatrixAmbiguousRef checks an unrecognised reference base
// (an IUPAC ambiguity code, lower-case, or a gap) maps to the N row.
func TestSubstMatrixAmbiguousRef(t *testing.T) {
	m := newSubstMatrix([]byte{0x1B, 0x1B, 0x1B, 0x1B, 0x1B})
	// 'R' is an IUPAC ambiguity code; it should resolve through the N row.
	if got := m.lookup('R', 0); got != m.lookup('N', 0) {
		t.Errorf("an ambiguous reference base must use the N row")
	}
}

// TestRefBaseIndex checks the ACGTN index mapping.
func TestRefBaseIndex(t *testing.T) {
	for i, b := range []byte{'A', 'C', 'G', 'T', 'N'} {
		if refBaseIndex(b) != i {
			t.Errorf("refBaseIndex(%c) = %d, want %d", b, refBaseIndex(b), i)
		}
	}
	if refBaseIndex('X') != 4 {
		t.Error("an unknown base must index the N slot")
	}
}
