package samtools

// Binary-free unit tests for the phase base-admission and
// variant-column decision helpers (admitPhaseBase / phaseHetLOD /
// isPhaseVariantColumn). These mirror the per-base body of phase.c's
// het-detection loop (phase.c:738-758) and pass with the upstream
// submodule UNPOPULATED — no samtools binary required.

import "testing"

// TestUnitAdmitPhaseBase exercises the per-base admission rule:
// the -Q/min-BQ quality gate, the non-ACGT drop, the
// q = clamp(min(baseQ, mapQ), 4, 63) clamp, and the packed-word
// layout q<<5 | rev<<4 | base.
func TestUnitAdmitPhaseBase(t *testing.T) {
	const minBQ = 13
	cases := []struct {
		name     string
		seq      byte
		baseQ    byte
		mapQ     byte
		rev      bool
		minBaseQ uint8
		wantOK   bool
		wantWord uint16
	}{
		{
			name: "G high quality forward", seq: 'G', baseQ: 40, mapQ: 60,
			minBaseQ: minBQ, wantOK: true,
			// q = min(40,60)=40 (<=63), base G=2: 40<<5 | 0 | 2
			wantWord: uint16(40)<<5 | 2,
		},
		{
			name: "T marginal Q15 admitted at default -Q", seq: 'T', baseQ: 15, mapQ: 60,
			minBaseQ: minBQ, wantOK: true,
			// q = min(15,60)=15, base T=3
			wantWord: uint16(15)<<5 | 3,
		},
		{
			name: "marginal Q12 dropped below default -Q 13", seq: 'T', baseQ: 12, mapQ: 60,
			minBaseQ: minBQ, wantOK: false,
		},
		{
			name: "marginal Q12 admitted when -Q lowered to 1", seq: 'T', baseQ: 12, mapQ: 60,
			minBaseQ: 1, wantOK: true,
			// q = min(12,60)=12, base T=3
			wantWord: uint16(12)<<5 | 3,
		},
		{
			name: "N base dropped (non-ACGT)", seq: 'N', baseQ: 40, mapQ: 60,
			minBaseQ: minBQ, wantOK: false,
		},
		{
			name: "low MAPQ caps q below baseQ", seq: 'A', baseQ: 40, mapQ: 7,
			minBaseQ: minBQ, wantOK: true,
			// q = min(40,7)=7, base A=0
			wantWord: uint16(7) << 5,
		},
		{
			name: "q clamped up to floor 4", seq: 'C', baseQ: 40, mapQ: 2,
			minBaseQ: 1, wantOK: true,
			// q = min(40,2)=2 -> clamp to 4, base C=1
			wantWord: uint16(4)<<5 | 1,
		},
		{
			name: "q clamped down to ceiling 63", seq: 'G', baseQ: 90, mapQ: 90,
			minBaseQ: minBQ, wantOK: true,
			// q = min(90,90)=90 -> clamp to 63, base G=2
			wantWord: uint16(63)<<5 | 2,
		},
		{
			name: "reverse-strand bit set", seq: 'A', baseQ: 30, mapQ: 60, rev: true,
			minBaseQ: minBQ, wantOK: true,
			// q=30, rev=1, base A=0: 30<<5 | 1<<4 | 0
			wantWord: uint16(30)<<5 | 1<<4,
		},
		{
			name: "baseQ exactly at threshold is admitted", seq: 'G', baseQ: 13, mapQ: 60,
			minBaseQ: minBQ, wantOK: true,
			wantWord: uint16(13)<<5 | 2,
		},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			word, ok := admitPhaseBase(tc.seq, tc.baseQ, tc.mapQ, tc.rev, tc.minBaseQ)
			if ok != tc.wantOK {
				t.Fatalf("admitPhaseBase ok=%v, want %v", ok, tc.wantOK)
			}
			if ok && word != tc.wantWord {
				t.Fatalf("admitPhaseBase word=%#x, want %#x", word, tc.wantWord)
			}
		})
	}
}

// TestUnitPhaseHetLOD checks the LOD extraction from a gl2cns word
// matches phase.c's (c&0xffff)>>2 expression.
func TestUnitPhaseHetLOD(t *testing.T) {
	// c = 0x00060057 is the observed gl2cns word for the marginal
	// 3G/3T Q15 column (het bit 1<<18, alleles G/T, LOD field).
	// (0x0057 & 0xffff) >> 2 = 0x57 >> 2 = 0x15 = 21.
	if got := phaseHetLOD(0x00060057); got != 21 {
		t.Fatalf("phaseHetLOD(0x00060057)=%d, want 21", got)
	}
	// Homozygous gl2cns returns 0 -> LOD 0.
	if got := phaseHetLOD(0); got != 0 {
		t.Fatalf("phaseHetLOD(0)=%d, want 0", got)
	}
}

// TestUnitIsPhaseVariantColumn pins the variant-column gate: a column
// is a variant iff it is forced by the -l site list (inSet) OR its LOD
// reaches the -q minVarLOD threshold. This is the exact decision that
// the misrouted-flag bug used to get wrong: with the default 37 the
// marginal LOD-21 column was dropped; once -q is honoured the same
// column passes at any threshold <= 21.
func TestUnitIsPhaseVariantColumn(t *testing.T) {
	const marginal = uint32(0x00060057) // LOD 21
	cases := []struct {
		name      string
		c         uint32
		minVarLOD int
		inSet     bool
		want      bool
	}{
		{"marginal dropped at default 37", marginal, 37, false, false},
		{"marginal admitted at q=21 (boundary)", marginal, 21, false, true},
		{"marginal admitted at q=1", marginal, 1, false, true},
		{"marginal admitted at q=0", marginal, 0, false, true},
		{"marginal dropped at q=22", marginal, 22, false, false},
		{"forced by site list passes despite high threshold", marginal, 37, true, true},
		{"homozygous never a variant unless forced", 0, 0, false, true},
		{"homozygous dropped at q=1", 0, 1, false, false},
		{"homozygous forced by site list", 0, 37, true, true},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			if got := isPhaseVariantColumn(tc.c, tc.minVarLOD, tc.inSet); got != tc.want {
				t.Fatalf("isPhaseVariantColumn(%#x, %d, %v)=%v, want %v",
					tc.c, tc.minVarLOD, tc.inSet, got, tc.want)
			}
		})
	}
}
