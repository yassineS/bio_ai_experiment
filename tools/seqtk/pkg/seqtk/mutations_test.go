package seqtk

import (
	"bytes"
	"math/rand"
	"strings"
	"testing"
)

func TestParseMutfile(t *testing.T) {
	tests := []struct {
		name    string
		input   string
		want    MutationSet
		wantLen int
	}{
		{
			name:  "three columns",
			input: "chr1\t10\tA\nchr1\t20\tC\nchr2\t5\tG\n",
			want: MutationSet{
				"chr1": {{Pos: 9, Base: 'A'}, {Pos: 19, Base: 'C'}},
				"chr2": {{Pos: 4, Base: 'G'}},
			},
			wantLen: 2,
		},
		{
			name:  "four columns (upstream format chrom pos ref alt)",
			input: "chr1\t10\tT\tA\nchr1\t20\tG\tC\n",
			want: MutationSet{
				"chr1": {{Pos: 9, Base: 'A'}, {Pos: 19, Base: 'C'}},
			},
			wantLen: 1,
		},
		{
			name:    "comments and blanks are skipped",
			input:   "# comment\n\nchr1\t1\tT\n",
			want:    MutationSet{"chr1": {{Pos: 0, Base: 'T'}}},
			wantLen: 1,
		},
		{
			name:    "invalid base is skipped",
			input:   "chr1\t1\tQQ\nchr1\t2\tA\n",
			want:    MutationSet{"chr1": {{Pos: 1, Base: 'A'}}},
			wantLen: 1,
		},
		{
			name:    "non-integer position is skipped",
			input:   "chr1\tfoo\tA\nchr1\t3\tC\n",
			want:    MutationSet{"chr1": {{Pos: 2, Base: 'C'}}},
			wantLen: 1,
		},
		{
			name:    "zero position is skipped",
			input:   "chr1\t0\tA\nchr1\t1\tG\n",
			want:    MutationSet{"chr1": {{Pos: 0, Base: 'G'}}},
			wantLen: 1,
		},
		{
			name:    "too few fields is skipped",
			input:   "chr1\t1\nchr1\t1\tA\n",
			want:    MutationSet{"chr1": {{Pos: 0, Base: 'A'}}},
			wantLen: 1,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := ParseMutfile(strings.NewReader(tt.input))
			if err != nil {
				t.Fatalf("ParseMutfile error: %v", err)
			}
			if len(got) != tt.wantLen {
				t.Errorf("got %d chromosomes, want %d", len(got), tt.wantLen)
			}
			for chrom, muts := range tt.want {
				gotMuts, ok := got[chrom]
				if !ok {
					t.Errorf("missing chrom %q", chrom)
					continue
				}
				if len(gotMuts) != len(muts) {
					t.Errorf("chrom %q: got %d muts, want %d", chrom, len(gotMuts), len(muts))
					continue
				}
				for i, m := range muts {
					if gotMuts[i] != m {
						t.Errorf("chrom %q mut[%d]: got %+v want %+v", chrom, i, gotMuts[i], m)
					}
				}
			}
		})
	}
}

func TestMutfa_BasicSubstitution(t *testing.T) {
	in := ">chr1\nACGTACGT\n>chr2\nTTTTTTTT\n"
	mut := "chr1\t1\tT\nchr1\t4\tA\nchr2\t5\tG\n"
	want := ">chr1\nTCGAACGT\n>chr2\nTTTTGTTT\n"

	var buf bytes.Buffer
	if err := Mutfa(strings.NewReader(in), strings.NewReader(mut), &buf); err != nil {
		t.Fatalf("Mutfa: %v", err)
	}
	if got := buf.String(); got != want {
		t.Errorf("Mutfa output:\ngot  %q\nwant %q", got, want)
	}
}

func TestMutfa_PreservesLineWidth(t *testing.T) {
	// Sequence split across multiple lines of width 4: positions 1..12.
	in := ">chr1\nACGT\nACGT\nACGT\n"
	// Mutate position 1 (-> X), 5 (-> Y), 12 (-> Z). Line layout must be kept.
	mut := "chr1\t1\tX\nchr1\t5\tY\nchr1\t12\tZ\n"
	want := ">chr1\nXCGT\nYCGT\nACGZ\n"

	var buf bytes.Buffer
	if err := Mutfa(strings.NewReader(in), strings.NewReader(mut), &buf); err != nil {
		t.Fatalf("Mutfa: %v", err)
	}
	if got := buf.String(); got != want {
		t.Errorf("Mutfa output:\ngot  %q\nwant %q", got, want)
	}
}

func TestMutfa_MissingChromIgnored(t *testing.T) {
	in := ">chr1\nACGT\n"
	mut := "chr2\t1\tA\nchr1\t2\tT\n"
	want := ">chr1\nATGT\n"

	var buf bytes.Buffer
	if err := Mutfa(strings.NewReader(in), strings.NewReader(mut), &buf); err != nil {
		t.Fatalf("Mutfa: %v", err)
	}
	if got := buf.String(); got != want {
		t.Errorf("Mutfa output:\ngot  %q\nwant %q", got, want)
	}
}

func TestMutfa_OutOfRangePositionSkipped(t *testing.T) {
	in := ">chr1\nACGT\n"
	mut := "chr1\t999\tA\nchr1\t1\tG\n" // 999 is out of range; 1 should apply
	want := ">chr1\nGCGT\n"

	var buf bytes.Buffer
	if err := Mutfa(strings.NewReader(in), strings.NewReader(mut), &buf); err != nil {
		t.Fatalf("Mutfa: %v", err)
	}
	if got := buf.String(); got != want {
		t.Errorf("Mutfa output:\ngot  %q\nwant %q", got, want)
	}
}

func TestMutfa_HeaderWithDescription(t *testing.T) {
	// Mutfa keys mutations by the first whitespace-delimited token after '>'.
	in := ">chr1 some description\nACGT\n"
	mut := "chr1\t2\tT\n"
	want := ">chr1 some description\nATGT\n"

	var buf bytes.Buffer
	if err := Mutfa(strings.NewReader(in), strings.NewReader(mut), &buf); err != nil {
		t.Fatalf("Mutfa: %v", err)
	}
	if got := buf.String(); got != want {
		t.Errorf("Mutfa output:\ngot  %q\nwant %q", got, want)
	}
}

func TestMutfa_MultiRecordOutOfRangeWarning(t *testing.T) {
	// Two records, the second has a mutation past its end. The mutation file
	// also lists a third sequence that is not in the input.
	in := ">chr1\nAAAA\n>chr2\nCCCC\n"
	mut := "chr1\t2\tT\nchr2\t99\tG\nghost\t1\tA\n"
	want := ">chr1\nATAA\n>chr2\nCCCC\n"

	var buf bytes.Buffer
	if err := Mutfa(strings.NewReader(in), strings.NewReader(mut), &buf); err != nil {
		t.Fatalf("Mutfa: %v", err)
	}
	if got := buf.String(); got != want {
		t.Errorf("Mutfa output:\ngot  %q\nwant %q", got, want)
	}
}

func TestMutfa_NoLeadingHeaderPassesThrough(t *testing.T) {
	// Stray content before any header should be emitted verbatim so we don't
	// silently drop input.
	in := "stray\n>chr1\nACGT\n"
	mut := "chr1\t1\tT\n"
	want := "stray\n>chr1\nTCGT\n"

	var buf bytes.Buffer
	if err := Mutfa(strings.NewReader(in), strings.NewReader(mut), &buf); err != nil {
		t.Fatalf("Mutfa: %v", err)
	}
	if got := buf.String(); got != want {
		t.Errorf("Mutfa output:\ngot  %q\nwant %q", got, want)
	}
}

func TestMutfa_EmptyMutfile(t *testing.T) {
	in := ">chr1\nACGT\n"
	want := ">chr1\nACGT\n"

	var buf bytes.Buffer
	if err := Mutfa(strings.NewReader(in), strings.NewReader(""), &buf); err != nil {
		t.Fatalf("Mutfa: %v", err)
	}
	if got := buf.String(); got != want {
		t.Errorf("Mutfa output:\ngot  %q\nwant %q", got, want)
	}
}

func TestPickIUPAC_Table(t *testing.T) {
	// For each IUPAC code, after many draws every member of its expansion
	// should have been chosen at least once, and no non-member should appear.
	rng := rand.New(rand.NewSource(1))
	cases := []struct {
		code     byte
		expected string
	}{
		{'R', "AG"},
		{'Y', "CT"},
		{'S', "GC"},
		{'W', "AT"},
		{'K', "GT"},
		{'M', "AC"},
		{'B', "CGT"},
		{'D', "AGT"},
		{'H', "ACT"},
		{'V', "ACG"},
		{'N', "ACGT"},
	}
	for _, c := range cases {
		seen := make(map[byte]bool)
		for i := 0; i < 500; i++ {
			b := pickIUPAC(c.code, rng)
			if !strings.ContainsRune(c.expected, rune(b)) {
				t.Errorf("pickIUPAC(%c) returned %c, not in expansion %q", c.code, b, c.expected)
			}
			seen[b] = true
		}
		for _, e := range []byte(c.expected) {
			if !seen[e] {
				t.Errorf("pickIUPAC(%c): %c never sampled", c.code, e)
			}
		}
	}
}

func TestPickIUPAC_CasePreserved(t *testing.T) {
	rng := rand.New(rand.NewSource(2))
	for i := 0; i < 50; i++ {
		b := pickIUPAC('r', rng)
		if b != 'a' && b != 'g' {
			t.Errorf("pickIUPAC('r') returned %c, want 'a' or 'g'", b)
		}
	}
}

func TestPickIUPAC_NonAmbiguousUnchanged(t *testing.T) {
	rng := rand.New(rand.NewSource(3))
	for _, b := range []byte{'A', 'C', 'G', 'T', 'a', 'c', 'g', 't', '-', 'X'} {
		if got := pickIUPAC(b, rng); got != b {
			t.Errorf("pickIUPAC(%c) = %c, want unchanged", b, got)
		}
	}
}

func TestRandbase_Deterministic(t *testing.T) {
	in := ">s\nNNNNNNNN\n"
	var a, b bytes.Buffer
	if err := Randbase(strings.NewReader(in), &a, 42); err != nil {
		t.Fatalf("Randbase a: %v", err)
	}
	if err := Randbase(strings.NewReader(in), &b, 42); err != nil {
		t.Fatalf("Randbase b: %v", err)
	}
	if a.String() != b.String() {
		t.Errorf("Randbase not deterministic for same seed:\n a=%q\n b=%q", a.String(), b.String())
	}
}

func TestRandbase_PreservesLineWidth(t *testing.T) {
	in := ">s\nACGT\nNNNN\n"
	var buf bytes.Buffer
	if err := Randbase(strings.NewReader(in), &buf, 99); err != nil {
		t.Fatalf("Randbase: %v", err)
	}
	out := buf.String()
	// Two sequence lines, each 4 chars, the second all from {A,C,G,T}.
	lines := strings.Split(out, "\n")
	if len(lines) < 3 || lines[0] != ">s" || len(lines[1]) != 4 || len(lines[2]) != 4 {
		t.Fatalf("Randbase did not preserve line layout: %q", out)
	}
	// Line 1 should be unchanged (no ambiguity bases).
	if lines[1] != "ACGT" {
		t.Errorf("ACGT line was changed: got %q", lines[1])
	}
	for _, b := range []byte(lines[2]) {
		if b != 'A' && b != 'C' && b != 'G' && b != 'T' {
			t.Errorf("N-line contains non-ACGT byte %c after Randbase: %q", b, lines[2])
		}
	}
}

func TestRandbase_OnlyAmbiguityChanges(t *testing.T) {
	in := ">s\nACGTRYSWKMBDHVN-acgtryswkmbdhvn\n"
	var buf bytes.Buffer
	if err := Randbase(strings.NewReader(in), &buf, 7); err != nil {
		t.Fatalf("Randbase: %v", err)
	}
	out := buf.String()
	lines := strings.Split(out, "\n")
	if len(lines) < 2 {
		t.Fatalf("unexpected output: %q", out)
	}
	seq := lines[1]
	for i, b := range []byte(seq) {
		orig := byte(in[len(">s\n")+i])
		// Unambiguous and non-letter bytes should be unchanged.
		switch orig {
		case 'A', 'C', 'G', 'T', 'a', 'c', 'g', 't', '-':
			if b != orig {
				t.Errorf("byte %d changed unexpectedly: %c -> %c", i, orig, b)
			}
		}
		// Case must be preserved for IUPAC codes.
		if orig >= 'a' && orig <= 'z' {
			if b >= 'A' && b <= 'Z' {
				t.Errorf("lowercase IUPAC %c became uppercase %c", orig, b)
			}
		}
		if orig >= 'A' && orig <= 'Z' && orig != '-' {
			if b >= 'a' && b <= 'z' {
				t.Errorf("uppercase IUPAC %c became lowercase %c", orig, b)
			}
		}
	}
}

func TestRandbase_MultipleRecords(t *testing.T) {
	in := ">a\nNN\n>b\nRR\n"
	var buf bytes.Buffer
	if err := Randbase(strings.NewReader(in), &buf, 5); err != nil {
		t.Fatalf("Randbase: %v", err)
	}
	out := buf.String()
	lines := strings.Split(strings.TrimSuffix(out, "\n"), "\n")
	if len(lines) != 4 {
		t.Fatalf("want 4 lines, got %d: %q", len(lines), out)
	}
	if lines[0] != ">a" || lines[2] != ">b" {
		t.Errorf("record headers not preserved: %q", out)
	}
	for _, b := range []byte(lines[1]) {
		if b != 'A' && b != 'C' && b != 'G' && b != 'T' {
			t.Errorf("record a contains non-ACGT %c: %q", b, lines[1])
		}
	}
	for _, b := range []byte(lines[3]) {
		if b != 'A' && b != 'G' {
			t.Errorf("record b contains non-AG %c: %q", b, lines[3])
		}
	}
}
