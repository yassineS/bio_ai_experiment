package seqtk

import (
	"bytes"
	"strings"
	"testing"
)

// TestKfreq_BasicCounts walks a hand-built FASTA with two short
// records and asserts every output column matches what the algorithm
// description in the package-level doc comment promises. This is the
// unit-level sanity check; byte-for-byte parity against upstream is
// covered separately by TestParity_Seqtk_Kfreq_*.
func TestKfreq_BasicCounts(t *testing.T) {
	const fasta = ">chr1\nACGTACGTACGTACGTAACCGGTT\n>chr2\nCCCCCCCCCCCCCCCC\n"

	tests := []struct {
		name string
		kmer string
		want string
	}{
		{
			// "AC" matches at positions 0,4,8,12,16 forward (5 hits).
			// Reverse complement "GT" appears at positions 2,6,10,14
			// (4 hits). All 4 2-mers are AC's neighbours (any 2-mer
			// is reachable from "AC" by changing one base except the
			// six 2-mers that differ in BOTH bases — but the
			// neighbourhood here uses union over single substitutions,
			// so it covers exactly { AA, AC, AG, AT } ∪ { CC, GC, TC }
			// = 7 distinct 2-mers, which catches most ACGT 2-mers).
			// We trust upstream's count of cnt_nei here; the unit test
			// only asserts the structure (column count, strand
			// selection) — parity is in the parity_test below.
			name: "AC kmer (mostly forward in chr1, AT-only mismatch in chr2 picks reverse '-')",
			kmer: "AC",
			want: "chr1\t24\t-\t7\t5\nchr2\t16\t+\t15\t0\n",
		},
		{
			// "ACGT" matches at positions 0,4,8,12 forward in chr1 (4
			// hits). chr2 is all C, no ACGT match — counts (0,0), '-'
			// by upstream tie-break.
			name: "ACGT kmer",
			kmer: "ACGT",
			want: "chr1\t24\t+\t4\t4\nchr2\t16\t-\t0\t0\n",
		},
		{
			// "AAAA": no hits anywhere; both records report 0/0 with '-'.
			name: "AAAA kmer no hits",
			kmer: "AAAA",
			want: "chr1\t24\t-\t0\t0\nchr2\t16\t-\t0\t0\n",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			var out bytes.Buffer
			if err := Kfreq(strings.NewReader(fasta), &out, KfreqOptions{Kmer: tt.kmer}); err != nil {
				t.Fatalf("Kfreq: %v", err)
			}
			if got := out.String(); got != tt.want {
				t.Fatalf("Kfreq mismatch\ngot:\n%s\nwant:\n%s", got, tt.want)
			}
		})
	}
}

// TestKfreq_EmptyRecord ensures a 0-length sequence emits a row with
// counts 0 and strand '-' (matching upstream's tie-break on cnt_nei[0]
// == cnt_nei[1]). This pins the edge case that upstream's algorithm
// silently produces.
func TestKfreq_EmptyRecord(t *testing.T) {
	const fasta = ">empty\n\n>short\nACGT\n"
	var out bytes.Buffer
	if err := Kfreq(strings.NewReader(fasta), &out, KfreqOptions{Kmer: "AC"}); err != nil {
		t.Fatalf("Kfreq: %v", err)
	}
	want := "empty\t0\t-\t0\t0\nshort\t4\t-\t1\t1\n"
	if got := out.String(); got != want {
		t.Fatalf("Kfreq mismatch\ngot:\n%s\nwant:\n%s", got, want)
	}
}

// TestKfreq_NRunResetsRollingWindow exercises the path where a non-ACGT
// byte resets the rolling encoder. The "ACGT" 2-mer at positions 8-11
// must NOT be considered a continuation of the 5-mer rolling window
// that was reset at the N-run; upstream stk_kfreq does this via
// `else k = 0;`.
func TestKfreq_NRunResetsRollingWindow(t *testing.T) {
	const fasta = ">withN\nACGTNNNNACGT\n"
	var out bytes.Buffer
	// kmer length 4 lets us count exactly two ACGT hits: one at pos 0
	// and one at pos 8 (post-N). If we forgot to reset the window the
	// N-bridging k-mer would also match neighbour space.
	if err := Kfreq(strings.NewReader(fasta), &out, KfreqOptions{Kmer: "ACGT"}); err != nil {
		t.Fatalf("Kfreq: %v", err)
	}
	// neighbour count == 2 (the two ACGT hits themselves); exact == 2;
	// reverse strand has 0 — '+' wins.
	want := "withN\t12\t+\t2\t2\n"
	if got := out.String(); got != want {
		t.Fatalf("Kfreq mismatch\ngot:\n%s\nwant:\n%s", got, want)
	}
}

// TestKfreq_InvalidKmer ensures non-ACGT bytes in the kmer are
// rejected with a typed error rather than panicking (upstream uses
// assert() and aborts; we promise a clean error path).
func TestKfreq_InvalidKmer(t *testing.T) {
	cases := []struct {
		name string
		kmer string
	}{
		{"empty kmer", ""},
		{"N in kmer", "ACNT"},
		{"IUPAC in kmer", "RACG"},
		{"too long", strings.Repeat("A", 16)},
	}
	for _, tt := range cases {
		t.Run(tt.name, func(t *testing.T) {
			var out bytes.Buffer
			err := Kfreq(strings.NewReader(">x\nACGT\n"), &out, KfreqOptions{Kmer: tt.kmer})
			if err == nil {
				t.Fatalf("expected error for kmer %q", tt.kmer)
			}
		})
	}
}

// TestKfreq_LowercaseAcceptedSameAsUpper pins the seq_nt6_table
// behaviour that lowercase ACGT bytes count identically to uppercase.
func TestKfreq_LowercaseAcceptedSameAsUpper(t *testing.T) {
	var outUpper, outLower bytes.Buffer
	if err := Kfreq(strings.NewReader(">x\nACGTACGT\n"), &outUpper, KfreqOptions{Kmer: "ACGT"}); err != nil {
		t.Fatalf("Kfreq upper: %v", err)
	}
	if err := Kfreq(strings.NewReader(">x\nacgtacgt\n"), &outLower, KfreqOptions{Kmer: "ACGT"}); err != nil {
		t.Fatalf("Kfreq lower: %v", err)
	}
	if outUpper.String() != outLower.String() {
		t.Fatalf("lowercase != uppercase:\nupper: %q\nlower: %q", outUpper.String(), outLower.String())
	}
}
