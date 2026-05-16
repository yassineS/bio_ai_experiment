package seqtk

import (
	"bytes"
	"strings"
	"testing"
)

func TestGap_TableDriven(t *testing.T) {
	// Fixture used across cases: chr1 has two N-runs (lengths 5 and 10),
	// chr2 is all-N (length 10), chr3 has a 3-byte IUPAC run (RRR).
	const fixture = ">chr1\nACGTNNNNNACGTNNACGTNNNNNNNNNNGGAA\n" +
		">chr2\nNNNNNNNNNNACGTACGT\n" +
		">chr3\nACGTRRRACGT\n" +
		">empty\n\n"

	tests := []struct {
		name    string
		minSize int
		want    string
	}{
		{
			name:    "default-ish min 5 keeps long Ns",
			minSize: 5,
			want:    "chr1\t4\t9\nchr1\t19\t29\nchr2\t0\t10\n",
		},
		{
			name:    "min 3 picks up the IUPAC R run",
			minSize: 3,
			want:    "chr1\t4\t9\nchr1\t19\t29\nchr2\t0\t10\nchr3\t4\t7\n",
		},
		{
			name:    "min 11 filters everything out",
			minSize: 11,
			want:    "",
		},
		{
			name:    "min 1 catches every non-ACGT byte including 2-byte runs",
			minSize: 1,
			want:    "chr1\t4\t9\nchr1\t13\t15\nchr1\t19\t29\nchr2\t0\t10\nchr3\t4\t7\n",
		},
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			var got bytes.Buffer
			if err := Gap(strings.NewReader(fixture), &got, GapOptions{MinSize: tc.minSize}); err != nil {
				t.Fatalf("Gap returned error: %v", err)
			}
			if got.String() != tc.want {
				t.Fatalf("output mismatch.\nwant:\n%q\ngot:\n%q", tc.want, got.String())
			}
		})
	}
}

func TestGap_RejectsBadMinSize(t *testing.T) {
	for _, n := range []int{0, -1, -100} {
		var out bytes.Buffer
		err := Gap(strings.NewReader(">x\nACGT\n"), &out, GapOptions{MinSize: n})
		if err == nil {
			t.Fatalf("MinSize=%d: expected error, got none", n)
		}
		if !strings.Contains(err.Error(), "min-size") {
			t.Fatalf("MinSize=%d: error message should mention min-size, got %v", n, err)
		}
	}
}

func TestGap_TrailingGapFlushed(t *testing.T) {
	// A gap that runs right to end-of-sequence must still be reported —
	// upstream's `i == len` sentinel branch handles this.
	const fixture = ">chr1\nACGTNNNNN\n"
	var got bytes.Buffer
	if err := Gap(strings.NewReader(fixture), &got, GapOptions{MinSize: 3}); err != nil {
		t.Fatalf("Gap returned error: %v", err)
	}
	if got.String() != "chr1\t4\t9\n" {
		t.Fatalf("trailing gap not flushed: %q", got.String())
	}
}

// TestGap_ParityWithUpstream byte-compares against fixtures pre-computed by
// running upstream `seqtk gap` on the same inputs (see
// reference_code/seqtk/seqtk v1.5). Fixtures live under
// tools/seqtk/testdata/parity/ alongside the input FASTA.
func TestGap_ParityWithUpstream(t *testing.T) {
	cases := []struct {
		name    string
		input   string
		want    string
		minSize int
	}{
		{
			name:    "gap_small l=5",
			input:   "gap_small.fa",
			want:    "gap_small_l5.expected.bed",
			minSize: 5,
		},
		{
			name:    "gap_small l=3 (picks up RRR)",
			input:   "gap_small.fa",
			want:    "gap_small_l3.expected.bed",
			minSize: 3,
		},
		{
			name:    "nruns.fa l=1",
			input:   "nruns.fa",
			want:    "gap_nruns_l1.expected.bed",
			minSize: 1,
		},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			in := readParityFile(t, tc.input)
			want := readParityFile(t, tc.want)
			var got bytes.Buffer
			if err := Gap(bytes.NewReader(in), &got, GapOptions{MinSize: tc.minSize}); err != nil {
				t.Fatalf("Gap: %v", err)
			}
			mustEqualBytes(t, "gap parity "+tc.name, got.Bytes(), want)
		})
	}
}
