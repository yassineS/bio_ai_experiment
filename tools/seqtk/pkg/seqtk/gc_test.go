package seqtk

import (
	"bytes"
	"strings"
	"testing"
)

func TestGC_TableDriven(t *testing.T) {
	// chr1: 20 A's + 30 alternating GC + 20 A's -> obvious GC island at 20-50.
	// chr2: small AT/GC alternations -> shorter region.
	// chr3: all-N -> nothing.
	const fixture = ">chr1\nAAAAAAAAAAAAAAAAAAAAGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCAAAAAAAAAAAAAAAAAAAA\n" +
		">chr2\nATATATATGCGCGCGCATATATATGCGCGCGCATATATAT\n" +
		">chr3\nNNNNNNNN\n"

	tests := []struct {
		name string
		opts GCOptions
		want string
	}{
		{
			name: "defaults find the chr1 island",
			opts: GCOptions{
				MinLength: DefaultGCMinLength,
				MinFrac:   DefaultGCMinFrac,
				XDropoff:  DefaultGCXDropoff,
			},
			want: "chr1\t20\t50\t30\n",
		},
		{
			name: "min-length too large filters chr1 out",
			opts: GCOptions{
				MinLength: 100,
				MinFrac:   DefaultGCMinFrac,
				XDropoff:  DefaultGCXDropoff,
			},
			want: "",
		},
		{
			name: "AT mode (-w) finds the A-flanks on chr1",
			opts: GCOptions{
				MinLength: 10,
				MinFrac:   0.7,
				XDropoff:  DefaultGCXDropoff,
				IsAT:      true,
			},
			want: "chr1\t0\t20\t20\nchr1\t50\t70\t20\n",
		},
	}

	for _, tc := range tests {
		t.Run(tc.name, func(t *testing.T) {
			var got bytes.Buffer
			if err := GC(strings.NewReader(fixture), &got, tc.opts); err != nil {
				t.Fatalf("GC returned error: %v", err)
			}
			if got.String() != tc.want {
				t.Fatalf("output mismatch.\nwant:\n%q\ngot:\n%q", tc.want, got.String())
			}
		})
	}
}

func TestGC_RejectsBadOptions(t *testing.T) {
	cases := []struct {
		name string
		opts GCOptions
		want string
	}{
		{"negative length", GCOptions{MinLength: 0, MinFrac: 0.6, XDropoff: 10}, "min-length"},
		{"frac 0", GCOptions{MinLength: 5, MinFrac: 0, XDropoff: 10}, "min-frac"},
		{"frac 1", GCOptions{MinLength: 5, MinFrac: 1.0, XDropoff: 10}, "min-frac"},
		{"frac >1", GCOptions{MinLength: 5, MinFrac: 1.5, XDropoff: 10}, "min-frac"},
		{"frac negative", GCOptions{MinLength: 5, MinFrac: -0.1, XDropoff: 10}, "min-frac"},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			var out bytes.Buffer
			err := GC(strings.NewReader(">x\nGCGC\n"), &out, tc.opts)
			if err == nil {
				t.Fatalf("expected error, got none")
			}
			if !strings.Contains(err.Error(), tc.want) {
				t.Fatalf("error %q should contain %q", err.Error(), tc.want)
			}
		})
	}
}

func TestGC_EmptyAndSingleRecord(t *testing.T) {
	// Empty input -> no output, no error.
	var got bytes.Buffer
	if err := GC(strings.NewReader(""), &got, GCOptions{
		MinLength: 5, MinFrac: 0.6, XDropoff: 10,
	}); err != nil {
		t.Fatalf("empty input: %v", err)
	}
	if got.Len() != 0 {
		t.Fatalf("empty input produced %q", got.String())
	}

	// A pure-GC sequence that ends without ever triggering the dropoff —
	// the final flush block must still emit it.
	got.Reset()
	if err := GC(strings.NewReader(">all\nGCGCGCGCGCGCGCGCGCGC\n"), &got, GCOptions{
		MinLength: 5, MinFrac: 0.6, XDropoff: 10,
	}); err != nil {
		t.Fatalf("pure GC: %v", err)
	}
	if got.String() != "all\t0\t20\t20\n" {
		t.Fatalf("pure-GC trailing flush wrong: %q", got.String())
	}
}

// TestGC_ParityWithUpstream byte-compares against fixtures pre-computed with
// upstream `seqtk gc` (v1.5). See gap_test.go for the parity testing pattern.
func TestGC_ParityWithUpstream(t *testing.T) {
	cases := []struct {
		name  string
		input string
		want  string
		opts  GCOptions
	}{
		{
			name:  "gc_small defaults",
			input: "gc_small.fa",
			want:  "gc_small_default.expected.bed",
			opts: GCOptions{
				MinLength: DefaultGCMinLength,
				MinFrac:   DefaultGCMinFrac,
				XDropoff:  DefaultGCXDropoff,
			},
		},
		{
			name:  "gc_small f=0.7 l=10",
			input: "gc_small.fa",
			want:  "gc_small_f07_l10.expected.bed",
			opts: GCOptions{
				MinLength: 10,
				MinFrac:   0.7,
				XDropoff:  DefaultGCXDropoff,
			},
		},
		{
			name:  "gc_small -w f=0.7 l=10",
			input: "gc_small.fa",
			want:  "gc_small_w_f07_l10.expected.bed",
			opts: GCOptions{
				MinLength: 10,
				MinFrac:   0.7,
				XDropoff:  DefaultGCXDropoff,
				IsAT:      true,
			},
		},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			in := readParityFile(t, tc.input)
			want := readParityFile(t, tc.want)
			var got bytes.Buffer
			if err := GC(bytes.NewReader(in), &got, tc.opts); err != nil {
				t.Fatalf("GC: %v", err)
			}
			mustEqualBytes(t, "gc parity "+tc.name, got.Bytes(), want)
		})
	}
}
