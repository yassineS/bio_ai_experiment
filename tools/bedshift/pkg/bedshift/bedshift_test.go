package bedshift

import (
	"bytes"
	"strings"
	"testing"
)

func sizes(m map[string]int64) ChromSizes { return ChromSizes(m) }

func TestReadChromSizes(t *testing.T) {
	in := "chr1\t1000\nchr2 2000 extra\n# comment\n\nchr3\t3000\n"
	cs, err := ReadChromSizes(strings.NewReader(in))
	if err != nil {
		t.Fatalf("ReadChromSizes: %v", err)
	}
	want := map[string]int64{"chr1": 1000, "chr2": 2000, "chr3": 3000}
	for k, v := range want {
		if cs[k] != v {
			t.Fatalf("%s = %d, want %d", k, cs[k], v)
		}
	}
}

func runShift(t *testing.T, in string, cs ChromSizes, opts Options, header bool) string {
	t.Helper()
	var out bytes.Buffer
	if _, err := Shift(strings.NewReader(in), &out, cs, opts, header); err != nil {
		t.Fatalf("Shift: %v", err)
	}
	return out.String()
}

func TestShift(t *testing.T) {
	tiny := sizes(map[string]int64{"chr1": 1000})
	a := "chr1\t100\t200\ta1\t1\t+\nchr1\t100\t200\ta2\t2\t-\n"

	tests := []struct {
		name   string
		in     string
		cs     ChromSizes
		opts   Options
		header bool
		want   string
	}{
		{
			name: "s5-forward",
			in:   a, cs: tiny, opts: Options{ShiftPlus: 5, ShiftMinus: 5},
			want: "chr1\t105\t205\ta1\t1\t+\nchr1\t105\t205\ta2\t2\t-\n",
		},
		{
			name: "s-5-backward",
			in:   a, cs: tiny, opts: Options{ShiftPlus: -5, ShiftMinus: -5},
			want: "chr1\t95\t195\ta1\t1\t+\nchr1\t95\t195\ta2\t2\t-\n",
		},
		{
			name: "m-5-p0",
			in:   a, cs: tiny, opts: Options{ShiftPlus: 0, ShiftMinus: -5},
			want: "chr1\t100\t200\ta1\t1\t+\nchr1\t95\t195\ta2\t2\t-\n",
		},
		{
			name: "p5-m0",
			in:   a, cs: tiny, opts: Options{ShiftPlus: 5, ShiftMinus: 0},
			want: "chr1\t105\t205\ta1\t1\t+\nchr1\t100\t200\ta2\t2\t-\n",
		},
		{
			name: "past-start",
			in:   a, cs: tiny, opts: Options{ShiftPlus: -200, ShiftMinus: -200},
			want: "chr1\t0\t1\ta1\t1\t+\nchr1\t0\t1\ta2\t2\t-\n",
		},
		{
			name: "past-end",
			in:   a, cs: tiny, opts: Options{ShiftPlus: 1000, ShiftMinus: 1000},
			want: "chr1\t999\t1000\ta1\t1\t+\nchr1\t999\t1000\ta2\t2\t-\n",
		},
		{
			name: "fractional-truncates-toward-zero",
			in:   "chr1\t100\t200\tx\t0\t+\n", cs: tiny,
			opts: Options{ShiftPlus: -0.123, ShiftMinus: -0.123, Fractional: true},
			want: "chr1\t87\t187\tx\t0\t+\n",
		},
		{
			name: "header-passthrough",
			in:   "# header line\ntrack name=foo\nchr1\t100\t200\tx\t0\t+\n", cs: tiny,
			opts: Options{ShiftPlus: 5, ShiftMinus: 5}, header: true,
			want: "# header line\ntrack name=foo\nchr1\t105\t205\tx\t0\t+\n",
		},
		{
			name: "header-dropped-without-flag",
			in:   "# header line\nchr1\t100\t200\tx\t0\t+\n", cs: tiny,
			opts: Options{ShiftPlus: 5, ShiftMinus: 5},
			want: "chr1\t105\t205\tx\t0\t+\n",
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := runShift(t, tt.in, tt.cs, tt.opts, tt.header)
			if got != tt.want {
				t.Fatalf("got %q, want %q", got, tt.want)
			}
		})
	}
}

func TestShiftMissingChromMatchesUpstreamBug(t *testing.T) {
	// Upstream getChromSize returns (uint64_t)-1 reinterpreted as -1 for an
	// unknown chromosome, producing negative output coordinates. We replicate
	// that for byte parity.
	got := runShift(t, "chrZ\t10\t20\ty\t0\t+\n", sizes(map[string]int64{"chr1": 1000}),
		Options{ShiftPlus: 5, ShiftMinus: 5}, false)
	want := "chrZ\t-2\t-1\ty\t0\t+\n"
	if got != want {
		t.Fatalf("got %q, want %q", got, want)
	}
}

func TestShiftHugeShift(t *testing.T) {
	// -s larger than a 32-bit int still clamps to the chrom bounds.
	got := runShift(t, "chr1\t100\t200\ta1\t1\t+\n", sizes(map[string]int64{"chr1": 1000}),
		Options{ShiftPlus: 3000000000, ShiftMinus: 3000000000}, false)
	want := "chr1\t999\t1000\ta1\t1\t+\n"
	if got != want {
		t.Fatalf("got %q, want %q", got, want)
	}
}
