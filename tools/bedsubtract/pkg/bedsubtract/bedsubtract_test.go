package bedsubtract

import (
	"bytes"
	"strings"
	"testing"
)

// runSubtract is a small helper that runs Subtract on two input strings and
// returns the resulting output (or an error).
func runSubtract(t *testing.T, a, b string, opts Options) (string, int, error) {
	t.Helper()
	var buf bytes.Buffer
	n, err := Subtract(strings.NewReader(a), strings.NewReader(b), &buf, opts)
	return buf.String(), n, err
}

func TestSubtractBasic(t *testing.T) {
	tests := []struct {
		name    string
		a, b    string
		opts    Options
		want    string
		wantN   int
		wantErr bool
	}{
		{
			name:  "no_overlap",
			a:     "chr1\t10\t20\nchr1\t30\t40\n",
			b:     "chr1\t50\t60\n",
			want:  "chr1\t10\t20\nchr1\t30\t40\n",
			wantN: 2,
		},
		{
			name:  "full_cover_drops_a",
			a:     "chr1\t10\t20\n",
			b:     "chr1\t5\t30\n",
			want:  "",
			wantN: 0,
		},
		{
			name:  "left_trim",
			a:     "chr1\t10\t20\n",
			b:     "chr1\t5\t15\n",
			want:  "chr1\t15\t20\n",
			wantN: 1,
		},
		{
			name:  "right_trim",
			a:     "chr1\t10\t20\n",
			b:     "chr1\t15\t25\n",
			want:  "chr1\t10\t15\n",
			wantN: 1,
		},
		{
			name:  "split_middle",
			a:     "chr1\t10\t30\n",
			b:     "chr1\t15\t20\n",
			want:  "chr1\t10\t15\nchr1\t20\t30\n",
			wantN: 2,
		},
		{
			name:  "multiple_holes",
			a:     "chr1\t10\t50\n",
			b:     "chr1\t15\t20\nchr1\t30\t35\n",
			want:  "chr1\t10\t15\nchr1\t20\t30\nchr1\t35\t50\n",
			wantN: 3,
		},
		{
			name:  "different_chrom_no_subtract",
			a:     "chr1\t10\t20\n",
			b:     "chr2\t10\t20\n",
			want:  "chr1\t10\t20\n",
			wantN: 1,
		},
		{
			name:  "preserve_bed6_columns",
			a:     "chr1\t10\t30\tnameA\t100\t+\n",
			b:     "chr1\t15\t20\n",
			want:  "chr1\t10\t15\tnameA\t100\t+\nchr1\t20\t30\tnameA\t100\t+\n",
			wantN: 2,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, n, err := runSubtract(t, tt.a, tt.b, tt.opts)
			if (err != nil) != tt.wantErr {
				t.Fatalf("err = %v, wantErr = %v", err, tt.wantErr)
			}
			if got != tt.want {
				t.Errorf("output mismatch.\nwant:\n%q\ngot:\n%q", tt.want, got)
			}
			if n != tt.wantN {
				t.Errorf("count = %d, want %d", n, tt.wantN)
			}
		})
	}
}

func TestSubtractRemoveEntire(t *testing.T) {
	a := "chr1\t10\t30\nchr1\t100\t200\n"
	b := "chr1\t20\t25\n"
	got, n, err := runSubtract(t, a, b, Options{RemoveEntire: true})
	if err != nil {
		t.Fatal(err)
	}
	want := "chr1\t100\t200\n"
	if got != want {
		t.Errorf("want %q got %q", want, got)
	}
	if n != 1 {
		t.Errorf("count = %d, want 1", n)
	}
}

func TestSubtractMinFraction(t *testing.T) {
	a := "chr1\t0\t100\n"
	// b covers 10/100 = 10% of A. Default 0 -> subtract; 0.5 -> ignore.
	b := "chr1\t10\t20\n"

	got, _, err := runSubtract(t, a, b, Options{})
	if err != nil {
		t.Fatal(err)
	}
	if got != "chr1\t0\t10\nchr1\t20\t100\n" {
		t.Errorf("got %q", got)
	}

	got2, _, err := runSubtract(t, a, b, Options{MinFraction: 0.5})
	if err != nil {
		t.Fatal(err)
	}
	if got2 != "chr1\t0\t100\n" {
		t.Errorf("min-fraction 0.5 should keep A intact, got %q", got2)
	}
}

func TestSubtractStrandSame(t *testing.T) {
	a := "chr1\t10\t30\ta1\t0\t+\n"
	b := "chr1\t15\t20\tb1\t0\t-\nchr1\t22\t25\tb2\t0\t+\n"
	got, _, err := runSubtract(t, a, b, Options{SameStrand: true})
	if err != nil {
		t.Fatal(err)
	}
	// Only the '+'-strand B (22..25) is subtracted; '-'-strand B is ignored.
	want := "chr1\t10\t22\ta1\t0\t+\nchr1\t25\t30\ta1\t0\t+\n"
	if got != want {
		t.Errorf("want %q got %q", want, got)
	}
}

func TestSubtractStrandOpposite(t *testing.T) {
	a := "chr1\t10\t30\ta1\t0\t+\n"
	b := "chr1\t15\t20\tb1\t0\t-\nchr1\t22\t25\tb2\t0\t+\n"
	got, _, err := runSubtract(t, a, b, Options{OppositeStrand: true})
	if err != nil {
		t.Fatal(err)
	}
	// Only the '-'-strand B is subtracted.
	want := "chr1\t10\t15\ta1\t0\t+\nchr1\t20\t30\ta1\t0\t+\n"
	if got != want {
		t.Errorf("want %q got %q", want, got)
	}
}

func TestSubtractStrandMissingDropped(t *testing.T) {
	a := "chr1\t10\t30\n"
	b := "chr1\t15\t20\n"
	// BED3 has no strand; with -s nothing should be subtracted.
	got, _, err := runSubtract(t, a, b, Options{SameStrand: true})
	if err != nil {
		t.Fatal(err)
	}
	if got != "chr1\t10\t30\n" {
		t.Errorf("expected A unchanged, got %q", got)
	}
}

func TestSubtractValidate(t *testing.T) {
	bad := []Options{
		{SameStrand: true, OppositeStrand: true},
		{MinFraction: -0.1},
		{MinFraction: 1.5},
	}
	for _, o := range bad {
		if err := o.Validate(); err == nil {
			t.Errorf("expected validation error for %+v", o)
		}
	}
	if err := (Options{MinFraction: 0.5}.Validate()); err != nil {
		t.Errorf("unexpected validation error: %v", err)
	}
}

func TestSubtractErrors(t *testing.T) {
	// Invalid start in A.
	_, _, err := runSubtract(t, "chr1\tNOPE\t10\n", "", Options{})
	if err == nil {
		t.Error("expected error for invalid start")
	}
	// Too few fields.
	_, _, err = runSubtract(t, "chr1\t10\n", "", Options{})
	if err == nil {
		t.Error("expected error for too few fields")
	}
	// Invalid end.
	_, _, err = runSubtract(t, "chr1\t0\tBAD\n", "", Options{})
	if err == nil {
		t.Error("expected error for invalid end")
	}
	// end < start.
	_, _, err = runSubtract(t, "chr1\t10\t5\n", "", Options{})
	if err == nil {
		t.Error("expected error for end < start")
	}
	// Bad options surface via Subtract too.
	_, _, err = runSubtract(t, "", "", Options{MinFraction: 2})
	if err == nil {
		t.Error("expected validation error from Subtract")
	}
}

func TestSubtractSkipsCommentsAndBlankLines(t *testing.T) {
	a := "# comment\n\ntrack name=foo\nbrowser hide all\nchr1\t10\t20\n"
	b := ""
	got, n, err := runSubtract(t, a, b, Options{})
	if err != nil {
		t.Fatal(err)
	}
	if got != "chr1\t10\t20\n" || n != 1 {
		t.Errorf("got %q n=%d", got, n)
	}
}

func TestSubtractSortsBByStart(t *testing.T) {
	// B comes in unsorted; subtraction should still be correct.
	a := "chr1\t0\t100\n"
	b := "chr1\t60\t70\nchr1\t10\t20\nchr1\t30\t40\n"
	got, n, err := runSubtract(t, a, b, Options{})
	if err != nil {
		t.Fatal(err)
	}
	want := "chr1\t0\t10\nchr1\t20\t30\nchr1\t40\t60\nchr1\t70\t100\n"
	if got != want || n != 4 {
		t.Errorf("got %q n=%d", got, n)
	}
}

func TestSubtractRemoveEntireDoesNothingWhenNoOverlap(t *testing.T) {
	a := "chr1\t0\t10\n"
	b := "chr1\t20\t30\n"
	got, _, err := runSubtract(t, a, b, Options{RemoveEntire: true})
	if err != nil {
		t.Fatal(err)
	}
	if got != "chr1\t0\t10\n" {
		t.Errorf("got %q", got)
	}
}
