package bedoverlap

import (
	"bytes"
	"strings"
	"testing"
)

func TestParseCols(t *testing.T) {
	tests := []struct {
		name    string
		spec    string
		want    Cols
		wantErr bool
	}{
		{"basic", "2,3,6,7", Cols{2, 3, 6, 7}, false},
		{"spaces", " 1 , 2 , 3 , 4 ", Cols{1, 2, 3, 4}, false},
		{"too-few", "2,3,6", Cols{}, true},
		{"too-many", "2,3,6,7,8", Cols{}, true},
		{"non-int", "a,b,c,d", Cols{}, true},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := ParseCols(tt.spec)
			if tt.wantErr {
				if err == nil {
					t.Fatalf("expected error, got %v", got)
				}
				return
			}
			if err != nil {
				t.Fatalf("unexpected error: %v", err)
			}
			if got != tt.want {
				t.Fatalf("got %+v, want %+v", got, tt.want)
			}
		})
	}
}

func TestOverlap(t *testing.T) {
	cols := Cols{2, 3, 6, 7}
	tests := []struct {
		name  string
		input string
		want  string
	}{
		{
			name:  "overlap-and-distance",
			input: "chr1\t10\t20\tA\tchr1\t15\t25\tB\nchr1\t10\t20\tC\tchr1\t25\t35\tD\n",
			want:  "chr1\t10\t20\tA\tchr1\t15\t25\tB\t5\nchr1\t10\t20\tC\tchr1\t25\t35\tD\t-5\n",
		},
		{
			name:  "touching-zero",
			input: "chr1\t10\t20\tx\tchr1\t20\t30\ty\n",
			want:  "chr1\t10\t20\tx\tchr1\t20\t30\ty\t0\n",
		},
		{
			name:  "nested",
			input: "chr1\t10\t100\tx\tchr1\t40\t50\ty\n",
			want:  "chr1\t10\t100\tx\tchr1\t40\t50\ty\t10\n",
		},
		{
			name:  "disjoint-negative",
			input: "chr1\t10\t20\tx\tchr1\t50\t60\ty\n",
			want:  "chr1\t10\t20\tx\tchr1\t50\t60\ty\t-30\n",
		},
		{
			name:  "single-field-skipped",
			input: "loneword\nchr1\t10\t20\tx\tchr1\t15\t25\ty\n",
			want:  "chr1\t10\t20\tx\tchr1\t15\t25\ty\t5\n",
		},
		{
			name:  "negative-coordinate",
			input: "chr1\t-10\t20\tx\tchr1\t5\t25\ty\n",
			want:  "chr1\t-10\t20\tx\tchr1\t5\t25\ty\t15\n",
		},
		{
			name:  "empty",
			input: "",
			want:  "",
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			var out bytes.Buffer
			if err := Overlap(strings.NewReader(tt.input), &out, cols); err != nil {
				t.Fatalf("Overlap error: %v", err)
			}
			if out.String() != tt.want {
				t.Fatalf("got %q, want %q", out.String(), tt.want)
			}
		})
	}
}

func TestOverlapNonNumeric(t *testing.T) {
	var out bytes.Buffer
	err := Overlap(strings.NewReader("chr1\tfoo\t20\tx\tchr1\t50\t60\ty\n"), &out, Cols{2, 3, 6, 7})
	if err == nil {
		t.Fatalf("expected error for non-numeric column")
	}
	if !strings.Contains(err.Error(), "non-numeric at line 1") {
		t.Fatalf("unexpected error text: %v", err)
	}
}
