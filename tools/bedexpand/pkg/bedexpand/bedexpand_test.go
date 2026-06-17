package bedexpand

import (
	"bytes"
	"reflect"
	"strings"
	"testing"
)

func TestParseColumns(t *testing.T) {
	cases := []struct {
		in      string
		want    []int
		wantErr bool
	}{
		{"4", []int{4}, false},
		{"4,5", []int{4, 5}, false},
		{"5,4", []int{5, 4}, false},
		{" 5 , 4 ", []int{5, 4}, false},
		{"", nil, true},
		{",", nil, true},
		{"4,", nil, true},
		{"x", nil, true},
		{"0", nil, true},
		{"-3", nil, true},
	}
	for _, tc := range cases {
		got, err := ParseColumns(tc.in)
		if (err != nil) != tc.wantErr {
			t.Errorf("ParseColumns(%q): err=%v wantErr=%v", tc.in, err, tc.wantErr)
			continue
		}
		if err == nil && !reflect.DeepEqual(got, tc.want) {
			t.Errorf("ParseColumns(%q) = %v, want %v", tc.in, got, tc.want)
		}
	}
}

func TestExpand_SingleColumn(t *testing.T) {
	in := "chr1\t10\t20\t1,2,3\t10,20,30\nchr1\t40\t50\t4,5,6\t40,50,60\n"
	want := "chr1\t10\t20\t1\t10,20,30\n" +
		"chr1\t10\t20\t2\t10,20,30\n" +
		"chr1\t10\t20\t3\t10,20,30\n" +
		"chr1\t40\t50\t4\t40,50,60\n" +
		"chr1\t40\t50\t5\t40,50,60\n" +
		"chr1\t40\t50\t6\t40,50,60\n"
	var buf bytes.Buffer
	n, err := Expand(strings.NewReader(in), &buf, Options{Columns: []int{4}})
	if err != nil {
		t.Fatalf("Expand: %v", err)
	}
	if n != 6 {
		t.Errorf("written = %d, want 6", n)
	}
	if buf.String() != want {
		t.Errorf("output mismatch.\nwant:\n%s\ngot:\n%s", want, buf.String())
	}
}

func TestExpand_MultipleColumns(t *testing.T) {
	in := "chr1\t10\t20\t1,2,3\t10,20,30\nchr1\t40\t50\t4,5,6\t40,50,60\n"
	want := "chr1\t10\t20\t1\t10\n" +
		"chr1\t10\t20\t2\t20\n" +
		"chr1\t10\t20\t3\t30\n" +
		"chr1\t40\t50\t4\t40\n" +
		"chr1\t40\t50\t5\t50\n" +
		"chr1\t40\t50\t6\t60\n"
	var buf bytes.Buffer
	if _, err := Expand(strings.NewReader(in), &buf, Options{Columns: []int{4, 5}}); err != nil {
		t.Fatalf("Expand: %v", err)
	}
	if buf.String() != want {
		t.Errorf("output mismatch.\nwant:\n%s\ngot:\n%s", want, buf.String())
	}
}

func TestExpand_SwappedColumns(t *testing.T) {
	// Upstream expand.t3: -c 5,4 substitutes column-5 elements at position 4
	// and column-4 elements at position 5.
	in := "chr1\t10\t20\t1,2,3\t10,20,30\nchr1\t40\t50\t4,5,6\t40,50,60\n"
	want := "chr1\t10\t20\t10\t1\n" +
		"chr1\t10\t20\t20\t2\n" +
		"chr1\t10\t20\t30\t3\n" +
		"chr1\t40\t50\t40\t4\n" +
		"chr1\t40\t50\t50\t5\n" +
		"chr1\t40\t50\t60\t6\n"
	var buf bytes.Buffer
	if _, err := Expand(strings.NewReader(in), &buf, Options{Columns: []int{5, 4}}); err != nil {
		t.Fatalf("Expand: %v", err)
	}
	if buf.String() != want {
		t.Errorf("output mismatch.\nwant:\n%s\ngot:\n%s", want, buf.String())
	}
}

func TestExpand_HeaderLinesPassthrough(t *testing.T) {
	in := "# header\ntrack name=foo\nchr1\t1\t2\ta,b\n"
	want := "# header\ntrack name=foo\nchr1\t1\t2\ta\nchr1\t1\t2\tb\n"
	var buf bytes.Buffer
	if _, err := Expand(strings.NewReader(in), &buf, Options{Columns: []int{4}}); err != nil {
		t.Fatalf("Expand: %v", err)
	}
	if buf.String() != want {
		t.Errorf("output mismatch.\nwant:\n%s\ngot:\n%s", want, buf.String())
	}
}

func TestExpand_EmptyLines(t *testing.T) {
	in := "\nchr1\t0\t1\tx,y\n\n"
	want := "\nchr1\t0\t1\tx\nchr1\t0\t1\ty\n\n"
	var buf bytes.Buffer
	if _, err := Expand(strings.NewReader(in), &buf, Options{Columns: []int{4}}); err != nil {
		t.Fatalf("Expand: %v", err)
	}
	if buf.String() != want {
		t.Errorf("output mismatch.\nwant:\n%s\ngot:\n%s", want, buf.String())
	}
}

func TestExpand_SingleElement(t *testing.T) {
	// No comma → one element per list → one output row per input row.
	in := "chr1\t1\t2\tx\ty\n"
	want := "chr1\t1\t2\tx\ty\n"
	var buf bytes.Buffer
	if _, err := Expand(strings.NewReader(in), &buf, Options{Columns: []int{4}}); err != nil {
		t.Fatalf("Expand: %v", err)
	}
	if buf.String() != want {
		t.Errorf("output mismatch.\nwant:\n%s\ngot:\n%s", want, buf.String())
	}
}

// TestUnit_tokenizeCSV checks the comma tokenizer matches C++ getline
// semantics: a single terminating empty (trailing comma) is dropped, while
// leading and interior empties are preserved, and an empty cell yields no
// elements.
func TestUnit_tokenizeCSV(t *testing.T) {
	cases := []struct {
		in   string
		want []string
	}{
		{"", nil},
		{"a", []string{"a"}},
		{"a,b,c", []string{"a", "b", "c"}},
		{"a,b,c,", []string{"a", "b", "c"}}, // trailing comma dropped
		{"10,20,30,", []string{"10", "20", "30"}},
		{",a", []string{"", "a"}},             // leading empty kept
		{"a,,b", []string{"a", "", "b"}},      // interior empty kept
		{",a,,b", []string{"", "a", "", "b"}}, // leading + interior kept
		{"a,,", []string{"a", ""}},            // only the final empty dropped
		{",", []string{""}},                   // single delimiter: one empty token
	}
	for _, c := range cases {
		got := tokenizeCSV(c.in)
		if !reflect.DeepEqual(got, c.want) {
			t.Errorf("tokenizeCSV(%q) = %#v, want %#v", c.in, got, c.want)
		}
	}
}

func TestExpand_EmptyInput(t *testing.T) {
	var buf bytes.Buffer
	n, err := Expand(strings.NewReader(""), &buf, Options{Columns: []int{4}})
	if err != nil {
		t.Fatalf("Expand: %v", err)
	}
	if n != 0 || buf.Len() != 0 {
		t.Errorf("empty input should produce no output: n=%d out=%q", n, buf.String())
	}
}

func TestExpand_NoColumns(t *testing.T) {
	var buf bytes.Buffer
	_, err := Expand(strings.NewReader("a\tb\tc\n"), &buf, Options{})
	if err == nil {
		t.Errorf("expected error for empty Columns")
	}
}

func TestExpand_BadColumnIndex(t *testing.T) {
	var buf bytes.Buffer
	_, err := Expand(strings.NewReader("a\tb\tc\n"), &buf, Options{Columns: []int{0}})
	if err == nil {
		t.Errorf("expected error for 0-based column")
	}
}

func TestExpand_OutOfBoundsColumn(t *testing.T) {
	var buf bytes.Buffer
	_, err := Expand(strings.NewReader("a\tb\tc\n"), &buf, Options{Columns: []int{99}})
	if err == nil {
		t.Errorf("expected error for column 99 in 3-col row")
	}
}

func TestExpand_ListLengthMismatch(t *testing.T) {
	var buf bytes.Buffer
	_, err := Expand(strings.NewReader("chr1\t0\t1\ta,b\tx,y,z\n"), &buf, Options{Columns: []int{4, 5}})
	if err == nil {
		t.Errorf("expected error for mismatched list lengths")
	}
}

// errReader returns an error on the first Read; used to exercise scanner.Err
// in Expand.
type errReader struct{}

func (errReader) Read(_ []byte) (int, error) { return 0, errReadFail }

var errReadFail = &readErr{"forced"}

type readErr struct{ s string }

func (e *readErr) Error() string { return e.s }

func TestExpand_ScannerError(t *testing.T) {
	var buf bytes.Buffer
	_, err := Expand(errReader{}, &buf, Options{Columns: []int{1}})
	if err == nil {
		t.Errorf("expected scanner error")
	}
}
