package bedmerge

import (
	"bytes"
	"strings"
	"testing"
)

func TestMergeBasic(t *testing.T) {
	input := `chr1	100	200
chr1	150	250
chr1	300	400`

	expected := `chr1	100	250
chr1	300	400
`

	reader := strings.NewReader(input)
	var buf bytes.Buffer

	count, err := Merge(reader, &buf, MergeOptions{})
	if err != nil {
		t.Fatalf("Merge failed: %v", err)
	}

	if count != 2 {
		t.Errorf("Expected 2 merged intervals, got %d", count)
	}

	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected:\n%s\nGot:\n%s", expected, buf.String())
	}
}

func TestMergeStrandFilter(t *testing.T) {
	// Mixed strands: only the requested strand survives, then a positional
	// merge runs over the survivors (BED3 output).
	input := "chr1\t10\t50\ta\t0\t+\n" +
		"chr1\t20\t60\tb\t0\t-\n" +
		"chr1\t40\t80\tc\t0\t+\n" +
		"chr1\t30\t70\td\t0\t.\n"

	plus, err := mergeToString(input, MergeOptions{StrandFilter: "+"})
	if err != nil {
		t.Fatalf("Merge(+): %v", err)
	}
	if plus != "chr1\t10\t80\n" {
		t.Errorf("-S +: want chr1 10 80, got %q", plus)
	}

	minus, err := mergeToString(input, MergeOptions{StrandFilter: "-"})
	if err != nil {
		t.Fatalf("Merge(-): %v", err)
	}
	if minus != "chr1\t20\t60\n" {
		t.Errorf("-S -: want chr1 20 60, got %q", minus)
	}
}

func TestMergeStrandFilterValidation(t *testing.T) {
	// -S with a bad argument is rejected.
	if _, err := Merge(strings.NewReader(""), &bytes.Buffer{}, MergeOptions{StrandFilter: "x"}); err == nil {
		t.Error("expected error for invalid -S argument")
	}
	// -s and -S together are rejected.
	if _, err := Merge(strings.NewReader(""), &bytes.Buffer{}, MergeOptions{StrandFilter: "+", StrandSpec: true}); err == nil {
		t.Error("expected error for -s and -S together")
	}
}

// mergeToString runs Merge on an input string and returns the output.
func mergeToString(input string, opts MergeOptions) (string, error) {
	var buf bytes.Buffer
	if _, err := Merge(strings.NewReader(input), &buf, opts); err != nil {
		return "", err
	}
	return buf.String(), nil
}

func TestMergeAdjacent(t *testing.T) {
	input := `chr1	100	200
chr1	200	300
chr1	300	400`

	expected := `chr1	100	400
`

	reader := strings.NewReader(input)
	var buf bytes.Buffer

	count, err := Merge(reader, &buf, MergeOptions{})
	if err != nil {
		t.Fatalf("Merge failed: %v", err)
	}

	if count != 1 {
		t.Errorf("Expected 1 merged interval, got %d", count)
	}

	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected:\n%s\nGot:\n%s", expected, buf.String())
	}
}

func TestMergeNonOverlapping(t *testing.T) {
	input := `chr1	100	200
chr1	300	400
chr1	500	600`

	expected := `chr1	100	200
chr1	300	400
chr1	500	600
`

	reader := strings.NewReader(input)
	var buf bytes.Buffer

	count, err := Merge(reader, &buf, MergeOptions{})
	if err != nil {
		t.Fatalf("Merge failed: %v", err)
	}

	if count != 3 {
		t.Errorf("Expected 3 intervals, got %d", count)
	}

	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected:\n%s\nGot:\n%s", expected, buf.String())
	}
}

func TestMergeMaxDistance(t *testing.T) {
	input := `chr1	100	200
chr1	250	300
chr1	500	600`

	expected := `chr1	100	300
chr1	500	600
`

	reader := strings.NewReader(input)
	var buf bytes.Buffer

	// Merge intervals within 50bp
	count, err := Merge(reader, &buf, MergeOptions{MaxDistance: 50})
	if err != nil {
		t.Fatalf("Merge failed: %v", err)
	}

	if count != 2 {
		t.Errorf("Expected 2 merged intervals, got %d", count)
	}

	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected:\n%s\nGot:\n%s", expected, buf.String())
	}
}

func TestMergeMultipleChromosomes(t *testing.T) {
	input := `chr1	100	200
chr1	150	250
chr2	100	200
chr2	150	250
chr3	100	200`

	expected := `chr1	100	250
chr2	100	250
chr3	100	200
`

	reader := strings.NewReader(input)
	var buf bytes.Buffer

	count, err := Merge(reader, &buf, MergeOptions{})
	if err != nil {
		t.Fatalf("Merge failed: %v", err)
	}

	if count != 3 {
		t.Errorf("Expected 3 intervals, got %d", count)
	}

	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected:\n%s\nGot:\n%s", expected, buf.String())
	}
}

func TestMergeEmpty(t *testing.T) {
	input := ``

	reader := strings.NewReader(input)
	var buf bytes.Buffer

	count, err := Merge(reader, &buf, MergeOptions{})
	if err != nil {
		t.Fatalf("Merge failed: %v", err)
	}

	if count != 0 {
		t.Errorf("Expected 0 intervals, got %d", count)
	}
}

func TestMergeWithStats(t *testing.T) {
	input := `chr1	100	200
chr1	150	250
chr1	300	400
chr1	350	450`

	reader := strings.NewReader(input)
	var buf bytes.Buffer

	stats, err := MergeWithStats(reader, &buf, MergeOptions{})
	if err != nil {
		t.Fatalf("MergeWithStats failed: %v", err)
	}

	if stats.InputIntervals != 4 {
		t.Errorf("Expected 4 input intervals, got %d", stats.InputIntervals)
	}

	if stats.OutputIntervals != 2 {
		t.Errorf("Expected 2 output intervals, got %d", stats.OutputIntervals)
	}

	if stats.MergedCount != 2 {
		t.Errorf("Expected 2 merged intervals, got %d", stats.MergedCount)
	}
}

func TestMergeUnsorted(t *testing.T) {
	// Input is not sorted - should still work
	input := `chr1	300	400
chr1	100	200
chr1	150	250`

	expected := `chr1	100	250
chr1	300	400
`

	reader := strings.NewReader(input)
	var buf bytes.Buffer

	count, err := Merge(reader, &buf, MergeOptions{})
	if err != nil {
		t.Fatalf("Merge failed: %v", err)
	}

	if count != 2 {
		t.Errorf("Expected 2 merged intervals, got %d", count)
	}

	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected:\n%s\nGot:\n%s", expected, buf.String())
	}
}

func TestMergeWithCount(t *testing.T) {
	input := `chr1	100	200
chr1	150	250
chr1	300	400`

	expected := `chr1	100	250	2
chr1	300	400	1
`

	reader := strings.NewReader(input)
	var buf bytes.Buffer

	opts := MergeOptions{
		OutputFields: OutputFields{
			Count: true,
		},
	}

	count, err := Merge(reader, &buf, opts)
	if err != nil {
		t.Fatalf("Merge failed: %v", err)
	}

	if count != 2 {
		t.Errorf("Expected 2 merged intervals, got %d", count)
	}

	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected:\n%s\nGot:\n%s", expected, buf.String())
	}
}

func TestMergeBedGraph(t *testing.T) {
	input := `chr1	100	200	10
chr1	150	250	20
chr1	300	400	30`

	expected := `chr1	100	250	10
chr1	300	400	30
`

	reader := strings.NewReader(input)
	var buf bytes.Buffer

	opts := MergeOptions{
		OutputFields: OutputFields{
			BedGraph: true,
		},
	}

	count, err := Merge(reader, &buf, opts)
	if err != nil {
		t.Fatalf("Merge failed: %v", err)
	}

	if count != 2 {
		t.Errorf("Expected 2 merged intervals, got %d", count)
	}

	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected:\n%s\nGot:\n%s", expected, buf.String())
	}
}

func TestMergeStreaming(t *testing.T) {
	input := `chr1	100	200
chr1	150	250
chr2	100	200
chr2	150	250
chr3	100	200`

	expected := `chr1	100	250
chr2	100	250
chr3	100	200
`

	reader := strings.NewReader(input)
	var buf bytes.Buffer

	opts := MergeOptions{
		Streaming: true,
	}

	count, err := Merge(reader, &buf, opts)
	if err != nil {
		t.Fatalf("Merge failed: %v", err)
	}

	if count != 3 {
		t.Errorf("Expected 3 intervals, got %d", count)
	}

	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected:\n%s\nGot:\n%s", expected, buf.String())
	}
}

func TestMergeStreamingWithCount(t *testing.T) {
	input := `chr1	100	200
chr1	150	250
chr2	100	200
chr2	150	250`

	expected := `chr1	100	250	2
chr2	100	250	2
`

	reader := strings.NewReader(input)
	var buf bytes.Buffer

	opts := MergeOptions{
		Streaming: true,
		OutputFields: OutputFields{
			Count: true,
		},
	}

	count, err := Merge(reader, &buf, opts)
	if err != nil {
		t.Fatalf("Merge failed: %v", err)
	}

	if count != 2 {
		t.Errorf("Expected 2 intervals, got %d", count)
	}

	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected:\n%s\nGot:\n%s", expected, buf.String())
	}
}

// columnOpsInput is a small set of intervals: the first three on chr1 merge
// into one group (10..30), the fourth on chr1 is separate (40..50), and chr2
// has one interval. Columns: 4 = name, 5 = score-ish number.
const columnOpsInput = `chr1	10	20	a	5
chr1	15	30	b	7
chr1	18	25	a	7
chr1	40	50	c	3
chr2	5	15	d	9`

func TestMergeColumnOps(t *testing.T) {
	tests := []struct {
		name     string
		cols     string
		ops      string
		input    string
		expected string
	}{
		{
			name:  "sum on numeric column",
			cols:  "5",
			ops:   "sum",
			input: columnOpsInput,
			expected: "chr1\t10\t30\t19\n" +
				"chr1\t40\t50\t3\n" +
				"chr2\t5\t15\t9\n",
		},
		{
			name:  "min and max",
			cols:  "5,5",
			ops:   "min,max",
			input: columnOpsInput,
			expected: "chr1\t10\t30\t5\t7\n" +
				"chr1\t40\t50\t3\t3\n" +
				"chr2\t5\t15\t9\t9\n",
		},
		{
			// Upstream formats with 10 significant digits (KeyListOps default
			// precision), so 19/3 prints as 6.333333333, not full float64 width.
			name:  "mean produces float",
			cols:  "5",
			ops:   "mean",
			input: columnOpsInput,
			expected: "chr1\t10\t30\t6.333333333\n" +
				"chr1\t40\t50\t3\n" +
				"chr2\t5\t15\t9\n",
		},
		{
			name:  "median odd count",
			cols:  "5",
			ops:   "median",
			input: columnOpsInput,
			expected: "chr1\t10\t30\t7\n" +
				"chr1\t40\t50\t3\n" +
				"chr2\t5\t15\t9\n",
		},
		{
			name:  "count counts merged intervals",
			cols:  "4",
			ops:   "count",
			input: columnOpsInput,
			expected: "chr1\t10\t30\t3\n" +
				"chr1\t40\t50\t1\n" +
				"chr2\t5\t15\t1\n",
		},
		{
			name:  "count_distinct counts distinct values",
			cols:  "4",
			ops:   "count_distinct",
			input: columnOpsInput,
			expected: "chr1\t10\t30\t2\n" +
				"chr1\t40\t50\t1\n" +
				"chr2\t5\t15\t1\n",
		},
		{
			name:  "distinct lists unique values in sorted order",
			cols:  "4",
			ops:   "distinct",
			input: columnOpsInput,
			expected: "chr1\t10\t30\ta,b\n" +
				"chr1\t40\t50\tc\n" +
				"chr2\t5\t15\td\n",
		},
		{
			name:  "collapse keeps all values with dups",
			cols:  "4",
			ops:   "collapse",
			input: columnOpsInput,
			expected: "chr1\t10\t30\ta,b,a\n" +
				"chr1\t40\t50\tc\n" +
				"chr2\t5\t15\td\n",
		},
		{
			name:  "first and last",
			cols:  "4,4",
			ops:   "first,last",
			input: columnOpsInput,
			expected: "chr1\t10\t30\ta\ta\n" +
				"chr1\t40\t50\tc\tc\n" +
				"chr2\t5\t15\td\td\n",
		},
		{
			// chr1 group col5 = [5,7,7]; most frequent is 7.
			name:  "mode picks most frequent value",
			cols:  "5",
			ops:   "mode",
			input: columnOpsInput,
			expected: "chr1\t10\t30\t7\n" +
				"chr1\t40\t50\t3\n" +
				"chr2\t5\t15\t9\n",
		},
		{
			// chr1 group col5 = [5,7,7]; least frequent is 5 (tie-break: 5 seen first).
			name:  "antimode picks least frequent value, tie -> first seen",
			cols:  "5",
			ops:   "antimode",
			input: columnOpsInput,
			expected: "chr1\t10\t30\t5\n" +
				"chr1\t40\t50\t3\n" +
				"chr2\t5\t15\t9\n",
		},
		{
			// All distinct values are equally (least) frequent; antimode tie-breaks
			// to the first-seen value.
			name: "antimode tie-break across all-unique values",
			cols: "5",
			ops:  "antimode",
			input: "chr1\t10\t20\ta\t8\n" +
				"chr1\t12\t22\tb\t3\n" +
				"chr1\t14\t24\tc\t9\n",
			expected: "chr1\t10\t24\t8\n",
		},
		{
			name:  "single op applied to multiple columns",
			cols:  "4,5",
			ops:   "distinct",
			input: columnOpsInput,
			expected: "chr1\t10\t30\ta,b\t5,7\n" +
				"chr1\t40\t50\tc\t3\n" +
				"chr2\t5\t15\td\t9\n",
		},
		{
			name: "distinct then sum (task smoke example)",
			cols: "4,5",
			ops:  "distinct,sum",
			input: "chr1\t10\t20\ta\t5\n" +
				"chr1\t15\t30\tb\t7\n" +
				"chr1\t40\t50\tc\t3\n",
			expected: "chr1\t10\t30\ta,b\t12\n" +
				"chr1\t40\t50\tc\t3\n",
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			co, err := ParseColumnOps(tt.cols, tt.ops)
			if err != nil {
				t.Fatalf("ParseColumnOps(%q, %q) failed: %v", tt.cols, tt.ops, err)
			}
			var buf bytes.Buffer
			if _, err := Merge(strings.NewReader(tt.input), &buf, MergeOptions{ColumnOps: co}); err != nil {
				t.Fatalf("Merge failed: %v", err)
			}
			if buf.String() != tt.expected {
				t.Errorf("Output mismatch.\nExpected:\n%q\nGot:\n%q", tt.expected, buf.String())
			}
		})
	}
}

func TestParseColumnOpsErrors(t *testing.T) {
	tests := []struct {
		name string
		cols string
		ops  string
	}{
		{"columns without operations", "4,5", ""},
		{"operations without columns", "", "sum"},
		{"mismatched lengths", "4,5,6", "sum,min"},
		{"unsupported operation", "4", "bogus"},
		{"non-numeric column number", "x", "sum"},
		{"zero column number", "0", "sum"},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if _, err := ParseColumnOps(tt.cols, tt.ops); err == nil {
				t.Errorf("ParseColumnOps(%q, %q) expected error, got nil", tt.cols, tt.ops)
			}
		})
	}
}

func TestParseColumnOpsNeitherGiven(t *testing.T) {
	co, err := ParseColumnOps("", "")
	if err != nil {
		t.Fatalf("ParseColumnOps(\"\", \"\") returned error: %v", err)
	}
	if co != nil {
		t.Errorf("ParseColumnOps(\"\", \"\") expected nil ColumnOps, got %+v", co)
	}
}

// TestMergeColumnOpsNonNumericWarns confirms a non-numeric value under a numeric
// op produces the null value "." and an upstream-formatted warning (parity with
// bedtools merge.t23a/t23b) rather than aborting.
func TestMergeColumnOpsNonNumericWarns(t *testing.T) {
	input := "chr1\t10\t20\ta\tx\nchr1\t15\t30\tb\t7\n"
	co, err := ParseColumnOps("5", "sum")
	if err != nil {
		t.Fatalf("ParseColumnOps failed: %v", err)
	}
	var buf, warn bytes.Buffer
	if _, err := Merge(strings.NewReader(input), &buf, MergeOptions{ColumnOps: co, Warn: &warn}); err != nil {
		t.Fatalf("Merge errored on non-numeric value: %v", err)
	}
	if got := buf.String(); got != "chr1\t10\t30\t.\n" {
		t.Errorf("expected null-value output, got %q", got)
	}
	if !strings.Contains(warn.String(), "Non numeric value x in 5") {
		t.Errorf("warning should name the offending value+column, got: %q", warn.String())
	}
}

func TestMergeBED3UnchangedWithoutColumnOps(t *testing.T) {
	input := "chr1\t10\t20\ta\t5\nchr1\t15\t30\tb\t7\nchr1\t40\t50\tc\t3\n"
	expected := "chr1\t10\t30\nchr1\t40\t50\n"
	var buf bytes.Buffer
	count, err := Merge(strings.NewReader(input), &buf, MergeOptions{})
	if err != nil {
		t.Fatalf("Merge failed: %v", err)
	}
	if count != 2 {
		t.Errorf("expected 2 intervals, got %d", count)
	}
	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected:\n%q\nGot:\n%q", expected, buf.String())
	}
}

func TestMergeOutputFields(t *testing.T) {
	input := `chr1	100	200	name1	100	+
chr1	150	250	name2	200	+
chr1	300	400	name3	300	-`

	// Test with name and strand
	expected := `chr1	100	250	name1	100	+
chr1	300	400	name3	300	-
`

	reader := strings.NewReader(input)
	var buf bytes.Buffer

	opts := MergeOptions{
		OutputFields: OutputFields{
			Name:   true,
			Score:  true,
			Strand: true,
		},
	}

	count, err := Merge(reader, &buf, opts)
	if err != nil {
		t.Fatalf("Merge failed: %v", err)
	}

	if count != 2 {
		t.Errorf("Expected 2 merged intervals, got %d", count)
	}

	if buf.String() != expected {
		t.Errorf("Output mismatch.\nExpected:\n%s\nGot:\n%s", expected, buf.String())
	}
}
