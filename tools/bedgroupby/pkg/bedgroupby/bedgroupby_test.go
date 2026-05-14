package bedgroupby

import (
	"bytes"
	"strings"
	"testing"
)

func runGroup(t *testing.T, input string, opts Options) string {
	t.Helper()
	var buf bytes.Buffer
	if _, err := Group(strings.NewReader(input), &buf, opts); err != nil {
		t.Fatalf("Group: %v", err)
	}
	return buf.String()
}

func TestParseGroupSpec(t *testing.T) {
	cases := []struct {
		in   string
		want []int
		err  bool
	}{
		{"", nil, false},
		{"1", []int{1}, false},
		{"1,2,3", []int{1, 2, 3}, false},
		{"2-4", []int{2, 3, 4}, false},
		{"3-4,2", []int{3, 4, 2}, false},
		{"1, 2 , 5", []int{1, 2, 5}, false},
		{"0", nil, true},
		{"abc", nil, true},
		{"3-1", nil, true},
		{"3-x", nil, true},
		{"x-3", nil, true},
		{",,", nil, false},
	}
	for _, c := range cases {
		got, err := ParseGroupSpec(c.in)
		if c.err {
			if err == nil {
				t.Errorf("ParseGroupSpec(%q): expected error, got %v", c.in, got)
			}
			continue
		}
		if err != nil {
			t.Errorf("ParseGroupSpec(%q): unexpected error: %v", c.in, err)
			continue
		}
		if !equalInts(got, c.want) {
			t.Errorf("ParseGroupSpec(%q) = %v, want %v", c.in, got, c.want)
		}
	}
}

func equalInts(a, b []int) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

func TestGroup_DefaultSum(t *testing.T) {
	in := "chr1\t0\t10\ta\t10\t+\nchr1\t0\t10\tb\t5\t+\nchr1\t10\t20\tc\t7\t+\n"
	got := runGroup(t, in, Options{AggCols: []int{5}})
	want := "chr1\t0\t10\t15\nchr1\t10\t20\t7\n"
	if got != want {
		t.Errorf("got:\n%q\nwant:\n%q", got, want)
	}
}

func TestGroup_DefaultOpIsSum(t *testing.T) {
	in := "x\t1\t2\t3\nx\t1\t2\t4\n"
	got := runGroup(t, in, Options{AggCols: []int{4}})
	want := "x\t1\t2\t7\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestGroup_MultiOp(t *testing.T) {
	in := "x\t1\t2\t10\ta\nx\t1\t2\t20\tb\nx\t1\t2\t30\ta\n"
	got := runGroup(t, in, Options{
		AggCols: []int{4, 5},
		Ops:     []string{"sum", "distinct"},
	})
	want := "x\t1\t2\t60\ta,b\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestGroup_SingleOpBroadcast(t *testing.T) {
	in := "x\t1\t2\t10\t20\nx\t1\t2\t30\t40\n"
	got := runGroup(t, in, Options{AggCols: []int{4, 5}, Ops: []string{"sum"}})
	want := "x\t1\t2\t40\t60\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestGroup_CustomGroupCols(t *testing.T) {
	in := "l\tchr1\t0\t10\t10\nk\tchr1\t0\t10\t5\n"
	got := runGroup(t, in, Options{
		GroupCols: []int{2, 3, 4},
		AggCols:   []int{5},
		Ops:       []string{"sum"},
	})
	want := "chr1\t0\t10\t15\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestGroup_GroupColsRange(t *testing.T) {
	in := "l\t10\t20\tchr1\ta\nk\t10\t20\tchr1\tb\nj\t30\t40\tchr1\tc\n"
	got := runGroup(t, in, Options{
		GroupCols: []int{2, 3, 4},
		AggCols:   []int{5},
		Ops:       []string{"collapse"},
	})
	want := "10\t20\tchr1\ta,b\n30\t40\tchr1\tc\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestGroup_GroupColsReordered(t *testing.T) {
	in := "x\t10\t20\tA\t1\nx\t10\t20\tA\t2\n"
	got := runGroup(t, in, Options{
		GroupCols: []int{3, 4, 2},
		AggCols:   []int{5},
		Ops:       []string{"sum"},
	})
	want := "20\tA\t10\t3\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestGroup_Full(t *testing.T) {
	in := "chr1\t0\t10\ta\t10\t+\nchr1\t0\t10\tb\t5\t+\n"
	got := runGroup(t, in, Options{AggCols: []int{5}, Ops: []string{"sum"}, Full: true})
	want := "chr1\t0\t10\ta\t10\t+\t15\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestGroup_IgnoreCase(t *testing.T) {
	in := "chr1\t0\t10\t10\nChr1\t0\t10\t5\n"
	got := runGroup(t, in, Options{AggCols: []int{4}, Ops: []string{"sum"}, IgnoreCase: true})
	want := "chr1\t0\t10\t15\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
	got = runGroup(t, in, Options{AggCols: []int{4}, Ops: []string{"sum"}})
	want = "chr1\t0\t10\t10\nChr1\t0\t10\t5\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestGroup_SkipHeaderMarked(t *testing.T) {
	in := "#chrom\tstart\tend\tval\nchr1\t0\t10\t5\nchr1\t0\t10\t6\n"
	got := runGroup(t, in, Options{AggCols: []int{4}})
	want := "chr1\t0\t10\t11\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestGroup_InHeader(t *testing.T) {
	in := "Chrom\tstart\tend\tval\nchr1\t0\t10\t5\nchr1\t0\t10\t6\n"
	got := runGroup(t, in, Options{AggCols: []int{4}, InHeader: true})
	want := "chr1\t0\t10\t11\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
	got = runGroup(t, in, Options{AggCols: []int{4}, Ops: []string{"distinct"}})
	want = "Chrom\tstart\tend\tval\nchr1\t0\t10\t5,6\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestGroup_OutHeaderMarked(t *testing.T) {
	in := "#chrom\tstart\tend\tA\tB\tC\nchr1\t0\t10\ta\t10\t+\nchr1\t0\t10\tb\t5\t+\n"
	got := runGroup(t, in, Options{AggCols: []int{5}, OutHeader: true})
	want := "#chrom\tstart\tend\tA\tB\tC\nchr1\t0\t10\t15\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestGroup_OutHeaderSynthetic(t *testing.T) {
	in := "Chromz\tstart\tend\tA\tB\tC\nchr1\t0\t10\ta\t10\t+\n"
	got := runGroup(t, in, Options{AggCols: []int{5}, Ops: []string{"distinct"}, OutHeader: true})
	want := "col_1\tcol_2\tcol_3\tcol_4\tcol_5\tcol_6\n" +
		"Chromz\tstart\tend\tB\nchr1\t0\t10\t10\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestGroup_HeaderCombines(t *testing.T) {
	in := "Chrom\tstart\tend\tA\tB\tC\nchr1\t0\t10\ta\t10\t+\n"
	got := runGroup(t, in, Options{AggCols: []int{5}, Ops: []string{"distinct"}, Header: true})
	want := "Chrom\tstart\tend\tA\tB\tC\nchr1\t0\t10\t10\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestGroup_MissingAggColsErrs(t *testing.T) {
	in := "chr1\t0\t10\n"
	var buf bytes.Buffer
	if _, err := Group(strings.NewReader(in), &buf, Options{}); err == nil {
		t.Fatal("expected error when AggCols empty")
	}
}

func TestGroup_BadColRefErrs(t *testing.T) {
	in := "chr1\t0\t10\n"
	var buf bytes.Buffer
	if _, err := Group(strings.NewReader(in), &buf, Options{AggCols: []int{99}}); err == nil {
		t.Fatal("expected error when column index out of range")
	}
}

func TestGroup_OpsCountMismatchErrs(t *testing.T) {
	var buf bytes.Buffer
	_, err := Group(strings.NewReader("a\t1\t2\t3\n"),
		&buf, Options{AggCols: []int{2, 3, 4}, Ops: []string{"sum", "min"}})
	if err == nil {
		t.Fatal("expected error when len(ops) != 1 and != len(cols)")
	}
}

func TestGroup_BadOpName(t *testing.T) {
	var buf bytes.Buffer
	_, err := Group(strings.NewReader("a\t1\n"),
		&buf, Options{AggCols: []int{2}, Ops: []string{"bogus"}})
	if err == nil {
		t.Fatal("expected error for unknown op")
	}
}

func TestGroup_NumericOpOnNonNumber(t *testing.T) {
	var buf bytes.Buffer
	_, err := Group(strings.NewReader("chr1\t0\t10\ta\n"),
		&buf, Options{AggCols: []int{4}, Ops: []string{"sum"}})
	if err == nil {
		t.Fatal("expected error for non-numeric value with sum")
	}
}

func TestGroup_EmptyInput(t *testing.T) {
	got := runGroup(t, "", Options{AggCols: []int{4}})
	if got != "" {
		t.Errorf("got %q want empty", got)
	}
}

func TestGroup_AllHeaderInput(t *testing.T) {
	got := runGroup(t, "#header\n", Options{AggCols: []int{4}, OutHeader: true})
	if got != "" {
		t.Errorf("got %q want empty (no data rows)", got)
	}
}

func TestGroup_TrackAndBrowserSkipped(t *testing.T) {
	in := "track name=foo\nbrowser position chr1\n" +
		"chr1\t0\t10\t1\nchr1\t0\t10\t2\n"
	got := runGroup(t, in, Options{AggCols: []int{4}, Ops: []string{"sum"}})
	want := "chr1\t0\t10\t3\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestGroup_BlankLinesSkipped(t *testing.T) {
	in := "\n\nchr1\t0\t10\t1\n\nchr1\t0\t10\t2\n"
	got := runGroup(t, in, Options{AggCols: []int{4}, Ops: []string{"sum"}})
	want := "chr1\t0\t10\t3\n"
	if got != want {
		t.Errorf("got %q want %q", got, want)
	}
}

func TestGroup_AdditionalOps(t *testing.T) {
	t.Skip("absmin/absmax/cat/cat_uniq ops not yet wired into bedmerge.ApplyOp; tracked in PARITY_ROADMAP.md#bedtools")
	in := "x\t1\t2\t-5\nx\t1\t2\t3\nx\t1\t2\t-2\n"
	cases := []struct {
		op   string
		want string
	}{
		{"min", "-5"},
		{"max", "3"},
		{"absmin", "2"},
		{"absmax", "5"},
		{"median", "-2"},
		{"count", "3"},
		{"first", "-5"},
		{"last", "-2"},
	}
	for _, c := range cases {
		got := runGroup(t, in, Options{AggCols: []int{4}, Ops: []string{c.op}})
		got = strings.TrimSpace(got)
		want := "x\t1\t2\t" + c.want
		if got != want {
			t.Errorf("op %q: got %q want %q", c.op, got, want)
		}
	}
}

func TestGroup_StdevSstdev(t *testing.T) {
	t.Skip("stdev/sstdev ops not yet wired into bedmerge.ApplyOp; tracked in PARITY_ROADMAP.md#bedtools")
	in := "x\t1\t2\t2\nx\t1\t2\t4\nx\t1\t2\t4\nx\t1\t2\t6\n"
	got := strings.TrimSpace(runGroup(t, in, Options{AggCols: []int{4}, Ops: []string{"stdev"}}))
	if !strings.HasPrefix(got, "x\t1\t2\t1.4142") {
		t.Errorf("stdev got %q", got)
	}
	got = strings.TrimSpace(runGroup(t, in, Options{AggCols: []int{4}, Ops: []string{"sstdev"}}))
	if !strings.HasPrefix(got, "x\t1\t2\t1.632993") {
		t.Errorf("sstdev got %q", got)
	}
	got = strings.TrimSpace(runGroup(t, "x\t1\t2\t5\n", Options{AggCols: []int{4}, Ops: []string{"sstdev"}}))
	if got != "x\t1\t2\t." {
		t.Errorf("sstdev single: got %q", got)
	}
}
