package bedunionbedg

import (
	"bytes"
	"io"
	"strings"
	"testing"
)

func readers(srcs ...string) []io.Reader {
	rs := make([]io.Reader, len(srcs))
	for i, s := range srcs {
		rs[i] = strings.NewReader(s)
	}
	return rs
}

func runUnion(t *testing.T, opts Options, srcs ...string) string {
	t.Helper()
	var out bytes.Buffer
	if err := Union(readers(srcs...), &out, opts); err != nil {
		t.Fatalf("Union: %v", err)
	}
	return out.String()
}

func TestUnionBasic(t *testing.T) {
	one := "chr1\t1000\t1500\t10\nchr1\t2000\t2100\t20\n"
	two := "chr1\t900\t1600\t60\nchr1\t1700\t2050\t50\n"
	three := "chr1\t1980\t2070\t80\nchr1\t2090\t2100\t20\n"
	got := runUnion(t, Options{}, one, two, three)
	want := "chr1\t900\t1000\t0\t60\t0\n" +
		"chr1\t1000\t1500\t10\t60\t0\n" +
		"chr1\t1500\t1600\t0\t60\t0\n" +
		"chr1\t1700\t1980\t0\t50\t0\n" +
		"chr1\t1980\t2000\t0\t50\t80\n" +
		"chr1\t2000\t2050\t20\t50\t80\n" +
		"chr1\t2050\t2070\t20\t0\t80\n" +
		"chr1\t2070\t2090\t20\t0\t0\n" +
		"chr1\t2090\t2100\t20\t0\t20\n"
	if got != want {
		t.Fatalf("got:\n%s\nwant:\n%s", got, want)
	}
}

func TestUnionHeader(t *testing.T) {
	got := runUnion(t, Options{PrintHeader: true},
		"chr1\t0\t10\t1\n", "chr1\t5\t15\t2\n")
	if !strings.HasPrefix(got, "chrom\tstart\tend\n") {
		t.Fatalf("expected header, got:\n%s", got)
	}
}

func TestUnionHeaderNames(t *testing.T) {
	got := runUnion(t, Options{PrintHeader: true, Names: []string{"A", "B"}},
		"chr1\t0\t10\t1\n", "chr1\t5\t15\t2\n")
	if !strings.HasPrefix(got, "chrom\tstart\tend\tA\tB\n") {
		t.Fatalf("expected named header, got:\n%s", got)
	}
}

func TestUnionEmpty(t *testing.T) {
	got := runUnion(t, Options{PrintEmpty: true, Sizes: map[string]int64{"chr1": 30}},
		"chr1\t10\t20\t5\n", "chr1\t15\t25\t7\n")
	want := "chr1\t0\t10\t0\t0\n" +
		"chr1\t10\t15\t5\t0\n" +
		"chr1\t15\t20\t5\t7\n" +
		"chr1\t20\t25\t0\t7\n" +
		"chr1\t25\t30\t0\t0\n"
	if got != want {
		t.Fatalf("got:\n%s\nwant:\n%s", got, want)
	}
}

func TestUnionFiller(t *testing.T) {
	got := runUnion(t, Options{PrintEmpty: true, Sizes: map[string]int64{"chr1": 30}, Filler: "N/A"},
		"chr1\t10\t20\t5\n", "chr1\t15\t25\t7\n")
	if !strings.Contains(got, "N/A") || !strings.Contains(got, "chr1\t0\t10\tN/A\tN/A\n") {
		t.Fatalf("filler not applied, got:\n%s", got)
	}
}

func TestUnionFloatDepths(t *testing.T) {
	got := runUnion(t, Options{},
		"chr1\t0\t10\t1.5\nchr1\t10\t20\t2.75\n", "chr1\t5\t15\t3.0\n")
	want := "chr1\t0\t5\t1.5\t0\n" +
		"chr1\t5\t10\t1.5\t3.0\n" +
		"chr1\t10\t15\t2.75\t3.0\n" +
		"chr1\t15\t20\t2.75\t0\n"
	if got != want {
		t.Fatalf("got:\n%s\nwant:\n%s", got, want)
	}
}

func TestUnionMultiChromUnsorted(t *testing.T) {
	// Files not globally chrom-sorted: matches upstream's per-front-item chrom
	// selection and queue sweep.
	got := runUnion(t, Options{},
		"chr2\t10\t20\t5\nchr1\t10\t20\t7\n", "chr1\t15\t25\t3\n")
	want := "chr1\t15\t25\t0\t3\n" +
		"chr2\t10\t20\t5\t0\n" +
		"chr1\t10\t20\t7\t0\n"
	if got != want {
		t.Fatalf("got:\n%s\nwant:\n%s", got, want)
	}
}

func TestUnionTrackLineSkipped(t *testing.T) {
	got := runUnion(t, Options{},
		"track type=bedGraph\nchr1\t0\t10\t1\n", "chr1\t5\t15\t2\n")
	want := "chr1\t0\t5\t1\t0\n" +
		"chr1\t5\t10\t1\t2\n" +
		"chr1\t10\t15\t0\t2\n"
	if got != want {
		t.Fatalf("got:\n%s\nwant:\n%s", got, want)
	}
}

func TestUnionTooFewFiles(t *testing.T) {
	var out bytes.Buffer
	if err := Union(readers("chr1\t0\t10\t1\n"), &out, Options{}); err == nil {
		t.Fatalf("expected error for single input")
	}
}

func TestReadChromSizes(t *testing.T) {
	cs, err := ReadChromSizes(strings.NewReader("chr1\t100\nchr2\t200\n"))
	if err != nil {
		t.Fatalf("ReadChromSizes: %v", err)
	}
	if cs["chr1"] != 100 || cs["chr2"] != 200 {
		t.Fatalf("unexpected sizes: %+v", cs)
	}
}
