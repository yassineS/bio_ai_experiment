package bedcluster

import (
	"bytes"
	"strings"
	"testing"
)

func TestCluster_Basic(t *testing.T) {
	// 3 records on chr1: 0-10, 5-15 overlap (cluster 1); 100-200 (cluster 2)
	in := strings.NewReader("chr1\t0\t10\nchr1\t5\t15\nchr1\t100\t200\n")
	var out bytes.Buffer
	n, err := Cluster(in, &out, Options{})
	if err != nil {
		t.Fatalf("Cluster: %v", err)
	}
	if n != 3 {
		t.Errorf("n = %d, want 3", n)
	}
	got := out.String()
	want := "chr1\t0\t10\t1\nchr1\t5\t15\t1\nchr1\t100\t200\t2\n"
	if got != want {
		t.Errorf("got:\n%s\nwant:\n%s", got, want)
	}
}

func TestCluster_GapDistance(t *testing.T) {
	// gap=5: book-end + 5bp gap should still cluster with -d 5
	in := strings.NewReader("chr1\t0\t10\nchr1\t15\t20\n")
	var out bytes.Buffer
	if _, err := Cluster(in, &out, Options{MaxDistance: 5}); err != nil {
		t.Fatalf("Cluster: %v", err)
	}
	want := "chr1\t0\t10\t1\nchr1\t15\t20\t1\n"
	if out.String() != want {
		t.Errorf("got:\n%s\nwant:\n%s", out.String(), want)
	}

	// gap=4: 5bp gap > 4 means new cluster
	in2 := strings.NewReader("chr1\t0\t10\nchr1\t15\t20\n")
	var out2 bytes.Buffer
	if _, err := Cluster(in2, &out2, Options{MaxDistance: 4}); err != nil {
		t.Fatalf("Cluster: %v", err)
	}
	want2 := "chr1\t0\t10\t1\nchr1\t15\t20\t2\n"
	if out2.String() != want2 {
		t.Errorf("got:\n%s\nwant:\n%s", out2.String(), want2)
	}
}

func TestCluster_BookEndIsSameCluster(t *testing.T) {
	// Default MaxDistance=0: an interval ending at 10 and one starting at 10
	// have a gap of 0, which is NOT > 0, so they cluster (matches upstream).
	in := strings.NewReader("chr1\t0\t10\nchr1\t10\t20\n")
	var out bytes.Buffer
	if _, err := Cluster(in, &out, Options{}); err != nil {
		t.Fatalf("Cluster: %v", err)
	}
	if !strings.Contains(out.String(), "chr1\t10\t20\t1\n") {
		t.Errorf("book-ended should cluster: %q", out.String())
	}
}

func TestCluster_StrandSpec(t *testing.T) {
	// Two overlapping intervals on different strands should be in different
	// clusters when -s is set; outputs grouped by strand alphabetically (- < +).
	in := strings.NewReader(
		"chr1\t0\t10\ta\t1\t+\n" +
			"chr1\t5\t15\tb\t2\t-\n",
	)
	var out bytes.Buffer
	if _, err := Cluster(in, &out, Options{StrandSpec: true}); err != nil {
		t.Fatalf("Cluster: %v", err)
	}
	got := out.String()
	// '+' (ASCII 43) sorts before '-' (ASCII 45), matching upstream's
	// per-strand grouping order in test cluster.t2.
	if !strings.HasPrefix(got, "chr1\t0\t10\ta\t1\t+\t1\n") {
		t.Errorf("expected '+' strand row first with cid=1, got:\n%s", got)
	}
	if !strings.Contains(got, "chr1\t5\t15\tb\t2\t-\t2\n") {
		t.Errorf("expected '-' strand row with cid=2, got:\n%s", got)
	}
}

func TestCluster_DifferentChromsRestartIDs(t *testing.T) {
	// Cluster IDs DO NOT restart across chroms — they are global, monotonic.
	in := strings.NewReader("chr1\t0\t10\nchr2\t0\t10\n")
	var out bytes.Buffer
	if _, err := Cluster(in, &out, Options{}); err != nil {
		t.Fatalf("Cluster: %v", err)
	}
	want := "chr1\t0\t10\t1\nchr2\t0\t10\t2\n"
	if out.String() != want {
		t.Errorf("got:\n%s\nwant:\n%s", out.String(), want)
	}
}

func TestCluster_SortsUnsortedInput(t *testing.T) {
	// Input is unsorted; cluster should sort internally.
	in := strings.NewReader("chr1\t100\t200\nchr1\t0\t10\nchr1\t5\t15\n")
	var out bytes.Buffer
	if _, err := Cluster(in, &out, Options{}); err != nil {
		t.Fatalf("Cluster: %v", err)
	}
	want := "chr1\t0\t10\t1\nchr1\t5\t15\t1\nchr1\t100\t200\t2\n"
	if out.String() != want {
		t.Errorf("got:\n%s\nwant:\n%s", out.String(), want)
	}
}

func TestCluster_EmptyInput(t *testing.T) {
	var out bytes.Buffer
	n, err := Cluster(strings.NewReader(""), &out, Options{})
	if err != nil {
		t.Fatalf("Cluster: %v", err)
	}
	if n != 0 || out.Len() != 0 {
		t.Errorf("n=%d outLen=%d, want 0/0", n, out.Len())
	}
}

func TestCluster_BadInput(t *testing.T) {
	var out bytes.Buffer
	if _, err := Cluster(strings.NewReader("chr1\toops\t10\n"), &out, Options{}); err == nil {
		t.Errorf("expected error on bad start")
	}
}

func TestCluster_PreservesAllColumns(t *testing.T) {
	// 6 columns in -> 7 columns out (cluster ID appended).
	in := strings.NewReader("chr1\t0\t10\tname\t100\t+\n")
	var out bytes.Buffer
	if _, err := Cluster(in, &out, Options{}); err != nil {
		t.Fatalf("Cluster: %v", err)
	}
	if out.String() != "chr1\t0\t10\tname\t100\t+\t1\n" {
		t.Errorf("expected all columns preserved, got: %q", out.String())
	}
}
