package bedigv

import (
	"bytes"
	"errors"
	"io"
	"strings"
	"testing"
)

func TestRun_Default(t *testing.T) {
	in := "chr1\t100\t200\nchr2\t1000\t1100\n"
	want := "snapshotDirectory ./\n" +
		"goto chr1:100-200\n" +
		"snapshot chr1_100_200.png\n" +
		"goto chr2:1000-1100\n" +
		"snapshot chr2_1000_1100.png\n"
	var buf bytes.Buffer
	n, err := Run(strings.NewReader(in), &buf, Options{})
	if err != nil {
		t.Fatalf("Run: %v", err)
	}
	if n != 2 {
		t.Errorf("snapshots = %d, want 2", n)
	}
	if got := buf.String(); got != want {
		t.Errorf("output mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

func TestRun_PathSessionSortCollapseSlopName(t *testing.T) {
	in := "chr1\t100\t200\tregionA\nchrX\t50\t60\tregionB\n"
	opts := Options{
		Path:      "/tmp/shots/",
		Session:   "/data/sess.xml",
		Sort:      SortPosition,
		Collapse:  true,
		UseNames:  true,
		Slop:      10,
		ImageType: ImageSVG,
	}
	want := "snapshotDirectory /tmp/shots/\n" +
		"load /data/sess.xml\n" +
		"goto chr1:90-210\n" +
		"sort position\n" +
		"collapse\n" +
		"snapshot chr1_100_200_regionA_slop10.svg\n" +
		"goto chrX:40-70\n" +
		"sort position\n" +
		"collapse\n" +
		"snapshot chrX_50_60_regionB_slop10.svg\n"
	var buf bytes.Buffer
	if _, err := Run(strings.NewReader(in), &buf, opts); err != nil {
		t.Fatalf("Run: %v", err)
	}
	if got := buf.String(); got != want {
		t.Errorf("output mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

func TestRun_NegativeSlopUnderflowMatchesUpstream(t *testing.T) {
	// Upstream emits the raw arithmetic result, including negative
	// start-slop. We preserve that to keep byte-for-byte parity.
	in := "chr1\t5\t10\n"
	opts := Options{Slop: 100}
	var buf bytes.Buffer
	if _, err := Run(strings.NewReader(in), &buf, opts); err != nil {
		t.Fatalf("Run: %v", err)
	}
	got := buf.String()
	if !strings.Contains(got, "goto chr1:-95-110\n") {
		t.Errorf("expected raw negative locus; got: %s", got)
	}
	if !strings.Contains(got, "snapshot chr1_5_10_slop100.png\n") {
		t.Errorf("expected slop suffix in filename; got: %s", got)
	}
}

func TestRun_UseNamesErrorsOnEmptyName(t *testing.T) {
	// Record has only 3 columns → Name is empty → UseNames must error.
	in := "chr1\t100\t200\n"
	var buf bytes.Buffer
	_, err := Run(strings.NewReader(in), &buf, Options{UseNames: true})
	if err == nil {
		t.Fatalf("expected error for empty name with UseNames")
	}
}

func TestRun_RejectsBadSlop(t *testing.T) {
	var buf bytes.Buffer
	_, err := Run(strings.NewReader("chr1\t1\t2\n"), &buf, Options{Slop: -1})
	if err == nil {
		t.Fatalf("expected error for negative slop")
	}
}

func TestRun_RejectsBadSort(t *testing.T) {
	var buf bytes.Buffer
	_, err := Run(strings.NewReader("chr1\t1\t2\n"), &buf, Options{Sort: SortType("bogus")})
	if err == nil {
		t.Fatalf("expected error for invalid sort type")
	}
}

func TestRun_EmptyInputStillEmitsHeader(t *testing.T) {
	var buf bytes.Buffer
	n, err := Run(strings.NewReader(""), &buf, Options{Path: "snaps/"})
	if err != nil {
		t.Fatalf("Run: %v", err)
	}
	if n != 0 {
		t.Errorf("snapshots = %d, want 0", n)
	}
	if got := buf.String(); got != "snapshotDirectory snaps/\n" {
		t.Errorf("output mismatch.\nwant header only, got: %s", got)
	}
}

func TestRun_JPGImageType(t *testing.T) {
	in := "chrM\t0\t100\n"
	var buf bytes.Buffer
	if _, err := Run(strings.NewReader(in), &buf, Options{ImageType: ImageJPG}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	if !strings.HasSuffix(strings.TrimRight(buf.String(), "\n"), "snapshot chrM_0_100.jpg") {
		t.Errorf("expected .jpg snapshot, got: %s", buf.String())
	}
}

func TestRun_PropagatesBedParseError(t *testing.T) {
	// chromStart is not an integer → bed.Reader returns an error.
	in := "chr1\tNOT_A_NUMBER\t200\n"
	var buf bytes.Buffer
	_, err := Run(strings.NewReader(in), &buf, Options{})
	if err == nil {
		t.Fatalf("expected parse error from bed.Reader")
	}
}

func TestIsValidSort(t *testing.T) {
	good := []SortType{SortNone, SortBase, SortPosition, SortStrand, SortQuality, SortSample, SortReadGroup}
	for _, s := range good {
		if !IsValidSort(s) {
			t.Errorf("IsValidSort(%q) = false, want true", s)
		}
	}
	for _, s := range []SortType{"random", "BASE", "foo"} {
		if IsValidSort(s) {
			t.Errorf("IsValidSort(%q) = true, want false", s)
		}
	}
}

// errWriter returns an error on Write, exercising the write-error paths.
type errWriter struct{ failAfter int }

func (e *errWriter) Write(p []byte) (int, error) {
	if e.failAfter <= 0 {
		return 0, errors.New("forced write failure")
	}
	e.failAfter--
	return len(p), nil
}

func TestRun_WriteErrorPropagates(t *testing.T) {
	// bufio.Writer buffers, so we need enough records to flush. Use a
	// very small failAfter so even the buffered flush triggers the
	// error path.
	in := strings.Repeat("chr1\t100\t200\n", 2000)
	w := &errWriter{failAfter: 0}
	_, err := Run(strings.NewReader(in), w, Options{})
	if err == nil {
		t.Fatalf("expected write error to propagate")
	}
	if errors.Is(err, io.EOF) {
		t.Errorf("io.EOF should not surface as the public error")
	}
}
