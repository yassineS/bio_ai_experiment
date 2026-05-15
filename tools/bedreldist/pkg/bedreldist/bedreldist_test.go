package bedreldist

import (
	"bytes"
	"strings"
	"testing"
)

// runReldist is a small helper that feeds string inputs through Run.
func runReldist(t *testing.T, a, b string, opts Options) (*Result, string) {
	t.Helper()
	var out bytes.Buffer
	res, err := Run(strings.NewReader(a), strings.NewReader(b), &out, opts)
	if err != nil {
		t.Fatalf("Run: %v", err)
	}
	return res, out.String()
}

// TestHandComputed_Symmetric:
// A: chr1 0-10 (mid=5), 90-100 (mid=95)
// B mids: chr1: [0, 50, 100] (records 0-1 -> mid 0, 49-51 -> mid 50, 99-101 -> mid 100).
//
//	For A1 (mid=5): lower_bound(5) -> idx=1 (50>=5), low_idx=0, high_idx=1,
//	  left=0, right=50, leftDist=5, rightDist=45, min=5, rel=5/50=0.10 -> bin 10.
//	For A2 (mid=95): lower_bound(95) -> idx=2 (100>=95), low_idx=1, high_idx=2,
//	  left=50, right=100, leftDist=45, rightDist=5, rel=5/50=0.10 -> bin 10.
//
// Total=2, bin 10: 2 -> "0.10\t2\t2\t1.000\n".
func TestHandComputed_Symmetric(t *testing.T) {
	a := "chr1\t0\t10\nchr1\t90\t100\n"
	b := "chr1\t0\t1\nchr1\t49\t51\nchr1\t99\t101\n"
	_, out := runReldist(t, a, b, Options{})
	want := "reldist\tcount\ttotal\tfraction\n0.10\t2\t2\t1.000\n"
	if out != want {
		t.Errorf("histogram mismatch.\nwant:\n%q\ngot:\n%q", want, out)
	}
}

// TestHandComputed_SelfIntersect: every A matches a B interval at the same
// position, so the minimum distance from each A's midpoint to a B-midpoint is
// 0, putting all queries in bin 0.00.
func TestHandComputed_SelfIntersect(t *testing.T) {
	a := "chr1\t10\t20\nchr1\t30\t40\nchr1\t50\t60\n"
	// At least 4 B intervals so that A2 and A3 (the middle queries) have a
	// midpoint that is between two distinct B midpoints, not at the rightmost
	// boundary (upstream skips that case).
	b := "chr1\t10\t20\nchr1\t30\t40\nchr1\t50\t60\nchr1\t70\t80\n"
	res, out := runReldist(t, a, b, Options{})
	want := "reldist\tcount\ttotal\tfraction\n0.00\t3\t3\t1.000\n"
	if out != want {
		t.Errorf("histogram mismatch.\nwant:\n%q\ngot:\n%q", want, out)
	}
	if res.Total != 3 {
		t.Errorf("total=%d want 3", res.Total)
	}
}

// TestHandComputed_Detail emits one row per A in detail mode and skips queries
// whose chromosome is not in B (chr2 here is absent from B).
func TestHandComputed_Detail(t *testing.T) {
	a := "chr1\t10\t30\tnameA\nchr2\t0\t10\tnameC\n" // chr2 absent in B
	b := "chr1\t0\t10\nchr1\t40\t60\nchr1\t80\t100\n"
	// A1 mid=20. B mids: [5,50,90]. lower_bound(20)=1 (50>=20), low_idx=0,
	// high_idx=1, left=5, right=50, leftDist=15, rightDist=30, rel=15/45 ~= 0.333.
	// chr2 query is dropped because chr2 is absent in B.
	res, out := runReldist(t, a, b, Options{Detail: true})
	if res.Total != 1 {
		t.Errorf("total=%d want 1", res.Total)
	}
	// Detail rows are followed by the histogram NOT being printed.
	wantLine := "chr1\t10\t30\tnameA\t0.333\n"
	if out != wantLine {
		t.Errorf("detail mismatch.\nwant:\n%q\ngot:\n%q", wantLine, out)
	}
}

// TestHandComputed_AtLastBin: A's midpoint is exactly between two B midpoints,
// so rel=0.5 and the query falls in bin 0.50.
func TestHandComputed_AtLastBin(t *testing.T) {
	a := "chr1\t1\t2\nchr1\t2\t3\n" // mid 1, 2
	b := "chr1\t1\t2\nchr1\t1\t2\nchr1\t3\t4\n"
	// Mirrors upstream issue_711 fixture: B midpoints [1,1,3].
	// A1 mid=1: low_idx=0, left=1, right=1, leftDist=0 -> bin 0.
	// A2 mid=2: lower_bound(2) -> 2 (3>=2), low_idx=1, left=1, right=3,
	//   leftDist=1, rightDist=1, rel=0.5 -> bin 50.
	_, out := runReldist(t, a, b, Options{})
	want := "reldist\tcount\ttotal\tfraction\n0.00\t1\t2\t0.500\n0.50\t1\t2\t0.500\n"
	if out != want {
		t.Errorf("histogram mismatch.\nwant:\n%q\ngot:\n%q", want, out)
	}
}

// TestEmptyB returns an empty histogram (total=0) without errors. Upstream
// prints just the header in this case.
func TestEmptyB(t *testing.T) {
	a := "chr1\t0\t10\n"
	_, out := runReldist(t, a, "", Options{})
	want := "reldist\tcount\ttotal\tfraction\n"
	if out != want {
		t.Errorf("histogram mismatch.\nwant:\n%q\ngot:\n%q", want, out)
	}
}

// TestDetail_BED6 exercises the BED4/BED5/BED6 branches of writeDetail.
func TestDetail_BED6(t *testing.T) {
	// A1 (BED6, strand +): mid=20. B mids: [5,50,90]. lower_bound(20)=1,
	// low_idx=0, left=5, right=50, leftDist=15, rightDist=30 -> rel=15/45=0.333.
	a := "chr1\t10\t30\tnameA\t100\t+\n"
	b := "chr1\t0\t10\nchr1\t40\t60\nchr1\t80\t100\n"
	_, out := runReldist(t, a, b, Options{Detail: true})
	want := "chr1\t10\t30\tnameA\t100\t+\t0.333\n"
	if out != want {
		t.Errorf("detail BED6 mismatch.\nwant:\n%q\ngot:\n%q", want, out)
	}
}

// TestDetail_BED5_NoStrand: name + score populated, no strand.
func TestDetail_BED5_NoStrand(t *testing.T) {
	a := "chr1\t10\t30\tnameA\t250\n"
	b := "chr1\t0\t10\nchr1\t40\t60\nchr1\t80\t100\n"
	_, out := runReldist(t, a, b, Options{Detail: true})
	want := "chr1\t10\t30\tnameA\t250\t0.333\n"
	if out != want {
		t.Errorf("detail BED5 mismatch.\nwant:\n%q\ngot:\n%q", want, out)
	}
}

// errWriter returns an error after `goodBytes` bytes have been written; used
// to exercise the error paths in writeHistogram / writeDetail.
type errWriter struct {
	good, written int
}

func (e *errWriter) Write(p []byte) (int, error) {
	remain := e.good - e.written
	if remain <= 0 {
		return 0, errBoom
	}
	if len(p) <= remain {
		e.written += len(p)
		return len(p), nil
	}
	e.written += remain
	return remain, errBoom
}

var errBoom = &boomErr{}

type boomErr struct{}

func (*boomErr) Error() string { return "boom" }

func TestRun_WriterErrorsBubbleUp(t *testing.T) {
	a := "chr1\t0\t10\n"
	b := "chr1\t0\t10\nchr1\t20\t30\n"
	w := &errWriter{good: 0}
	_, err := Run(strings.NewReader(a), strings.NewReader(b), w, Options{})
	if err == nil {
		t.Fatalf("expected write error, got nil")
	}
}

func TestRun_WriterErrorsBubbleUpDetail(t *testing.T) {
	a := "chr1\t0\t10\n"
	b := "chr1\t0\t10\nchr1\t20\t30\n"
	w := &errWriter{good: 0}
	_, err := Run(strings.NewReader(a), strings.NewReader(b), w, Options{Detail: true})
	if err == nil {
		t.Fatalf("expected write error, got nil")
	}
}

// TestRun_MalformedA / B return parser errors.
func TestRun_MalformedA(t *testing.T) {
	_, err := Run(strings.NewReader("chr1\tnotanumber\t10\n"), strings.NewReader(""), &bytesSink{}, Options{})
	if err == nil {
		t.Fatalf("expected error for malformed A")
	}
}

func TestRun_MalformedB(t *testing.T) {
	_, err := Run(strings.NewReader(""), strings.NewReader("chr1\tnotanumber\t10\n"), &bytesSink{}, Options{})
	if err == nil {
		t.Fatalf("expected error for malformed B")
	}
}

type bytesSink struct{}

func (*bytesSink) Write(p []byte) (int, error) { return len(p), nil }
