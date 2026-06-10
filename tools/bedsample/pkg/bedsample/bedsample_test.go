package bedsample

import (
	"bytes"
	"errors"
	"fmt"
	"strconv"
	"strings"
	"testing"
)

// buildInput emits n BED rows of the form "chr1\t<start>\t<end>\trec_<i>\n".
// Each row's first integer is monotonic so output-order assertions are easy.
func buildInput(n int) string {
	var sb strings.Builder
	for i := 0; i < n; i++ {
		fmt.Fprintf(&sb, "chr1\t%d\t%d\trec_%s\n", i*10, i*10+5, strconv.Itoa(i))
	}
	return sb.String()
}

func TestSample_BasicCount(t *testing.T) {
	in := buildInput(100)
	var buf bytes.Buffer
	n, err := Sample(strings.NewReader(in), &buf, Options{N: 10, Seed: 42})
	if err != nil {
		t.Fatalf("Sample: %v", err)
	}
	if n != 10 {
		t.Errorf("returned n=%d, want 10", n)
	}
	got := strings.Count(buf.String(), "\n")
	if got != 10 {
		t.Errorf("output line count = %d, want 10", got)
	}
}

func TestSample_Deterministic(t *testing.T) {
	in := buildInput(200)
	var a, b bytes.Buffer
	if _, err := Sample(strings.NewReader(in), &a, Options{N: 50, Seed: 4}); err != nil {
		t.Fatalf("Sample a: %v", err)
	}
	if _, err := Sample(strings.NewReader(in), &b, Options{N: 50, Seed: 4}); err != nil {
		t.Fatalf("Sample b: %v", err)
	}
	if a.String() != b.String() {
		t.Errorf("seeded runs disagree:\na=%s\nb=%s", a.String(), b.String())
	}
}

func TestSample_DifferentSeeds(t *testing.T) {
	in := buildInput(200)
	var a, b bytes.Buffer
	if _, err := Sample(strings.NewReader(in), &a, Options{N: 50, Seed: 1}); err != nil {
		t.Fatalf("Sample a: %v", err)
	}
	if _, err := Sample(strings.NewReader(in), &b, Options{N: 50, Seed: 2}); err != nil {
		t.Fatalf("Sample b: %v", err)
	}
	if a.String() == b.String() {
		t.Errorf("different seeds produced identical samples (improbable)")
	}
}

func TestSample_NoReplacement(t *testing.T) {
	in := buildInput(50)
	var buf bytes.Buffer
	if _, err := Sample(strings.NewReader(in), &buf, Options{N: 50, Seed: 7}); err != nil {
		t.Fatalf("Sample: %v", err)
	}
	// All 50 input records should be present exactly once.
	seen := make(map[string]int)
	for _, line := range strings.Split(strings.TrimRight(buf.String(), "\n"), "\n") {
		seen[line]++
	}
	for line, count := range seen {
		if count != 1 {
			t.Errorf("record %q appears %d times (want 1)", line, count)
		}
	}
	if len(seen) != 50 {
		t.Errorf("unique records = %d, want 50", len(seen))
	}
}

func TestSample_OutputIsSubsetOfInput(t *testing.T) {
	// Upstream emits sampled BED records in reservoir-slot order, NOT input
	// order (it only re-sorts when the output type is BAM). So the only
	// invariants we can assert are: every output line is a distinct input
	// record, and the count is exactly N.
	in := buildInput(50)
	inputSet := make(map[string]bool)
	for _, line := range strings.Split(strings.TrimRight(in, "\n"), "\n") {
		inputSet[line] = true
	}

	var buf bytes.Buffer
	if _, err := Sample(strings.NewReader(in), &buf, Options{N: 20, Seed: 99}); err != nil {
		t.Fatalf("Sample: %v", err)
	}
	lines := strings.Split(strings.TrimRight(buf.String(), "\n"), "\n")
	if len(lines) != 20 {
		t.Fatalf("got %d lines, want 20", len(lines))
	}
	seen := make(map[string]bool)
	for _, line := range lines {
		if !inputSet[line] {
			t.Errorf("output line %q is not an input record", line)
		}
		if seen[line] {
			t.Errorf("output line %q appears more than once", line)
		}
		seen[line] = true
	}
}

func TestSample_FillPhasePreservesInputOrder(t *testing.T) {
	// When N >= total records the reservoir never enters the replacement
	// phase, so the slot order equals input order. This documents upstream's
	// fill-phase behaviour (keepRecord without an RNG draw).
	in := buildInput(15)
	var buf bytes.Buffer
	if _, err := Sample(strings.NewReader(in), &buf, Options{N: 15, Seed: 99}); err != nil {
		t.Fatalf("Sample: %v", err)
	}
	lines := strings.Split(strings.TrimRight(buf.String(), "\n"), "\n")
	prev := -1
	for _, line := range lines {
		parts := strings.Split(line, "\t")
		if len(parts) < 2 {
			t.Fatalf("malformed line %q", line)
		}
		n, err := strconv.Atoi(parts[1])
		if err != nil {
			t.Fatalf("bad int %q: %v", parts[1], err)
		}
		if n <= prev {
			t.Errorf("fill-phase output not in input order: %d after %d", n, prev)
		}
		prev = n
	}
}

func TestSample_NEqualsTotal(t *testing.T) {
	in := buildInput(7)
	var buf bytes.Buffer
	n, err := Sample(strings.NewReader(in), &buf, Options{N: 7, Seed: 5})
	if err != nil {
		t.Fatalf("Sample: %v", err)
	}
	if n != 7 {
		t.Errorf("n=%d, want 7", n)
	}
	lines := strings.Split(strings.TrimRight(buf.String(), "\n"), "\n")
	if len(lines) != 7 {
		t.Errorf("lines=%d, want 7", len(lines))
	}
}

func TestSample_TooFewRecords(t *testing.T) {
	in := buildInput(3)
	var buf bytes.Buffer
	_, err := Sample(strings.NewReader(in), &buf, Options{N: 10, Seed: 1})
	if err == nil {
		t.Errorf("expected ErrTooFewRecords")
	}
	var tooFew *ErrTooFewRecords
	if !errors.As(err, &tooFew) {
		t.Errorf("err type = %T, want *ErrTooFewRecords", err)
	} else if tooFew.Have != 3 || tooFew.Want != 10 {
		t.Errorf("ErrTooFewRecords{Have=%d, Want=%d}, want {3,10}", tooFew.Have, tooFew.Want)
	}
}

func TestSample_NonPositiveN(t *testing.T) {
	var buf bytes.Buffer
	if _, err := Sample(strings.NewReader("chr1\t0\t1\n"), &buf, Options{N: 0}); err == nil {
		t.Errorf("expected error for N=0")
	}
	if _, err := Sample(strings.NewReader("chr1\t0\t1\n"), &buf, Options{N: -1}); err == nil {
		t.Errorf("expected error for N<0")
	}
}

func TestSample_HeaderForwardingOn(t *testing.T) {
	in := "#hdr1\ntrack name=foo\nchr1\t0\t1\nchr1\t2\t3\nchr1\t4\t5\n"
	var buf bytes.Buffer
	if _, err := Sample(strings.NewReader(in), &buf, Options{N: 1, Seed: 1, Header: true}); err != nil {
		t.Fatalf("Sample: %v", err)
	}
	out := buf.String()
	if !strings.HasPrefix(out, "#hdr1\ntrack name=foo\n") {
		t.Errorf("expected header lines first; got:\n%s", out)
	}
	if strings.Count(out, "\n") != 3 {
		t.Errorf("expected 3 lines (2 header + 1 sampled); got\n%s", out)
	}
}

func TestSample_HeaderForwardingOff(t *testing.T) {
	in := "#hdr1\ntrack name=foo\nchr1\t0\t1\nchr1\t2\t3\nchr1\t4\t5\n"
	var buf bytes.Buffer
	if _, err := Sample(strings.NewReader(in), &buf, Options{N: 1, Seed: 1}); err != nil {
		t.Fatalf("Sample: %v", err)
	}
	if strings.HasPrefix(buf.String(), "#hdr") {
		t.Errorf("header should not be forwarded by default; got:\n%s", buf.String())
	}
}

func TestSample_BlankLinesIgnored(t *testing.T) {
	in := "\nchr1\t0\t1\n\nchr1\t2\t3\n\n"
	var buf bytes.Buffer
	if _, err := Sample(strings.NewReader(in), &buf, Options{N: 2, Seed: 1}); err != nil {
		t.Fatalf("Sample: %v", err)
	}
	if strings.Count(buf.String(), "\n") != 2 {
		t.Errorf("blanks should be dropped; got\n%s", buf.String())
	}
}

func TestSample_TimeSeed(t *testing.T) {
	// Seed=0 → time-based seed. Just verify it doesn't error.
	in := buildInput(20)
	var buf bytes.Buffer
	if _, err := Sample(strings.NewReader(in), &buf, Options{N: 5, Seed: 0}); err != nil {
		t.Fatalf("Sample: %v", err)
	}
}

type errReader struct{}

var errForcedRead = errors.New("forced read error")

func (errReader) Read(_ []byte) (int, error) { return 0, errForcedRead }

func TestSample_ReaderError(t *testing.T) {
	var buf bytes.Buffer
	if _, err := Sample(errReader{}, &buf, Options{N: 3, Seed: 1}); err == nil {
		t.Errorf("expected reader error to propagate")
	}
}

func TestErrTooFewRecords_Error(t *testing.T) {
	e := &ErrTooFewRecords{Have: 1, Want: 5}
	if e.Error() == "" {
		t.Errorf("Error() empty")
	}
}
