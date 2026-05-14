package bedspacing

import (
	"bytes"
	"errors"
	"strings"
	"testing"
)

func runSpacing(t *testing.T, in string) string {
	t.Helper()
	var buf bytes.Buffer
	if _, err := Spacing(strings.NewReader(in), &buf); err != nil {
		t.Fatalf("Spacing: %v", err)
	}
	return buf.String()
}

func TestSpacing_BasicGaps(t *testing.T) {
	in := "chr1\t0\t10\nchr1\t20\t30\n"
	want := "chr1\t0\t10\t.\nchr1\t20\t30\t10\n"
	if got := runSpacing(t, in); got != want {
		t.Errorf("output mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

func TestSpacing_Overlap(t *testing.T) {
	in := "chr1\t0\t20\nchr1\t10\t30\n"
	want := "chr1\t0\t20\t.\nchr1\t10\t30\t-1\n"
	if got := runSpacing(t, in); got != want {
		t.Errorf("output mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

func TestSpacing_Abuts(t *testing.T) {
	in := "chr1\t0\t10\nchr1\t10\t20\n"
	want := "chr1\t0\t10\t.\nchr1\t10\t20\t0\n"
	if got := runSpacing(t, in); got != want {
		t.Errorf("output mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

func TestSpacing_FirstPerChrom(t *testing.T) {
	// Each chromosome's first row is ".".
	in := "chr1\t0\t10\nchr2\t100\t200\n"
	want := "chr1\t0\t10\t.\nchr2\t100\t200\t.\n"
	if got := runSpacing(t, in); got != want {
		t.Errorf("output mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

func TestSpacing_PreservesAllColumns(t *testing.T) {
	in := "chr1\t0\t10\tfoo\t100\t+\nchr1\t20\t30\tbar\t200\t-\n"
	want := "chr1\t0\t10\tfoo\t100\t+\t.\nchr1\t20\t30\tbar\t200\t-\t10\n"
	if got := runSpacing(t, in); got != want {
		t.Errorf("output mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

func TestSpacing_HeaderLinesPassthrough(t *testing.T) {
	in := "# comment\ntrack name=foo\nchr1\t0\t10\nchr1\t20\t30\n"
	want := "# comment\ntrack name=foo\nchr1\t0\t10\t.\nchr1\t20\t30\t10\n"
	if got := runSpacing(t, in); got != want {
		t.Errorf("output mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

func TestSpacing_BlankLinesPassthrough(t *testing.T) {
	in := "\nchr1\t0\t10\n\nchr1\t20\t30\n"
	want := "\nchr1\t0\t10\t.\n\nchr1\t20\t30\t10\n"
	if got := runSpacing(t, in); got != want {
		t.Errorf("output mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

func TestSpacing_EmptyInput(t *testing.T) {
	if got := runSpacing(t, ""); got != "" {
		t.Errorf("empty input should produce empty output, got %q", got)
	}
}

func TestSpacing_SingleRecord(t *testing.T) {
	if got := runSpacing(t, "chr1\t0\t10\n"); got != "chr1\t0\t10\t.\n" {
		t.Errorf("single record output = %q", got)
	}
}

func TestSpacing_PreviousPointerAdvances(t *testing.T) {
	// Per upstream, the "previous record" is only the immediately previous
	// one on the chrom. A contained interval rewinds the running end.
	in := "chr1\t60\t80\nchr1\t65\t70\nchr1\t75\t100\n"
	want := "chr1\t60\t80\t.\nchr1\t65\t70\t-1\nchr1\t75\t100\t5\n"
	if got := runSpacing(t, in); got != want {
		t.Errorf("output mismatch.\nwant:\n%s\ngot:\n%s", want, got)
	}
}

func TestSpacing_TooFewColumns(t *testing.T) {
	var buf bytes.Buffer
	_, err := Spacing(strings.NewReader("chr1\t0\n"), &buf)
	if err == nil {
		t.Errorf("expected error for 2-column row")
	}
}

func TestSpacing_BadStart(t *testing.T) {
	var buf bytes.Buffer
	_, err := Spacing(strings.NewReader("chr1\tabc\t10\n"), &buf)
	if err == nil {
		t.Errorf("expected error for non-numeric chromStart")
	}
}

func TestSpacing_BadEnd(t *testing.T) {
	var buf bytes.Buffer
	_, err := Spacing(strings.NewReader("chr1\t0\txyz\n"), &buf)
	if err == nil {
		t.Errorf("expected error for non-numeric chromEnd")
	}
}

// erroringReader returns an error on read; exercises the scanner.Err path.
type erroringReader struct{}

var errForcedRead = errors.New("forced read error")

func (erroringReader) Read(_ []byte) (int, error) { return 0, errForcedRead }

func TestSpacing_ReaderError(t *testing.T) {
	var buf bytes.Buffer
	_, err := Spacing(erroringReader{}, &buf)
	if err == nil {
		t.Errorf("expected reader error to propagate")
	}
}

// limitedWriter fails after writing N bytes; lets us exercise write errors.
type limitedWriter struct {
	n      int
	remain int
}

var errLimitedWriter = errors.New("write capacity exhausted")

func (w *limitedWriter) Write(p []byte) (int, error) {
	if w.remain <= 0 {
		return 0, errLimitedWriter
	}
	if len(p) > w.remain {
		w.remain = 0
		return 0, errLimitedWriter
	}
	w.remain -= len(p)
	w.n += len(p)
	return len(p), nil
}

func TestSpacing_WriterErrorPath(t *testing.T) {
	// Just ensure that Spacing doesn't panic if the underlying writer
	// errors; the buffered writer will hide errors until flush time. We
	// don't assert on the specific error because bufio swallows it until
	// flush — but at minimum the call must not panic.
	defer func() {
		if r := recover(); r != nil {
			t.Errorf("Spacing panicked on writer error: %v", r)
		}
	}()
	// A 0-byte writer never accepts anything: bufio will eventually flush
	// and surface a downstream error.
	_, _ = Spacing(strings.NewReader("chr1\t0\t10\nchr1\t20\t30\n"), &limitedWriter{remain: 0})
}
