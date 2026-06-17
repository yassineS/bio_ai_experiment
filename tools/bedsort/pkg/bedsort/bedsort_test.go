package bedsort

import (
	"bytes"
	"strings"
	"testing"
)

func runString(t *testing.T, input string, opts Options) string {
	t.Helper()
	var buf bytes.Buffer
	if err := Run(strings.NewReader(input), &buf, opts); err != nil {
		t.Fatalf("Run failed: %v", err)
	}
	return buf.String()
}

func TestSortDefault(t *testing.T) {
	// Records with equal (chrom, start) must keep their INPUT order, not be
	// reordered by chromEnd. Here chr1:100-200 precedes chr1:100-150 in the
	// input, so it must precede it in the output too — matching upstream
	// bedtools sort, whose start-only stable sort preserves input order on
	// (chrom, start) ties.
	input := "chr2\t100\t200\n" +
		"chr1\t300\t400\n" +
		"chr1\t100\t200\n" +
		"chr1\t100\t150\n"
	want := "chr1\t100\t200\n" +
		"chr1\t100\t150\n" +
		"chr1\t300\t400\n" +
		"chr2\t100\t200\n"
	if got := runString(t, input, Options{}); got != want {
		t.Errorf("default sort mismatch.\nGot:\n%q\nWant:\n%q", got, want)
	}
}

func TestSortLexicographicChromOrder(t *testing.T) {
	// 'chr10' sorts before 'chr2' lexicographically; this matches bedtools sort.
	input := "chr2\t10\t20\nchr10\t10\t20\nchr1\t10\t20\n"
	want := "chr1\t10\t20\nchr10\t10\t20\nchr2\t10\t20\n"
	if got := runString(t, input, Options{}); got != want {
		t.Errorf("lex order mismatch.\nGot:\n%q\nWant:\n%q", got, want)
	}
}

func TestSortPreservesAllColumns(t *testing.T) {
	// BED6 input: name/score/strand must round-trip through the sort.
	input := "chr1\t300\t400\tg1\t500\t+\n" +
		"chr1\t100\t200\tg2\t900\t-\textra\n"
	want := "chr1\t100\t200\tg2\t900\t-\textra\n" +
		"chr1\t300\t400\tg1\t500\t+\n"
	if got := runString(t, input, Options{}); got != want {
		t.Errorf("column round-trip mismatch.\nGot:\n%q\nWant:\n%q", got, want)
	}
}

func TestSortIgnoresHeadersAndBlanks(t *testing.T) {
	input := "track name=foo\n" +
		"# header comment\n" +
		"\n" +
		"chr1\t10\t20\n" +
		"browser position chr1\n" +
		"chr1\t5\t8\n"
	want := "chr1\t5\t8\nchr1\t10\t20\n"
	if got := runString(t, input, Options{}); got != want {
		t.Errorf("header filter mismatch.\nGot:\n%q\nWant:\n%q", got, want)
	}
}

func TestSortSizeAsc(t *testing.T) {
	input := "chr1\t0\t100\n" +
		"chr2\t0\t10\n" +
		"chr1\t0\t50\n"
	want := "chr2\t0\t10\nchr1\t0\t50\nchr1\t0\t100\n"
	if got := runString(t, input, Options{Mode: ModeSizeA}); got != want {
		t.Errorf("sizeA mismatch.\nGot:\n%q\nWant:\n%q", got, want)
	}
}

func TestSortSizeDesc(t *testing.T) {
	input := "chr1\t0\t100\n" +
		"chr2\t0\t10\n" +
		"chr1\t0\t50\n"
	want := "chr1\t0\t100\nchr1\t0\t50\nchr2\t0\t10\n"
	if got := runString(t, input, Options{Mode: ModeSizeD}); got != want {
		t.Errorf("sizeD mismatch.\nGot:\n%q\nWant:\n%q", got, want)
	}
}

func TestSortChrThenSize(t *testing.T) {
	input := "chr1\t0\t30\n" +
		"chr1\t0\t10\n" +
		"chr2\t0\t100\n" +
		"chr2\t0\t20\n"
	wantAsc := "chr1\t0\t10\nchr1\t0\t30\nchr2\t0\t20\nchr2\t0\t100\n"
	if got := runString(t, input, Options{Mode: ModeChrThenSizeA}); got != wantAsc {
		t.Errorf("chrThenSizeA mismatch.\nGot:\n%q\nWant:\n%q", got, wantAsc)
	}
	wantDesc := "chr1\t0\t30\nchr1\t0\t10\nchr2\t0\t100\nchr2\t0\t20\n"
	if got := runString(t, input, Options{Mode: ModeChrThenSizeD}); got != wantDesc {
		t.Errorf("chrThenSizeD mismatch.\nGot:\n%q\nWant:\n%q", got, wantDesc)
	}
}

func TestSortChrThenScore(t *testing.T) {
	input := "chr1\t0\t100\tg1\t30\t+\n" +
		"chr1\t10\t100\tg2\t10\t+\n" +
		"chr2\t0\t100\tg3\t900\t+\n" +
		"chr2\t0\t100\tg4\t200\t+\n"
	wantAsc := "chr1\t10\t100\tg2\t10\t+\n" +
		"chr1\t0\t100\tg1\t30\t+\n" +
		"chr2\t0\t100\tg4\t200\t+\n" +
		"chr2\t0\t100\tg3\t900\t+\n"
	if got := runString(t, input, Options{Mode: ModeChrThenScoreA}); got != wantAsc {
		t.Errorf("chrThenScoreA mismatch.\nGot:\n%q\nWant:\n%q", got, wantAsc)
	}
	wantDesc := "chr1\t0\t100\tg1\t30\t+\n" +
		"chr1\t10\t100\tg2\t10\t+\n" +
		"chr2\t0\t100\tg3\t900\t+\n" +
		"chr2\t0\t100\tg4\t200\t+\n"
	if got := runString(t, input, Options{Mode: ModeChrThenScoreD}); got != wantDesc {
		t.Errorf("chrThenScoreD mismatch.\nGot:\n%q\nWant:\n%q", got, wantDesc)
	}
}

func TestSortChrThenScoreMissingScore(t *testing.T) {
	// BED3 rows have no score and should sort with effective score=0,
	// preserving input order under stable sort.
	input := "chr1\t0\t30\n" +
		"chr1\t10\t40\n"
	want := "chr1\t0\t30\nchr1\t10\t40\n"
	if got := runString(t, input, Options{Mode: ModeChrThenScoreA}); got != want {
		t.Errorf("chrThenScoreA with no score mismatch.\nGot:\n%q\nWant:\n%q", got, want)
	}
}

func TestSortChromOrderFromFaidx(t *testing.T) {
	order, err := ReadFaidx(strings.NewReader("chrX\t100\nchr2\t100\nchr1\t100\n"))
	if err != nil {
		t.Fatalf("ReadFaidx: %v", err)
	}
	wantOrder := []string{"chrX", "chr2", "chr1"}
	if len(order) != len(wantOrder) {
		t.Fatalf("ReadFaidx length mismatch: %v", order)
	}
	for i, c := range wantOrder {
		if order[i] != c {
			t.Errorf("ReadFaidx[%d] = %q, want %q", i, order[i], c)
		}
	}
	input := "chr1\t0\t10\nchr2\t0\t10\nchrX\t0\t10\n"
	want := "chrX\t0\t10\nchr2\t0\t10\nchr1\t0\t10\n"
	if got := runString(t, input, Options{ChromOrder: order}); got != want {
		t.Errorf("custom chrom order mismatch.\nGot:\n%q\nWant:\n%q", got, want)
	}
}

func TestSortFaidxIgnoresCommentsAndBlanks(t *testing.T) {
	in := "# comment\n\nchr1 100\nchr2 200\nchr1 200\n"
	got, err := ReadFaidx(strings.NewReader(in))
	if err != nil {
		t.Fatalf("ReadFaidx: %v", err)
	}
	want := []string{"chr1", "chr2"}
	if len(got) != len(want) || got[0] != want[0] || got[1] != want[1] {
		t.Errorf("got %v, want %v", got, want)
	}
}

func TestSortChromOrderUnlistedSortLast(t *testing.T) {
	// chrUn is not in ChromOrder and should sort *after* chr1/chr2, in lex
	// order with any other unlisted chromosomes.
	input := "chrUn\t0\t10\nchr2\t0\t10\nchr1\t0\t10\nchrZ\t0\t10\n"
	want := "chr1\t0\t10\nchr2\t0\t10\nchrUn\t0\t10\nchrZ\t0\t10\n"
	got := runString(t, input, Options{ChromOrder: []string{"chr1", "chr2"}})
	if got != want {
		t.Errorf("unlisted-chrom order mismatch.\nGot:\n%q\nWant:\n%q", got, want)
	}
}

func TestSortInvalidStart(t *testing.T) {
	_, err := ReadAll(strings.NewReader("chr1\tNOT_A_NUMBER\t20\n"))
	if err == nil {
		t.Errorf("expected error on bad chromStart, got nil")
	}
}

func TestSortInvalidEnd(t *testing.T) {
	_, err := ReadAll(strings.NewReader("chr1\t10\tNOT_A_NUMBER\n"))
	if err == nil {
		t.Errorf("expected error on bad chromEnd, got nil")
	}
}

func TestSortTooFewFields(t *testing.T) {
	_, err := ReadAll(strings.NewReader("chr1 10 20\n")) // space-separated -> 1 field
	if err == nil {
		t.Errorf("expected error on too-few fields, got nil")
	}
}

func TestSortEmptyInput(t *testing.T) {
	got := runString(t, "", Options{})
	if got != "" {
		t.Errorf("empty input should produce empty output, got %q", got)
	}
}

func TestSortHeaderPreservesLeadingHeaderLines(t *testing.T) {
	// -header: leading comment/track/browser lines are kept verbatim before
	// the sorted body. Mid-file header-style lines are still dropped.
	input := "#Header line\n" +
		"track name=foo\n" +
		"browser position chr1\n" +
		"chr2\t10\t20\n" +
		"# mid-file comment\n" +
		"chr1\t5\t8\n"
	want := "#Header line\n" +
		"track name=foo\n" +
		"browser position chr1\n" +
		"chr1\t5\t8\n" +
		"chr2\t10\t20\n"
	if got := runString(t, input, Options{Header: true}); got != want {
		t.Errorf("-header mismatch.\nGot:\n%q\nWant:\n%q", got, want)
	}
}

func TestSortHeaderEmptyHeaderBlock(t *testing.T) {
	// -header on input with no header lines should be a no-op.
	input := "chr2\t10\t20\nchr1\t5\t8\n"
	want := "chr1\t5\t8\nchr2\t10\t20\n"
	if got := runString(t, input, Options{Header: true}); got != want {
		t.Errorf("-header empty-block mismatch.\nGot:\n%q\nWant:\n%q", got, want)
	}
}

func TestSortStableWithinTies(t *testing.T) {
	// Same chrom/start/end: ordering must be preserved (stable sort).
	input := "chr1\t10\t20\tfirst\n" +
		"chr1\t10\t20\tsecond\n" +
		"chr1\t10\t20\tthird\n"
	want := input
	if got := runString(t, input, Options{}); got != want {
		t.Errorf("stability mismatch.\nGot:\n%q\nWant:\n%q", got, want)
	}
}
