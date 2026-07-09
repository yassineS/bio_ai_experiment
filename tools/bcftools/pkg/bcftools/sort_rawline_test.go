package bcftools

import (
	"bytes"
	"strings"
	"testing"
)

// TestSortRawLineMatchesParsePath pins the sort -O v raw-line fast path: sorting
// through readAllVariantsRawLine (KeepRawLine, verbatim re-emit) must produce the
// exact same output as sorting through the full parse→re-encode path. For a
// well-formed record the two are byte-identical, so this compares the -O v output
// against a reference produced by the parse path over the same input.
func TestSortRawLineMatchesParsePath(t *testing.T) {
	in := makeSortVCF(
		"chr10\t100\trsX\tA\tT\t.\tPASS\tDP=1\tGT\t0/1",
		"chr2\t50\trsB\tA\tG,C\t.\tPASS\tDP=2\tGT\t1/2",
		"chr1\t300\trsC\tA\tT\t60\tPASS\tDP=1\tGT\t0/1",
		"chr1\t100\trsA\tA\tT\t.\tPASS\tDP=1\tGT\t0/1",
	)

	// Raw-line path (-O v is OutputVCF, the default).
	var raw bytes.Buffer
	if _, err := Sort(strings.NewReader(in), &raw, SortOptions{OutputFormat: OutputVCF}); err != nil {
		t.Fatalf("raw-line Sort: %v", err)
	}

	// Reference: read via the full parse path, sort, and re-encode, mirroring what
	// the pre-fix code did.
	hdr, recs, err := readAllVariants(strings.NewReader(in))
	if err != nil {
		t.Fatalf("readAllVariants: %v", err)
	}
	order := contigOrder(hdr)
	sortVariantsForSort(recs, order)
	var ref bytes.Buffer
	w, finish, err := openOutput(&ref, ViewOptions{OutputFormat: OutputVCF}, hdr)
	if err != nil {
		t.Fatalf("openOutput: %v", err)
	}
	if err := w.WriteHeader(); err != nil {
		t.Fatalf("WriteHeader: %v", err)
	}
	for _, v := range recs {
		if err := w.Write(v); err != nil {
			t.Fatalf("Write: %v", err)
		}
	}
	if err := w.Flush(); err != nil {
		t.Fatalf("Flush: %v", err)
	}
	finish()

	if raw.String() != ref.String() {
		t.Fatalf("raw-line -O v output differs from parse path:\n--raw-line--\n%s\n--parse--\n%s", raw.String(), ref.String())
	}
}

// TestSortRawLineIsVerbatim confirms the raw-line path re-emits the ORIGINAL data
// line byte-for-byte (not a re-serialised form), so any field the parser would
// normalise is preserved. It sorts records whose ID/INFO ordering is already in
// input form and checks each output data line is one of the input lines verbatim.
func TestSortRawLineIsVerbatim(t *testing.T) {
	lines := []string{
		"chr1\t200\trsB\tA\tT\t.\tPASS\tDP=5;AF=0.5\tGT\t0/1",
		"chr1\t100\trsA\tC\tG\t.\tPASS\tAF=0.1;DP=9\tGT\t1/1",
	}
	in := makeSortVCF(lines...)

	var out bytes.Buffer
	if _, err := Sort(strings.NewReader(in), &out, SortOptions{OutputFormat: OutputVCF}); err != nil {
		t.Fatalf("Sort: %v", err)
	}

	inputSet := map[string]bool{}
	for _, l := range lines {
		inputSet[l] = true
	}
	var dataLines []string
	for _, l := range strings.Split(strings.TrimRight(out.String(), "\n"), "\n") {
		if strings.HasPrefix(l, "#") {
			continue
		}
		dataLines = append(dataLines, l)
	}
	if len(dataLines) != len(lines) {
		t.Fatalf("expected %d data lines, got %d", len(lines), len(dataLines))
	}
	// After sorting by POS the order is 100 then 200.
	if dataLines[0] != lines[1] || dataLines[1] != lines[0] {
		t.Fatalf("data lines are not the verbatim input lines in sorted order:\n%q\n%q", dataLines[0], dataLines[1])
	}
}
