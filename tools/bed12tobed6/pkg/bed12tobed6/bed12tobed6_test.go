package bed12tobed6

import (
	"bytes"
	"errors"
	"strings"
	"testing"
)

// TestConvert_DifferingFieldCounts is the parity test for upstream bedtools'
// "Differing number of BED fields encountered" behaviour: the first data record
// fixes the expected column count and any later record with a different count
// aborts. Records emitted before the offending line are preserved (streaming),
// matching upstream's stdout on exit.
func TestConvert_DifferingFieldCounts(t *testing.T) {
	// First record BED12 (12 cols) → establishes 12; second record has 6 cols.
	in := "chr1\t100\t200\tn1\t0\t+\t100\t200\t0\t2\t10,20\t0,50\nchr1\t300\t400\tn2\t0\t+\n"
	wantOut := "chr1\t100\t110\tn1\t0\t+\nchr1\t150\t170\tn1\t0\t+\n"

	var out bytes.Buffer
	n, err := Convert(strings.NewReader(in), &out, Options{})
	if err == nil {
		t.Fatal("expected a FieldCountError for a mixed-field file, got nil")
	}
	var fcErr *FieldCountError
	if !errors.As(err, &fcErr) {
		t.Fatalf("error type: got %T (%v), want *FieldCountError", err, err)
	}
	if fcErr.Line != 2 {
		t.Errorf("FieldCountError.Line = %d, want 2", fcErr.Line)
	}
	if want := "Differing number of BED fields encountered at line: 2.  Exiting..."; fcErr.Error() != want {
		t.Errorf("message = %q, want %q", fcErr.Error(), want)
	}
	// Blocks from the first (valid) record must still have been written.
	if n != 2 {
		t.Errorf("written = %d, want 2 (first record's blocks)", n)
	}
	if got := out.String(); got != wantOut {
		t.Errorf("partial output = %q, want %q", got, wantOut)
	}
}

// TestConvert_UniformSixColUnchanged guards against the field-count check
// firing on a well-formed uniform file (all records share a column count).
func TestConvert_UniformSixColUnchanged(t *testing.T) {
	in := "chr1\t10\t20\ta\t0\t+\nchr2\t30\t40\tb\t0\t-\n"
	want := "chr1\t10\t20\ta\t0\t+\nchr2\t30\t40\tb\t0\t-\n"
	var out bytes.Buffer
	if _, err := Convert(strings.NewReader(in), &out, Options{}); err != nil {
		t.Fatalf("unexpected err on uniform file: %v", err)
	}
	if got := out.String(); got != want {
		t.Errorf("output = %q, want %q", got, want)
	}
}

func TestConvert_OneBlock(t *testing.T) {
	in := "chr1\t0\t50\tone_blocks_match\t0\t+\t0\t0\t0\t1\t50,\t0,\n"
	want := "chr1\t0\t50\tone_blocks_match\t0\t+\n"
	var out bytes.Buffer
	n, err := Convert(strings.NewReader(in), &out, Options{})
	if err != nil {
		t.Fatalf("unexpected err: %v", err)
	}
	if n != 1 {
		t.Fatalf("want 1 record, got %d", n)
	}
	if got := out.String(); got != want {
		t.Fatalf("output mismatch:\nwant=%q\ngot =%q", want, got)
	}
}

func TestConvert_TwoBlocks(t *testing.T) {
	in := "chr1\t0\t50\ttwo_blocks_match\t0\t+\t0\t0\t0\t2\t10,10,\t0,40,\n"
	want := "chr1\t0\t10\ttwo_blocks_match\t0\t+\nchr1\t40\t50\ttwo_blocks_match\t0\t+\n"
	var out bytes.Buffer
	n, err := Convert(strings.NewReader(in), &out, Options{})
	if err != nil {
		t.Fatalf("unexpected err: %v", err)
	}
	if n != 2 {
		t.Fatalf("want 2 records, got %d", n)
	}
	if got := out.String(); got != want {
		t.Fatalf("output mismatch:\nwant=%q\ngot =%q", want, got)
	}
}

func TestConvert_ThreeBlocks(t *testing.T) {
	in := "chr1\t0\t50\tthree_blocks_match\t0\t+\t0\t0\t0\t3\t10,10,10,\t0,20,40,\n"
	want := "chr1\t0\t10\tthree_blocks_match\t0\t+\n" +
		"chr1\t20\t30\tthree_blocks_match\t0\t+\n" +
		"chr1\t40\t50\tthree_blocks_match\t0\t+\n"
	var out bytes.Buffer
	if _, err := Convert(strings.NewReader(in), &out, Options{}); err != nil {
		t.Fatalf("err: %v", err)
	}
	if got := out.String(); got != want {
		t.Fatalf("mismatch:\nwant=%q\ngot =%q", want, got)
	}
}

func TestConvert_NumberBlocksForward(t *testing.T) {
	in := "chr1\t0\t50\tthree_blocks_match\t0\t+\t0\t0\t0\t3\t10,10,10,\t0,20,40,\n"
	want := "chr1\t0\t10\tthree_blocks_match\t1\t+\n" +
		"chr1\t20\t30\tthree_blocks_match\t2\t+\n" +
		"chr1\t40\t50\tthree_blocks_match\t3\t+\n"
	var out bytes.Buffer
	if _, err := Convert(strings.NewReader(in), &out, Options{NumberBlocks: true}); err != nil {
		t.Fatalf("err: %v", err)
	}
	if got := out.String(); got != want {
		t.Fatalf("mismatch:\nwant=%q\ngot =%q", want, got)
	}
}

func TestConvert_NumberBlocksReverseStrand(t *testing.T) {
	in := "chr1\t0\t50\tthree_blocks_match\t0\t-\t0\t0\t0\t3\t10,10,10,\t0,20,40,\n"
	want := "chr1\t0\t10\tthree_blocks_match\t3\t-\n" +
		"chr1\t20\t30\tthree_blocks_match\t2\t-\n" +
		"chr1\t40\t50\tthree_blocks_match\t1\t-\n"
	var out bytes.Buffer
	if _, err := Convert(strings.NewReader(in), &out, Options{NumberBlocks: true}); err != nil {
		t.Fatalf("err: %v", err)
	}
	if got := out.String(); got != want {
		t.Fatalf("mismatch:\nwant=%q\ngot =%q", want, got)
	}
}

// TestUnit_ScorePropagation checks (binary-free) that the parent record's
// score column is carried unchanged onto each emitted BED6 block, not zeroed.
func TestUnit_ScorePropagation(t *testing.T) {
	in := "chr1\t100\t300\tfeatA\t500\t+\t100\t300\t0,0,0\t2\t50,50\t0,150\n"
	want := "chr1\t100\t150\tfeatA\t500\t+\n" +
		"chr1\t250\t300\tfeatA\t500\t+\n"
	var out bytes.Buffer
	if _, err := Convert(strings.NewReader(in), &out, Options{}); err != nil {
		t.Fatalf("err: %v", err)
	}
	if got := out.String(); got != want {
		t.Fatalf("score not propagated:\nwant=%q\ngot =%q", want, got)
	}
}

// TestUnit_NumberBlocksDotStrand checks that -n reverses numbering for any
// non-"+" strand, including ".", matching upstream's `strand == "+"` test.
func TestUnit_NumberBlocksDotStrand(t *testing.T) {
	in := "chr3\t0\t120\tfeatC\t42\t.\t0\t120\t0,0,0\t2\t20,30\t0,90\n"
	want := "chr3\t0\t20\tfeatC\t2\t.\n" +
		"chr3\t90\t120\tfeatC\t1\t.\n"
	var out bytes.Buffer
	if _, err := Convert(strings.NewReader(in), &out, Options{NumberBlocks: true}); err != nil {
		t.Fatalf("err: %v", err)
	}
	if got := out.String(); got != want {
		t.Fatalf("dot-strand numbering mismatch:\nwant=%q\ngot =%q", want, got)
	}
}

func TestConvert_EmptyInput(t *testing.T) {
	var out bytes.Buffer
	n, err := Convert(strings.NewReader(""), &out, Options{})
	if err != nil {
		t.Fatalf("err: %v", err)
	}
	if n != 0 {
		t.Fatalf("want 0 records, got %d", n)
	}
}

func TestConvert_SkipsHeadersAndComments(t *testing.T) {
	in := "# header\ntrack name=foo\nbrowser pos chr1\n\n" +
		"chr1\t0\t50\tname\t0\t+\t0\t0\t0\t1\t50,\t0,\n"
	var out bytes.Buffer
	n, err := Convert(strings.NewReader(in), &out, Options{})
	if err != nil {
		t.Fatalf("err: %v", err)
	}
	if n != 1 {
		t.Fatalf("want 1 record, got %d", n)
	}
}

// TestConvert_Bed6RoundTrips confirms a full BED6 record (6 columns, not 12)
// is re-emitted identically as a normalised single-block BED6.
func TestConvert_Bed6RoundTrips(t *testing.T) {
	in := "chr1\t0\t50\tbed6\t1\t+\n"
	var out bytes.Buffer
	if _, err := Convert(strings.NewReader(in), &out, Options{}); err != nil {
		t.Fatalf("err: %v", err)
	}
	if got := out.String(); got != in {
		t.Fatalf("BED6 should round-trip unchanged; got %q", got)
	}
}

// TestConvert_NormaliseNonBed12 verifies that upstream's behaviour of treating
// any non-12-column record as a single whole-feature block, then always
// printing six columns with empty defaults for absent trailing fields, is
// matched byte-for-byte for BED3, BED4 and BED5 inputs — with and without -n.
func TestConvert_NormaliseNonBed12(t *testing.T) {
	cases := []struct {
		name  string
		in    string
		want  string
		wantN string // expected output with NumberBlocks
	}{
		{
			name:  "bed3",
			in:    "chr20\t81335\t101267\n",
			want:  "chr20\t81335\t101267\t\t\t\n",
			wantN: "chr20\t81335\t101267\t\t1\t\n",
		},
		{
			name:  "bed4",
			in:    "chr1\t10\t20\tfeat\n",
			want:  "chr1\t10\t20\tfeat\t\t\n",
			wantN: "chr1\t10\t20\tfeat\t1\t\n",
		},
		{
			name:  "bed5",
			in:    "chr1\t10\t20\tfeat\t500\n",
			want:  "chr1\t10\t20\tfeat\t500\t\n",
			wantN: "chr1\t10\t20\tfeat\t1\t\n",
		},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			var out bytes.Buffer
			if _, err := Convert(strings.NewReader(tc.in), &out, Options{}); err != nil {
				t.Fatalf("err: %v", err)
			}
			if got := out.String(); got != tc.want {
				t.Fatalf("normalise mismatch:\nwant=%q\ngot =%q", tc.want, got)
			}
			var outN bytes.Buffer
			if _, err := Convert(strings.NewReader(tc.in), &outN, Options{NumberBlocks: true}); err != nil {
				t.Fatalf("err (-n): %v", err)
			}
			if got := outN.String(); got != tc.wantN {
				t.Fatalf("normalise -n mismatch:\nwant=%q\ngot =%q", tc.wantN, got)
			}
		})
	}
}

// TestConvert_ThirteenColumns verifies that records with MORE than 12 columns
// are also normalised (fields.size() != 12) rather than block-split.
func TestConvert_ThirteenColumns(t *testing.T) {
	in := "chr1\t0\t50\tname\t0\t+\t0\t0\t0\t1\t50,\t0,\textra\n"
	want := "chr1\t0\t50\tname\t0\t+\n"
	var out bytes.Buffer
	if _, err := Convert(strings.NewReader(in), &out, Options{}); err != nil {
		t.Fatalf("err: %v", err)
	}
	if got := out.String(); got != want {
		t.Fatalf(">12-column record should normalise:\nwant=%q\ngot =%q", want, got)
	}
}

func TestConvert_PassThroughZeroBlocks(t *testing.T) {
	in := "chr1\t0\t50\tx\t0\t+\t0\t0\t0\t0\t\t\n"
	var out bytes.Buffer
	if _, err := Convert(strings.NewReader(in), &out, Options{}); err != nil {
		t.Fatalf("err: %v", err)
	}
	if got := out.String(); got != in {
		t.Fatalf("zero-block record should pass through; got %q", got)
	}
}

func TestConvert_BadChromStart(t *testing.T) {
	in := "chr1\tBAD\t50\tn\t0\t+\t0\t0\t0\t1\t50,\t0,\n"
	var out bytes.Buffer
	if _, err := Convert(strings.NewReader(in), &out, Options{}); err == nil {
		t.Fatalf("expected error for bad chromStart")
	}
}

func TestConvert_BadBlockCount(t *testing.T) {
	in := "chr1\t0\t50\tn\t0\t+\t0\t0\t0\tBAD\t50,\t0,\n"
	var out bytes.Buffer
	if _, err := Convert(strings.NewReader(in), &out, Options{}); err == nil {
		t.Fatalf("expected error for bad blockCount")
	}
}

func TestConvert_BadBlockSizes(t *testing.T) {
	in := "chr1\t0\t50\tn\t0\t+\t0\t0\t0\t1\tBAD,\t0,\n"
	var out bytes.Buffer
	if _, err := Convert(strings.NewReader(in), &out, Options{}); err == nil {
		t.Fatalf("expected error for bad blockSizes")
	}
}

func TestConvert_BadBlockStarts(t *testing.T) {
	in := "chr1\t0\t50\tn\t0\t+\t0\t0\t0\t1\t50,\tBAD,\n"
	var out bytes.Buffer
	if _, err := Convert(strings.NewReader(in), &out, Options{}); err == nil {
		t.Fatalf("expected error for bad blockStarts")
	}
}

func TestConvert_MismatchedBlockCountSizes(t *testing.T) {
	in := "chr1\t0\t50\tn\t0\t+\t0\t0\t0\t2\t10,\t0,\n"
	var out bytes.Buffer
	if _, err := Convert(strings.NewReader(in), &out, Options{}); err == nil {
		t.Fatalf("expected error for mismatched blockCount/blockSizes")
	}
}

func TestParseIntList(t *testing.T) {
	got, err := parseIntList("10,20,30,")
	if err != nil {
		t.Fatalf("err: %v", err)
	}
	want := []int{10, 20, 30}
	if len(got) != len(want) {
		t.Fatalf("len mismatch: want %v got %v", want, got)
	}
	for i := range want {
		if got[i] != want[i] {
			t.Fatalf("index %d: want %d got %d", i, want[i], got[i])
		}
	}
	empty, err := parseIntList("")
	if err != nil || empty != nil {
		t.Fatalf("empty case: %v %v", empty, err)
	}
	if _, err := parseIntList("a,b"); err == nil {
		t.Fatalf("expected error for non-numeric list")
	}
}
