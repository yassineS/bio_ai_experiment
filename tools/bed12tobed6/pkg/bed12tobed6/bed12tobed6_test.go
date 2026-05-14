package bed12tobed6

import (
	"bytes"
	"strings"
	"testing"
)

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

func TestConvert_PassThroughShortRecord(t *testing.T) {
	in := "chr1\t0\t50\tbed6\t1\t+\n"
	var out bytes.Buffer
	if _, err := Convert(strings.NewReader(in), &out, Options{}); err != nil {
		t.Fatalf("err: %v", err)
	}
	if got := out.String(); got != in {
		t.Fatalf("short record should pass through unchanged; got %q", got)
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
