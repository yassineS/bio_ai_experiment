package samtools

import (
	"bytes"
	"fmt"
	"os"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/cram"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// buildMultiContainerCRAMBytes encodes a reference-free CRAM with n records.
// The CRAM writer caps each container at its per-slice record limit, so an n
// well above that cap spans several containers — the structure that exercises
// the container/slice-parallel decode wired to -@. The records mix two contigs
// and a fraction of unmapped reads so the decode covers placed and unplaced
// paths.
func buildMultiContainerCRAMBytes(t *testing.T, n int) []byte {
	t.Helper()
	text := "@HD\tVN:1.6\tSO:unsorted\n" +
		"@SQ\tSN:chr1\tLN:100000\n" +
		"@SQ\tSN:chr2\tLN:50000\n" +
		"@RG\tID:rg1\tSM:sample1\n"
	h, err := sam.ParseHeaderText(text)
	if err != nil {
		t.Fatalf("ParseHeaderText: %v", err)
	}
	cig, err := sam.ParseCigar("10M")
	if err != nil {
		t.Fatalf("ParseCigar: %v", err)
	}
	recs := make([]*sam.Record, 0, n)
	for i := 0; i < n; i++ {
		rec := &sam.Record{
			QName: fmt.Sprintf("read%07d", i),
			RName: "chr1", Pos: int64(i%90000 + 1), MapQ: 40, Cigar: cig,
			Seq: "ACGTACGTAC", Qual: []byte{30, 31, 32, 33, 34, 35, 36, 37, 38, 39}, RNext: "*",
		}
		switch i % 5 {
		case 1:
			rec.RName = "chr2"
			rec.Pos = int64(i%40000 + 1)
		case 4:
			rec.Flag = 0x4
			rec.RName = "*"
			rec.Pos = 0
			rec.MapQ = 0
			rec.Cigar = nil
		}
		recs = append(recs, rec)
	}
	var buf bytes.Buffer
	if err := cram.WriteCRAM(&buf, h, recs); err != nil {
		t.Fatalf("WriteCRAM: %v", err)
	}
	return buf.Bytes()
}

// TestViewCRAMThreadsByteIdentical asserts that `samtools view` of a
// multi-container CRAM produces byte-identical SAM text for every -@ thread
// count. The single-threaded output is the oracle; -@ {2,4,8} must reproduce it
// exactly. Only decode throughput changes with threads, never the bytes.
func TestViewCRAMThreadsByteIdentical(t *testing.T) {
	data := buildMultiContainerCRAMBytes(t, 25000)

	run := func(threads int) []byte {
		var out bytes.Buffer
		if _, err := View(bytes.NewReader(data), &out, ViewOptions{
			WithHeader: true,
			Threads:    threads,
		}); err != nil {
			t.Fatalf("View(threads=%d): %v", threads, err)
		}
		return out.Bytes()
	}

	want := run(1)
	if len(want) == 0 {
		t.Fatal("single-threaded view produced no output")
	}
	for _, threads := range []int{2, 4, 8} {
		got := run(threads)
		if !bytes.Equal(got, want) {
			t.Fatalf("view -@ %d SAM text differs from single-threaded output (len got=%d want=%d)",
				threads, len(got), len(want))
		}
	}
}

// TestFlagstatCRAMThreadsByteIdentical asserts flagstat's text output over a
// multi-container CRAM is identical for every thread count. flagstat reaches
// CRAM input through the same alnio.NewReaderThreaded path as view.
func TestFlagstatCRAMThreadsByteIdentical(t *testing.T) {
	data := buildMultiContainerCRAMBytes(t, 25000)

	run := func(threads int) []byte {
		var out bytes.Buffer
		if err := FlagstatThreaded(bytes.NewReader(data), &out, threads); err != nil {
			t.Fatalf("FlagstatThreaded(threads=%d): %v", threads, err)
		}
		return out.Bytes()
	}

	want := run(1)
	for _, threads := range []int{2, 4, 8} {
		got := run(threads)
		if !bytes.Equal(got, want) {
			t.Fatalf("flagstat -@ %d output differs from single-threaded:\n got=%s\nwant=%s",
				threads, got, want)
		}
	}
}

// TestViewVendoredCRAMThreadsIdentical runs the vendored upstream CRAM fixture
// (a real samtools test file, not synthetic data) through view at several
// thread counts and asserts byte-identical SAM text. It complements the
// synthetic multi-container test by covering a real on-disk CRAM.
func TestViewVendoredCRAMThreadsIdentical(t *testing.T) {
	path := openCRAMFixture(t)
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("read fixture: %v", err)
	}
	run := func(threads int) []byte {
		var out bytes.Buffer
		if _, err := View(bytes.NewReader(data), &out, ViewOptions{WithHeader: true, Threads: threads}); err != nil {
			t.Fatalf("View(threads=%d): %v", threads, err)
		}
		return out.Bytes()
	}
	want := run(1)
	for _, threads := range []int{2, 4, 8} {
		if !bytes.Equal(run(threads), want) {
			t.Fatalf("vendored CRAM view -@ %d differs from single-threaded output", threads)
		}
	}
}
