package cram

import (
	"bytes"
	"fmt"
	"reflect"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// buildMultiContainerCRAM encodes a reference-free CRAM with n records. The
// writer caps each container at defaultRecordsPerSlice records, so an n well
// above that cap spans several containers — exactly the structure the parallel
// driver decodes across its worker pool. Mapped reads on two contigs plus a few
// unmapped reads give the decode a realistic mix.
func buildMultiContainerCRAM(t *testing.T, n int) ([]byte, *sam.Header) {
	t.Helper()
	h := writerTestHeader()
	cig, err := sam.ParseCigar("10M")
	if err != nil {
		t.Fatalf("ParseCigar: %v", err)
	}
	recs := make([]*sam.Record, 0, n)
	for i := 0; i < n; i++ {
		rec := &sam.Record{
			QName: fmt.Sprintf("read%07d", i),
			Flag:  0,
			RName: "chr1",
			Pos:   int64(i%90000 + 1),
			MapQ:  40,
			Cigar: cig,
			Seq:   "ACGTACGTAC",
			Qual:  []byte{30, 31, 32, 33, 34, 35, 36, 37, 38, 39},
			RNext: "*",
		}
		switch i % 5 {
		case 1:
			rec.RName = "chr2"
			rec.Pos = int64(i%40000 + 1)
		case 4:
			// An unmapped read: no placement, no CIGAR.
			rec.Flag = 0x4
			rec.RName = "*"
			rec.Pos = 0
			rec.MapQ = 0
			rec.Cigar = nil
		}
		recs = append(recs, rec)
	}
	var buf bytes.Buffer
	if err := WriteCRAM(&buf, h, recs); err != nil {
		t.Fatalf("WriteCRAM: %v", err)
	}
	return buf.Bytes(), h
}

// decodeAllThreaded decodes every record of a CRAM byte stream with the given
// thread count, returning the reconstructed records.
func decodeAllThreaded(t *testing.T, data []byte, threads int) []*sam.Record {
	t.Helper()
	rr, err := NewRecordReader(bytes.NewReader(data))
	if err != nil {
		t.Fatalf("NewRecordReader(threads=%d): %v", threads, err)
	}
	defer rr.Close()
	rr.SetThreads(threads)
	recs, err := rr.ReadAll()
	if err != nil {
		t.Fatalf("ReadAll(threads=%d): %v", threads, err)
	}
	return recs
}

// TestParallelCRAMRecordsIdentical asserts that container/slice-parallel CRAM
// decode yields exactly the same records, in the same order, for any thread
// count. The single-threaded decode is the oracle; -@ {2,4,8} must reproduce it
// record-for-record. This is the non-negotiable cross-thread byte-identity
// invariant at the record level.
func TestParallelCRAMRecordsIdentical(t *testing.T) {
	data, _ := buildMultiContainerCRAM(t, 25000)
	want := decodeAllThreaded(t, data, 1)
	if len(want) != 25000 {
		t.Fatalf("oracle decoded %d records, want 25000", len(want))
	}
	for _, threads := range []int{2, 4, 8} {
		got := decodeAllThreaded(t, data, threads)
		if len(got) != len(want) {
			t.Fatalf("threads=%d: decoded %d records, want %d", threads, len(got), len(want))
		}
		for i := range want {
			if !reflect.DeepEqual(got[i], want[i]) {
				t.Fatalf("threads=%d: record %d differs from single-threaded decode:\n got=%+v\nwant=%+v",
					threads, i, got[i], want[i])
			}
		}
	}
}

// TestParallelCRAMSAMTextIdentical asserts that the full SAM text WriteSAM
// produces is byte-for-byte identical for every thread count. This is the
// invariant at the output-byte level: header plus every record line must match
// the single-threaded rendering exactly.
func TestParallelCRAMSAMTextIdentical(t *testing.T) {
	data, _ := buildMultiContainerCRAM(t, 25000)

	render := func(threads int) []byte {
		rr, err := NewRecordReader(bytes.NewReader(data))
		if err != nil {
			t.Fatalf("NewRecordReader(threads=%d): %v", threads, err)
		}
		defer rr.Close()
		rr.SetThreads(threads)
		var out bytes.Buffer
		if err := rr.WriteSAM(&out); err != nil {
			t.Fatalf("WriteSAM(threads=%d): %v", threads, err)
		}
		return out.Bytes()
	}

	want := render(1)
	for _, threads := range []int{2, 4, 8} {
		got := render(threads)
		if !bytes.Equal(got, want) {
			t.Fatalf("threads=%d: SAM text differs from single-threaded output (len got=%d want=%d)",
				threads, len(got), len(want))
		}
	}
}

// TestParallelCRAMSingleContainer checks the small-file case: a CRAM with a
// single data container still decodes identically under the parallel path,
// where the worker pool sees exactly one job.
func TestParallelCRAMSingleContainer(t *testing.T) {
	data, _ := buildMultiContainerCRAM(t, 50)
	want := decodeAllThreaded(t, data, 1)
	got := decodeAllThreaded(t, data, 4)
	if len(got) != len(want) {
		t.Fatalf("decoded %d records, want %d", len(got), len(want))
	}
	for i := range want {
		if !reflect.DeepEqual(got[i], want[i]) {
			t.Fatalf("record %d differs under parallel decode", i)
		}
	}
}

// TestParallelCRAMSetThreadsAfterReadNoop verifies SetThreads is a no-op once
// decoding has started, so a late call cannot corrupt an in-flight read.
func TestParallelCRAMSetThreadsAfterReadNoop(t *testing.T) {
	data, _ := buildMultiContainerCRAM(t, 50)
	rr, err := NewRecordReader(bytes.NewReader(data))
	if err != nil {
		t.Fatalf("NewRecordReader: %v", err)
	}
	defer rr.Close()
	if _, err := rr.Read(); err != nil {
		t.Fatalf("first Read: %v", err)
	}
	rr.SetThreads(8) // too late: must not engage the parallel driver.
	if rr.par != nil || rr.threads != 0 {
		t.Fatalf("SetThreads after first Read engaged parallel decode (par=%v threads=%d)", rr.par, rr.threads)
	}
	// The rest of the stream still decodes cleanly on the single-threaded path.
	if _, err := rr.ReadAll(); err != nil {
		t.Fatalf("ReadAll after late SetThreads: %v", err)
	}
}

// TestParallelCRAMEarlyClose closes a parallel reader before draining it,
// exercising the driver teardown path (the consumer abandons the stream). Run
// under -race this catches a goroutine left blocked on a channel.
func TestParallelCRAMEarlyClose(t *testing.T) {
	data, _ := buildMultiContainerCRAM(t, 25000)
	rr, err := NewRecordReader(bytes.NewReader(data))
	if err != nil {
		t.Fatalf("NewRecordReader: %v", err)
	}
	rr.SetThreads(4)
	// Read only a handful of records, then close mid-stream.
	for i := 0; i < 5; i++ {
		if _, err := rr.Read(); err != nil {
			t.Fatalf("Read %d: %v", i, err)
		}
	}
	if err := rr.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
}

// BenchmarkParallelCRAMDecode measures decode throughput across thread counts
// on a multi-container CRAM. It reports the speedup of the parallel path over
// the single-threaded baseline.
func BenchmarkParallelCRAMDecode(b *testing.B) {
	h := writerTestHeader()
	cig, _ := sam.ParseCigar("10M")
	const n = 50000
	recs := make([]*sam.Record, 0, n)
	for i := 0; i < n; i++ {
		recs = append(recs, &sam.Record{
			QName: fmt.Sprintf("read%07d", i),
			RName: "chr1", Pos: int64(i%90000 + 1), MapQ: 40, Cigar: cig,
			Seq: "ACGTACGTAC", Qual: []byte{30, 31, 32, 33, 34, 35, 36, 37, 38, 39}, RNext: "*",
		})
	}
	var buf bytes.Buffer
	if err := WriteCRAM(&buf, h, recs); err != nil {
		b.Fatalf("WriteCRAM: %v", err)
	}
	data := buf.Bytes()

	for _, threads := range []int{1, 2, 4, 8} {
		b.Run(fmt.Sprintf("threads=%d", threads), func(b *testing.B) {
			b.ReportAllocs()
			for i := 0; i < b.N; i++ {
				rr, err := NewRecordReader(bytes.NewReader(data))
				if err != nil {
					b.Fatalf("NewRecordReader: %v", err)
				}
				rr.SetThreads(threads)
				if _, err := rr.ReadAll(); err != nil {
					b.Fatalf("ReadAll: %v", err)
				}
				rr.Close()
			}
		})
	}
}
