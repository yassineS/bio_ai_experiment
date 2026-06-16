package alnio

import (
	"bytes"
	"fmt"
	"io"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// buildMultiBlockBAM encodes a BGZF-wrapped BAM with n records, large enough to
// span many BGZF blocks so the parallel decode path is genuinely exercised.
func buildMultiBlockBAM(t *testing.T, n int) []byte {
	t.Helper()
	var buf bytes.Buffer
	bw := sam.NewBAMWriter(&buf)
	hdr := &sam.Header{Refs: []sam.Reference{{Name: "chr1", Length: 100000000}}}
	if err := bw.WriteHeader(hdr); err != nil {
		t.Fatalf("WriteHeader: %v", err)
	}
	cig, err := sam.ParseCigar("100M")
	if err != nil {
		t.Fatalf("ParseCigar: %v", err)
	}
	seq := make([]byte, 100)
	for i := range seq {
		seq[i] = "ACGT"[i%4]
	}
	qual := make([]byte, 100)
	for i := range qual {
		qual[i] = 40
	}
	for i := 0; i < n; i++ {
		rec := &sam.Record{
			QName: fmt.Sprintf("read%06d", i),
			Flag:  0, RName: "chr1", Pos: int32((i*37)%99999000 + 1), MapQ: 60,
			Cigar: cig, Seq: string(seq), Qual: qual,
			RNext: "*", PNext: 0, TLen: 0,
		}
		if err := bw.Write(rec); err != nil {
			t.Fatalf("Write: %v", err)
		}
	}
	if err := bw.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	return buf.Bytes()
}

// readAllKeys decodes a stream via NewReaderThreaded at the given thread count
// and returns one key per record.
func readAllKeys(t *testing.T, bam []byte, threads int) []string {
	t.Helper()
	rd, err := NewReaderThreaded(bytes.NewReader(bam), "", threads)
	if err != nil {
		t.Fatalf("NewReaderThreaded(%d): %v", threads, err)
	}
	if rc, ok := rd.(io.Closer); ok {
		defer rc.Close()
	}
	var keys []string
	for {
		rec, err := rd.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("threads=%d read: %v", threads, err)
		}
		keys = append(keys, fmt.Sprintf("%s|%d|%s|%d", rec.QName, rec.Flag, rec.RName, rec.Pos))
	}
	return keys
}

// TestNewReaderThreaded_RecordsIdenticalAcrossThreads proves the thread-aware
// BGZF input reader yields the exact same record stream for -@ {1,2,4,8}: the
// cross-thread byte-identity invariant at the alnio layer.
func TestNewReaderThreaded_RecordsIdenticalAcrossThreads(t *testing.T) {
	bam := buildMultiBlockBAM(t, 20000)
	base := readAllKeys(t, bam, 1)
	if len(base) != 20000 {
		t.Fatalf("serial decode got %d records, want 20000", len(base))
	}
	for _, threads := range []int{2, 4, 8} {
		got := readAllKeys(t, bam, threads)
		if len(got) != len(base) {
			t.Fatalf("threads=%d record count %d != %d", threads, len(got), len(base))
		}
		for i := range base {
			if got[i] != base[i] {
				t.Fatalf("threads=%d record %d mismatch: %q vs %q", threads, i, got[i], base[i])
			}
		}
	}
}

// TestNewReaderThreaded_SAMPassthrough confirms a plain-text SAM stream is read
// unchanged regardless of thread count (nothing BGZF to parallelise).
func TestNewReaderThreaded_SAMPassthrough(t *testing.T) {
	const samText = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n" +
		"r1\t0\tchr1\t10\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
	for _, threads := range []int{1, 2, 4} {
		rd, err := NewReaderThreaded(bytes.NewReader([]byte(samText)), "", threads)
		if err != nil {
			t.Fatalf("threads=%d: %v", threads, err)
		}
		rec, err := rd.Read()
		if err != nil {
			t.Fatalf("threads=%d read: %v", threads, err)
		}
		if rec.QName != "r1" || rec.Pos != 10 {
			t.Fatalf("threads=%d unexpected record %+v", threads, rec)
		}
		if rc, ok := rd.(io.Closer); ok {
			rc.Close()
		}
	}
}
