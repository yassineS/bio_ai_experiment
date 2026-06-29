package codec

import (
	"fmt"
	"testing"
)

// benchNameBlock builds a method-8 name-token block from n synthetic
// Illumina-style read names — the shape that dominates real CRAM RN
// series — so BenchmarkNameTokDecode measures decode allocations against
// a realistic workload.
func benchNameBlock(tb testing.TB, n int) []byte {
	tb.Helper()
	var raw []byte
	for i := 0; i < n; i++ {
		name := fmt.Sprintf("HISEQ1:11:H8GV6ADXX:1:1101:%d:%d", 1200+i, 2000+i*3)
		raw = append(raw, name...)
		raw = append(raw, 0)
	}
	enc, err := NameTokEncode(raw, 1)
	if err != nil {
		tb.Fatalf("encode: %v", err)
	}
	return enc
}

// BenchmarkNameTokDecode reports allocations per decode of a name-token
// block. It is the direct measure of the codec churn the buffer pool
// targets.
func BenchmarkNameTokDecode(b *testing.B) {
	enc := benchNameBlock(b, 5000)
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		out, err := NameTokDecode(enc)
		if err != nil {
			b.Fatalf("decode: %v", err)
		}
		if len(out) == 0 {
			b.Fatal("empty decode")
		}
	}
}
