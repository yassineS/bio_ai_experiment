package bcftools

import (
	"bytes"
	"fmt"
	"io"
	"strings"
	"testing"
)

// BenchmarkMcall exercises the multiallelic caller over a synthetic mpileup-PL
// stream. It is self-contained (no external fixture) so it runs in CI and, via
// -benchmem, guards the call hot path against allocation regressions — notably
// the per-record VCF reuse that vcfRecordSource.Read relies on.
func BenchmarkMcall(b *testing.B) {
	const n = 4000
	recs := make([]string, n)
	for i := 0; i < n; i++ {
		recs[i] = fmt.Sprintf(
			"17\t%d\t.\tA\t<*>\t0\t.\tDP=5;I16=5,0,0,0,202,8170,0,0,145,4205,0,0,107,2459,0,0;QS=1,0;MQ0F=0\tPL:AD\t0,15,100:5,0",
			i+1)
	}
	in := []byte(mcallVCFHeader + strings.Join(recs, "\n") + "\n")
	b.SetBytes(int64(len(in)))
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, err := Call(bytes.NewReader(in), io.Discard, CallOptions{Model: CallModelMultiallelic}); err != nil {
			b.Fatal(err)
		}
	}
}
