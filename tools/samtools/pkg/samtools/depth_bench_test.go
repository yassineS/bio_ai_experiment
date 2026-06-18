package samtools

import (
	"io"
	"os"
	"testing"
)

func BenchmarkDepth(b *testing.B) {
	bam := "../../../../pipeline/.fixtures/medium/reads.bam"
	if _, err := os.Stat(bam); err != nil {
		b.Skip()
	}
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		f, _ := os.Open(bam)
		if err := Depth([]io.Reader{f}, io.Discard, DepthOptions{}); err != nil {
			b.Fatal(err)
		}
		f.Close()
	}
}
