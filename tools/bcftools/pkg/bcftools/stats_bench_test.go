package bcftools

import (
	"bytes"
	"io"
	"os"
	"testing"
)

func BenchmarkStats(b *testing.B) {
	data, err := os.ReadFile("/home/user/bio_ai_experiment/pipeline/.fixtures/medium/variants.vcf")
	if err != nil {
		b.Skip("fixture not available")
	}
	b.SetBytes(int64(len(data)))
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, err := Stats(bytes.NewReader(data), io.Discard, StatsOptions{InputFile: "-"}); err != nil {
			b.Fatal(err)
		}
	}
}

func BenchmarkStatsSamples(b *testing.B) {
	data, err := os.ReadFile("/home/user/bio_ai_experiment/pipeline/.fixtures/medium/variants.vcf")
	if err != nil {
		b.Skip("fixture not available")
	}
	b.SetBytes(int64(len(data)))
	b.ReportAllocs()
	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		if _, err := Stats(bytes.NewReader(data), io.Discard, StatsOptions{InputFile: "-", Samples: []string{"-"}}); err != nil {
			b.Fatal(err)
		}
	}
}
