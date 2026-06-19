package bcftools

import (
	"io"
	"os"
	"testing"
)

func BenchmarkViewVCF(b *testing.B) {
	vcf := "../../../../pipeline/.fixtures/medium/variants.vcf"
	if _, err := os.Stat(vcf); err != nil {
		b.Skip()
	}
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		f, _ := os.Open(vcf)
		if _, err := View(f, io.Discard, ViewOptions{NoHeader: true}); err != nil {
			b.Fatal(err)
		}
		f.Close()
	}
}
