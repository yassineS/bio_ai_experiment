package bcftools

import (
	"io"
	"os"
	"testing"
)

func BenchmarkQueryVCF(b *testing.B) {
	vcf := "../../../../pipeline/.fixtures/medium/variants.vcf"
	if _, err := os.Stat(vcf); err != nil {
		b.Skip()
	}
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		f, _ := os.Open(vcf)
		if _, err := Query(f, io.Discard, QueryOptions{Format: `%CHROM\t%POS\t%REF\t%ALT\n`}); err != nil {
			b.Fatal(err)
		}
		f.Close()
	}
}
