package samtools

import (
	"io"
	"os"
	"testing"
)

const benchStatsBAM = "../../../../pipeline/.fixtures/medium/reads.bam"

func BenchmarkStats(b *testing.B) {
	if _, err := os.Stat(benchStatsBAM); err != nil {
		b.Skipf("fixture absent")
	}
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		f, _ := os.Open(benchStatsBAM)
		if err := Stats(f, io.Discard, StatsOptions{}); err != nil {
			b.Fatal(err)
		}
		f.Close()
	}
}
