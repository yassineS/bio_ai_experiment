package samtools

import (
	"io"
	"os"
	"testing"
)

// benchBAM points at the real chr20 fixture when present; the region/BED
// view benchmarks skip when it is absent (CI without the large fixture).
const benchBAM = "../../../../pipeline/.fixtures/realchr20/chr20.bam"
const benchBED = "../../../../pipeline/.fixtures/realchr20/chr20.bed"

func benchFixturePresent(tb testing.TB) {
	if _, err := os.Stat(benchBAM); err != nil {
		tb.Skip("realchr20 fixture absent")
	}
}

// BenchmarkViewRegionSAM measures the whole-chromosome region->SAM fast path.
func BenchmarkViewRegionSAM(b *testing.B) {
	benchFixturePresent(b)
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		if _, err := ViewFile(benchBAM, io.Discard, ViewOptions{Regions: []string{"chr20"}}, io.Discard); err != nil {
			b.Fatal(err)
		}
	}
}

// BenchmarkViewBEDSAM measures the -L BED-filtered linear scan fast path.
func BenchmarkViewBEDSAM(b *testing.B) {
	benchFixturePresent(b)
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		if _, err := ViewFile(benchBAM, io.Discard, ViewOptions{BedPath: benchBED}, io.Discard); err != nil {
			b.Fatal(err)
		}
	}
}
