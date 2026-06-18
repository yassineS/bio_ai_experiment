package sam

import (
	"os"
	"testing"
)

// benchBAM is the medium parity fixture; skipped if absent.
const benchBAM = "../../../pipeline/.fixtures/medium/reads.bam"

func BenchmarkBAMReadDecode(b *testing.B) {
	if _, err := os.Stat(benchBAM); err != nil {
		b.Skipf("fixture %s not present", benchBAM)
	}
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		f, err := os.Open(benchBAM)
		if err != nil {
			b.Fatal(err)
		}
		br, err := NewBAMReader(f)
		if err != nil {
			f.Close()
			b.Fatal(err)
		}
		var n int
		var flagAcc uint64
		for {
			rec, err := br.Read()
			if err != nil {
				break
			}
			flagAcc += uint64(rec.Flag) // flagstat-like: touch only Flag
			n++
		}
		br.Close()
		f.Close()
		if i == 0 {
			b.Logf("records=%d flagAcc=%d", n, flagAcc)
		}
	}
}

func BenchmarkBAMReadInto(b *testing.B) {
	if _, err := os.Stat(benchBAM); err != nil {
		b.Skipf("fixture %s not present", benchBAM)
	}
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		f, err := os.Open(benchBAM)
		if err != nil {
			b.Fatal(err)
		}
		br, err := NewBAMReader(f)
		if err != nil {
			f.Close()
			b.Fatal(err)
		}
		var n int
		var flagAcc uint64
		var rec Record
		for {
			if err := br.ReadInto(&rec); err != nil {
				break
			}
			flagAcc += uint64(rec.Flag)
			n++
		}
		br.Close()
		f.Close()
		if i == 0 {
			b.Logf("records=%d flagAcc=%d", n, flagAcc)
		}
	}
}
