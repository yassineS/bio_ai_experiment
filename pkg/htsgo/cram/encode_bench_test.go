package cram

import (
	"io"
	"os"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

func loadBenchRecords(b *testing.B) (*sam.Header, []*sam.Record, map[string][]byte) {
	bamPath := "../../../pipeline/.fixtures/medium/reads.bam"
	faPath := "../../../pipeline/.fixtures/medium/ref.fa"
	if _, err := os.Stat(bamPath); err != nil {
		b.Skip("fixture absent")
	}
	f, err := os.Open(bamPath)
	if err != nil {
		b.Fatal(err)
	}
	defer f.Close()
	br, err := sam.NewBAMReader(f)
	if err != nil {
		b.Fatal(err)
	}
	hdr := br.Header()
	var recs []*sam.Record
	for {
		r, err := br.Read()
		if err != nil {
			break
		}
		recs = append(recs, r)
	}
	ra, err := fasta.OpenRandomAccess(faPath)
	if err != nil {
		b.Fatal(err)
	}
	defer ra.Close()
	ref := make(map[string][]byte)
	for _, sq := range hdr.Refs {
		n := ra.Length(sq.Name)
		if n > 0 {
			seq, _ := ra.Fetch(sq.Name, 0, n)
			ref[sq.Name] = seq
		}
	}
	return hdr, recs, ref
}

func BenchmarkCRAMEncode(b *testing.B) {
	hdr, recs, ref := loadBenchRecords(b)
	b.ResetTimer()
	b.ReportAllocs()
	for i := 0; i < b.N; i++ {
		rw, err := NewRecordWriterOpts(io.Discard, hdr, WriterOptions{Reference: ref})
		if err != nil {
			b.Fatal(err)
		}
		for _, r := range recs {
			if err := rw.Write(r); err != nil {
				b.Fatal(err)
			}
		}
		if err := rw.Close(); err != nil {
			b.Fatal(err)
		}
	}
}
