package cram

import (
	"bytes"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// fuzzRecordsFromBytes derives a deterministic batch of sam.Records from
// arbitrary fuzz bytes. The construction is total — every input yields a
// valid (if odd) record set — so the fuzz target exercises the writer's
// encoding paths rather than its argument validation.
func fuzzRecordsFromBytes(data []byte) []*sam.Record {
	if len(data) == 0 {
		return nil
	}
	refNames := []string{"chr1", "chr2"}
	bases := "ACGT"
	var recs []*sam.Record
	i := 0
	for i < len(data) && len(recs) < 64 {
		// Each record consumes a small, fixed window of fuzz bytes.
		b := func(off int) byte {
			if i+off < len(data) {
				return data[i+off]
			}
			return 0
		}
		ctrl := b(0)
		unmapped := ctrl&1 != 0
		seqLen := int(b(1)%16) + 1
		var seq bytes.Buffer
		for j := 0; j < seqLen; j++ {
			seq.WriteByte(bases[int(b(2+j%4))%4])
		}
		rec := &sam.Record{
			QName: "fuzz" + string(rune('0'+len(recs)%10)),
			Seq:   seq.String(),
		}
		qual := make([]byte, seqLen)
		for j := range qual {
			qual[j] = b(1) % 60
		}
		rec.Qual = qual
		if unmapped {
			rec.Flag = sam.FlagUnmapped
			rec.RName = "*"
			rec.RNext = "*"
		} else {
			rec.RName = refNames[int(b(0)>>1)%2]
			rec.Pos = int64(b(1))%1000 + 1
			rec.MapQ = b(0) % 60
			// A plain full-length match CIGAR keeps the record encodable.
			cig, _ := sam.ParseCigar(itoa(seqLen) + "M")
			rec.Cigar = cig
			rec.RNext = "*"
		}
		recs = append(recs, rec)
		i += 8
	}
	return recs
}

// itoa is a tiny integer-to-string helper for the fuzz CIGAR builder.
func itoa(n int) string {
	if n == 0 {
		return "0"
	}
	var buf [12]byte
	pos := len(buf)
	for n > 0 {
		pos--
		buf[pos] = byte('0' + n%10)
		n /= 10
	}
	return string(buf[pos:])
}

// fuzzWriteRoundTrip writes records as the given CRAM version and reads
// them back, asserting an error-free write always yields a file that
// re-reads to the same records. It is the shared body of the v3.0 and
// v3.1 fuzz targets.
func fuzzWriteRoundTrip(t *testing.T, records []*sam.Record, version Version) {
	h := writerTestHeader()
	var buf bytes.Buffer
	if err := WriteCRAMVersion(&buf, h, records, version); err != nil {
		// A record the simple writer rejects is a clean error, not a
		// panic; nothing more to check.
		return
	}
	rr, err := NewRecordReader(bytes.NewReader(buf.Bytes()))
	if err != nil {
		t.Fatalf("writer produced an unreadable %s file: %v", version, err)
	}
	out, err := rr.ReadAll()
	if err != nil {
		t.Fatalf("writer produced a %s file that fails to decode: %v", version, err)
	}
	if len(out) != len(records) {
		t.Fatalf("round-trip record count = %d, want %d", len(out), len(records))
	}
	// A written CRAM must not just re-read — it must re-read to the
	// same records. Check every field of every record.
	for i := range records {
		assertRecordEqual(t, i, out[i], records[i])
	}
}

// fuzzWriterSeeds adds the shared corpus seeds to a fuzz target.
func fuzzWriterSeeds(f *testing.F) {
	f.Add([]byte{})
	f.Add([]byte{0x00, 0x05, 0x41, 0x42, 0x43, 0x44, 0x10, 0x20})
	f.Add([]byte{0x01, 0x08, 0x47, 0x47, 0x47, 0x47, 0x05, 0x05})
	f.Add(bytes.Repeat([]byte{0x02, 0x0a, 0x54, 0x41, 0x43, 0x47, 0x1f, 0x3c}, 8))
}

// FuzzRecordWriter drives the CRAM v3.0 writer with record batches
// derived from arbitrary bytes, then reads the result back. The writer
// must never panic, and any file it produces without error must re-read
// without error: a written CRAM is always a valid CRAM.
func FuzzRecordWriter(f *testing.F) {
	fuzzWriterSeeds(f)
	f.Fuzz(func(t *testing.T, data []byte) {
		fuzzWriteRoundTrip(t, fuzzRecordsFromBytes(data), VersionV30)
	})
}

// FuzzRecordWriterV31 is the v3.1 counterpart of FuzzRecordWriter: it
// drives the rANS 4x16-capable write path with the same derived record
// batches, holding it to the same "a written CRAM is always a valid
// CRAM" invariant.
func FuzzRecordWriterV31(f *testing.F) {
	fuzzWriterSeeds(f)
	f.Fuzz(func(t *testing.T, data []byte) {
		fuzzWriteRoundTrip(t, fuzzRecordsFromBytes(data), VersionV31)
	})
}
