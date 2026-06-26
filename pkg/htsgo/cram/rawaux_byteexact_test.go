package cram

import (
	"bytes"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// decodeCRAMToBAM decodes the reference-backed CRAM in cramBytes to a BAM byte
// stream, feeding the supplied reference. rawAux selects the memory-lean
// passthrough decode mode (rec.RawAux) when true, or the eager []sam.Aux path
// when false. The two must produce byte-identical BAM output — that equality is
// the local pre-bioval proof that the raw-aux passthrough emits exactly what the
// eager path does.
func decodeCRAMToBAM(t *testing.T, cramBytes []byte, ref *stubReference, hdr *sam.Header, rawAux bool) []byte {
	return decodeCRAMToBAMThreads(t, cramBytes, ref, hdr, rawAux, 0)
}

func decodeCRAMToBAMThreads(t *testing.T, cramBytes []byte, ref *stubReference, hdr *sam.Header, rawAux bool, threads int) []byte {
	t.Helper()
	rr, err := NewRecordReader(bytes.NewReader(cramBytes))
	if err != nil {
		t.Fatalf("NewRecordReader: %v", err)
	}
	rr.SetReference(ref)
	if threads >= 2 {
		rr.SetThreads(threads)
	}
	if rawAux {
		rr.SetRawAuxBAMSink(true)
	}
	var out bytes.Buffer
	bw := sam.NewBAMWriter(&out)
	if err := bw.WriteHeader(hdr); err != nil {
		t.Fatalf("WriteHeader: %v", err)
	}
	for {
		rec, err := rr.Read()
		if err != nil {
			break
		}
		if err := bw.Write(rec); err != nil {
			t.Fatalf("BAM Write: %v", err)
		}
	}
	if err := bw.Close(); err != nil {
		t.Fatalf("BAM Close: %v", err)
	}
	return out.Bytes()
}

// TestRawAuxBAMSinkMatchesDecode is the byte-exact gate for the CRAM→BAM
// memory-lean raw-aux passthrough (task #43): a reference-backed CRAM is decoded
// twice — once on the eager []sam.Aux path (mode OFF) and once on the RawAux
// passthrough (mode ON) — to BAM, and the two BAM byte streams must be
// identical. The eager output is already byte-identical to upstream
// `samtools view -b -T ref` per the #40 baseline, so equality here proves the
// passthrough is byte-exact too. The fixture covers records that need MD/NM
// regeneration (mismatched/indel reads against the reference) and records with a
// trailing data-series RG tag, so the MD<NM<RG raw splice is exercised.
func TestRawAuxBAMSinkMatchesDecode(t *testing.T) {
	h := writerTestHeader() // chr1/chr2 + @RG ID:rg1

	ref := []byte(strings.Repeat("ACGTACGTGGCCATGCTAGC", 20)) // 400 bp
	refMap := map[string][]byte{"chr1": ref}

	mut := func(pos1, n int, edits map[int]byte) string {
		b := append([]byte(nil), ref[pos1-1:pos1-1+n]...)
		for off, base := range edits {
			b[off] = base
		}
		return string(b)
	}
	// withAux returns rec with the given aux tags attached (RG last, mimicking
	// the data-series RG position upstream emits). A perfect-match read yields
	// MD:Z:<len> NM:i:0 on decode; a mismatched/indel read yields a non-trivial
	// MD string and NM > 0 — both regenerated because the decode has a reference.
	withAux := func(rec *sam.Record, aux ...sam.Aux) *sam.Record {
		rec.Aux = aux
		return rec
	}

	in := []*sam.Record{
		// Perfect match, trailing RG only → MD:Z:30 NM:i:0 spliced before RG.
		withAux(mkRec("r_exact", "chr1", 50, "30M", mut(50, 30, nil)),
			sam.Aux{Tag: "RG", Type: 'Z', Value: "rg1"}),
		// Substitutions, a non-RG dictionary tag then a trailing RG → MD/NM land
		// between the dictionary tag and the RG.
		withAux(mkRec("r_subs", "chr1", 100, "40M", mut(100, 40, map[int]byte{5: 'T', 22: 'A', 31: 'C'})),
			sam.Aux{Tag: "NH", Type: 'i', Value: int64(1)},
			sam.Aux{Tag: "RG", Type: 'Z', Value: "rg1"}),
		// Deletion + subs, no aux at all → MD/NM appended at the end.
		mkRec("r_indel", "chr1", 200, "20M5D20M",
			mut(200, 45, map[int]byte{3: 'G', 38: 'T'})[:20]+mut(225, 20, nil)),
		// A record carrying a B-array aux and an 'A' aux (exercises the wider aux
		// serialisers in the raw builder) plus trailing RG.
		withAux(mkRec("r_arr", "chr1", 300, "25M", mut(300, 25, map[int]byte{10: 'A'})),
			sam.Aux{Tag: "BC", Type: 'A', Value: "P"},
			sam.Aux{Tag: "ZB", Type: 'B', ArrayType: 'S', ArrayValues: []interface{}{int64(1), int64(258), int64(3)}},
			sam.Aux{Tag: "RG", Type: 'Z', Value: "rg1"}),
		// An unmapped read (no MD/NM regen) with an aux tag — must pass through
		// unchanged on both paths.
		func() *sam.Record {
			r := mkRec("r_unmap", "*", 0, "*", "ACGTACGTAC")
			r.Flag = sam.FlagUnmapped
			r.MapQ = 0
			return withAux(r, sam.Aux{Tag: "RG", Type: 'Z', Value: "rg1"})
		}(),
	}

	var cramBuf bytes.Buffer
	rw, err := NewRecordWriterOpts(&cramBuf, h, WriterOptions{Reference: refMap})
	if err != nil {
		t.Fatalf("NewRecordWriterOpts(reference): %v", err)
	}
	for _, rec := range in {
		if err := rw.Write(rec); err != nil {
			t.Fatalf("Write: %v", err)
		}
	}
	if err := rw.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	stub := &stubReference{name: "chr1", seq: string(ref)}
	eager := decodeCRAMToBAM(t, cramBuf.Bytes(), stub, h, false)
	raw := decodeCRAMToBAM(t, cramBuf.Bytes(), stub, h, true)

	if !bytes.Equal(eager, raw) {
		t.Fatalf("CRAM→BAM mode-OFF (eager, %d B) and mode-ON (RawAux, %d B) BAM output differ",
			len(eager), len(raw))
	}

	// The parallel decoder threads the same rawAuxBAMSink flag through its
	// workers, so its mode-ON output must also equal the eager output.
	rawParallel := decodeCRAMToBAMThreads(t, cramBuf.Bytes(), stub, h, true, 4)
	if !bytes.Equal(eager, rawParallel) {
		t.Fatalf("CRAM→BAM mode-ON parallel (%d B) differs from eager (%d B)",
			len(rawParallel), len(eager))
	}

	// Sanity: the eager decode must actually have produced MD/NM and RG aux so
	// the equality above is not vacuous (e.g. both empty).
	rr, err := NewRecordReader(bytes.NewReader(cramBuf.Bytes()))
	if err != nil {
		t.Fatalf("NewRecordReader (verify): %v", err)
	}
	rr.SetReference(stub)
	recs, err := rr.ReadAll()
	if err != nil {
		t.Fatalf("ReadAll (verify): %v", err)
	}
	first := recs[0]
	if _, ok := first.GetAux("MD"); !ok {
		t.Errorf("expected regenerated MD on r_exact, got none (aux=%v)", first.Aux)
	}
	if _, ok := first.GetAux("NM"); !ok {
		t.Errorf("expected regenerated NM on r_exact, got none")
	}
	if _, ok := first.GetAux("RG"); !ok {
		t.Errorf("expected trailing RG on r_exact, got none")
	}
}

// TestRawAuxRoundTripsThroughDecode asserts the raw-aux byte block a record
// carries on the passthrough path round-trips losslessly:
// AppendBAMAux(decodeBAMAux(RawAux)) == RawAux for every decoded record. This
// pins the raw bytes to exactly the canonical BAM aux encoding the writer would
// emit for the equivalent []sam.Aux.
func TestRawAuxRoundTripsThroughDecode(t *testing.T) {
	h := writerTestHeader()
	ref := []byte(strings.Repeat("ACGTACGTGGCCATGCTAGC", 20))
	refMap := map[string][]byte{"chr1": ref}
	mut := func(pos1, n int, edits map[int]byte) string {
		b := append([]byte(nil), ref[pos1-1:pos1-1+n]...)
		for off, base := range edits {
			b[off] = base
		}
		return string(b)
	}
	in := []*sam.Record{
		func() *sam.Record {
			r := mkRec("a", "chr1", 50, "30M", mut(50, 30, map[int]byte{4: 'T'}))
			r.Aux = []sam.Aux{
				{Tag: "NH", Type: 'i', Value: int64(1)},
				{Tag: "ZB", Type: 'B', ArrayType: 'i', ArrayValues: []interface{}{int64(-1), int64(70000)}},
				{Tag: "RG", Type: 'Z', Value: "rg1"},
			}
			return r
		}(),
	}
	var cramBuf bytes.Buffer
	rw, err := NewRecordWriterOpts(&cramBuf, h, WriterOptions{Reference: refMap})
	if err != nil {
		t.Fatalf("NewRecordWriterOpts: %v", err)
	}
	for _, rec := range in {
		if err := rw.Write(rec); err != nil {
			t.Fatalf("Write: %v", err)
		}
	}
	if err := rw.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}

	rr, err := NewRecordReader(bytes.NewReader(cramBuf.Bytes()))
	if err != nil {
		t.Fatalf("NewRecordReader: %v", err)
	}
	rr.SetReference(&stubReference{name: "chr1", seq: string(ref)})
	rr.SetRawAuxBAMSink(true)
	recs, err := rr.ReadAll()
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	checked := 0
	for i, rec := range recs {
		if rec.RawAux == nil {
			continue
		}
		decoded, derr := sam.DecodeBAMAux(rec.RawAux)
		if derr != nil {
			t.Fatalf("record %d decodeBAMAux(RawAux): %v", i, derr)
		}
		var reenc []byte
		for _, a := range decoded {
			reenc, err = sam.AppendBAMAux(reenc, a)
			if err != nil {
				t.Fatalf("record %d AppendBAMAux: %v", i, err)
			}
		}
		if !bytes.Equal(reenc, rec.RawAux) {
			t.Errorf("record %d raw-aux round-trip mismatch:\n got %x\nwant %x", i, reenc, rec.RawAux)
		}
		checked++
	}
	if checked == 0 {
		t.Fatal("no records carried RawAux; raw-aux mode did not engage")
	}
}
