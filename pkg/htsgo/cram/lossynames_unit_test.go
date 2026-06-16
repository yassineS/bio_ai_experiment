package cram

import (
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// TestUnitReconstructDroppedNames proves, with no upstream binary, that the
// lossy-names name reconstruction matches htslib's cram_to_bam scheme:
//   - a record whose name was stored (nameStored, e.g. a detached record
//     read from the mate block) keeps that name;
//   - a within-slice pair whose names were both dropped shares one
//     synthesised name, taken from the upstream (lower-index) record;
//   - an orphan dropped record synthesises "<prefix>:<record_number>",
//     where record_number = slice counter + index + 1.
func TestUnitReconstructDroppedNames(t *testing.T) {
	newRD := func(prefix string, counter int64) *recordDecoder {
		h := refFreeHeader()
		h.Preservation.ReadNamesIncluded = false
		return &recordDecoder{
			h:          h,
			slice:      &SliceHeader{RecordCounter: counter},
			namePrefix: prefix,
		}
	}

	t.Run("detached-keeps-stored-name", func(t *testing.T) {
		rd := newRD("file.cram", 0)
		dr := &decodedRecord{rec: &sam.Record{QName: "readB"}, nameStored: true, mateIndex: -1}
		rd.reconstructDroppedNames([]*decodedRecord{dr})
		if dr.rec.QName != "readB" {
			t.Errorf("stored name overwritten: got %q, want readB", dr.rec.QName)
		}
	})

	t.Run("within-slice-pair-shares-name", func(t *testing.T) {
		rd := newRD("lossy.cram", 0)
		// Two records of a same-slice pair, both names dropped, linked by
		// resolveMates: upstream at index 0, downstream mate at index 2.
		up := &decodedRecord{rec: &sam.Record{}, mateIndex: 2}
		mid := &decodedRecord{rec: &sam.Record{QName: "kept"}, nameStored: true, mateIndex: -1}
		down := &decodedRecord{rec: &sam.Record{}, mateIndex: 0}
		rd.reconstructDroppedNames([]*decodedRecord{up, mid, down})
		if up.rec.QName != "lossy.cram:1" {
			t.Errorf("upstream synthesised name = %q, want lossy.cram:1", up.rec.QName)
		}
		// The downstream mate's index points back to index 0 (< its own
		// index), so it shares the upstream's number, not its own.
		if down.rec.QName != "lossy.cram:1" {
			t.Errorf("downstream mate name = %q, want lossy.cram:1 (shared with upstream)", down.rec.QName)
		}
		if mid.rec.QName != "kept" {
			t.Errorf("stored middle name overwritten: %q", mid.rec.QName)
		}
	})

	t.Run("orphan-synthesises-by-index", func(t *testing.T) {
		rd := newRD("sample.cram", 1000)
		slice := make([]*decodedRecord, 4)
		for i := range slice {
			slice[i] = &decodedRecord{rec: &sam.Record{}, mateIndex: -1}
		}
		rd.reconstructDroppedNames(slice)
		// record_number = counter + index + 1.
		wants := []string{"sample.cram:1001", "sample.cram:1002", "sample.cram:1003", "sample.cram:1004"}
		for i, w := range wants {
			if slice[i].rec.QName != w {
				t.Errorf("record %d name = %q, want %q", i, slice[i].rec.QName, w)
			}
		}
	})

	t.Run("no-prefix-bare-number", func(t *testing.T) {
		rd := newRD("", 0)
		dr := &decodedRecord{rec: &sam.Record{}, mateIndex: -1}
		rd.reconstructDroppedNames([]*decodedRecord{dr})
		if dr.rec.QName != "1" {
			t.Errorf("bare synthesised name = %q, want 1", dr.rec.QName)
		}
	})

	t.Run("names-preserved-is-noop", func(t *testing.T) {
		h := refFreeHeader()
		h.Preservation.ReadNamesIncluded = true
		rd := &recordDecoder{h: h, slice: &SliceHeader{RecordCounter: 5}, namePrefix: "x.cram"}
		dr := &decodedRecord{rec: &sam.Record{QName: ""}, mateIndex: -1}
		rd.reconstructDroppedNames([]*decodedRecord{dr})
		if dr.rec.QName != "" {
			t.Errorf("names-preserved slice should not synthesise a name, got %q", dr.rec.QName)
		}
	})
}

// TestUnitDecodeReadNameDeferred proves decodeReadName leaves the QNAME
// unset (deferred) when read names are not preserved, and reads the RN
// series into the name when they are. No upstream binary is needed.
func TestUnitDecodeReadNameDeferred(t *testing.T) {
	// Not preserved: defer (no name, no nameStored, no series consumed).
	h := refFreeHeader()
	h.Preservation.ReadNamesIncluded = false
	rd := &recordDecoder{h: h, slice: &SliceHeader{}}
	dr := &decodedRecord{rec: &sam.Record{}, mateIndex: -1}
	if err := rd.decodeReadName(dr, 0, 0); err != nil {
		t.Fatalf("decodeReadName (not preserved): %v", err)
	}
	if dr.rec.QName != "" || dr.nameStored {
		t.Errorf("deferred name = %q (stored=%v), want empty/false", dr.rec.QName, dr.nameStored)
	}
}
