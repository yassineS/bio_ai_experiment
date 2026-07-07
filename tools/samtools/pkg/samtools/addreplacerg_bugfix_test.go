package samtools

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// TestAddReplaceRG_ThreadsByteIdentical verifies -@/--threads only parallelises
// the BGZF I/O: the emitted BAM is byte-for-byte identical for any worker count
// (the RG-tagging pass is single-threaded, as in upstream bam_addrprg.c). It
// covers both the in-memory reader path (AddReplaceRG) and the path-based
// AddReplaceRGFile path, whose -@ >= 2 raw opener engages the block-parallel
// input inflate.
func TestAddReplaceRG_ThreadsByteIdentical(t *testing.T) {
	bamBytes := localSAMToBAM(t,
		"@HD\tVN:1.6\tSO:coordinate\n"+
			"@SQ\tSN:chr1\tLN:1000\n"+
			"@RG\tID:a\tSM:s1\n"+
			"@RG\tID:b\tSM:s2\n"+
			"r1\t0\tchr1\t10\t60\t5M\t*\t0\t0\tACGTA\tIIIII\tRG:Z:a\n"+
			"r2\t0\tchr1\t20\t60\t5M\t*\t0\t0\tACGTA\tIIIII\tRG:Z:b\n"+
			"r3\t0\tchr1\t30\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n"+
			"r4\t0\tchr1\t40\t60\t5M\t*\t0\t0\tACGTA\tIIIII\tRG:Z:a\n")

	opts := func(threads int) AddReplaceRGOptions {
		return AddReplaceRGOptions{RGLine: "ID:rgz\tSM:sz", Mode: AddReplaceRGOverwriteAll, Threads: threads}
	}

	// In-memory reader path.
	run := func(threads int) []byte {
		var out bytes.Buffer
		if err := AddReplaceRG(bytes.NewReader(bamBytes), &out, opts(threads)); err != nil {
			t.Fatalf("AddReplaceRG -@%d: %v", threads, err)
		}
		return out.Bytes()
	}
	if one, many := run(1), run(4); !bytes.Equal(one, many) {
		t.Errorf("AddReplaceRG -@1 (%d bytes) vs -@4 (%d bytes) differ", len(one), len(many))
	}

	// Path-based AddReplaceRGFile path: the -@ >= 2 raw opener drives the
	// block-parallel BGZF input inflate over the on-disk BAM.
	bamPath := filepath.Join(t.TempDir(), "in.bam")
	if err := os.WriteFile(bamPath, bamBytes, 0o644); err != nil {
		t.Fatalf("write bam: %v", err)
	}
	runFile := func(threads int) []byte {
		var out bytes.Buffer
		if err := AddReplaceRGFile(bamPath, &out, opts(threads)); err != nil {
			t.Fatalf("AddReplaceRGFile -@%d: %v", threads, err)
		}
		return out.Bytes()
	}
	if one, many := runFile(1), runFile(4); !bytes.Equal(one, many) {
		t.Errorf("AddReplaceRGFile -@1 (%d bytes) vs -@4 (%d bytes) differ", len(one), len(many))
	}
}

// TestUnitAddReplaceRGHeaderAndRecord pins bug #5's fix with binary-free unit
// tests over the pure helpers: setRecordRG (the per-record RG:Z: setter for
// each mode) and the header add/replace + prune helpers (removeRG /
// keepOnlyRG) that back overwrite_all / -w.
func TestUnitAddReplaceRGHeaderAndRecord(t *testing.T) {
	// --- setRecordRG: overwrite_all replaces any existing RG ---
	rec := &sam.Record{QName: "r", Aux: []sam.Aux{{Tag: "RG", Type: 'Z', Value: "old"}}}
	setRecordRG(rec, "new", AddReplaceRGOverwriteAll)
	if got, _ := rec.GetAux("RG"); got.Value != "new" {
		t.Errorf("overwrite_all: RG = %v, want new", got.Value)
	}

	// --- setRecordRG: orphan_only leaves an existing RG untouched ---
	rec2 := &sam.Record{QName: "r", Aux: []sam.Aux{{Tag: "RG", Type: 'Z', Value: "old"}}}
	setRecordRG(rec2, "new", AddReplaceRGOrphanOnly)
	if got, _ := rec2.GetAux("RG"); got.Value != "old" {
		t.Errorf("orphan_only: RG = %v, want old (unchanged)", got.Value)
	}

	// --- setRecordRG: orphan_only ADDS RG when the record lacks one ---
	rec3 := &sam.Record{QName: "r"}
	setRecordRG(rec3, "new", AddReplaceRGOrphanOnly)
	if got, ok := rec3.GetAux("RG"); !ok || got.Value != "new" {
		t.Errorf("orphan_only on orphan: RG = %v ok=%v, want new", got.Value, ok)
	}

	// --- header prune helpers ---
	makeHdr := func() *sam.Header {
		return &sam.Header{
			Lines: []sam.HeaderLine{
				{Tag: "HD", Fields: []sam.HeaderField{{Tag: "VN", Value: "1.6"}}},
				{Tag: "RG", Fields: []sam.HeaderField{{Tag: "ID", Value: "old"}, {Tag: "SM", Value: "s1"}}},
				{Tag: "RG", Fields: []sam.HeaderField{{Tag: "ID", Value: "new"}, {Tag: "SM", Value: "s2"}}},
			},
			ReadGroups: []sam.ReadGroup{{ID: "old"}, {ID: "new"}},
		}
	}

	// keepOnlyRG drops every @RG except the named ID (overwrite_all on -r).
	h := makeHdr()
	keepOnlyRG(h, "new")
	if findRG(h, "old") != -1 {
		t.Error("keepOnlyRG: @RG old should be removed")
	}
	if findRG(h, "new") < 0 {
		t.Error("keepOnlyRG: @RG new should be retained")
	}
	rgLines := 0
	for _, l := range h.Lines {
		if l.Tag == "RG" {
			rgLines++
		}
	}
	if rgLines != 1 {
		t.Errorf("keepOnlyRG: %d @RG lines remain, want 1", rgLines)
	}

	// removeRG drops only the named ID (the -w replace-existing path).
	h2 := makeHdr()
	removeRG(h2, "old")
	if findRG(h2, "old") != -1 {
		t.Error("removeRG: @RG old should be gone")
	}
	if findRG(h2, "new") < 0 {
		t.Error("removeRG: @RG new should remain")
	}
}

// TestUnitAddReplaceRGDefaultModeOverwriteAll guards the default-mode parity
// detail: a freshly built header with two @RG lines, after AddReplaceRG with
// a brand-new -r line at the default (overwrite_all) mode, must end with only
// the new @RG line and every record tagged RG:Z:<new>.
func TestUnitAddReplaceRGRLinePrunesHeader(t *testing.T) {
	bamBytes := localSAMToBAM(t,
		"@HD\tVN:1.6\n"+
			"@SQ\tSN:chr1\tLN:100\n"+
			"@RG\tID:old\tSM:s1\n"+
			"r1\t0\tchr1\t10\t60\t5M\t*\t0\t0\tACGTA\tIIIII\tRG:Z:old\n"+
			"r2\t0\tchr1\t20\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n")

	var out bytes.Buffer
	if err := AddReplaceRG(bytes.NewReader(bamBytes), &out, AddReplaceRGOptions{
		RGLine: "ID:new\tSM:s2",
		Mode:   AddReplaceRGOverwriteAll,
	}); err != nil {
		t.Fatalf("AddReplaceRG: %v", err)
	}

	rd, err := newBAMReader(out.Bytes())
	if err != nil {
		t.Fatalf("re-read: %v", err)
	}
	hdr := rd.Header()
	if findRG(hdr, "old") != -1 {
		t.Error("default mode -r should prune @RG old from the header")
	}
	if findRG(hdr, "new") < 0 {
		t.Error("default mode -r should keep @RG new in the header")
	}
	for {
		rec, err := rd.Read()
		if err != nil {
			break
		}
		if got, _ := rec.GetAux("RG"); got.Value != "new" {
			t.Errorf("record %s: RG = %v, want new", rec.QName, got.Value)
		}
	}
}
