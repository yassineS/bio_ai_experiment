package samtools

import (
	"bytes"
	"io"
	"os"
	"path/filepath"
	"testing"
)

// writeIndexedBAM writes the given SAM text to a coordinate-sorted BAM on disk
// at dir/name and builds its sibling .bai, returning the BAM path. It mirrors
// the on-disk fixture pattern used by the BAI region-query parity tests.
func writeIndexedBAM(t *testing.T, dir, name, samText string) string {
	t.Helper()
	bamPath := filepath.Join(dir, name)
	bf, err := os.Create(bamPath)
	if err != nil {
		t.Fatal(err)
	}
	// Sort consumes the SAM text directly and emits a coordinate-sorted BAM.
	if err := Sort(bytes.NewReader([]byte(samText)), bf, SortOptions{Order: SortCoordinate, OutputBAM: true}); err != nil {
		bf.Close()
		t.Fatalf("Sort: %v", err)
	}
	if err := bf.Close(); err != nil {
		t.Fatal(err)
	}
	if err := IndexFile(bamPath, "", IndexOptions{}); err != nil {
		t.Fatalf("IndexFile(%s): %v", name, err)
	}
	return bamPath
}

// TestMergeRegionRestrictsAndRenames exercises `samtools merge -R STR`: it
// merges two indexed inputs that share an @RG ID, restricting to a region, and
// asserts that (1) only reads overlapping the region survive, (2) reads outside
// the region are dropped, and (3) the colliding @RG from the second input is
// still seeded-renamed and its records retagged exactly as in a whole-file
// merge. This pins the indexed-query path added for -R against the existing
// merge semantics.
func TestMergeRegionRestrictsAndRenames(t *testing.T) {
	dir := t.TempDir()
	hdr := "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000\n@SQ\tSN:chr2\tLN:100000\n@RG\tID:grp\tSM:s1\n"
	// Input A: one read inside the target region (chr1:1000-2000), one well
	// outside it (chr1:50000), one on a different reference (chr2).
	aSAM := hdr +
		"a_in\t0\tchr1\t1500\t60\t5M\t*\t0\t0\tACGTA\tIIIII\tRG:Z:grp\n" +
		"a_out\t0\tchr1\t50000\t60\t5M\t*\t0\t0\tACGTA\tIIIII\tRG:Z:grp\n" +
		"a_chr2\t0\tchr2\t1500\t60\t5M\t*\t0\t0\tACGTA\tIIIII\tRG:Z:grp\n"
	// Input B: one read inside the region, one outside.
	bSAM := hdr +
		"b_in\t0\tchr1\t1600\t60\t5M\t*\t0\t0\tACGTA\tIIIII\tRG:Z:grp\n" +
		"b_out\t0\tchr1\t60000\t60\t5M\t*\t0\t0\tACGTA\tIIIII\tRG:Z:grp\n"

	aPath := writeIndexedBAM(t, dir, "a.bam", aSAM)
	bPath := writeIndexedBAM(t, dir, "b.bam", bSAM)

	var out bytes.Buffer
	if err := MergeFiles([]string{aPath, bPath}, &out, MergeOptions{
		Region:     "chr1:1000-2000",
		RandomSeed: 1,
		SeedSet:    true,
	}); err != nil {
		t.Fatalf("MergeFiles(-R): %v", err)
	}

	rd, err := newBAMReader(out.Bytes())
	if err != nil {
		t.Fatalf("re-read merged BAM: %v", err)
	}
	// The @RG collision must still be renamed under the seed.
	ids := map[string]bool{}
	for _, rg := range rd.Header().ReadGroups {
		ids[rg.ID] = true
	}
	if !ids["grp"] || !ids["grp-055424A4"] {
		t.Errorf("merged @RG IDs = %v, want grp and grp-055424A4", ids)
	}

	// Only the in-region reads survive; the second input's read is retagged.
	wantRG := map[string]string{"a_in": "grp", "b_in": "grp-055424A4"}
	got := map[string]string{}
	for {
		rec, err := rd.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("Read: %v", err)
		}
		got[rec.QName] = recordRGID(rec)
	}
	if len(got) != len(wantRG) {
		t.Fatalf("merged reads = %v, want exactly %v", got, wantRG)
	}
	for q, rg := range wantRG {
		if got[q] != rg {
			t.Errorf("read %s RG = %q, want %q (present=%v)", q, got[q], rg, got)
		}
	}
	// Explicitly confirm the out-of-region and other-reference reads are gone.
	for _, q := range []string{"a_out", "a_chr2", "b_out"} {
		if _, ok := got[q]; ok {
			t.Errorf("read %s should have been excluded by -R chr1:1000-2000", q)
		}
	}
}

// TestMergeRegionUnknownChrom checks that a region naming a reference present in
// the header but with no overlapping reads yields a header-only, record-free
// merge rather than an error — matching `samtools merge -R` on an empty region.
func TestMergeRegionUnknownChrom(t *testing.T) {
	dir := t.TempDir()
	hdr := "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000\n@RG\tID:grp\tSM:s1\n"
	aSAM := hdr + "a1\t0\tchr1\t1500\t60\t5M\t*\t0\t0\tACGTA\tIIIII\tRG:Z:grp\n"
	aPath := writeIndexedBAM(t, dir, "a.bam", aSAM)

	var out bytes.Buffer
	if err := MergeFiles([]string{aPath}, &out, MergeOptions{
		Region: "chr1:90000-95000",
	}); err != nil {
		t.Fatalf("MergeFiles(-R empty region): %v", err)
	}
	rd, err := newBAMReader(out.Bytes())
	if err != nil {
		t.Fatalf("re-read: %v", err)
	}
	if _, err := rd.Read(); err != io.EOF {
		t.Errorf("expected no records for an empty region, got err=%v", err)
	}
}
