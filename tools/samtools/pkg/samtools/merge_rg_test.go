package samtools

import (
	"bytes"
	"io"
	"strings"
	"testing"
)

// TestRand48Sequence pins the drand48-compatible PRNG to glibc's lrand48
// stream: seed 1's first draw is 0x055424A4, which is exactly the suffix
// upstream samtools merge appends to a colliding @RG ID under `-s 1`.
func TestRand48Sequence(t *testing.T) {
	r := newRand48(1)
	if got := r.lrand48(); got != 0x055424A4 {
		t.Errorf("lrand48(seed=1)#1 = %08X, want 055424A4", got)
	}
	// genUniqueID formats the same draw as the upstream "-%08X" suffix.
	if id := genUniqueID("rg1", map[string]bool{"rg1": true}, newRand48(1)); id != "rg1-055424A4" {
		t.Errorf("genUniqueID = %q, want rg1-055424A4", id)
	}
}

// TestMergeRenamesCollidingRG verifies that merging two inputs that share an
// @RG ID renames the second one (seeded, so reproducible) and retags that
// input's records, mirroring `samtools merge` (no -c).
func TestMergeRenamesCollidingRG(t *testing.T) {
	hdr := "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@RG\tID:grp\tSM:s1\n"
	a := localSAMToBAM(t, hdr+"r1\t0\tchr1\t100\t60\t5M\t*\t0\t0\tACGTA\tIIIII\tRG:Z:grp\n")
	b := localSAMToBAM(t, hdr+"r2\t0\tchr1\t200\t60\t5M\t*\t0\t0\tACGTA\tIIIII\tRG:Z:grp\n")

	var out bytes.Buffer
	if err := Merge([]io.Reader{bytes.NewReader(a), bytes.NewReader(b)}, &out,
		MergeOptions{RandomSeed: 1, SeedSet: true}); err != nil {
		t.Fatalf("Merge: %v", err)
	}
	rd, err := newBAMReader(out.Bytes())
	if err != nil {
		t.Fatalf("re-read: %v", err)
	}
	h := rd.Header()
	ids := map[string]bool{}
	for _, rg := range h.ReadGroups {
		ids[rg.ID] = true
	}
	if !ids["grp"] || !ids["grp-055424A4"] {
		t.Errorf("merged @RG IDs = %v, want grp and grp-055424A4", ids)
	}
	// r1 keeps grp; r2 (from the second input) is retagged to grp-055424A4.
	want := map[string]string{"r1": "grp", "r2": "grp-055424A4"}
	for {
		rec, err := rd.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("Read: %v", err)
		}
		if got := recordRGID(rec); got != want[rec.QName] {
			t.Errorf("%s RG = %q, want %q", rec.QName, got, want[rec.QName])
		}
	}
}

// TestMergeCombineRGKeepsID checks that -c (CombineRG) collapses the colliding
// @RG into one ID and leaves records untouched.
func TestMergeCombineRGKeepsID(t *testing.T) {
	hdr := "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@RG\tID:grp\tSM:s1\n"
	a := localSAMToBAM(t, hdr+"r1\t0\tchr1\t100\t60\t5M\t*\t0\t0\tACGTA\tIIIII\tRG:Z:grp\n")
	b := localSAMToBAM(t, hdr+"r2\t0\tchr1\t200\t60\t5M\t*\t0\t0\tACGTA\tIIIII\tRG:Z:grp\n")

	var out bytes.Buffer
	if err := Merge([]io.Reader{bytes.NewReader(a), bytes.NewReader(b)}, &out,
		MergeOptions{CombineRG: true, RandomSeed: 1, SeedSet: true}); err != nil {
		t.Fatalf("Merge: %v", err)
	}
	rd, err := newBAMReader(out.Bytes())
	if err != nil {
		t.Fatalf("re-read: %v", err)
	}
	var rgLines int
	for _, rg := range rd.Header().ReadGroups {
		if strings.HasPrefix(rg.ID, "grp") {
			rgLines++
		}
		if strings.Contains(rg.ID, "-") {
			t.Errorf("-c should not rename, got %q", rg.ID)
		}
	}
	if rgLines != 1 {
		t.Errorf("-c should collapse to one @RG, got %d", rgLines)
	}
}
