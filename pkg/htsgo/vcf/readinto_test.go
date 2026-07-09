package vcf

import (
	"io"
	"strings"
	"testing"
)

const readIntoVCF = `##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1	s2
chr1	100	rs1	A	G	30	PASS	DP=60	GT	0/0	0/1
chr1	200	.	C	T,A	25	PASS	DP=45	GT	0/1	1/2
chr2	100	rs3	AT	A	40	q10	DP=90	GT	1/1	0/1
`

// TestReadIntoMatchesRead checks that ReadInto parses each record identically to
// Read for the value-copyable fields, so the reuse path produces the same stats
// input as the allocating path.
func TestReadIntoMatchesRead(t *testing.T) {
	rr := NewReader(strings.NewReader(readIntoVCF))
	if _, err := rr.ReadHeader(); err != nil {
		t.Fatalf("ReadHeader: %v", err)
	}
	ri := NewReader(strings.NewReader(readIntoVCF))
	if _, err := ri.ReadHeader(); err != nil {
		t.Fatalf("ReadHeader: %v", err)
	}
	var reuse Variant
	for i := 0; ; i++ {
		want, errR := rr.Read()
		errI := ri.ReadInto(&reuse)
		if errR == io.EOF && errI == io.EOF {
			break
		}
		if errR != nil || errI != nil {
			t.Fatalf("record %d: Read err=%v ReadInto err=%v", i, errR, errI)
		}
		if reuse.Chrom != want.Chrom || reuse.Pos != want.Pos || reuse.ID != want.ID ||
			reuse.Ref != want.Ref || reuse.Qual != want.Qual {
			t.Fatalf("record %d scalar mismatch: ReadInto=%+v Read=%+v", i, reuse, want)
		}
		if strings.Join(reuse.Alt, ",") != strings.Join(want.Alt, ",") {
			t.Fatalf("record %d Alt mismatch: ReadInto=%v Read=%v", i, reuse.Alt, want.Alt)
		}
		if strings.Join(reuse.Filter, ";") != strings.Join(want.Filter, ";") {
			t.Fatalf("record %d Filter mismatch: ReadInto=%v Read=%v", i, reuse.Filter, want.Filter)
		}
	}
}

// TestReadIntoInternStable is the aliasing regression guard: it snapshots the
// retained-across-reads string fields (Chrom/Ref and the Alt alleles) from
// every ReadInto and asserts they still hold the correct value after later reads
// overwrote the shared line buffer. Without interning these would alias r.lineBuf
// and read back as a later record's bytes. ID is deliberately excluded: it is not
// interned (no consumer retains it across sites and rsIDs are ~unique), so it must
// not be snapshotted here — a retained v.ID legitimately aliases the buffer.
func TestReadIntoInternStable(t *testing.T) {
	r := NewReader(strings.NewReader(readIntoVCF))
	if _, err := r.ReadHeader(); err != nil {
		t.Fatalf("ReadHeader: %v", err)
	}
	type snap struct {
		chrom, ref string
		alt        []string
	}
	var snaps []snap
	var reuse Variant
	for {
		if err := r.ReadInto(&reuse); err == io.EOF {
			break
		} else if err != nil {
			t.Fatalf("ReadInto: %v", err)
		}
		// Retain the interned strings; the Alt slice backing is reused, so copy it.
		altCopy := make([]string, len(reuse.Alt))
		copy(altCopy, reuse.Alt)
		snaps = append(snaps, snap{reuse.Chrom, reuse.Ref, altCopy})
	}
	want := []snap{
		{"chr1", "A", []string{"G"}},
		{"chr1", "C", []string{"T", "A"}},
		{"chr2", "AT", []string{"A"}},
	}
	if len(snaps) != len(want) {
		t.Fatalf("got %d records, want %d", len(snaps), len(want))
	}
	for i, w := range want {
		g := snaps[i]
		if g.chrom != w.chrom || g.ref != w.ref ||
			strings.Join(g.alt, ",") != strings.Join(w.alt, ",") {
			t.Fatalf("record %d retained value drifted (aliasing): got %+v want %+v", i, g, w)
		}
	}
	// The two chr1 rows must share one interned backing string (same pointer via
	// value equality is all Go exposes; assert equality holds after reuse).
	if snaps[0].chrom != snaps[1].chrom {
		t.Fatalf("interned chr1 strings not equal across reads: %q vs %q", snaps[0].chrom, snaps[1].chrom)
	}
}
