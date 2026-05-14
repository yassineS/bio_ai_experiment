package bcftools

import (
	"path/filepath"
	"testing"
)

// TestReadBCFRegions exercises the BCF region-query path: build a BCF with
// three records on two contigs, build its CSI, and query each chrom.
func TestReadBCFRegions(t *testing.T) {
	const vcfText = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##contig=<ID=chr2,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="DP">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	a	A	T	30	PASS	DP=10
chr1	500	b	C	G	30	PASS	DP=20
chr2	200	c	G	A	30	PASS	DP=30
`
	bcfPath := writeBCFForIndex(t, vcfText)
	if _, err := BuildIndex(bcfPath, IndexOptions{Format: IndexCSI, Force: true}); err != nil {
		t.Fatal(err)
	}

	hdr, recs, err := ReadBCFRegions(bcfPath, []region{{chrom: "chr1", beg: 50, end: 200}})
	if err != nil {
		t.Fatal(err)
	}
	if hdr == nil {
		t.Fatal("nil header")
	}
	if len(recs) != 1 {
		t.Errorf("got %d records, want 1", len(recs))
	}
	if recs[0].ID != "a" {
		t.Errorf("got id %q want a", recs[0].ID)
	}

	hdr2, recs2, err := ReadBCFRegions(bcfPath, []region{{chrom: "chr2", beg: 1, end: 1000}})
	if err != nil {
		t.Fatal(err)
	}
	_ = hdr2
	if len(recs2) != 1 || recs2[0].ID != "c" {
		t.Errorf("chr2: got %d recs (first id %q)", len(recs2), func() string {
			if len(recs2) > 0 {
				return recs2[0].ID
			}
			return ""
		}())
	}
}

// TestReadBCFRegionsUnknownChrom returns no records for an unindexed chrom.
func TestReadBCFRegionsUnknownChrom(t *testing.T) {
	const vcfText = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	a	A	T	30	PASS	.
`
	bcfPath := writeBCFForIndex(t, vcfText)
	if _, err := BuildIndex(bcfPath, IndexOptions{Format: IndexCSI, Force: true}); err != nil {
		t.Fatal(err)
	}
	_, recs, err := ReadBCFRegions(bcfPath, []region{{chrom: "chrUNKNOWN", beg: 1, end: 100}})
	if err != nil {
		t.Fatal(err)
	}
	if len(recs) != 0 {
		t.Errorf("expected no records, got %d", len(recs))
	}
}

// TestHasCSI checks the sibling-file probe.
func TestHasCSI(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "x.bcf")
	if HasCSI(path) {
		t.Error("HasCSI returned true for missing file")
	}
}

// TestReadBCFRegionsMissingIndex surfaces the index-load error path.
func TestReadBCFRegionsMissingIndex(t *testing.T) {
	const vcfText = `##fileformat=VCFv4.2
##contig=<ID=chr1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	100	a	A	T	30	PASS	.
`
	bcfPath := writeBCFForIndex(t, vcfText)
	if _, _, err := ReadBCFRegions(bcfPath, []region{{chrom: "chr1", beg: 1, end: 200}}); err == nil {
		t.Fatal("expected error for missing .csi")
	}
}
