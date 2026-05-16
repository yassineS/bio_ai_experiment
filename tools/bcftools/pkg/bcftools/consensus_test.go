package bcftools

import (
	"bytes"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fasta"
)

// TestConsensusSNPApplyAllAlts verifies the "no -s" default of applying
// the first ALT at each record's position.
func TestConsensusSNPApplyAllAlts(t *testing.T) {
	ref := []*fasta.Record{
		{ID: "chr1", Sequence: []byte("AAAACCCCGGGGTTTT")},
	}
	vcf := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=16>\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chr1\t1\t.\tA\tG\t.\tPASS\t.\n" +
		"chr1\t9\t.\tG\tT\t.\tPASS\t.\n"
	var out bytes.Buffer
	n, err := Consensus(strings.NewReader(vcf), &out, ConsensusOptions{
		Reference: ref,
		LineWidth: 80,
	})
	if err != nil {
		t.Fatalf("Consensus: %v", err)
	}
	if n != 2 {
		t.Fatalf("applied %d variants, want 2", n)
	}
	want := ">chr1\nGAAACCCCTGGGTTTT\n"
	if out.String() != want {
		t.Errorf("output:\n got %q\n want %q", out.String(), want)
	}
}

// TestConsensusInsertion verifies a simple insertion (REF=A, ALT=AC).
func TestConsensusInsertion(t *testing.T) {
	ref := []*fasta.Record{
		{ID: "chr1", Sequence: []byte("AAAACCCC")},
	}
	vcf := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=8>\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chr1\t4\t.\tA\tAGGG\t.\tPASS\t.\n"
	var out bytes.Buffer
	if _, err := Consensus(strings.NewReader(vcf), &out, ConsensusOptions{
		Reference: ref,
		LineWidth: 80,
	}); err != nil {
		t.Fatalf("Consensus: %v", err)
	}
	want := ">chr1\nAAAAGGGCCCC\n"
	if out.String() != want {
		t.Errorf("output:\n got %q\n want %q", out.String(), want)
	}
}

// TestConsensusDeletion verifies a simple deletion (REF=ACC, ALT=A).
func TestConsensusDeletion(t *testing.T) {
	ref := []*fasta.Record{
		{ID: "chr1", Sequence: []byte("AAAACCCCGGGG")},
	}
	vcf := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=12>\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chr1\t4\t.\tACC\tA\t.\tPASS\t.\n"
	var out bytes.Buffer
	if _, err := Consensus(strings.NewReader(vcf), &out, ConsensusOptions{
		Reference: ref,
		LineWidth: 80,
	}); err != nil {
		t.Fatalf("Consensus: %v", err)
	}
	want := ">chr1\nAAAACCGGGG\n"
	if out.String() != want {
		t.Errorf("output:\n got %q\n want %q", out.String(), want)
	}
}

// TestConsensusMarkDel pads deletions with the --mark-del character so
// downstream coordinates stay aligned.
func TestConsensusMarkDel(t *testing.T) {
	ref := []*fasta.Record{
		{ID: "chr1", Sequence: []byte("AAAACCCCGGGG")},
	}
	vcf := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=12>\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chr1\t4\t.\tACC\tA\t.\tPASS\t.\n"
	var out bytes.Buffer
	if _, err := Consensus(strings.NewReader(vcf), &out, ConsensusOptions{
		Reference: ref,
		LineWidth: 80,
		MarkDel:   MarkSpec{Mode: MarkChar, Char: '*'},
	}); err != nil {
		t.Fatalf("Consensus: %v", err)
	}
	// REF=ACC at POS=4 means positions 4,5,6. ALT=A keeps the 4th base
	// and pads the 5th-6th with '*'.
	want := ">chr1\nAAAA**CCGGGG\n"
	if out.String() != want {
		t.Errorf("output:\n got %q\n want %q", out.String(), want)
	}
}

// TestConsensusSamplePicksGT picks the right ALT based on a sample's GT.
func TestConsensusSamplePicksGT(t *testing.T) {
	ref := []*fasta.Record{
		{ID: "chr1", Sequence: []byte("AAAACCCC")},
	}
	vcf := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=8>\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n" +
		// S1 has 0/0 (ref) at POS 1; S2 has 1/1 (alt). With -s S2 we
		// expect the ALT to be applied.
		"chr1\t1\t.\tA\tG\t.\tPASS\t.\tGT\t0/0\t1/1\n"
	var out bytes.Buffer
	if _, err := Consensus(strings.NewReader(vcf), &out, ConsensusOptions{
		Reference: ref,
		Sample:    "S2",
		LineWidth: 80,
	}); err != nil {
		t.Fatalf("Consensus: %v", err)
	}
	if !strings.HasPrefix(out.String(), ">chr1\nGAAACCCC\n") {
		t.Errorf("S2 hom-alt: expected GAAACCCC, got %q", out.String())
	}

	// Re-run with S1 (hom-ref) — should be unchanged.
	out.Reset()
	if _, err := Consensus(strings.NewReader(vcf), &out, ConsensusOptions{
		Reference: ref,
		Sample:    "S1",
		LineWidth: 80,
	}); err != nil {
		t.Fatalf("Consensus: %v", err)
	}
	if !strings.HasPrefix(out.String(), ">chr1\nAAAACCCC\n") {
		t.Errorf("S1 hom-ref: expected AAAACCCC, got %q", out.String())
	}
}

// TestConsensusHaplotypeR picks REF in het sites.
func TestConsensusHaplotypeR(t *testing.T) {
	ref := []*fasta.Record{
		{ID: "chr1", Sequence: []byte("AAAA")},
	}
	vcf := "##fileformat=VCFv4.2\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n" +
		"chr1\t2\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\n"
	var out bytes.Buffer
	if _, err := Consensus(strings.NewReader(vcf), &out, ConsensusOptions{
		Reference: ref,
		Sample:    "S1",
		Haplotype: HapRef,
		LineWidth: 80,
	}); err != nil {
		t.Fatalf("Consensus: %v", err)
	}
	if !strings.HasPrefix(out.String(), ">chr1\nAAAA\n") {
		t.Errorf("hap R het: expected unchanged AAAA, got %q", out.String())
	}
}

// TestConsensusHaplotypeA picks ALT in het sites.
func TestConsensusHaplotypeA(t *testing.T) {
	ref := []*fasta.Record{
		{ID: "chr1", Sequence: []byte("AAAA")},
	}
	vcf := "##fileformat=VCFv4.2\n" +
		"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n" +
		"chr1\t2\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\n"
	var out bytes.Buffer
	if _, err := Consensus(strings.NewReader(vcf), &out, ConsensusOptions{
		Reference: ref,
		Sample:    "S1",
		Haplotype: HapAlt,
		LineWidth: 80,
	}); err != nil {
		t.Fatalf("Consensus: %v", err)
	}
	if !strings.HasPrefix(out.String(), ">chr1\nATAA\n") {
		t.Errorf("hap A het: expected ATAA, got %q", out.String())
	}
}

// TestConsensusPrefix prepends the prefix to the FASTA header.
func TestConsensusPrefix(t *testing.T) {
	ref := []*fasta.Record{
		{ID: "chr1", Sequence: []byte("AAAA")},
	}
	vcf := "##fileformat=VCFv4.2\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	var out bytes.Buffer
	if _, err := Consensus(strings.NewReader(vcf), &out, ConsensusOptions{
		Reference: ref,
		Prefix:    "sample_",
		LineWidth: 80,
	}); err != nil {
		t.Fatalf("Consensus: %v", err)
	}
	if !strings.HasPrefix(out.String(), ">sample_chr1\n") {
		t.Errorf("prefix: got %q", out.String())
	}
}

// TestConsensusMarkSnvUC highlights substituted bases.
func TestConsensusMarkSnvUC(t *testing.T) {
	ref := []*fasta.Record{
		{ID: "chr1", Sequence: []byte("aaaa")},
	}
	vcf := "##fileformat=VCFv4.2\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chr1\t2\t.\ta\tg\t.\tPASS\t.\n"
	var out bytes.Buffer
	if _, err := Consensus(strings.NewReader(vcf), &out, ConsensusOptions{
		Reference: ref,
		MarkSnv:   MarkSpec{Mode: MarkUpper},
		LineWidth: 80,
	}); err != nil {
		t.Fatalf("Consensus: %v", err)
	}
	if !strings.HasPrefix(out.String(), ">chr1\naGaa\n") {
		t.Errorf("mark-snv uc: got %q", out.String())
	}
}

// TestConsensusMaskN replaces masked positions with N.
func TestConsensusMaskN(t *testing.T) {
	ref := []*fasta.Record{
		{ID: "chr1", Sequence: []byte("AAAACCCC")},
	}
	vcf := "##fileformat=VCFv4.2\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
	var out bytes.Buffer
	if _, err := Consensus(strings.NewReader(vcf), &out, ConsensusOptions{
		Reference: ref,
		LineWidth: 80,
		Mask:      []MaskRegion{{Chrom: "chr1", Beg: 3, End: 5}},
	}); err != nil {
		t.Fatalf("Consensus: %v", err)
	}
	if !strings.HasPrefix(out.String(), ">chr1\nAANNNCCC\n") {
		t.Errorf("mask N: got %q", out.String())
	}
}

func TestParseMarkSpec(t *testing.T) {
	cases := []struct {
		in   string
		want MarkSpec
		err  bool
	}{
		{"", MarkSpec{Mode: MarkNone}, false},
		{"uc", MarkSpec{Mode: MarkUpper}, false},
		{"lc", MarkSpec{Mode: MarkLower}, false},
		{"*", MarkSpec{Mode: MarkChar, Char: '*'}, false},
		{"abc", MarkSpec{}, true},
	}
	for _, tc := range cases {
		got, err := ParseMarkSpec(tc.in)
		if tc.err && err == nil {
			t.Errorf("ParseMarkSpec(%q) expected error", tc.in)
		}
		if !tc.err && (got.Mode != tc.want.Mode || got.Char != tc.want.Char) {
			t.Errorf("ParseMarkSpec(%q) = %+v, want %+v", tc.in, got, tc.want)
		}
	}
}

func TestParseHaplotypeSelector(t *testing.T) {
	cases := []struct {
		in       string
		wantSel  HaplotypeSelector
		wantIdx  int
		wantsErr bool
	}{
		{"", HapAuto, 0, false},
		{"R", HapRef, 0, false},
		{"A", HapAlt, 0, false},
		{"I", HapIUPAC, 0, false},
		{"LR", HapLongRef, 0, false},
		{"LA", HapLongAlt, 0, false},
		{"SR", HapShortRef, 0, false},
		{"SA", HapShortAlt, 0, false},
		{"1", HapIndex, 1, false},
		{"3", HapIndex, 3, false},
		{"0", HapAuto, 0, true},
		{"banana", HapAuto, 0, true},
	}
	for _, tc := range cases {
		gotSel, gotIdx, err := ParseHaplotypeSelector(tc.in)
		if tc.wantsErr && err == nil {
			t.Errorf("ParseHaplotypeSelector(%q) expected error", tc.in)
		}
		if !tc.wantsErr {
			if gotSel != tc.wantSel || gotIdx != tc.wantIdx {
				t.Errorf("ParseHaplotypeSelector(%q) = (%v,%d), want (%v,%d)", tc.in, gotSel, gotIdx, tc.wantSel, tc.wantIdx)
			}
		}
	}
}

// TestConsensusOverlappingFirstWins documents the v1 "first wins" rule
// for overlapping variants. A clean upstream port would emit a warning
// and apply the more left-aligned variant.
func TestConsensusOverlappingFirstWins(t *testing.T) {
	ref := []*fasta.Record{
		{ID: "chr1", Sequence: []byte("AAAACCCCGGGG")},
	}
	vcf := "##fileformat=VCFv4.2\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		// First variant deletes AACC (POS 3, REF=AACC ALT=A).
		"chr1\t3\t.\tAACC\tA\t.\tPASS\t.\n" +
		// Second variant overlaps -> skipped.
		"chr1\t5\t.\tC\tT\t.\tPASS\t.\n"
	var out bytes.Buffer
	n, err := Consensus(strings.NewReader(vcf), &out, ConsensusOptions{
		Reference: ref,
		LineWidth: 80,
	})
	if err != nil {
		t.Fatalf("Consensus: %v", err)
	}
	if n != 1 {
		t.Errorf("applied %d, want exactly 1 (overlap skipped)", n)
	}
	// First variant: AAAACCCCGGGG -> AA + A + CCGGGG = AAACCGGGG (3 bases
	// removed by REF=AACC -> ALT=A). Second variant (POS=5) overlaps and
	// is skipped.
	if !strings.HasPrefix(out.String(), ">chr1\nAAACCGGGG\n") {
		t.Errorf("overlap-first-wins: got %q", out.String())
	}
}
