package bcftools

import (
	"bytes"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// TestUnitIUPACAlleleAcrossSamples checks the default iupac_GTs allele folding
// without invoking any binary: alleles present across the given samples are
// OR-ed into one IUPAC ambiguity code (C+A->M, A+G->R), and a site with no set
// allele is skipped.
func TestUnitIUPACAlleleAcrossSamples(t *testing.T) {
	mkVar := func(ref, alt string, gts ...string) *vcf.Variant {
		v := &vcf.Variant{Ref: ref, Alt: []string{alt}}
		for i, g := range gts {
			v.Samples = append(v.Samples, vcf.Sample{
				Name: "S" + string(rune('1'+i)),
				Data: map[string]string{"GT": g},
			})
		}
		return v
	}
	// C ref + A alt, samples 0/0 and 1/1 -> alleles {C,A} -> M.
	got, ok := iupacAlleleAcrossSamples(mkVar("C", "A", "0/0", "1/1"), []int{0, 1}, ConsensusOptions{})
	if !ok || got != "M" {
		t.Errorf("C/A across samples: got %q ok=%v, want M", got, ok)
	}
	// A ref + G alt, one het sample -> {A,G} -> R.
	got, ok = iupacAlleleAcrossSamples(mkVar("A", "G", "0/0", "0/1"), []int{0, 1}, ConsensusOptions{})
	if !ok || got != "R" {
		t.Errorf("A/G across samples: got %q ok=%v, want R", got, ok)
	}
	// All-missing genotypes -> skip (ok=false) without -M.
	if _, ok := iupacAlleleAcrossSamples(mkVar("A", "G", "./.", "./."), []int{0, 1}, ConsensusOptions{}); ok {
		t.Errorf("all-missing should be skipped")
	}
}

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

// TestConsensusMarkDelDoesNotShiftDownstream pins PR #108 review
// finding: with --mark-del padding active, the post-deletion offset
// MUST NOT shift, otherwise the next variant lands at the wrong
// position. The prior code always did `offset += len(alt) - len(ref)`
// even when the padding restored the emitted length to len(ref).
func TestConsensusMarkDelDoesNotShiftDownstream(t *testing.T) {
	ref := []*fasta.Record{
		{ID: "chr1", Sequence: []byte("AAAACCCCGGGG")},
	}
	vcf := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=12>\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chr1\t4\t.\tACC\tA\t.\tPASS\t.\n" +
		"chr1\t10\t.\tG\tT\t.\tPASS\t.\n"
	var out bytes.Buffer
	if _, err := Consensus(strings.NewReader(vcf), &out, ConsensusOptions{
		Reference: ref,
		LineWidth: 80,
		MarkDel:   MarkSpec{Mode: MarkChar, Char: '*'},
	}); err != nil {
		t.Fatalf("Consensus: %v", err)
	}
	// Variant #1 (chr1:4 ACC>A with mark-del) emits 'A**' at positions 4,5,6.
	// Variant #2 (chr1:10 G>T) must land at the ORIGINAL position 10
	// because the padding preserved length. Prior buggy offset shift
	// would land it at position 8 (offset = 1 - 3 = -2).
	want := ">chr1\nAAAA**CCGTGG\n"
	if out.String() != want {
		t.Errorf("downstream variant landed at wrong position:\n got %q\n want %q", out.String(), want)
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

// TestConsensusOverlappingAnchoredDeletions pins the anchor-preservation
// behaviour for two overlapping anchored deletions in a tandem repeat.
// Upstream consensus.c writes a deletion's alt starting at i=trim_beg
// (consensus.c:1014), so an anchored deletion never rewrites its leading
// anchor base — it keeps whatever is already in the buffer. When a second
// anchored deletion overlaps the first exactly on the base the first
// deletion preserved, the first variant's anchor must survive.
//
// Reference:      AAGTCTCTGTGAAAA (1-based positions)
//
//	var1 POS=3 GTC>G : keeps G at pos3, deletes T(4),C(5). frz_pos = pos5.
//	var2 POS=5 CTCTGTG>C : anchored deletion landing exactly on frz_pos
//	  (last base var1 consumed). Its anchor (C) must NOT overwrite the G
//	  that var1 preserved.
//
// Verified byte-for-byte against `bcftools consensus` (upstream): AAGAAAA.
// Before the fix ours emitted AACAAAA (var2's C clobbered var1's G).
func TestConsensusOverlappingAnchoredDeletions(t *testing.T) {
	ref := []*fasta.Record{
		{ID: "chr1", Sequence: []byte("AAGTCTCTGTGAAAA")},
	}
	vcf := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=15>\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		// First anchored deletion: GTC -> G (keeps G at pos 3).
		"chr1\t3\t.\tGTC\tG\t.\tPASS\t.\n" +
		// Second anchored deletion overlaps exactly on the frozen base and
		// is applied (clean anchored indel landing on frz_pos).
		"chr1\t5\t.\tCTCTGTG\tC\t.\tPASS\t.\n"
	var out bytes.Buffer
	n, err := Consensus(strings.NewReader(vcf), &out, ConsensusOptions{
		Reference: ref,
		LineWidth: 80,
	})
	if err != nil {
		t.Fatalf("Consensus: %v", err)
	}
	if n != 2 {
		t.Errorf("applied %d, want 2 (both anchored deletions applied)", n)
	}
	// AA + G (var1 anchor, preserved) + tail AAAA = AAGAAAA. The var2 anchor
	// 'C' must not overwrite the preserved G.
	want := ">chr1\nAAGAAAA\n"
	if out.String() != want {
		t.Errorf("anchor not preserved:\n got %q\n want %q", out.String(), want)
	}
}

// TestParseHaplotypeAliases pins the upstream consensus.c:1312-1313
// shortcuts: "L" === "LR", "S" === "SR".
func TestParseHaplotypeAliases(t *testing.T) {
	cases := []struct {
		in   string
		want HaplotypeSelector
	}{
		{"L", HapLongRef},
		{"LR", HapLongRef},
		{"S", HapShortRef},
		{"SR", HapShortRef},
	}
	for _, c := range cases {
		got, _, err := ParseHaplotypeSelector(c.in)
		if err != nil {
			t.Errorf("ParseHaplotypeSelector(%q): unexpected error %v", c.in, err)
			continue
		}
		if got != c.want {
			t.Errorf("ParseHaplotypeSelector(%q): got %v, want %v", c.in, got, c.want)
		}
	}
}
