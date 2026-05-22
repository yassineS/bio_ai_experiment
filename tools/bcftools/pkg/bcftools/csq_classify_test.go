package bcftools

import (
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/gff"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// buildClassifyIndex builds a forward-strand transcript with an
// explicit exon/CDS/UTR layout so every consequence class can be
// exercised with hand-derived expectations.
//
// Genomic layout (1-based, inclusive):
//
//	transcript span : 1 .. 120
//	exon 1          : 1 .. 30      (5'UTR 1..10, CDS 11..30)
//	intron 1        : 31 .. 50
//	exon 2          : 51 .. 80     (all CDS)
//	intron 2        : 81 .. 100
//	exon 3          : 101 .. 120   (CDS 101..110, 3'UTR 111..120)
//
// CDS coding sequence (forward) starts at genomic 11 with ATG so the
// transcript is a complete annotation (Trim5 == false).
func buildClassifyIndex() *CSQIndex {
	// 120-base reference. Position 11 must begin "ATG".
	// We craft the CDS so codon 1 = ATG (start), and the very last
	// coding codon (genomic 108..110) = TAA (stop).
	ref := make([]byte, 120)
	for i := range ref {
		ref[i] = 'A'
	}
	put := func(pos int, s string) {
		for i := 0; i < len(s); i++ {
			ref[pos-1+i] = s[i]
		}
	}
	// CDS exon1 11..30 (20 bp): ATG GCG TAC AAC TAC GTG C  (last base
	// of codon 7 spills into exon2).
	put(11, "ATGGCGTACAACTACGTGCA")
	// CDS exon2 51..80 (30 bp).
	put(51, "AAGCTGACGACGACGACGACGACGACGACG")
	// CDS exon3 101..110 (10 bp). Make 108..110 = TAA (stop codon).
	put(101, "ACGACGGTAA")

	tx := &CSQTranscript{
		ID:      "txF",
		Gene:    "GENEF",
		Biotype: "protein_coding",
		Chrom:   "chr1",
		Strand:  gff.StrandForward,
		CDSExons: []CSQExon{
			{Start: 11, End: 30},
			{Start: 51, End: 80},
			{Start: 101, End: 110},
		},
		Exons: []CSQExon{
			{Start: 1, End: 30},
			{Start: 51, End: 80},
			{Start: 101, End: 120},
		},
		UTRs: []CSQUTR{
			{Start: 1, End: 10, Prime5: true},
			{Start: 111, End: 120, Prime5: false},
		},
		Beg:    1,
		End:    120,
		Coding: true,
	}
	idx := &CSQIndex{
		Refs:        map[string][]byte{"chr1": ref},
		Transcripts: map[string]*CSQTranscript{"txF": tx},
		ByChrom:     map[string][]*CSQTranscript{"chr1": {tx}},
	}
	return idx
}

// classifyOne is a test helper: classify a single (pos, ref, alt)
// variant against the index and return the joined BCSQ entries.
func classifyOne(idx *CSQIndex, pos int, ref, alt string) []string {
	v := &vcf.Variant{Chrom: "chr1", Pos: pos, Ref: ref, Alt: []string{alt}}
	return classifyCSQRecord(v, idx)
}

// firstConsequence extracts the leading SO term ("missense" from
// "missense|GENE|...").
func firstConsequence(entry string) string {
	if i := strings.IndexByte(entry, '|'); i >= 0 {
		return entry[:i]
	}
	return entry
}

// TestClassifyRegionClasses exercises every region/SO class with a
// SNP whose consequence is hand-derived from csq.c.
func TestClassifyRegionClasses(t *testing.T) {
	idx := buildClassifyIndex()
	cases := []struct {
		name string
		pos  int
		ref  string
		alt  string
		want string // expected leading SO term
	}{
		// 5'UTR: genomic 5 is inside exon1's 5'UTR (1..10).
		{"utr5", 5, "A", "C", "5_prime_utr"},
		// 3'UTR: genomic 115 is inside exon3's 3'UTR (111..120).
		{"utr3", 115, "A", "C", "3_prime_utr"},
		// Intron: genomic 40 sits in intron 1 (31..50), >8bp from
		// either exon edge, so a pure intron consequence.
		{"intron", 40, "A", "C", "intron"},
		// Splice donor: intron 1 base 31 is the first intron base
		// after exon1 (fwd strand => donor side).
		{"splice_donor", 31, "A", "C", "splice_donor"},
		// Splice acceptor: intron 1 base 50 is the last intron base
		// before exon2 (fwd strand => acceptor side).
		{"splice_acceptor", 50, "A", "C", "splice_acceptor"},
		// Splice region (intronic): intron base 35 is within 3..8bp
		// of exon1's end (30) -> splice_region only.
		{"splice_region_intron", 35, "A", "C", "splice_region"},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			got := classifyOne(idx, tc.pos, tc.ref, tc.alt)
			if len(got) != 1 {
				t.Fatalf("got %d entries %v, want 1", len(got), got)
			}
			if c := firstConsequence(got[0]); c != tc.want {
				t.Errorf("pos %d: leading consequence %q, want %q (full %q)", tc.pos, c, tc.want, got[0])
			}
		})
	}
}

// TestClassifySNPCodonClasses pins the CDS codon-level SO terms.
func TestClassifySNPCodonClasses(t *testing.T) {
	idx := buildClassifyIndex()
	// CDS coding: exon1 11..30 holds "ATGGCGTACAACTACGTGCA".
	// codon 1 (11..13) ATG -> M (start)
	// codon 2 (14..16) GCG -> A
	// codon 3 (17..19) TAC -> Y
	// codon 7 last coding codon is genomic 108..110 = TAA -> *.
	cases := []struct {
		name string
		pos  int
		ref  string
		alt  string
		want string
	}{
		// codon 3 TAC: pos 17 T>A => AAC (Asn) -> missense.
		{"missense", 17, "T", "A", "missense"},
		// codon 2 GCG: pos 16 G>A => GCA (Ala) -> synonymous.
		{"synonymous", 16, "G", "A", "synonymous"},
		// codon 3 TAC: pos 19 C>A => TAA (stop) -> stop_gained.
		{"stop_gained", 19, "C", "A", "stop_gained"},
		// start codon ATG: pos 11 A>C => CTG (Leu) -> start_lost.
		{"start_lost", 11, "A", "C", "start_lost"},
		// stop codon TAA at 108..110: pos 110 A>C => TAC (Tyr)
		// -> stop_lost.
		{"stop_lost", 110, "A", "C", "stop_lost"},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			got := classifyOne(idx, tc.pos, tc.ref, tc.alt)
			if len(got) != 1 {
				t.Fatalf("got %d entries %v, want 1", len(got), got)
			}
			if c := firstConsequence(got[0]); c != tc.want {
				t.Errorf("pos %d: leading consequence %q, want %q (full %q)", tc.pos, c, tc.want, got[0])
			}
		})
	}
}

// TestClassifyIndelConsequences exercises the splice_csq indel arms:
// frameshift vs inframe insertion/deletion inside the CDS.
func TestClassifyIndelConsequences(t *testing.T) {
	idx := buildClassifyIndex()
	cases := []struct {
		name string
		pos  int
		ref  string
		alt  string
		want string
	}{
		// 1bp insertion inside CDS exon2 -> frameshift (1 %% 3 != 0).
		{"frameshift_ins", 60, "A", "AT", "frameshift"},
		// 3bp insertion inside CDS exon2 -> inframe_insertion.
		{"inframe_ins", 60, "A", "ATTT", "inframe_insertion"},
		// 1bp deletion inside CDS exon2 -> frameshift.
		{"frameshift_del", 60, "AC", "A", "frameshift"},
		// 3bp deletion inside CDS exon2 -> inframe_deletion.
		{"inframe_del", 60, "ACGA", "A", "inframe_deletion"},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			got := classifyOne(idx, tc.pos, tc.ref, tc.alt)
			if len(got) != 1 {
				t.Fatalf("got %d entries %v, want 1", len(got), got)
			}
			if c := firstConsequence(got[0]); c != tc.want {
				t.Errorf("%s: leading consequence %q, want %q (full %q)", tc.name, c, tc.want, got[0])
			}
		})
	}
}

// TestClassifySpliceSiteIndel pins an indel that lands on a splice
// site: a deletion spanning the last exon base and first intron bases.
// Upstream runs test_cds, test_utr AND test_splice unconditionally
// (csq.c:3733-3736 `hit += ...`), so a CDS-exon-terminal indel stages
// BOTH the CDS-arm entry (which carries the codon-level frame bit) and
// an independent test_splice splice_region entry — two BCSQ items.
func TestClassifySpliceSiteIndel(t *testing.T) {
	idx := buildClassifyIndex()
	// Deletion at genomic 30..32: CDS exon1 ends at 30, intron 1 is
	// 31..50. The variant sits in the terminal 3bp of CDS exon1.
	got := classifyOne(idx, 29, "AGA", "A")
	if len(got) != 2 {
		t.Fatalf("got %d entries %v, want 2 (CDS arm + test_splice arm)", len(got), got)
	}
	var sawSpliceArm bool
	for _, e := range got {
		if !strings.Contains(e, "splice_") {
			t.Errorf("entry %q has no splice consequence", e)
		}
		if firstConsequence(e) == "splice_region" {
			sawSpliceArm = true
		}
	}
	if !sawSpliceArm {
		t.Errorf("expected a standalone splice_region entry from test_splice, got %v", got)
	}
}

// TestClassifyUTRTerminalAdditive pins the additive UTR + splice
// behaviour: a SNP in the terminal 3bp of a UTR exon (at an INTERNAL
// exon edge, i.e. a splice site) stages BOTH a 5'/3'-UTR entry (from
// test_utr) AND a splice_region entry (from test_splice). Upstream
// csq.c:3733-3736 runs both tests unconditionally (`hit += ...`); the
// per-record port must not drop the splice entry.
//
// Fixture: a 3-exon forward transcript whose first exon is entirely
// 5'UTR and whose internal edge (exon1.End) is a real splice site.
//
//	transcript : 1 .. 200
//	exon 1     : 1 .. 40     (all 5'UTR; internal 3' edge at 40)
//	intron 1   : 41 .. 80
//	exon 2     : 81 .. 130   (CDS 81..130, complete: ATG..stop)
//	intron 2   : 131 .. 160
//	exon 3     : 161 .. 200  (CDS 161..170, 3'UTR 171..200)
func TestClassifyUTRTerminalAdditive(t *testing.T) {
	ref := make([]byte, 200)
	for i := range ref {
		ref[i] = 'A'
	}
	put := func(pos int, s string) {
		for i := 0; i < len(s); i++ {
			ref[pos-1+i] = s[i]
		}
	}
	// CDS exon2 81..130 begins with ATG; CDS exon3 161..170 ends TAA.
	put(81, "ATG")
	put(168, "TAA")
	tx := &CSQTranscript{
		ID:      "txU",
		Gene:    "GENEU",
		Biotype: "protein_coding",
		Chrom:   "chr1",
		Strand:  gff.StrandForward,
		CDSExons: []CSQExon{
			{Start: 81, End: 130},
			{Start: 161, End: 170},
		},
		Exons: []CSQExon{
			{Start: 1, End: 40},
			{Start: 81, End: 130},
			{Start: 161, End: 200},
		},
		UTRs: []CSQUTR{
			{Start: 1, End: 40, Prime5: true},
			{Start: 171, End: 200, Prime5: false},
		},
		Beg:    1,
		End:    200,
		Coding: true,
	}
	idx := &CSQIndex{
		Refs:        map[string][]byte{"chr1": ref},
		Transcripts: map[string]*CSQTranscript{"txU": tx},
		ByChrom:     map[string][]*CSQTranscript{"chr1": {tx}},
	}
	// SNP at genomic 39: inside the 5'UTR (1..40) AND within the
	// terminal 3bp (38..40) of exon1, whose 3' edge (40) is an internal
	// splice site. test_utr stages 5_prime_utr; test_splice stages
	// splice_region. Both entries must appear.
	got := classifyOne(idx, 39, "A", "C")
	if len(got) != 2 {
		t.Fatalf("UTR-exon-terminal SNP at pos 39: got %d entries %v, want 2", len(got), got)
	}
	var sawUTR, sawSplice bool
	for _, e := range got {
		switch firstConsequence(e) {
		case "5_prime_utr":
			sawUTR = true
		case "splice_region":
			sawSplice = true
		}
	}
	if !sawUTR || !sawSplice {
		t.Errorf("want both 5_prime_utr and splice_region, got %v", got)
	}
}

// TestClassifyNonCoding verifies a variant in a non-coding transcript
// gets the non_coding SO term, not intron.
func TestClassifyNonCoding(t *testing.T) {
	tx := &CSQTranscript{
		ID:      "txNC",
		Gene:    "LNC1",
		Biotype: "lncRNA",
		Chrom:   "chr1",
		Strand:  gff.StrandForward,
		Exons:   []CSQExon{{Start: 1, End: 50}},
		Beg:     1,
		End:     50,
		Coding:  false,
	}
	idx := &CSQIndex{
		Refs:        map[string][]byte{"chr1": make([]byte, 50)},
		Transcripts: map[string]*CSQTranscript{"txNC": tx},
		ByChrom:     map[string][]*CSQTranscript{"chr1": {tx}},
	}
	idx.Refs["chr1"][24] = 'A'
	got := classifyOne(idx, 25, "A", "C")
	if len(got) != 1 {
		t.Fatalf("got %d entries %v, want 1", len(got), got)
	}
	if c := firstConsequence(got[0]); c != "non_coding" {
		t.Errorf("non-coding transcript: leading consequence %q, want non_coding (full %q)", c, got[0])
	}
}

// TestClassifyReverseStrandSplice verifies that on a reverse-strand
// transcript the donor/acceptor sides flip: the first intron base
// after the genomically-higher exon edge is the acceptor side.
func TestClassifyReverseStrandSplice(t *testing.T) {
	ref := make([]byte, 120)
	for i := range ref {
		ref[i] = 'A'
	}
	// Reverse-strand transcript: CDS exon at the genomically-higher
	// coordinates is the 5' (first) coding exon.
	tx := &CSQTranscript{
		ID:      "txR",
		Gene:    "GENER",
		Biotype: "protein_coding",
		Chrom:   "chr1",
		Strand:  gff.StrandReverse,
		CDSExons: []CSQExon{
			{Start: 11, End: 30},
			{Start: 51, End: 70},
		},
		Exons: []CSQExon{
			{Start: 11, End: 30},
			{Start: 51, End: 70},
		},
		Beg:    11,
		End:    70,
		Coding: true,
		Trim5:  true, // skip start/stop checks for this minimal fixture
		Trim3:  true,
	}
	idx := &CSQIndex{
		Refs:        map[string][]byte{"chr1": ref},
		Transcripts: map[string]*CSQTranscript{"txR": tx},
		ByChrom:     map[string][]*CSQTranscript{"chr1": {tx}},
	}
	// Intron is 31..50. On the reverse strand, intron base 31 (just
	// after exon 11..30) is the acceptor side; base 50 is the donor.
	gotAcc := classifyOne(idx, 31, "A", "C")
	if len(gotAcc) != 1 || !strings.Contains(gotAcc[0], "splice_acceptor") {
		t.Errorf("reverse strand pos 31: want splice_acceptor, got %v", gotAcc)
	}
	gotDon := classifyOne(idx, 50, "A", "C")
	if len(gotDon) != 1 || !strings.Contains(gotDon[0], "splice_donor") {
		t.Errorf("reverse strand pos 50: want splice_donor, got %v", gotDon)
	}
}

// TestConsequencePrecedence pins the kput_vcsq SO-term ordering: when
// several bits are set, the lowest-index term leads and the rest are
// joined with '&'. It also pins kput_vcsq's field gating
// (csq.c:2176-2186): the `|±strand` suffix is printed only when
// CSQ_PRN_STRAND holds (a CSQ_COMPOUND bit is set AND no
// splice/elongation/truncation bit), and the transcript field is empty
// unless CSQ_PRN_TSCRIPT holds (some bit other than intron/non_coding).
func TestConsequencePrecedence(t *testing.T) {
	tx := &CSQTranscript{ID: "t", Gene: "G", Biotype: "protein_coding", Strand: gff.StrandForward}
	cases := []struct {
		csq  uint32
		want string
	}{
		// Pure missense: compound bit, no splice -> strand printed.
		{csqMissense, "missense|G|t|protein_coding|+"},
		// missense&splice_region: a splice bit is set, so CSQ_PRN_STRAND
		// is false -> no strand suffix.
		{csqMissense | csqSpliceRegion, "missense&splice_region|G|t|protein_coding"},
		// splice_donor&splice_region: no compound bit -> no strand.
		{csqSpliceDonor | csqSpliceRegion, "splice_donor&splice_region|G|t|protein_coding"},
		// missense is dropped when a start/stop change co-occurs;
		// stop_gained is compound and has no splice bit -> strand printed.
		{csqMissense | csqStopGained, "stop_gained|G|t|protein_coding|+"},
		// frameshift is compound but co-occurs with a splice bit -> no
		// strand suffix.
		{csqFrameshift | csqSpliceRegion, "frameshift&splice_region|G|t|protein_coding"},
		// Empty mask -> coding_sequence; no compound bit (no strand) and
		// CSQ_PRN_TSCRIPT false (empty transcript field).
		{0, "coding_sequence|G||protein_coding"},
	}
	for _, tc := range cases {
		got := formatConsequence(tx, tc.csq)
		if got != tc.want {
			t.Errorf("formatConsequence(%#x) = %q, want %q", tc.csq, got, tc.want)
		}
	}
}

// TestConsequenceFieldGating pins the per-SO-term field gating of
// formatConsequence against kput_vcsq (csq.c:2176-2186): intron /
// non_coding records print an EMPTY transcript field, and the strand
// suffix is gated on CSQ_PRN_STRAND.
func TestConsequenceFieldGating(t *testing.T) {
	fwd := &CSQTranscript{ID: "txF", Gene: "GENEF", Biotype: "protein_coding", Strand: gff.StrandForward}
	nc := &CSQTranscript{ID: "txNC", Gene: "LNC1", Biotype: "lncRNA", Strand: gff.StrandReverse}
	cases := []struct {
		name string
		tx   *CSQTranscript
		csq  uint32
		want string
	}{
		// intron: CSQ_PRN_TSCRIPT false -> empty transcript; no compound
		// bit -> no strand.
		{"intron", fwd, csqIntron, "intron|GENEF||protein_coding"},
		// non_coding: same gating, empty transcript, no strand.
		{"non_coding", nc, csqNonCoding, "non_coding|LNC1||lncRNA"},
		// splice_region: transcript printed (not intron/non_coding), but
		// no compound bit -> no strand.
		{"splice_region", fwd, csqSpliceRegion, "splice_region|GENEF|txF|protein_coding"},
		// synonymous: compound, no splice -> strand printed.
		{"synonymous", fwd, csqSynonymous, "synonymous|GENEF|txF|protein_coding|+"},
		// 5_prime_utr: transcript printed, no compound bit -> no strand.
		{"utr5", fwd, csqUTR5, "5_prime_utr|GENEF|txF|protein_coding"},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			if got := formatConsequence(tc.tx, tc.csq); got != tc.want {
				t.Errorf("formatConsequence = %q, want %q", got, tc.want)
			}
		})
	}
}

// TestSpliceCSQTrim verifies the allele-trimming step of spliceCSQ
// matches upstream: identical bases are trimmed from the right then
// the left, and a no-op variant (ref==alt) returns spliceVarRef.
func TestSpliceCSQTrim(t *testing.T) {
	tx := &CSQTranscript{Strand: gff.StrandForward}
	// ACGT > ACGT inside an exon: not a real variant.
	sc := newSpliceCtx(tx, 20, "ACGT", "ACGT")
	if ret := spliceCSQ(&sc, 1, 100); ret != spliceVarRef {
		t.Errorf("ACGT>ACGT: ret = %d, want spliceVarRef", ret)
	}
	// CAT > CGT: an MNP; the shared C and T trim, leaving a 1bp
	// substitution at the middle base.
	sc = newSpliceCtx(tx, 20, "CAT", "CGT")
	spliceCSQ(&sc, 1, 100)
	if sc.tbeg != 1 || sc.tend != 1 {
		t.Errorf("CAT>CGT: tbeg=%d tend=%d, want 1,1", sc.tbeg, sc.tend)
	}
}

// TestClassifyMultiContigDispatch verifies that transcripts on other
// contigs are not consulted.
func TestClassifyMultiContigDispatch(t *testing.T) {
	idx := buildClassifyIndex()
	v := &vcf.Variant{Chrom: "chrX", Pos: 17, Ref: "T", Alt: []string{"A"}}
	if got := classifyCSQRecord(v, idx); len(got) != 0 {
		t.Errorf("variant on chrX should not match chr1 transcript, got %v", got)
	}
}

// TestClassifyOutsideTranscript verifies a variant well past the
// transcript span produces no consequence.
func TestClassifyOutsideTranscript(t *testing.T) {
	idx := buildClassifyIndex()
	if got := classifyOne(idx, 200, "A", "C"); len(got) != 0 {
		t.Errorf("variant at pos 200 (transcript ends at 120) should be unannotated, got %v", got)
	}
}
