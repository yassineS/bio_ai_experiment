package bcftools

import (
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// TestEstSeqQ_OpenExt covers the q = openQ + extQ*(|l|-1) branch from
// bam2bcf_indel.c:85. Picking lRun<3 forces qh = 1000 so q wins.
func TestEstSeqQ_OpenExt(t *testing.T) {
	bca := &bcfCallauxIndel{OpenQ: 40, ExtQ: 20, TandemQ: 500}
	// |l|=1, lRun=1 => qh=1000, q = 40 + 20*0 = 40.
	if got := estSeqQ(bca, 1, 1); got != 40 {
		t.Errorf("estSeqQ(l=1,lRun=1) = %d, want 40", got)
	}
	// |l|=3, lRun=2 => qh=1000, q = 40 + 20*2 = 80.
	if got := estSeqQ(bca, 3, 2); got != 80 {
		t.Errorf("estSeqQ(l=3,lRun=2) = %d, want 80", got)
	}
	// Deletion (negative l) takes |l|: |-2|=2, lRun=1 => q = 60.
	if got := estSeqQ(bca, -2, 1); got != 60 {
		t.Errorf("estSeqQ(l=-2,lRun=1) = %d, want 60", got)
	}
}

// TestEstSeqQ_TandemCap covers the qh = tandemQ * |l|/lRun branch
// (bam2bcf_indel.c:86) — i.e. lRun>=3 with a low |l|. Picking
// tandemQ=500, lRun=10, |l|=1 gives qh = round(500*1/10+0.499) = 50,
// which beats q = 40+0 = 40, so q wins. Bumping lRun to 50 gives
// qh = round(500/50+0.499) = 10, so qh wins.
func TestEstSeqQ_TandemCap(t *testing.T) {
	bca := &bcfCallauxIndel{OpenQ: 40, ExtQ: 20, TandemQ: 500}
	if got := estSeqQ(bca, 1, 10); got != 40 {
		t.Errorf("estSeqQ(l=1,lRun=10) = %d, want 40", got)
	}
	if got := estSeqQ(bca, 1, 50); got != 10 {
		t.Errorf("estSeqQ(l=1,lRun=50) = %d, want 10", got)
	}
}

// TestBcfCgpLRun_Homopolymer: a clean homopolymer "AAAAAA" centred on
// pos=2 yields lRun = length of the run. ref[pos+1] is the homopolymer
// nucleotide (bam2bcf_indel.c:417).
func TestBcfCgpLRun_Homopolymer(t *testing.T) {
	ref := []byte("AAAAAA")
	if got := bcfCgpLRun(ref, 2); got != 6 {
		t.Errorf("bcfCgpLRun(AAAAAA, pos=2) = %d, want 6", got)
	}
}

// TestBcfCgpLRun_AdjacentMismatch: when pos+1 sits at the start of a
// homopolymer the run extends only to the right.
func TestBcfCgpLRun_AdjacentMismatch(t *testing.T) {
	// ref[pos+1] = 'C' at index 1; run "CCC" extends to index 3.
	ref := []byte("ACCCG")
	if got := bcfCgpLRun(ref, 0); got != 3 {
		t.Errorf("bcfCgpLRun(ACCCG, pos=0) = %d, want 3", got)
	}
}

// TestBcfCgpLRun_NBase: when ref[pos+1] is N upstream returns 1
// (bam2bcf_indel.c:419-420).
func TestBcfCgpLRun_NBase(t *testing.T) {
	ref := []byte("ANAA")
	if got := bcfCgpLRun(ref, 0); got != 1 {
		t.Errorf("bcfCgpLRun(ANAA, pos=0) = %d, want 1", got)
	}
}

// TestEstIndelreg_DeletionMatchesRef: when ins4 is nil, est_indelreg
// (bam2bcf_indel.c:90) compares ref[i] to ref[pos+1+j%l]. For a clean
// homopolymer this matches throughout — every score increment is +1
// so max_i tracks i to the end and the returned distance is the
// remaining-length minus 1 from pos.
func TestEstIndelreg_DeletionMatchesRef(t *testing.T) {
	ref := []byte("AAAAAA")
	// pos=0, l=-1 (single-base deletion). ref[1..] is all A so score
	// increments each step. Return value should equal len(ref)-1.
	got := estIndelreg(0, ref, -1, nil)
	if got != len(ref)-1 {
		t.Errorf("estIndelreg(AAAAAA, pos=0, l=-1) = %d, want %d", got, len(ref)-1)
	}
}

// TestEstIndelreg_StopsOnMismatch verifies the score -= 10 / break
// branch (bam2bcf_indel.c:97). With a mismatch right after the
// homopolymer ends, the score goes negative and the loop exits; the
// returned max_i is the last position with score > 0.
func TestEstIndelreg_StopsOnMismatch(t *testing.T) {
	ref := []byte("AAAACG")
	// pos=0, l=-1: walk i=1..5 expecting A. i=1..3 score=1..3; i=4 'C'
	// score=3-10=-7 -> break. max_i = 3 -> return 3-0 = 3.
	if got := estIndelreg(0, ref, -1, nil); got != 3 {
		t.Errorf("estIndelreg(AAAACG, pos=0, l=-1) = %d, want 3", got)
	}
}

// TestEstIndelreg_InsertionConsensus: when ins4 is provided the
// helper compares against the encoded insertion consensus rather than
// the reference itself. Encode an "A" repeat as a single byte 0
// (= 'A' in upstream's "ACGTN" table) and verify the result.
func TestEstIndelreg_InsertionConsensus(t *testing.T) {
	ref := []byte("XAAAAA") // first base does not matter (pos=0 starts scan at 1)
	ins4 := []byte{0}       // 'A'
	if got := estIndelreg(0, ref, 1, ins4); got != 5 {
		t.Errorf("estIndelreg(insertion 'A', AAAAA) = %d, want 5", got)
	}
}

// TestTpos2qpos_SimpleMatch: a 10M read at recPos=100 should map
// tpos=105 to qpos=5 (bam2bcf_indel.c:62).
func TestTpos2qpos_SimpleMatch(t *testing.T) {
	cig := sam.Cigar{sam.CigarOp(uint32(10<<4 | sam.CigarMatch))}
	var tposOut int
	q := tpos2qpos(cig, 100, 105, true, &tposOut)
	if q != 5 || tposOut != 105 {
		t.Errorf("tpos2qpos simple match q=%d tposOut=%d, want q=5 tposOut=105", q, tposOut)
	}
}

// TestTpos2qpos_DeletionIsLeft: 5M2D5M at recPos=100. tpos=106 falls
// in the deletion. isLeft=true clamps to the deletion's left edge.
func TestTpos2qpos_DeletionIsLeft(t *testing.T) {
	cig := sam.Cigar{
		sam.CigarOp(uint32(5<<4 | sam.CigarMatch)),
		sam.CigarOp(uint32(2<<4 | sam.CigarDeletion)),
		sam.CigarOp(uint32(5<<4 | sam.CigarMatch)),
	}
	var tposOut int
	q := tpos2qpos(cig, 100, 106, true, &tposOut)
	// left edge of D op is x=105 (after the 5M).
	if tposOut != 105 || q != 5 {
		t.Errorf("tpos2qpos D isLeft q=%d tposOut=%d, want q=5 tposOut=105", q, tposOut)
	}
}

// TestTpos2qpos_DeletionIsRight: same setup, isLeft=false clamps to
// the deletion's right edge.
func TestTpos2qpos_DeletionIsRight(t *testing.T) {
	cig := sam.Cigar{
		sam.CigarOp(uint32(5<<4 | sam.CigarMatch)),
		sam.CigarOp(uint32(2<<4 | sam.CigarDeletion)),
		sam.CigarOp(uint32(5<<4 | sam.CigarMatch)),
	}
	var tposOut int
	q := tpos2qpos(cig, 100, 106, false, &tposOut)
	if tposOut != 107 || q != 5 {
		t.Errorf("tpos2qpos D !isLeft q=%d tposOut=%d, want q=5 tposOut=107", q, tposOut)
	}
}

// TestTpos2qpos_AfterInsertion: 5M3I5M at recPos=100. tpos=108 falls
// in the second match run; qpos accounts for the +3 insertion offset.
func TestTpos2qpos_AfterInsertion(t *testing.T) {
	cig := sam.Cigar{
		sam.CigarOp(uint32(5<<4 | sam.CigarMatch)),
		sam.CigarOp(uint32(3<<4 | sam.CigarInsertion)),
		sam.CigarOp(uint32(5<<4 | sam.CigarMatch)),
	}
	var tposOut int
	q := tpos2qpos(cig, 100, 108, true, &tposOut)
	// x after 5M = 105 = tpos. The check `if x+l > tpos` with x=105,
	// l=5 gives 110>108 true, so qpos = y(8) + (108-105) = 11.
	if q != 11 || tposOut != 108 {
		t.Errorf("tpos2qpos after I q=%d tposOut=%d, want q=11 tposOut=108", q, tposOut)
	}
}

// TestGetPos_NoSoftClip: a 10M read at qpos=5, qlen=10, no indel: the
// epos bin should be round(5/(10+1) * 100) = 45; scLen=0; end=-1.
func TestGetPos_NoSoftClip(t *testing.T) {
	cig := sam.Cigar{sam.CigarOp(uint32(10<<4 | sam.CigarMatch))}
	r := getPos(cig, 0, 5, 10, 0)
	if r.ScLen != 0 || r.End != -1 || r.SLen != 10 {
		t.Errorf("getPos no-SC: %+v, want ScLen=0 End=-1 SLen=10", r)
	}
	f := float64(5) / float64(11) * 100
	want := int(f)
	if r.EPos != want {
		t.Errorf("getPos EPos=%d, want %d", r.EPos, want)
	}
}

// TestGetPos_LeftSoftClip: 3S7M, qpos=5 (i.e. the 3rd matched base).
// slen drops to 7; epos drops by sc_len(3) to 2. scLen bin computed via
// 15*sc_len/(sc_dist+1).
func TestGetPos_LeftSoftClip(t *testing.T) {
	cig := sam.Cigar{
		sam.CigarOp(uint32(3<<4 | sam.CigarSoftClip)),
		sam.CigarOp(uint32(7<<4 | sam.CigarMatch)),
	}
	r := getPos(cig, 0, 5, 10, 0)
	if r.SLen != 7 {
		t.Errorf("getPos SLen=%d, want 7", r.SLen)
	}
	if r.End != 0 {
		t.Errorf("getPos End=%d, want 0 (left clip)", r.End)
	}
	// sc_len=3, sc_dist=epos(5-3)=2; scLen bin = 15*3/(2+1) = 15.
	if r.ScLen != 15 {
		t.Errorf("getPos ScLen=%d, want 15", r.ScLen)
	}
}

// TestBcfCgpFindTypes_NoIndels: a pile of pure-match reads has no
// indel observations, so bcf_cgp_find_types returns nil (n_types==1).
func TestBcfCgpFindTypes_NoIndels(t *testing.T) {
	bca := &bcfCallauxIndel{
		MinSupport:   2,
		MinFrac:      0.05,
		IndelWinSize: 110,
	}
	ref := []byte("ACGTACGTACGTACGTACGTACGTACGT")
	pile := []pileupBase{
		{qlen: 10},
		{qlen: 10},
		{qlen: 10},
	}
	got := bcfCgpFindTypes([][]pileupBase{pile}, 0, bca, ref)
	if got != nil {
		t.Errorf("bcfCgpFindTypes (no indels) = %+v, want nil", got)
	}
}

// TestBcfCgpFindTypes_WithIndels: a pile with two reads carrying +1
// insertions (out of three reads total) passes the per-site filter
// (na/nt = 2/3 > min_frac=0.05; na=2 >= min_support=2). Returns a
// types slice [0, +1] with REF at index 0.
func TestBcfCgpFindTypes_WithIndels(t *testing.T) {
	bca := &bcfCallauxIndel{
		MinSupport:   2,
		MinFrac:      0.05,
		IndelWinSize: 110,
	}
	ref := []byte("ACGTACGTACGTACGTACGTACGTACGT")
	pile := []pileupBase{
		{qlen: 10, indel: 1},
		{qlen: 10, indel: 1},
		{qlen: 10},
	}
	got := bcfCgpFindTypes([][]pileupBase{pile}, 0, bca, ref)
	if got == nil {
		t.Fatal("bcfCgpFindTypes with 2/3 indels = nil, want a result")
	}
	if len(got.Types) != 2 {
		t.Fatalf("bcfCgpFindTypes types = %+v, want 2 entries", got.Types)
	}
	if got.Types[0] != 0 || got.Types[1] != 1 {
		t.Errorf("bcfCgpFindTypes types = %v, want [0, 1]", got.Types)
	}
	if got.RefType != 0 {
		t.Errorf("bcfCgpFindTypes RefType = %d, want 0", got.RefType)
	}
	if got.N != 3 {
		t.Errorf("bcfCgpFindTypes N = %d, want 3", got.N)
	}
	// MaxSupport must equal 2 and MaxFrac == 2/3.
	if bca.MaxSupport != 2 {
		t.Errorf("MaxSupport = %d, want 2", bca.MaxSupport)
	}
}

// TestBcfCgpFindTypes_BelowSupport: na=1 < min_support=2, so the
// site-level filter rejects.
func TestBcfCgpFindTypes_BelowSupport(t *testing.T) {
	bca := &bcfCallauxIndel{
		MinSupport:   2,
		MinFrac:      0.05,
		IndelWinSize: 110,
	}
	ref := []byte("ACGTACGTACGTACGTACGTACGTACGT")
	pile := []pileupBase{
		{qlen: 10, indel: 1},
		{qlen: 10},
		{qlen: 10},
	}
	if got := bcfCgpFindTypes([][]pileupBase{pile}, 0, bca, ref); got != nil {
		t.Errorf("bcfCgpFindTypes (below support) = %+v, want nil", got)
	}
}

// TestBcfCgpFindTypes_NRichRefSkipped: when the reference window has
// majority N's the position is dropped (bam2bcf_indel.c:241-247).
func TestBcfCgpFindTypes_NRichRefSkipped(t *testing.T) {
	bca := &bcfCallauxIndel{
		MinSupport:   1,
		MinFrac:      0.05,
		IndelWinSize: 5, // small window so N's dominate
	}
	ref := []byte("ANNNNNNNNNNNNNN")
	pile := []pileupBase{
		{qlen: 10, indel: 1},
		{qlen: 10, indel: 1},
		{qlen: 10},
	}
	if got := bcfCgpFindTypes([][]pileupBase{pile}, 0, bca, ref); got != nil {
		t.Errorf("bcfCgpFindTypes (N-rich) = %+v, want nil", got)
	}
}

// TestNewBcfCallauxIndel_Defaults asserts that an empty MpileupOptions
// gives the upstream defaults (mpileup.c:1381-1383 and the indel
// getopt long table) when fed through newBcfCallauxIndel.
func TestNewBcfCallauxIndel_Defaults(t *testing.T) {
	bca := newBcfCallauxIndel(MpileupOptions{})
	if bca.OpenQ != DefaultMpileupOpenProb {
		t.Errorf("OpenQ = %d, want %d", bca.OpenQ, DefaultMpileupOpenProb)
	}
	if bca.ExtQ != DefaultMpileupExtProb {
		t.Errorf("ExtQ = %d, want %d", bca.ExtQ, DefaultMpileupExtProb)
	}
	if bca.TandemQ != DefaultMpileupTandemQual {
		t.Errorf("TandemQ = %d, want %d", bca.TandemQ, DefaultMpileupTandemQual)
	}
	if bca.MinSupport != DefaultMpileupMinIReads {
		t.Errorf("MinSupport = %d, want %d", bca.MinSupport, DefaultMpileupMinIReads)
	}
	if bca.MinFrac != DefaultMpileupGapFrac {
		t.Errorf("MinFrac = %g, want %g", bca.MinFrac, DefaultMpileupGapFrac)
	}
	if bca.IndelBias != DefaultMpileupIndelBias {
		t.Errorf("IndelBias = %g, want %g", bca.IndelBias, DefaultMpileupIndelBias)
	}
	if bca.IndelWinSize != DefaultMpileupIndelSize {
		t.Errorf("IndelWinSize = %d, want %d", bca.IndelWinSize, DefaultMpileupIndelSize)
	}
}

// TestNewBcfCallauxIndel_Overrides: CLI knobs override the defaults
// (sub-slice 4e will then feed bca into the indel-aware emitter).
func TestNewBcfCallauxIndel_Overrides(t *testing.T) {
	bca := newBcfCallauxIndel(MpileupOptions{
		OpenProb:   45,
		ExtProb:    25,
		TandemQual: 250,
		MinIReads:  3,
		GapFrac:    0.1,
		IndelBias:  1.5,
		IndelSize:  90,
	})
	if bca.OpenQ != 45 || bca.ExtQ != 25 || bca.TandemQ != 250 ||
		bca.MinSupport != 3 || bca.MinFrac != 0.1 ||
		bca.IndelBias != 1.5 || bca.IndelWinSize != 90 {
		t.Errorf("override-knob plumb failed: %+v", bca)
	}
}

// TestAccumulateMpileupBases_IndelSet covers the wiring added in
// mpileup.go: an alignment with CIGAR 3M2D3M sets the indel field to
// -2 on the LAST base of the first match run, and 0 elsewhere. The
// rec back-pointer is set on every column (4e.2: every read in the
// pile needs to be available to bcfCallGapPrep for the per-read
// iref/ialt bias accumulation, not just the indel-bearing ones).
func TestAccumulateMpileupBases_IndelSet(t *testing.T) {
	rec := &sam.Record{
		QName: "r",
		RName: "chr1",
		Pos:   1,
		MapQ:  60,
		Cigar: sam.Cigar{
			sam.CigarOp(uint32(3<<4 | sam.CigarMatch)),
			sam.CigarOp(uint32(2<<4 | sam.CigarDeletion)),
			sam.CigarOp(uint32(3<<4 | sam.CigarMatch)),
		},
		Seq:  "AAATTT",
		Qual: []byte{30, 30, 30, 30, 30, 30},
	}
	events := make([][]pileupBase, 10)
	accumulateMpileupBases(rec, events, 0, nil, nil, nil, nil, nil)
	// Ref positions 0..2 are M, 3..4 are D (no bases emitted), 5..7 are M.
	if len(events[2]) != 1 || events[2][0].indel != -2 {
		t.Errorf("events[2].indel = %v, want a single base with indel=-2", events[2])
	}
	if events[2][0].rec != rec {
		t.Errorf("events[2].rec = %v, want rec pointer set", events[2][0].rec)
	}
	if len(events[0]) != 1 || events[0][0].indel != 0 || events[0][0].rec != rec {
		t.Errorf("events[0] indel/rec wrong: %+v", events[0])
	}
	if len(events[5]) != 1 || events[5][0].indel != 0 {
		t.Errorf("events[5].indel = %d, want 0", events[5][0].indel)
	}
}

// TestAccumulateMpileupBases_InsertionIndel: CIGAR 4M3I4M sets indel
// to +3 on the 4th base (last base of first match run) and 0 elsewhere.
func TestAccumulateMpileupBases_InsertionIndel(t *testing.T) {
	rec := &sam.Record{
		QName: "r",
		RName: "chr1",
		Pos:   1,
		MapQ:  60,
		Cigar: sam.Cigar{
			sam.CigarOp(uint32(4<<4 | sam.CigarMatch)),
			sam.CigarOp(uint32(3<<4 | sam.CigarInsertion)),
			sam.CigarOp(uint32(4<<4 | sam.CigarMatch)),
		},
		Seq:  "AAAAGGGTTTT",
		Qual: []byte{30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30},
	}
	events := make([][]pileupBase, 10)
	accumulateMpileupBases(rec, events, 0, nil, nil, nil, nil, nil)
	if events[3][0].indel != 3 {
		t.Errorf("events[3].indel = %d, want 3", events[3][0].indel)
	}
	if events[3][0].rec != rec {
		t.Errorf("events[3].rec not set")
	}
	if events[2][0].indel != 0 {
		t.Errorf("events[2].indel = %d, want 0", events[2][0].indel)
	}
}

// TestAccumulateMpileupBases_MergedDeletions mirrors htslib sam.c:5466-5475:
// consecutive same-type D ops are merged on the column preceding them
// (CIGAR 3M 1D 2D 3M -> events[2].indel == -3, not -1).
func TestAccumulateMpileupBases_MergedDeletions(t *testing.T) {
	rec := &sam.Record{
		QName: "r",
		RName: "chr1",
		Pos:   1,
		MapQ:  60,
		Cigar: sam.Cigar{
			sam.CigarOp(uint32(3<<4 | sam.CigarMatch)),
			sam.CigarOp(uint32(1<<4 | sam.CigarDeletion)),
			sam.CigarOp(uint32(2<<4 | sam.CigarDeletion)),
			sam.CigarOp(uint32(3<<4 | sam.CigarMatch)),
		},
		Seq:  "AAATTT",
		Qual: []byte{30, 30, 30, 30, 30, 30},
	}
	events := make([][]pileupBase, 12)
	accumulateMpileupBases(rec, events, 0, nil, nil, nil, nil, nil)
	if len(events[2]) != 1 || events[2][0].indel != -3 {
		t.Errorf("events[2].indel = %v, want -3 (1D+2D merged)", events[2])
	}
}

// TestAccumulateMpileupBases_MergedInsertions mirrors htslib sam.c:5476-5482:
// consecutive same-type I ops are accumulated (CIGAR 4M 2I 1I 4M ->
// events[3].indel == +3).
func TestAccumulateMpileupBases_MergedInsertions(t *testing.T) {
	rec := &sam.Record{
		QName: "r",
		RName: "chr1",
		Pos:   1,
		MapQ:  60,
		Cigar: sam.Cigar{
			sam.CigarOp(uint32(4<<4 | sam.CigarMatch)),
			sam.CigarOp(uint32(2<<4 | sam.CigarInsertion)),
			sam.CigarOp(uint32(1<<4 | sam.CigarInsertion)),
			sam.CigarOp(uint32(4<<4 | sam.CigarMatch)),
		},
		Seq:  "AAAAGGGTTTT",
		Qual: []byte{30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30},
	}
	events := make([][]pileupBase, 12)
	accumulateMpileupBases(rec, events, 0, nil, nil, nil, nil, nil)
	if events[3][0].indel != 3 {
		t.Errorf("events[3].indel = %d, want 3 (2I+1I merged)", events[3][0].indel)
	}
}

// TestAccumulateMpileupBases_InsertionsAcrossPad mirrors htslib
// sam.c:5476-5482: in an I-led run, intervening CPAD ops are skipped
// while same-type I ops continue to accumulate. CIGAR 4M 2I 1P 1I 4M
// -> events[3].indel == +3.
func TestAccumulateMpileupBases_InsertionsAcrossPad(t *testing.T) {
	rec := &sam.Record{
		QName: "r",
		RName: "chr1",
		Pos:   1,
		MapQ:  60,
		Cigar: sam.Cigar{
			sam.CigarOp(uint32(4<<4 | sam.CigarMatch)),
			sam.CigarOp(uint32(2<<4 | sam.CigarInsertion)),
			sam.CigarOp(uint32(1<<4 | sam.CigarPadding)),
			sam.CigarOp(uint32(1<<4 | sam.CigarInsertion)),
			sam.CigarOp(uint32(4<<4 | sam.CigarMatch)),
		},
		Seq:  "AAAAGGGTTTT",
		Qual: []byte{30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30},
	}
	events := make([][]pileupBase, 12)
	accumulateMpileupBases(rec, events, 0, nil, nil, nil, nil, nil)
	if events[3][0].indel != 3 {
		t.Errorf("events[3].indel = %d, want 3 (2I across CPAD + 1I)", events[3][0].indel)
	}
}
