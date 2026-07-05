package samtools

import (
	"strings"
	"testing"
)

// insPadRunningMinSAM is the synthetic regression fixture for task #58: the
// bayesian insertion-column PAD branch must run a minimum of the carried
// quality against the base AFTER the insertion point, mirroring upstream's
// consensus_pileup.c:188-189
//
//	if (p->seq_offset < b->core.l_qseq)
//	    p->qual = MIN(p->qual, p->b_qual[p->seq_offset+1]);
//
// Three reads insert a 'T' after reference position 2 (CIGAR 2M1I4M); a fourth
// read carries no insertion (CIGAR 6M) and therefore contributes a '*' pad to
// the nth==1 insertion column. The pad read's base AT the reference position
// (seq_offset == 1, the 'C') is high quality ('I' == Q40), but the base
// immediately AFTER it (seq_offset+1 == 2, the 'G') is Q0 ('!'). Upstream's
// running minimum pulls the pad's quality down to Q0; the pre-fix port used the
// raw post-column base quality (e.qual == Q40), leaving the '*' pad far too
// confident.
//
// The pad-read quals are "II!IIII" so the discriminating byte is the Q0 'G'
// that follows the insertion point.
const insPadRunningMinSAM = "@HD\tVN:1.6\tSO:coordinate\n" +
	"@SQ\tSN:chr1\tLN:20\n" +
	"ins1\t0\tchr1\t1\t60\t2M1I4M\t*\t0\t0\tACTGTAC\tIIIIIII\n" +
	"ins2\t0\tchr1\t1\t60\t2M1I4M\t*\t0\t0\tACTGTAC\tIIIIIII\n" +
	"ins3\t0\tchr1\t1\t60\t2M1I4M\t*\t0\t0\tACTGTAC\tIIIIIII\n" +
	"pad1\t0\tchr1\t1\t60\t6M\t*\t0\t0\tACGTAC\tII!IIII\n"

// TestConsensus_Bayesian_InsertionPadQual_RunningMin is the task #58 regression.
// It asserts that the bayesian caller applies the insertion-column PAD
// running-min so the '*' pad quality follows upstream, AND that the inserted
// base call itself is unchanged (the fix is quality-only).
//
// Before the fix the pad read kept its raw post-column quality (Q40), so the
// per-read pileup pad-qual byte was 'I' and the bayesian consensus quality at
// the insertion column was depressed (the over-confident pad fought the
// inserted 'T'). After the fix the pad qual byte is '!' (Q0, the running-min
// against the post-insertion 'G') and the consensus quality rises, exactly as
// upstream emits. The test FAILS before the fix and PASSES after.
func TestConsensus_Bayesian_InsertionPadQual_RunningMin(t *testing.T) {
	// --- Pileup: the per-read pad-qual byte is the direct discriminator. ---
	pile := runConsensusOnSAM(t, insPadRunningMinSAM, ConsensusOptions{
		Format: ConsensusPileup,
		Mode:   ConsensusModeBayesian,
		// Insertions are shown by default (NoShowIns == false), reproducing
		// upstream's show_ins=1.
	})
	var insLine string
	for _, ln := range strings.Split(strings.TrimRight(pile, "\n"), "\n") {
		f := strings.Split(ln, "\t")
		// pos==2, nth==1 is the inserted column following reference position 2.
		if len(f) >= 8 && f[1] == "2" && f[2] == "1" {
			insLine = ln
			break
		}
	}
	if insLine == "" {
		t.Fatalf("no nth==1 insertion column emitted at pos 2:\n%s", pile)
	}
	f := strings.Split(insLine, "\t")
	call, seq, qual := f[4], f[6], f[7]

	// The base call for the insertion column must be 'T' (the inserted base),
	// unchanged by the quality-only fix.
	if call != "T" {
		t.Errorf("insertion-column call = %q, want \"T\" (the fix must not "+
			"change the base call)\nline: %s", call, insLine)
	}
	// The pad read is the 4th read in sorted order; its '*' pad must carry the
	// running-min quality (Q0 == '!'), not the raw post-column quality (Q40 ==
	// 'I'). The seq column is "TTT*" (three inserts then the pad).
	if !strings.HasSuffix(seq, "*") {
		t.Fatalf("insertion-column seq = %q, want a trailing '*' pad\nline: %s", seq, insLine)
	}
	padQual := qual[len(qual)-1]
	if padQual != '!' {
		t.Errorf("insertion-column pad qual = %q, want '!' (Q0, the running-min "+
			"MIN(e.qual, rec.Qual[seqOff+1]) against the post-insertion base). "+
			"Got the raw post-column quality instead — the #58 running-min was "+
			"not applied.\nline: %s", string(padQual), insLine)
	}

	// --- FASTQ: the inserted base's consensus quality must rise once the
	// over-confident pad is correctly down-weighted. ---
	fq := runConsensusOnSAM(t, insPadRunningMinSAM, ConsensusOptions{
		Format: ConsensusFASTQ,
		Mode:   ConsensusModeBayesian,
	})
	lines := strings.Split(strings.TrimRight(fq, "\n"), "\n")
	if len(lines) != 4 || !strings.HasPrefix(lines[0], "@") || lines[2] != "+" {
		t.Fatalf("unexpected FASTQ shape:\n%s", fq)
	}
	seqLine, qualLine := lines[1], lines[3]
	// The emitted SEQ embeds the inserted base: "ACTGTAC". This must be
	// unchanged by the quality-only fix (the gate: 0 base-call diffs).
	const wantSeq = "ACTGTAC"
	if seqLine != wantSeq {
		t.Errorf("FASTQ seq = %q, want %q (the fix is quality-only and must not "+
			"change the consensus sequence)", seqLine, wantSeq)
	}
	if len(qualLine) != len(seqLine) {
		t.Fatalf("FASTQ qual/seq length mismatch: %d vs %d", len(qualLine), len(seqLine))
	}
	// Index 2 is the inserted 'T'. Before the fix the over-confident pad
	// depressed its quality ('.' == Q13); after the fix it rises ('J' == Q41).
	insQ := qualLine[2]
	if insQ <= '.' {
		t.Errorf("inserted-base FASTQ quality = %q (Q%d), want it raised by the "+
			"running-min down-weighting the over-confident pad (pre-fix it sits "+
			"at '.' == Q13).\nseq=%q qual=%q", string(insQ), int(insQ)-33, seqLine, qualLine)
	}
}

// insDelMultiNthSAM exercises the STATEFUL running-min pre-pass across a
// MULTI-nth insertion column with a DELETION-spanning read (task #58, bayesian
// port). Three reads insert two bases ("TT") after reference position 3
// (CIGAR 3M2I3M), producing two insertion columns (nth==1 and nth==2). A
// fourth read deletes reference position 3 (CIGAR 2M1D3M) and therefore
// contributes a '*' pad to BOTH insertion columns.
//
// The deletion read's pad quality must follow upstream's stateful engine
// (consensus_pileup.c:180-228): entering the insertion it is already inside a
// deletion run, so seq_offset is PINNED at the PRE-gap base (query index 1,
// the 'C' at reference position 2) and p->qual already holds the deletion
// running minimum MIN(pre-gap, post-gap). Its per-base quals are "I5III", so
// the pre-gap 'C' is Q20 ('5') while the post-gap 'G' (query index 2) is Q40
// ('I'). The running minimum against b_qual[seq_offset+1] (== the post-gap
// Q40) keeps the pad at Q20 ('5'), CARRIED unchanged across BOTH nth columns.
//
// The pre-fix bayesian builder seeded each pad independently from the raw
// post-gap quality (seq_offset = readBP-1, the POST-gap index) and applied a
// single non-carried MIN, so it would have emitted Q40 ('I') — matching
// neither the simple builder nor upstream. This fixture pins the fix: the
// bayesian pad qual must be Q20 ('5') at both nth columns and equal to the
// simple builder byte-for-byte.
const insDelMultiNthSAM = "@HD\tVN:1.6\tSO:coordinate\n" +
	"@SQ\tSN:chr1\tLN:20\n" +
	"ins1\t0\tchr1\t1\t60\t3M2I3M\t*\t0\t0\tACGTTTAC\tIIIIIIII\n" +
	"ins2\t0\tchr1\t1\t60\t3M2I3M\t*\t0\t0\tACGTTTAC\tIIIIIIII\n" +
	"ins3\t0\tchr1\t1\t60\t3M2I3M\t*\t0\t0\tACGTTTAC\tIIIIIIII\n" +
	"del1\t0\tchr1\t1\t60\t2M1D3M\t*\t0\t0\tACGTA\tI5III\n"

// insertionColumnLine returns the pileup row for reference position pos with
// insertion index nth (both as decimal strings), or "" if absent.
func insertionColumnLine(pile, pos, nth string) string {
	for _, ln := range strings.Split(strings.TrimRight(pile, "\n"), "\n") {
		f := strings.Split(ln, "\t")
		if len(f) >= 8 && f[1] == pos && f[2] == nth {
			return ln
		}
	}
	return ""
}

// TestConsensus_Bayesian_InsertionPadQual_StatefulMultiNth asserts the bayesian
// insertion builder replays the SAME stateful running-min pre-pass as the simple
// builder / upstream across a multi-nth insertion column that a deletion read
// pads: (1) the deletion-spanning pad qual is Q20 ('5', the PRE-gap base pulled
// in via the deletion running-min + pinned pre-gap seq_offset) at BOTH nth
// columns, and (2) the bayesian per-read seq/qual bytes equal the simple
// builder's byte-for-byte, confirming both consume the shared pre-pass.
func TestConsensus_Bayesian_InsertionPadQual_StatefulMultiNth(t *testing.T) {
	simple := runConsensusOnSAM(t, insDelMultiNthSAM, ConsensusOptions{
		Format: ConsensusPileup,
		Mode:   ConsensusModeSimple,
	})
	bayes := runConsensusOnSAM(t, insDelMultiNthSAM, ConsensusOptions{
		Format: ConsensusPileup,
		Mode:   ConsensusModeBayesian,
	})

	for _, nth := range []string{"1", "2"} {
		sLine := insertionColumnLine(simple, "3", nth)
		bLine := insertionColumnLine(bayes, "3", nth)
		if sLine == "" || bLine == "" {
			t.Fatalf("missing pos 3 nth %s insertion column\nsimple:\n%s\nbayes:\n%s",
				nth, simple, bayes)
		}
		sf := strings.Split(sLine, "\t")
		bf := strings.Split(bLine, "\t")
		sSeq, sQual := sf[6], sf[7]
		bSeq, bQual := bf[6], bf[7]

		// The deletion read is the 4th (last) in sorted order: seq "TTT*".
		if !strings.HasSuffix(bSeq, "*") {
			t.Fatalf("nth %s bayesian seq = %q, want a trailing '*' pad", nth, bSeq)
		}
		// (1) The pad qual must be the pre-gap Q20 ('5'), carried across nth.
		if got := bQual[len(bQual)-1]; got != '5' {
			t.Errorf("nth %s bayesian pad qual = %q (Q%d), want '5' (Q20). The "+
				"deletion-spanning pad must use the stateful running-min seeded "+
				"from the PRE-gap base (delPileupQual, seq_offset=readBP-2) and "+
				"carried across nth columns, not the raw post-gap quality.\n"+
				"line: %s", nth, string(got), int(got)-33, bLine)
		}
		// (2) Bayesian per-read seq/qual must equal the simple builder's — both
		// now consume buildInsertionColumnCells, so they cannot diverge.
		if bSeq != sSeq || bQual != sQual {
			t.Errorf("nth %s bayesian per-read columns differ from simple:\n"+
				"  bayes seq=%q qual=%q\n  simple seq=%q qual=%q", nth,
				bSeq, bQual, sSeq, sQual)
		}
	}
}
