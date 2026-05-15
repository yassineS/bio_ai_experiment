package samtools

import (
	"bytes"
	"strings"
	"testing"
)

// TestTargetcut_Basic exercises soft-clipped flanks, insertions
// (kept), deletions (no query bases consumed), and unmapped /
// secondary record skipping.
func TestTargetcut_Basic(t *testing.T) {
	// Reference is irrelevant for targetcut; we only need a header.
	// All quals are 'I' (Phred 40) so MinBaseQ default (13) keeps
	// every base unless explicitly tested.
	samText := strings.Join([]string{
		"@HD\tVN:1.6\tSO:coordinate",
		"@SQ\tSN:chr1\tLN:100",
		// r_basic: 10M with all-quality bases, full slice kept.
		"r_basic\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII",
		// r_softclip: 2S6M2S — soft-clipped flanks removed, 6 middle bases kept.
		"r_softclip\t0\tchr1\t5\t60\t2S6M2S\t*\t0\t0\tNNACGTACNN\tIIIIIIIIII",
		// r_insertion: 3M2I3M — insertions kept (query-consuming), 8 bases total.
		"r_insertion\t0\tchr1\t10\t60\t3M2I3M\t*\t0\t0\tACGTTACG\tIIIIIIII",
		// r_deletion: 3M2D3M — deletions consume no query, output is 6 bases (the 6 M bases).
		"r_deletion\t0\tchr1\t20\t60\t3M2D3M\t*\t0\t0\tACGTAC\tIIIIII",
		// r_unmapped: flag 0x4 — skipped.
		"r_unmapped\t4\t*\t0\t0\t*\t*\t0\t0\tACGTAC\tIIIIII",
		// r_secondary: flag 0x100 — skipped.
		"r_secondary\t256\tchr1\t30\t60\t6M\t*\t0\t0\tACGTAC\tIIIIII",
		// r_lowQ: a 6M record where all bases have qual '!'=0 — all dropped.
		"r_lowQ\t0\tchr1\t40\t60\t6M\t*\t0\t0\tACGTAC\t!!!!!!",
		// r_starseq: SEQ='*' — skipped (no sequence to emit).
		"r_starseq\t0\tchr1\t50\t60\t6M\t*\t0\t0\t*\t*",
	}, "\n") + "\n"

	var buf bytes.Buffer
	n, err := Targetcut(strings.NewReader(samText), &buf, TargetcutOptions{MinBaseQ: DefaultTargetcutMinBaseQ})
	if err != nil {
		t.Fatalf("Targetcut: %v", err)
	}

	// r_basic + r_softclip + r_insertion + r_deletion = 4 emitted records.
	// r_unmapped, r_secondary, r_lowQ, r_starseq → no output (lowQ becomes empty).
	if n != 4 {
		t.Fatalf("emitted = %d, want 4 (records emitted incl. r_basic/r_softclip/r_insertion/r_deletion). full output:\n%s", n, buf.String())
	}

	want := strings.Join([]string{
		">r_basic",
		"ACGTACGTAC",
		">r_softclip",
		"ACGTAC",
		">r_insertion",
		"ACGTTACG",
		">r_deletion",
		"ACGTAC",
		"",
	}, "\n")
	if got := buf.String(); got != want {
		t.Errorf("Targetcut output mismatch\ngot:\n%s\nwant:\n%s", got, want)
	}
}

// TestTargetcut_MinBaseQDropsBases checks per-base quality filtering.
// '!' is Phred 0, '"' is Phred 1, ..., 'I' is Phred 40. With MinBaseQ
// = 30, qual chars below '?'=30 are dropped. Mixed-quality reads
// should emit only the high-q subset.
func TestTargetcut_MinBaseQDropsBases(t *testing.T) {
	samText := strings.Join([]string{
		"@HD\tVN:1.6\tSO:coordinate",
		"@SQ\tSN:chr1\tLN:100",
		// 6M with alternating high/low quality.
		// Bases:   A C G T A C
		// Quals:   I ! I ! I !   (40 0 40 0 40 0)
		// At MinBaseQ=30 we keep A G A (positions 0,2,4).
		"r_mixed\t0\tchr1\t1\t60\t6M\t*\t0\t0\tACGTAC\tI!I!I!",
	}, "\n") + "\n"

	var buf bytes.Buffer
	if _, err := Targetcut(strings.NewReader(samText), &buf, TargetcutOptions{MinBaseQ: 30}); err != nil {
		t.Fatalf("Targetcut: %v", err)
	}
	want := ">r_mixed\nAGA\n"
	if got := buf.String(); got != want {
		t.Errorf("Targetcut min-baseq output:\ngot:  %q\nwant: %q", got, want)
	}
}

// TestTargetcut_DefaultExportedConst guards against drift in the
// upstream-default base-quality cutoff.
func TestTargetcut_DefaultExportedConst(t *testing.T) {
	if DefaultTargetcutMinBaseQ != 13 {
		t.Errorf("DefaultTargetcutMinBaseQ = %d, want 13 (upstream default)", DefaultTargetcutMinBaseQ)
	}
}
