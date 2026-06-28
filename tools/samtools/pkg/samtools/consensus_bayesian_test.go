package samtools

import (
	"strings"
	"testing"
)

// TestApplyMDCosts_ConsecutiveSubsDoNotSlide is the unit-level regression for
// task #53. Upstream's nm_init MD loop (bam_consensus.c:1174-1201) advances
// only the md-string cursor after a substitution (md++); it does NOT advance
// the query offset `pos`. Consecutive substitutions therefore stack their
// halo bands at one centre rather than sliding one base per substitution.
//
// The Go port previously carried a spurious `pos++` after each substitution,
// which slid every band right by one base per extra substitution. For the
// real-data pattern MD:Z:0T0G0A0A106 (four substitutions sharing query offset
// 0) with the default halo of 50, that over-counted local_nm in the outer
// band: the flat plateau (offsets 50..99 == +20) acquired a sliding ramp
// (offsets 50..53 == 35,30,25,20) and bled +5/+10/+15 into offsets 100..102
// that must be zero.
//
// This test pins the byte-exact band shape upstream produces. It FAILS with
// the spurious `pos++` (e.g. localNM[50] becomes 35, localNM[100] becomes 15)
// and PASSES once it is removed.
func TestApplyMDCosts_ConsecutiveSubsDoNotSlide(t *testing.T) {
	const halo = 50
	const qlen = 160
	nm := make([]int32, qlen)
	// Four consecutive substitutions all at query offset 0, then 106 matches.
	applyMDCosts(nm, "0T0G0A0A106", halo, qlen)

	low24 := func(i int) int { return int(nm[i] & ((1 << 24) - 1)) }

	// Inner band (offsets 0..49): four subs each contribute +10 -> 40, FLAT.
	for i := 0; i < halo; i++ {
		if got := low24(i); got != 40 {
			t.Fatalf("inner band: localNM[%d] = %d, want 40 (flat). A sliding "+
				"ramp here means the spurious pos++ regressed.", i, got)
		}
	}
	// Outer band (offsets 50..99): four subs each contribute +5 -> 20, FLAT.
	// With the bug these slide to 35,30,25,20 at offsets 50,51,52,53.
	for i := halo; i < 2*halo; i++ {
		if got := low24(i); got != 20 {
			t.Fatalf("outer band: localNM[%d] = %d, want 20 (flat plateau). "+
				"The spurious pos++ slides this to 35/30/25/... at the left edge.", i, got)
		}
	}
	// Beyond 2*halo (offsets >= 100): no band reaches here, must be zero.
	// With the bug, offsets 100/101/102 acquire 15/10/5.
	for _, i := range []int{100, 101, 102, 103, 150} {
		if got := low24(i); got != 0 {
			t.Fatalf("localNM[%d] = %d, want 0. The spurious pos++ leaks the "+
				"outer band past 2*halo (15/10/5 at 100/101/102).", i, got)
		}
	}
}

// islandEdgeMDSAM is the synthetic end-to-end regression fixture for task #53.
// A single 110bp read aligns flush at chr1:1 carrying MD:Z:0T0G0A0A106 — four
// consecutive substitutions at its left edge, exactly the real contig-20
// pattern (the proxy 20:160157 read) that the diagnose pinned. With per-base
// quality Q11 the bayesian caller sits right on the consensus cutoff (default
// -C 10) across the read, so the MD halo's extra cost is decisive: with the
// spurious pos++ the over-counted local_nm depresses the adjusted MAPQ enough
// to mask the coverage-island's first callable bases to N; without it those
// bases are called.
//
// The read sequence is "ACGT" repeated to 110bp; with the fix the island
// becomes callable starting at query offset 49 ("...NNNNCGTACG..."), whereas
// the buggy pos++ pushes the first callable base three positions right
// ("...NNNNNNNACG..."), masking C/G/T at offsets 49/50/51 to N.
const islandEdgeMDSAM = "@HD\tVN:1.6\tSO:coordinate\n" +
	"@SQ\tSN:chr1\tLN:200\n" +
	"r0\t0\tchr1\t1\t60\t110M\t*\t0\t0\t" +
	"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\t" +
	",,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,\t" +
	"MD:Z:0T0G0A0A106\n"

// TestConsensus_Bayesian_MDIsland_CallsEdgeBase is the end-to-end regression
// for task #53: the bayesian -f fasta caller must call the coverage-island's
// first callable bases rather than masking them to N. The MD halo over-count
// from the spurious pos++ slid the island edge three bases to the right; with
// the fix the edge base (query offset 49, a 'C') is called. This test FAILS
// with the spurious pos++ (offsets 49/50/51 render N) and PASSES without it.
func TestConsensus_Bayesian_MDIsland_CallsEdgeBase(t *testing.T) {
	out := runConsensusOnSAM(t, islandEdgeMDSAM, ConsensusOptions{
		Format: ConsensusFASTA,
		Mode:   ConsensusModeBayesian,
		// UseMQual / ConsCutoff / nm_halo left at their upstream defaults
		// (true / 10 / 50): the defaults are what exposed the bug in the
		// whole-contig-20 run.
	})
	// Unwrap the FASTA body (strip the header and join wrapped lines).
	lines := strings.Split(strings.TrimRight(out, "\n"), "\n")
	if len(lines) < 2 || !strings.HasPrefix(lines[0], ">") {
		t.Fatalf("unexpected FASTA shape:\n%s", out)
	}
	body := strings.Join(lines[1:], "")
	if len(body) < 52 {
		t.Fatalf("consensus body too short (%d):\n%s", len(body), out)
	}
	// The island's first callable bases are at query offsets 49,50,51 == C,G,T.
	// With the spurious pos++ these are masked to N.
	for off, want := range map[int]byte{49: 'C', 50: 'G', 51: 'T'} {
		if got := body[off]; got != want {
			t.Errorf("island-edge offset %d = %q, want %q (a base, not N): the "+
				"MD-halo off-by-one (spurious pos++) masked the coverage-island "+
				"edge.\nbody[40:55]=%q", off, string(got), string(want), body[40:55])
		}
	}
	// Guard: the buggy output starts the island at offset 52, so offset 51
	// being a base (not N) is the load-bearing discriminator.
	if body[51] == 'N' {
		t.Errorf("offset 51 is N: the island edge is still slid right by the "+
			"MD-halo off-by-one.\nbody[40:55]=%q", body[40:55])
	}
}
