package samtools

// Complex-input live-oracle test for `samtools phase`.
//
// TestLivePhase (phase_live_oracle_test.go) pins byte-parity on a tiny
// two-het fixture. That fixture does not exercise the deterministic
// phasing engine hard enough: it has no chimeric reads, no FL-masked
// region, and too few reads at any one het index to trigger a khash
// table grow. This file builds a deliberately COMPLEX fixture — ~10
// heterozygous SNP sites, many reads spanning several sites, and a
// batch of chimeric reads that switch haplotype mid-read — and asserts
// that the full `samtools phase` text stream (CC banner + PS / FL / M /
// EV / //) matches the upstream binary BYTE-FOR-BYTE.
//
// Upstream phase is fully deterministic on this path: phase() computes
// an int8 `path` via the dynaprog Viterbi recurrence, fragphase()
// repairs chimeras, genmask() emits FL regions, and EV lines are dumped
// in ks_introsort_rseq order over the khash buckets. There is no MCMC
// and no RNG in this path (the only drand48 use is the -b BAM split,
// not exercised here). The port reproduces every one of those steps,
// including the in-place Cuckoo-style khash rehash that fixes the
// equal-vpos EV tie order after the fragment table grows.

import (
	"bytes"
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// complexPhaseRefAlt is the per-position (ref, alt) allele table used to
// synthesise the complex fixture. Positions are 0-based on chr1.
var complexPhaseRefAlt = []struct {
	pos      int
	ref, alt byte
}{
	{9, 'A', 'C'}, {19, 'G', 'T'}, {29, 'A', 'G'}, {39, 'C', 'T'},
	{49, 'A', 'T'}, {59, 'G', 'C'}, {69, 'A', 'C'}, {79, 'G', 'A'},
	{89, 'C', 'A'}, {99, 'T', 'G'},
}

// buildComplexPhaseSAM returns a coordinate-sorted SAM string with a
// reference length of 120, ten het sites, sixteen clean reads (eight on
// each haplotype) and ten strong chimeric reads that switch haplotype
// in the middle. The chimeras carry >= 3 supporting hets on each side
// so they trip fragphase's FLIP_THRES repair, and the read pile is dense
// enough to grow the khash fragment table past its initial 16 buckets —
// the exact regime where a clean (non-kick-out) rehash would reorder the
// equal-vpos EV lines relative to upstream.
func buildComplexPhaseSAM() string {
	const L = 120
	refbase := make([]byte, L)
	for i := range refbase {
		refbase[i] = 'A'
	}
	for _, ra := range complexPhaseRefAlt {
		refbase[ra.pos] = ra.ref
	}
	refseq := string(refbase)

	altAt := func(gp int) (byte, byte, bool) {
		for _, ra := range complexPhaseRefAlt {
			if ra.pos == gp {
				return ra.ref, ra.alt, true
			}
		}
		return 0, 0, false
	}
	// readSeq builds the read substring [start,start+length) on the given
	// haplotype; when switch >= 0 the read flips to the other haplotype
	// at reference positions >= switch (a chimera).
	readSeq := func(start, length, hap, switchPos int) string {
		s := []byte(refseq[start : start+length])
		for j := 0; j < length; j++ {
			gp := start + j
			if r, a, ok := altAt(gp); ok {
				h := hap
				if switchPos >= 0 && gp >= switchPos {
					h = 1 - hap
				}
				if h == 0 {
					s[j] = r
				} else {
					s[j] = a
				}
			}
		}
		return string(s)
	}

	type spec struct{ start, length, hap, switchPos int }
	var specs []spec
	// Clean reads at staggered starts, both haplotypes.
	for _, st := range []int{0, 0, 0, 5, 5, 10, 10, 15, 20, 20, 25, 30, 30, 35, 40, 40} {
		specs = append(specs, spec{st, 80, 0, -1})
		specs = append(specs, spec{st, 80, 1, -1})
	}
	// Strong chimeras at the same locus, switching mid-read.
	for i := 0; i < 5; i++ {
		specs = append(specs, spec{0, 90, 0, 50})
		specs = append(specs, spec{0, 90, 1, 40})
	}

	type rec struct {
		start int
		line  string
	}
	recs := make([]rec, 0, len(specs))
	for i, sp := range specs {
		length := sp.length
		if sp.start+length > L {
			length = L - sp.start
		}
		seq := readSeq(sp.start, length, sp.hap, sp.switchPos)
		qual := strings.Repeat("I", length)
		recs = append(recs, rec{sp.start, fmt.Sprintf(
			"r%d\t0\tchr1\t%d\t60\t%dM\t*\t0\t0\t%s\t%s", i, sp.start+1, length, seq, qual)})
	}
	// Stable sort by start (coordinate order).
	for i := 1; i < len(recs); i++ {
		for j := i; j > 0 && recs[j].start < recs[j-1].start; j-- {
			recs[j], recs[j-1] = recs[j-1], recs[j]
		}
	}

	var b strings.Builder
	b.WriteString("@HD\tVN:1.6\tSO:coordinate\n")
	fmt.Fprintf(&b, "@SQ\tSN:chr1\tLN:%d\n", L)
	for _, r := range recs {
		b.WriteString(r.line)
		b.WriteByte('\n')
	}
	return b.String()
}

// TestLivePhaseComplex asserts byte-for-byte equality of the full
// `samtools phase --no-PG` stream between the upstream binary and our
// port on the complex (chimera + FL + table-grow) fixture. This is the
// parity gate that pins the deterministic-DP phasing engine — including
// the FL-masked regions, the fragphase chimera flips (YF:i:1 reads), and
// the EV emit order over equal-vpos fragments — against real upstream.
func TestLivePhaseComplex(t *testing.T) {
	live := upstreamSamtools(t)
	ours := ourSamtoolsBinary(t)

	dir := t.TempDir()
	samPath := filepath.Join(dir, "complex.sam")
	if err := os.WriteFile(samPath, []byte(buildComplexPhaseSAM()), 0o644); err != nil {
		t.Fatal(err)
	}
	bamPath := filepath.Join(dir, "complex.bam")
	if err := os.WriteFile(bamPath, runSamtools(t, ours, "view", "-b", "--no-PG", samPath), 0o644); err != nil {
		t.Fatal(err)
	}

	up := runSamtools(t, live, "phase", "--no-PG", bamPath)
	up2 := runSamtools(t, live, "phase", "--no-PG", bamPath)
	if !bytes.Equal(up, up2) {
		t.Fatalf("upstream phase output non-deterministic across runs")
	}

	// Sanity: the fixture must actually exercise the hard paths, else the
	// test is not pinning what it claims. Require at least one FL line, at
	// least one flipped (YF:i:1) chimera, and a het index high enough to
	// have grown the khash table.
	if !bytes.Contains(up, []byte("\nFL\t")) {
		t.Fatalf("fixture did not produce an FL (masked) region:\n%s", up)
	}
	if !bytes.Contains(up, []byte("YF:i:1")) {
		t.Fatalf("fixture did not produce a fragphase chimera flip (YF:i:1):\n%s", up)
	}

	gp := runSamtools(t, ours, "phase", "--no-PG", bamPath)
	if !bytes.Equal(up, gp) {
		t.Errorf("DIVERGENCE: complex phase byte-stream differs\nupstream (%d bytes):\n%s\nours (%d bytes):\n%s",
			len(up), up, len(gp), gp)
	}
}
