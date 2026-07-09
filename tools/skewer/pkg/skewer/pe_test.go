package skewer

import (
	"bufio"
	"bytes"
	"fmt"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fastq"
)

// The default Illumina pair-1/pair-2 3' adapter prefixes upstream applies for
// `-m pe` (parameter.cpp:48,50). Kept here so the PE regression tests exercise
// the same adapters the drop-in CLI configures.
const (
	testPEAdapter1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
	testPEAdapter2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA"
)

func peOpts() TrimOptions {
	return TrimOptions{
		Adapter3:      testPEAdapter1,
		Adapter3Pair2: testPEAdapter2,
		MinLength:     18,
		MinOverlap:    3,
		ErrorRate:     0.1,
		IndelRate:     0.03,
		PEMatrixMode:  true,
	}
}

func runPEPair(t *testing.T, r1, r2 string) (out1, out2 string) {
	t.Helper()
	var o1, o2 bytes.Buffer
	if _, err := TrimPairedEnd(bufio.NewReader(bytes.NewReader([]byte(r1))),
		bufio.NewReader(bytes.NewReader([]byte(r2))), &o1, &o2, nil, fastq.Phred33, peOpts()); err != nil {
		t.Fatalf("TrimPairedEnd failed: %v", err)
	}
	return o1.String(), o2.String()
}

// TestPE_OverlapCorrection_RealRecord is a byte-for-byte regression against
// upstream skewer 0.2.2 for one real GIAB exome read pair
// (NIST7035 NA12878, read 1101:12564:45475). The pair overlaps with a short
// insert, so upstream:
//   - cuts both mates from 101 bp to 99 bp at the overlap-derived position, and
//   - error-corrects read 1 over the overlap (combinePairSeqs): the mid-read
//     base A→C and the 3' half of read 1's quality string is rewritten with the
//     higher-quality complement from read 2.
//
// Read 2's sequence is only truncated (never corrected). The expected outputs
// below were produced by the upstream binary (reference_code/skewer/skewer
// -m pe) on exactly this pair.
func TestPE_OverlapCorrection_RealRecord(t *testing.T) {
	r1 := "@HWI-D00119:50:H7AP8ADXX:1:1101:12564:45475 1:N:0:TAAGGCGA\n" +
		"CGCCTGCACTGCGTCTGGTGTTTCATTTCCTCCCGGCCTCCTCGCCTGCAATGGGTCTGGTGTTTCATTTCCTCCCGGCCTCCTCGCCTGCACTGCGTCTG\n" +
		"+\n" +
		"CCCFFFFFHHHHHJJJJJHHHIJJJJJJJJJJJJIJJJJJJJJJJJJJJI;CHHHHHHFFFFFCEDEEDDDEDDDBBBDDBDDDDDDDDAB?CCDCBBDBB\n"
	r2 := "@HWI-D00119:50:H7AP8ADXX:1:1101:12564:45475 2:N:0:TAAGGCGA\n" +
		"GACGCAGTGCAGGCGAGGAGGCCGGGAGGAAATGAAACACCAGACCCAGTGCAGGCGAGGAGGCCGGGAGGAAATGAAACACCAGACCCAGTGCAGGCGAG\n" +
		"+\n" +
		"CCCFFFFFHHHHHJJJJJJJJJIIJJIIJGJJJIJJJIJHHHHFFBEDDEEEDDDDD;;<7?@DB;;;@>>0?>:>:@>>9A<8<<288A2@>>CC?####\n"

	wantR1 := "@HWI-D00119:50:H7AP8ADXX:1:1101:12564:45475 1:N:0:TAAGGCGA\n" +
		"CGCCTGCACTGCGTCTGGTGTTTCATTTCCTCCCGGCCTCCTCGCCTGCACTGGGTCTGGTGTTTCATTTCCTCCCGGCCTCCTCGCCTGCACTGCGTC\n" +
		"+\n" +
		"CCCFFFFFHHHHHJJJJJHHHIJJJJJJJJJJJJIJJJJJJJJJJJJJJIDDHHHHHHHHJIJJJIJJJGJIIJJIIJJJJJJJJJHHHHHFFFFFCCD\n"
	wantR2 := "@HWI-D00119:50:H7AP8ADXX:1:1101:12564:45475 2:N:0:TAAGGCGA\n" +
		"GACGCAGTGCAGGCGAGGAGGCCGGGAGGAAATGAAACACCAGACCCAGTGCAGGCGAGGAGGCCGGGAGGAAATGAAACACCAGACCCAGTGCAGGCG\n" +
		"+\n" +
		"CCCFFFFFHHHHHJJJJJJJJJIIJJIIJGJJJIJJJIJHHHHFFBEDDEEEDDDDD;;<7?@DB;;;@>>0?>:>:@>>9A<8<<288A2@>>CC?##\n"

	gotR1, gotR2 := runPEPair(t, r1, r2)
	if gotR1 != wantR1 {
		t.Errorf("pair1 mismatch.\nwant:\n%s\ngot:\n%s", wantR1, gotR1)
	}
	if gotR2 != wantR2 {
		t.Errorf("pair2 mismatch.\nwant:\n%s\ngot:\n%s", wantR2, gotR2)
	}
}

// TestPE_NonOverlapping_PassThrough confirms a genuinely non-overlapping pair
// (no adapter, prefixes are not reverse-complements) passes through untrimmed
// and uncorrected — the PE overlap analysis must not invent a cut.
func TestPE_NonOverlapping_PassThrough(t *testing.T) {
	r1 := "@p1 1\n" +
		"GAATAAGGGAATCTGCTGGAAGTCTTTATTTTTAAAAAACACCAAAAGTGGGAAGAAATCAGTGATGAGGGATTTATGAACCAGGATTTTGTGGATTTAAG\n" +
		"+\n" +
		"@@@FFFFFGFHHGIIIIJJGIJIJHHJJJJJJJJJJJIJJJJIJJJJIFHIJJIJJIIIJIICHIHHHHHHFEFFEEEDDEDBDD?ADDDDACBCDDDDEC\n"
	r2 := "@p1 2\n" +
		"GTTACAAATTATCACAATGCTTAAATCCACAAAATCCTGGTTCATAAATCCCTCATCACTGATTTCTTCCCACTTTTGGTGTTTTTTAAAAATAAAGACTT\n" +
		"+\n" +
		"CCCFFFFFGFGHHJJJFIIIJJJIJJIIIIJJJJJJHGGIGHIJIGIJIIIIIJIHIIIJJJJIIIGCHIJJIIJJIIHFHCDDFDDDDDDCBDCDCCCDD\n"

	gotR1, gotR2 := runPEPair(t, r1, r2)
	if gotR1 != r1 {
		t.Errorf("pair1 should pass through unchanged.\nwant:\n%s\ngot:\n%s", r1, gotR1)
	}
	if gotR2 != r2 {
		t.Errorf("pair2 should pass through unchanged.\nwant:\n%s\ngot:\n%s", r2, gotR2)
	}
}

// synthPEPair builds a deterministic n-pair PE FASTQ. About half the pairs
// overlap (so they exercise the trim + error-correction path) and half pass
// through, giving the determinism test a mix of outcomes.
func synthPEPair(n int) (r1, r2 []byte) {
	var r1b, r2b bytes.Buffer
	// One overlapping template (trims 101->99, corrects R1) reused from the
	// real-record regression, plus a non-overlapping template.
	ovR1seq := "CGCCTGCACTGCGTCTGGTGTTTCATTTCCTCCCGGCCTCCTCGCCTGCAATGGGTCTGGTGTTTCATTTCCTCCCGGCCTCCTCGCCTGCACTGCGTCTG"
	ovR1q := "CCCFFFFFHHHHHJJJJJHHHIJJJJJJJJJJJJIJJJJJJJJJJJJJJI;CHHHHHHFFFFFCEDEEDDDEDDDBBBDDBDDDDDDDDAB?CCDCBBDBB"
	ovR2seq := "GACGCAGTGCAGGCGAGGAGGCCGGGAGGAAATGAAACACCAGACCCAGTGCAGGCGAGGAGGCCGGGAGGAAATGAAACACCAGACCCAGTGCAGGCGAG"
	ovR2q := "CCCFFFFFHHHHHJJJJJJJJJIIJJIIJGJJJIJJJIJHHHHFFBEDDEEEDDDDD;;<7?@DB;;;@>>0?>:>:@>>9A<8<<288A2@>>CC?####"
	npR1seq := "GAATAAGGGAATCTGCTGGAAGTCTTTATTTTTAAAAAACACCAAAAGTGGGAAGAAATCAGTGATGAGGGATTTATGAACCAGGATTTTGTGGATTTAAG"
	npR1q := "@@@FFFFFGFHHGIIIIJJGIJIJHHJJJJJJJJJJJIJJJJIJJJJIFHIJJIJJIIIJIICHIHHHHHHFEFFEEEDDEDBDD?ADDDDACBCDDDDEC"
	npR2seq := "GTTACAAATTATCACAATGCTTAAATCCACAAAATCCTGGTTCATAAATCCCTCATCACTGATTTCTTCCCACTTTTGGTGTTTTTTAAAAATAAAGACTT"
	npR2q := "CCCFFFFFGFGHHJJJFIIIJJJIJJIIIIJJJJJJHGGIGHIJIGIJIIIIIJIHIIIJJJJIIIGCHIJJIIJJIIHFHCDDFDDDDDDCBDCDCCCDD"
	for i := 0; i < n; i++ {
		tag := string(rune('A'+i%26)) + string(rune('a'+(i/26)%26))
		if i%2 == 0 {
			fmt.Fprintf(&r1b, "@ov%s 1\n%s\n+\n%s\n", tag, ovR1seq, ovR1q)
			fmt.Fprintf(&r2b, "@ov%s 2\n%s\n+\n%s\n", tag, ovR2seq, ovR2q)
		} else {
			fmt.Fprintf(&r1b, "@np%s 1\n%s\n+\n%s\n", tag, npR1seq, npR1q)
			fmt.Fprintf(&r2b, "@np%s 2\n%s\n+\n%s\n", tag, npR2seq, npR2q)
		}
	}
	return r1b.Bytes(), r2b.Bytes()
}

func runPE(t *testing.T, r1, r2 []byte, threads int) (out1, out2 string, st *TrimStats) {
	t.Helper()
	opts := peOpts()
	opts.Threads = threads
	var o1, o2 bytes.Buffer
	st, err := TrimPairedEnd(bufio.NewReader(bytes.NewReader(r1)),
		bufio.NewReader(bytes.NewReader(r2)), &o1, &o2, nil, fastq.Phred33, opts)
	if err != nil {
		t.Fatalf("TrimPairedEnd(threads=%d) failed: %v", threads, err)
	}
	return o1.String(), o2.String(), st
}

// TestPE_KeepsAllPairs confirms FIX A: every input pair is emitted (no
// final-block drop). With 5 non-overlapping, passing pairs all 5 survive.
func TestPE_KeepsAllPairs(t *testing.T) {
	var r1b, r2b bytes.Buffer
	for i := 0; i < 5; i++ {
		r1b.WriteString("@r" + string(rune('A'+i)) + " 1\n")
		r1b.WriteString("ACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n")
		r2b.WriteString("@r" + string(rune('A'+i)) + " 2\n")
		r2b.WriteString("TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAA\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n")
	}
	o1, o2, _ := runPE(t, r1b.Bytes(), r2b.Bytes(), 1)
	count1 := strings.Count(o1, "\n") / 4
	count2 := strings.Count(o2, "\n") / 4
	if count1 != 5 || count2 != 5 {
		t.Errorf("keep-all: want 5 records per mate, got pair1=%d pair2=%d", count1, count2)
	}
}

// TestPE_ThreadInvariant is the determinism contract: TrimPairedEnd output must
// be byte-identical across thread counts (and equal the sequential path), and
// the merged stats must be order-independent.
func TestPE_ThreadInvariant(t *testing.T) {
	r1, r2 := synthPEPair(500)
	base1, base2, baseStats := runPE(t, r1, r2, 1)
	for _, threads := range []int{2, 4, 8} {
		o1, o2, st := runPE(t, r1, r2, threads)
		if o1 != base1 {
			t.Errorf("threads=%d: R1 output differs from -t1", threads)
		}
		if o2 != base2 {
			t.Errorf("threads=%d: R2 output differs from -t1", threads)
		}
		if st.TotalReads != baseStats.TotalReads || st.TrimmedReads != baseStats.TrimmedReads ||
			st.TrimmedBases != baseStats.TrimmedBases || st.DiscardedReads != baseStats.DiscardedReads ||
			st.TotalBases != baseStats.TotalBases || st.AdapterFound3 != baseStats.AdapterFound3 {
			t.Errorf("threads=%d: stats differ from -t1.\nwant %+v\ngot  %+v", threads, baseStats, st)
		}
	}
}
