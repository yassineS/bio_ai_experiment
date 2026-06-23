package skewer

import (
	"bufio"
	"bytes"
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
		// PEBlockSize left 0: keep every pair (the final-block drop is a
		// drop-in-CLI-only behaviour, validated separately in the real-data run).
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

// TestPE_BlockTruncation reproduces upstream's final-block drop: with a block
// size of 3 and 5 input pairs, only the first 3 (one complete block) survive;
// the trailing partial block of 2 is silently dropped.
func TestPE_BlockTruncation(t *testing.T) {
	var r1b, r2b bytes.Buffer
	for i := 0; i < 5; i++ {
		// Non-overlapping reads that pass through untrimmed, so the only effect
		// under test is the block-boundary drop.
		r1b.WriteString("@r" + string(rune('A'+i)) + " 1\n")
		r1b.WriteString("ACGTACGTACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n")
		r2b.WriteString("@r" + string(rune('A'+i)) + " 2\n")
		r2b.WriteString("TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAA\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n")
	}

	opts := peOpts()
	opts.PEBlockSize = 3

	var o1, o2 bytes.Buffer
	if _, err := TrimPairedEnd(bufio.NewReader(bytes.NewReader(r1b.Bytes())),
		bufio.NewReader(bytes.NewReader(r2b.Bytes())), &o1, &o2, nil, fastq.Phred33, opts); err != nil {
		t.Fatalf("TrimPairedEnd failed: %v", err)
	}

	count1 := bytes.Count(o1.Bytes(), []byte("\n")) / 4
	count2 := bytes.Count(o2.Bytes(), []byte("\n")) / 4
	if count1 != 3 || count2 != 3 {
		t.Errorf("block truncation: want 3 records per mate, got pair1=%d pair2=%d", count1, count2)
	}
}
