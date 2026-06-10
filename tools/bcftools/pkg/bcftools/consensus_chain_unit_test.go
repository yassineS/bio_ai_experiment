package bcftools

import (
	"bytes"
	"os"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
)

// TestChainPushAndPrint exercises the pushGap / print routines directly,
// checking the UCSC chain rendering for a single insertion gap. The
// reference is 10 bp; an insertion of 2 bp at ref position 3 splits the
// alignment into two ungapped blocks.
func TestChainPushAndPrint(t *testing.T) {
	c := NewChain(0)
	// Insertion: 0 ref bases consumed, 2 alt bases inserted at ref pos 3.
	c.pushGap(3, 0, 3, 2)
	var buf bytes.Buffer
	if err := c.print(&buf, "chr1", 10, 1); err != nil {
		t.Fatalf("print: %v", err)
	}
	// score = block0 (3) + last block (7) = 10; alt length = 12.
	want := "chain 10 chr1 10 + 0 10 chr1 12 + 0 12 1\n" +
		"3 0 2\n" +
		"7\n\n"
	if buf.String() != want {
		t.Fatalf("chain mismatch:\n got %q\n want %q", buf.String(), want)
	}
}

// TestChainBackToBackMerge checks that two abutting gaps are merged into a
// single block, mirroring push_chain_gap's back-to-back branch.
func TestChainBackToBackMerge(t *testing.T) {
	c := NewChain(0)
	c.pushGap(3, 1, 3, 0) // deletion of 1 ref base at pos 3
	c.pushGap(4, 1, 3, 0) // immediately abutting deletion of 1 ref base
	if len(c.blockLengths) != 1 {
		t.Fatalf("expected back-to-back merge into 1 block, got %d blocks", len(c.blockLengths))
	}
	if c.refGaps[0] != 2 || c.altGaps[0] != 0 {
		t.Fatalf("merged gap = (ref %d, alt %d), want (2, 0)", c.refGaps[0], c.altGaps[0])
	}
}

// TestConsensusChainInsertionDeletion drives the full Consensus path with a
// chain writer and verifies both the FASTA and the in-memory chain content
// for a deterministic insertion + deletion case (no upstream binary).
func TestConsensusChainInsertionDeletion(t *testing.T) {
	dir := t.TempDir()
	chainPath := dir + "/out.chain"
	ref := []*fasta.Record{{ID: "chr1", Sequence: []byte("AAAACCCCGGGGTTTT")}}
	body := "##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=16>\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chr1\t5\t.\tC\tCAA\t.\tPASS\t.\n" // insertion of AA after pos 5
	var out bytes.Buffer
	if _, err := Consensus(strings.NewReader(body), &out, ConsensusOptions{
		Reference: ref,
		LineWidth: 80,
		ChainFile: chainPath,
	}); err != nil {
		t.Fatalf("Consensus: %v", err)
	}
	if got, want := out.String(), ">chr1\nAAAACAACCCGGGGTTTT\n"; got != want {
		t.Fatalf("fasta:\n got %q\n want %q", got, want)
	}
	chainBytes, err := os.ReadFile(chainPath)
	if err != nil {
		t.Fatalf("read chain: %v", err)
	}
	chain := string(chainBytes)
	// Insertion shares the base before the event, so the block extends by 1
	// (push at ref pos 5, ref gap 0, alt gap 2). ref length 16, alt 18.
	want := "chain 16 chr1 16 + 0 16 chr1 18 + 0 18 1\n" +
		"5 0 2\n" +
		"11\n\n"
	if chain != want {
		t.Fatalf("chain:\n got %q\n want %q", chain, want)
	}
}

// TestConsensusNoChainWhenUnset confirms no chain file is created when
// ChainFile is empty.
func TestConsensusNoChainWhenUnset(t *testing.T) {
	ref := []*fasta.Record{{ID: "chr1", Sequence: []byte("ACGT")}}
	body := "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t1\t.\tA\tT\t.\t.\t.\n"
	var out bytes.Buffer
	if _, err := Consensus(strings.NewReader(body), &out, ConsensusOptions{Reference: ref, LineWidth: 80}); err != nil {
		t.Fatalf("Consensus: %v", err)
	}
	if out.String() != ">chr1\nTCGT\n" {
		t.Fatalf("unexpected fasta: %q", out.String())
	}
}

// TestPhasedIUPACSelection covers the NpIu logic at the library level: a
// phased genotype takes the N-th haplotype index; an unphased genotype
// collapses to an IUPAC ambiguity code.
func TestPhasedIUPACSelection(t *testing.T) {
	ref := []*fasta.Record{{ID: "chr1", Sequence: []byte("ACGT")}}
	cases := []struct {
		name string
		gt   string
		idx  int
		want string // resulting consensus base at pos 1
	}{
		{"phased-hap1", "0|1", 1, "A"},     // phased, index 1 -> REF
		{"phased-hap2", "0|1", 2, "C"},     // phased, index 2 -> ALT
		{"unphased-iupac", "0/1", 1, "M"},  // unphased A/C -> M
		{"unphased-iupac2", "0/1", 2, "M"}, // index irrelevant when unphased
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			body := "##fileformat=VCFv4.2\n" +
				"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
				"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n" +
				"chr1\t1\t.\tA\tC\t.\t.\t.\tGT\t" + tc.gt + "\n"
			var out bytes.Buffer
			if _, err := Consensus(strings.NewReader(body), &out, ConsensusOptions{
				Reference:      ref,
				Sample:         "S1",
				Haplotype:      HapPhasedIUPAC,
				HaplotypeIndex: tc.idx,
				LineWidth:      80,
			}); err != nil {
				t.Fatalf("Consensus: %v", err)
			}
			got := strings.TrimSpace(strings.SplitN(out.String(), "\n", 2)[1])
			if got[:1] != tc.want {
				t.Fatalf("pos1 base = %q, want %q (full %q)", got[:1], tc.want, got)
			}
		})
	}
}

// TestIUPACBitmaskRoundTrip checks the IUPAC encode/decode helpers against
// the canonical nucleotide ambiguity table.
func TestIUPACBitmaskRoundTrip(t *testing.T) {
	cases := []struct {
		ch   byte
		mask byte
	}{
		{'A', 1}, {'C', 2}, {'G', 4}, {'T', 8},
		{'M', 3}, {'R', 5}, {'W', 9}, {'S', 6}, {'Y', 10}, {'K', 12},
		{'V', 7}, {'H', 11}, {'D', 13}, {'B', 14}, {'N', 15},
	}
	for _, tc := range cases {
		if got := iupac2bitmask(tc.ch); got != int(tc.mask) {
			t.Errorf("iupac2bitmask(%c) = %d, want %d", tc.ch, got, tc.mask)
		}
		// Lowercase must map identically.
		if got := iupac2bitmask(tc.ch + 32); got != int(tc.mask) {
			t.Errorf("iupac2bitmask(%c) = %d, want %d", tc.ch+32, got, tc.mask)
		}
		if got := bitmask2iupac(tc.mask); got != tc.ch {
			t.Errorf("bitmask2iupac(%d) = %c, want %c", tc.mask, got, tc.ch)
		}
	}
	if iupac2bitmask('-') != -1 {
		t.Errorf("iupac2bitmask('-') should be -1")
	}
	if bitmask2iupac(0) != 0 || bitmask2iupac(16) != 0 {
		t.Errorf("bitmask2iupac out-of-range should be 0")
	}
}

// TestGTIsPhased exercises the phasing detector.
func TestGTIsPhased(t *testing.T) {
	cases := []struct {
		gt   string
		want bool
	}{
		{"0|1", true},
		{"1|1", true},
		{"0/1", false},
		{"1", true}, // haploid is trivially phased
		{".", true}, // haploid missing
		{"./.", false},
	}
	for _, tc := range cases {
		if got := gtIsPhased(tc.gt); got != tc.want {
			t.Errorf("gtIsPhased(%q) = %v, want %v", tc.gt, got, tc.want)
		}
	}
}
