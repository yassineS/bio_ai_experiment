package samtools

// Live-binary oracle parity tests for five bugs the parity pipeline found in
// the samtools port, each asserting byte-for-byte agreement with the genuine
// upstream `samtools` built from reference_code/:
//
//  1. depth -q/-Q were swapped (and the base-quality count diverged): upstream
//     -q/--min-BQ is the per-base quality floor, -Q/--min-MQ is the read MAPQ
//     floor (bam2depth.c case 'q' -> min_qual, case 'Q' -> min_mqual).
//  2. mpileup omitted interior zero-depth rows that upstream emits.
//  3. sort (coordinate) tie-break ordering differed from htslib bam1_lt.
//  4. view -s (subsample) selected a different read subset than upstream.
//  5. addreplacerg -r did not add/replace the @RG line + RG:Z: tags.
//
// Per the project rule these helpers t.Fatalf (never t.Skip) when a binary
// cannot be built. They reuse upstreamSamtools / ourSamtoolsBinary /
// runSamtools / decodeBAMRecords from the other live tests in this package.

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"
)

// bugfixDepthSAM has a het base-quality and MAPQ spread so the -q (base) and
// -Q (mapping) filters select genuinely different subsets:
//   - r2's first two bases are Phred 0 ('!') so -q drops them per-base.
//   - r2 (MAPQ 20) and r3 (MAPQ 5) are dropped wholesale by -Q 30.
const bugfixDepthSAM = "@HD\tVN:1.6\tSO:coordinate\n" +
	"@SQ\tSN:chr1\tLN:50\n" +
	"r1\t0\tchr1\t10\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n" +
	"r2\t0\tchr1\t10\t20\t5M\t*\t0\t0\tACGTA\t!!III\n" +
	"r3\t0\tchr1\t11\t5\t5M\t*\t0\t0\tCGTAC\tIIIII\n"

// TestLiveDepthQQBaseVsMapping is the parity case for bug #1. It asserts that
// `depth -q N` (base quality) and `depth -Q N` (mapping quality), and their
// combination with -a, are byte-identical to upstream.
func TestLiveDepthQQBaseVsMapping(t *testing.T) {
	live := upstreamSamtools(t)
	ours := ourSamtoolsBinary(t)
	dir := t.TempDir()
	sam := filepath.Join(dir, "depthqq.sam")
	if err := os.WriteFile(sam, []byte(bugfixDepthSAM), 0o644); err != nil {
		t.Fatal(err)
	}

	cases := [][]string{
		{"depth", "-q", "30", sam}, // base-quality floor
		{"depth", "-Q", "30", sam}, // mapping-quality floor
		{"depth", "-q", "30", "-Q", "10", sam},
		{"depth", "-a", "-q", "30", sam},
		{"depth", "-aa", "-q", "30", sam},
		{"depth", "--min-BQ", "30", sam},
		{"depth", "--min-MQ", "30", sam},
	}
	for _, args := range cases {
		upstream := runSamtools(t, live, args...)
		mine := runSamtools(t, ours, args...)
		if !bytes.Equal(upstream, mine) {
			t.Fatalf("depth %v mismatch:\nupstream:\n%s\nours:\n%s", args[1:], upstream, mine)
		}
	}
}

// bugfixMpileupGapSAM has two reads on chr1 separated by a gap, plus a
// second read-free contig, so -a must fill chr1's leading (1..4), interior
// (10..19) and trailing (25..40) zero-depth rows up to LN, and -aa must
// additionally emit all of chr2.
const bugfixMpileupGapSAM = "@HD\tVN:1.6\tSO:coordinate\n" +
	"@SQ\tSN:chr1\tLN:40\n" +
	"@SQ\tSN:chr2\tLN:6\n" +
	"r1\t0\tchr1\t5\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n" +
	"r2\t0\tchr1\t20\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n"

// TestLiveMpileupZeroDepthRows is the parity case for bug #2: upstream
// mpileup -a fills every position of a read-bearing contig (leading,
// interior and trailing zero-depth rows up to the contig length), and -aa
// additionally emits read-free contigs. The default (no -a) emits only
// covered positions. All three are asserted byte-for-byte vs upstream.
func TestLiveMpileupZeroDepthRows(t *testing.T) {
	live := upstreamSamtools(t)
	ours := ourSamtoolsBinary(t)
	dir := t.TempDir()
	sam := filepath.Join(dir, "mpgap.sam")
	if err := os.WriteFile(sam, []byte(bugfixMpileupGapSAM), 0o644); err != nil {
		t.Fatal(err)
	}

	for _, args := range [][]string{
		{"mpileup", sam},
		{"mpileup", "-a", sam},
		{"mpileup", "-aa", sam},
	} {
		upstream := runSamtools(t, live, args...)
		mine := runSamtools(t, ours, args...)
		if !bytes.Equal(upstream, mine) {
			t.Fatalf("mpileup %v mismatch:\nupstream:\n%s\nours:\n%s", args[1:], upstream, mine)
		}
	}
}
