package samtools

// Live-binary oracle tests for the `samtools phase -b STR` per-haplotype
// BAM split (phase.c dump_aln). These run the genuine upstream samtools
// binary AND our locally-built port over the SAME sorted BAM, then decode
// each of the three output BAMs (<prefix>.0.bam, .1.bam, .chimera.bam)
// with the upstream `samtools view` and assert the decoded SAM records and
// headers are byte-identical, modulo @PG (which our port omits by design;
// both sides are invoked with --no-PG so no @PG is written either way).
//
// Determinism note: upstream phase.c routes reads through drand48() in
// dump_aln but never calls srand48(), so it runs from glibc's default seed
// state. Our drand48 port (phase_drand48.go) reproduces that exact
// sequence, which is what makes the .0/.1/.chimera assignment agree
// record-for-record rather than merely up to a 0<->1 relabelling. Per the
// project's testing rules the upstream check actually executes: the
// helpers t.Fatalf rather than t.Skip when the binaries cannot be produced.

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"
)

// phaseSplitFixture describes a phase -b oracle scenario: a SAM body and
// the phase CLI flags (beyond --no-PG -b PREFIX) to exercise.
type phaseSplitFixture struct {
	name string
	sam  string
	args []string // extra flags, e.g. {"-A"} or {"-F"}
}

// mkHetSeq builds a length-`total` read string of all 'A' with `alleles`
// placed at het positions 2, 7, 12, ... (0-based; every 5 bases).
func mkHetSeq(alleles string, total int) string {
	b := []byte(strings.Repeat("A", total))
	for i := 0; i < len(alleles); i++ {
		p := 2 + 5*i
		if p < total {
			b[p] = alleles[i]
		}
	}
	return string(b)
}

// buildPhaseSplitSAM assembles a coordinate-sorted SAM with two clean
// haplotype cohorts (all-G vs all-T at `nHet` het sites) plus one extra
// read whose allele string is `special` (e.g. a chimera or an ambiguous
// read). readLen must cover the het sites.
func buildPhaseSplitSAM(nHet, readLen, cohort int, special string) string {
	hapG := mkHetSeq(strings.Repeat("G", nHet), readLen)
	hapT := mkHetSeq(strings.Repeat("T", nHet), readLen)
	qual := strings.Repeat("I", readLen)
	cig := itoaM(readLen)
	var sb strings.Builder
	sb.WriteString("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:200\n")
	for i := 0; i < cohort; i++ {
		sb.WriteString("r_h0_" + string(rune('a'+i)) + "\t0\tchr1\t1\t60\t" + cig + "\t*\t0\t0\t" + hapG + "\t" + qual + "\n")
	}
	for i := 0; i < cohort; i++ {
		sb.WriteString("r_h1_" + string(rune('a'+i)) + "\t0\tchr1\t1\t60\t" + cig + "\t*\t0\t0\t" + hapT + "\t" + qual + "\n")
	}
	if special != "" {
		sp := mkHetSeq(special, readLen)
		sb.WriteString("r_special\t0\tchr1\t1\t60\t" + cig + "\t*\t0\t0\t" + sp + "\t" + qual + "\n")
	}
	return sb.String()
}

// itoaM renders "<n>M" without importing strconv at call sites.
func itoaM(n int) string {
	if n == 0 {
		return "0M"
	}
	var digits []byte
	for n > 0 {
		digits = append([]byte{byte('0' + n%10)}, digits...)
		n /= 10
	}
	return string(digits) + "M"
}

// TestLivePhaseBamSplit drives several phase -b scenarios against the real
// upstream binary and asserts byte-identical decoded output for all three
// buckets. It is the parity gate for the -b / -A / -F BAM-split path.
func TestLivePhaseBamSplit(t *testing.T) {
	live := upstreamSamtools(t)
	ours := ourSamtoolsBinary(t)

	fixtures := []phaseSplitFixture{
		{
			// Two clean cohorts; no extra read. Exercises the confident
			// haplotype routing (phase.c which=0/1 plus the is_flip
			// shuffle) — the canonical phase fixture.
			name: "clean",
			sam:  buildPhaseSplitSAM(2, 10, 3, ""),
		},
		{
			// A symmetric 4G/4T chimera over 8 hets. fragphase flips it
			// (flip=1) but in==out keeps phased=0, so dump_aln routes it
			// via which=3 (random) — NOT the chimera bucket. Stresses the
			// phased==0 random path interacting with the flip rewrite.
			name: "chimera_symmetric",
			sam:  buildPhaseSplitSAM(8, 40, 4, strings.Repeat("G", 4)+strings.Repeat("T", 4)),
		},
		{
			// A 4G/5T chimera over 9 hets: fragphase flips it AND keeps
			// in!=out, so phased=1 && flip=1 → dump_aln which=2 → the
			// chimera bucket fires. This is the case that distinguishes
			// the default from -F.
			name: "chimera_repaired",
			sam:  buildPhaseSplitSAM(9, 100, 6, strings.Repeat("G", 4)+strings.Repeat("T", 5)),
		},
		{
			// Same chimera_repaired input but with -F: chimera repair is
			// off, so the flip never happens, which=2 never fires and the
			// read lands in a haplotype bucket. Proves -F changes -b
			// output.
			name: "chimera_repaired_F",
			sam:  buildPhaseSplitSAM(9, 100, 6, strings.Repeat("G", 4)+strings.Repeat("T", 5)),
			args: []string{"-F"},
		},
		{
			// An ambiguous read: 3G/2T over 5 hets gives in=3, out=2 →
			// f.ambig=1. Default routes it via which=3 (random); with -A
			// it is dropped to the chimera bucket (which=2). This pair of
			// scenarios proves -A changes -b output.
			name: "ambiguous_default",
			sam:  buildPhaseSplitSAM(5, 40, 4, "GGTGT"),
		},
		{
			name: "ambiguous_drop",
			sam:  buildPhaseSplitSAM(5, 40, 4, "GGTGT"),
			args: []string{"-A"},
		},
	}

	for _, fx := range fixtures {
		fx := fx
		t.Run(fx.name, func(t *testing.T) {
			dir := t.TempDir()
			samPath := filepath.Join(dir, "in.sam")
			if err := os.WriteFile(samPath, []byte(fx.sam), 0o644); err != nil {
				t.Fatal(err)
			}
			// Build the BAM with OUR view so both binaries read identical
			// input bytes.
			bamPath := filepath.Join(dir, "in.bam")
			if err := os.WriteFile(bamPath, runSamtools(t, ours, "view", "-b", "--no-PG", samPath), 0o644); err != nil {
				t.Fatal(err)
			}

			upPrefix := filepath.Join(dir, "up")
			ourPrefix := filepath.Join(dir, "ours")
			upArgs := append([]string{"phase", "--no-PG", "-b", upPrefix}, fx.args...)
			ourArgs := append([]string{"phase", "--no-PG", "-b", ourPrefix}, fx.args...)
			upArgs = append(upArgs, bamPath)
			ourArgs = append(ourArgs, bamPath)
			// Run both (TSV stdout is discarded here; the TSV stream is
			// covered by TestLivePhase).
			runSamtools(t, live, upArgs...)
			runSamtools(t, ours, ourArgs...)

			for _, suffix := range []string{"0", "1", "chimera"} {
				upBam := upPrefix + "." + suffix + ".bam"
				ourBam := ourPrefix + "." + suffix + ".bam"
				upRecs := decodeSAMRecords(t, live, upBam)
				ourRecs := decodeSAMRecords(t, live, ourBam)
				if !bytes.Equal(upRecs, ourRecs) {
					t.Errorf("bucket %s.bam records differ:\nupstream:\n%s\nours:\n%s",
						suffix, upRecs, ourRecs)
				}
				upHdr := decodeSAMHeaderNoPG(t, live, upBam)
				ourHdr := decodeSAMHeaderNoPG(t, live, ourBam)
				if !bytes.Equal(upHdr, ourHdr) {
					t.Errorf("bucket %s.bam headers differ (modulo @PG):\nupstream:\n%s\nours:\n%s",
						suffix, upHdr, ourHdr)
				}
			}
		})
	}
}

// decodeSAMRecords returns the decoded SAM alignment records (no header)
// for a BAM file, via the upstream `samtools view`.
func decodeSAMRecords(t *testing.T, bin, bam string) []byte {
	t.Helper()
	return runSamtools(t, bin, "view", bam)
}

// decodeSAMHeaderNoPG returns the decoded SAM header for a BAM file with
// @PG lines stripped (our port never injects @PG; the comparison must be
// modulo @PG to be meaningful).
func decodeSAMHeaderNoPG(t *testing.T, bin, bam string) []byte {
	t.Helper()
	hdr := runSamtools(t, bin, "view", "-H", bam)
	var keep [][]byte
	for _, line := range bytes.Split(hdr, []byte("\n")) {
		if bytes.HasPrefix(line, []byte("@PG")) {
			continue
		}
		keep = append(keep, line)
	}
	return bytes.Join(keep, []byte("\n"))
}
