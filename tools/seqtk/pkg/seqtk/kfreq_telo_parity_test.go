package seqtk

// Byte-for-byte parity tests for `seqtk kfreq` and `seqtk telo`
// against the upstream C reference implementation (reference_code/seqtk
// v1.5-r133). Fixtures live under tools/seqtk/testdata/parity/ and were
// produced by running the pinned upstream binary with the same input
// files and flags as the Go invocations below:
//
//	reference_code/seqtk/seqtk kfreq <kmer> kfreq_small.fa > kfreq_small_<kmer>.expected.txt
//	reference_code/seqtk/seqtk telo  [flags] telo_basic.fa  > telo_basic_<tag>.stdout.expected.txt 2> ...stderr...
//
// For telo we capture stdout and stderr separately because upstream
// writes the BED rows to stdout and the summary line to stderr (and
// the Go port preserves that split).

import (
	"bytes"
	"testing"
)

// runKfreq drives Kfreq and returns its bytes.
func runKfreq(t *testing.T, inName, kmer string) []byte {
	t.Helper()
	in := readParityFile(t, inName)
	var out bytes.Buffer
	if err := Kfreq(bytes.NewReader(in), &out, KfreqOptions{Kmer: kmer}); err != nil {
		t.Fatalf("Kfreq(%s, %s): %v", inName, kmer, err)
	}
	return out.Bytes()
}

// TestParity_Seqtk_Kfreq_Small_AC covers the canonical case from the
// upstream usage docs: a 2-mer on a short fixture where forward and
// reverse-complement neighbour counts diverge between the two records.
func TestParity_Seqtk_Kfreq_Small_AC(t *testing.T) {
	got := runKfreq(t, "kfreq_small.fa", "AC")
	want := readParityFile(t, "kfreq_small_AC.expected.txt")
	mustEqualBytes(t, "kfreq AC (kfreq_small.fa)", got, want)
}

// TestParity_Seqtk_Kfreq_Small_ACGT exercises a 4-mer (longer rolling
// window) and confirms the "exact-match-of-target wins over neighbour"
// branch upstream's `else if` chain encodes.
func TestParity_Seqtk_Kfreq_Small_ACGT(t *testing.T) {
	got := runKfreq(t, "kfreq_small.fa", "ACGT")
	want := readParityFile(t, "kfreq_small_ACGT.expected.txt")
	mustEqualBytes(t, "kfreq ACGT (kfreq_small.fa)", got, want)
}

// TestParity_Seqtk_Kfreq_Small_AAAA tests the all-zero-counts path,
// pinning the upstream tie-break that picks '-' when both neighbour
// counts are equal (including the all-zero case).
func TestParity_Seqtk_Kfreq_Small_AAAA(t *testing.T) {
	got := runKfreq(t, "kfreq_small.fa", "AAAA")
	want := readParityFile(t, "kfreq_small_AAAA.expected.txt")
	mustEqualBytes(t, "kfreq AAAA (kfreq_small.fa)", got, want)
}

// TestParity_Seqtk_Kfreq_Edge exercises a zero-length record + a
// length-4 record where the kmer is longer than the sequence on the
// short side; verifies upstream's silent emit of the zero-counts row.
func TestParity_Seqtk_Kfreq_Edge(t *testing.T) {
	got := runKfreq(t, "kfreq_edge.fa", "AC")
	want := readParityFile(t, "kfreq_edge_AC.expected.txt")
	mustEqualBytes(t, "kfreq AC (kfreq_edge.fa)", got, want)
}

// TestParity_Seqtk_Kfreq_Mixed_AA covers a multi-record fixture that
// includes a repetitive "AA" record (max forward hits), a random
// record, a mixed-case record (lowercase ACGT must count via
// seq_nt6_table), and an N-containing record (the rolling-window
// reset path).
func TestParity_Seqtk_Kfreq_Mixed_AA(t *testing.T) {
	got := runKfreq(t, "kfreq_mixed.fa", "AA")
	want := readParityFile(t, "kfreq_mixed_AA.expected.txt")
	mustEqualBytes(t, "kfreq AA (kfreq_mixed.fa)", got, want)
}

// TestParity_Seqtk_Kfreq_Mixed_ACGT pins the 4-mer behaviour on the
// mixed fixture, in particular the N-reset path (kfreq must not span
// the NNNN gap).
func TestParity_Seqtk_Kfreq_Mixed_ACGT(t *testing.T) {
	got := runKfreq(t, "kfreq_mixed.fa", "ACGT")
	want := readParityFile(t, "kfreq_mixed_ACGT.expected.txt")
	mustEqualBytes(t, "kfreq ACGT (kfreq_mixed.fa)", got, want)
}

// TestParity_Seqtk_Kfreq_Mixed_CCGG hits a kmer whose only matches are
// inside the mixed-case record — pins case-folding parity with
// upstream's seq_nt6_table.
func TestParity_Seqtk_Kfreq_Mixed_CCGG(t *testing.T) {
	got := runKfreq(t, "kfreq_mixed.fa", "CCGG")
	want := readParityFile(t, "kfreq_mixed_CCGG.expected.txt")
	mustEqualBytes(t, "kfreq CCGG (kfreq_mixed.fa)", got, want)
}

// TestParity_Seqtk_Kfreq_Mixed_CCCTAA tests a 6-mer; the larger rolling
// window also exercises the neighbour-set sizing (1 << 12 = 4096
// entries).
func TestParity_Seqtk_Kfreq_Mixed_CCCTAA(t *testing.T) {
	got := runKfreq(t, "kfreq_mixed.fa", "CCCTAA")
	want := readParityFile(t, "kfreq_mixed_CCCTAA.expected.txt")
	mustEqualBytes(t, "kfreq CCCTAA (kfreq_mixed.fa)", got, want)
}

// runTelo drives Telo and returns (stdoutBytes, stderrBytes).
func runTelo(t *testing.T, inName string, opts TeloOptions) ([]byte, []byte) {
	t.Helper()
	in := readParityFile(t, inName)
	var out, errb bytes.Buffer
	if err := Telo(bytes.NewReader(in), &out, &errb, opts); err != nil {
		t.Fatalf("Telo(%s): %v", inName, err)
	}
	return out.Bytes(), errb.Bytes()
}

// teloDefaults returns the upstream-default TeloOptions; helper to
// keep the parity-test call sites readable.
func teloDefaults() TeloOptions {
	return TeloOptions{
		Motif:    DefaultTeloMotif,
		Penalty:  DefaultTeloPenalty,
		MaxDrop:  DefaultTeloMaxDrop,
		MinScore: DefaultTeloMinScore,
	}
}

// TestParity_Seqtk_Telo_BasicDefault covers the headline case: 5' +
// 3' BED rows on a synthetic chromosome, with the summary line on
// stderr.
func TestParity_Seqtk_Telo_BasicDefault(t *testing.T) {
	gotStdout, gotStderr := runTelo(t, "telo_basic.fa", teloDefaults())
	mustEqualBytes(t, "telo default stdout (telo_basic.fa)",
		gotStdout, readParityFile(t, "telo_basic_default.stdout.expected.txt"))
	mustEqualBytes(t, "telo default stderr (telo_basic.fa)",
		gotStderr, readParityFile(t, "telo_basic_default.stderr.expected.txt"))
}

// TestParity_Seqtk_Telo_BasicMotifTTAGGG flips the motif. Because the
// telomere we built has CCCTAA at 5' and TTAGGG at 3', passing
// -m TTAGGG should produce no BED rows from the 5' scan; the 3' scan
// reverse-encodes so it also fails. We pin stdout to its (empty)
// expected file and stderr to the bytes-with-zero-telomere summary.
func TestParity_Seqtk_Telo_BasicMotifTTAGGG(t *testing.T) {
	opts := teloDefaults()
	opts.Motif = "TTAGGG"
	gotStdout, gotStderr := runTelo(t, "telo_basic.fa", opts)
	mustEqualBytes(t, "telo -m TTAGGG stdout (telo_basic.fa)",
		gotStdout, readParityFile(t, "telo_basic_mTTAGGG.stdout.expected.txt"))
	mustEqualBytes(t, "telo -m TTAGGG stderr (telo_basic.fa)",
		gotStderr, readParityFile(t, "telo_basic_mTTAGGG.stderr.expected.txt"))
}

// TestParity_Seqtk_Telo_BasicProfile drives the -P show-profile path
// and asserts the byte-for-byte match of every P / Q row. With
// -s 0 every record contributes (the BED-row guard is the
// `if (max >= min_score)` check after the loop).
func TestParity_Seqtk_Telo_BasicProfile(t *testing.T) {
	opts := teloDefaults()
	opts.ShowProfile = true
	opts.MinScore = 0
	gotStdout, gotStderr := runTelo(t, "telo_basic.fa", opts)
	mustEqualBytes(t, "telo -P -s 0 stdout (telo_basic.fa)",
		gotStdout, readParityFile(t, "telo_basic_P_s0.stdout.expected.txt"))
	mustEqualBytes(t, "telo -P -s 0 stderr (telo_basic.fa)",
		gotStderr, readParityFile(t, "telo_basic_P_s0.stderr.expected.txt"))
}

// TestParity_Seqtk_Telo_ComplexDefault exercises a multi-record
// fixture with non-telomere records, N-runs in the middle, and a
// mixed-case lowercase-only record. Upstream's seq_nt6_table treats
// lowercase identically to uppercase, so the lowercase record's
// behaviour is pinned by the expected fixture.
func TestParity_Seqtk_Telo_ComplexDefault(t *testing.T) {
	gotStdout, gotStderr := runTelo(t, "telo_complex.fa", teloDefaults())
	mustEqualBytes(t, "telo default stdout (telo_complex.fa)",
		gotStdout, readParityFile(t, "telo_complex_default.stdout.expected.txt"))
	mustEqualBytes(t, "telo default stderr (telo_complex.fa)",
		gotStderr, readParityFile(t, "telo_complex_default.stderr.expected.txt"))
}

// TestParity_Seqtk_Telo_ComplexMinScore100 verifies the -s flag by
// lowering the min-score threshold. This pins the
// `if (max >= min_score)` guard against upstream.
func TestParity_Seqtk_Telo_ComplexMinScore100(t *testing.T) {
	opts := teloDefaults()
	opts.MinScore = 100
	gotStdout, gotStderr := runTelo(t, "telo_complex.fa", opts)
	mustEqualBytes(t, "telo -s 100 stdout (telo_complex.fa)",
		gotStdout, readParityFile(t, "telo_complex_s100.stdout.expected.txt"))
	mustEqualBytes(t, "telo -s 100 stderr (telo_complex.fa)",
		gotStderr, readParityFile(t, "telo_complex_s100.stderr.expected.txt"))
}

// TestParity_Seqtk_Telo_ComplexPenaltyMaxDrop tunes both -p and -d.
// Penalty changes the per-miss score adjustment; max-drop changes the
// abort threshold of the scan. Together they shift the BED-row
// boundaries vs the default run.
func TestParity_Seqtk_Telo_ComplexPenaltyMaxDrop(t *testing.T) {
	opts := teloDefaults()
	opts.Penalty = 2
	opts.MaxDrop = 500
	gotStdout, gotStderr := runTelo(t, "telo_complex.fa", opts)
	mustEqualBytes(t, "telo -p 2 -d 500 stdout (telo_complex.fa)",
		gotStdout, readParityFile(t, "telo_complex_p2d500.stdout.expected.txt"))
	mustEqualBytes(t, "telo -p 2 -d 500 stderr (telo_complex.fa)",
		gotStderr, readParityFile(t, "telo_complex_p2d500.stderr.expected.txt"))
}

// TestParity_Seqtk_Telo_EdgeCases covers the zero-length-record /
// length-4 record path. No BED rows expected; stderr summary is
// "0\t4\n".
func TestParity_Seqtk_Telo_EdgeCases(t *testing.T) {
	gotStdout, gotStderr := runTelo(t, "telo_edge.fa", teloDefaults())
	mustEqualBytes(t, "telo default stdout (telo_edge.fa)",
		gotStdout, readParityFile(t, "telo_edge_default.stdout.expected.txt"))
	mustEqualBytes(t, "telo default stderr (telo_edge.fa)",
		gotStderr, readParityFile(t, "telo_edge_default.stderr.expected.txt"))
}
