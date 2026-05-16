package sickle

// Byte-for-byte parity tests for sickle's Go port against the upstream C
// reference implementation (najoshi/sickle v1.33). The upstream binary is
// built from reference_code/sickle (a git submodule pinned to v1.33 +
// upstream's post-1.33 commits) by running `make` in that directory. Each
// fixture under tools/sickle/testdata/parity/ was generated once by piping
// the matching input file through the upstream binary with the exact same
// command-line flags as the Go invocation below. The Go invocation drives
// our in-process library directly so the test runs without spawning a
// subprocess.
//
// Goal: every test that doesn't t.Skip should byte-match upstream's output.
// Any divergence is either:
//   - documented in tools/PARITY_VALIDATION.md (which the test points to),
//   - listed as an upstream bug in docs/UPSTREAM_BUGS.md and t.Skip'd here.

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/fastq"
)

// readParityFile reads a fixture from tools/sickle/testdata/parity/.
func readParityFile(t *testing.T, name string) []byte {
	t.Helper()
	path := filepath.Join("..", "..", "testdata", "parity", name)
	data, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("read parity fixture %s: %v", name, err)
	}
	return data
}

// runSickleSE runs the SE trimmer against an input fixture and returns the
// trimmed FASTQ bytes. Opens through bufio so the algorithm sees the same
// reader shape as the production CLI.
func runSickleSE(t *testing.T, inputName string, enc fastq.QualityEncoding, opts TrimOptions) []byte {
	t.Helper()
	in := readParityFile(t, inputName)
	var out bytes.Buffer
	if _, err := TrimSingleEnd(bufio.NewReader(bytes.NewReader(in)), &out, enc, opts); err != nil {
		t.Fatalf("TrimSingleEnd(%s) failed: %v", inputName, err)
	}
	return out.Bytes()
}

// runSicklePE runs the PE trimmer and returns (R1, R2, singletons) bytes.
// The singletons buffer is always provided so that the test exercises the
// upstream `-s` code path (which is the common case in real workflows).
func runSicklePE(t *testing.T, r1Name, r2Name string, enc fastq.QualityEncoding, opts TrimOptions) (r1Out, r2Out, sOut []byte) {
	t.Helper()
	r1 := readParityFile(t, r1Name)
	r2 := readParityFile(t, r2Name)
	var o1, o2, single bytes.Buffer
	if _, err := TrimPairedEnd(bufio.NewReader(bytes.NewReader(r1)), bufio.NewReader(bytes.NewReader(r2)),
		&o1, &o2, &single, enc, opts); err != nil {
		t.Fatalf("TrimPairedEnd failed: %v", err)
	}
	return o1.Bytes(), o2.Bytes(), single.Bytes()
}

func mustEqualBytes(t *testing.T, label string, got, want []byte) {
	t.Helper()
	if !bytes.Equal(got, want) {
		t.Fatalf("%s mismatch.\nwant:\n%s\ngot:\n%s", label, want, got)
	}
}

// case01 — SE basic -q 20 -l 20 sanger encoding.
// Mix of high-Q, gradient (high then low), and all-low-Q reads.
// Tests the core sliding-window trim + the no-5'-cut-found discard rule.
func TestParity_Sickle_Case01_SEBasic(t *testing.T) {
	opts := TrimOptions{QualThreshold: 20, LengthThreshold: 20}
	got := runSickleSE(t, "case01_se_basic.fq", fastq.Phred33, opts)
	want := readParityFile(t, "case01_se_basic.expected.fq")
	mustEqualBytes(t, "case01 SE basic", got, want)
}

// case02 — SE -n (truncate at first N).
// Reads with internal N, leading N, and clean reads.
func TestParity_Sickle_Case02_SETruncN(t *testing.T) {
	opts := TrimOptions{QualThreshold: 20, LengthThreshold: 10, TruncateN: true}
	got := runSickleSE(t, "case02_se_truncn.fq", fastq.Phred33, opts)
	want := readParityFile(t, "case02_se_truncn.expected.fq")
	mustEqualBytes(t, "case02 SE -n", got, want)
}

// case03 — SE -x (no 5' trimming).
// 5'-low / 3'-high read: with -x, the algorithm scans for 3' cut from the
// start; if the first window is below threshold, the 3' cut lands at 0 and
// the whole read is discarded — matching upstream behaviour.
func TestParity_Sickle_Case03_SENoFivePrime(t *testing.T) {
	opts := TrimOptions{QualThreshold: 20, LengthThreshold: 10, NoFivePrime: true}
	got := runSickleSE(t, "case03_se_nofiveprime.fq", fastq.Phred33, opts)
	want := readParityFile(t, "case03_se_nofiveprime.expected.fq")
	mustEqualBytes(t, "case03 SE -x", got, want)
}

// case04 — PE basic with -s singletons output.
func TestParity_Sickle_Case04_PEBasic(t *testing.T) {
	opts := TrimOptions{QualThreshold: 20, LengthThreshold: 10}
	r1, r2, s := runSicklePE(t, "case04_pe_r1.fq", "case04_pe_r2.fq", fastq.Phred33, opts)
	mustEqualBytes(t, "case04 PE R1", r1, readParityFile(t, "case04_pe_r1.expected.fq"))
	mustEqualBytes(t, "case04 PE R2", r2, readParityFile(t, "case04_pe_r2.expected.fq"))
	mustEqualBytes(t, "case04 PE singles", s, readParityFile(t, "case04_pe_s.expected.fq"))
}

// case05 — PE with singletons (one mate passes, the other fails).
func TestParity_Sickle_Case05_PESingles(t *testing.T) {
	opts := TrimOptions{QualThreshold: 20, LengthThreshold: 10}
	r1, r2, s := runSicklePE(t, "case05_pe_r1.fq", "case05_pe_r2.fq", fastq.Phred33, opts)
	mustEqualBytes(t, "case05 PE R1", r1, readParityFile(t, "case05_pe_r1.expected.fq"))
	mustEqualBytes(t, "case05 PE R2", r2, readParityFile(t, "case05_pe_r2.expected.fq"))
	mustEqualBytes(t, "case05 PE singles", s, readParityFile(t, "case05_pe_s.expected.fq"))
}

// case06 — SE illumina (Phred+64) encoding. Same trim semantics with the
// per-encoding quality-byte decode.
func TestParity_Sickle_Case06_SEIllumina(t *testing.T) {
	opts := TrimOptions{QualThreshold: 20, LengthThreshold: 10}
	got := runSickleSE(t, "case06_se_illumina.fq", fastq.Phred64, opts)
	want := readParityFile(t, "case06_se_illumina.expected.fq")
	mustEqualBytes(t, "case06 SE -t illumina", got, want)
}

// case07 — empty input. Upstream exits cleanly with an empty output; the
// Go port must do the same (no crash, no spurious header lines).
func TestParity_Sickle_Case07_SEEmpty(t *testing.T) {
	opts := TrimOptions{QualThreshold: 20, LengthThreshold: 20}
	got := runSickleSE(t, "case07_se_empty.fq", fastq.Phred33, opts)
	want := readParityFile(t, "case07_se_empty.expected.fq")
	mustEqualBytes(t, "case07 SE empty", got, want)
}

// case08 — all reads have quality strictly below threshold; upstream
// discards every one of them (no 5'-cut found in any).
func TestParity_Sickle_Case08_SEAllLowQual(t *testing.T) {
	opts := TrimOptions{QualThreshold: 20, LengthThreshold: 10}
	got := runSickleSE(t, "case08_se_alllow.fq", fastq.Phred33, opts)
	want := readParityFile(t, "case08_se_alllow.expected.fq")
	mustEqualBytes(t, "case08 SE all-low", got, want)
}

// case09 — boundary read at exactly Q=threshold and exactly Q=threshold-1.
// The exact-threshold read must pass; the just-below read must be discarded.
// Confirms the `>=` vs `<` boundary matches upstream.
func TestParity_Sickle_Case09_SEThreshold(t *testing.T) {
	opts := TrimOptions{QualThreshold: 20, LengthThreshold: 10}
	got := runSickleSE(t, "case09_se_threshold.fq", fastq.Phred33, opts)
	want := readParityFile(t, "case09_se_threshold.expected.fq")
	mustEqualBytes(t, "case09 SE threshold-boundary", got, want)
}

// case10 — gzipped input (upstream reads via zlib; our port via
// pkg/bioformats/iohelper). The TrimSingleEnd library entry point gets
// an io.Reader, so we manually wrap the gzip stream here to exercise the
// same code path.
func TestParity_Sickle_Case10_SEGzip(t *testing.T) {
	t.Helper()
	gzData := readParityFile(t, "case10_se_gz.fq.gz")
	gzr, err := gzip.NewReader(bytes.NewReader(gzData))
	if err != nil {
		t.Fatalf("gzip.NewReader: %v", err)
	}
	defer gzr.Close()
	var got bytes.Buffer
	opts := TrimOptions{QualThreshold: 20, LengthThreshold: 10}
	if _, err := TrimSingleEnd(bufio.NewReader(gzr), &got, fastq.Phred33, opts); err != nil {
		t.Fatalf("TrimSingleEnd(gz) failed: %v", err)
	}
	want := readParityFile(t, "case10_se_gz.expected.fq")
	mustEqualBytes(t, "case10 SE gzip", got.Bytes(), want)
}

// case11 — short reads below the length threshold. Upstream discards
// them up-front; the Go port replicates this.
func TestParity_Sickle_Case11_SEShort(t *testing.T) {
	opts := TrimOptions{QualThreshold: 20, LengthThreshold: 20}
	got := runSickleSE(t, "case11_se_short.fq", fastq.Phred33, opts)
	want := readParityFile(t, "case11_se_short.expected.fq")
	mustEqualBytes(t, "case11 SE short-read filter", got, want)
}

// case12 — PE with both mates passing/failing in a synced layout. The
// upstream PE writer emits R1 and R2 in lockstep and never crosses the
// streams; we should match.
func TestParity_Sickle_Case12_PESynced(t *testing.T) {
	opts := TrimOptions{QualThreshold: 20, LengthThreshold: 10}
	r1, r2, s := runSicklePE(t, "case12_pe_synced_r1.fq", "case12_pe_synced_r2.fq", fastq.Phred33, opts)
	mustEqualBytes(t, "case12 PE R1", r1, readParityFile(t, "case12_pe_synced_r1.expected.fq"))
	mustEqualBytes(t, "case12 PE R2", r2, readParityFile(t, "case12_pe_synced_r2.expected.fq"))
	mustEqualBytes(t, "case12 PE singles", s, readParityFile(t, "case12_pe_synced_s.expected.fq"))
}

// case13 — SE strict thresholds (-q 30 -l 5). High quality threshold forces
// the 3' tail of r1 to be cut at the first sub-Q30 window; r2 (uniform Q30)
// passes through whole; r3's leading N+homopolymer should pass with strict
// q30 because the average quality stays above 30.
func TestParity_Sickle_Case13_SEStrict(t *testing.T) {
	opts := TrimOptions{QualThreshold: 30, LengthThreshold: 5}
	got := runSickleSE(t, "case13_se_strict.fq", fastq.Phred33, opts)
	want := readParityFile(t, "case13_se_strict.expected.fq")
	mustEqualBytes(t, "case13 SE strict q30 l5", got, want)
}

// case14 — SE lax thresholds (-q 0 -l 0). All reads must pass through
// unmodified including the all-low-quality and the very short ones.
// Tests the boundary where the threshold is at 0 and the length filter
// is disabled.
func TestParity_Sickle_Case14_SELax(t *testing.T) {
	opts := TrimOptions{QualThreshold: 0, LengthThreshold: 0}
	got := runSickleSE(t, "case14_se_lax.fq", fastq.Phred33, opts)
	want := readParityFile(t, "case14_se_lax.expected.fq")
	mustEqualBytes(t, "case14 SE lax q0 l0", got, want)
}

// case15 — PE strict thresholds (-q 30 -l 10). r1 mate p2 has a low-quality
// 5' run that fails the threshold; the mate p2/2 still passes and lands in
// the singletons file. Tests strict-threshold PE singletons routing.
func TestParity_Sickle_Case15_PEStrict(t *testing.T) {
	opts := TrimOptions{QualThreshold: 30, LengthThreshold: 10}
	r1, r2, s := runSicklePE(t, "case15_pe_strict_r1.fq", "case15_pe_strict_r2.fq", fastq.Phred33, opts)
	mustEqualBytes(t, "case15 PE strict R1", r1, readParityFile(t, "case15_pe_strict_r1.expected.fq"))
	mustEqualBytes(t, "case15 PE strict R2", r2, readParityFile(t, "case15_pe_strict_r2.expected.fq"))
	mustEqualBytes(t, "case15 PE strict singles", s, readParityFile(t, "case15_pe_strict_s.expected.fq"))
}

// Smoke test: confirm fixtures exist under testdata/parity/ at the path the
// other tests expect; gives a single clear failure if the submodule wasn't
// initialised when the test corpus was regenerated.
func TestParity_Sickle_FixturesPresent(t *testing.T) {
	required := []string{
		"case01_se_basic.fq", "case01_se_basic.expected.fq",
		"case04_pe_r1.fq", "case04_pe_r1.expected.fq",
		"case10_se_gz.fq.gz",
		"case13_se_strict.fq", "case13_se_strict.expected.fq",
		"case14_se_lax.fq", "case14_se_lax.expected.fq",
		"case15_pe_strict_r1.fq", "case15_pe_strict_r1.expected.fq",
	}
	missing := []string{}
	for _, name := range required {
		p := filepath.Join("..", "..", "testdata", "parity", name)
		if _, err := os.Stat(p); err != nil {
			missing = append(missing, name)
		}
	}
	if len(missing) > 0 {
		t.Fatalf("missing parity fixtures (regenerate via reference_code/sickle/sickle):\n  %s", strings.Join(missing, "\n  "))
	}
}
