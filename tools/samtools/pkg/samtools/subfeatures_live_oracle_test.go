package samtools

// Live-binary oracle tests for the samtools per-subcommand feature gaps closed
// in this change set:
//
//   - coverage -A / -m / -D / -w : the ASCII / UTF histogram output mode.
//   - markdup  -s / -f          : the duplicate-statistics report.
//   - markdup  -d               : optical-duplicate detection (dt:Z:SQ/LB).
//   - markdup  -S               : supplementary / secondary duplicate marking.
//   - calmd    -e / -u / -C     : the previously-deferred calmd knobs.
//
// Each test runs the genuine upstream `samtools` (built from the vendored
// reference_code submodule, shared via the upstreamSamtools helper) alongside
// the Go port on the same fixture and asserts byte-for-byte agreement. BAM
// output is decoded back to SAM text (with the upstream binary) so the @PG
// provenance line the Go port intentionally omits does not perturb the
// comparison. Per the project rules these helpers t.Fatalf — never t.Skip —
// when the upstream binary cannot be produced.

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"
)

// runSamtoolsSubfeatures runs bin with args and returns (stdout, stderr). It
// does not fail the test on a non-zero exit so the caller can compare warning
// output and partial results; callers that require success assert on the
// returned streams. Uniquely named to avoid clashing with the shared
// runSamtools helper, which fails on any non-zero exit.
func runSamtoolsSubfeatures(t *testing.T, bin string, args ...string) (stdout, stderr []byte) {
	t.Helper()
	cmd := exec.Command(bin, args...)
	var out, errb bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = &errb
	_ = cmd.Run()
	return out.Bytes(), errb.Bytes()
}

// writeSubfeatureBAM writes sam to dir/in.sam and converts it to dir/in.bam
// with the upstream binary, returning the BAM path.
func writeSubfeatureBAM(t *testing.T, live, dir, sam string) string {
	t.Helper()
	samPath := filepath.Join(dir, "in.sam")
	if err := os.WriteFile(samPath, []byte(sam), 0o644); err != nil {
		t.Fatal(err)
	}
	bamPath := filepath.Join(dir, "in.bam")
	if out, err := exec.Command(live, "view", "-b", "-o", bamPath, samPath).CombinedOutput(); err != nil {
		t.Fatalf("upstream view -b failed: %v\n%s", err, out)
	}
	return bamPath
}

// viewSAMRecords decodes a BAM byte slice to SAM text using the upstream
// binary and returns only the alignment lines (header dropped) so @PG / @HD
// provenance differences are ignored.
func viewSAMRecords(t *testing.T, live string, bam []byte) string {
	t.Helper()
	cmd := exec.Command(live, "view", "-")
	cmd.Stdin = bytes.NewReader(bam)
	out, err := cmd.CombinedOutput()
	if err != nil {
		t.Fatalf("upstream view decode failed: %v\n%s", err, out)
	}
	return string(out)
}

// coverageHistFixtureSAM has two references with a couple of overlapping reads,
// one reverse-strand read, and an unmapped read, exercising the histogram
// breadth/depth binning, the per-bin block selection, and the side-panel
// statistics.
const coverageHistFixtureSAM = "@HD\tVN:1.6\tSO:coordinate\n" +
	"@SQ\tSN:chr1\tLN:120\n" +
	"@SQ\tSN:chr2\tLN:50\n" +
	"r1\t0\tchr1\t10\t60\t20M\t*\t0\t0\tACGTACGTACGTACGTACGT\tIIIIIIIIIIIIIIIIIIII\n" +
	"r2\t0\tchr1\t15\t40\t20M\t*\t0\t0\tACGTACGTACGTACGTACGT\tIIIIIIIIIIIIIIIIIIII\n" +
	"r3\t16\tchr1\t50\t30\t15M\t*\t0\t0\tCGTACGTACGTACGT\tIIIIIIIIIIIIIII\n" +
	"r4\t0\tchr2\t5\t25\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n" +
	"r5\t4\t*\t0\t0\t*\t*\t0\t0\tACGTA\tIIIII\n"

// TestLiveCoverageHistogram asserts byte equality of every histogram output
// mode (-m, -A, -D, -w N) plus the now-upstream-matching tabular form.
func TestLiveCoverageHistogram(t *testing.T) {
	live := upstreamSamtools(t)
	ours := ourSamtoolsBinary(t)
	dir := t.TempDir()
	bam := writeSubfeatureBAM(t, live, dir, coverageHistFixtureSAM)

	cases := [][]string{
		{"coverage", bam},             // tabular
		{"coverage", "-H", bam},       // tabular, no header
		{"coverage", "-m", bam},       // UTF histogram
		{"coverage", "-A", bam},       // ASCII histogram
		{"coverage", "-D", bam},       // depth-plot histogram
		{"coverage", "-D", "-A", bam}, // ASCII depth-plot
		{"coverage", "-w", "20", bam}, // 20-bin histogram
		{"coverage", "-Q", "20", bam}, // base-quality filter (tabular)
		{"coverage", "-q", "30", bam}, // mapq filter (tabular)
		{"coverage", "--min-depth", "2", "-m", bam},
		// -D plot-depth with a >1 min-depth: upstream sums depth at every
		// position with depth>=1 (before the mindepth gate), so the depth
		// bars differ from the breadth gate.
		{"coverage", "-D", "--min-depth", "2", bam},
	}
	for _, args := range cases {
		// COLUMNS unset so upstream's terminal-width probe yields the 40-bin
		// fallback the Go port assumes.
		up, _ := runSamtoolsSubfeatures(t, live, args...)
		mine, _ := runSamtoolsSubfeatures(t, ours, args...)
		if !bytes.Equal(up, mine) {
			t.Fatalf("coverage %v mismatch:\nupstream:\n%s\nours:\n%s", args[1:], up, mine)
		}
	}
}

// coverageEmptyRefFixtureSAM lists two references but only aligns reads to the
// first; upstream's histogram mode visits only references with pileup columns,
// so chr2 must be omitted from the histogram (but still listed in the table).
const coverageEmptyRefFixtureSAM = "@HD\tVN:1.6\tSO:coordinate\n" +
	"@SQ\tSN:chr1\tLN:60\n" +
	"@SQ\tSN:chr2\tLN:50\n" +
	"r1\t0\tchr1\t10\t60\t20M\t*\t0\t0\tACGTACGTACGTACGTACGT\tIIIIIIIIIIIIIIIIIIII\n"

// TestLiveCoverageEmptyRefHistogram pins that the histogram skips references
// with no selected reads while the tabular form still lists them.
func TestLiveCoverageEmptyRefHistogram(t *testing.T) {
	live := upstreamSamtools(t)
	ours := ourSamtoolsBinary(t)
	dir := t.TempDir()
	bam := writeSubfeatureBAM(t, live, dir, coverageEmptyRefFixtureSAM)

	for _, args := range [][]string{
		{"coverage", "-m", bam},
		{"coverage", "-A", bam},
		{"coverage", bam}, // tabular: chr2 still listed
	} {
		up, _ := runSamtoolsSubfeatures(t, live, args...)
		mine, _ := runSamtoolsSubfeatures(t, ours, args...)
		if !bytes.Equal(up, mine) {
			t.Fatalf("coverage %v empty-ref mismatch:\nupstream:\n%s\nours:\n%s", args[1:], up, mine)
		}
	}
}

// markdupExcludedFixtureSAM mixes a normal duplicate pair, a primary record
// that already carries the duplicate flag (FDUP, must be EXAMINED not
// EXCLUDED), and a trailing unmapped read (must be EXCLUDED). It pins the
// EXCLUDED / EXAMINED / SINGLE counters.
const markdupExcludedFixtureSAM = "@HD\tVN:1.6\tSO:coordinate\n" +
	"@SQ\tSN:chr1\tLN:1000\n" +
	"a\t99\tchr1\t10\t60\t10M\t=\t100\t100\tACGTACGTAC\tIIIIIIIIII\tms:i:200\tMC:Z:10M\n" +
	"p1\t1024\tchr1\t50\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n" +
	"a\t147\tchr1\t100\t60\t10M\t=\t10\t-100\tACGTACGTAC\tIIIIIIIIII\tms:i:200\tMC:Z:10M\n" +
	"u1\t4\t*\t0\t0\t*\t*\t0\t0\tACGTA\tIIIII\n"

// TestLiveMarkdupExcludedCounters asserts the EXCLUDED/EXAMINED accounting
// matches upstream when the input has a pre-marked duplicate primary and an
// unmapped read (upstream's exclude mask is SECONDARY|SUPP|UNMAP|QCFAIL).
func TestLiveMarkdupExcludedCounters(t *testing.T) {
	live := upstreamSamtools(t)
	ours := ourSamtoolsBinary(t)
	dir := t.TempDir()
	bam := writeSubfeatureBAM(t, live, dir, markdupExcludedFixtureSAM)

	args := []string{"markdup", "-s", bam, filepath.Join(dir, "out.bam")}
	_, upErr := runSamtoolsSubfeatures(t, live, args...)
	_, myErr := runSamtoolsSubfeatures(t, ours, args...)
	if statsBody(upErr) != statsBody(myErr) {
		t.Fatalf("markdup -s EXCLUDED/EXAMINED mismatch:\nupstream:\n%s\nours:\n%s", statsBody(upErr), statsBody(myErr))
	}
}

// markdupStatsFixtureSAM has a duplicate pair (a / b at the same coordinates),
// a unique pair (c), and a supplementary alignment of the duplicate that
// carries an SA tag so the -S path has a record to mark.
const markdupStatsFixtureSAM = "@HD\tVN:1.6\tSO:coordinate\n" +
	"@SQ\tSN:chr1\tLN:1000\n" +
	"a\t99\tchr1\t10\t60\t10M\t=\t100\t100\tACGTACGTAC\tIIIIIIIIII\tms:i:200\tMC:Z:10M\n" +
	"b\t99\tchr1\t10\t60\t10M\t=\t100\t100\tACGTACGTAC\tIIIIIIIIII\tms:i:200\tMC:Z:10M\tSA:Z:chr1,500,+,10M,60,0;\n" +
	"c\t99\tchr1\t30\t60\t10M\t=\t130\t100\tACGTACGTAC\tIIIIIIIIII\tms:i:200\tMC:Z:10M\n" +
	"a\t147\tchr1\t100\t60\t10M\t=\t10\t-100\tACGTACGTAC\tIIIIIIIIII\tms:i:200\tMC:Z:10M\n" +
	"b\t147\tchr1\t100\t60\t10M\t=\t10\t-100\tACGTACGTAC\tIIIIIIIIII\tms:i:200\tMC:Z:10M\n" +
	"c\t147\tchr1\t130\t60\t10M\t=\t30\t-100\tACGTACGTAC\tIIIIIIIIII\tms:i:200\tMC:Z:10M\n" +
	"b\t2049\tchr1\t500\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\tSA:Z:chr1,10,+,10M,60,0;\n"

// statsBody keeps only the comparable counter lines from a markdup stats
// report: it drops the "COMMAND:" line (which echoes the verbatim argv and so
// differs between the two binaries by design), blank lines, and any
// "samtools markdup: warning/error" diagnostics upstream prints to stderr but
// the Go port (by documented design) does not reproduce.
func statsBody(b []byte) string {
	var keep []string
	for _, l := range strings.Split(string(b), "\n") {
		if l == "" || strings.HasPrefix(l, "COMMAND:") || strings.HasPrefix(l, "samtools markdup:") {
			continue
		}
		keep = append(keep, l)
	}
	return strings.Join(keep, "\n")
}

// TestLiveMarkdupStats compares the markdup -s statistics block (sans the
// COMMAND line) and the marked records for the plain, -S and -d modes.
func TestLiveMarkdupStats(t *testing.T) {
	live := upstreamSamtools(t)
	ours := ourSamtoolsBinary(t)
	dir := t.TempDir()
	bam := writeSubfeatureBAM(t, live, dir, markdupStatsFixtureSAM)

	for _, flags := range [][]string{{"-s"}, {"-S", "-s"}} {
		args := append(append([]string{"markdup"}, flags...), bam, filepath.Join(dir, "out.bam"))
		_, upErr := runSamtoolsSubfeatures(t, live, args...)
		_, myErr := runSamtoolsSubfeatures(t, ours, args...)
		if statsBody(upErr) != statsBody(myErr) {
			t.Fatalf("markdup %v stats mismatch:\nupstream:\n%s\nours:\n%s", flags, statsBody(upErr), statsBody(myErr))
		}
	}

	// Marked records must agree for plain and -S modes.
	for _, flags := range [][]string{nil, {"-S"}} {
		args := append(append([]string{"markdup"}, flags...), bam, "-")
		up, _ := runSamtoolsSubfeatures(t, live, args...)
		mine, _ := runSamtoolsSubfeatures(t, ours, args...)
		if viewSAMRecords(t, live, up) != viewSAMRecords(t, live, mine) {
			t.Fatalf("markdup %v records mismatch:\nupstream:\n%s\nours:\n%s",
				flags, viewSAMRecords(t, live, up), viewSAMRecords(t, live, mine))
		}
	}
}

// markdupOpticalFixtureSAM uses Illumina-style colon-delimited read names so
// the -d optical detector can parse tile coordinates: a / b share a tile and
// are 5 units apart (optical at -d 100), e is on the same tile but far away
// (library duplicate).
const markdupOpticalFixtureSAM = "@HD\tVN:1.6\tSO:coordinate\n" +
	"@SQ\tSN:chr1\tLN:1000\n" +
	"M:1:FC:1:1:100:200\t99\tchr1\t10\t60\t10M\t=\t100\t100\tACGTACGTAC\tIIIIIIIIII\tms:i:200\tMC:Z:10M\n" +
	"M:1:FC:1:1:105:205\t99\tchr1\t10\t60\t10M\t=\t100\t100\tACGTACGTAC\tIIIIIIIIII\tms:i:200\tMC:Z:10M\n" +
	"M:1:FC:1:1:900:900\t99\tchr1\t10\t60\t10M\t=\t100\t100\tACGTACGTAC\tIIIIIIIIII\tms:i:200\tMC:Z:10M\n" +
	"M:1:FC:1:1:100:200\t147\tchr1\t100\t60\t10M\t=\t10\t-100\tACGTACGTAC\tIIIIIIIIII\tms:i:200\tMC:Z:10M\n" +
	"M:1:FC:1:1:105:205\t147\tchr1\t100\t60\t10M\t=\t10\t-100\tACGTACGTAC\tIIIIIIIIII\tms:i:200\tMC:Z:10M\n" +
	"M:1:FC:1:1:900:900\t147\tchr1\t100\t60\t10M\t=\t10\t-100\tACGTACGTAC\tIIIIIIIIII\tms:i:200\tMC:Z:10M\n"

// TestLiveMarkdupOptical asserts that -d marks optical duplicates with dt:Z:SQ
// and library duplicates with dt:Z:LB, byte-for-byte against upstream, and
// that the optical counters in the stats block agree.
func TestLiveMarkdupOptical(t *testing.T) {
	live := upstreamSamtools(t)
	ours := ourSamtoolsBinary(t)
	dir := t.TempDir()
	bam := writeSubfeatureBAM(t, live, dir, markdupOpticalFixtureSAM)

	up, _ := runSamtoolsSubfeatures(t, live, "markdup", "-d", "100", bam, "-")
	mine, _ := runSamtoolsSubfeatures(t, ours, "markdup", "-d", "100", bam, "-")
	if viewSAMRecords(t, live, up) != viewSAMRecords(t, live, mine) {
		t.Fatalf("markdup -d records mismatch:\nupstream:\n%s\nours:\n%s",
			viewSAMRecords(t, live, up), viewSAMRecords(t, live, mine))
	}

	args := []string{"markdup", "-d", "100", "-s", bam, filepath.Join(dir, "out.bam")}
	_, upErr := runSamtoolsSubfeatures(t, live, args...)
	_, myErr := runSamtoolsSubfeatures(t, ours, args...)
	if statsBody(upErr) != statsBody(myErr) {
		t.Fatalf("markdup -d stats mismatch:\nupstream:\n%s\nours:\n%s", statsBody(upErr), statsBody(myErr))
	}
}

// calmdSubfeatureFixtureSAM has two short reads with a single mismatch each,
// against a 32bp reference, for the -e / -u / -C knobs.
const calmdSubfeatureRef = ">chr1\nACGTACGTACGTACGTACGTACGTACGTACGT\n"
const calmdSubfeatureSAM = "@HD\tVN:1.6\tSO:coordinate\n" +
	"@SQ\tSN:chr1\tLN:32\n" +
	"r1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGAAC\tIIIIIIIIII\n" +
	"r2\t0\tchr1\t5\t60\t8M\t*\t0\t0\tACGTACGT\tIIIIIIII\n"

// TestLiveCalmdKnobs verifies the -e (= SEQ), -u (uncompressed BAM) and -C
// (mapQ cap) flags produce records identical to upstream (header @PG aside).
func TestLiveCalmdKnobs(t *testing.T) {
	live := upstreamSamtools(t)
	ours := ourSamtoolsBinary(t)
	dir := t.TempDir()

	refPath := filepath.Join(dir, "ref.fa")
	if err := os.WriteFile(refPath, []byte(calmdSubfeatureRef), 0o644); err != nil {
		t.Fatal(err)
	}
	if out, err := exec.Command(live, "faidx", refPath).CombinedOutput(); err != nil {
		t.Fatalf("faidx failed: %v\n%s", err, out)
	}
	bam := writeSubfeatureBAM(t, live, dir, calmdSubfeatureSAM)

	// Text-SAM knobs: -e and -C compared after dropping the @PG line.
	for _, flags := range [][]string{{"-e"}, {"-C", "50"}, {"-e", "-C", "50"}} {
		args := append(append([]string{"calmd"}, flags...), bam, refPath)
		up, _ := runSamtoolsSubfeatures(t, live, args...)
		mine, _ := runSamtoolsSubfeatures(t, ours, args...)
		if dropPGLines(up) != dropPGLines(mine) {
			t.Fatalf("calmd %v mismatch:\nupstream:\n%s\nours:\n%s", flags, dropPGLines(up), dropPGLines(mine))
		}
	}

	// -u uncompressed BAM: decode both and compare records.
	up, _ := runSamtoolsSubfeatures(t, live, "calmd", "-u", bam, refPath)
	mine, _ := runSamtoolsSubfeatures(t, ours, "calmd", "-u", bam, refPath)
	if viewSAMRecords(t, live, up) != viewSAMRecords(t, live, mine) {
		t.Fatalf("calmd -u records mismatch:\nupstream:\n%s\nours:\n%s",
			viewSAMRecords(t, live, up), viewSAMRecords(t, live, mine))
	}
}

// dropPGLines removes @PG header lines so the Go port's intentional omission
// of @PG provenance does not perturb a text comparison.
func dropPGLines(b []byte) string {
	var keep []string
	for _, l := range strings.Split(string(b), "\n") {
		if strings.HasPrefix(l, "@PG") {
			continue
		}
		keep = append(keep, l)
	}
	return strings.Join(keep, "\n")
}
