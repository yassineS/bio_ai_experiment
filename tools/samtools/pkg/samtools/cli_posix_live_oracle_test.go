package samtools

// Live-binary oracle tests for the POSIX getopt short-flag bundling rolled
// out across the remaining samtools subcommands (flagstat, sort, idxstats,
// mpileup, depth, ...), now that each subcommand's argument parse is routed
// through cliflag.Parse.
//
// Each test runs the genuine upstream `samtools` and the Go port on the same
// fixture, using a representative BUNDLED command line, and asserts:
//   1. the bundled form parses and the two binaries' outputs agree (decoded
//      for binary BAM output, compared verbatim for text output); and
//   2. within our own binary, the bundled spelling is equivalent to the
//      canonical spelled-out form.
//
// Per the project's testing rules the helpers t.Fatalf rather than t.Skip
// when the binaries cannot be built. The upstreamSamtools / ourSamtoolsBinary
// / runSamtools / decodeBAMRecords helpers are shared with the other live
// tests in this package.

import (
	"bytes"
	"os"
	"path/filepath"
	"testing"
)

// cliPosixFixtureSAM is a small coordinate-sorted SAM with a spread of MAPQs
// and a het column, suitable for flagstat / depth / sort / mpileup.
const cliPosixFixtureSAM = "@HD\tVN:1.6\tSO:coordinate\n" +
	"@SQ\tSN:chr1\tLN:1000\n" +
	"@SQ\tSN:chr2\tLN:500\n" +
	"r1\t0\tchr1\t100\t60\t5M\t*\t0\t0\tACGTA\tIIIII\n" +
	"r2\t0\tchr1\t100\t40\t5M\t*\t0\t0\tACGTA\tIIIII\n" +
	"r3\t16\tchr1\t102\t30\t5M\t*\t0\t0\tCGTAC\tIIIII\n" +
	"r4\t0\tchr2\t10\t25\t4M\t*\t0\t0\tACGT\tIIII\n" +
	"r5\t4\t*\t0\t0\t*\t*\t0\t0\tACGTA\tIIIII\n"

// writeFixtureSAM writes the shared fixture SAM into dir and returns its path.
func writeFixtureSAM(t *testing.T, dir string) string {
	t.Helper()
	p := filepath.Join(dir, "in.sam")
	if err := os.WriteFile(p, []byte(cliPosixFixtureSAM), 0o644); err != nil {
		t.Fatal(err)
	}
	return p
}

// TestLiveFlagstatPosixBundled exercises `flagstat -@1` (a bundled
// value-concatenated thread count) and asserts the text report matches
// upstream and that the bundled thread spelling is a no-op equal to the bare
// invocation in our binary.
func TestLiveFlagstatPosixBundled(t *testing.T) {
	live := upstreamSamtools(t)
	ours := ourSamtoolsBinary(t)
	dir := t.TempDir()
	sam := writeFixtureSAM(t, dir)

	// `-@1`: threads value concatenated onto the short flag — must parse.
	upstream := runSamtools(t, live, "flagstat", "-@1", sam)
	mine := runSamtools(t, ours, "flagstat", "-@1", sam)
	if !bytes.Equal(upstream, mine) {
		t.Fatalf("flagstat -@1 mismatch:\nupstream:\n%s\nours:\n%s", upstream, mine)
	}

	// Bundled thread count is a no-op vs the bare invocation in our binary.
	bare := runSamtools(t, ours, "flagstat", sam)
	if !bytes.Equal(bare, mine) {
		t.Fatalf("flagstat -@1 differs from bare:\n-@1:\n%s\nbare:\n%s", mine, bare)
	}
}

// TestLiveSortPosixBundled exercises `sort -nu` (name-sort, uncompressed BAM)
// as a single bundled cluster: -n -u. It decodes both binaries' BAM output
// and compares the record order/identity.
func TestLiveSortPosixBundled(t *testing.T) {
	live := upstreamSamtools(t)
	ours := ourSamtoolsBinary(t)
	dir := t.TempDir()
	sam := writeFixtureSAM(t, dir)

	// `-nu`: -n (name sort) + -u (uncompressed BAM) bundled into one cluster.
	upstream := decodeBAMRecords(t, runSamtools(t, live, "sort", "-nu", sam))
	mine := decodeBAMRecords(t, runSamtools(t, ours, "sort", "-nu", sam))
	if len(upstream) != len(mine) {
		t.Fatalf("sort -nu record count: upstream=%d ours=%d", len(upstream), len(mine))
	}
	for i := range upstream {
		if upstream[i] != mine[i] {
			t.Errorf("sort -nu record %d differs:\nupstream=%+v\nours=%+v", i, upstream[i], mine[i])
		}
	}

	// Within our binary, the bundled `-nu` equals the spelled-out `-n -u`.
	bundled := decodeBAMRecords(t, runSamtools(t, ours, "sort", "-nu", sam))
	canonical := decodeBAMRecords(t, runSamtools(t, ours, "sort", "-n", "-u", sam))
	if len(bundled) != len(canonical) {
		t.Fatalf("sort -nu vs -n -u count: %d vs %d", len(bundled), len(canonical))
	}
	for i := range bundled {
		if bundled[i] != canonical[i] {
			t.Errorf("sort -nu != -n -u at record %d:\n-nu=%+v\n-n -u=%+v", i, bundled[i], canonical[i])
		}
	}
}

// TestLiveMpileupPosixFusedAA exercises the headline fused `-aa` form ("all
// positions, all chromosomes"), which cliflag.Parse expands to the repeatable
// `-a -a` via the count flag. It asserts our text pileup matches upstream and
// that `-aa` is equivalent to `-a -a` within our binary.
func TestLiveMpileupPosixFusedAA(t *testing.T) {
	live := upstreamSamtools(t)
	ours := ourSamtoolsBinary(t)
	dir := t.TempDir()
	sam := writeFixtureSAM(t, dir)

	upstream := runSamtools(t, live, "mpileup", "-aa", sam)
	mine := runSamtools(t, ours, "mpileup", "-aa", sam)
	if !bytes.Equal(upstream, mine) {
		t.Fatalf("mpileup -aa mismatch:\nupstream:\n%s\nours:\n%s", upstream, mine)
	}

	// `-aa` == `-a -a` within our binary (the count flag resolves both to
	// all-positions-all-chroms).
	doubled := runSamtools(t, ours, "mpileup", "-a", "-a", sam)
	if !bytes.Equal(mine, doubled) {
		t.Fatalf("mpileup -aa differs from -a -a:\n-aa:\n%s\n-a -a:\n%s", mine, doubled)
	}
}

// TestLiveDepthPosixBundled exercises `depth -aa` (all positions, all
// chromosomes — repeatable -a), asserting the text output matches upstream
// and the bundled form equals `-a -a` in our binary.
func TestLiveDepthPosixBundled(t *testing.T) {
	live := upstreamSamtools(t)
	ours := ourSamtoolsBinary(t)
	dir := t.TempDir()
	sam := writeFixtureSAM(t, dir)

	upstream := runSamtools(t, live, "depth", "-aa", sam)
	mine := runSamtools(t, ours, "depth", "-aa", sam)
	if !bytes.Equal(upstream, mine) {
		t.Fatalf("depth -aa mismatch:\nupstream:\n%s\nours:\n%s", upstream, mine)
	}

	doubled := runSamtools(t, ours, "depth", "-a", "-a", sam)
	if !bytes.Equal(mine, doubled) {
		t.Fatalf("depth -aa differs from -a -a:\n-aa:\n%s\n-a -a:\n%s", mine, doubled)
	}
}
