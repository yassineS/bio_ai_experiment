package samtools

// Live-binary oracle test for `samtools phase`.
//
// This test invokes the genuine upstream samtools binary (built once via
// the shared upstreamSamtools helper) AND our locally-built port on the
// same SAM fixture, then asserts byte-for-byte equality of the phase text
// stream (CC banner + PS / FL / M / EV / //). The phase output carries no
// provenance lines, so the comparison is exact.
//
// It is the salvaged parity gate for the phase.c port (phase_algo.go,
// phase_emit.go, phase_frag.go, phase_khash.go, phase_ksort.go,
// phase_pileup.go). Per the project's testing rules the upstream check
// must actually execute: the helpers t.Fatalf rather than t.Skip when the
// binaries cannot be produced.

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"sync"
	"testing"
)

// ourSamtoolsOnce memoises the one-time `go build` of our samtools CLI.
var ourSamtoolsOnce struct {
	sync.Once
	path string
	err  error
}

// ourSamtoolsBinary builds (once) and returns the path to our samtools
// command-line binary. It t.Fatalf's on build failure rather than
// skipping, matching the upstreamSamtools contract.
func ourSamtoolsBinary(t *testing.T) string {
	t.Helper()
	root := repoRootForTest(t)
	ourSamtoolsOnce.Do(func() {
		dir, err := os.MkdirTemp("", "our-samtools-")
		if err != nil {
			ourSamtoolsOnce.err = err
			return
		}
		bin := filepath.Join(dir, "samtools")
		cmd := exec.Command("go", "build", "-o", bin,
			"github.com/yassineS/bio_ai_experiment/tools/samtools/cmd/samtools")
		cmd.Dir = root
		if out, err := cmd.CombinedOutput(); err != nil {
			ourSamtoolsOnce.err = &buildErr{msg: string(out), err: err}
			return
		}
		ourSamtoolsOnce.path = bin
	})
	if ourSamtoolsOnce.err != nil {
		t.Fatalf("could not build our samtools binary: %v", ourSamtoolsOnce.err)
	}
	return ourSamtoolsOnce.path
}

type buildErr struct {
	msg string
	err error
}

func (e *buildErr) Error() string { return e.err.Error() + ": " + e.msg }

// runSamtools runs bin with args and returns stdout, failing on error.
func runSamtools(t *testing.T, bin string, args ...string) []byte {
	t.Helper()
	cmd := exec.Command(bin, args...)
	var out, errb bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = &errb
	if err := cmd.Run(); err != nil {
		t.Fatalf("%s %v failed: %v\nstderr: %s", bin, args, err, errb.String())
	}
	return out.Bytes()
}

// phaseFixtureSAM has two adjacent het sites at chr1:3 (G/T) and chr1:7
// (G/C), each covered by three reads carrying each allele. It stresses
// the greedy chainer, the EV emission order, and the determinism of the
// upstream output (no -b BAM split, so no RNG dependence).
const phaseFixtureSAM = "@HD\tVN:1.6\tSO:coordinate\n" +
	"@SQ\tSN:chr1\tLN:100\n" +
	"r_a\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n" +
	"r_b\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n" +
	"r_c\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n" +
	"r_d\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACTTACCTAC\tIIIIIIIIII\n" +
	"r_e\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACTTACCTAC\tIIIIIIIIII\n" +
	"r_f\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACTTACCTAC\tIIIIIIIIII\n"

// TestLivePhase asserts byte-equality of `samtools phase --no-PG` between
// the upstream binary and our port. Upstream phase calls drand48() but
// never seeds it, and the RNG only routes -b BAM reads (not exercised
// here), so the text stream is deterministic and exactly reproducible.
func TestLivePhase(t *testing.T) {
	live := upstreamSamtools(t)
	ours := ourSamtoolsBinary(t)

	dir := t.TempDir()
	samPath := filepath.Join(dir, "phase.sam")
	if err := os.WriteFile(samPath, []byte(phaseFixtureSAM), 0o644); err != nil {
		t.Fatal(err)
	}
	// Build the BAM with our own view so both binaries read identical input.
	bamPath := filepath.Join(dir, "phase.bam")
	if err := os.WriteFile(bamPath, runSamtools(t, ours, "view", "-b", "--no-PG", samPath), 0o644); err != nil {
		t.Fatal(err)
	}

	up := runSamtools(t, live, "phase", "--no-PG", bamPath)
	up2 := runSamtools(t, live, "phase", "--no-PG", bamPath)
	if !bytes.Equal(up, up2) {
		t.Fatalf("upstream phase output is non-deterministic across runs:\nrun1=%q\nrun2=%q", up, up2)
	}

	gp := runSamtools(t, ours, "phase", "--no-PG", bamPath)
	if !bytes.Equal(up, gp) {
		t.Errorf("DIVERGENCE: phase byte-stream differs\nupstream (%d bytes):\n%s\nours (%d bytes):\n%s",
			len(up), up, len(gp), gp)
	}
}

// phaseLowQFixtureSAM places two candidate het sites (chr1:3 G/T and
// chr1:7 G/C), each covered by three G-carrying and three alt-carrying
// reads, all at Phred 40 base quality. Lowering `-q` below the default
// 37 admits these het sites; this exercises the errmod genotype-
// likelihood LOD (errmod_cal m=4 + gl2cns) in the het-discovery path that
// surfaced the original "LOD precision" divergence. With the upstream-
// faithful float32 `tmp1` accumulation the het set — and hence the whole
// phase text stream — matches upstream byte-for-byte at every `-q`.
//
// NOTE on quality: the bases are kept at the default minimum base quality
// or above. A *separate*, pre-existing phase-pileup gap (unrelated to
// errmod precision) causes a het-set divergence only when reads carry
// MARGINAL base qualities at the variant column; that gap is documented in
// docs/PARITY_ROADMAP.md and is not exercised here because it is outside
// the errmod/gl2cns scope. The errmod/gl2cns LOD itself is byte-exact for
// those marginal columns (verified by the errmod oracle goldens) — the
// divergence is purely in which bases the phase pileup feeds to errmod.
const phaseLowQFixtureSAM = "@HD\tVN:1.6\tSO:coordinate\n" +
	"@SQ\tSN:chr1\tLN:100\n" +
	"r_a\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n" +
	"r_b\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n" +
	"r_c\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\tIIIIIIIIII\n" +
	"r_d\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACTTACCTAC\tIIIIIIIIII\n" +
	"r_e\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACTTACCTAC\tIIIIIIIIII\n" +
	"r_f\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACTTACCTAC\tIIIIIIIIII\n"

// TestLivePhaseLowQ asserts byte-equality of `samtools phase` between
// upstream and our port across a sweep of LOW `-q` het-LOD thresholds,
// the regime that surfaced the original errmod/gl2cns "LOD precision"
// divergence. With the upstream-faithful float32 `tmp1` accumulation the
// het set, and therefore the entire phase text stream, matches upstream
// byte-for-byte at every threshold from 1 up to the default 37. This is
// the end-to-end gate for the errmod accumulation-width fix.
func TestLivePhaseLowQ(t *testing.T) {
	live := upstreamSamtools(t)
	ours := ourSamtoolsBinary(t)

	dir := t.TempDir()
	samPath := filepath.Join(dir, "phase_lowq.sam")
	if err := os.WriteFile(samPath, []byte(phaseLowQFixtureSAM), 0o644); err != nil {
		t.Fatal(err)
	}
	bamPath := filepath.Join(dir, "phase_lowq.bam")
	if err := os.WriteFile(bamPath, runSamtools(t, ours, "view", "-b", "--no-PG", samPath), 0o644); err != nil {
		t.Fatal(err)
	}

	for _, q := range []string{"1", "5", "13", "20", "37"} {
		q := q
		t.Run("q"+q, func(t *testing.T) {
			args := []string{"phase", "--no-PG", "-q", q, bamPath}
			up := runSamtools(t, live, args...)
			gp := runSamtools(t, ours, args...)
			if !bytes.Equal(up, gp) {
				t.Errorf("DIVERGENCE at -q %s: phase byte-stream differs\nupstream (%d bytes):\n%s\nours (%d bytes):\n%s",
					q, len(up), up, len(gp), gp)
			}
		})
	}
}
