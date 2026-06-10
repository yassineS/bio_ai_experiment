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
