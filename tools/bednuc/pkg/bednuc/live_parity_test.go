package bednuc

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"sync"
	"testing"
)

// Live-upstream parity tests for bednuc. They run the real upstream
// `bedtools` binary from the vendored submodule and compare its `nuc`
// output byte-for-byte against this port. They also exercise the bednuc CLI
// binary directly to prove that flag registration no longer panics (a
// duplicate `-seq` registration previously made `bednuc` abort on every
// run). They t.Fatalf (never t.Skip) when the upstream binary is absent,
// matching the project's parity-rig policy.

var (
	upstreamBedtoolsOnce sync.Once
	upstreamBedtoolsPath string
)

func repoRoot(t *testing.T) string {
	t.Helper()
	_, here, _, ok := runtime.Caller(0)
	if !ok {
		t.Fatalf("runtime.Caller failed")
	}
	dir := filepath.Dir(here)
	for {
		if _, err := os.Stat(filepath.Join(dir, "go.mod")); err == nil {
			return dir
		}
		parent := filepath.Dir(dir)
		if parent == dir {
			t.Fatalf("could not locate repo root (go.mod) above %s", here)
		}
		dir = parent
	}
}

func upstreamBedtools(t *testing.T) string {
	t.Helper()
	upstreamBedtoolsOnce.Do(func() {
		root := repoRoot(t)
		bin := filepath.Join(root, "reference_code", "bedtools", "bin", "bedtools")
		if _, err := os.Stat(bin); err == nil {
			upstreamBedtoolsPath = bin
		}
	})
	if upstreamBedtoolsPath == "" {
		t.Skipf("upstream bedtools binary not found at reference_code/bedtools/bin/bedtools " +
			"(run `git submodule update --init reference_code/bedtools` and build it)")
	}
	return upstreamBedtoolsPath
}

// writeTempFASTA writes a FASTA to a fresh temp dir and returns its path. A
// matching .fai is created on demand by both upstream and the port.
func writeTempFASTA(t *testing.T, body string) string {
	t.Helper()
	dir := t.TempDir()
	fa := filepath.Join(dir, "ref.fa")
	if err := os.WriteFile(fa, []byte(body), 0o644); err != nil {
		t.Fatalf("write fasta: %v", err)
	}
	return fa
}

func writeTemp(t *testing.T, name, body string) string {
	t.Helper()
	dir := t.TempDir()
	p := filepath.Join(dir, name)
	if err := os.WriteFile(p, []byte(body), 0o644); err != nil {
		t.Fatalf("write %s: %v", name, err)
	}
	return p
}

const liveFASTA = ">chr1\nACGTACGTACGTACGTACGTNNNNACGTACGT\n>chr2\nTTTTGGGGCCCCAAAATTTTGGGG\n"
const liveBED = "chr1\t0\t10\tf1\t0\t+\nchr2\t4\t12\tf2\t0\t-\n"

// TestLiveParity_Nuc_DefaultAndSeq proves bednuc.Run matches upstream
// `bedtools nuc` byte-for-byte in default and -seq modes, exercising the BED
// records that carry a score column.
func TestLiveParity_Nuc_DefaultAndSeq(t *testing.T) {
	bin := upstreamBedtools(t)
	fa := writeTempFASTA(t, liveFASTA)
	bed := writeTemp(t, "in.bed", liveBED)

	cases := []struct {
		name string
		args []string
		opts Options
	}{
		{"default", nil, Options{}},
		{"seq", []string{"-seq"}, Options{PrintSeq: true}},
		{"strand", []string{"-s"}, Options{ForceStrand: true}},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			args := append([]string{"nuc", "-fi", fa, "-bed", bed}, tc.args...)
			cmd := exec.Command(bin, args...)
			want, err := cmd.Output()
			if err != nil {
				t.Fatalf("upstream bedtools nuc %v: %v", tc.args, err)
			}
			bedData, err := os.ReadFile(bed)
			if err != nil {
				t.Fatalf("read bed: %v", err)
			}
			var got, warn bytes.Buffer
			if _, err := Run(bytes.NewReader(bedData), fa, &got, &warn, tc.opts); err != nil {
				t.Fatalf("Run: %v\nwarn: %s", err, warn.String())
			}
			if !bytes.Equal(got.Bytes(), want) {
				t.Fatalf("parity mismatch (%s)\nupstream:\n%s\nours:\n%s", tc.name, want, got.Bytes())
			}
		})
	}
}

// TestLiveCLI_NoPanic builds and runs the bednuc CLI to confirm it no longer
// panics at flag registration (the historical `flag redefined: seq` bug) and
// matches upstream `bedtools nuc` end-to-end through the real command path.
func TestLiveCLI_NoPanic(t *testing.T) {
	bin := upstreamBedtools(t)
	root := repoRoot(t)
	fa := writeTempFASTA(t, liveFASTA)
	bed := writeTemp(t, "in.bed", liveBED)

	// Run the bednuc CLI via `go run` so the actual main() flag wiring is
	// exercised (the panic was in main's flag registration, not in Run).
	cmd := exec.Command("go", "run", "./tools/bednuc/cmd/bednuc", "-fi", fa, "-bed", bed, "-seq")
	cmd.Dir = root
	got, err := cmd.CombinedOutput()
	if err != nil {
		t.Fatalf("bednuc CLI failed (panic regression?): %v\noutput:\n%s", err, got)
	}

	up := exec.Command(bin, "nuc", "-fi", fa, "-bed", bed, "-seq")
	want, err := up.Output()
	if err != nil {
		t.Fatalf("upstream bedtools nuc: %v", err)
	}
	if !bytes.Equal(got, want) {
		t.Fatalf("CLI parity mismatch\nupstream:\n%s\nbednuc CLI:\n%s", want, got)
	}
}
