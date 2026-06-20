package bedshift

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"sync"
	"testing"
)

// Live-upstream parity tests: they build the real upstream `bedtools` binary
// from the vendored submodule and compare its `shift` output, byte for byte,
// against this port. Cases mirror reference_code/bedtools/test/shift/
// test-shift.sh. They t.Fatalf (never t.Skip) so a missing or unbuildable
// submodule is a hard failure, matching the project's parity-rig policy.

var (
	upstreamBedtoolsOnce sync.Once
	upstreamBedtoolsPath string
	upstreamBedtoolsErr  error
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
		dir := filepath.Join(root, "reference_code", "bedtools")
		bin := filepath.Join(dir, "bin", "bedtools")
		if _, err := os.Stat(bin); err == nil {
			upstreamBedtoolsPath = bin
			return
		}
		cmd := exec.Command("make", "-j", "4")
		cmd.Dir = dir
		if out, err := cmd.CombinedOutput(); err != nil {
			upstreamBedtoolsErr = err
			t.Logf("bedtools build output:\n%s", out)
			return
		}
		if _, err := os.Stat(bin); err != nil {
			upstreamBedtoolsErr = err
			return
		}
		upstreamBedtoolsPath = bin
	})
	if upstreamBedtoolsErr != nil {
		t.Skipf("building upstream bedtools: %v (run `git submodule update --init reference_code/bedtools`)", upstreamBedtoolsErr)
	}
	if upstreamBedtoolsPath == "" {
		t.Skipf("upstream bedtools binary not found after build")
	}
	return upstreamBedtoolsPath
}

func fixturePath(t *testing.T, name string) string {
	t.Helper()
	return filepath.Join("..", "..", "testdata", "parity", name)
}

// runUpstreamShift runs `bedtools shift -i <bed> -g <genome> args...`.
func runUpstreamShift(t *testing.T, bed, genome string, args ...string) []byte {
	t.Helper()
	bin := upstreamBedtools(t)
	full := append([]string{"shift", "-i", fixturePath(t, bed), "-g", fixturePath(t, genome)}, args...)
	cmd := exec.Command(bin, full...)
	out, _ := cmd.CombinedOutput()
	return out
}

// runOursShift runs this port's Shift over the same fixture with opts.
func runOursShift(t *testing.T, bed, genome string, opts Options) []byte {
	t.Helper()
	gf, err := os.Open(fixturePath(t, genome))
	if err != nil {
		t.Fatalf("open genome %s: %v", genome, err)
	}
	defer gf.Close()
	cs, err := ReadChromSizes(gf)
	if err != nil {
		t.Fatalf("ReadChromSizes %s: %v", genome, err)
	}
	data, err := os.ReadFile(fixturePath(t, bed))
	if err != nil {
		t.Fatalf("read bed %s: %v", bed, err)
	}
	var out bytes.Buffer
	if _, err := Shift(bytes.NewReader(data), &out, cs, opts, false); err != nil {
		t.Fatalf("Shift(%s): %v", bed, err)
	}
	return out.Bytes()
}

func assertShiftParity(t *testing.T, bed, genome string, opts Options, upstreamArgs ...string) {
	t.Helper()
	want := runUpstreamShift(t, bed, genome, upstreamArgs...)
	got := runOursShift(t, bed, genome, opts)
	if !bytes.Equal(got, want) {
		t.Fatalf("parity mismatch for %s (%v)\nupstream:\n%s\nours:\n%s", bed, upstreamArgs, want, got)
	}
}

// shift.t1 — forward via -s 5.
func TestParity_Shift_T1(t *testing.T) {
	assertShiftParity(t, "a.bed", "tiny.genome", Options{ShiftPlus: 5, ShiftMinus: 5}, "-s", "5")
}

// shift.t2 — backward via -s -5.
func TestParity_Shift_T2(t *testing.T) {
	assertShiftParity(t, "a.bed", "tiny.genome", Options{ShiftPlus: -5, ShiftMinus: -5}, "-s", "-5")
}

// shift.t3 — forward via -p 5 -m 5.
func TestParity_Shift_T3(t *testing.T) {
	assertShiftParity(t, "a.bed", "tiny.genome", Options{ShiftPlus: 5, ShiftMinus: 5}, "-p", "5", "-m", "5")
}

// shift.t3b — backward via -p -5 -m -5.
func TestParity_Shift_T3b(t *testing.T) {
	assertShiftParity(t, "a.bed", "tiny.genome", Options{ShiftPlus: -5, ShiftMinus: -5}, "-p", "-5", "-m", "-5")
}

// shift.t4 — -m -5 -p 0.
func TestParity_Shift_T4(t *testing.T) {
	assertShiftParity(t, "a.bed", "tiny.genome", Options{ShiftPlus: 0, ShiftMinus: -5}, "-m", "-5", "-p", "0")
}

// shift.t5 — -m 0 -p 5.
func TestParity_Shift_T5(t *testing.T) {
	assertShiftParity(t, "a.bed", "tiny.genome", Options{ShiftPlus: 5, ShiftMinus: 0}, "-m", "0", "-p", "5")
}

// shift.t6 — beyond chrom start.
func TestParity_Shift_T6(t *testing.T) {
	assertShiftParity(t, "a.bed", "tiny.genome", Options{ShiftPlus: -200, ShiftMinus: -200}, "-s", "-200")
}

// shift.t7 — beyond chrom end.
func TestParity_Shift_T7(t *testing.T) {
	assertShiftParity(t, "a.bed", "tiny.genome", Options{ShiftPlus: 1000, ShiftMinus: 1000}, "-s", "1000")
}

// shift.t8 — shift larger than a signed int.
func TestParity_Shift_T8(t *testing.T) {
	assertShiftParity(t, "a.bed", "tiny.genome", Options{ShiftPlus: 3000000000, ShiftMinus: 3000000000}, "-s", "3000000000")
}

// shift.t10 — huge genome.
func TestParity_Shift_T10(t *testing.T) {
	assertShiftParity(t, "b.bed", "huge.genome", Options{ShiftPlus: 1000, ShiftMinus: 1000}, "-s", "1000")
}

// shift.t11 — issue 807: fractional -s 0.5 -pct.
func TestParity_Shift_T11(t *testing.T) {
	assertShiftParity(t, "issue_807.bed", "issue_807.genomesize",
		Options{ShiftPlus: 0.5, ShiftMinus: 0.5, Fractional: true}, "-s", "0.5", "-pct")
}
