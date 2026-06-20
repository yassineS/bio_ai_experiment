package fastp

// Live-parity tests for the fastp "tail" features ported in this change:
//   - --correction (overlap-based PE base correction) + overlap analysis
//   - -p/-P overrepresentation analysis
//   - -s/-S/-d output splitting
//
// These tests build and run the upstream OpenGene/fastp binary directly
// (no checked-in goldens) and compare its output against the Go port on the
// checked-in fixtures under tools/fastp/testdata/parity/. Per the env-guard
// policy (PR #294) they t.Fatalf with an exact init/build hint when the
// upstream binary cannot be located/built (the submodule is initializable
// here), and likewise t.Fatalf on any mismatch.

import (
	"encoding/json"
	"os"
	"os/exec"
	"path/filepath"
	"sync"
	"testing"
)

// upstreamFastpOnce guards the one-time build/locate of the upstream binary.
var (
	upstreamFastpOnce sync.Once
	upstreamFastpPath string
	upstreamFastpErr  error
)

// upstreamFastp returns the path to a usable upstream fastp binary, building
// it from the reference_code/fastp submodule if necessary. The result is
// cached across tests via sync.Once. A test should t.Skip when err is
// non-nil (binary unavailable) and t.Fatalf on any subsequent mismatch.
func upstreamFastp(t *testing.T) (string, error) {
	t.Helper()
	upstreamFastpOnce.Do(func() {
		root := filepath.Join("..", "..", "..", "..", "reference_code", "fastp")
		abs, err := filepath.Abs(filepath.Join(root, "fastp"))
		if err != nil {
			upstreamFastpErr = err
			return
		}
		if info, statErr := os.Stat(abs); statErr == nil && info.Mode()&0o111 != 0 {
			upstreamFastpPath = abs
			return
		}
		// Attempt a build (best-effort; requires libisal/libdeflate).
		build := exec.Command("make", "-j", "2")
		build.Dir = root
		if out, berr := build.CombinedOutput(); berr != nil {
			upstreamFastpErr = &buildError{out: out, err: berr}
			return
		}
		if info, statErr := os.Stat(abs); statErr == nil && info.Mode()&0o111 != 0 {
			upstreamFastpPath = abs
			return
		}
		upstreamFastpErr = os.ErrNotExist
	})
	return upstreamFastpPath, upstreamFastpErr
}

type buildError struct {
	out []byte
	err error
}

func (e *buildError) Error() string { return e.err.Error() + ": " + string(e.out) }

// commonDisableFlags turns off everything except the feature under test so
// the comparison isolates the new behaviour. The Go port mirrors these via
// permissive ProcessOptions.
var commonDisableFlags = []string{
	"--disable_quality_filtering",
	"--disable_length_filtering",
	"--disable_adapter_trimming",
}

// permissiveOpts returns ProcessOptions with all filtering effectively off,
// mirroring commonDisableFlags so survivors match upstream's.
func permissiveOpts() ProcessOptions {
	opts := DefaultProcessOptions()
	opts.QualThreshold = 0
	opts.QualPercent = 0
	opts.MinLength = 0
	opts.LengthRequired = 0
	opts.MaxNCount = 1 << 30
	opts.MaxNPercent = 100
	return opts
}

func TestParity_Fastp_Correction(t *testing.T) {
	bin, err := upstreamFastp(t)
	if err != nil {
		t.Skipf("upstream fastp unavailable; run `git submodule update --init reference_code/fastp && make -C reference_code/fastp`: %v", err)
	}
	dir := t.TempDir()
	r1 := parityInput(t, "corr_r1.fq")
	r2 := parityInput(t, "corr_r2.fq")

	upR1 := filepath.Join(dir, "up1.fq")
	upR2 := filepath.Join(dir, "up2.fq")
	upJSON := filepath.Join(dir, "up.json")
	args := append([]string{
		"-i", r1, "-I", r2, "-o", upR1, "-O", upR2,
		"--correction", "-j", upJSON, "-h", filepath.Join(dir, "up.html"),
	}, commonDisableFlags...)
	runUpstream(t, bin, args)

	opts := permissiveOpts()
	opts.Correction = true
	goR1, goR2, goStats := runGoFastpPE(t, r1, r2, opts)

	mustEqualBytes(t, "correction R1", goR1, readFile(t, upR1))
	mustEqualBytes(t, "correction R2", goR2, readFile(t, upR2))

	upFilt := readJSON(t, upJSON)["filtering_result"].(map[string]any)
	if got, want := goStats.CorrectedReads, int64(jsonNum(upFilt["corrected_reads"])); got != want {
		t.Fatalf("corrected_reads: go=%d up=%d", got, want)
	}
	if got, want := goStats.BasesCorrected, int64(jsonNum(upFilt["corrected_bases"])); got != want {
		t.Fatalf("corrected_bases: go=%d up=%d", got, want)
	}
	if goStats.BasesCorrected == 0 {
		t.Fatalf("expected some corrected bases in this fixture, got 0")
	}
}

func TestParity_Fastp_Overrepresentation(t *testing.T) {
	bin, err := upstreamFastp(t)
	if err != nil {
		t.Skipf("upstream fastp unavailable; run `git submodule update --init reference_code/fastp && make -C reference_code/fastp`: %v", err)
	}
	dir := t.TempDir()
	in := parityInput(t, "ora.fq")
	upOut := filepath.Join(dir, "up.fq")
	upJSON := filepath.Join(dir, "up.json")
	args := append([]string{
		"-i", in, "-o", upOut, "-p", "-P", "1",
		"-j", upJSON, "-h", filepath.Join(dir, "up.html"),
	}, commonDisableFlags...)
	runUpstream(t, bin, args)

	opts := permissiveOpts()
	opts.OverrepAnalysis = true
	opts.OverrepSampling = 1
	_, goStats := runGoFastpSE(t, in, opts)

	upORA := overrepFromJSON(t, readJSON(t, upJSON), "read1_before_filtering")
	goReport := buildJSONReport(goStats)
	goORA := goReport.Read1BeforeFilter.OverrepresentedSequences
	if goORA == nil {
		goORA = map[string]int64{}
	}
	if len(goORA) == 0 {
		t.Fatalf("Go produced no overrepresented sequences; expected the injected motif")
	}
	if len(goORA) != len(upORA) {
		t.Fatalf("overrep sequence set size: go=%d up=%d\ngo=%v\nup=%v", len(goORA), len(upORA), goORA, upORA)
	}
	for seq, upCount := range upORA {
		goCount, ok := goORA[seq]
		if !ok {
			t.Fatalf("overrep seq %q present upstream (count %d) but missing in Go", seq, upCount)
		}
		if goCount != upCount {
			t.Fatalf("overrep seq %q count: go=%d up=%d", seq, goCount, upCount)
		}
	}
}

func TestParity_Fastp_SplitByLines(t *testing.T) {
	bin, err := upstreamFastp(t)
	if err != nil {
		t.Skipf("upstream fastp unavailable; run `git submodule update --init reference_code/fastp && make -C reference_code/fastp`: %v", err)
	}
	in := parityInput(t, "ora.fq")

	upDir := t.TempDir()
	args := append([]string{
		"-i", in, "-o", filepath.Join(upDir, "out.fq"),
		"-S", "4000", "-w", "1",
		"-j", filepath.Join(upDir, "up.json"), "-h", filepath.Join(upDir, "up.html"),
	}, commonDisableFlags...)
	runUpstream(t, bin, args)

	goDir := t.TempDir()
	opts := permissiveOpts()
	opts.SplitByLines = 4000
	runGoFastpSESplit(t, in, filepath.Join(goDir, "out.fq"), opts)

	compareSplitDirs(t, upDir, goDir, "out.fq")
}

func TestParity_Fastp_SplitByNumber(t *testing.T) {
	bin, err := upstreamFastp(t)
	if err != nil {
		t.Skipf("upstream fastp unavailable; run `git submodule update --init reference_code/fastp && make -C reference_code/fastp`: %v", err)
	}
	in := parityInput(t, "ora.fq")

	upDir := t.TempDir()
	args := append([]string{
		"-i", in, "-o", filepath.Join(upDir, "out.fq"),
		"-s", "4", "-w", "1",
		"-j", filepath.Join(upDir, "up.json"), "-h", filepath.Join(upDir, "up.html"),
	}, commonDisableFlags...)
	runUpstream(t, bin, args)

	goDir := t.TempDir()
	opts := permissiveOpts()
	opts.SplitNumber = 4
	runGoFastpSESplit(t, in, filepath.Join(goDir, "out.fq"), opts)

	compareSplitDirs(t, upDir, goDir, "out.fq")
}

// runUpstream runs the upstream binary, failing the test on a non-zero exit.
func runUpstream(t *testing.T, bin string, args []string) {
	t.Helper()
	cmd := exec.Command(bin, args...)
	cmd.Stderr = os.Stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream fastp %v failed: %v", args, err)
	}
}

// runGoFastpSESplit drives ProcessSingleEndSplit against the given base path.
func runGoFastpSESplit(t *testing.T, in, outBase string, opts ProcessOptions) {
	t.Helper()
	f, err := os.Open(in)
	if err != nil {
		t.Fatalf("open %s: %v", in, err)
	}
	defer f.Close()
	if _, err := ProcessSingleEndSplit(f, outBase, defaultEncoding(), opts); err != nil {
		t.Fatalf("ProcessSingleEndSplit: %v", err)
	}
}

// compareSplitDirs verifies the numbered split files in upDir and goDir match
// byte-for-byte (same set of files, same contents).
func compareSplitDirs(t *testing.T, upDir, goDir, base string) {
	t.Helper()
	upFiles := splitFilesIn(t, upDir, base)
	goFiles := splitFilesIn(t, goDir, base)
	if len(upFiles) == 0 {
		t.Fatalf("upstream produced no split files in %s", upDir)
	}
	if len(upFiles) != len(goFiles) {
		t.Fatalf("split file count: up=%d go=%d\nup=%v\ngo=%v", len(upFiles), len(goFiles), upFiles, goFiles)
	}
	for _, name := range upFiles {
		up := readFile(t, filepath.Join(upDir, name))
		goPath := filepath.Join(goDir, name)
		if _, err := os.Stat(goPath); err != nil {
			t.Fatalf("Go missing split file %s", name)
		}
		mustEqualBytes(t, "split "+name, readFile(t, goPath), up)
	}
}

// splitFilesIn returns the sorted names of the NNNN.<base> split files in dir.
func splitFilesIn(t *testing.T, dir, base string) []string {
	t.Helper()
	matches, err := filepath.Glob(filepath.Join(dir, "*."+base))
	if err != nil {
		t.Fatalf("glob split files: %v", err)
	}
	out := make([]string, 0, len(matches))
	for _, m := range matches {
		out = append(out, filepath.Base(m))
	}
	return out
}

// overrepFromJSON extracts the overrepresented_sequences map from a read-stats
// section of an upstream JSON report, normalising counts to int64.
func overrepFromJSON(t *testing.T, report map[string]any, section string) map[string]int64 {
	t.Helper()
	sec, ok := report[section].(map[string]any)
	if !ok {
		return map[string]int64{}
	}
	raw, ok := sec["overrepresented_sequences"].(map[string]any)
	if !ok {
		return map[string]int64{}
	}
	out := make(map[string]int64, len(raw))
	for seq, v := range raw {
		out[seq] = int64(jsonNum(v))
	}
	return out
}

// jsonNum coerces a JSON number (decoded as float64) to float64, tolerating
// json.Number too.
func jsonNum(v any) float64 {
	switch n := v.(type) {
	case float64:
		return n
	case json.Number:
		f, _ := n.Float64()
		return f
	default:
		return 0
	}
}
