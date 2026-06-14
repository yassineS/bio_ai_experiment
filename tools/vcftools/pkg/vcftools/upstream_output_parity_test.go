package vcftools

// Live-binary parity tests for the per-individual missingness (.imiss)
// output and the --max-indv individual-thinning logic.
//
// These tests build (once per package test run) the upstream vcftools binary
// from the reference_code submodule and run BOTH it and the Go port over the
// same VCF fixture, comparing the relevant output file in-process. They do not
// rely on checked-in goldens: the expected output IS upstream's output.
//
// For --max-indv, upstream uses srand(time(NULL)) + random_shuffle (see
// reference_code/vcftools/src/cpp/variant_file_filters.cpp:105), so the
// kept-sample identity is non-deterministic and byte-parity is impossible.
// Per docs/PARITY_ROADMAP.md#rng--stochastic-output-policy we therefore assert
// the deterministic CONTRACT instead: the correct N individuals are kept, the
// selection is a subset of upstream's candidate set, and the Go port's choice
// is stable and documented (first N in header order).

import (
	"bufio"
	"errors"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"strconv"
	"strings"
	"sync"
	"testing"
)

// upstreamVcftools locates the upstream vcftools binary, building it from the
// reference_code/vcftools submodule on demand. The result is memoised so the
// (slow) autotools build runs at most once per package test process. The
// helper is intentionally self-contained in this file so it does not collide
// with sibling work touching other vcftools output writers.
var (
	upstreamVcftoolsOnce sync.Once
	upstreamVcftoolsPath string
	upstreamVcftoolsErr  error
)

// errUpstreamNotInitialised signals that the reference_code/vcftools submodule
// is not checked out, so upstream-comparison tests should skip rather than
// fail (distinct from a genuine build failure).
var errUpstreamNotInitialised = errors.New("reference_code/vcftools submodule not initialised")

func upstreamVcftools(t *testing.T) string {
	t.Helper()
	upstreamVcftoolsOnce.Do(func() {
		upstreamVcftoolsPath, upstreamVcftoolsErr = buildUpstreamVcftools()
	})
	if upstreamVcftoolsErr != nil {
		// The reference_code/vcftools submodule is the live-parity oracle.
		// Per the project's parity-rig policy (see PR #294), a missing
		// submodule is a hard failure with an init hint, not a silent skip;
		// a genuine build failure (submodule present but won't compile) also
		// fails hard.
		if errors.Is(upstreamVcftoolsErr, errUpstreamNotInitialised) {
			t.Fatalf("upstream vcftools submodule not initialised; run `git submodule update --init reference_code/vcftools` to enable live parity: %v", upstreamVcftoolsErr)
		}
		t.Fatalf("upstream vcftools unavailable: %v", upstreamVcftoolsErr)
	}
	return upstreamVcftoolsPath
}

// upstreamVcftoolsRepoRoot walks up from this test file to the repository root
// (the directory containing go.mod) and returns the vcftools submodule path.
func upstreamVcftoolsRepoRoot() (string, error) {
	_, thisFile, _, ok := runtime.Caller(0)
	if !ok {
		return "", os.ErrNotExist
	}
	dir := filepath.Dir(thisFile)
	for {
		if _, err := os.Stat(filepath.Join(dir, "go.mod")); err == nil {
			return dir, nil
		}
		parent := filepath.Dir(dir)
		if parent == dir {
			return "", os.ErrNotExist
		}
		dir = parent
	}
}

// buildUpstreamVcftools returns the path to a built upstream vcftools binary,
// running the autotools build if it is not already present.
func buildUpstreamVcftools() (string, error) {
	root, err := upstreamVcftoolsRepoRoot()
	if err != nil {
		return "", err
	}
	submodule := filepath.Join(root, "reference_code", "vcftools")
	binary := filepath.Join(submodule, "src", "cpp", "vcftools")
	if _, err := os.Stat(binary); err == nil {
		return binary, nil
	}

	// The submodule must already be initialised; the validation loop does
	// `git submodule update --init reference_code/vcftools` before running.
	if _, err := os.Stat(filepath.Join(submodule, "autogen.sh")); err != nil {
		return "", errUpstreamNotInitialised
	}
	steps := [][]string{
		{"./autogen.sh"},
		{"./configure"},
		{"make"},
	}
	for _, step := range steps {
		cmd := exec.Command(step[0], step[1:]...)
		cmd.Dir = submodule
		if out, err := cmd.CombinedOutput(); err != nil {
			return "", &buildError{step: strings.Join(step, " "), output: string(out), err: err}
		}
	}
	if _, err := os.Stat(binary); err != nil {
		return "", err
	}
	return binary, nil
}

type buildError struct {
	step   string
	output string
	err    error
}

func (e *buildError) Error() string {
	tail := e.output
	if len(tail) > 2000 {
		tail = tail[len(tail)-2000:]
	}
	return "build step '" + e.step + "' failed: " + e.err.Error() + "\n" + tail
}

// runUpstream runs the upstream binary with the given args plus --vcf/--out
// wiring and returns the output prefix used.
func runUpstream(t *testing.T, vcfPath string, extraArgs ...string) string {
	t.Helper()
	bin := upstreamVcftools(t)
	prefix := filepath.Join(t.TempDir(), "up")
	args := append([]string{"--vcf", vcfPath, "--out", prefix}, extraArgs...)
	cmd := exec.Command(bin, args...)
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("upstream run failed (%v):\n%s", err, out)
	}
	return prefix
}

// runGo runs the Go port over the same fixture file as upstream and returns
// the output prefix used.
func runGo(t *testing.T, vcfPath string, params *Params) string {
	t.Helper()
	prefix := filepath.Join(t.TempDir(), "go")
	params.OutPrefix = prefix
	in, err := os.Open(vcfPath)
	if err != nil {
		t.Fatalf("open vcf: %v", err)
	}
	defer in.Close()
	if err := Run(in, params); err != nil {
		t.Fatalf("Run: %v", err)
	}
	return prefix
}

// readLinesTrim reads a file into a slice of lines (no trailing empties).
func readLinesTrim(t *testing.T, path string) []string {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("open %s: %v", path, err)
	}
	defer f.Close()
	var lines []string
	sc := bufio.NewScanner(f)
	sc.Buffer(make([]byte, 0, 64*1024), 1<<20)
	for sc.Scan() {
		lines = append(lines, sc.Text())
	}
	if err := sc.Err(); err != nil {
		t.Fatalf("scan %s: %v", path, err)
	}
	return lines
}

// fixtureVCF returns the absolute path to a parity fixture VCF.
func fixtureVCF(t *testing.T, name string) string {
	t.Helper()
	return filepath.Join(vcftoolsFixtureDir(t), name)
}

// TestVcftools_MissingIndvUpstreamParity compares the .imiss output of the Go
// port against the live upstream binary across several genotype-filter
// configurations. The N_GENOTYPES_FILTERED column (and its effect on N_DATA /
// N_MISS / F_MISS) is exercised by the --minGQ and --minDP cases.
func TestVcftools_MissingIndvUpstreamParity(t *testing.T) {
	vcf := fixtureVCF(t, "sample.vcf")

	cases := []struct {
		name   string
		upArgs []string
		params *Params
	}{
		{
			name:   "no_geno_filter",
			upArgs: []string{"--missing-indv"},
			params: &Params{MissingIndv: true},
		},
		{
			name:   "min_gq_30",
			upArgs: []string{"--missing-indv", "--minGQ", "30"},
			params: &Params{MissingIndv: true, MinGQ: 30},
		},
		{
			name:   "min_dp_4",
			upArgs: []string{"--missing-indv", "--minDP", "4"},
			params: &Params{MissingIndv: true, MinDP: 4},
		},
		{
			name:   "min_gq_50",
			upArgs: []string{"--missing-indv", "--minGQ", "50"},
			params: &Params{MissingIndv: true, MinGQ: 50},
		},
	}

	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			upPrefix := runUpstream(t, vcf, tc.upArgs...)
			goPrefix := runGo(t, vcf, tc.params)

			want := readLinesTrim(t, upPrefix+".imiss")
			got := readLinesTrim(t, goPrefix+".imiss")
			if len(want) != len(got) {
				t.Fatalf(".imiss line count: want %d, got %d\nwant=%q\ngot=%q",
					len(want), len(got), want, got)
			}
			for i := range want {
				if want[i] != got[i] {
					t.Errorf(".imiss line %d mismatch:\nwant: %q\ngot:  %q", i, want[i], got[i])
				}
			}
		})
	}
}

// TestVcftools_MaxIndvUpstream validates --max-indv against the live upstream
// binary at the contract level mandated by the RNG policy: the same number of
// individuals must survive in both tools, and the Go port's deterministic
// selection (first N in header order) must be a subset of the input samples
// upstream draws from. Byte-parity is impossible because upstream shuffles
// with a time-seeded PRNG, so we do NOT assert identical sample identity.
func TestVcftools_MaxIndvUpstream(t *testing.T) {
	vcf := fixtureVCF(t, "sample.vcf")
	allSamples := recodeSampleNames(t, runUpstream(t, vcf, "--recode")+".recode.vcf")

	for _, n := range []int{1, 2, 3, 4} {
		n := n
		t.Run(maxIndvName(n), func(t *testing.T) {
			upPrefix := runUpstream(t, vcf, "--recode", "--max-indv", strconv.Itoa(n))
			goPrefix := runGo(t, vcf, &Params{MaxIndv: n, MaxIndvSet: true, Recode: true})

			upSamples := recodeSampleNames(t, upPrefix+".recode.vcf")
			goSamples := recodeSampleNames(t, goPrefix+".recode.vcf")

			wantCount := n
			if wantCount > len(allSamples) {
				wantCount = len(allSamples)
			}
			// COUNT contract: both tools keep min(N, |all|) individuals.
			if len(upSamples) != wantCount {
				t.Fatalf("upstream kept %d samples for --max-indv %d, want %d",
					len(upSamples), n, wantCount)
			}
			if len(goSamples) != wantCount {
				t.Errorf("Go kept %d samples for --max-indv %d, want %d (upstream kept %d)",
					len(goSamples), n, wantCount, len(upSamples))
			}
			// SUBSET contract: every Go-kept sample is a real input sample.
			set := make(map[string]bool, len(allSamples))
			for _, s := range allSamples {
				set[s] = true
			}
			for _, s := range goSamples {
				if !set[s] {
					t.Errorf("Go kept sample %q not present in input set %v", s, allSamples)
				}
			}
			// DETERMINISM contract: the Go selection is the first N in header
			// order (stable + documented, see PARITY_ROADMAP.md).
			wantGo := allSamples
			if n < len(allSamples) {
				wantGo = allSamples[:n]
			}
			if strings.Join(goSamples, ",") != strings.Join(wantGo, ",") {
				t.Errorf("Go --max-indv %d kept %v, want deterministic prefix %v", n, goSamples, wantGo)
			}
		})
	}
}

// recodeSampleNames extracts the sample column names from a recode VCF's
// #CHROM header line.
func recodeSampleNames(t *testing.T, path string) []string {
	t.Helper()
	for _, ln := range readLinesTrim(t, path) {
		if strings.HasPrefix(ln, "#CHROM") {
			cols := strings.Split(ln, "\t")
			if len(cols) <= 9 {
				return nil
			}
			return cols[9:]
		}
	}
	t.Fatalf("no #CHROM header in %s", path)
	return nil
}

func maxIndvName(n int) string { return "max_indv_" + strconv.Itoa(n) }

// TestVcftools_FreqUpstreamParity compares the --freq (.frq) and --counts
// (.frq.count) outputs of the Go port against the live upstream binary,
// byte-for-byte, on a purely-biallelic fixture.
//
// Upstream's output_frequency (variant_file_output.cpp:131) writes each
// allele frequency straight to a default-configured C++ ostream, i.e.
// `defaultfloat` with precision 6: six significant digits, trailing zeros
// stripped. So 6/12 prints "0.5" (not "0.500000"), 1/12 prints "0.0833333",
// and 11/12 prints "0.916667". The Go port previously emitted "%.6f"
// ("0.500000"); this test pins the corrected formatting.
//
// The fixture is deliberately all-biallelic so the Go port (which restricts
// --freq/--counts to biallelic loci — the multi-allelic gap tracked in
// docs/PARITY_ROADMAP.md#vcftools) emits exactly the same rows upstream does,
// making a byte-for-byte file comparison valid.
func TestVcftools_FreqUpstreamParity(t *testing.T) {
	vcf := fixtureVCF(t, "freq_fmt_fixture.vcf")

	cases := []struct {
		name   string
		upArg  string
		suffix string
		params *Params
	}{
		{"freq", "--freq", ".frq", &Params{Freq: true}},
		{"counts", "--counts", ".frq.count", &Params{Counts: true}},
		// --freq2 / --counts2 write the SAME .frq / .frq.count file as the
		// plain variants but strip the allele label (header {FREQ}/{COUNT},
		// bare values) via upstream's suppress_allele_output flag.
		{"freq2", "--freq2", ".frq", &Params{Freq2: true}},
		{"counts2", "--counts2", ".frq.count", &Params{Counts2: true}},
	}

	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			upPrefix := runUpstream(t, vcf, tc.upArg)
			goPrefix := runGo(t, vcf, tc.params)

			want := readLinesTrim(t, upPrefix+tc.suffix)
			got := readLinesTrim(t, goPrefix+tc.suffix)
			if len(want) != len(got) {
				t.Fatalf("%s line count: want %d, got %d\nwant=%q\ngot=%q",
					tc.suffix, len(want), len(got), want, got)
			}
			for i := range want {
				if want[i] != got[i] {
					t.Errorf("%s line %d mismatch:\nwant: %q\ngot:  %q",
						tc.suffix, i, want[i], got[i])
				}
			}
		})
	}
}
