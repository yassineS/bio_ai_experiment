package bcftools

import (
	"bytes"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"reflect"
	"sync"
	"testing"
	"time"
)

// statsUSRSection extracts the USR (user-tstv) section from a full stats
// report: the `# USR:` comment/header line plus every `USR:` data row. This
// is what the in-process upstream-parity comparison checks (no committed
// snapshot — both outputs are produced live in the same test run).
func statsUSRSection(stats []byte) []byte {
	var out bytes.Buffer
	for _, line := range bytes.Split(stats, []byte("\n")) {
		if bytes.HasPrefix(line, []byte("# USR:")) || bytes.HasPrefix(line, []byte("USR:")) {
			out.Write(line)
			out.WriteByte('\n')
		}
	}
	return out.Bytes()
}

func runUserTSTVStats(t *testing.T, in []byte, spec UserTSTVSpec) []byte {
	t.Helper()
	var out bytes.Buffer
	if _, err := Stats(bytes.NewReader(in), &out, StatsOptions{UserTSTV: []UserTSTVSpec{spec}}); err != nil {
		t.Fatalf("Stats: %v", err)
	}
	return out.Bytes()
}

// =====================================================================
// Unit tests (always run — no upstream binary or fixtures required)
// =====================================================================

// TestParseUserTSTV covers the TAG[idx][:min:max:n] grammar including
// defaults (min:max:n = 0:1:100), an explicit value index, partial binning
// specs, and error cases.
func TestParseUserTSTV(t *testing.T) {
	tests := []struct {
		name    string
		in      string
		want    UserTSTVSpec
		wantErr bool
	}{
		{"tag only defaults", "DP", UserTSTVSpec{Tag: "DP", Idx: 0, Min: 0, Max: 1, NBins: 100}, false},
		{"full spec", "BQB:0:1:10", UserTSTVSpec{Tag: "BQB", Idx: 0, Min: 0, Max: 1, NBins: 10}, false},
		{"int range", "MQ:0:100:5", UserTSTVSpec{Tag: "MQ", Idx: 0, Min: 0, Max: 100, NBins: 5}, false},
		{"indexed tag", "PV4[1]:0:1:4", UserTSTVSpec{Tag: "PV4", Idx: 1, Min: 0, Max: 1, NBins: 4}, false},
		{"indexed tag no binning", "PV4[2]", UserTSTVSpec{Tag: "PV4", Idx: 2, Min: 0, Max: 1, NBins: 100}, false},
		{"partial min only", "DP:5", UserTSTVSpec{Tag: "DP", Idx: 0, Min: 5, Max: 1, NBins: 100}, false},
		{"negative min", "DP:-2:2:4", UserTSTVSpec{Tag: "DP", Idx: 0, Min: -2, Max: 2, NBins: 4}, false},
		{"empty", "", UserTSTVSpec{}, true},
		{"bad min", "DP:x:1:4", UserTSTVSpec{}, true},
		{"bad n", "DP:0:1:x", UserTSTVSpec{}, true},
		{"zero bins", "DP:0:1:0", UserTSTVSpec{}, true},
		{"negative index", "PV4[-1]", UserTSTVSpec{}, true},
		{"too many fields", "DP:0:1:4:9", UserTSTVSpec{}, true},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := ParseUserTSTV(tt.in)
			if tt.wantErr {
				if err == nil {
					t.Fatalf("ParseUserTSTV(%q) = %+v, want error", tt.in, got)
				}
				return
			}
			if err != nil {
				t.Fatalf("ParseUserTSTV(%q) unexpected error: %v", tt.in, err)
			}
			if !reflect.DeepEqual(got, tt.want) {
				t.Fatalf("ParseUserTSTV(%q) = %+v, want %+v", tt.in, got, tt.want)
			}
		})
	}
}

// TestUserTSTVBin checks the binning rule against the upstream
// (val-min)/(max-min)*(nbins-1) formula, including the clamped endpoints.
func TestUserTSTVBin(t *testing.T) {
	spec := UserTSTVSpec{Min: 0, Max: 1, NBins: 10}
	tests := []struct {
		val  float64
		want int
	}{
		{-0.5, 0}, // <= min clamps to 0
		{0, 0},    // == min
		{0.05, 0}, // 0.05 * 9 = 0.45 -> 0
		{0.15, 1}, // 0.15 * 9 = 1.35 -> 1
		{0.95, 8}, // 0.95 * 9 = 8.55 -> 8
		{1.0, 9},  // >= max clamps to last bin
		{1.5, 9},  // > max clamps to last bin
	}
	for _, tt := range tests {
		if got := userTSTVBin(spec, tt.val); got != tt.want {
			t.Errorf("userTSTVBin(%v) = %d, want %d", tt.val, got, tt.want)
		}
	}
}

// TestUserTSTVBin_IntRange exercises a wider integer-style range so the
// scaling factor isn't a unit interval (mirrors the MQ:0:100:5 fixture
// case used by the live upstream parity test).
func TestUserTSTVBin_IntRange(t *testing.T) {
	spec := UserTSTVSpec{Min: 0, Max: 100, NBins: 5}
	tests := []struct {
		val  float64
		want int
	}{
		{5, 0},   // 5/100*4 = 0.2 -> 0
		{10, 0},  // 10/100*4 = 0.4 -> 0
		{20, 0},  // 20/100*4 = 0.8 -> 0
		{30, 1},  // 30/100*4 = 1.2 -> 1
		{50, 2},  // 50/100*4 = 2.0 -> 2
		{70, 2},  // 70/100*4 = 2.8 -> 2
		{95, 3},  // 95/100*4 = 3.8 -> 3
		{120, 4}, // >= max clamps to last bin
		{-1, 0},  // <= min clamps to 0
	}
	for _, tt := range tests {
		if got := userTSTVBin(spec, tt.val); got != tt.want {
			t.Errorf("userTSTVBin(%v) = %d, want %d", tt.val, got, tt.want)
		}
	}
}

// TestUserTSTV_TsTvClassification pins the transition/transversion split
// that drives the [5]ts vs [6]tv columns: A<->G and C<->T are transitions,
// everything else is a transversion (case-insensitive).
func TestUserTSTV_TsTvClassification(t *testing.T) {
	tests := []struct {
		ref, alt string
		want     string
	}{
		{"A", "G", "ts"},
		{"G", "A", "ts"},
		{"C", "T", "ts"},
		{"T", "C", "ts"},
		{"a", "g", "ts"}, // case-insensitive
		{"A", "C", "tv"},
		{"A", "T", "tv"},
		{"C", "G", "tv"},
		{"G", "T", "tv"},
		{"G", "C", "tv"},
		{"T", "A", "tv"},
	}
	for _, tt := range tests {
		if got := transitionType(tt.ref, tt.alt); got != tt.want {
			t.Errorf("transitionType(%q,%q) = %q, want %q", tt.ref, tt.alt, got, tt.want)
		}
	}
}

// TestStatsUserTSTV_MissingTag confirms that records without the requested
// INFO tag are silently skipped (no panic, empty USR data section).
func TestStatsUserTSTV_MissingTag(t *testing.T) {
	in := []byte("##fileformat=VCFv4.2\n" +
		"##contig=<ID=chr1,length=1000>\n" +
		"##INFO=<ID=AB,Number=1,Type=Float,Description=\"x\">\n" +
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
		"chr1\t100\t.\tA\tG\t50\t.\t.\n")
	got := runUserTSTVStats(t, in, UserTSTVSpec{Tag: "AB", Min: 0, Max: 1, NBins: 4})
	for _, line := range bytes.Split(got, []byte("\n")) {
		if bytes.HasPrefix(line, []byte("USR:")) {
			t.Fatalf("expected no USR data rows for missing tag, got:\n%s", got)
		}
	}
}

// TestStatsUserTSTV_IntFloatFormat is a self-contained format check: the
// Float-typed accumulator prints bin values in `%e` (scientific) form while
// the Integer-typed one prints `%.0f`. This holds without any upstream
// binary, so it always runs.
func TestStatsUserTSTV_IntFloatFormat(t *testing.T) {
	in := readStatsFixtureVCF(t)

	float := statsUSRSection(runUserTSTVStats(t, in, UserTSTVSpec{Tag: "BQB", Min: 0, Max: 1, NBins: 10}))
	if !bytes.Contains(float, []byte("USR:BQB/0\t0\t0.000000e+00\t")) {
		t.Errorf("Float USR section missing scientific-format bin label:\n%s", float)
	}

	integer := statsUSRSection(runUserTSTVStats(t, in, UserTSTVSpec{Tag: "MQ", Min: 0, Max: 100, NBins: 5}))
	if !bytes.Contains(integer, []byte("USR:MQ/0\t0\t0\t")) || !bytes.Contains(integer, []byte("USR:MQ/0\t0\t100\t")) {
		t.Errorf("Integer USR section missing plain-integer bin labels:\n%s", integer)
	}
}

// =====================================================================
// Live upstream-binary parity (no committed golden snapshots)
// =====================================================================

// TestStats_UserTSTVUpstreamParity runs BOTH the upstream C `bcftools stats
// -u TAG:min:max:n` and our Go implementation on the same fixture, then
// asserts the USR:TAG/idx section matches byte-for-byte IN-PROCESS. The
// upstream binary is located (or built from the vendored submodules) at
// test time; if it genuinely cannot be built the test fails (it never
// skips). A Float-typed tag (BQB) and an Integer-typed tag (MQ) are both
// exercised so the `%e` vs `%.0f` bin-label divergence is covered live.
func TestStats_UserTSTVUpstreamParity(t *testing.T) {
	bin := ensureUpstreamBcftools(t)
	fixture := statsFixturePath(t)
	in := readStatsFixtureVCF(t)

	cases := []struct {
		name string
		spec UserTSTVSpec
		arg  string // upstream -u argument
	}{
		{"Float_BQB", UserTSTVSpec{Tag: "BQB", Min: 0, Max: 1, NBins: 10}, "BQB:0:1:10"},
		{"Int_MQ", UserTSTVSpec{Tag: "MQ", Min: 0, Max: 100, NBins: 5}, "MQ:0:100:5"},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			cmd := exec.Command(bin, "stats", "-u", c.arg, fixture)
			cmd.Stderr = os.Stderr
			upOut, err := cmd.Output()
			if err != nil {
				t.Fatalf("upstream bcftools stats -u %s failed: %v", c.arg, err)
			}
			want := statsUSRSection(upOut)
			got := statsUSRSection(runUserTSTVStats(t, in, c.spec))
			if len(want) == 0 {
				t.Fatalf("upstream produced no USR section for -u %s", c.arg)
			}
			equalBytes(t, got, want, fmt.Sprintf("stats -u %s (live upstream parity)", c.arg))
		})
	}
}

// statsFixturePath returns the absolute path of the shared user-tstv input
// VCF fixture.
func statsFixturePath(t *testing.T) string {
	t.Helper()
	abs, err := filepath.Abs(filepath.Join("..", "..", "testdata", "parity", "stats", "user_tstv.vcf"))
	if err != nil {
		t.Fatalf("abs path: %v", err)
	}
	return abs
}

// readStatsFixtureVCF reads the shared user-tstv input fixture. Unlike the
// removed golden helper it never skips: the fixture is committed, so a read
// failure is a real test error.
func readStatsFixtureVCF(t *testing.T) []byte {
	t.Helper()
	data, err := os.ReadFile(statsFixturePath(t))
	if err != nil {
		t.Fatalf("read user-tstv fixture: %v", err)
	}
	return data
}

var (
	upstreamBcftoolsOnce sync.Once
	upstreamBcftoolsPath string
	upstreamBcftoolsErr  error
)

// ensureUpstreamBcftools returns the path to a working upstream bcftools
// binary, building it from the vendored submodules on first use. The build
// is memoised with sync.Once so a single `go test ./tools/bcftools/...`
// invocation pays for it at most once. On unrecoverable failure the calling
// test is FAILED (t.Fatalf), never skipped — the owner's rule is that the
// live upstream parity test must actually execute.
func ensureUpstreamBcftools(t *testing.T) string {
	t.Helper()
	upstreamBcftoolsOnce.Do(func() {
		upstreamBcftoolsPath, upstreamBcftoolsErr = buildUpstreamBcftools(t)
	})
	if upstreamBcftoolsErr != nil {
		t.Fatalf("upstream bcftools unavailable: %v", upstreamBcftoolsErr)
	}
	return upstreamBcftoolsPath
}

// repoRootDir resolves the absolute path of the repository root (four levels
// above the package directory tools/bcftools/pkg/bcftools).
func repoRootDir(t *testing.T) string {
	t.Helper()
	abs, err := filepath.Abs(filepath.Join("..", "..", "..", ".."))
	if err != nil {
		t.Fatalf("abs repo root: %v", err)
	}
	return abs
}

// referenceRoot resolves the absolute path of reference_code/<name>.
func referenceRoot(t *testing.T, name string) string {
	t.Helper()
	return filepath.Join(repoRootDir(t), "reference_code", name)
}

// buildUpstreamBcftools locates an already-built `bcftools` under
// reference_code/bcftools, and otherwise initialises the htslib + bcftools
// submodules (recursively, so htslib's nested htscodecs is present) and
// builds them. It returns the absolute path to the resulting binary.
func buildUpstreamBcftools(t *testing.T) (string, error) {
	t.Helper()
	bcftoolsDir := referenceRoot(t, "bcftools")
	htslibDir := referenceRoot(t, "htslib")
	binary := filepath.Join(bcftoolsDir, "bcftools")

	if isExecutable(binary) {
		return binary, nil
	}

	// Initialise the submodules recursively from the repo root. htslib's
	// ./configure aborts without its nested htscodecs submodule, so a
	// non-recursive init is not enough.
	if err := runBuildStep(t, repoRootDir(t), "git", "submodule", "update", "--init", "--recursive",
		"reference_code/htslib", "reference_code/bcftools"); err != nil {
		return "", fmt.Errorf("submodule init: %w", err)
	}

	// Build htslib first (autoreconf + configure + make). bcftools links
	// against this in-tree htslib build.
	if err := runBuildStep(t, htslibDir, "autoreconf", "-i"); err != nil {
		return "", fmt.Errorf("htslib autoreconf: %w", err)
	}
	if err := runBuildStep(t, htslibDir, "./configure"); err != nil {
		return "", fmt.Errorf("htslib configure: %w", err)
	}
	if err := runBuildStep(t, htslibDir, "make", "-j"); err != nil {
		return "", fmt.Errorf("htslib make: %w", err)
	}

	// Build bcftools with plain make. Do NOT run bcftools' ./configure — it
	// can clobber htslib's config.mk and break the in-tree link.
	if err := runBuildStep(t, bcftoolsDir, "make", "-j"); err != nil {
		return "", fmt.Errorf("bcftools make: %w", err)
	}

	if !isExecutable(binary) {
		return "", fmt.Errorf("bcftools build finished but %s is missing/not executable", binary)
	}
	return binary, nil
}

// isExecutable reports whether path exists and has an execute bit set.
func isExecutable(path string) bool {
	info, err := os.Stat(path)
	return err == nil && !info.IsDir() && info.Mode()&0o111 != 0
}

// runBuildStep runs a build/setup command in dir (the calling process's cwd
// when dir is empty). Network-touching steps (git submodule update) are
// retried with exponential backoff (2/4/8/16s) so transient fetch failures
// don't fail the build outright.
func runBuildStep(t *testing.T, dir, name string, args ...string) error {
	t.Helper()
	network := name == "git"
	backoff := []time.Duration{2 * time.Second, 4 * time.Second, 8 * time.Second, 16 * time.Second}

	var lastErr error
	attempts := 1
	if network {
		attempts = len(backoff) + 1
	}
	for attempt := 0; attempt < attempts; attempt++ {
		if attempt > 0 {
			t.Logf("retrying %s %v (attempt %d) after %s", name, args, attempt+1, backoff[attempt-1])
			time.Sleep(backoff[attempt-1])
		}
		cmd := exec.Command(name, args...)
		if dir != "" {
			cmd.Dir = dir
		}
		out, err := cmd.CombinedOutput()
		if err == nil {
			return nil
		}
		lastErr = fmt.Errorf("%s %v: %v\n%s", name, args, err, out)
		if !network {
			break
		}
	}
	return lastErr
}
