// Package upstream locates the vendored upstream tool binaries under
// reference_code/ and our own freshly built tool binaries. Both the fixture
// generator and the parity runner depend on it.
//
// Upstream binaries are expected to already be built inside the (possibly
// submodule) reference_code/ tree, exactly as the existing per-tool live parity
// tests expect. In an isolated git worktree the submodules may be unpopulated;
// the documented workaround is to symlink the binaries from the main checkout
// (see pipeline/README.md). A missing binary produces a clear, actionable
// error rather than a silent skip.
package upstream

import (
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"sync"
)

// RepoRoot walks up from this source file to the module root (the directory
// containing go.mod).
func RepoRoot() (string, error) {
	_, file, _, _ := runtime.Caller(0)
	dir := filepath.Dir(file)
	for {
		if _, err := os.Stat(filepath.Join(dir, "go.mod")); err == nil {
			return dir, nil
		}
		parent := filepath.Dir(dir)
		if parent == dir {
			return "", fmt.Errorf("could not locate go.mod above %s", filepath.Dir(file))
		}
		dir = parent
	}
}

// relPaths maps an upstream binary key to its path(s) relative to repo root.
// The first existing candidate wins.
var relPaths = map[string][]string{
	"samtools": {"reference_code/samtools/samtools"},
	"bcftools": {"reference_code/bcftools/bcftools"},
	"bgzip":    {"reference_code/htslib/bgzip"},
	"tabix":    {"reference_code/htslib/tabix"},
	"htsfile":  {"reference_code/htslib/htsfile"},
	"bedtools": {"reference_code/bedtools/bin/bedtools"},
	"seqtk":    {"reference_code/seqtk/seqtk"},
	"sickle":   {"reference_code/sickle/sickle"},
	"skewer":   {"reference_code/skewer/skewer"},
	"fastp":    {"reference_code/fastp/fastp"},
	"vcftools": {"reference_code/vcftools/src/cpp/vcftools", "reference_code/vcftools/bin/vcftools"},
	// prinseq is a Perl script, not a compiled binary; the runner detects the
	// ".pl" suffix and invokes it through `perl` (see runner.timedRun).
	"prinseq": {"reference_code/prinseq/prinseq-lite.pl"},
}

// hint maps a key to the command that produces the missing binary.
var hint = map[string]string{
	"samtools": "git submodule update --init reference_code/samtools && (cd reference_code/samtools && make)",
	"bcftools": "git submodule update --init reference_code/bcftools && (cd reference_code/bcftools && make)",
	"bgzip":    "git submodule update --init reference_code/htslib && (cd reference_code/htslib && make)",
	"tabix":    "git submodule update --init reference_code/htslib && (cd reference_code/htslib && make)",
	"htsfile":  "git submodule update --init reference_code/htslib && (cd reference_code/htslib && make)",
	"bedtools": "git submodule update --init reference_code/bedtools && (cd reference_code/bedtools && make -j)",
	"seqtk":    "git submodule update --init reference_code/seqtk && (cd reference_code/seqtk && make)",
	"sickle":   "git submodule update --init reference_code/sickle && (cd reference_code/sickle && make)",
	"skewer":   "git submodule update --init reference_code/skewer && (cd reference_code/skewer && make)",
	"fastp":    "git submodule update --init reference_code/fastp && (cd reference_code/fastp && make)",
	"vcftools": "git submodule update --init reference_code/vcftools && (cd reference_code/vcftools && ./autogen.sh && ./configure && make)",
	"prinseq":  "git submodule update --init reference_code/prinseq (prinseq-lite.pl is a Perl script; ensure `perl` is on PATH)",
}

// MosdepthEnv is the environment variable a caller can set to point the runner
// at a prebuilt upstream mosdepth binary (mirroring the per-tool parity test's
// MOSDEPTH_BIN). mosdepth ships only as a linux/amd64 release asset, so the
// pipeline does not build it from source.
const MosdepthEnv = "MOSDEPTH_BIN"

// mosdepthCacheNames are the temp-dir cache file names the per-tool mosdepth
// parity tests download the upstream release binary into. The pipeline reuses
// that cache rather than building mosdepth (a Nim project) from source.
var mosdepthCacheNames = []string{"mosdepth_v0.3.14", "mosdepth"}

// Binary returns the absolute path to the vendored upstream binary for key
// (one of samtools, bcftools, bgzip, tabix, htsfile, bedtools, seqtk, sickle,
// skewer, fastp, vcftools, prinseq, mosdepth), or an actionable error.
//
// mosdepth is special: it ships only as a linux/amd64 GitHub release asset, so
// it is resolved from the MOSDEPTH_BIN environment variable or the temp-dir
// cache the per-tool parity tests populate, never from reference_code/.
func Binary(key string) (string, error) {
	if key == "mosdepth" {
		return mosdepthBinary()
	}
	root, err := RepoRoot()
	if err != nil {
		return "", err
	}
	cands, ok := relPaths[key]
	if !ok {
		return "", fmt.Errorf("unknown upstream binary key %q", key)
	}
	for _, rel := range cands {
		p := filepath.Join(root, rel)
		if st, err := os.Stat(p); err == nil && !st.IsDir() {
			return p, nil
		}
	}
	h := hint[key]
	return "", fmt.Errorf("upstream binary %q not found under reference_code/.\n"+
		"In an isolated worktree, symlink it from the main checkout, e.g.:\n"+
		"  ln -sf /path/to/main/%s %s\n"+
		"Or build it: %s", key, cands[0], cands[0], h)
}

// mosdepthBinary resolves the upstream mosdepth binary from MOSDEPTH_BIN or the
// temp-dir release cache. It returns an actionable error (never a silent skip)
// describing how to populate the cache; callers that want to skip on an
// unsupported platform check runtime.GOOS/GOARCH themselves.
func mosdepthBinary() (string, error) {
	if p := os.Getenv(MosdepthEnv); p != "" {
		if st, err := os.Stat(p); err == nil && !st.IsDir() {
			return p, nil
		}
		return "", fmt.Errorf("%s=%q does not point at an existing file", MosdepthEnv, p)
	}
	for _, name := range mosdepthCacheNames {
		p := filepath.Join(os.TempDir(), name)
		if st, err := os.Stat(p); err == nil && !st.IsDir() && st.Size() > 0 {
			return p, nil
		}
	}
	return "", fmt.Errorf("upstream mosdepth binary not found.\n"+
		"mosdepth ships only as a linux/amd64 GitHub release asset (it is a Nim\n"+
		"project, not built from source here). Provide it via one of:\n"+
		"  - set %s=/path/to/mosdepth, or\n"+
		"  - run the per-tool parity test once to populate the cache:\n"+
		"      go test ./tools/mosdepth/... -run Upstream\n"+
		"    (it downloads %s)", MosdepthEnv, mosdepthCacheNames[0])
}

// MosdepthSupported reports whether the current platform has a published
// upstream mosdepth release binary (linux/amd64 only).
func MosdepthSupported() bool {
	return runtime.GOOS == "linux" && runtime.GOARCH == "amd64"
}

var (
	ourBins   = map[string]string{}
	ourBinsMu sync.Mutex
)

// OurBinary builds tools/<tool>/cmd/<tool> once into a cache dir and returns
// the absolute path. Repeated calls for the same tool reuse the build.
func OurBinary(tool, cacheDir string) (string, error) {
	ourBinsMu.Lock()
	defer ourBinsMu.Unlock()
	if p, ok := ourBins[tool]; ok {
		return p, nil
	}
	if err := os.MkdirAll(cacheDir, 0o755); err != nil {
		return "", err
	}
	out := filepath.Join(cacheDir, tool)
	pkg := fmt.Sprintf("github.com/yassineS/bio_ai_experiment/tools/%s/cmd/%s", tool, tool)
	cmd := exec.Command("go", "build", "-o", out, pkg)
	if b, err := cmd.CombinedOutput(); err != nil {
		return "", fmt.Errorf("building our %s (%s): %v\n%s", tool, pkg, err, b)
	}
	ourBins[tool] = out
	return out, nil
}
