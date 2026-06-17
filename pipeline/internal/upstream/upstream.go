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
	"bedtools": {"reference_code/bedtools/bin/bedtools"},
}

// hint maps a key to the command that produces the missing binary.
var hint = map[string]string{
	"samtools": "git submodule update --init reference_code/samtools && (cd reference_code/samtools && make)",
	"bcftools": "git submodule update --init reference_code/bcftools && (cd reference_code/bcftools && make)",
	"bgzip":    "git submodule update --init reference_code/htslib && (cd reference_code/htslib && make)",
	"tabix":    "git submodule update --init reference_code/htslib && (cd reference_code/htslib && make)",
	"bedtools": "git submodule update --init reference_code/bedtools && (cd reference_code/bedtools && make -j)",
}

// Binary returns the absolute path to the vendored upstream binary for key
// (one of samtools, bcftools, bgzip, tabix, bedtools), or an actionable error.
func Binary(key string) (string, error) {
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
