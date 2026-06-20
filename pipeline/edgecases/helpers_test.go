package edgecases

import (
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pipeline/internal/upstream"
)

// cacheDir returns a stable per-package build cache so OurBinary rebuilds each
// tool at most once across the whole test binary.
func cacheDir(t testing.TB) string {
	t.Helper()
	d := filepath.Join(os.TempDir(), "bioai-edgecases-cache")
	if err := os.MkdirAll(d, 0o755); err != nil {
		t.Fatalf("mkdir cache: %v", err)
	}
	return d
}

// ourBin builds (once) and returns the path to our binary for tool.
func ourBin(t testing.TB, tool string) string {
	t.Helper()
	bin, err := upstream.OurBinary(tool, cacheDir(t))
	if err != nil {
		t.Fatalf("building our %s: %v", tool, err)
	}
	return bin
}

// upBin resolves an upstream binary or SKIPs with an actionable message.
func upBin(t testing.TB, key string) string {
	t.Helper()
	bin, err := upstream.Binary(key)
	if err != nil {
		t.Skipf("upstream %s unavailable: %v", key, err)
	}
	return bin
}

// run executes bin with args and returns (stdout, stderr, err).
func run(t testing.TB, bin string, args ...string) (string, string, error) {
	t.Helper()
	cmd := exec.Command(bin, args...)
	var out, errb strings.Builder
	cmd.Stdout = &out
	cmd.Stderr = &errb
	err := cmd.Run()
	return out.String(), errb.String(), err
}

// mustRun runs bin and fails the test on a non-zero exit.
func mustRun(t testing.TB, bin string, args ...string) string {
	t.Helper()
	out, errOut, err := run(t, bin, args...)
	if err != nil {
		t.Fatalf("%s %s failed: %v\n%s", filepath.Base(bin), strings.Join(args, " "), err, errOut)
	}
	return out
}

// writeFile writes content to dir/name and returns the full path.
func writeFile(t testing.TB, dir, name, content string) string {
	t.Helper()
	p := filepath.Join(dir, name)
	if err := os.WriteFile(p, []byte(content), 0o644); err != nil {
		t.Fatalf("write %s: %v", p, err)
	}
	return p
}

// smallFixtureDir returns pipeline/.fixtures/small if it exists, else "".
func smallFixtureDir(t testing.TB) string {
	t.Helper()
	root, err := upstream.RepoRoot()
	if err != nil {
		return ""
	}
	d := filepath.Join(root, "pipeline", ".fixtures", "small")
	if st, err := os.Stat(d); err == nil && st.IsDir() {
		return d
	}
	return ""
}

// dropVCFHeaderNoise returns only the data lines of a VCF (no ## header), so two
// VCFs are compared on records, not on tool-version/##bcftools_* provenance
// lines that legitimately differ between producers.
func dropVCFHeaderNoise(vcf string) string {
	var out []string
	for _, ln := range strings.Split(vcf, "\n") {
		if strings.HasPrefix(ln, "##") {
			continue
		}
		if ln == "" {
			continue
		}
		out = append(out, ln)
	}
	return strings.Join(out, "\n")
}
