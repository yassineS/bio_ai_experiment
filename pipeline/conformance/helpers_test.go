package conformance

import (
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pipeline/internal/upstream"
)

// cacheDir returns a stable per-package build cache so that OurBinary only
// rebuilds each tool once across the whole test binary.
func cacheDir(t testing.TB) string {
	t.Helper()
	d := filepath.Join(os.TempDir(), "bioai-conformance-cache")
	if err := os.MkdirAll(d, 0o755); err != nil {
		t.Fatalf("mkdir cache: %v", err)
	}
	return d
}

// ourSamtools builds (once) and returns the path to our samtools binary.
func ourSamtools(t testing.TB) string {
	t.Helper()
	bin, err := upstream.OurBinary("samtools", cacheDir(t))
	if err != nil {
		t.Fatalf("building our samtools: %v", err)
	}
	return bin
}

// upSamtools resolves the upstream samtools binary or SKIPs.
func upSamtools(t testing.TB) string {
	t.Helper()
	bin, err := upstream.Binary("samtools")
	if err != nil {
		t.Skipf("upstream samtools unavailable: %v", err)
	}
	return bin
}

// htslibTest resolves reference_code/htslib/test or SKIPs with a pointer to the
// submodule-init command documented in docs/CONFORMANCE.md.
func htslibTest(t testing.TB) string {
	t.Helper()
	dir, ok := upstream.HtslibTestDir()
	if !ok {
		t.Skipf("htslib test corpus not initialised (%s missing); run:\n"+
			"  git submodule update --init reference_code/htslib\n"+
			"see docs/CONFORMANCE.md", dir)
	}
	return dir
}

// runCapture runs bin with args, returns stdout, and reports failure with the
// captured stderr when the command exits non-zero.
func runCapture(t testing.TB, bin string, args ...string) (stdout string, stderr string, err error) {
	t.Helper()
	cmd := exec.Command(bin, args...)
	var out, errb strings.Builder
	cmd.Stdout = &out
	cmd.Stderr = &errb
	err = cmd.Run()
	return out.String(), errb.String(), err
}
