package bcftools

// TestMain for the live-oracle suite. Builds the local bcftools port
// binary into a per-suite temp dir exactly once so each individual
// TestLive* case can shell out to it cheaply.
//
// If the build fails we record the failure into ourBinPath = "" and
// every live-oracle test will t.Skip. The non-live tests (parity,
// stats_live, filter_live, mpileup_golden, etc.) continue to run.

import (
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"testing"
)

func TestMain(m *testing.M) {
	tmp, err := os.MkdirTemp("", "bcftools-live-oracle-*")
	if err != nil {
		fmt.Fprintln(os.Stderr, "live-oracle: failed to make tempdir:", err)
		os.Exit(m.Run())
	}
	bin := filepath.Join(tmp, "bcftools")
	// Build from the cmd dir using a path relative to repo root. The
	// test binary's cwd is the package dir, so we walk back to the
	// module root.
	pkg := "../../cmd/bcftools"
	cmd := exec.Command("go", "build", "-o", bin, pkg)
	cmd.Stderr = os.Stderr
	if err := cmd.Run(); err != nil {
		// Leave ourBinPath empty so requireLive will Skip.
		fmt.Fprintln(os.Stderr, "live-oracle: go build failed:", err)
	} else {
		ourBinPath = bin
	}
	code := m.Run()
	_ = os.RemoveAll(tmp)
	os.Exit(code)
}
