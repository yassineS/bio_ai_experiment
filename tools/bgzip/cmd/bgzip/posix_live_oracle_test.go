package main

// Live-binary oracle test for POSIX getopt-style short-flag bundling in
// bgzip, now that the CLI is routed through cliflag.Parse.
//
// BGZF is a deterministic container but the deflate bit-stream our Go
// compressor emits is not guaranteed byte-identical to upstream's zlib, so
// the cross-binary assertion is on the *decompressed* round-trip: our
// bundled `bgzip -cl6` output, fed to the upstream `bgzip -d`, must
// reproduce the original bytes exactly. The intra-binary assertion proves a
// bundled cluster (`-cl6`) parses and behaves identically to the spelled-out
// canonical form (`-c -l 6`).
//
// Per the project's testing rules the helpers t.Fatalf rather than t.Skip
// when the upstream binary cannot be built.

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"sync"
	"testing"
)

var (
	bgzipOurOnce sync.Once
	bgzipOurPath string
	bgzipOurErr  error
)

func buildOurBgzipBinary(t *testing.T) string {
	t.Helper()
	bgzipOurOnce.Do(func() {
		dir, err := os.MkdirTemp("", "our-bgzip-")
		if err != nil {
			bgzipOurErr = err
			return
		}
		bin := filepath.Join(dir, "bgzip")
		cmd := exec.Command("go", "build", "-o", bin, ".")
		if out, err := cmd.CombinedOutput(); err != nil {
			bgzipOurErr = bgzipBuildErr{err: err, out: out}
			return
		}
		bgzipOurPath = bin
	})
	if bgzipOurErr != nil {
		t.Fatalf("build our bgzip: %v", bgzipOurErr)
	}
	return bgzipOurPath
}

type bgzipBuildErr struct {
	err error
	out []byte
}

func (e bgzipBuildErr) Error() string { return e.err.Error() + ": " + string(e.out) }

func upstreamBgzipBinary(t *testing.T) string {
	t.Helper()
	root := bgzipRepoRoot(t)
	bin := filepath.Join(root, "reference_code", "htslib", "bgzip")
	if _, err := os.Stat(bin); err != nil {
		t.Fatalf("upstream bgzip not found at %s (build reference_code/htslib first): %v", bin, err)
	}
	return bin
}

func bgzipRepoRoot(t *testing.T) string {
	t.Helper()
	dir, err := os.Getwd()
	if err != nil {
		t.Fatal(err)
	}
	for {
		if _, err := os.Stat(filepath.Join(dir, "go.mod")); err == nil {
			return dir
		}
		parent := filepath.Dir(dir)
		if parent == dir {
			t.Fatal("could not locate repo root (go.mod)")
		}
		dir = parent
	}
}

// runBgzipStdin runs bin with args, feeding stdin, returning stdout bytes.
func runBgzipStdin(t *testing.T, bin string, stdin []byte, args ...string) []byte {
	t.Helper()
	cmd := exec.Command(bin, args...)
	cmd.Stdin = bytes.NewReader(stdin)
	var out, errb bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = &errb
	if err := cmd.Run(); err != nil {
		t.Fatalf("%s %v failed: %v\nstderr: %s", bin, args, err, errb.String())
	}
	return out.Bytes()
}

// TestLiveBgzipPosixBundledRoundTrip is the upstream-parity gate: our bundled
// `bgzip -cl6` (stdout + level 6) output must decompress, via the genuine
// upstream `bgzip -d`, back to the original payload.
func TestLiveBgzipPosixBundledRoundTrip(t *testing.T) {
	ours := buildOurBgzipBinary(t)
	up := upstreamBgzipBinary(t)

	payload := []byte("hello bgzip world\nsecond line\n" + makeFiller(5000))

	// Bundled `-cl6`: -c (stdout) and -l 6 (level) fused into one cluster —
	// exactly the form upstream getopt accepts and the form that must now
	// parse in our port.
	compressed := runBgzipStdin(t, ours, payload, "-cl6")
	decoded := runBgzipStdin(t, up, compressed, "-d")
	if !bytes.Equal(decoded, payload) {
		t.Fatalf("round-trip mismatch: decoded %d bytes, want %d", len(decoded), len(payload))
	}
}

// TestBgzipPosixBundlingEquivalentToCanonical proves, within our binary, that
// bundled / value-concatenated spellings behave identically to the canonical
// spelled-out forms (compared on the decompressed bytes, since the deflate
// stream is sensitive to internal buffering).
func TestBgzipPosixBundlingEquivalentToCanonical(t *testing.T) {
	ours := buildOurBgzipBinary(t)
	up := upstreamBgzipBinary(t)

	payload := []byte("payload for equivalence test\n" + makeFiller(3000))

	cases := []struct {
		name      string
		bundled   []string
		canonical []string
	}{
		{"cl6 == c l 6", []string{"-cl6"}, []string{"-c", "-l", "6"}},
		{"cl1 == c l 1", []string{"-cl1"}, []string{"-c", "-l", "1"}},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			bundled := runBgzipStdin(t, ours, payload, tc.bundled...)
			canonical := runBgzipStdin(t, ours, payload, tc.canonical...)
			// Decode both through upstream and compare payloads.
			db := runBgzipStdin(t, up, bundled, "-d")
			dc := runBgzipStdin(t, up, canonical, "-d")
			if !bytes.Equal(db, dc) || !bytes.Equal(db, payload) {
				t.Fatalf("bundled %v not equivalent to canonical %v", tc.bundled, tc.canonical)
			}
		})
	}
}

func makeFiller(n int) string {
	b := make([]byte, n)
	for i := range b {
		b[i] = byte('A' + (i % 26))
	}
	return string(b)
}
