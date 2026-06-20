package bgzf

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"sync"
	"testing"
)

// upstreamBgzipOnce memoises the build of htslib's bgzip so the parity test
// pays the configure/make cost at most once per `go test` process.
var (
	upstreamBgzipOnce sync.Once
	upstreamBgzipPath string
	upstreamBgzipErr  error
)

// upstreamHtslibBgzip returns the path to a built htslib bgzip binary, building
// it from the reference_code/htslib submodule on first use. It is self-contained
// to this test: it inits the submodule (recursively, for the nested htscodecs
// dependency) when missing and runs autoreconf/configure/make for the bgzip
// target only. The result is memoised via upstreamBgzipOnce.
func upstreamHtslibBgzip(t *testing.T) string {
	t.Helper()
	upstreamBgzipOnce.Do(func() {
		root := repoRoot(t)
		htsDir := filepath.Join(root, "reference_code", "htslib")
		bin := filepath.Join(htsDir, "bgzip")

		if _, err := os.Stat(bin); err == nil {
			upstreamBgzipPath = bin
			return
		}

		// Init the submodule (recursive: htslib nests htscodecs).
		if _, err := os.Stat(filepath.Join(htsDir, "configure.ac")); err != nil {
			run(t, root, "git", "submodule", "update", "--init", "--recursive", "reference_code/htslib")
		}

		run(t, htsDir, "autoreconf", "-i")
		run(t, htsDir, "./configure")
		run(t, htsDir, "make", "-j4", "bgzip")

		if _, err := os.Stat(bin); err != nil {
			upstreamBgzipErr = err
			return
		}
		upstreamBgzipPath = bin
	})
	if upstreamBgzipErr != nil {
		t.Skipf("building upstream bgzip: %v", upstreamBgzipErr)
	}
	return upstreamBgzipPath
}

// repoRoot returns the module root by walking up from this test file's
// directory until it finds go.mod.
func repoRoot(t *testing.T) string {
	t.Helper()
	_, file, _, ok := runtime.Caller(0)
	if !ok {
		t.Fatalf("cannot determine caller path")
	}
	dir := filepath.Dir(file)
	for {
		if _, err := os.Stat(filepath.Join(dir, "go.mod")); err == nil {
			return dir
		}
		parent := filepath.Dir(dir)
		if parent == dir {
			t.Fatalf("go.mod not found walking up from %s", filepath.Dir(file))
		}
		dir = parent
	}
}

// run executes a command in dir and fails the test on error, surfacing output.
func run(t *testing.T, dir string, name string, args ...string) {
	t.Helper()
	cmd := exec.Command(name, args...)
	cmd.Dir = dir
	out, err := cmd.CombinedOutput()
	if err != nil {
		t.Fatalf("%s %v failed: %v\n%s", name, args, err, out)
	}
}

// upstreamDecompress runs `bgzip -d -c` (decompress to stdout) on data.
func upstreamDecompress(t *testing.T, bin string, data []byte) []byte {
	t.Helper()
	cmd := exec.Command(bin, "-d", "-c")
	cmd.Stdin = bytes.NewReader(data)
	var out, errBuf bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = &errBuf
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream bgzip -d failed: %v\n%s", err, errBuf.String())
	}
	return out.Bytes()
}

// upstreamCompress runs `bgzip -c -@ threads -l level` on data.
func upstreamCompress(t *testing.T, bin string, data []byte, threads, level int) []byte {
	t.Helper()
	cmd := exec.Command(bin, "-c", "-@", itoa(threads), "-l", itoa(level))
	cmd.Stdin = bytes.NewReader(data)
	var out, errBuf bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = &errBuf
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream bgzip compress failed: %v\n%s", err, errBuf.String())
	}
	return out.Bytes()
}

func itoa(n int) string {
	if n == 0 {
		return "0"
	}
	neg := n < 0
	if neg {
		n = -n
	}
	var b [20]byte
	i := len(b)
	for n > 0 {
		i--
		b[i] = byte('0' + n%10)
		n /= 10
	}
	if neg {
		i--
		b[i] = '-'
	}
	return string(b[i:])
}

// ourCompress compresses data with our MultiWriter at the given thread count.
func ourCompress(t *testing.T, data []byte, threads, level int) []byte {
	t.Helper()
	var buf bytes.Buffer
	w, err := NewMultiWriter(&buf, level, threads)
	if err != nil {
		t.Fatalf("NewMultiWriter: %v", err)
	}
	if _, err := w.Write(data); err != nil {
		t.Fatalf("Write: %v", err)
	}
	if err := w.Close(); err != nil {
		t.Fatalf("Close: %v", err)
	}
	return buf.Bytes()
}

// ourDecompress decompresses a BGZF stream with our Reader.
func ourDecompress(t *testing.T, data []byte) []byte {
	t.Helper()
	r, err := NewReader(bytes.NewReader(data))
	if err != nil {
		t.Fatalf("NewReader: %v", err)
	}
	var out bytes.Buffer
	if _, err := out.ReadFrom(r); err != nil {
		t.Fatalf("read: %v", err)
	}
	return out.Bytes()
}

// TestBgzip_ThreadsUpstreamParity is a live parity test against htslib's bgzip.
//
// Multi-threaded BGZF output is NOT guaranteed byte-identical to upstream
// (block boundaries and per-block deflate output differ between implementations
// and thread counts), so we validate the PLAINTEXT and structural invariants
// instead:
//
//	(a) our -@N output decompresses to the original plaintext via BOTH upstream
//	    `bgzip -d` and our own Reader;
//	(b) upstream `bgzip -@N` output decompresses to the same plaintext via our
//	    Reader (and trivially via upstream);
//	(c) our single-thread and multi-thread outputs decompress-equal;
//	(d) the stream is structurally valid (ends in the canonical EOF block, all
//	    blocks Scan cleanly).
func TestBgzip_ThreadsUpstreamParity(t *testing.T) {
	bin := upstreamHtslibBgzip(t)

	big := make([]byte, MaxBlockSize*6+4321)
	fillPseudoRandom(big)

	inputs := []struct {
		name string
		data []byte
	}{
		{"empty", nil},
		{"small", []byte("a single short line of text\n")},
		{"one-block", bytes.Repeat([]byte("Q"), MaxBlockSize)},
		{"multi-block-repetitive", bytes.Repeat([]byte("ACGTACGTNN\n"), 80000)},
		{"multi-block-random", big},
	}

	for _, in := range inputs {
		in := in
		t.Run(in.name, func(t *testing.T) {
			for _, threads := range []int{1, 2, 4} {
				level := DefaultCompression

				ours := ourCompress(t, in.data, threads, level)

				// (d) structural validity.
				if !bytes.HasSuffix(ours, EOFBlock) {
					t.Fatalf("threads=%d: missing EOF block", threads)
				}
				if _, err := Scan(bytes.NewReader(ours)); err != nil {
					t.Fatalf("threads=%d: Scan: %v", threads, err)
				}

				// (a) our output -> upstream decompress == original.
				if got := upstreamDecompress(t, bin, ours); !bytes.Equal(got, in.data) {
					t.Fatalf("threads=%d: upstream decode of our output != original", threads)
				}
				// (a) our output -> our decompress == original.
				if got := ourDecompress(t, ours); !bytes.Equal(got, in.data) {
					t.Fatalf("threads=%d: our decode of our output != original", threads)
				}

				// (b) upstream output -> our decompress == original.
				up := upstreamCompress(t, bin, in.data, threads, level)
				if got := ourDecompress(t, up); !bytes.Equal(got, in.data) {
					t.Fatalf("threads=%d: our decode of upstream output != original", threads)
				}
			}

			// (c) single-thread vs multi-thread decompress-equal.
			single := ourDecompress(t, ourCompress(t, in.data, 1, DefaultCompression))
			multi := ourDecompress(t, ourCompress(t, in.data, 4, DefaultCompression))
			if !bytes.Equal(single, multi) {
				t.Fatalf("single- and multi-thread plaintext differ")
			}
		})
	}
}

// fillPseudoRandom fills b with a deterministic pseudo-random byte pattern so
// the test is reproducible without depending on crypto/rand.
func fillPseudoRandom(b []byte) {
	var x uint64 = 0x9e3779b97f4a7c15
	for i := range b {
		x ^= x << 13
		x ^= x >> 7
		x ^= x << 17
		b[i] = byte(x)
	}
}
