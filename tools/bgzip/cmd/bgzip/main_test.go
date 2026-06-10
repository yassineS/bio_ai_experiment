package main

import (
	"bytes"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"testing"

	bgzip "github.com/yassineS/bio_ai_experiment/pkg/htsgo/bgzf"
)

// runCLI invokes run with the given args and stdin, returning the exit code and
// captured stdout/stderr.
func runCLI(t *testing.T, args []string, stdin []byte) (code int, stdout, stderr []byte) {
	t.Helper()
	var out, errBuf bytes.Buffer
	code = run(args, bytes.NewReader(stdin), &out, &errBuf)
	return code, out.Bytes(), errBuf.Bytes()
}

// decodeBGZF decompresses a BGZF stream into plaintext.
func decodeBGZF(t *testing.T, data []byte) []byte {
	t.Helper()
	r, err := bgzip.NewReader(bytes.NewReader(data))
	if err != nil {
		t.Fatalf("NewReader: %v", err)
	}
	var plain bytes.Buffer
	if _, err := plain.ReadFrom(r); err != nil {
		t.Fatalf("read: %v", err)
	}
	return plain.Bytes()
}

func TestCompressStdinToStdout(t *testing.T) {
	in := []byte("the quick brown fox\n")
	code, stdout, stderr := runCLI(t, []string{"-c"}, in)
	if code != 0 {
		t.Fatalf("exit=%d stderr=%s", code, stderr)
	}
	if got := decodeBGZF(t, stdout); !bytes.Equal(got, in) {
		t.Fatalf("roundtrip mismatch: %q", got)
	}
}

func TestCompressStdinThreadsToStdout(t *testing.T) {
	in := bytes.Repeat([]byte("ACGT\n"), 50000)
	for _, threads := range []string{"1", "2", "4"} {
		code, stdout, stderr := runCLI(t, []string{"-c", "-@", threads}, in)
		if code != 0 {
			t.Fatalf("threads=%s exit=%d stderr=%s", threads, code, stderr)
		}
		if got := decodeBGZF(t, stdout); !bytes.Equal(got, in) {
			t.Fatalf("threads=%s roundtrip mismatch", threads)
		}
		if !bytes.HasSuffix(stdout, bgzip.EOFBlock) {
			t.Fatalf("threads=%s missing EOF block", threads)
		}
	}
}

// TestStdinOutputRename verifies the -o/--output naming convention: reading
// stdin and naming the output file explicitly writes BGZF to that file rather
// than stdout.
func TestStdinOutputRename(t *testing.T) {
	dir := t.TempDir()
	outPath := filepath.Join(dir, "renamed.gz")
	in := []byte("named output from stdin\n")

	code, stdout, stderr := runCLI(t, []string{"-o", outPath}, in)
	if code != 0 {
		t.Fatalf("exit=%d stderr=%s", code, stderr)
	}
	if len(stdout) != 0 {
		t.Fatalf("expected no stdout, got %d bytes", len(stdout))
	}
	data, err := os.ReadFile(outPath)
	if err != nil {
		t.Fatalf("read output: %v", err)
	}
	if got := decodeBGZF(t, data); !bytes.Equal(got, in) {
		t.Fatalf("roundtrip mismatch: %q", got)
	}
}

// TestOutputDashMeansStdout verifies `-o -` is treated as -c.
func TestOutputDashMeansStdout(t *testing.T) {
	in := []byte("dash means stdout\n")
	code, stdout, stderr := runCLI(t, []string{"-o", "-"}, in)
	if code != 0 {
		t.Fatalf("exit=%d stderr=%s", code, stderr)
	}
	if got := decodeBGZF(t, stdout); !bytes.Equal(got, in) {
		t.Fatalf("roundtrip mismatch: %q", got)
	}
}

func TestOutputAndStdoutConflict(t *testing.T) {
	code, _, stderr := runCLI(t, []string{"-c", "-o", "x.gz"}, []byte("hi"))
	if code != 2 {
		t.Fatalf("expected exit 2, got %d", code)
	}
	if !strings.Contains(string(stderr), "stdout") {
		t.Fatalf("expected conflict message, got %q", stderr)
	}
}

// TestCompressNamedFileWithThreads exercises in-place compression of a real
// file with the multi-threaded path, then decompresses with -o to a new file.
func TestCompressNamedFileWithThreads(t *testing.T) {
	dir := t.TempDir()
	src := filepath.Join(dir, "data.txt")
	content := bytes.Repeat([]byte("line of text\n"), 20000)
	if err := os.WriteFile(src, content, 0o644); err != nil {
		t.Fatalf("write src: %v", err)
	}

	code, _, stderr := runCLI(t, []string{"-@", "4", "-k", src}, nil)
	if code != 0 {
		t.Fatalf("compress exit=%d stderr=%s", code, stderr)
	}
	// -k keeps the source.
	if _, err := os.Stat(src); err != nil {
		t.Fatalf("source removed despite -k: %v", err)
	}
	gz := src + ".gz"
	data, err := os.ReadFile(gz)
	if err != nil {
		t.Fatalf("read gz: %v", err)
	}
	if got := decodeBGZF(t, data); !bytes.Equal(got, content) {
		t.Fatalf("roundtrip mismatch")
	}

	// Decompress to an explicit output path.
	dst := filepath.Join(dir, "restored.txt")
	code, _, stderr = runCLI(t, []string{"-d", "-k", "-o", dst, gz}, nil)
	if code != 0 {
		t.Fatalf("decompress exit=%d stderr=%s", code, stderr)
	}
	restored, err := os.ReadFile(dst)
	if err != nil {
		t.Fatalf("read restored: %v", err)
	}
	if !bytes.Equal(restored, content) {
		t.Fatalf("decompressed content mismatch")
	}
}

func TestDecompressStdinToStdout(t *testing.T) {
	// First compress to get a BGZF stream.
	orig := []byte("decompress me\n")
	_, compressed, _ := runCLI(t, []string{"-c"}, orig)

	code, stdout, stderr := runCLI(t, []string{"-d", "-c"}, compressed)
	if code != 0 {
		t.Fatalf("exit=%d stderr=%s", code, stderr)
	}
	if !bytes.Equal(stdout, orig) {
		t.Fatalf("got %q, want %q", stdout, orig)
	}
}

func TestThreadsMustBePositive(t *testing.T) {
	code, _, stderr := runCLI(t, []string{"-@", "0", "-c"}, []byte("x"))
	if code != 2 {
		t.Fatalf("expected exit 2, got %d", code)
	}
	if !strings.Contains(string(stderr), "threads") {
		t.Fatalf("expected threads error, got %q", stderr)
	}
}

// TestQuerySubcommands exercises -s (size), -r (reindex), and -b (offset) on a
// real compressed file.
func TestQuerySubcommands(t *testing.T) {
	dir := t.TempDir()
	src := filepath.Join(dir, "q.txt")
	content := bytes.Repeat([]byte("payload\n"), 30000)
	if err := os.WriteFile(src, content, 0o644); err != nil {
		t.Fatalf("write: %v", err)
	}
	gz := filepath.Join(dir, "q.txt.gz")
	// Compress to a named file via -o with the multi-threaded path.
	if code, _, stderr := runCLI(t, []string{"-k", "-@", "2", "-o", gz, src}, nil); code != 0 {
		t.Fatalf("compress -o exit=%d %s", code, stderr)
	}

	// -s prints the decompressed size.
	code, stdout, stderr := runCLI(t, []string{"-s", gz}, nil)
	if code != 0 {
		t.Fatalf("size exit=%d %s", code, stderr)
	}
	if strings.TrimSpace(string(stdout)) != itoa(len(content)) {
		t.Fatalf("size = %q, want %d", stdout, len(content))
	}

	// -r writes a .gzi index.
	if code, _, stderr := runCLI(t, []string{"-r", gz}, nil); code != 0 {
		t.Fatalf("reindex exit=%d %s", code, stderr)
	}
	if _, err := os.Stat(gz + ".gzi"); err != nil {
		t.Fatalf("gzi not written: %v", err)
	}

	// -b 0 maps compressed offset 0 to uncompressed offset 0.
	code, stdout, stderr = runCLI(t, []string{"-b", "0", gz}, nil)
	if code != 0 {
		t.Fatalf("offset exit=%d %s", code, stderr)
	}
	if strings.TrimSpace(string(stdout)) != "0" {
		t.Fatalf("offset 0 = %q, want 0", stdout)
	}
}

func itoa(n int) string {
	return strconv.Itoa(n)
}

func TestInvalidCompressLevel(t *testing.T) {
	code, _, stderr := runCLI(t, []string{"-l", "99", "-c"}, []byte("x"))
	if code != 2 {
		t.Fatalf("expected exit 2, got %d (%s)", code, stderr)
	}
}

func TestTooManyInputs(t *testing.T) {
	code, _, stderr := runCLI(t, []string{"a", "b"}, nil)
	if code != 2 {
		t.Fatalf("expected exit 2, got %d (%s)", code, stderr)
	}
}

func TestHelpAndVersion(t *testing.T) {
	code, stdout, _ := runCLI(t, []string{"-h"}, nil)
	if code != 0 || !strings.Contains(string(stdout), "Block-gzip") {
		t.Fatalf("help: code=%d out=%q", code, stdout)
	}
	code, stdout, _ = runCLI(t, []string{"-v"}, nil)
	if code != 0 || !strings.Contains(string(stdout), "bgzip") {
		t.Fatalf("version: code=%d out=%q", code, stdout)
	}
}
