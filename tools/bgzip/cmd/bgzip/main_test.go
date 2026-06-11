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

// TestIntegrityCheckValidStream verifies that `bgzip --test` on a well-formed
// BGZF file decompresses cleanly, writes nothing to stdout, and exits 0
// without touching the input.
func TestIntegrityCheckValidStream(t *testing.T) {
	dir := t.TempDir()
	plain := bytes.Repeat([]byte("ACGTACGTNN\n"), 20000)

	// Produce a real multi-block BGZF file via the compress path.
	gzPath := filepath.Join(dir, "data.gz")
	if code, _, stderr := runCLI(t, []string{"-@", "2", "-o", gzPath}, plain); code != 0 {
		t.Fatalf("compress exit=%d stderr=%s", code, stderr)
	}
	before, err := os.ReadFile(gzPath)
	if err != nil {
		t.Fatalf("read gz: %v", err)
	}

	code, stdout, stderr := runCLI(t, []string{"--test", gzPath}, nil)
	if code != 0 {
		t.Fatalf("--test on valid stream exit=%d stderr=%s", code, stderr)
	}
	if len(stdout) != 0 {
		t.Fatalf("--test wrote %d bytes to stdout; want none", len(stdout))
	}
	after, err := os.ReadFile(gzPath)
	if err != nil {
		t.Fatalf("re-read gz: %v", err)
	}
	if !bytes.Equal(before, after) {
		t.Fatalf("--test modified the input file")
	}
}

// TestIntegrityCheckCorruptStream verifies that `bgzip --test` exits non-zero
// when the compressed payload is corrupted (CRC / inflate failure), and via
// stdin too (no on-disk file needed).
func TestIntegrityCheckCorruptStream(t *testing.T) {
	plain := bytes.Repeat([]byte("the quick brown fox\n"), 5000)
	code, stdout, _ := runCLI(t, []string{"-c"}, plain)
	if code != 0 {
		t.Fatalf("compress exit=%d", code)
	}
	good := stdout

	// Flip bytes in the middle of the first block's deflate payload (past the
	// 18-byte BGZF header) so the block decodes to wrong data / bad CRC.
	corrupt := append([]byte(nil), good...)
	if len(corrupt) < 40 {
		t.Fatalf("compressed output unexpectedly short: %d bytes", len(corrupt))
	}
	for i := 20; i < 30; i++ {
		corrupt[i] ^= 0xff
	}

	code, _, _ = runCLI(t, []string{"--test", "-"}, corrupt)
	if code == 0 {
		t.Fatalf("--test accepted a corrupted stream; want non-zero exit")
	}
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
