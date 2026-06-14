package bcftools

import (
	"bytes"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"strings"
	"testing"
)

// pluginVCF is a small VCF whose second record carries a soft FILTER ("q10").
const pluginVCF = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=200>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FILTER=<ID=q10,Description="below 10">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	100	rs1	A	T	30	PASS	DP=50	GT	0/1
chr1	150	rs2	C	G	10	q10	DP=15	GT	0/0
`

// buildExamplePlugin compiles the reference example plugin into dir under the
// given name and returns nothing; it fails the test on a build error.
func buildExamplePlugin(t *testing.T, dir, name string) {
	t.Helper()
	src, err := filepath.Abs("../../plugins/example")
	if err != nil {
		t.Fatalf("resolving example plugin source: %v", err)
	}
	bin := filepath.Join(dir, name)
	cmd := exec.Command("go", "build", "-o", bin, src)
	var stderr bytes.Buffer
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("building example plugin: %v\n%s", err, stderr.String())
	}
}

// writeScriptPlugin writes a tiny shell-script plugin (a plain `cat` filter)
// and marks it executable. Callers gate on requirePOSIXShell first.
func writeScriptPlugin(t *testing.T, dir, name, body string) {
	t.Helper()
	path := filepath.Join(dir, name)
	if err := os.WriteFile(path, []byte(body), 0o755); err != nil {
		t.Fatalf("writing script plugin: %v", err)
	}
}

// requirePOSIXShell asserts a POSIX `/bin/sh` is reachable so the
// shell-script plugin fixtures can run. Per the env-guard policy
// (PR #294) the absence of the dependency is a loud t.Fatalf rather than
// a silent skip everywhere a POSIX shell legitimately exists. On Windows
// — where there is genuinely no `/bin/sh` — this is a true platform
// guard and is allowed to skip.
func requirePOSIXShell(t *testing.T) {
	t.Helper()
	if runtime.GOOS == "windows" {
		t.Skip("shell-script plugin fixture needs a POSIX shell; not available on windows")
	}
	if _, err := exec.LookPath("sh"); err != nil {
		t.Fatalf("POSIX shell `sh` not found on PATH; the shell-script plugin parity fixtures require it: %v", err)
	}
}

func TestRunPluginRoundTrip(t *testing.T) {
	dir := t.TempDir()
	buildExamplePlugin(t, dir, "example")
	t.Setenv(pluginEnvVar, dir)

	inFile := filepath.Join(dir, "in.vcf")
	if err := os.WriteFile(inFile, []byte(pluginVCF), 0o644); err != nil {
		t.Fatalf("writing input: %v", err)
	}

	var out, stderr bytes.Buffer
	err := RunPlugin(PluginOptions{
		Name:         "example",
		InputFile:    inFile,
		OutputFormat: OutputVCF,
	}, &out, &stderr)
	if err != nil {
		t.Fatalf("RunPlugin: %v", err)
	}

	got := out.String()
	// The example plugin resets FILTER to PASS on every record.
	if strings.Contains(got, "\tq10\t") {
		t.Errorf("expected q10 FILTER to be cleared, got:\n%s", got)
	}
	if n := strings.Count(got, "\tPASS\t"); n != 2 {
		t.Errorf("expected 2 PASS records, got %d:\n%s", n, got)
	}
	if !strings.Contains(got, "rs1") || !strings.Contains(got, "rs2") {
		t.Errorf("records lost in round-trip:\n%s", got)
	}
	if !strings.Contains(stderr.String(), "processed 2 record") {
		t.Errorf("plugin stderr not forwarded, got: %q", stderr.String())
	}
}

// TestRunPluginScriptPlugin verifies that a plugin written in another
// language (here, a POSIX shell `cat` filter) round-trips identically: this
// is the language-agnostic guarantee of the subprocess protocol.
func TestRunPluginScriptPlugin(t *testing.T) {
	requirePOSIXShell(t)
	dir := t.TempDir()
	writeScriptPlugin(t, dir, "passthru", "#!/bin/sh\nexec cat\n")
	t.Setenv(pluginEnvVar, dir)

	inFile := filepath.Join(dir, "in.vcf")
	if err := os.WriteFile(inFile, []byte(pluginVCF), 0o644); err != nil {
		t.Fatalf("writing input: %v", err)
	}

	var out, stderr bytes.Buffer
	err := RunPlugin(PluginOptions{
		Name:         "passthru",
		InputFile:    inFile,
		OutputFormat: OutputVCF,
	}, &out, &stderr)
	if err != nil {
		t.Fatalf("RunPlugin: %v", err)
	}
	if !strings.Contains(out.String(), "rs1") || !strings.Contains(out.String(), "rs2") {
		t.Errorf("pass-through plugin lost records:\n%s", out.String())
	}
}

// TestRunPluginLargeInput streams a plugin over ~50k VCF records and confirms
// the data round-trips intact. It is a deadlock canary: the host currently
// buffers the whole input before spawning the child, but a future streaming
// rewrite that pipes stdin and stdout concurrently could reintroduce a
// classic pipe deadlock (child blocked writing stdout while the host is
// still blocked writing stdin). An input large enough to overflow the OS
// pipe buffer (64 KiB on Linux) would hang such a rewrite; this test would
// then time out instead of silently regressing.
func TestRunPluginLargeInput(t *testing.T) {
	requirePOSIXShell(t)
	dir := t.TempDir()
	// A plain `cat` filter: streams every byte of stdin straight to stdout.
	writeScriptPlugin(t, dir, "passthru", "#!/bin/sh\nexec cat\n")
	t.Setenv(pluginEnvVar, dir)

	const nRecords = 50000
	var buf bytes.Buffer
	buf.WriteString("##fileformat=VCFv4.2\n")
	buf.WriteString("##contig=<ID=chr1,length=300000000>\n")
	buf.WriteString("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">\n")
	buf.WriteString("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GT\">\n")
	buf.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
	for i := 0; i < nRecords; i++ {
		pos := i + 1
		fmt.Fprintf(&buf, "chr1\t%d\trs%d\tA\tT\t30\tPASS\tDP=%d\tGT\t0/1\n", pos, i, i%100)
	}

	inFile := filepath.Join(dir, "big.vcf")
	if err := os.WriteFile(inFile, buf.Bytes(), 0o644); err != nil {
		t.Fatalf("writing input: %v", err)
	}

	var out, stderr bytes.Buffer
	err := RunPlugin(PluginOptions{
		Name:         "passthru",
		InputFile:    inFile,
		OutputFormat: OutputVCF,
	}, &out, &stderr)
	if err != nil {
		t.Fatalf("RunPlugin: %v", err)
	}

	if n := strings.Count(out.String(), "\tPASS\t"); n != nRecords {
		t.Errorf("expected %d records to round-trip, got %d", nRecords, n)
	}
	if !strings.Contains(out.String(), "rs0\t") ||
		!strings.Contains(out.String(), fmt.Sprintf("rs%d\t", nRecords-1)) {
		t.Errorf("first/last record lost in large-input round-trip")
	}
}

func TestRunPluginMissing(t *testing.T) {
	t.Setenv(pluginEnvVar, t.TempDir())
	var out, stderr bytes.Buffer
	err := RunPlugin(PluginOptions{Name: "no-such-plugin"}, &out, &stderr)
	if err == nil {
		t.Fatal("expected an error for a missing plugin")
	}
	var nf *PluginNotFoundError
	if !asPluginNotFound(err, &nf) {
		t.Fatalf("expected *PluginNotFoundError, got %T: %v", err, err)
	}
}

func TestRunPluginNonZeroExit(t *testing.T) {
	requirePOSIXShell(t)
	dir := t.TempDir()
	writeScriptPlugin(t, dir, "boom", "#!/bin/sh\necho 'plugin blew up' >&2\nexit 3\n")
	t.Setenv(pluginEnvVar, dir)

	inFile := filepath.Join(dir, "in.vcf")
	if err := os.WriteFile(inFile, []byte(pluginVCF), 0o644); err != nil {
		t.Fatalf("writing input: %v", err)
	}

	var out, stderr bytes.Buffer
	err := RunPlugin(PluginOptions{
		Name:         "boom",
		InputFile:    inFile,
		OutputFormat: OutputVCF,
	}, &out, &stderr)
	if err == nil {
		t.Fatal("expected an error for a non-zero plugin exit")
	}
	var ee *PluginExecError
	if !asPluginExec(err, &ee) {
		t.Fatalf("expected *PluginExecError, got %T: %v", err, err)
	}
	// Plugin stderr is streamed through to the host stderr.
	if !strings.Contains(stderr.String(), "plugin blew up") {
		t.Errorf("plugin stderr not forwarded, got: %q", stderr.String())
	}
}

func TestListPlugins(t *testing.T) {
	dir := t.TempDir()
	buildExamplePlugin(t, dir, "example")
	// A non-executable file in the directory must be ignored.
	if err := os.WriteFile(filepath.Join(dir, "readme.txt"), []byte("hi"), 0o644); err != nil {
		t.Fatalf("writing non-exec file: %v", err)
	}
	t.Setenv(pluginEnvVar, dir)

	plugins := ListPlugins(false)
	if len(plugins) != 1 || plugins[0].Name != "example" {
		t.Fatalf("expected just the example plugin, got %+v", plugins)
	}
	if got := FormatPluginList(plugins, false); strings.TrimSpace(got) != "example" {
		t.Errorf("plain listing = %q, want %q", got, "example")
	}

	// Verbose probing must populate About via the optional --about flag.
	verbose := ListPlugins(true)
	if len(verbose) != 1 || !strings.Contains(verbose[0].About, "reference plugin") {
		t.Fatalf("verbose listing missing --about line: %+v", verbose)
	}
	if out := FormatPluginList(verbose, true); !strings.Contains(out, verbose[0].Path) {
		t.Errorf("verbose listing should include the path: %q", out)
	}
}

func TestPluginDirsEmpty(t *testing.T) {
	t.Setenv(pluginEnvVar, "")
	if dirs := PluginDirs(); dirs != nil {
		t.Errorf("expected nil dirs when %s is unset, got %v", pluginEnvVar, dirs)
	}
	if _, err := ResolvePlugin("anything"); err == nil {
		t.Error("expected ResolvePlugin to fail when no dirs are configured")
	}
}

// asPluginNotFound and asPluginExec are tiny errors.As wrappers kept local so
// the test file does not need to import the errors package directly.
func asPluginNotFound(err error, target **PluginNotFoundError) bool {
	for err != nil {
		if e, ok := err.(*PluginNotFoundError); ok {
			*target = e
			return true
		}
		u, ok := err.(interface{ Unwrap() error })
		if !ok {
			return false
		}
		err = u.Unwrap()
	}
	return false
}

func asPluginExec(err error, target **PluginExecError) bool {
	for err != nil {
		if e, ok := err.(*PluginExecError); ok {
			*target = e
			return true
		}
		u, ok := err.(interface{ Unwrap() error })
		if !ok {
			return false
		}
		err = u.Unwrap()
	}
	return false
}
