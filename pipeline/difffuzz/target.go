package difffuzz

import (
	"bytes"
	"context"
	"errors"
	"os"
	"os/exec"
	"path/filepath"
	"time"

	"github.com/yassineS/bio_ai_experiment/pipeline/internal/upstream"
)

// Target is one tool/subcommand invocation the fuzzer differentially tests.
//
// Both binaries are run with the SAME argument template. The token "{in}" in
// Args is replaced with the path of a temp file holding the fuzzed input; if no
// "{in}" token is present the input is piped on stdin instead. Tools whose CLI
// shape differs between our port and upstream can override per side with
// OurArgs / UpstreamArgs (same {in} substitution applies).
type Target struct {
	// Name is a unique label for reports and -target filtering.
	Name string

	// Tool is our tool binary name (tools/<Tool>/cmd/<Tool>) and, unless
	// UpstreamKey overrides it, the upstream binary key.
	Tool string

	// UpstreamKey overrides the upstream binary lookup key when our tool name
	// differs from the upstream binary (e.g. our "bedintersect" → "bedtools").
	UpstreamKey string

	// Subcommand is the upstream subcommand (e.g. "view", "flagstat"); prepended
	// to the upstream args always, and to our args when UsesSubcommand is true.
	Subcommand string

	// UsesSubcommand reports whether OUR binary also takes Subcommand as its
	// first argument.
	UsesSubcommand bool

	// Args is the shared argument template (with "{in}"). OurArgs / UpstreamArgs
	// override it per side when non-nil.
	Args         []string
	OurArgs      []string
	UpstreamArgs []string

	// Format is the input format; selects the structured generator and the seed
	// fixture key.
	Format Format

	// SeedFixture is the manifest fixture key whose bytes seed the mutation
	// strategy (e.g. "vcf_plain", "bam"). Empty means no seed (mutation falls
	// back to perturbing raw random bytes, fine for RawBytes targets).
	SeedFixture string
}

func (t Target) upstreamKey() string {
	if t.UpstreamKey != "" {
		return t.UpstreamKey
	}
	return t.Tool
}

// resolvedTarget is a Target with its binaries located.
type resolvedTarget struct {
	Target
	ourBin string
	upBin  string
}

// resolve locates both binaries for the target, returning an error (used to
// SKIP the target) if either is unavailable.
func (t Target) resolve(cacheDir string) (resolvedTarget, error) {
	ourBin, err := upstream.OurBinary(t.Tool, cacheDir)
	if err != nil {
		return resolvedTarget{}, err
	}
	upBin, err := upstream.Binary(t.upstreamKey())
	if err != nil {
		return resolvedTarget{}, err
	}
	return resolvedTarget{Target: t, ourBin: ourBin, upBin: upBin}, nil
}

// run executes one side of the target on input, returning its outcome. The
// input is written to a temp file (if the template uses "{in}") or piped on
// stdin. timeout bounds the run; exceeding it sets TimedOut.
func runOneSide(bin string, template []string, subcommand string, prependSub bool,
	input []byte, dir string, timeout time.Duration, env []string) RunOutcome {
	args, usesFile, inputPath := resolveArgs(template, dir, input)
	if prependSub && subcommand != "" {
		args = append([]string{subcommand}, args...)
	}

	ctx, cancel := context.WithTimeout(context.Background(), timeout)
	defer cancel()
	cmd := exec.CommandContext(ctx, bin, args...)
	if env != nil {
		cmd.Env = env
	}
	if !usesFile {
		cmd.Stdin = bytes.NewReader(input)
	}
	var out, errb bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = &errb
	err := cmd.Run()

	res := RunOutcome{Stdout: out.Bytes(), Stderr: errb.Bytes(), inputPath: inputPath}
	if ctx.Err() == context.DeadlineExceeded {
		res.TimedOut = true
		return res
	}
	if err == nil {
		res.ExitCode = 0
		return res
	}
	var ee *exec.ExitError
	if errors.As(err, &ee) {
		res.ExitCode = ee.ExitCode()
		// A negative ExitCode (-1) means the process was terminated by a signal
		// (e.g. SIGSEGV / SIGABRT) rather than exiting normally: a crash.
		if res.ExitCode < 0 || (ee.ProcessState != nil && ee.ProcessState.ExitCode() < 0) {
			res.Crashed = true
		}
		return res
	}
	// Could not start the process at all (e.g. ENOENT). Treat as a crash so it
	// is surfaced rather than silently dropped.
	res.Crashed = true
	res.ExitCode = -1
	return res
}

// resolveArgs substitutes "{in}" in template with a temp-file path holding
// input (written under dir). It reports whether a file was used (usesFile) and
// the substituted temp-file path (inputPath, empty for stdin runs). The input
// is written once even when "{in}" appears multiple times (e.g. -a {in} -b
// {in}), so both occurrences point at the same file.
func resolveArgs(template []string, dir string, input []byte) (args []string, usesFile bool, inputPath string) {
	args = make([]string, len(template))
	for i, a := range template {
		if a == "{in}" {
			if inputPath == "" {
				inputPath = writeTemp(dir, input)
			}
			args[i] = inputPath
			usesFile = true
			continue
		}
		args[i] = a
	}
	return args, usesFile, inputPath
}

// writeTemp writes input to a uniquely named file under dir and returns its
// path. Errors are unlikely (dir is a fresh temp dir) and are swallowed into an
// empty path, which makes the run fail loudly downstream rather than panic.
func writeTemp(dir string, input []byte) string {
	f, err := os.CreateTemp(dir, "in-*")
	if err != nil {
		return filepath.Join(dir, "in")
	}
	defer f.Close()
	_, _ = f.Write(input)
	return f.Name()
}

// execute runs BOTH binaries on input and returns their outcomes plus the
// divergence classification. ourEnv, when non-nil, is the environment our
// (possibly coverage-instrumented) binary runs with; upstream always uses the
// inherited environment.
func (rt resolvedTarget) execute(input []byte, dir string, timeout time.Duration, ourEnv []string) (ours, up RunOutcome, class DivergenceClass, detail string) {
	ours = runOneSide(rt.ourBin, rt.ourTemplate(), rt.Subcommand, rt.UsesSubcommand, input, dir, timeout, ourEnv)
	up = runOneSide(rt.upBin, rt.upTemplate(), rt.Subcommand, true, input, dir, timeout, nil)
	// Each side wrote the input to its own temp file; a tool that echoes that
	// path in a diagnostic would otherwise diverge on the path token alone.
	// Rewrite each side's own path to a common placeholder before classifying.
	ours.Stdout = normalizePath(ours.Stdout, ours.inputPath)
	ours.Stderr = normalizePath(ours.Stderr, ours.inputPath)
	up.Stdout = normalizePath(up.Stdout, up.inputPath)
	up.Stderr = normalizePath(up.Stderr, up.inputPath)
	class, detail = Diff(ours, up)
	return ours, up, class, detail
}

func (rt resolvedTarget) ourTemplate() []string {
	if rt.OurArgs != nil {
		return rt.OurArgs
	}
	return rt.Args
}

func (rt resolvedTarget) upTemplate() []string {
	if rt.UpstreamArgs != nil {
		return rt.UpstreamArgs
	}
	return rt.Args
}
