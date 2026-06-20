package difffuzz

import (
	"bytes"
	"fmt"

	"github.com/yassineS/bio_ai_experiment/pipeline/runner"
)

// DivergenceClass categorizes how two runs of the same input differed. The
// classes are ordered by severity for reporting: a crash mismatch outranks an
// exit-code-only difference, which outranks a stderr/stdout text difference.
type DivergenceClass string

// The divergence classes. None means the two runs agreed (after normalization).
const (
	ClassNone          DivergenceClass = "none"
	ClassStdoutDiffers DivergenceClass = "stdout-differs"
	ClassStderrDiffers DivergenceClass = "stderr-differs"
	ClassExitDiffers   DivergenceClass = "exitcode-differs"
	ClassOneCrashed    DivergenceClass = "one-crashed"
	ClassBothCrashed   DivergenceClass = "both-crashed"
)

// RunOutcome captures one binary's result on one input.
type RunOutcome struct {
	Stdout   []byte
	Stderr   []byte
	ExitCode int
	// Crashed reports an abnormal termination (signal / panic), distinct from a
	// clean non-zero exit. For drop-in parity a clean "exit 1 with an error
	// message" is fine; a SIGSEGV is never acceptable.
	Crashed bool
	// TimedOut reports the run exceeded the per-invocation deadline; treated like
	// a crash for classification (it indicates a hang on adversarial input).
	TimedOut bool
	// inputPath is the temp file the input was written to for this side (empty
	// for stdin runs). Each side gets a distinct temp file, so tools that echo
	// their -i argument in an error message (e.g. bedtools merge's sort/field
	// diagnostics) would otherwise diverge purely on the path token. The
	// per-target driver normalizes this path to a fixed placeholder before
	// classifying, so only a genuine wording difference is compared.
	inputPath string
}

// inputPathPlaceholder is the token both sides' input file paths are rewritten
// to before stderr/stdout comparison, so a diagnostic that legitimately echoes
// the input file name does not register as a divergence on the path alone.
const inputPathPlaceholder = "<INPUT>"

// normalizePath returns a copy of out with every occurrence of path (the side's
// own temp input file) replaced by inputPathPlaceholder. A nil/empty path or
// output is returned unchanged.
func normalizePath(out []byte, path string) []byte {
	if len(out) == 0 || path == "" {
		return out
	}
	return bytes.ReplaceAll(out, []byte(path), []byte(inputPathPlaceholder))
}

// terminal reports whether the outcome ended abnormally (crash or timeout).
func (o RunOutcome) terminal() bool { return o.Crashed || o.TimedOut }

// Diff classifies the divergence between our outcome and upstream's.
//
// stdout and stderr are compared after applying runner.StripProvenance (the
// SAME normalization the parity harness uses) so version stamps and @PG/##*
// command lines do not count as divergences. Exit codes are compared exactly
// for the clean-exit case. A normalized-equal pair classifies as ClassNone.
//
// Precedence (most severe first): both-crashed, one-crashed, exitcode-differs,
// stdout-differs, stderr-differs. stdout is ranked above stderr because a
// stdout content divergence is a data bug, whereas a stderr-only divergence is
// usually a diagnostics-wording difference (still reported, just lower).
func Diff(ours, up RunOutcome) (DivergenceClass, string) {
	switch {
	case ours.terminal() && up.terminal():
		// Both blew up. This is the "both-crashed" bucket: not a parity FAILURE
		// in the drop-in sense (both reject the garbage), but tracked because a
		// crash (vs clean error) is itself a robustness concern.
		return ClassBothCrashed, fmt.Sprintf("both terminated abnormally (ours crashed=%v timeout=%v, upstream crashed=%v timeout=%v)",
			ours.Crashed, ours.TimedOut, up.Crashed, up.TimedOut)
	case ours.terminal() != up.terminal():
		who := "ours"
		if up.terminal() {
			who = "upstream"
		}
		return ClassOneCrashed, fmt.Sprintf("%s terminated abnormally while the other did not (ours: crashed=%v timeout=%v exit=%d; upstream: crashed=%v timeout=%v exit=%d)",
			who, ours.Crashed, ours.TimedOut, ours.ExitCode, up.Crashed, up.TimedOut, up.ExitCode)
	}

	// Neither crashed: compare the clean results. Exit-code parity first.
	if ours.ExitCode != up.ExitCode {
		return ClassExitDiffers, fmt.Sprintf("exit code differs: ours=%d upstream=%d\n  ours stderr: %s\n  upstream stderr: %s",
			ours.ExitCode, up.ExitCode, truncStr(ours.Stderr), truncStr(up.Stderr))
	}

	// stdout content (normalized).
	oOut := runner.StripProvenance(ours.Stdout)
	uOut := runner.StripProvenance(up.Stdout)
	if !bytes.Equal(oOut, uOut) {
		return ClassStdoutDiffers, fmt.Sprintf("stdout differs (normalized): %s", firstLineDiff(oOut, uOut))
	}

	// stderr content (normalized).
	oErr := runner.StripProvenance(ours.Stderr)
	uErr := runner.StripProvenance(up.Stderr)
	if !bytes.Equal(oErr, uErr) {
		return ClassStderrDiffers, fmt.Sprintf("stderr differs (normalized): %s", firstLineDiff(oErr, uErr))
	}

	return ClassNone, ""
}

// IsDivergence reports whether a class represents a reportable divergence.
// ClassBothCrashed is included: while it is not a drop-in parity failure (both
// reject the input), a crash on adversarial input is a robustness finding worth
// surfacing. Callers that only care about strict parity can filter it out.
func IsDivergence(c DivergenceClass) bool { return c != ClassNone }

// firstLineDiff returns a short human description of the first differing line
// between two normalized streams.
func firstLineDiff(a, b []byte) string {
	al := bytes.Split(a, []byte("\n"))
	bl := bytes.Split(b, []byte("\n"))
	n := len(al)
	if len(bl) < n {
		n = len(bl)
	}
	for i := 0; i < n; i++ {
		if !bytes.Equal(al[i], bl[i]) {
			return fmt.Sprintf("first diff at line %d:\n    ours:     %s\n    upstream: %s",
				i+1, truncStr(al[i]), truncStr(bl[i]))
		}
	}
	if len(al) != len(bl) {
		return fmt.Sprintf("line count differs: ours=%d upstream=%d", len(al), len(bl))
	}
	return "streams differ"
}

func truncStr(b []byte) string {
	const max = 160
	s := string(b)
	if len(s) > max {
		return s[:max] + "..."
	}
	return s
}
