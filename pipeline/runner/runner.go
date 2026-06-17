// Package runner executes the parity matrix: for each Entry it resolves our
// freshly built tool binary and the vendored upstream binary, runs both on the
// generated fixtures, compares their output per the entry's compare mode, and
// records timing. Results aggregate into a Report.
package runner

import (
	"bytes"
	"fmt"
	"os/exec"
	"strings"
	"time"

	"github.com/yassineS/bio_ai_experiment/pipeline/fixtures"
	"github.com/yassineS/bio_ai_experiment/pipeline/internal/upstream"
	"github.com/yassineS/bio_ai_experiment/pipeline/matrix"
)

// Status is the outcome category for one entry.
type Status string

const (
	StatusPass    Status = "PASS"    // byte-exact (or within tolerance) match
	StatusSimilar Status = "SIMILAR" // similarity mode matched within tolerance
	StatusDiverge Status = "DIVERGE" // outputs differ (a real failure)
	StatusSkip    Status = "SKIP"    // entry was skipped (Entry.Skip set)
	StatusError   Status = "ERROR"   // a binary failed to run / setup error
)

// Result is the outcome of one matrix entry.
type Result struct {
	Tool         string        `json:"tool"`
	Name         string        `json:"name"`
	Args         []string      `json:"args"`
	Input        string        `json:"input"`
	Compare      string        `json:"compare"`
	Heavy        bool          `json:"heavy"`
	Status       Status        `json:"status"`
	Detail       string        `json:"detail,omitempty"`
	MaxDeviation float64       `json:"max_deviation,omitempty"`
	OursMillis   int64         `json:"ours_ms"`
	UpstreamMs   int64         `json:"upstream_ms"`
	TimingRatio  float64       `json:"timing_ratio,omitempty"` // ours/upstream
	oursDur      time.Duration `json:"-"`
	upDur        time.Duration `json:"-"`
}

// Config controls a run.
type Config struct {
	Manifest *fixtures.Manifest
	CacheDir string // where to build our binaries
	Logf     func(format string, args ...any)
}

func (c Config) log(format string, args ...any) {
	if c.Logf != nil {
		c.Logf(format, args...)
	}
}

// RunEntry executes one entry and returns its result.
func RunEntry(cfg Config, e matrix.Entry) Result {
	res := Result{
		Tool:    e.Tool,
		Name:    e.Name,
		Args:    e.Args,
		Input:   string(e.Input),
		Compare: string(e.CompareModeOrDefault()),
		Heavy:   e.Heavy,
	}
	if e.Skip != "" {
		res.Status = StatusSkip
		res.Detail = e.Skip
		return res
	}

	oursBin, err := upstream.OurBinary(e.Tool, cfg.CacheDir)
	if err != nil {
		res.Status = StatusError
		res.Detail = err.Error()
		return res
	}
	upBin, err := upstream.Binary(e.UpstreamKey())
	if err != nil {
		res.Status = StatusError
		res.Detail = err.Error()
		return res
	}

	args, err := resolvePlaceholders(e.Args, cfg.Manifest)
	if err != nil {
		res.Status = StatusError
		res.Detail = err.Error()
		return res
	}

	// Our invocation: prepend subcommand only when our binary uses it.
	ourArgs := args
	if e.UsesSubcommand && e.Subcommand != "" {
		ourArgs = append([]string{e.Subcommand}, args...)
	}
	// Upstream invocation: prepend the subcommand if present.
	upArgs := args
	if e.Subcommand != "" {
		upArgs = append([]string{e.Subcommand}, args...)
	}

	ourOut, ourErr, ourDur, ourRunErr := timedRun(oursBin, ourArgs)
	upOut, upErr, upDur, upRunErr := timedRun(upBin, upArgs)

	res.oursDur, res.upDur = ourDur, upDur
	res.OursMillis = ourDur.Milliseconds()
	res.UpstreamMs = upDur.Milliseconds()
	if upDur > 0 {
		res.TimingRatio = float64(ourDur) / float64(upDur)
	}

	// A run error on either side that the other did not share is a divergence.
	if (ourRunErr == nil) != (upRunErr == nil) {
		res.Status = StatusDiverge
		res.Detail = fmt.Sprintf("exit mismatch: ours_err=%v upstream_err=%v\nours stderr: %s\nupstream stderr: %s",
			ourRunErr, upRunErr, trunc(string(ourErr)), trunc(string(upErr)))
		return res
	}

	var cmp CompareResult
	switch e.CompareModeOrDefault() {
	case matrix.Similarity:
		cmp = CompareSimilarity(ourOut, upOut)
		if cmp.Equal {
			res.Status = StatusSimilar
		} else {
			res.Status = StatusDiverge
		}
	default: // ByteExact, DirContents (DirContents falls back to byte-exact of stdout for now)
		cmp = CompareByteExact(ourOut, upOut)
		if cmp.Equal {
			res.Status = StatusPass
		} else {
			res.Status = StatusDiverge
		}
	}
	res.MaxDeviation = cmp.MaxDeviation
	if !cmp.Equal {
		res.Detail = cmp.Detail
	}
	return res
}

// timedRun runs a binary with args, returning stdout, stderr, wall-clock, and a
// run error (non-zero exit).
func timedRun(bin string, args []string) (stdout, stderr []byte, dur time.Duration, err error) {
	cmd := exec.Command(bin, args...)
	var out, errb bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = &errb
	start := time.Now()
	err = cmd.Run()
	dur = time.Since(start)
	return out.Bytes(), errb.Bytes(), dur, err
}

// resolvePlaceholders substitutes {bam}, {cram}, {vcf}, ... tokens in args with
// the corresponding fixture path from the manifest.
func resolvePlaceholders(args []string, m *fixtures.Manifest) ([]string, error) {
	out := make([]string, len(args))
	for i, a := range args {
		if !strings.Contains(a, "{") {
			out[i] = a
			continue
		}
		r := a
		for _, key := range []string{"bam", "cram", "vcf", "vcf_plain", "bed", "bed12", "fasta", "genome"} {
			tok := "{" + key + "}"
			if strings.Contains(r, tok) {
				p := m.Path(key)
				if p == "" {
					return nil, fmt.Errorf("fixture %q not in manifest (needed by arg %q)", key, a)
				}
				r = strings.ReplaceAll(r, tok, p)
			}
		}
		out[i] = r
	}
	return out, nil
}
