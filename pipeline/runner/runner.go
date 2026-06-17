// Package runner executes the parity matrix: for each Entry it resolves our
// freshly built tool binary and the vendored upstream binary, runs both on the
// generated fixtures, compares their output per the entry's compare mode, and
// records timing. Results aggregate into a Report.
package runner

import (
	"bytes"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"
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

	// When the entry writes to an output prefix (vcftools/mosdepth), each side
	// gets its own temp dir and {out} resolves to "<dir>/out". Otherwise both
	// sides share the resolved args and the runner compares stdout.
	usesOutPrefix := len(e.OutputFiles) > 0
	var ourDir, upDir string
	if usesOutPrefix {
		ourDir, upDir, err = mkOutDirs(cfg.CacheDir)
		if err != nil {
			res.Status = StatusError
			res.Detail = err.Error()
			return res
		}
		defer os.RemoveAll(ourDir)
		defer os.RemoveAll(upDir)
	}

	// Per-side argument templates: OurArgs / UpstreamArgs override the shared
	// Args for tools whose CLI shape differs from upstream's.
	ourTemplate := e.Args
	if e.OurArgs != nil {
		ourTemplate = e.OurArgs
	}
	upTemplate := e.Args
	if e.UpstreamArgs != nil {
		upTemplate = e.UpstreamArgs
	}

	ourArgs, err := resolvePlaceholders(ourTemplate, cfg.Manifest, filepath.Join(ourDir, "out"))
	if err != nil {
		res.Status = StatusError
		res.Detail = err.Error()
		return res
	}
	upArgs, err := resolvePlaceholders(upTemplate, cfg.Manifest, filepath.Join(upDir, "out"))
	if err != nil {
		res.Status = StatusError
		res.Detail = err.Error()
		return res
	}

	// Our invocation: prepend subcommand only when our binary uses it.
	if e.UsesSubcommand && e.Subcommand != "" {
		ourArgs = append([]string{e.Subcommand}, ourArgs...)
	}
	// Upstream invocation: prepend the subcommand if present.
	if e.Subcommand != "" {
		upArgs = append([]string{e.Subcommand}, upArgs...)
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
	switch {
	case usesOutPrefix:
		// Compare the named output files between the two prefixes (decompressing
		// .gz, stripping provenance). The compare mode still selects byte-exact
		// vs similarity for the file contents.
		cmp = CompareOutputFiles(filepath.Join(ourDir, "out"), filepath.Join(upDir, "out"),
			e.OutputFiles, e.CompareModeOrDefault())
		if cmp.Equal {
			if e.CompareModeOrDefault() == matrix.Similarity {
				res.Status = StatusSimilar
			} else {
				res.Status = StatusPass
			}
		} else {
			res.Status = StatusDiverge
		}
	case e.CompareModeOrDefault() == matrix.Similarity:
		cmp = CompareSimilarity(ourOut, upOut)
		if cmp.Equal {
			res.Status = StatusSimilar
		} else {
			res.Status = StatusDiverge
		}
	default: // ByteExact on stdout
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

// mkOutDirs creates a fresh pair of per-entry output directories (one for our
// tool, one for upstream) under the cache dir and returns their absolute paths.
func mkOutDirs(cacheDir string) (ourDir, upDir string, err error) {
	base, err := os.MkdirTemp(cacheDir, "out-")
	if err != nil {
		return "", "", err
	}
	ourDir = filepath.Join(base, "ours")
	upDir = filepath.Join(base, "up")
	if err := os.MkdirAll(ourDir, 0o755); err != nil {
		return "", "", err
	}
	if err := os.MkdirAll(upDir, 0o755); err != nil {
		return "", "", err
	}
	return ourDir, upDir, nil
}

// timedRun runs a binary with args, returning stdout, stderr, wall-clock, and a
// run error (non-zero exit). A binary path ending in ".pl" is invoked through
// `perl` (prinseq ships as a Perl script, not a compiled binary).
func timedRun(bin string, args []string) (stdout, stderr []byte, dur time.Duration, err error) {
	var cmd *exec.Cmd
	if strings.HasSuffix(bin, ".pl") {
		cmd = exec.Command("perl", append([]string{bin}, args...)...)
	} else {
		cmd = exec.Command(bin, args...)
	}
	var out, errb bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = &errb
	start := time.Now()
	err = cmd.Run()
	dur = time.Since(start)
	return out.Bytes(), errb.Bytes(), dur, err
}

// placeholderKeys are the manifest-backed fixture tokens resolvePlaceholders
// substitutes. {out} is handled separately because it is per-invocation.
var placeholderKeys = []string{
	"bam", "cram", "vcf", "vcf_plain", "vcf_multi", "bed", "bed12", "fasta", "genome",
	"fastq", "fastq_gz", "fastq1", "fastq2", "gff",
}

// resolvePlaceholders substitutes {bam}, {cram}, {vcf}, {fastq}, {gff}, ...
// tokens in args with the corresponding fixture path from the manifest, and the
// {out} token with outPrefix (the per-invocation output prefix; "" when the
// entry does not use one).
func resolvePlaceholders(args []string, m *fixtures.Manifest, outPrefix string) ([]string, error) {
	out := make([]string, len(args))
	for i, a := range args {
		if !strings.Contains(a, "{") {
			out[i] = a
			continue
		}
		r := a
		if strings.Contains(r, "{out}") {
			if outPrefix == "" {
				return nil, fmt.Errorf("arg %q uses {out} but the entry declares no OutputFiles", a)
			}
			r = strings.ReplaceAll(r, "{out}", outPrefix)
		}
		for _, key := range placeholderKeys {
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
