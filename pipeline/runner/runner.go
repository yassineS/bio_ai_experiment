// Package runner executes the parity matrix: for each Entry it resolves our
// freshly built tool binary and the vendored upstream binary, runs both on the
// generated fixtures, compares their output per the entry's compare mode, and
// records timing. Results aggregate into a Report.
package runner

import (
	"bytes"
	"compress/gzip"
	"fmt"
	"io"
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
	// PIPELINE_NO_SKIP=1 force-runs entries that carry a documented Skip, so a
	// maintainer can re-triage whether a skip is still warranted (e.g. after a
	// tool fix lands) without editing the matrix. Off by default.
	if e.Skip != "" && os.Getenv("PIPELINE_NO_SKIP") == "" {
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

	// CopyToOut: stage fixture copies into each side's output dir so tools that
	// write alongside their input (e.g. bgzip -r → <file>.gzi) operate on an
	// isolated per-side copy.
	if len(e.CopyToOut) > 0 {
		for tok, suffix := range e.CopyToOut {
			src := cfg.Manifest.Path(tok)
			if src == "" {
				res.Status = StatusError
				res.Detail = fmt.Sprintf("CopyToOut: fixture %q not in manifest", tok)
				return res
			}
			for _, dir := range []string{ourDir, upDir} {
				if err := copyFile(src, filepath.Join(dir, "out"+suffix)); err != nil {
					res.Status = StatusError
					res.Detail = fmt.Sprintf("CopyToOut %s: %v", tok, err)
					return res
				}
			}
		}
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

	// Decode binary stdout before comparing, for tools that emit BGZF/BAM whose
	// framing is not byte-comparable but whose decoded content must match.
	switch e.CompareModeOrDefault() {
	case matrix.BGZFDecoded:
		o, oerr := gunzipAll(ourOut)
		u, uerr := gunzipAll(upOut)
		if oerr != nil || uerr != nil {
			res.Status = StatusError
			res.Detail = fmt.Sprintf("BGZF decode failed: ours=%v upstream=%v", oerr, uerr)
			return res
		}
		ourOut, upOut = o, u
	case matrix.BAMDecoded:
		// Decode BOTH sides with the upstream samtools (the canonical decoder),
		// regardless of which tool produced the BAM (e.g. bedtobam's upstream is
		// bedtools), so only the record content is compared.
		samBin, serr := upstream.Binary("samtools")
		if serr != nil {
			res.Status = StatusError
			res.Detail = "BAM decode needs the upstream samtools binary: " + serr.Error()
			return res
		}
		o, oerr := decodeBAM(samBin, ourOut)
		u, uerr := decodeBAM(samBin, upOut)
		if oerr != nil || uerr != nil {
			res.Status = StatusError
			res.Detail = fmt.Sprintf("BAM decode failed: ours=%v upstream=%v", oerr, uerr)
			return res
		}
		ourOut, upOut = o, u
	}

	var cmp CompareResult
	switch {
	case usesOutPrefix:
		// Compare the named output files between the two prefixes (decompressing
		// .gz, stripping provenance). The compare mode still selects byte-exact
		// vs similarity for the file contents.
		cmp = CompareOutputFiles(filepath.Join(ourDir, "out"), filepath.Join(upDir, "out"),
			e.OutputFiles, e.CompareModeOrDefault(), resolveEpsilon(e.Tolerance))
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
		cmp = CompareSimilarity(ourOut, upOut, resolveEpsilon(e.Tolerance))
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

// copyFile copies src to dst, creating dst's parent directory if needed.
func copyFile(src, dst string) error {
	if err := os.MkdirAll(filepath.Dir(dst), 0o755); err != nil {
		return err
	}
	in, err := os.Open(src)
	if err != nil {
		return err
	}
	defer in.Close()
	out, err := os.Create(dst)
	if err != nil {
		return err
	}
	if _, err := io.Copy(out, in); err != nil {
		out.Close()
		return err
	}
	return out.Close()
}

// gunzipAll decompresses a (possibly multi-member, e.g. BGZF) gzip stream in
// full. BGZF is a series of concatenated gzip members, which compress/gzip
// reads transparently via Reader.Multistream.
func gunzipAll(b []byte) ([]byte, error) {
	if len(b) == 0 {
		return b, nil
	}
	zr, err := gzip.NewReader(bytes.NewReader(b))
	if err != nil {
		return nil, err
	}
	defer zr.Close()
	return io.ReadAll(zr)
}

// decodeBAM pipes BAM bytes through `samtools view -h -` and returns the SAM
// text, so two BAMs with different BGZF framing can be compared by their
// decoded records (the @PG/@CO provenance is stripped by CompareByteExact).
func decodeBAM(samtoolsBin string, bam []byte) ([]byte, error) {
	if len(bam) == 0 {
		return bam, nil
	}
	cmd := exec.Command(samtoolsBin, "view", "-h", "-")
	cmd.Stdin = bytes.NewReader(bam)
	var out, errb bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = &errb
	if err := cmd.Run(); err != nil {
		return nil, fmt.Errorf("%v: %s", err, trunc(errb.String()))
	}
	return out.Bytes(), nil
}

// placeholderKeys are the manifest-backed fixture tokens resolvePlaceholders
// substitutes. {out} is handled separately because it is per-invocation.
var placeholderKeys = []string{
	"bam", "bam_namesorted", "cram", "vcf", "vcf_plain", "vcf_multi", "vcf_pl", "vcf_samples", "vcf_phased_plain", "bed", "bed12", "bedpe", "bedgraph1", "bedgraph2", "fasta", "genome",
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
