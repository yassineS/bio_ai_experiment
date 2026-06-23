package main

import (
	"bytes"
	"fmt"
	"os"
	"path/filepath"
	"strings"
	"time"

	"github.com/yassineS/bio_ai_experiment/pipeline/runner"
)

// binset holds the resolved binary paths for one tool on each side.
type binset struct {
	oursSamtools, upSamtools string
	oursBcftools, upBcftools string
}

func (b binset) ours(tool string) string {
	if tool == "bcftools" {
		return b.oursBcftools
	}
	return b.oursSamtools
}

func (b binset) upstream(tool string) string {
	if tool == "bcftools" {
		return b.upBcftools
	}
	return b.upSamtools
}

// inputs holds the resolved input file paths (empty when not provided).
type inputs struct {
	ref, bam, vcf string
}

// config is the fully-resolved run configuration.
type config struct {
	bins    binset
	in      inputs
	region  string
	reps    int
	outDir  string
	tmpDir  string
	verbose bool
}

// runBattery executes every cell in the battery, returning the assembled
// report. It never returns an error for an individual cell failure (those
// become ERROR/DIVERGE/SKIP rows); it only returns an error for a setup
// problem that prevents producing a report at all.
func runBattery(cfg config) (*report, error) {
	rep := &report{
		Generated: time.Now().UTC(),
		Region:    cfg.region,
		Reps:      cfg.reps,
		OurBin:    binDir(cfg.bins.oursSamtools, cfg.bins.oursBcftools),
		UpBin:     binDir(cfg.bins.upSamtools, cfg.bins.upBcftools),
	}
	rep.Inputs = describeInputs(cfg)

	for _, spec := range battery() {
		rep.Cells = append(rep.Cells, runCell(cfg, spec))
	}
	rep.finalize()
	return rep, nil
}

// runCell runs one cell: SKIP if its required input/binaries are absent,
// otherwise time and compare both sides. The returned cellResult is fully
// populated (status, measurements, ratios via finalize, diff snippet).
func runCell(cfg config, spec cellSpec) cellResult {
	res := cellResult{Name: spec.Name, Tool: spec.Tool, Multi: spec.Multi}

	ourBin := cfg.bins.ours(spec.Tool)
	upBin := cfg.bins.upstream(spec.Tool)
	if ourBin == "" || upBin == "" {
		res.Status = "SKIP"
		res.Detail = fmt.Sprintf("missing %s binary (ours=%q upstream=%q)", spec.Tool, ourBin, upBin)
		return res
	}
	if reason, ok := inputAvailable(cfg.in, spec.Need); !ok {
		res.Status = "SKIP"
		res.Detail = reason
		return res
	}

	ourSum, ourHead, ourMeas, ourErr := runSide(cfg, spec, ourBin)
	upSum, upHead, upMeas, upErr := runSide(cfg, spec, upBin)

	if ourErr != nil || upErr != nil {
		// A quickcheck-style postNone cell encodes pass/fail as exit status, so a
		// non-nil error there is meaningful data, compared below. Any other
		// command failing is a hard cell error.
		if spec.Post != postNone {
			res.Status = "ERROR"
			res.Detail = fmt.Sprintf("execution failed: ours_err=%v upstream_err=%v", ourErr, upErr)
			if ourMeas != nil {
				res.Ours = ptrSide(measToSide(*ourMeas))
			}
			if upMeas != nil {
				res.Upstream = ptrSide(measToSide(*upMeas))
			}
			return res
		}
	}

	if ourMeas != nil {
		res.Ours = ptrSide(measToSide(*ourMeas))
	}
	if upMeas != nil {
		res.Upstream = ptrSide(measToSide(*upMeas))
	}

	// For postNone (quickcheck) parity is "same exit verdict".
	if spec.Post == postNone {
		ourOK := ourErr == nil
		upOK := upErr == nil
		if ourOK == upOK {
			res.Status = "PASS"
			res.Detail = fmt.Sprintf("both %s", verdictWord(ourOK))
		} else {
			res.Status = "DIVERGE"
			res.Detail = fmt.Sprintf("exit verdict differs: ours_ok=%v upstream_ok=%v", ourOK, upOK)
		}
		return res
	}

	// Parity is exactly the provenance-stripped byte equality of CompareByteExact,
	// computed streaming: equal iff the two md5 digests of the provenance-stripped
	// streams match. The heads (first 64 KiB of stripped output) drive the diff
	// snippet on mismatch.
	if ourSum == upSum {
		res.Status = "PASS"
		return res
	}
	res.Status = "DIVERGE"
	res.Detail, res.DiffSnippet = headDiff(ourHead, upHead, ourSum, upSum)
	return res
}

// runSide builds the argv for one binary, runs it reps times streaming stdout
// through a provenance-filtering digester, and returns the digest of the
// comparable stream (stdout, or the re-decoded SAM for file-producing cells),
// the first 64 KiB of that stripped stream (head, for diff snippets), the
// reduced measurement, and any execution error. No rep ever holds the full
// output in memory.
func runSide(cfg config, spec cellSpec, bin string) (sum [16]byte, head []byte, meas *Measurement, err error) {
	args, outPath, cleanup, berr := buildArgs(cfg, spec)
	if cleanup != nil {
		defer cleanup()
	}
	if berr != nil {
		return sum, nil, nil, berr
	}

	var env []string
	if spec.Need == needBAMRef && spec.WriteOut == ".cram" {
		env = cramEnv(os.Environ())
	}

	switch spec.Post {
	case postViewSAM:
		// The command itself writes a file; its stdout is not the comparison
		// stream, so the timed reps discard stdout. The comparison stream is the
		// re-decoded SAM produced afterward (not timed).
		m, rerr := repeatRun(cfg.reps, bin, args, "", env, nil)
		if rerr != nil {
			return sum, nil, &m, rerr
		}
		// Re-decode the written BAM/CRAM through the SAME binary's `view -h` so the
		// two outputs are compared by decoded records, not BGZF/CRAM framing (the
		// repo-documented caveat), STREAMING the SAM through the digester so the
		// multi-GB decode never buffers. The decode itself is not timed.
		s, h, derr := decodeAlignment(cfg, bin, outPath, spec)
		if derr != nil {
			return sum, nil, &m, fmt.Errorf("re-decoding %s output: %w", spec.Name, derr)
		}
		return s, h, &m, nil
	default: // postStdout, postNone
		dig := runner.NewStreamDigester()
		m, rerr := repeatRun(cfg.reps, bin, args, "", env, dig)
		if rerr != nil {
			return sum, nil, &m, rerr
		}
		if cerr := dig.Close(); cerr != nil {
			return sum, nil, &m, cerr
		}
		return dig.Sum(), dig.Head(), &m, nil
	}
}

// decodeAlignment runs `bin view -h <file>` (with -T ref for CRAM), streaming
// its stdout through a provenance-filtering digester so two alignment files are
// compared by their decoded SAM records WITHOUT buffering the (potentially
// multi-GB) decoded SAM. It returns the digest and the 64 KiB head.
func decodeAlignment(cfg config, bin, file string, spec cellSpec) ([16]byte, []byte, error) {
	args := []string{"view", "-h"}
	if spec.WriteOut == ".cram" {
		args = append(args, "-T", cfg.in.ref)
	}
	args = append(args, file)
	var env []string
	if spec.WriteOut == ".cram" {
		env = cramEnv(os.Environ())
	}
	dig := runner.NewStreamDigester()
	_, err := runOnce(bin, args, "", env, dig)
	if err != nil {
		return [16]byte{}, dig.Head(), err
	}
	if cerr := dig.Close(); cerr != nil {
		return [16]byte{}, dig.Head(), cerr
	}
	return dig.Sum(), dig.Head(), nil
}

// buildArgs substitutes placeholders into the cell's argv, allocating a temp
// output path for file-producing cells and appending -region where the command
// accepts it. It returns the argv, the {out} path (empty if none), and a
// cleanup func for any temp file.
func buildArgs(cfg config, spec cellSpec) (args []string, outPath string, cleanup func(), err error) {
	repl := strings.NewReplacer(
		"{bam}", cfg.in.bam,
		"{vcf}", cfg.in.vcf,
		"{ref}", cfg.in.ref,
	)
	if spec.WriteOut != "" {
		f, ferr := os.CreateTemp(cfg.tmpDir, "realparity-*"+spec.WriteOut)
		if ferr != nil {
			return nil, "", nil, ferr
		}
		outPath = f.Name()
		f.Close()
		cleanup = func() { os.Remove(outPath) }
	}
	for _, a := range spec.Args {
		a = repl.Replace(a)
		a = strings.ReplaceAll(a, "{out}", outPath)
		args = append(args, a)
	}
	// Apply -region (when requested) using the mechanism each command actually
	// accepts. samtools/bcftools do NOT take a region uniformly: `samtools view`
	// and `bcftools view` take a trailing POSITIONAL region, but `samtools depth`,
	// `bcftools query` and `bcftools stats` take it via `-r REGION` (a bare
	// trailing "20" is parsed as a second input FILE and fails — and `samtools
	// view -r` would mean read-group, not region). norm/sort/flagstat/idxstats/
	// stats(sam)/quickcheck take no region at all.
	if cfg.region != "" {
		switch regionStyle(spec) {
		case regionPositional:
			args = append(args, cfg.region)
		case regionRFlag:
			args = append(args, "-r", cfg.region)
		}
	}
	return args, outPath, cleanup, nil
}

// regionMode names how a cell accepts an operator-supplied -region.
type regionMode int

const (
	regionNone       regionMode = iota // command takes no region
	regionPositional                   // trailing positional region (samtools/bcftools view)
	regionRFlag                        // -r REGION (samtools depth, bcftools query/stats)
)

// regionStyle reports how this cell's argv accepts a region, so buildArgs adds
// it in the form both binaries actually parse.
func regionStyle(spec cellSpec) regionMode {
	switch spec.Name {
	case "samtools_view_sam", "samtools_view_sam_header",
		"bcftools_view", "bcftools_view_body":
		return regionPositional
	case "samtools_depth_a", "bcftools_stats", "bcftools_query":
		return regionRFlag
	}
	return regionNone
}

// inputAvailable reports whether the required inputs for a cell are present.
func inputAvailable(in inputs, need inputKind) (string, bool) {
	switch need {
	case needBAM:
		if in.bam == "" {
			return "no -bam provided", false
		}
	case needBAMRef:
		if in.bam == "" {
			return "no -bam provided", false
		}
		if in.ref == "" {
			return "no -ref provided (CRAM path)", false
		}
	case needVCF:
		if in.vcf == "" {
			return "no -vcf provided", false
		}
	case needVCFRef:
		if in.vcf == "" {
			return "no -vcf provided", false
		}
		if in.ref == "" {
			return "no -ref provided (norm -f)", false
		}
	}
	return "", true
}

func ptrSide(s sideMeas) *sideMeas { return &s }

func verdictWord(ok bool) string {
	if ok {
		return "OK"
	}
	return "FAIL"
}

// headDiff builds the DIVERGE Detail + DiffSnippet from the two 64 KiB heads of
// provenance-stripped output. When the first difference lies WITHIN the head
// window, it reuses snippetFromStripped to render the labelled first-diff
// excerpt. When the heads are byte-identical but the full-stream digests differ,
// the divergence is beyond the captured window, so it reports that fact and the
// two digests (the heads are already provenance-stripped, so an in-window diff
// would have been found).
func headDiff(ourHead, upHead []byte, ourSum, upSum [16]byte) (detail, snippet string) {
	if !bytes.Equal(ourHead, upHead) {
		snippet = snippetFromStripped(ourHead, upHead)
		return "outputs differ (first diff within first 64KiB)", snippet
	}
	return fmt.Sprintf("differ beyond first 64KiB; digests %x vs %x", ourSum, upSum), ""
}

// snippetFromStripped renders a short first-diff excerpt from two ALREADY
// provenance-stripped byte windows (the digester heads). It does NOT re-strip;
// the bytes it receives are exactly the normalized stream parity compared.
func snippetFromStripped(ours, upstream []byte) string {
	a := strings.Split(string(ours), "\n")
	b := strings.Split(string(upstream), "\n")
	n := len(a)
	if len(b) < n {
		n = len(b)
	}
	for i := 0; i < n; i++ {
		if a[i] != b[i] {
			lo := i - 1
			if lo < 0 {
				lo = 0
			}
			var sb strings.Builder
			for j := lo; j <= i && j < n; j++ {
				fmt.Fprintf(&sb, "  ours[%d]: %s\n", j+1, trunc(a[j]))
				fmt.Fprintf(&sb, "  upst[%d]: %s\n", j+1, trunc(b[j]))
			}
			return strings.TrimRight(sb.String(), "\n")
		}
	}
	if len(a) != len(b) {
		return fmt.Sprintf("line count differs (within head): ours=%d upstream=%d", len(a), len(b))
	}
	return "streams differ (no line-level diff located within head)"
}

// diffSnippet returns a short provenance-stripped first-diff excerpt for the
// report from two RAW (un-stripped) outputs. It reuses runner.StripProvenance so
// the embedded snippet shows the SAME normalized bytes the parity comparison
// saw. It is retained for tests and any caller holding full outputs; the live
// streaming path uses snippetFromStripped on the already-stripped heads.
func diffSnippet(ours, upstream []byte) string {
	return snippetFromStripped(runner.StripProvenance(ours), runner.StripProvenance(upstream))
}

func trunc(s string) string {
	if len(s) > 160 {
		return s[:160] + "…"
	}
	return s
}

// binDir returns the common parent directory of the two resolved binaries, for
// the report header. Falls back to the samtools binary's dir.
func binDir(sam, bcf string) string {
	if sam != "" {
		return filepath.Dir(sam)
	}
	if bcf != "" {
		return filepath.Dir(bcf)
	}
	return ""
}
