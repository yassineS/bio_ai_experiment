package realbench

import (
	"compress/gzip"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strings"
	"time"

	"github.com/yassineS/bio_ai_experiment/pipeline/runner"
)

// Config is the fully-resolved run configuration for one tier.
type Config struct {
	Tier     string
	Inputs   Inputs
	Reps     int
	OutDir   string
	TmpDir   string
	Verbose  bool
	resolver *BinResolver
}

// ResolveInputs builds an Inputs from the raw CLI paths, absolutising each.
func ResolveInputs(ref, bam, cram, vcf, fastq1, fastq2, bed, gff string) Inputs {
	return Inputs{
		Ref:    abs(ref),
		BAM:    abs(bam),
		CRAM:   abs(cram),
		VCF:    abs(vcf),
		Fastq1: abs(fastq1),
		Fastq2: abs(fastq2),
		BED:    abs(bed),
		GFF:    abs(gff),
	}
}

// WithResolver returns cfg with the binary resolver attached (the resolver is an
// unexported field so the CLI sets it through this helper).
func WithResolver(cfg Config, r *BinResolver) Config {
	cfg.resolver = r
	return cfg
}

// Run executes the whole matrix for the configured tier and returns the report.
// It never returns an error for an individual cell (those become
// PASS/DIFF/SKIP/ERROR rows); it only errors on a setup problem that prevents
// producing a report at all.
func Run(cfg Config) (*Report, error) {
	rep := &Report{
		Tier:      cfg.Tier,
		Generated: time.Now().UTC(),
		Reps:      cfg.Reps,
		Machine:   machineInfo(),
	}
	for _, spec := range Matrix(cfg.Tier) {
		rep.Cells = append(rep.Cells, runCell(cfg, spec))
	}
	rep.finalize()
	return rep, nil
}

// runCell runs one cell: SKIP when its required input or binary is absent,
// otherwise time and compare both sides. The returned CellRecord is fully
// populated (parity, measurements, note); finalize() derives the ratios.
func runCell(cfg Config, spec CellSpec) CellRecord {
	res := CellRecord{
		Tool:       spec.Tool,
		Name:       spec.Name,
		Subcommand: spec.Subcommand,
		Tier:       cfg.Tier,
	}

	if miss := cfg.Inputs.missing(spec.Need); miss != "" {
		res.Parity = "SKIP"
		res.Note = "missing required input " + miss
		return res
	}

	ourBin := cfg.resolver.ourBinary(spec.Tool)
	if ourBin == "" {
		res.Parity = "SKIP"
		res.Note = "our " + spec.Tool + " binary not resolved"
		return res
	}

	// Ours always runs (perf cell at minimum).
	ourSum, ourMeas, ourErr := runSide(cfg, spec, ourBin, nil, false)
	if ourMeas != nil {
		res.Ours = measToSide(*ourMeas)
	}

	// Ours-only cells (no upstream pair) are perf cells: record SKIP parity.
	if spec.Post == PostOursOnly {
		res.Parity = "SKIP"
		if ourErr != nil {
			res.Parity = "ERROR"
			res.Note = "ours-only run failed: " + ourErr.Error()
		} else {
			res.Note = "ours-only perf cell (no upstream pair)"
		}
		return res
	}

	up := cfg.resolver.upstreamBinary(spec.Tool)
	if up.Path == "" {
		// No upstream binary available — degrade to an ours-only perf cell.
		res.Parity = "SKIP"
		note := up.NoteWhen
		if note == "" {
			note = "upstream " + spec.Tool + " not resolved"
		}
		res.Note = note + " (ran ours only)"
		if ourErr != nil {
			res.Parity = "ERROR"
			res.Note = "ours run failed: " + ourErr.Error()
		}
		return res
	}

	upSum, upMeas, upErr := runSide(cfg, spec, up.Path, up.UpStub, up.IsPerl)
	if upMeas != nil {
		res.Up = measToSide(*upMeas)
	}

	// PostNone (quickcheck): parity is "same exit verdict".
	if spec.Post == PostNone {
		ourOK := ourErr == nil
		upOK := upErr == nil
		if ourOK == upOK {
			res.Parity = "PASS"
			res.Note = fmt.Sprintf("both %s", verdictWord(ourOK))
		} else {
			res.Parity = "DIFF"
			res.Note = fmt.Sprintf("exit verdict differs: ours_ok=%v upstream_ok=%v", ourOK, upOK)
		}
		return res
	}

	if ourErr != nil || upErr != nil {
		res.Parity = "ERROR"
		res.Note = fmt.Sprintf("execution failed: ours_err=%v upstream_err=%v", ourErr, upErr)
		return res
	}

	if ourSum == upSum {
		res.Parity = "PASS"
		return res
	}
	res.Parity = "DIFF"
	res.Note = fmt.Sprintf("provenance-stripped outputs differ (digests %x vs %x)", ourSum, upSum)
	return res
}

// runSide builds the argv for one binary, runs it reps times streaming stdout
// through a provenance-filtering digester, and returns the digest of the
// comparable stream, the reduced measurement, and any execution error. No rep
// ever holds the full output in memory.
//
// upStub is the leading upstream argv (e.g. the bedtools subcommand) prepended
// before the cell's args; it is nil/empty for our binary and for upstream tools
// that take no sub. perl=true runs the binary through `perl` (prinseq).
func runSide(cfg Config, spec CellSpec, bin string, upStub []string, perl bool) (sum [16]byte, meas *Measurement, err error) {
	args, runBin, outPath, workDir, cleanup, berr := buildSide(cfg, spec, bin, upStub, perl)
	if cleanup != nil {
		defer cleanup()
	}
	if berr != nil {
		return sum, nil, berr
	}

	var env []string
	if spec.NeedsCRAMEnv {
		env = cramEnv(os.Environ())
	}

	switch spec.Post {
	case PostStdout:
		dig := runner.NewStreamDigester()
		m, rerr := repeatRun(cfg.Reps, runBin, args, "", workDir, env, dig)
		if rerr != nil {
			return sum, &m, rerr
		}
		if cerr := dig.Close(); cerr != nil {
			return sum, &m, cerr
		}
		return dig.Sum(), &m, nil

	case PostStdoutGzip:
		// The command's stdout is a gzip/BGZF stream. The timed reps discard it;
		// one extra, untimed rep streams stdout through a gzip reader into the
		// digester so the decompressed payload is compared (framing-independent).
		m, rerr := repeatRun(cfg.Reps, runBin, args, "", workDir, env, nil)
		if rerr != nil {
			return sum, &m, rerr
		}
		s, derr := digestGzipStdout(runBin, args, workDir, env)
		if derr != nil {
			return sum, &m, derr
		}
		return s, &m, nil

	case PostNone, PostOursOnly:
		// Exit status is the signal; output discarded. PostNone (e.g. quickcheck)
		// compares the two sides' exit verdicts; PostOursOnly has no upstream pair
		// to compare against, so it is a pure perf/exit measurement. Either way we
		// only need to run the command reps times and surface any execution error.
		m, rerr := repeatRun(cfg.Reps, runBin, args, "", workDir, env, nil)
		return sum, &m, rerr

	case PostViewSAM, PostBgzipD, PostFile:
		// The command writes a file; stdout is not the comparison stream, so the
		// timed reps discard it. The comparison stream is the re-decoded file
		// (not timed).
		m, rerr := repeatRun(cfg.Reps, runBin, args, "", workDir, env, nil)
		if rerr != nil {
			return sum, &m, rerr
		}
		cmpPath := outPath
		if spec.WorkDirOut {
			cmpPath = filepath.Join(workDir, spec.Compare)
		}
		s, derr := decodeOutput(cfg, spec, bin, cmpPath)
		if derr != nil {
			return sum, &m, fmt.Errorf("post-processing %s output: %w", spec.Name, derr)
		}
		return s, &m, nil
	}
	return sum, nil, fmt.Errorf("unhandled post kind for %s", spec.Name)
}

// buildSide assembles the argv (with placeholders substituted), the binary to
// exec (handling the `perl <script>` case), the {out} path, the per-side work
// dir, and a cleanup func. For file-producing cells a unique temp output path
// (or work dir) is allocated per side so ours and upstream never clobber each
// other.
func buildSide(cfg Config, spec CellSpec, bin string, upStub []string, perl bool) (args []string, runBin, outPath, workDir string, cleanup func(), err error) {
	var cleanups []func()
	cleanup = func() {
		for _, c := range cleanups {
			c()
		}
	}

	if spec.WorkDirOut {
		d, derr := os.MkdirTemp(cfg.TmpDir, "rb-wd-*")
		if derr != nil {
			return nil, "", "", "", cleanup, derr
		}
		workDir = d
		cleanups = append(cleanups, func() { os.RemoveAll(d) })
	} else if spec.WriteOut != "" {
		f, ferr := os.CreateTemp(cfg.TmpDir, "rb-out-*"+spec.WriteOut)
		if ferr != nil {
			return nil, "", "", "", cleanup, ferr
		}
		outPath = f.Name()
		f.Close()
		// bgzip refuses to overwrite an existing file by default; remove the
		// placeholder so the command can create it cleanly, and clean up after.
		os.Remove(outPath)
		cleanups = append(cleanups, func() { os.Remove(outPath) })
	}

	// Choose the cell's argv for this side.
	cellArgs := spec.OurArgs
	if len(upStub) > 0 && spec.UpArgs != nil {
		cellArgs = spec.UpArgs
	}
	subst := substituteArgs(cellArgs, cfg.Inputs, outPath, workDir)
	args = append(append([]string{}, upStub...), subst...)

	runBin = bin
	if perl {
		args = append([]string{bin}, args...)
		runBin = "perl"
	}
	return args, runBin, outPath, workDir, cleanup, nil
}

// decodeOutput turns a file-producing cell's written output into the comparable
// digest, per PostKind. For PostViewSAM it re-decodes through `<bin> view -h`
// (so a BAM/CRAM/BCF is compared by decoded records, framing-independent — bin
// is the side's own samtools/bcftools); for PostBgzipD it gzip-decompresses;
// for PostFile it reads the bytes directly. The decode streams through a
// provenance-filtering digester so a multi-GB decode never buffers.
func decodeOutput(cfg Config, spec CellSpec, bin, file string) ([16]byte, error) {
	switch spec.Post {
	case PostViewSAM:
		return decodeAlignment(cfg, spec, bin, file)
	case PostBgzipD:
		return decodeGzip(file)
	case PostFile:
		return digestFile(file)
	}
	return [16]byte{}, fmt.Errorf("decodeOutput: unhandled post kind for %s", spec.Name)
}

// decodeAlignment runs `<bin> view -h <file>` (with -T ref for CRAM), streaming
// its stdout through a provenance-filtering digester so two alignment files are
// compared by their decoded records WITHOUT buffering the (potentially multi-GB)
// decoded output. bin is the SAME binary that wrote the file (our samtools/
// bcftools for our side, upstream for the other), so `view -h` decodes BAM,
// CRAM and BCF alike.
func decodeAlignment(cfg Config, spec CellSpec, bin, file string) ([16]byte, error) {
	args := []string{"view", "-h"}
	var env []string
	if strings.HasSuffix(file, ".cram") || spec.NeedsCRAMEnv {
		args = append(args, "-T", cfg.Inputs.Ref)
		env = cramEnv(os.Environ())
	}
	args = append(args, file)
	dig := runner.NewStreamDigester()
	if _, err := runOnce(bin, args, "", "", env, dig); err != nil {
		return [16]byte{}, err
	}
	if cerr := dig.Close(); cerr != nil {
		return [16]byte{}, cerr
	}
	return dig.Sum(), nil
}

// digestGzipStdout runs the command once (untimed), streaming its gzip/BGZF
// stdout through a multistream gzip reader into a provenance-filtering digester
// via an io.Pipe, so the decompressed payload is digested without buffering the
// (potentially multi-GB) stream.
func digestGzipStdout(bin string, args []string, workDir string, env []string) ([16]byte, error) {
	pr, pw := io.Pipe()
	type result struct {
		sum [16]byte
		err error
	}
	done := make(chan result, 1)
	go func() {
		gr, err := gzip.NewReader(pr)
		if err != nil {
			pr.CloseWithError(err)
			done <- result{err: err}
			return
		}
		gr.Multistream(true)
		sum, derr := digestReader(gr)
		gr.Close()
		// Drain any remaining bytes so the writer never blocks on a full pipe.
		io.Copy(io.Discard, pr)
		done <- result{sum: sum, err: derr}
	}()
	_, runErr := runOnce(bin, args, "", workDir, env, pw)
	pw.Close()
	res := <-done
	if runErr != nil {
		return [16]byte{}, runErr
	}
	return res.sum, res.err
}

// decodeGzip gzip-decompresses file (BGZF is gzip-compatible, multistream) and
// digests its provenance-stripped payload.
func decodeGzip(file string) ([16]byte, error) {
	f, err := os.Open(file)
	if err != nil {
		return [16]byte{}, err
	}
	defer f.Close()
	gr, err := gzip.NewReader(f)
	if err != nil {
		return [16]byte{}, err
	}
	defer gr.Close()
	gr.Multistream(true)
	return digestReader(gr)
}

// digestFile digests a plain file's provenance-stripped bytes.
func digestFile(file string) ([16]byte, error) {
	f, err := os.Open(file)
	if err != nil {
		return [16]byte{}, err
	}
	defer f.Close()
	return digestReader(f)
}

// digestReader streams r through a provenance-filtering digester.
func digestReader(r io.Reader) ([16]byte, error) {
	dig := runner.NewStreamDigester()
	if _, err := io.Copy(dig, r); err != nil {
		return [16]byte{}, err
	}
	if err := dig.Close(); err != nil {
		return [16]byte{}, err
	}
	return dig.Sum(), nil
}

func verdictWord(ok bool) string {
	if ok {
		return "OK"
	}
	return "FAIL"
}

// cramEnv returns the child environment used for CRAM operations so a missing
// MD5 in the local cache never triggers a network fetch. Pointing REF_PATH at a
// non-existent path plus passing -T <ref> on the command line makes CRAM fully
// local and offline (mirrors realparity).
func cramEnv(base []string) []string {
	out := make([]string, 0, len(base)+2)
	for _, e := range base {
		if strings.HasPrefix(e, "REF_PATH=") || strings.HasPrefix(e, "REF_CACHE=") {
			continue
		}
		out = append(out, e)
	}
	return append(out, "REF_PATH=/dev/null", "REF_CACHE=")
}
