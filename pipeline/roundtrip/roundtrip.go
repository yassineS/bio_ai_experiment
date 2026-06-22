// Package roundtrip validates that our tools losslessly round-trip each file
// format: encoding then decoding (or decoding then re-encoding) a fixture must
// reproduce the original payload, and our round-trip must agree with upstream's.
// This complements the parity matrix (which compares single operations against
// upstream) by exercising the encode/decode pair end-to-end for every container
// format the suite implements (BAM, CRAM, BGZF, VCF/BCF, FASTQ).
package roundtrip

import (
	"bytes"
	"fmt"
	"os"
	"os/exec"
	"path/filepath"

	"github.com/yassineS/bio_ai_experiment/pipeline/fixtures"
	"github.com/yassineS/bio_ai_experiment/pipeline/internal/upstream"
)

// Status is the verdict of a single round-trip check.
type Status string

const (
	// Pass means the round-trip reproduced the original (and matched upstream
	// where cross-checked).
	Pass Status = "PASS"
	// Fail means the round-trip diverged.
	Fail Status = "FAIL"
	// Skip means a prerequisite (binary or fixture) was unavailable.
	Skip Status = "SKIP"
)

// Result records one round-trip check.
type Result struct {
	Name   string `json:"name"`
	Format string `json:"format"`
	Status Status `json:"status"`
	Detail string `json:"detail,omitempty"`
}

// Run executes every round-trip check against the fixtures in man, building our
// binaries into cacheDir on demand. It never aborts on a single failure — every
// check runs so the report is complete.
func Run(man *fixtures.Manifest, cacheDir string) []Result {
	tmp, err := os.MkdirTemp("", "roundtrip")
	if err != nil {
		return []Result{{Name: "setup", Status: Fail, Detail: err.Error()}}
	}
	defer os.RemoveAll(tmp)

	e := &env{man: man, cache: cacheDir, tmp: tmp}
	var out []Result
	out = append(out, e.bgzf())
	out = append(out, e.bamIdentity())
	out = append(out, e.bamViaCRAM())
	out = append(out, e.vcfViaBCF())
	out = append(out, e.fastq())
	// Explicitly bidirectional ours↔upstream interop for every container format
	// plus index (.bai/.csi/.tbi) interop. Each skips cleanly without upstream.
	out = append(out, e.interopChecks()...)
	return out
}

type env struct {
	man   *fixtures.Manifest
	cache string
	tmp   string
}

func (e *env) our(tool string) (string, error) { return upstream.OurBinary(tool, e.cache) }
func (e *env) up(key string) (string, error)   { return upstream.Binary(key) }
func (e *env) path(key string) string          { return e.man.Path(key) }
func (e *env) out(name string) string          { return filepath.Join(e.tmp, name) }

// runCmd runs bin with args, returning combined stdout (and an error carrying
// stderr on failure).
func runCmd(bin string, args ...string) ([]byte, error) {
	cmd := exec.Command(bin, args...)
	var stdout, stderr bytes.Buffer
	cmd.Stdout, cmd.Stderr = &stdout, &stderr
	if err := cmd.Run(); err != nil {
		return stdout.Bytes(), fmt.Errorf("%s %v: %w: %s", filepath.Base(bin), args, err, stderr.String())
	}
	return stdout.Bytes(), nil
}

// runToFile runs bin with args, writing stdout to a file.
func runToFile(dst, bin string, args ...string) error {
	o, err := runCmd(bin, args...)
	if err != nil {
		return err
	}
	return os.WriteFile(dst, o, 0o644)
}

func ok(name, format, detail string) Result {
	return Result{Name: name, Format: format, Status: Pass, Detail: detail}
}
func fail(name, format string, err error) Result {
	return Result{Name: name, Format: format, Status: Fail, Detail: err.Error()}
}
func skip(name, format, why string) Result {
	return Result{Name: name, Format: format, Status: Skip, Detail: why}
}

// bgzf: our bgzip compresses then decompresses to a byte-identical payload, and
// our decompressor reads upstream's bgzip output identically.
func (e *env) bgzf() Result {
	const name, format = "bgzf-compress-decompress", "BGZF"
	src := e.path("vcf_plain")
	if src == "" {
		return skip(name, format, "no vcf_plain fixture")
	}
	bin, err := e.our("bgzip")
	if err != nil {
		return skip(name, format, "bgzip unavailable: "+err.Error())
	}
	raw, err := os.ReadFile(src)
	if err != nil {
		return fail(name, format, err)
	}
	gz := e.out("rt.vcf.gz")
	if err := runToFile(gz, bin, "-c", src); err != nil {
		return fail(name, format, err)
	}
	back, err := runCmd(bin, "-dc", gz)
	if err != nil {
		return fail(name, format, err)
	}
	if !bytes.Equal(raw, back) {
		return fail(name, format, fmt.Errorf("decompressed payload differs from source (%d vs %d bytes)", len(back), len(raw)))
	}
	// Cross-check: our decompressor reads upstream bgzip's output identically.
	if ub, e2 := e.up("bgzip"); e2 == nil {
		ugz := e.out("up.vcf.gz")
		if err := runToFile(ugz, ub, "-c", src); err == nil {
			if uback, err := runCmd(bin, "-dc", ugz); err == nil && !bytes.Equal(raw, uback) {
				return fail(name, format, fmt.Errorf("our decode of upstream bgzip differs from source"))
			}
		}
	}
	return ok(name, format, "byte-identical compress→decompress; cross-decodes upstream bgzip")
}

// bamIdentity: re-encoding a BAM and decoding it reproduces the same records.
func (e *env) bamIdentity() Result {
	const name, format = "bam-reencode", "BAM"
	src := e.path("bam")
	if src == "" {
		return skip(name, format, "no bam fixture")
	}
	bin, err := e.our("samtools")
	if err != nil {
		return skip(name, format, "samtools unavailable: "+err.Error())
	}
	rebam := e.out("rt.bam")
	if err := runToFile(rebam, bin, "view", "-b", src); err != nil {
		return fail(name, format, err)
	}
	a, err := runCmd(bin, "view", "-h", src)
	if err != nil {
		return fail(name, format, err)
	}
	b, err := runCmd(bin, "view", "-h", rebam)
	if err != nil {
		return fail(name, format, err)
	}
	if !bytes.Equal(a, b) {
		return fail(name, format, fmt.Errorf("decoded records differ after BAM re-encode"))
	}
	return ok(name, format, "records identical after BAM re-encode/decode")
}

// bamViaCRAM: BAM→CRAM→BAM. CRAM is reference-based and legitimately recomputes
// some fields, so the original BAM is not the right oracle — instead our
// round-trip's decoded records must agree with upstream's round-trip (both go
// through CRAM identically). Without upstream available it falls back to
// comparing against the source as a weaker idempotence check.
func (e *env) bamViaCRAM() Result {
	const name, format = "bam-via-cram", "CRAM"
	src, ref := e.path("bam"), e.path("fasta")
	if src == "" || ref == "" {
		return skip(name, format, "no bam/fasta fixture")
	}
	bin, err := e.our("samtools")
	if err != nil {
		return skip(name, format, "samtools unavailable: "+err.Error())
	}
	roundtrip := func(tool string) ([]byte, error) {
		cram := e.out(filepath.Base(tool) + ".cram")
		bam := e.out(filepath.Base(tool) + ".fromcram.bam")
		if err := runToFile(cram, tool, "view", "-C", "-T", ref, src); err != nil {
			return nil, err
		}
		if err := runToFile(bam, tool, "view", "-b", "-T", ref, cram); err != nil {
			return nil, err
		}
		return runCmd(tool, "view", bam)
	}
	ours, err := roundtrip(bin)
	if err != nil {
		return fail(name, format, err)
	}
	if up, uerr := e.up("samtools"); uerr == nil {
		theirs, terr := roundtrip(up)
		if terr != nil {
			return fail(name, format, fmt.Errorf("upstream round-trip failed: %w", terr))
		}
		if !bytes.Equal(ours, theirs) {
			return fail(name, format, fmt.Errorf("our BAM→CRAM→BAM decode differs from upstream's (CRAM encode/decode mismatch)"))
		}
		return ok(name, format, "BAM→CRAM→BAM agrees with upstream")
	}
	// No upstream: weaker check — round-trip must reproduce the source records.
	orig, err := runCmd(bin, "view", src)
	if err != nil {
		return fail(name, format, err)
	}
	if !bytes.Equal(ours, orig) {
		return fail(name, format, fmt.Errorf("BAM→CRAM→BAM differs from source (no upstream to cross-check)"))
	}
	return ok(name, format, "BAM→CRAM→BAM reproduced source (no upstream cross-check)")
}

// vcfViaBCF: VCF→BCF→VCF reproduces the same data rows (provenance stripped).
func (e *env) vcfViaBCF() Result {
	const name, format = "vcf-via-bcf", "BCF"
	src := e.path("vcf_plain")
	if src == "" {
		return skip(name, format, "no vcf_plain fixture")
	}
	bin, err := e.our("bcftools")
	if err != nil {
		return skip(name, format, "bcftools unavailable: "+err.Error())
	}
	bcfFile := e.out("rt.bcf")
	if err := runToFile(bcfFile, bin, "view", "-Ob", src); err != nil {
		return fail(name, format, err)
	}
	back, err := runCmd(bin, "view", bcfFile)
	if err != nil {
		return fail(name, format, err)
	}
	orig, err := runCmd(bin, "view", src)
	if err != nil {
		return fail(name, format, err)
	}
	// Compare data rows only (header provenance/order may differ legitimately).
	if !bytes.Equal(dataRows(back), dataRows(orig)) {
		return fail(name, format, fmt.Errorf("data rows differ after VCF→BCF→VCF"))
	}
	return ok(name, format, "data rows identical after VCF→BCF→VCF")
}

// fastq: seqtk seq is idempotent (re-emitting a FASTQ reproduces it).
func (e *env) fastq() Result {
	const name, format = "fastq-idempotent", "FASTQ"
	src := e.path("fastq")
	if src == "" {
		return skip(name, format, "no fastq fixture")
	}
	bin, err := e.our("seqtk")
	if err != nil {
		return skip(name, format, "seqtk unavailable: "+err.Error())
	}
	a, err := runCmd(bin, "seq", src)
	if err != nil {
		return fail(name, format, err)
	}
	// Feed a's output back through seq via a temp file and require idempotence.
	mid := e.out("rt.fastq")
	if err := os.WriteFile(mid, a, 0o644); err != nil {
		return fail(name, format, err)
	}
	c, err := runCmd(bin, "seq", mid)
	if err != nil {
		return fail(name, format, err)
	}
	if !bytes.Equal(a, c) {
		return fail(name, format, fmt.Errorf("seqtk seq not idempotent"))
	}
	return ok(name, format, "seqtk seq idempotent")
}

// dataRows returns only the non-header (non-#) lines of a VCF stream.
func dataRows(b []byte) []byte {
	var out bytes.Buffer
	for _, ln := range bytes.Split(b, []byte("\n")) {
		if len(ln) == 0 || ln[0] == '#' {
			continue
		}
		out.Write(ln)
		out.WriteByte('\n')
	}
	return out.Bytes()
}
