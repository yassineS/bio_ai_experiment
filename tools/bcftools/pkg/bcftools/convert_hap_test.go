package bcftools

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"os"
	"os/exec"
	"path/filepath"
	"runtime"
	"strings"
	"sync"
	"testing"
)

// ---------------------------------------------------------------------------
// Live upstream parity harness.
//
// These tests shell out to the freshly built upstream `bcftools` binary
// under reference_code/bcftools and compare its HAP/legend/sample output
// (and the VCF round-trip) against this port byte-for-byte. There are no
// recorded golden files: the harness always exercises the live binary.
// ---------------------------------------------------------------------------

var (
	upstreamConvertHapPath string
	upstreamConvertHapErr  error
	upstreamConvertHapOnce sync.Once
)

// upstreamBcftoolsConvertHap returns the path to the upstream bcftools
// binary, building it (and htslib) from the reference_code submodules if
// it is not already present. The result is cached across tests.
func upstreamBcftoolsConvertHap(t *testing.T) string {
	t.Helper()
	upstreamConvertHapOnce.Do(func() {
		root, err := repoRootConvertHap()
		if err != nil {
			upstreamConvertHapErr = err
			return
		}
		bin := filepath.Join(root, "reference_code", "bcftools", "bcftools")
		if _, statErr := os.Stat(bin); statErr == nil {
			upstreamConvertHapPath = bin
			return
		}
		// Build htslib then bcftools.
		htsDir := filepath.Join(root, "reference_code", "htslib")
		bcfDir := filepath.Join(root, "reference_code", "bcftools")
		for _, step := range []struct {
			dir  string
			name string
			args []string
		}{
			{htsDir, "make", []string{"-j", numCPUConvertHap()}},
			{bcfDir, "make", []string{"-j", numCPUConvertHap()}},
		} {
			cmd := exec.Command(step.name, step.args...)
			cmd.Dir = step.dir
			if out, err := cmd.CombinedOutput(); err != nil {
				upstreamConvertHapErr = &buildErr{step.dir, err, out}
				return
			}
		}
		if _, statErr := os.Stat(bin); statErr != nil {
			upstreamConvertHapErr = statErr
			return
		}
		upstreamConvertHapPath = bin
	})
	if upstreamConvertHapErr != nil {
		t.Fatalf("could not obtain upstream bcftools: %v", upstreamConvertHapErr)
	}
	return upstreamConvertHapPath
}

type buildErr struct {
	dir string
	err error
	out []byte
}

func (b *buildErr) Error() string {
	return "build failed in " + b.dir + ": " + b.err.Error() + "\n" + string(b.out)
}

func numCPUConvertHap() string {
	n := runtime.NumCPU()
	if n < 1 {
		n = 1
	}
	return itoaConvertHap(n)
}

func itoaConvertHap(n int) string {
	if n == 0 {
		return "0"
	}
	var b [20]byte
	i := len(b)
	for n > 0 {
		i--
		b[i] = byte('0' + n%10)
		n /= 10
	}
	return string(b[i:])
}

func repoRootConvertHap() (string, error) {
	dir, err := os.Getwd()
	if err != nil {
		return "", err
	}
	for {
		if _, err := os.Stat(filepath.Join(dir, "go.mod")); err == nil {
			return dir, nil
		}
		parent := filepath.Dir(dir)
		if parent == dir {
			return "", os.ErrNotExist
		}
		dir = parent
	}
}

// readMaybeGz reads a file, transparently gunzipping it when it has the
// gzip magic bytes (covers both BGZF and plain gzip outputs).
func readMaybeGz(t *testing.T, path string) []byte {
	t.Helper()
	raw, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("read %s: %v", path, err)
	}
	if len(raw) >= 2 && raw[0] == 0x1f && raw[1] == 0x8b {
		gr, err := gzip.NewReader(bytes.NewReader(raw))
		if err != nil {
			t.Fatalf("gunzip %s: %v", path, err)
		}
		gr.Multistream(true)
		var out bytes.Buffer
		if _, err := out.ReadFrom(gr); err != nil {
			t.Fatalf("gunzip %s: %v", path, err)
		}
		return out.Bytes()
	}
	return raw
}

// hapSampleVCF is a small biallelic mixed-phase corpus shared by parity tests.
const hapSampleVCF = `##fileformat=VCFv4.2
##contig=<ID=20>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
20	100	.	A	G	.	.	.	GT	0|0	0|1	1|1
20	200	rs5	C	T	.	.	.	GT	1|0	0/1	.|.
20	300	.	G	A	.	.	.	GT	0/0	1/1	0|1
`

// hapHaploidVCF exercises the haploid-allele branches of the hap encoder.
const hapHaploidVCF = `##fileformat=VCFv4.2
##contig=<ID=X>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	M1	M2
X	500	.	A	T	.	.	.	GT	0	1
X	600	.	C	G	.	.	.	GT	.	0
`

func writeTempVCF(t *testing.T, content string) string {
	t.Helper()
	dir := t.TempDir()
	p := filepath.Join(dir, "in.vcf")
	if err := os.WriteFile(p, []byte(content), 0o644); err != nil {
		t.Fatalf("write vcf: %v", err)
	}
	return p
}

func TestVCFToHapSampleParity(t *testing.T) {
	bin := upstreamBcftoolsConvertHap(t)
	in := writeTempVCF(t, hapSampleVCF)

	for _, hap2dip := range []bool{false, true} {
		hap2dip := hap2dip
		name := "plain"
		if hap2dip {
			name = "hap2dip"
		}
		t.Run(name, func(t *testing.T) {
			dir := t.TempDir()
			upPrefix := filepath.Join(dir, "up")
			goPrefix := filepath.Join(dir, "go")

			upArgs := []string{"convert"}
			if hap2dip {
				upArgs = append(upArgs, "--haploid2diploid")
			}
			upArgs = append(upArgs, "--hapsample", upPrefix, in)
			runUpstream(t, bin, upArgs...)

			opts := HapConvertOptions{Prefix: goPrefix, Hap2Dip: hap2dip}
			if _, err := VCFToHapSample(in, opts, devNull{}); err != nil {
				t.Fatalf("VCFToHapSample: %v", err)
			}

			assertFileEqual(t, "hap", readMaybeGz(t, upPrefix+".hap.gz"), readMaybeGz(t, goPrefix+".hap.gz"))
			assertFileEqual(t, "samples", readMaybeGz(t, upPrefix+".samples"), readMaybeGz(t, goPrefix+".samples"))
		})
	}
}

func TestVCFToHapLegendSampleParity(t *testing.T) {
	bin := upstreamBcftoolsConvertHap(t)
	in := writeTempVCF(t, hapSampleVCF)
	dir := t.TempDir()
	upPrefix := filepath.Join(dir, "up")
	goPrefix := filepath.Join(dir, "go")

	runUpstream(t, bin, "convert", "--haplegendsample", upPrefix, in)
	if _, err := VCFToHapLegendSample(in, HapConvertOptions{Prefix: goPrefix}, devNull{}); err != nil {
		t.Fatalf("VCFToHapLegendSample: %v", err)
	}

	assertFileEqual(t, "hap", readMaybeGz(t, upPrefix+".hap.gz"), readMaybeGz(t, goPrefix+".hap.gz"))
	assertFileEqual(t, "legend", readMaybeGz(t, upPrefix+".legend.gz"), readMaybeGz(t, goPrefix+".legend.gz"))
	assertFileEqual(t, "samples", readMaybeGz(t, upPrefix+".samples"), readMaybeGz(t, goPrefix+".samples"))
}

func TestVCFToHapSampleHaploidParity(t *testing.T) {
	bin := upstreamBcftoolsConvertHap(t)
	in := writeTempVCF(t, hapHaploidVCF)
	dir := t.TempDir()
	upPrefix := filepath.Join(dir, "up")
	goPrefix := filepath.Join(dir, "go")

	runUpstream(t, bin, "convert", "--hapsample", upPrefix, in)
	if _, err := VCFToHapSample(in, HapConvertOptions{Prefix: goPrefix}, devNull{}); err != nil {
		t.Fatalf("VCFToHapSample: %v", err)
	}
	assertFileEqual(t, "hap", readMaybeGz(t, upPrefix+".hap.gz"), readMaybeGz(t, goPrefix+".hap.gz"))
}

func TestHapSampleToVCFParity(t *testing.T) {
	bin := upstreamBcftoolsConvertHap(t)
	in := writeTempVCF(t, hapSampleVCF)
	dir := t.TempDir()
	prefix := filepath.Join(dir, "h")

	// Produce hap/sample with upstream, then round-trip both ways back to VCF.
	runUpstream(t, bin, "convert", "--hapsample", prefix, in)

	upVCF := runUpstreamStdout(t, bin, "convert", "--hapsample2vcf", prefix, "-O", "v")

	var goBuf bytes.Buffer
	if _, err := HapSampleToVCF(prefix, &goBuf, HapConvertOptions{OutputFormat: OutputVCF}, devNull{}); err != nil {
		t.Fatalf("HapSampleToVCF: %v", err)
	}
	assertVCFBodyEqual(t, upVCF, goBuf.Bytes())
}

func TestHapLegendSampleToVCFParity(t *testing.T) {
	bin := upstreamBcftoolsConvertHap(t)
	in := writeTempVCF(t, hapSampleVCF)
	dir := t.TempDir()
	prefix := filepath.Join(dir, "h")

	runUpstream(t, bin, "convert", "--haplegendsample", prefix, in)
	upVCF := runUpstreamStdout(t, bin, "convert", "--haplegendsample2vcf", prefix, "-O", "v")

	var goBuf bytes.Buffer
	if _, err := HapLegendSampleToVCF(prefix, &goBuf, HapConvertOptions{OutputFormat: OutputVCF}, devNull{}); err != nil {
		t.Fatalf("HapLegendSampleToVCF: %v", err)
	}
	assertVCFBodyEqual(t, upVCF, goBuf.Bytes())
}

// TestHapRoundTripParity verifies VCF -> hap -> VCF is loss-free relative
// to upstream for the body (CHROM..GT) of the records.
func TestHapRoundTripParity(t *testing.T) {
	bin := upstreamBcftoolsConvertHap(t)
	in := writeTempVCF(t, hapSampleVCF)
	dir := t.TempDir()
	prefix := filepath.Join(dir, "rt")

	if _, err := VCFToHapLegendSample(in, HapConvertOptions{Prefix: prefix}, devNull{}); err != nil {
		t.Fatalf("VCFToHapLegendSample: %v", err)
	}
	upVCF := runUpstreamStdout(t, bin, "convert", "--haplegendsample2vcf", prefix, "-O", "v")

	var goBuf bytes.Buffer
	if _, err := HapLegendSampleToVCF(prefix, &goBuf, HapConvertOptions{OutputFormat: OutputVCF}, devNull{}); err != nil {
		t.Fatalf("HapLegendSampleToVCF: %v", err)
	}
	assertVCFBodyEqual(t, upVCF, goBuf.Bytes())
}

// TestHapLegendSampleToVCFReversedAllelesParity exercises the rev_als path
// where the legend's a0/a1 columns are swapped relative to the
// CHROM:POS_REF_ALT id, so the haplotype 0/1 labels must be flipped.
func TestHapLegendSampleToVCFReversedAllelesParity(t *testing.T) {
	bin := upstreamBcftoolsConvertHap(t)
	dir := t.TempDir()
	prefix := filepath.Join(dir, "rev")

	writeGzFile(t, prefix+".legend.gz", "id position a0 a1\n20:100_A_G 100 G A\n")
	writeGzFile(t, prefix+".hap.gz", "0 1 1 0\n")
	if err := os.WriteFile(prefix+".samples", []byte("sample population group sex\ns1 s1 s1 2\ns2 s2 s2 2\n"), 0o644); err != nil {
		t.Fatalf("write samples: %v", err)
	}

	upVCF := runUpstreamStdout(t, bin, "convert", "--haplegendsample2vcf", prefix, "-O", "v")
	var goBuf bytes.Buffer
	if _, err := HapLegendSampleToVCF(prefix, &goBuf, HapConvertOptions{OutputFormat: OutputVCF}, devNull{}); err != nil {
		t.Fatalf("HapLegendSampleToVCF: %v", err)
	}
	assertVCFBodyEqual(t, upVCF, goBuf.Bytes())
}

func writeGzFile(t *testing.T, path, content string) {
	t.Helper()
	f, err := os.Create(path)
	if err != nil {
		t.Fatalf("create %s: %v", path, err)
	}
	gw := gzip.NewWriter(f)
	if _, err := gw.Write([]byte(content)); err != nil {
		t.Fatalf("write %s: %v", path, err)
	}
	if err := gw.Close(); err != nil {
		t.Fatalf("close gz %s: %v", path, err)
	}
	if err := f.Close(); err != nil {
		t.Fatalf("close %s: %v", path, err)
	}
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

type devNull struct{}

func (devNull) Write(p []byte) (int, error) { return len(p), nil }

func runUpstream(t *testing.T, bin string, args ...string) {
	t.Helper()
	cmd := exec.Command(bin, args...)
	if out, err := cmd.CombinedOutput(); err != nil {
		t.Fatalf("upstream %v failed: %v\n%s", args, err, out)
	}
}

func runUpstreamStdout(t *testing.T, bin string, args ...string) []byte {
	t.Helper()
	cmd := exec.Command(bin, args...)
	var stdout, stderr bytes.Buffer
	cmd.Stdout = &stdout
	cmd.Stderr = &stderr
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream %v failed: %v\n%s", args, err, stderr.String())
	}
	return stdout.Bytes()
}

func assertFileEqual(t *testing.T, label string, want, got []byte) {
	t.Helper()
	if !bytes.Equal(want, got) {
		t.Fatalf("%s mismatch:\n--- upstream ---\n%s\n--- port ---\n%s", label, want, got)
	}
}

// assertVCFBodyEqual compares two VCFs ignoring the provenance header
// lines that this port intentionally omits (##fileformat is kept, but the
// ##bcftools_convert* lines are upstream-only).
func assertVCFBodyEqual(t *testing.T, want, got []byte) {
	t.Helper()
	w := stripProvenance(want)
	g := stripProvenance(got)
	if w != g {
		t.Fatalf("VCF body mismatch:\n--- upstream ---\n%s\n--- port ---\n%s", w, g)
	}
}

func stripProvenance(b []byte) string {
	var sb strings.Builder
	sc := bufio.NewScanner(bytes.NewReader(b))
	sc.Buffer(make([]byte, 0, 64<<10), 16<<20)
	for sc.Scan() {
		line := sc.Text()
		if strings.HasPrefix(line, "##bcftools_") {
			continue
		}
		sb.WriteString(line)
		sb.WriteByte('\n')
	}
	return sb.String()
}

// ---------------------------------------------------------------------------
// Unit tests (no upstream binary)
// ---------------------------------------------------------------------------

func TestGtToHap(t *testing.T) {
	cases := []struct {
		gt      string
		hap2dip bool
		want    string
	}{
		{"0|0", false, "0 0"},
		{"0|1", false, "0 1"},
		{"1|0", false, "1 0"},
		{"1|1", false, "1 1"},
		{"0/0", false, "0* 0*"},
		{"0/1", false, "0* 1*"},
		{"1/1", false, "1* 1*"},
		{".|.", false, "? ?"},
		{"0|.", false, "? ?"},
		{".", false, "-1 -"},
		{"0", false, "0 -"},
		{"1", false, "1 -"},
		// hap2dip: haploid becomes diploid homozygote.
		{"0", true, "0 0"},
		{"1", true, "1 1"},
		{".", true, "? ?"},
		{"0|0", true, "0 0"},
		{"0/1", true, "0* 1*"},
	}
	for _, c := range cases {
		got := gtToHap(c.gt, c.hap2dip)
		if got != c.want {
			t.Errorf("gtToHap(%q, %v) = %q, want %q", c.gt, c.hap2dip, got, c.want)
		}
	}
}

func TestHapOutputNames(t *testing.T) {
	hap, sample, err := hapOutputNames("pfx")
	if err != nil || hap != "pfx.hap.gz" || sample != "pfx.samples" {
		t.Fatalf("hapOutputNames(prefix): %q %q %v", hap, sample, err)
	}
	hap, sample, err = hapOutputNames("a.hap.gz,b.samples")
	if err != nil || hap != "a.hap.gz" || sample != "b.samples" {
		t.Fatalf("hapOutputNames(list): %q %q %v", hap, sample, err)
	}
	hap, sample, err = hapOutputNames("a.hap.gz,.")
	if err != nil || hap != "a.hap.gz" || sample != "" {
		t.Fatalf("hapOutputNames(skip): %q %q %v", hap, sample, err)
	}
	if _, _, err := hapOutputNames("a,b,c"); err == nil {
		t.Fatalf("hapOutputNames: expected error for 3-element list")
	}
}

func TestHapLegendOutputNames(t *testing.T) {
	h, l, s, err := hapLegendOutputNames("pfx")
	if err != nil || h != "pfx.hap.gz" || l != "pfx.legend.gz" || s != "pfx.samples" {
		t.Fatalf("hapLegendOutputNames(prefix): %q %q %q %v", h, l, s, err)
	}
	h, l, s, err = hapLegendOutputNames("x,y,z")
	if err != nil || h != "x" || l != "y" || s != "z" {
		t.Fatalf("hapLegendOutputNames(list): %q %q %q %v", h, l, s, err)
	}
	if _, _, _, err := hapLegendOutputNames("a,b"); err == nil {
		t.Fatalf("hapLegendOutputNames: expected error for 2-element list")
	}
}

func TestParseHapChromPosRefAlt(t *testing.T) {
	chrom, pos, ref, alt, end, err := parseHapChromPosRefAlt("20:100_A_G")
	if err != nil || chrom != "20" || pos != 100 || ref != "A" || alt != "G" || end != 0 {
		t.Fatalf("parse: %q %d %q %q %d %v", chrom, pos, ref, alt, end, err)
	}
	chrom, pos, ref, alt, end, err = parseHapChromPosRefAlt("X:5_AT_A_42")
	if err != nil || chrom != "X" || pos != 5 || ref != "AT" || alt != "A" || end != 42 {
		t.Fatalf("parse END: %q %d %q %q %d %v", chrom, pos, ref, alt, end, err)
	}
	if _, _, _, _, _, err := parseHapChromPosRefAlt("nogthere"); err == nil {
		t.Fatalf("expected error for malformed token")
	}
}

func TestHapInputNames(t *testing.T) {
	h, s := hapInputNames("pfx")
	if h != "pfx.hap.gz" || s != "pfx.samples" {
		t.Fatalf("hapInputNames(prefix): %q %q", h, s)
	}
	h, s = hapInputNames("foo.hap.gz,foo.samples")
	if h != "foo.hap.gz" || s != "foo.samples" {
		t.Fatalf("hapInputNames(list): %q %q", h, s)
	}
}

func TestFillHapGenotypes(t *testing.T) {
	v := newHapVariant("20", 100, "A", "G", 0)
	if err := fillHapGenotypes(v, []string{"0", "1", "0*", "0*", "?", "?"}, []string{"a", "b", "c"}, false); err != nil {
		t.Fatalf("fillHapGenotypes: %v", err)
	}
	want := []string{"0|1", "0/0", ".|."}
	for i, w := range want {
		if got := v.Samples[i].Data["GT"]; got != w {
			t.Errorf("sample %d GT = %q, want %q", i, got, w)
		}
	}

	// Reversed alleles swap 0<->1.
	v = newHapVariant("20", 100, "A", "G", 0)
	if err := fillHapGenotypes(v, []string{"0", "1"}, []string{"a"}, true); err != nil {
		t.Fatalf("fillHapGenotypes rev: %v", err)
	}
	if got := v.Samples[0].Data["GT"]; got != "1|0" {
		t.Errorf("reversed GT = %q, want 1|0", got)
	}
}

func TestRefAltOrientation(t *testing.T) {
	if rev, err := refAltOrientation("A", "G", "A", "G"); err != nil || rev {
		t.Fatalf("matched orientation: rev=%v err=%v", rev, err)
	}
	if rev, err := refAltOrientation("A", "G", "G", "A"); err != nil || !rev {
		t.Fatalf("reversed orientation: rev=%v err=%v", rev, err)
	}
	if _, err := refAltOrientation("A", "G", "A", "T"); err == nil {
		t.Fatalf("expected mismatch error")
	}
}

func TestHapModeRejectsMultiple(t *testing.T) {
	// A wrong field count should be reported, not silently truncated.
	v := newHapVariant("20", 100, "A", "G", 0)
	if err := fillHapGenotypes(v, []string{"0", "1", "0"}, []string{"a", "b"}, false); err == nil {
		t.Fatalf("expected error for wrong field count")
	}
}
