package bcftools

// Live-binary parity tests for the per-subcommand feature gaps closed in
// the "bcftools subfeatures" change:
//
//   - query format tokens: bare %INFO (whole INFO column) and %SAMPLE
//   - roh -O z (BGZF-compressed output)
//   - filter -M/--mask-file (and --mask) soft-filter
//   - cnv --AF-file (per-site allele frequency + targets filter)
//
// Each test runs the genuine upstream bcftools binary vendored under
// reference_code/bcftools AND our Go port on the same fixture, then
// asserts byte-equality of the relevant output. Unlike the live_oracle
// suite these tests never t.Skip: when the upstream binary is missing
// they build it (htslib + bcftools) via the uniquely-named
// upstreamBcftoolsSubfeatures builder, and any failure is t.Fatalf.

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"sync"
	"testing"
)

// upstreamBcftoolsSubfeatures returns the path to the vendored upstream
// bcftools binary for the subfeatures parity suite. It delegates to the
// shared upstreamBcftools builder (upstream_test.go), which builds htslib
// (recursive submodule init for htscodecs) and bcftools under a
// cross-process build lock exactly once and t.Fatalf's on failure — never
// t.Skip. The distinct name keeps this suite's intent self-documenting
// while reusing the robust build path.
func upstreamBcftoolsSubfeatures(t *testing.T) string {
	t.Helper()
	return upstreamBcftools(t)
}

// subfeatOurBin builds the local bcftools port binary once and returns
// its path. It reuses the TestMain-built binary (ourBinPath) when
// available and otherwise builds into t.TempDir() so the artefact is
// cleaned up with the test.
var (
	subfeatOurOnce sync.Once
	subfeatOurPath string
	subfeatOurErr  error
)

func subfeatOursBin(t *testing.T) string {
	t.Helper()
	if ourBinPath != "" {
		return ourBinPath
	}
	subfeatOurOnce.Do(func() {
		bin := filepath.Join(t.TempDir(), "bcftools")
		cmd := exec.Command("go", "build", "-o", bin, "../../cmd/bcftools")
		if out, err := cmd.CombinedOutput(); err != nil {
			subfeatOurErr = err
			t.Logf("go build output: %s", out)
			return
		}
		subfeatOurPath = bin
	})
	if subfeatOurPath == "" {
		t.Fatalf("failed to build local bcftools port: %v", subfeatOurErr)
	}
	return subfeatOurPath
}

// runCapture runs a binary and returns stdout, failing on a non-zero
// exit (stderr is captured for the failure message).
func runCapture(t *testing.T, bin string, args ...string) []byte {
	t.Helper()
	cmd := exec.Command(bin, args...)
	var out, errb bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = &errb
	if err := cmd.Run(); err != nil {
		t.Fatalf("%s %v: %v\nstderr: %s", bin, args, err, errb.String())
	}
	return out.Bytes()
}

// writeTempFile writes content to a uniquely-named file in t.TempDir()
// and returns its path.
func writeTempFile(t *testing.T, dir, name, content string) string {
	t.Helper()
	p := filepath.Join(dir, name)
	if err := os.WriteFile(p, []byte(content), 0o644); err != nil {
		t.Fatalf("write %s: %v", p, err)
	}
	return p
}

const subfeatQueryVCF = `##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="dp">
##INFO=<ID=AA,Number=1,Type=String,Description="aa">
##INFO=<ID=FLAG1,Number=0,Type=Flag,Description="flag">
##FORMAT=<ID=GT,Number=1,Type=String,Description="gt">
##FORMAT=<ID=CH,Number=1,Type=Character,Description="char">
##contig=<ID=1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
1	100	.	A	T,C	50	PASS	DP=5;AA=x;FLAG1	GT:CH	0/1:a	1/2:b
1	200	.	G	A	30	PASS	DP=8	GT:CH	0/0:c	1/1:d
`

// TestSubfeatQueryTokensParity checks the bare %INFO and %SAMPLE query
// format tokens byte-for-byte against upstream.
func TestSubfeatQueryTokensParity(t *testing.T) {
	up := upstreamBcftoolsSubfeatures(t)
	ours := subfeatOursBin(t)
	dir := t.TempDir()
	vcf := writeTempFile(t, dir, "q.vcf", subfeatQueryVCF)

	formats := []string{
		`%CHROM\t%INFO\n`,
		`[%SAMPLE=%GT\n]`,
		`%CHROM\t%INFO\t[%SAMPLE:%CH ]\n`,
		`[%SAMPLE\t%CH\n]`,
	}
	for _, f := range formats {
		f := f
		t.Run(f, func(t *testing.T) {
			want := runCapture(t, up, "query", "-f", f, vcf)
			got := runCapture(t, ours, "query", "-f", f, vcf)
			if !bytes.Equal(want, got) {
				t.Fatalf("query -f %q mismatch:\nupstream=%q\nours    =%q", f, want, got)
			}
		})
	}
}

const subfeatRohVCF = `##fileformat=VCFv4.2
##INFO=<ID=AF,Number=A,Type=Float,Description="af">
##FORMAT=<ID=GT,Number=1,Type=String,Description="gt">
##contig=<ID=1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
1	100	.	A	T	50	PASS	AF=0.2	GT	0/0
1	200	.	G	C	50	PASS	AF=0.3	GT	0/0
1	300	.	A	T	50	PASS	AF=0.2	GT	1/1
1	400	.	G	C	50	PASS	AF=0.3	GT	1/1
1	500	.	A	T	50	PASS	AF=0.2	GT	0/1
`

// TestSubfeatRohOzParity checks roh -O z BGZF output decodes to the same
// data rows as upstream. (The provenance comment lines that upstream
// writes into the output file are stripped, as elsewhere in the oracle.)
func TestSubfeatRohOzParity(t *testing.T) {
	up := upstreamBcftoolsSubfeatures(t)
	ours := subfeatOursBin(t)
	dir := t.TempDir()
	vcf := writeTempFile(t, dir, "roh.vcf", subfeatRohVCF)

	// Cover both -Osrz (s,r,z explicit) and bare -Oz (compression only):
	// upstream defaults the sections to s+r when only z is given, so bare
	// -Oz must still emit the full report compressed.
	for _, ot := range []string{"-Osrz", "-Oz"} {
		ot := ot
		t.Run(ot, func(t *testing.T) {
			upOut := filepath.Join(dir, "up"+ot+".roh.gz")
			ourOut := filepath.Join(dir, "our"+ot+".roh.gz")
			runCapture(t, up, "roh", "-G30", "--AF-tag", "AF", ot, "-o", upOut, vcf)
			runCapture(t, ours, "roh", "-G30", "--AF-tag", "AF", ot, "-o", ourOut, vcf)

			upData := dataRows(t, gunzipFile(t, upOut))
			ourData := dataRows(t, gunzipFile(t, ourOut))
			if !bytes.Equal(upData, ourData) {
				t.Fatalf("roh %s data mismatch:\nupstream=%q\nours    =%q", ot, upData, ourData)
			}
			if len(ourData) == 0 {
				t.Fatalf("roh %s produced no data rows", ot)
			}
		})
	}
}

// dataRows keeps only ST/RG data lines (dropping the '#'-prefixed header
// and provenance lines that legitimately differ).
func dataRows(t *testing.T, b []byte) []byte {
	t.Helper()
	var out bytes.Buffer
	for _, line := range strings.Split(string(b), "\n") {
		if strings.HasPrefix(line, "ST\t") || strings.HasPrefix(line, "RG\t") {
			out.WriteString(line)
			out.WriteByte('\n')
		}
	}
	return out.Bytes()
}

const subfeatFilterVCF = `##fileformat=VCFv4.2
##contig=<ID=1>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	100	.	A	T	50	.	.
1	200	.	G	C	50	.	.
1	300	.	A	T	50	.	.
1	400	.	C	G	50	.	.
`

// TestSubfeatFilterMaskParity checks filter --mask-file (BED + tab),
// --mask region strings, negation, and --mask-overlap against upstream.
func TestSubfeatFilterMaskParity(t *testing.T) {
	up := upstreamBcftoolsSubfeatures(t)
	ours := subfeatOursBin(t)
	dir := t.TempDir()
	vcf := writeTempFile(t, dir, "f.vcf", subfeatFilterVCF)
	bed := writeTempFile(t, dir, "mask.bed", "1\t150\t250\n1\t380\t420\n")
	tab := writeTempFile(t, dir, "mask.txt", "1\t200\t300\n")

	cases := [][]string{
		{"filter", "-M", bed, "-s", "MASK", "--no-version"},
		{"filter", "-M", tab, "-s", "MASK", "--no-version"},
		{"filter", "-M", "^" + bed, "-s", "MASK", "--no-version"},
		{"filter", "--mask", "1:150-250", "-s", "M2", "--no-version"},
		{"filter", "-M", bed, "--mask-overlap", "0", "-s", "M", "-m", "x", "--no-version"},
	}
	for _, base := range cases {
		base := base
		t.Run(strings.Join(base[1:], " "), func(t *testing.T) {
			want := stripHeaders(runCapture(t, up, append(append([]string{}, base...), vcf)...))
			got := stripHeaders(runCapture(t, ours, append(append([]string{}, base...), vcf)...))
			if !bytes.Equal(want, got) {
				t.Fatalf("filter mask mismatch:\nupstream=%q\nours    =%q", want, got)
			}
		})
	}
}

// stripHeaders drops every '##'-prefixed header line, leaving the
// #CHROM line and data records for comparison.
func stripHeaders(b []byte) []byte {
	var out bytes.Buffer
	for _, line := range strings.Split(string(b), "\n") {
		if strings.HasPrefix(line, "##") {
			continue
		}
		out.WriteString(line)
		out.WriteByte('\n')
	}
	return out.Bytes()
}

// TestSubfeatCNVAFFileParity checks cnv --AF-file: the per-site genotype
// frequencies and the targets-filter behaviour. The AF file must be
// BGZF-compressed and tabix-indexed for upstream, so we build that with
// the vendored bgzip/tabix.
func TestSubfeatCNVAFFileParity(t *testing.T) {
	up := upstreamBcftoolsSubfeatures(t)
	ours := subfeatOursBin(t)
	dir := t.TempDir()

	vcf, afTab := buildCNVFixture(t, dir)
	afGz := bgzipTabix(t, afTab, dir)

	upDir := filepath.Join(dir, "up")
	ourDir := filepath.Join(dir, "our")
	if err := os.MkdirAll(upDir, 0o755); err != nil {
		t.Fatal(err)
	}
	if err := os.MkdirAll(ourDir, 0o755); err != nil {
		t.Fatal(err)
	}

	runCapture(t, up, "cnv", "-s", "S1", "--AF-file", afGz, "-o", upDir, vcf)
	ourOut := runCapture(t, ours, "cnv", "-s", "S1", "--AF-file", afGz, "-o", ourDir, vcf)

	upRG := rgRows(t, mustRead(t, filepath.Join(upDir, "summary.S1.tab")))
	ourRG := rgRows(t, ourOut)
	if len(upRG) == 0 {
		t.Fatalf("upstream cnv produced no RG rows")
	}
	// Upstream encodes the copy number as a bare integer (2) in column 5;
	// our port spells it "CN2". Normalise that single column before
	// comparing the rest (chrom, start, end, quality, nSites, nHETs).
	if !equalRGRows(upRG, ourRG) {
		t.Fatalf("cnv --AF-file RG mismatch:\nupstream=%q\nours    =%q", upRG, ourRG)
	}
}

// rgRows extracts the RG data rows from a cnv summary, splitting each on
// tabs.
func rgRows(t *testing.T, b []byte) [][]string {
	t.Helper()
	var out [][]string
	for _, line := range strings.Split(string(b), "\n") {
		if strings.HasPrefix(line, "RG\t") {
			out = append(out, strings.Split(line, "\t"))
		}
	}
	return out
}

// equalRGRows compares two RG row sets, mapping upstream's integer copy
// number (column index 4) onto our "CNn" spelling.
func equalRGRows(up, ours [][]string) bool {
	if len(up) != len(ours) {
		return false
	}
	for i := range up {
		a, b := up[i], ours[i]
		if len(a) != len(b) {
			return false
		}
		for j := range a {
			av, bv := a[j], b[j]
			if j == 4 { // copy-number column
				av = "CN" + av
			}
			if av != bv {
				return false
			}
		}
	}
	return true
}

func mustRead(t *testing.T, path string) []byte {
	t.Helper()
	b, err := os.ReadFile(path)
	if err != nil {
		t.Fatalf("read %s: %v", path, err)
	}
	return b
}

// gunzipFile reads a BGZF/gzip file and returns its decompressed bytes.
// Go's gzip.Reader transparently reads BGZF's concatenated blocks.
func gunzipFile(t *testing.T, path string) []byte {
	t.Helper()
	return gunzipBytes(t, mustRead(t, path))
}

// buildCNVFixture writes a synthetic single-sample CNV VCF (with FORMAT
// BAF/LRR) and a matching tab AF file into dir, returning their paths.
// The data is fixed (no RNG) so the parity comparison is deterministic.
func buildCNVFixture(t *testing.T, dir string) (vcfPath, afPath string) {
	t.Helper()
	var vcf, af strings.Builder
	vcf.WriteString("##fileformat=VCFv4.2\n")
	vcf.WriteString(`##FORMAT=<ID=GT,Number=1,Type=String,Description="gt">` + "\n")
	vcf.WriteString(`##FORMAT=<ID=BAF,Number=1,Type=Float,Description="baf">` + "\n")
	vcf.WriteString(`##FORMAT=<ID=LRR,Number=1,Type=Float,Description="lrr">` + "\n")
	vcf.WriteString("##contig=<ID=1>\n")
	vcf.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
	// Deterministic BAF pattern: cycle 0.0/0.5/1.0 with a het-heavy
	// middle block, so the HMM has structure to call.
	bafCycle := []string{"0.000", "0.502", "1.000"}
	for i := 0; i < 60; i++ {
		pos := 1000 + i*1000
		baf := bafCycle[i%3]
		lrr := "0.010"
		if i >= 20 && i < 40 {
			baf = "0.030"
			lrr = "-0.480"
		}
		vcf.WriteString(stringsJoinCNV("1", pos, baf, lrr))
		af.WriteString("1\t" + itoa(pos) + "\tA,G\t0." + itoa(100+i%300) + "\n")
	}
	return writeTempFile(t, dir, "cnv.vcf", vcf.String()),
		writeTempFile(t, dir, "cnvaf.tab", af.String())
}

// stringsJoinCNV builds one CNV VCF data line.
func stringsJoinCNV(chrom string, pos int, baf, lrr string) string {
	return chrom + "\t" + itoa(pos) + "\t.\tA\tG\t50\tPASS\t.\tGT:BAF:LRR\t0/1:" + baf + ":" + lrr + "\n"
}

// bgzipTabix bgzip-compresses tabPath and tabix-indexes it (the upstream
// targets-file requirement), returning the .gz path. It uses the
// vendored bgzip/tabix binaries built alongside the bcftools binary.
func bgzipTabix(t *testing.T, tabPath, dir string) string {
	t.Helper()
	htsBin := filepath.Join(filepath.Dir(filepath.Dir(upstreamBcftoolsSubfeatures(t))), "htslib")
	bgzip := filepath.Join(htsBin, "bgzip")
	tabix := filepath.Join(htsBin, "tabix")
	gz := tabPath + ".gz"
	raw, err := os.ReadFile(tabPath)
	if err != nil {
		t.Fatal(err)
	}
	cmd := exec.Command(bgzip, "-c")
	cmd.Stdin = bytes.NewReader(raw)
	out, err := cmd.Output()
	if err != nil {
		t.Fatalf("bgzip: %v", err)
	}
	if err := os.WriteFile(gz, out, 0o644); err != nil {
		t.Fatal(err)
	}
	if err := exec.Command(tabix, "-s1", "-b2", "-e2", gz).Run(); err != nil {
		t.Fatalf("tabix: %v", err)
	}
	return gz
}
