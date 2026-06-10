package bcftools

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"strings"
	"sync"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/iohelper"
)

// --- live upstream harness ---------------------------------------------------

var (
	upstreamGenOnce sync.Once
	upstreamGenPath string
	upstreamGenErr  error
)

// upstreamBcftoolsConvertGen locates (building if necessary) the upstream
// bcftools binary under reference_code/bcftools. It returns the absolute path
// or "" when the submodule is not checked out / cannot be built. The build is
// performed at most once per test process.
func upstreamBcftoolsConvertGen(t *testing.T) string {
	t.Helper()
	upstreamGenOnce.Do(func() {
		root, err := filepath.Abs(filepath.Join("..", "..", "..", ".."))
		if err != nil {
			upstreamGenErr = err
			return
		}
		bcfDir := filepath.Join(root, "reference_code", "bcftools")
		if _, err := os.Stat(filepath.Join(bcfDir, "vcfconvert.c")); err != nil {
			// Submodule not checked out: leave path empty so callers skip.
			return
		}
		bin := filepath.Join(bcfDir, "bcftools")
		if fi, err := os.Stat(bin); err == nil && fi.Mode()&0o111 != 0 {
			upstreamGenPath = bin
			return
		}
		// Try to build (htslib must already be built).
		cmd := exec.Command("make", "-j4")
		cmd.Dir = bcfDir
		if out, err := cmd.CombinedOutput(); err != nil {
			upstreamGenErr = errWithOutput(err, out)
			return
		}
		if fi, err := os.Stat(bin); err == nil && fi.Mode()&0o111 != 0 {
			upstreamGenPath = bin
		}
	})
	if upstreamGenErr != nil {
		t.Fatalf("building upstream bcftools failed: %v", upstreamGenErr)
	}
	return upstreamGenPath
}

type outputErr struct {
	err error
	out []byte
}

func (e *outputErr) Error() string { return e.err.Error() + ": " + string(e.out) }
func errWithOutput(err error, out []byte) error {
	return &outputErr{err: err, out: out}
}

// gunzipAll fully decompresses gzip/bgzf data via the project iohelper so the
// semantic .gen content is compared rather than the (gzip vs bgzf) framing.
func gunzipAll(t *testing.T, path string) []byte {
	t.Helper()
	r, err := iohelper.OpenReader(path)
	if err != nil {
		t.Fatalf("open %s: %v", path, err)
	}
	defer r.Close()
	var buf bytes.Buffer
	if _, err := buf.ReadFrom(r); err != nil {
		t.Fatalf("read %s: %v", path, err)
	}
	return buf.Bytes()
}

const genTestVCF = `##fileformat=VCFv4.2
##contig=<ID=20,length=63025520>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
20	100	rs1	C	T	.	.	.	GT	0/0	0/1	1/1
20	200	rs2	A	G	.	.	.	GT	1/1	./.	0/0
`

const genTestVCFProbs = `##fileformat=VCFv4.2
##contig=<ID=20,length=63025520>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GP,Number=G,Type=Float,Description="GP">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="PL">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
20	100	rs1	C	T	.	.	.	GT:GP:PL	0/0:0.9,0.1,0:0,10,40	0/1:0.1,0.8,0.1:20,0,15
`

func writeTemp(t *testing.T, dir, name, content string) string {
	t.Helper()
	p := filepath.Join(dir, name)
	if err := os.WriteFile(p, []byte(content), 0o644); err != nil {
		t.Fatalf("write %s: %v", name, err)
	}
	return p
}

// --- VCF -> GEN unit tests ---------------------------------------------------

func TestVCFToGenSample_GT(t *testing.T) {
	dir := t.TempDir()
	in := writeTemp(t, dir, "in.vcf", genTestVCF)
	prefix := filepath.Join(dir, "out")
	if err := VCFToGenSampleFile(in, prefix, GenSampleOptions{}); err != nil {
		t.Fatalf("VCFToGenSampleFile: %v", err)
	}
	gen := string(gunzipAll(t, prefix+".gen.gz"))
	want := "20:100_C_T 20:100_C_T 100 C T 1 0 0 0 1 0 0 0 1\n" +
		"20:200_A_G 20:200_A_G 200 A G 0 0 1 0.33 0.33 0.33 1 0 0\n"
	if gen != want {
		t.Fatalf("gen mismatch:\n got: %q\nwant: %q", gen, want)
	}
	smpl, err := os.ReadFile(prefix + ".samples")
	if err != nil {
		t.Fatalf("read samples: %v", err)
	}
	wantS := "ID_1 ID_2 missing\n0 0 0\nS1 S1 0\nS2 S2 0\nS3 S3 0\n"
	if string(smpl) != wantS {
		t.Fatalf("samples mismatch:\n got: %q\nwant: %q", smpl, wantS)
	}
}

func TestVCFToGenSample_3N6_VCFIDs(t *testing.T) {
	dir := t.TempDir()
	in := writeTemp(t, dir, "in.vcf", genTestVCF)
	prefix := filepath.Join(dir, "out")
	if err := VCFToGenSampleFile(in, prefix, GenSampleOptions{ThreeN6: true, VCFIDs: true}); err != nil {
		t.Fatalf("VCFToGenSampleFile: %v", err)
	}
	gen := string(gunzipAll(t, prefix+".gen.gz"))
	want := "20 20:100_C_T rs1 100 C T 1 0 0 0 1 0 0 0 1\n" +
		"20 20:200_A_G rs2 200 A G 0 0 1 0.33 0.33 0.33 1 0 0\n"
	if gen != want {
		t.Fatalf("gen mismatch:\n got: %q\nwant: %q", gen, want)
	}
}

func TestVCFToGenSample_TagGP(t *testing.T) {
	dir := t.TempDir()
	in := writeTemp(t, dir, "in.vcf", genTestVCFProbs)
	prefix := filepath.Join(dir, "out")
	if err := VCFToGenSampleFile(in, prefix, GenSampleOptions{Tag: "GP"}); err != nil {
		t.Fatalf("VCFToGenSampleFile: %v", err)
	}
	gen := string(gunzipAll(t, prefix+".gen.gz"))
	want := "20:100_C_T 20:100_C_T 100 C T 0.900000 0.100000 0.000000 0.100000 0.800000 0.100000\n"
	if gen != want {
		t.Fatalf("gen mismatch:\n got: %q\nwant: %q", gen, want)
	}
}

func TestVCFToGenSample_TagPL(t *testing.T) {
	dir := t.TempDir()
	in := writeTemp(t, dir, "in.vcf", genTestVCFProbs)
	prefix := filepath.Join(dir, "out")
	if err := VCFToGenSampleFile(in, prefix, GenSampleOptions{Tag: "PL"}); err != nil {
		t.Fatalf("VCFToGenSampleFile: %v", err)
	}
	gen := string(gunzipAll(t, prefix+".gen.gz"))
	// linear probs derived from phred PLs, normalised.
	want := "20:100_C_T 20:100_C_T 100 C T 0.909008 0.090901 0.000091 0.009600 0.960040 0.030359\n"
	if gen != want {
		t.Fatalf("gen mismatch:\n got: %q\nwant: %q", gen, want)
	}
}

func TestVCFToGenSample_Sex(t *testing.T) {
	dir := t.TempDir()
	in := writeTemp(t, dir, "in.vcf", genTestVCF)
	sex := writeTemp(t, dir, "sex.txt", "S1\tM\nS2\tF\nS3\tM\n")
	prefix := filepath.Join(dir, "out")
	if err := VCFToGenSampleFile(in, prefix, GenSampleOptions{SexFile: sex}); err != nil {
		t.Fatalf("VCFToGenSampleFile: %v", err)
	}
	smpl, err := os.ReadFile(prefix + ".samples")
	if err != nil {
		t.Fatalf("read samples: %v", err)
	}
	want := "ID_1 ID_2 missing sex\n0 0 0 0\nS1 S1 0 1\nS2 S2 0 2\nS3 S3 0 1\n"
	if string(smpl) != want {
		t.Fatalf("samples mismatch:\n got: %q\nwant: %q", smpl, want)
	}
}

func TestVCFToGenSample_SexMissing(t *testing.T) {
	dir := t.TempDir()
	in := writeTemp(t, dir, "in.vcf", genTestVCF)
	sex := writeTemp(t, dir, "sex.txt", "S1\tM\nS2\tF\n") // S3 missing
	prefix := filepath.Join(dir, "out")
	err := VCFToGenSampleFile(in, prefix, GenSampleOptions{SexFile: sex})
	if err == nil || !strings.Contains(err.Error(), "missing sex for sample S3") {
		t.Fatalf("expected missing-sex error, got %v", err)
	}
}

func TestVCFToGenSample_ExplicitNames(t *testing.T) {
	dir := t.TempDir()
	in := writeTemp(t, dir, "in.vcf", genTestVCF)
	genF := filepath.Join(dir, "g.gen")
	sampleF := filepath.Join(dir, "s.sample")
	if err := VCFToGenSampleFile(in, genF+","+sampleF, GenSampleOptions{}); err != nil {
		t.Fatalf("VCFToGenSampleFile: %v", err)
	}
	// .gen without .gz is written uncompressed.
	gen, err := os.ReadFile(genF)
	if err != nil {
		t.Fatalf("read gen: %v", err)
	}
	if !strings.HasPrefix(string(gen), "20:100_C_T") {
		t.Fatalf("unexpected gen content: %q", gen)
	}
	if _, err := os.Stat(sampleF); err != nil {
		t.Fatalf("sample file not written: %v", err)
	}
}

// --- GEN -> VCF unit tests ---------------------------------------------------

func TestGenSampleRoundTrip(t *testing.T) {
	dir := t.TempDir()
	in := writeTemp(t, dir, "in.vcf", genTestVCF)
	prefix := filepath.Join(dir, "rt")
	if err := VCFToGenSampleFile(in, prefix, GenSampleOptions{VCFIDs: true}); err != nil {
		t.Fatalf("export: %v", err)
	}
	var buf bytes.Buffer
	if err := GenSampleToVCFFile(prefix, &buf, GenSampleOptions{VCFIDs: true, OutputFormat: OutputVCF}); err != nil {
		t.Fatalf("import: %v", err)
	}
	got := buf.String()
	for _, want := range []string{
		"##contig=<ID=20,length=2147483647>",
		"20\t100\trs1\tC\tT\t.\t.\t.\tGT:GP\t0/0:1,0,0\t0/1:0,1,0\t1/1:0,0,1",
		"20\t200\trs2\tA\tG\t.\t.\t.\tGT:GP\t1/1:0,0,1\t0/0:0.33,0.33,0.33\t0/0:1,0,0",
	} {
		if !strings.Contains(got, want) {
			t.Fatalf("round-trip output missing %q:\n%s", want, got)
		}
	}
}

// TestGenSampleToVCF_LabelInSecondColumn covers the IMPUTE2 reference-panel
// case where the first column is "--" and the CHROM:POS_REF_ALT label lives in
// the second column. With --vcf-ids the (empty-ish) first column becomes the
// VCF ID, mirroring upstream's two-stage tsv setter.
func TestGenSampleToVCF_LabelInSecondColumn(t *testing.T) {
	dir := t.TempDir()
	genF := writeTemp(t, dir, "ref.gen",
		"rs99 20:100_C_T 100 C T 1 0 0 0 1 0\n")
	sampleF := writeTemp(t, dir, "ref.samples",
		"ID_1 ID_2 missing\n0 0 0\nS1 S1 0\nS2 S2 0\n")

	var buf bytes.Buffer
	if err := GenSampleToVCFFile(genF+","+sampleF, &buf, GenSampleOptions{VCFIDs: true, OutputFormat: OutputVCF}); err != nil {
		t.Fatalf("import: %v", err)
	}
	want := "20\t100\trs99\tC\tT\t.\t.\t.\tGT:GP\t0/0:1,0,0\t0/1:0,1,0"
	if !strings.Contains(buf.String(), want) {
		t.Fatalf("missing %q in:\n%s", want, buf.String())
	}
}

func TestParseChromPosRefAlt(t *testing.T) {
	cases := []struct {
		in              string
		chrom, ref, alt string
		pos             int
		ok              bool
	}{
		{"20:100_C_T", "20", "C", "T", 100, true},
		{"chrX:5_A_G_120", "chrX", "A", "G", 5, true}, // _END suffix dropped
		{"no_colon", "", "", "", 0, false},
		{"20:abc_C_T", "", "", "", 0, false},
	}
	for _, c := range cases {
		chrom, pos, ref, alt, ok := parseChromPosRefAlt(c.in)
		if ok != c.ok || chrom != c.chrom || pos != c.pos || ref != c.ref || alt != c.alt {
			t.Errorf("parseChromPosRefAlt(%q) = (%q,%d,%q,%q,%v), want (%q,%d,%q,%q,%v)",
				c.in, chrom, pos, ref, alt, ok, c.chrom, c.pos, c.ref, c.alt, c.ok)
		}
	}
}

func TestGtFromProb3(t *testing.T) {
	cases := []struct {
		aa, ab, bb float64
		want       string
	}{
		{1, 0, 0, "0/0"},
		{0, 1, 0, "0/1"},
		{0, 0, 1, "1/1"},
		{0.33, 0.33, 0.33, "0/0"}, // ties favour REF hom
		{0.2, 0.3, 0.5, "1/1"},
	}
	for _, c := range cases {
		if got := gtFromProb3(c.aa, c.ab, c.bb); got != c.want {
			t.Errorf("gtFromProb3(%v,%v,%v) = %q, want %q", c.aa, c.ab, c.bb, got, c.want)
		}
	}
}

func TestNormalizeGenTag(t *testing.T) {
	for in, want := range map[string]string{"": "GT", "GT": "GT", "PL": "PL", "GP": "GP"} {
		got, err := normalizeGenTag(in)
		if err != nil || got != want {
			t.Errorf("normalizeGenTag(%q) = (%q,%v), want %q", in, got, err, want)
		}
	}
	if _, err := normalizeGenTag("XX"); err == nil {
		t.Errorf("normalizeGenTag(XX) should error")
	}
}

func TestParseGenSamplePrefix(t *testing.T) {
	g, s, err := parseGenSamplePrefix("pfx", true)
	if err != nil || g != "pfx.gen.gz" || s != "pfx.samples" {
		t.Fatalf("prefix expand: %q %q %v", g, s, err)
	}
	g, s, err = parseGenSamplePrefix("a.gen,b.sample", true)
	if err != nil || g != "a.gen" || s != "b.sample" {
		t.Fatalf("explicit: %q %q %v", g, s, err)
	}
	g, s, _ = parseGenSamplePrefix("a.gen,.", true)
	if g != "a.gen" || s != "" {
		t.Fatalf("skip-sample: %q %q", g, s)
	}
	if _, _, err := parseGenSamplePrefix("a.gen,", false); err == nil {
		t.Fatalf("import with empty sample side should error")
	}
}

// --- live parity tests -------------------------------------------------------

func TestUpstreamParity_VCFToGen(t *testing.T) {
	bin := upstreamBcftoolsConvertGen(t)
	if bin == "" {
		t.Skip("upstream bcftools submodule not checked out")
	}
	dir := t.TempDir()
	in := writeTemp(t, dir, "in.vcf", genTestVCF)
	inP := writeTemp(t, dir, "inp.vcf", genTestVCFProbs)
	sex := writeTemp(t, dir, "sex.txt", "S1\tM\nS2\tF\nS3\tM\n")

	cases := []struct {
		name string
		src  string
		opts GenSampleOptions
		args []string
	}{
		{"plain", in, GenSampleOptions{}, nil},
		{"3N6", in, GenSampleOptions{ThreeN6: true}, []string{"--3N6"}},
		{"vcfids", in, GenSampleOptions{VCFIDs: true}, []string{"--vcf-ids"}},
		{"3N6_vcfids", in, GenSampleOptions{ThreeN6: true, VCFIDs: true}, []string{"--3N6", "--vcf-ids"}},
		{"tagGP", inP, GenSampleOptions{Tag: "GP"}, []string{"--tag", "GP"}},
		{"tagPL", inP, GenSampleOptions{Tag: "PL"}, []string{"--tag", "PL"}},
		{"sex", in, GenSampleOptions{SexFile: sex}, []string{"--sex", sex}},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			// upstream
			upPfx := filepath.Join(dir, "up_"+c.name)
			args := append([]string{"convert", "-g", upPfx}, c.args...)
			args = append(args, c.src)
			cmd := exec.Command(bin, args...)
			if out, err := cmd.CombinedOutput(); err != nil {
				t.Fatalf("upstream convert -g: %v\n%s", err, out)
			}
			// go
			goPfx := filepath.Join(dir, "go_"+c.name)
			if err := VCFToGenSampleFile(c.src, goPfx, c.opts); err != nil {
				t.Fatalf("go convert -g: %v", err)
			}

			upGen := gunzipAll(t, upPfx+".gen.gz")
			goGen := gunzipAll(t, goPfx+".gen.gz")
			if !bytes.Equal(upGen, goGen) {
				t.Fatalf(".gen mismatch:\nup: %q\ngo: %q", upGen, goGen)
			}
			upS, _ := os.ReadFile(upPfx + ".samples")
			goS, _ := os.ReadFile(goPfx + ".samples")
			if !bytes.Equal(upS, goS) {
				t.Fatalf(".samples mismatch:\nup: %q\ngo: %q", upS, goS)
			}
		})
	}
}

func TestUpstreamParity_GenToVCF(t *testing.T) {
	bin := upstreamBcftoolsConvertGen(t)
	if bin == "" {
		t.Skip("upstream bcftools submodule not checked out")
	}
	dir := t.TempDir()
	in := writeTemp(t, dir, "in.vcf", genTestVCF)
	inP := writeTemp(t, dir, "inp.vcf", genTestVCFProbs)

	cases := []struct {
		name string
		src  string
		opts GenSampleOptions
		args []string
	}{
		{"plain", in, GenSampleOptions{}, nil},
		{"vcfids", in, GenSampleOptions{VCFIDs: true}, []string{"--vcf-ids"}},
		{"3N6", in, GenSampleOptions{ThreeN6: true}, []string{"--3N6"}},
		{"gp", inP, GenSampleOptions{Tag: "GP"}, []string{"--tag", "GP"}},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			// Produce the .gen/.sample with the upstream binary so both sides
			// import identical input.
			pfx := filepath.Join(dir, "src_"+c.name)
			genArgs := append([]string{"convert", "-g", pfx}, c.args...)
			genArgs = append(genArgs, c.src)
			if out, err := exec.Command(bin, genArgs...).CombinedOutput(); err != nil {
				t.Fatalf("upstream convert -g: %v\n%s", err, out)
			}

			// upstream -G
			upArgs := append([]string{"convert", "-G", pfx, "--no-version", "-Ov"}, c.args...)
			upOut, err := exec.Command(bin, upArgs...).Output()
			if err != nil {
				t.Fatalf("upstream convert -G: %v", err)
			}
			// go -G
			var goBuf bytes.Buffer
			impOpts := c.opts
			impOpts.OutputFormat = OutputVCF
			impOpts.NoVersion = true
			if err := GenSampleToVCFFile(pfx, &goBuf, impOpts); err != nil {
				t.Fatalf("go convert -G: %v", err)
			}
			if !bytes.Equal(upOut, goBuf.Bytes()) {
				t.Fatalf("-G output mismatch:\n--- up ---\n%s\n--- go ---\n%s", upOut, goBuf.String())
			}
		})
	}
}
