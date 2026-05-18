package vcftools

// Tests for --recode-bcf (the BCF output sibling of --recode).
//
// BGZF byte parity against upstream is not a sensible target — deflate
// output is implementation-dependent and upstream's vcftools 0.1.18 has
// known bugs in its BCF emission (phantom unnamed INFO field, see
// `docs/UPSTREAM_BUGS.md`). Instead we exercise:
//
//  1. The port can emit a BCF that round-trips through its own reader.
//  2. The emitted BCF has the upstream-standard `,IDX=N` annotations in
//     the text header so an htslib-compatible reader picks the right
//     dictionary numbering.
//  3. The on-wire FORMAT field uses per-sample dimension as the
//     descriptor's `size` (the spec-correct interpretation; htslib reads
//     it as such).
//
// `TestParity_RecodeBCF_Roundtrip_Upstream` is gated behind the
// `/tmp/vcftools_install/bin/vcftools` binary and round-trips through
// upstream too — skipped when the binary isn't built.

import (
	"bytes"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/bcf"
	"github.com/yassineS/bio_ai_experiment/tools/bgzip/pkg/bgzip"
)

const recodeBCFFixture = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000>
##contig=<ID=chr2,length=2000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1	s2
chr1	100	.	A	G	30	PASS	DP=20	GT	0/0	0/1
chr1	200	.	C	T	30	PASS	DP=25	GT	0/1	1/1
chr2	150	.	G	A	30	PASS	DP=18	GT	0/1	0/0
`

// readBCFViaPort opens path, BGZF-decompresses, decodes the BCF through
// our own reader, and returns the resulting variants.
func readBCFViaPort(t *testing.T, path string) (samples []string, variants []decodedVariant) {
	t.Helper()
	f, err := os.Open(path)
	if err != nil {
		t.Fatalf("open bcf: %v", err)
	}
	defer f.Close()
	bz, err := bgzip.NewReader(f)
	if err != nil {
		t.Fatalf("bgzip reader: %v", err)
	}
	r, err := bcf.NewReader(bz)
	if err != nil {
		t.Fatalf("bcf reader: %v", err)
	}
	samples = append([]string(nil), r.Header().Samples...)
	for {
		rec, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("bcf read: %v", err)
		}
		v := rec.ToVariant(r.Header())
		dv := decodedVariant{Chrom: v.Chrom, Pos: v.Pos, Ref: v.Ref, Alt: append([]string(nil), v.Alt...)}
		for _, s := range v.Samples {
			dv.Samples = append(dv.Samples, sampleGT{Name: s.Name, GT: s.Data["GT"]})
		}
		variants = append(variants, dv)
	}
	return samples, variants
}

type decodedVariant struct {
	Chrom   string
	Pos     int
	Ref     string
	Alt     []string
	Samples []sampleGT
}

type sampleGT struct {
	Name string
	GT   string
}

// TestRun_RecodeBCF_Roundtrip writes a fresh BCF via --recode-bcf and
// reads it back through our own decoder. The decoded variants must match
// the input verbatim on chrom/pos/ref/alt and per-sample GT.
func TestRun_RecodeBCF_Roundtrip(t *testing.T) {
	tmp := t.TempDir()
	prefix := filepath.Join(tmp, "out")
	if err := Run(strings.NewReader(recodeBCFFixture), &Params{OutPrefix: prefix, RecodeBCF: true}); err != nil {
		t.Fatalf("Run --recode-bcf: %v", err)
	}
	samples, vs := readBCFViaPort(t, prefix+".recode.bcf")
	if got, want := samples, []string{"s1", "s2"}; !equalSS(got, want) {
		t.Errorf("samples: %v, want %v", got, want)
	}
	want := []decodedVariant{
		{Chrom: "chr1", Pos: 100, Ref: "A", Alt: []string{"G"}, Samples: []sampleGT{{Name: "s1", GT: "0/0"}, {Name: "s2", GT: "0/1"}}},
		{Chrom: "chr1", Pos: 200, Ref: "C", Alt: []string{"T"}, Samples: []sampleGT{{Name: "s1", GT: "0/1"}, {Name: "s2", GT: "1/1"}}},
		{Chrom: "chr2", Pos: 150, Ref: "G", Alt: []string{"A"}, Samples: []sampleGT{{Name: "s1", GT: "0/1"}, {Name: "s2", GT: "0/0"}}},
	}
	if len(vs) != len(want) {
		t.Fatalf("variant count: got %d want %d", len(vs), len(want))
	}
	for i := range want {
		if vs[i].Chrom != want[i].Chrom || vs[i].Pos != want[i].Pos ||
			vs[i].Ref != want[i].Ref || !equalSS(vs[i].Alt, want[i].Alt) {
			t.Errorf("variant %d shape mismatch: got %+v want %+v", i, vs[i], want[i])
			continue
		}
		if len(vs[i].Samples) != len(want[i].Samples) {
			t.Errorf("variant %d sample count: got %d want %d", i, len(vs[i].Samples), len(want[i].Samples))
			continue
		}
		for j := range want[i].Samples {
			if vs[i].Samples[j] != want[i].Samples[j] {
				t.Errorf("variant %d sample %d: got %+v want %+v", i, j, vs[i].Samples[j], want[i].Samples[j])
			}
		}
	}
}

// TestRun_RecodeBCF_HeaderHasIDXAnnotations checks the emitted text header
// carries `,IDX=N` on its INFO/FORMAT lines so htslib (and any other
// downstream consumer) picks the same unified dictionary numbering the
// port used for wire indices.
func TestRun_RecodeBCF_HeaderHasIDXAnnotations(t *testing.T) {
	tmp := t.TempDir()
	prefix := filepath.Join(tmp, "out")
	if err := Run(strings.NewReader(recodeBCFFixture), &Params{OutPrefix: prefix, RecodeBCF: true}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	f, err := os.Open(prefix + ".recode.bcf")
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	bz, err := bgzip.NewReader(f)
	if err != nil {
		t.Fatal(err)
	}
	body, err := io.ReadAll(bz)
	if err != nil {
		t.Fatal(err)
	}
	// Skip the BCF magic + length-prefixed text header and inspect the
	// text bytes. Magic is 5 bytes; next 4 bytes are uint32 LE text
	// length; then `length` bytes of text including a trailing NUL.
	if !bytes.HasPrefix(body, bcf.Magic[:]) {
		t.Fatal("missing BCF magic")
	}
	text := string(body[9:])
	idx := strings.IndexByte(text, '\x00')
	if idx >= 0 {
		text = text[:idx]
	}
	for _, want := range []string{",IDX=", "##INFO=<ID=DP", "##FORMAT=<ID=GT"} {
		if !strings.Contains(text, want) {
			t.Errorf("header missing %q in:\n%s", want, text)
		}
	}
}

func equalSS(a, b []string) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

// multiFormatBCFFixture exercises the FORMAT decoder beyond the GT
// fast-path: integer scalars (DP, GQ), the special "Number=G" likelihood
// vector (PL), an explicit ##FILTER=<ID=PASS,...> declaration (the
// implicit-IDX=0 invariant breaker the wave-21 review flagged), and
// per-sample multi-value INFO (AF as Number=A).
const multiFormatBCFFixture = `##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##contig=<ID=chr1,length=10000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Quality">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="PL likelihoods">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	s1	s2
chr1	100	.	A	G	30	PASS	DP=20;AF=0.5	GT:DP:GQ:PL	0/0:8:60:0,15,200	0/1:12:80:50,0,180
chr1	200	.	C	T	25	PASS	DP=15;AF=0.33	GT:DP:GQ:PL	0/1:7:55:30,0,150	1/1:9:70:200,15,0
`

// TestRun_RecodeBCF_MultiFormatRoundtrip writes a BCF with several
// non-GT FORMAT fields, decodes it through our reader, and asserts the
// per-sample data survives the GT/integer/likelihood-vector encoders.
// Regression for the wave-21 review's "encodeRecord uses unified IDX as
// slice position → non-GT FORMAT fields lost" bug.
func TestRun_RecodeBCF_MultiFormatRoundtrip(t *testing.T) {
	tmp := t.TempDir()
	prefix := filepath.Join(tmp, "out")
	if err := Run(strings.NewReader(multiFormatBCFFixture), &Params{OutPrefix: prefix, RecodeBCF: true, RecodeInfoAll: true}); err != nil {
		t.Fatalf("Run --recode-bcf: %v", err)
	}

	// Use the package's own bcf reader directly so we can inspect the
	// per-FORMAT-key data (the existing readBCFViaPort helper only
	// surfaces GT). Round-trip through the bgzip layer.
	f, err := os.Open(prefix + ".recode.bcf")
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	bz, err := bgzip.NewReader(f)
	if err != nil {
		t.Fatal(err)
	}
	r, err := bcf.NewReader(bz)
	if err != nil {
		t.Fatal(err)
	}
	var got []*decodedFullVariant
	for {
		rec, err := r.Read()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatalf("Read: %v", err)
		}
		v := rec.ToVariant(r.Header())
		dv := &decodedFullVariant{
			Chrom: v.Chrom, Pos: v.Pos, Ref: v.Ref, Alt: append([]string(nil), v.Alt...),
			Format: append([]string(nil), v.Format...),
		}
		for _, s := range v.Samples {
			row := map[string]string{}
			for k, val := range s.Data {
				row[k] = val
			}
			dv.Samples = append(dv.Samples, fullSample{Name: s.Name, Data: row})
		}
		got = append(got, dv)
	}
	if len(got) != 2 {
		t.Fatalf("variant count: got %d want 2", len(got))
	}

	// Site 1 assertions
	v := got[0]
	if v.Chrom != "chr1" || v.Pos != 100 {
		t.Errorf("v1: chrom/pos: %s:%d", v.Chrom, v.Pos)
	}
	checks := []struct{ sample, key, want string }{
		{"s1", "GT", "0/0"}, {"s1", "DP", "8"}, {"s1", "GQ", "60"}, {"s1", "PL", "0,15,200"},
		{"s2", "GT", "0/1"}, {"s2", "DP", "12"}, {"s2", "GQ", "80"}, {"s2", "PL", "50,0,180"},
	}
	for _, c := range checks {
		sample := findFullSample(v.Samples, c.sample)
		if sample == nil {
			t.Errorf("site 1 missing sample %s", c.sample)
			continue
		}
		if got := sample.Data[c.key]; got != c.want {
			t.Errorf("site 1 %s.%s: got %q want %q", c.sample, c.key, got, c.want)
		}
	}

	// Site 2 assertions
	v = got[1]
	checks = []struct{ sample, key, want string }{
		{"s1", "GT", "0/1"}, {"s1", "DP", "7"}, {"s1", "GQ", "55"}, {"s1", "PL", "30,0,150"},
		{"s2", "GT", "1/1"}, {"s2", "DP", "9"}, {"s2", "GQ", "70"}, {"s2", "PL", "200,15,0"},
	}
	for _, c := range checks {
		sample := findFullSample(v.Samples, c.sample)
		if sample == nil {
			t.Errorf("site 2 missing sample %s", c.sample)
			continue
		}
		if got := sample.Data[c.key]; got != c.want {
			t.Errorf("site 2 %s.%s: got %q want %q", c.sample, c.key, got, c.want)
		}
	}
}

type decodedFullVariant struct {
	Chrom   string
	Pos     int
	Ref     string
	Alt     []string
	Format  []string
	Samples []fullSample
}

type fullSample struct {
	Name string
	Data map[string]string
}

func findFullSample(samples []fullSample, name string) *fullSample {
	for i := range samples {
		if samples[i].Name == name {
			return &samples[i]
		}
	}
	return nil
}

// TestRun_RecodeBCF_ExplicitPASS_DoesNotShiftIDX verifies that an input
// VCF with an explicit `##FILTER=<ID=PASS,...>` declaration produces a
// BCF whose dictionary numbering still places PASS at IDX=0 (matching
// the htslib implicit-PASS convention). Regression for the wave-21
// review's "explicit PASS shifts every subsequent IDX by one" bug.
func TestRun_RecodeBCF_ExplicitPASS_DoesNotShiftIDX(t *testing.T) {
	tmp := t.TempDir()
	prefix := filepath.Join(tmp, "out")
	if err := Run(strings.NewReader(multiFormatBCFFixture), &Params{OutPrefix: prefix, RecodeBCF: true, RecodeInfoAll: true}); err != nil {
		t.Fatalf("Run: %v", err)
	}
	f, err := os.Open(prefix + ".recode.bcf")
	if err != nil {
		t.Fatal(err)
	}
	defer f.Close()
	bz, err := bgzip.NewReader(f)
	if err != nil {
		t.Fatal(err)
	}
	body, err := io.ReadAll(bz)
	if err != nil {
		t.Fatal(err)
	}
	text := string(body[9:])
	if i := strings.IndexByte(text, '\x00'); i >= 0 {
		text = text[:i]
	}
	// The emitted header must NOT contain `##FILTER=<ID=PASS` — the
	// implicit-IDX=0 PASS is enough and the explicit line would otherwise
	// duplicate the entry under a non-zero IDX.
	if strings.Contains(text, "##FILTER=<ID=PASS") {
		t.Errorf("emitted header should drop explicit PASS line; got:\n%s", text)
	}
	// First non-PASS FILTER/INFO/FORMAT line should be IDX=1, not IDX=2.
	// We expect INFO/DP at IDX=1 because the source declares it first.
	if !strings.Contains(text, "##INFO=<ID=DP,") || !strings.Contains(text, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\",IDX=1>") {
		t.Errorf("INFO/DP should be IDX=1; got:\n%s", text)
	}
}
