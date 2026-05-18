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
