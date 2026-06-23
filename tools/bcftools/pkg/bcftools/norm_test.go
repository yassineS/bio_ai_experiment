package bcftools

import (
	"bytes"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/bcf"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/fasta"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// testBCFWriter is a thin wrapper around bcf.Writer used by the BCF input
// test. Pulling the helper out keeps the test body focused.
type testBCFWriter struct {
	w *bcf.Writer
}

func newTestBCFWriter(out io.Writer, vh *vcf.Header) (*testBCFWriter, error) {
	w, err := bcf.NewWriterFromVCFHeader(out, vh)
	if err != nil {
		return nil, err
	}
	if err := w.WriteHeader(); err != nil {
		return nil, err
	}
	return &testBCFWriter{w: w}, nil
}

// MustWriteVariants encodes every variant; fatal on error.
func (w *testBCFWriter) MustWriteVariants(t *testing.T, variants []*vcf.Variant) {
	t.Helper()
	for _, v := range variants {
		if err := w.w.Write(v); err != nil {
			t.Fatalf("bcf write: %v", err)
		}
	}
	if err := w.w.Flush(); err != nil {
		t.Fatalf("bcf flush: %v", err)
	}
}

// writeRefFasta writes a tiny FASTA + .fai pair into a temp dir and
// returns its path.
func writeRefFasta(t *testing.T, contigs map[string]string) string {
	t.Helper()
	dir := t.TempDir()
	path := filepath.Join(dir, "ref.fa")
	var buf bytes.Buffer
	// Iterate contigs in deterministic order.
	names := make([]string, 0, len(contigs))
	for name := range contigs {
		names = append(names, name)
	}
	// Sort for determinism.
	for i := 0; i < len(names); i++ {
		for j := i + 1; j < len(names); j++ {
			if names[j] < names[i] {
				names[i], names[j] = names[j], names[i]
			}
		}
	}
	for _, name := range names {
		buf.WriteString(">")
		buf.WriteString(name)
		buf.WriteByte('\n')
		seq := contigs[name]
		// Write 60 bases per line; the builder also accepts a single
		// long line so we keep it simple here.
		for len(seq) > 0 {
			n := 60
			if n > len(seq) {
				n = len(seq)
			}
			buf.WriteString(seq[:n])
			buf.WriteByte('\n')
			seq = seq[n:]
		}
	}
	if err := os.WriteFile(path, buf.Bytes(), 0o644); err != nil {
		t.Fatalf("write ref: %v", err)
	}
	idx, err := fasta.BuildIndex(path)
	if err != nil {
		t.Fatalf("build idx: %v", err)
	}
	if err := idx.Save(path + ".fai"); err != nil {
		t.Fatalf("save idx: %v", err)
	}
	return path
}

// vcfDoc concatenates a minimal header with body lines into a VCF string.
func vcfDoc(samples []string, body ...string) string {
	var b strings.Builder
	b.WriteString("##fileformat=VCFv4.2\n")
	b.WriteString("##INFO=<ID=AC,Number=A,Type=Integer,Description=\"AC\">\n")
	b.WriteString("##INFO=<ID=AN,Number=1,Type=Integer,Description=\"AN\">\n")
	b.WriteString("##INFO=<ID=AF,Number=A,Type=Float,Description=\"AF\">\n")
	b.WriteString("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
	b.WriteString("##contig=<ID=chr1>\n")
	b.WriteString("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
	if len(samples) > 0 {
		b.WriteString("\tFORMAT")
		for _, s := range samples {
			b.WriteByte('\t')
			b.WriteString(s)
		}
	}
	b.WriteByte('\n')
	for _, line := range body {
		b.WriteString(line)
		b.WriteByte('\n')
	}
	return b.String()
}

// runNorm is a small wrapper that calls Norm on a VCF string.
func runNorm(t *testing.T, input string, opts NormOptions) (string, string, int) {
	t.Helper()
	var out, stderr bytes.Buffer
	n, err := Norm(strings.NewReader(input), &out, opts, &stderr)
	if err != nil {
		t.Fatalf("Norm: %v\nstderr: %s", err, stderr.String())
	}
	return out.String(), stderr.String(), n
}

func TestNormPassthroughNoFlags(t *testing.T) {
	input := vcfDoc(nil, "chr1\t10\t.\tA\tG\t.\tPASS\t.")
	out, _, n := runNorm(t, input, NormOptions{})
	if n != 1 {
		t.Fatalf("n=%d", n)
	}
	if !strings.Contains(out, "chr1\t10\t.\tA\tG") {
		t.Fatalf("output missing the record:\n%s", out)
	}
}

// TestNormOrdersByHeaderContigNotLexical is a regression test for the
// multi-contig ordering bug the 16-contig validation tier surfaced: norm must
// order output records by the VCF header's ##contig declaration order (rid),
// not by contig name as a string. With >=10 contigs, "chr10" sorts before
// "chr2" lexically, so a string sort wrongly emits chr10 before chr2; upstream
// bcftools preserves header order (chr1, chr2, ..., chr10). The records are fed
// out of header order to prove the sort actively reorders to header order.
func TestNormOrdersByHeaderContigNotLexical(t *testing.T) {
	input := strings.Join([]string{
		"##fileformat=VCFv4.2",
		"##contig=<ID=chr1>",
		"##contig=<ID=chr2>",
		"##contig=<ID=chr10>",
		"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
		"chr10\t30\t.\tG\tA\t.\tPASS\t.",
		"chr2\t20\t.\tC\tT\t.\tPASS\t.",
		"chr1\t10\t.\tA\tG\t.\tPASS\t.",
		"",
	}, "\n")
	out, _, _ := runNorm(t, input, NormOptions{})
	var got []string
	for _, line := range strings.Split(out, "\n") {
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}
		got = append(got, strings.SplitN(line, "\t", 2)[0])
	}
	if order := strings.Join(got, ","); order != "chr1,chr2,chr10" {
		t.Fatalf("contig order = %q, want \"chr1,chr2,chr10\" (header order, not lexical)", order)
	}
}

func TestNormLeftAlignSingleIndel(t *testing.T) {
	// Reference ATTTTG (1-based positions 1..6). A T-deletion can be
	// represented anywhere in the TTTT run; the canonical (left-aligned)
	// form puts it at the leftmost flanking-base position.
	//   variant: REF=TT ALT=T at pos 4 (the last two Ts).
	//   standard Tan-Abecasis-Durbin trimming:
	//     [TT,T] -> last bytes match, both not all len>1 -> trim?
	//     Our impl requires both len>1; len(T)==1 so we drop into the
	//     "any empty" branch only after stripping. Use the alternate
	//     stripping: shared trailing "T" -> [T,""], prepend ref[3]=T ->
	//     [TT,T] at pos 3. Repeat -> pos 2. Repeat -> [AT,A] at pos 1.
	ref := writeRefFasta(t, map[string]string{"chr1": "ATTTTG"})
	input := vcfDoc(nil, "chr1\t4\t.\tTT\tT\t.\tPASS\t.")
	out, _, _ := runNorm(t, input, NormOptions{FastaRef: ref})
	want := "chr1\t1\t.\tAT\tA"
	if !strings.Contains(out, want) {
		t.Fatalf("expected left-aligned %q in output, got:\n%s", want, out)
	}
}

func TestNormLeftAlignDoubleStep(t *testing.T) {
	// CC deletion in a CCC run: ACCCG. REF=CC ALT=C at pos 3 should
	// left-align to pos 1 (REF=AC ALT=A).
	ref := writeRefFasta(t, map[string]string{"chr1": "ACCCG"})
	input := vcfDoc(nil, "chr1\t3\t.\tCC\tC\t.\tPASS\t.")
	out, _, _ := runNorm(t, input, NormOptions{FastaRef: ref})
	want := "chr1\t1\t.\tAC\tA"
	if !strings.Contains(out, want) {
		t.Fatalf("expected left-aligned %q in output, got:\n%s", want, out)
	}
}

func TestNormLeftAlignInsertion(t *testing.T) {
	// AAACCCAAA (positions 1..9 are A,A,A,C,C,C,A,A,A). An insertion of
	// "C" in the CCC run can be written as REF=C ALT=CC at any of pos
	// 4, 5, or 6. Left-aligning should land at the leftmost flanking A:
	// REF=A ALT=AC at pos 3 (the A immediately before the run, used as
	// the VCF anchor base).
	ref := writeRefFasta(t, map[string]string{"chr1": "AAACCCAAA"})
	input := vcfDoc(nil, "chr1\t5\t.\tC\tCC\t.\tPASS\t.")
	out, _, _ := runNorm(t, input, NormOptions{FastaRef: ref})
	want := "chr1\t3\t.\tA\tAC"
	if !strings.Contains(out, want) {
		t.Fatalf("expected left-aligned %q in output, got:\n%s", want, out)
	}
}

func TestNormCheckRefModes(t *testing.T) {
	ref := writeRefFasta(t, map[string]string{"chr1": "AAAAAAAA"})

	// REF says "C" but FASTA says "A" at pos 4.
	input := vcfDoc(nil, "chr1\t4\t.\tC\tT\t.\tPASS\t.")

	// error mode → Norm returns an error.
	var out, stderr bytes.Buffer
	_, err := Norm(strings.NewReader(input), &out, NormOptions{
		FastaRef: ref,
		CheckRef: CheckRefError,
	}, &stderr)
	if err == nil {
		t.Fatalf("CheckRefError should fail")
	}

	// warn mode → record passes, stderr mentions mismatch.
	out.Reset()
	stderr.Reset()
	n, err := Norm(strings.NewReader(input), &out, NormOptions{
		FastaRef: ref,
		CheckRef: CheckRefWarn,
	}, &stderr)
	if err != nil {
		t.Fatalf("warn: %v", err)
	}
	if n != 1 {
		t.Fatalf("warn n=%d", n)
	}
	if !strings.Contains(stderr.String(), "REF_MISMATCH") {
		t.Fatalf("warn stderr missing message:\n%s", stderr.String())
	}

	// skip mode → record is dropped silently.
	out.Reset()
	stderr.Reset()
	n, err = Norm(strings.NewReader(input), &out, NormOptions{
		FastaRef: ref,
		CheckRef: CheckRefSkip,
	}, &stderr)
	if err != nil {
		t.Fatalf("skip: %v", err)
	}
	if n != 0 {
		t.Fatalf("skip n=%d, want 0", n)
	}
	if stderr.String() != "" {
		t.Fatalf("skip stderr should be empty: %q", stderr.String())
	}
}

func TestNormSplitMultiallelicSNPs(t *testing.T) {
	body := []string{
		"chr1\t10\t.\tA\tG,C\t.\tPASS\tAC=2,1;AN=4\tGT\t1/2",
	}
	input := vcfDoc([]string{"S1"}, body...)
	opts := NormOptions{
		Multiallelics: MultiallelicMode{Active: true, Split: true, Snps: true},
	}
	out, _, n := runNorm(t, input, opts)
	if n != 2 {
		t.Fatalf("split n=%d, want 2", n)
	}
	if !strings.Contains(out, "chr1\t10\t.\tA\tG") {
		t.Fatalf("missing G split:\n%s", out)
	}
	if !strings.Contains(out, "chr1\t10\t.\tA\tC") {
		t.Fatalf("missing C split:\n%s", out)
	}
	// Per-allele AC for the first record (allele G) should be 2.
	lines := strings.Split(strings.TrimSpace(out), "\n")
	for _, line := range lines {
		if strings.HasPrefix(line, "chr1\t10") && strings.Contains(line, "\tG\t") {
			if !strings.Contains(line, "AC=2") {
				t.Errorf("G record should carry AC=2: %s", line)
			}
			// FORMAT/GT: original 1/2 with G chosen → 1/0
			if !strings.HasSuffix(line, "\t1/0") {
				t.Errorf("G record GT should be 1/0: %s", line)
			}
		}
		if strings.HasPrefix(line, "chr1\t10") && strings.Contains(line, "\tC\t") {
			if !strings.Contains(line, "AC=1") {
				t.Errorf("C record should carry AC=1: %s", line)
			}
			// FORMAT/GT: original 1/2 with C chosen → 0/1
			if !strings.HasSuffix(line, "\t0/1") {
				t.Errorf("C record GT should be 0/1: %s", line)
			}
		}
	}
}

func TestNormJoinMultiallelicSNPs(t *testing.T) {
	body := []string{
		"chr1\t10\t.\tA\tG\t.\tPASS\tAC=2;AN=4\tGT\t0/1",
		"chr1\t10\t.\tA\tC\t.\tPASS\tAC=1;AN=4\tGT\t0/1",
	}
	input := vcfDoc([]string{"S1"}, body...)
	opts := NormOptions{
		Multiallelics: MultiallelicMode{Active: true, Split: false, Snps: true},
	}
	out, _, n := runNorm(t, input, opts)
	if n != 1 {
		t.Fatalf("join n=%d, want 1", n)
	}
	// Should have a single multiallelic with ALT=G,C and AC=2,1
	if !strings.Contains(out, "\tA\tG,C\t") {
		t.Fatalf("expected joined ALT 'G,C':\n%s", out)
	}
	if !strings.Contains(out, "AC=2,1") {
		t.Errorf("expected AC=2,1:\n%s", out)
	}
	// GT for the donor records: 0/1 (A>G) + 0/1 (A>C). Upstream's
	// merge_format_genotype keeps the first record's non-ref allele (the G at
	// slot 1) in place and writes the second record's allele (C, merged index
	// 2) into the first free strand, which is the leading ref slot — yielding
	// "2/1", not "0/2". (Verified byte-for-byte against bcftools 1.23.1.)
	if !strings.HasSuffix(strings.TrimRight(out, "\n"), "\t2/1") {
		t.Errorf("expected joined GT 2/1:\n%s", out)
	}
}

func TestNormRmDup(t *testing.T) {
	body := []string{
		"chr1\t10\t.\tA\tG\t.\tPASS\t.",
		"chr1\t10\t.\tA\tG\t.\tPASS\t.",
		"chr1\t11\t.\tA\tC\t.\tPASS\t.",
	}
	input := vcfDoc(nil, body...)
	// none = no-op
	_, _, n := runNorm(t, input, NormOptions{RmDup: RmDupNone})
	if n != 3 {
		t.Errorf("none n=%d", n)
	}
	// exact = drop second 10:A>G
	_, _, n = runNorm(t, input, NormOptions{RmDup: RmDupExact})
	if n != 2 {
		t.Errorf("exact n=%d", n)
	}
}

func TestNormAtomize(t *testing.T) {
	// REF=ACG ALT=AGT — first base matches; positions 2 and 3 mismatch.
	// Atomize → two single-bp variants at pos 2 (C>G) and pos 3 (G>T).
	body := []string{"chr1\t10\t.\tACG\tAGT\t.\tPASS\t."}
	input := vcfDoc(nil, body...)
	out, _, n := runNorm(t, input, NormOptions{Atomize: true})
	if n != 2 {
		t.Fatalf("atomize n=%d, want 2:\n%s", n, out)
	}
	if !strings.Contains(out, "chr1\t11\t.\tC\tG") {
		t.Errorf("missing C>G:\n%s", out)
	}
	if !strings.Contains(out, "chr1\t12\t.\tG\tT") {
		t.Errorf("missing G>T:\n%s", out)
	}
}

func TestNormRegionFilter(t *testing.T) {
	body := []string{
		"chr1\t50\t.\tA\tG\t.\tPASS\t.",
		"chr1\t150\t.\tA\tG\t.\tPASS\t.",
		"chr1\t250\t.\tA\tG\t.\tPASS\t.",
	}
	input := vcfDoc(nil, body...)
	opts := NormOptions{Regions: []string{"chr1:100-200"}}
	out, _, n := runNorm(t, input, opts)
	if n != 1 {
		t.Fatalf("region n=%d, want 1", n)
	}
	if !strings.Contains(out, "\t150\t") {
		t.Errorf("expected 150 to survive:\n%s", out)
	}
}

func TestParseCheckRefMode(t *testing.T) {
	cases := []struct {
		in   string
		want CheckRefMode
		err  bool
	}{
		{"", CheckRefExit, false},
		{"e", CheckRefExit, false},
		{"w", CheckRefWarn, false},
		{"x", CheckRefSkip, false},
		{"s", CheckRefFix, false},
		{"ws", CheckRefWarn | CheckRefFix, false},
		{"wx", CheckRefWarn | CheckRefSkip, false},
		{"q", 0, true},
	}
	for _, c := range cases {
		got, err := ParseCheckRefMode(c.in)
		if c.err {
			if err == nil {
				t.Errorf("%q want err", c.in)
			}
			continue
		}
		if err != nil {
			t.Errorf("%q unexpected err: %v", c.in, err)
		}
		if got != c.want {
			t.Errorf("%q got %v want %v", c.in, got, c.want)
		}
	}
}

func TestParseMultiallelicMode(t *testing.T) {
	cases := []struct {
		in   string
		want MultiallelicMode
		err  bool
	}{
		{"", MultiallelicMode{}, false},
		{"-snps", MultiallelicMode{Active: true, Split: true, Snps: true}, false},
		{"+indels", MultiallelicMode{Active: true, Split: false, Indels: true}, false},
		{"-both", MultiallelicMode{Active: true, Split: true, Snps: true, Indels: true}, false},
		// "+any" (COLLAPSE_ANY) sets Any so the joiner merges every variant
		// type into one record, distinct from "+both" which buckets by type.
		{"+any", MultiallelicMode{Active: true, Split: false, Snps: true, Indels: true, Any: true}, false},
		{"-any", MultiallelicMode{Active: true, Split: true, Snps: true, Indels: true, Any: true}, false},
		{"snps", MultiallelicMode{}, true},
		{"-other", MultiallelicMode{}, true},
	}
	for _, c := range cases {
		got, err := ParseMultiallelicMode(c.in)
		if c.err {
			if err == nil {
				t.Errorf("%q want err", c.in)
			}
			continue
		}
		if err != nil {
			t.Errorf("%q unexpected err: %v", c.in, err)
			continue
		}
		if got != c.want {
			t.Errorf("%q got %+v want %+v", c.in, got, c.want)
		}
	}
}

func TestParseRmDupMode(t *testing.T) {
	cases := []struct {
		in   string
		want RmDupMode
		err  bool
	}{
		{"", RmDupNone, false},
		{"none", RmDupNone, false},
		{"snps", RmDupSnps, false},
		{"indels", RmDupIndels, false},
		{"both", RmDupBoth, false},
		{"all", RmDupAll, false},
		{"exact", RmDupExact, false},
		{"bogus", RmDupNone, true},
	}
	for _, c := range cases {
		got, err := ParseRmDupMode(c.in)
		if c.err {
			if err == nil {
				t.Errorf("%q want err", c.in)
			}
			continue
		}
		if err != nil {
			t.Errorf("%q unexpected err: %v", c.in, err)
		}
		if got != c.want {
			t.Errorf("%q got %v want %v", c.in, got, c.want)
		}
	}
}

func TestNormStrictFilter(t *testing.T) {
	body := []string{
		"chr1\t10\t.\tA\tG,C\t.\tLowQual\t.",
	}
	input := vcfDoc(nil, body...)
	// With strict-filter ON and apply-filters=PASS we should drop the
	// record before any splitting.
	opts := NormOptions{
		StrictFilter:  true,
		ApplyFilters:  []string{"PASS"},
		Multiallelics: MultiallelicMode{Active: true, Split: true, Snps: true},
	}
	_, _, n := runNorm(t, input, opts)
	if n != 0 {
		t.Fatalf("strict filter should have dropped record; n=%d", n)
	}
}

func TestNormLaxFilterAfterSplit(t *testing.T) {
	body := []string{
		"chr1\t10\t.\tA\tG,C\t.\tPASS\t.",
		"chr1\t20\t.\tA\tT\t.\tLowQual\t.",
	}
	input := vcfDoc(nil, body...)
	opts := NormOptions{
		ApplyFilters:  []string{"PASS"},
		Multiallelics: MultiallelicMode{Active: true, Split: true, Snps: true},
	}
	_, _, n := runNorm(t, input, opts)
	// After splitting we have 3 records (2 from the multiallelic + 1
	// LowQual). The PASS filter drops the LowQual one. Expect 2.
	if n != 2 {
		t.Fatalf("lax filter n=%d, want 2", n)
	}
}

func TestNormDoNotNormalizeKeepsCoordinates(t *testing.T) {
	// ATTTTG (pos1=A,2=T,3=T,4=T,5=T,6=G); REF=TT at pos 4 matches bases
	// 4,5 = T,T. Left-align would push to pos 1 (REF=AT ALT=A); with -N
	// we expect the coordinates to stay put.
	ref := writeRefFasta(t, map[string]string{"chr1": "ATTTTG"})
	body := []string{"chr1\t4\t.\tTT\tT\t.\tPASS\t."}
	input := vcfDoc(nil, body...)
	out, _, _ := runNorm(t, input, NormOptions{
		FastaRef:       ref,
		DoNotNormalize: true,
	})
	if !strings.Contains(out, "chr1\t4\t.\tTT\tT") {
		t.Fatalf("-N should keep pos 4:\n%s", out)
	}
}

func TestRemapGT(t *testing.T) {
	cases := []struct {
		gt, wanted, want string
	}{
		{"1/2", "2", "0/1"},
		{"1|2", "2", "0|1"},
		{"0/0", "1", "0/0"},
		{"0/1", "1", "0/1"},
		{"./.", "1", "./."},
		{"", "1", ""},
		{"2|3", "3", "0|1"},
	}
	for _, c := range cases {
		got := remapGT(c.gt, c.wanted)
		if got != c.want {
			t.Errorf("remapGT(%q, %q) = %q want %q", c.gt, c.wanted, got, c.want)
		}
	}
}

func TestLeftAlignDoesNotRunPastChromStart(t *testing.T) {
	// Reference "CCG" puts a CC at pos 1-2 followed by G. REF=CC ALT=C
	// at pos 1 would normally shift left, but there's no upstream base
	// to consume — we expect an error.
	ref := writeRefFasta(t, map[string]string{"chr1": "CCG"})
	body := []string{"chr1\t1\t.\tCC\tC\t.\tPASS\t."}
	input := vcfDoc(nil, body...)
	var out, stderr bytes.Buffer
	_, err := Norm(strings.NewReader(input), &out, NormOptions{FastaRef: ref}, &stderr)
	if err == nil {
		t.Fatalf("expected error when aligning past chrom start")
	}
}

func TestVariantIsSnpAndIndel(t *testing.T) {
	v := &vcf.Variant{Ref: "A", Alt: []string{"G"}}
	if !variantIsSnp(v) {
		t.Fatalf("expected SNP")
	}
	if variantIsIndel(v) {
		t.Fatalf("expected not indel")
	}
	v2 := &vcf.Variant{Ref: "A", Alt: []string{"AT"}}
	if variantIsSnp(v2) {
		t.Fatalf("expected not SNP")
	}
	if !variantIsIndel(v2) {
		t.Fatalf("expected indel")
	}
	v3 := &vcf.Variant{Ref: "A", Alt: nil}
	if variantIsSnp(v3) {
		t.Fatalf("empty ALT should not be SNP")
	}
}

func TestRmDupSnps(t *testing.T) {
	body := []string{
		"chr1\t10\t.\tA\tG\t.\tPASS\t.",
		"chr1\t10\t.\tA\tT\t.\tPASS\t.",    // same pos, SNP — dropped
		"chr1\t10\t.\tACGT\tA\t.\tPASS\t.", // same pos, indel — kept
	}
	input := vcfDoc(nil, body...)
	_, _, n := runNorm(t, input, NormOptions{RmDup: RmDupSnps})
	if n != 2 {
		t.Fatalf("rm-dup snps n=%d, want 2", n)
	}
}

func TestRmDupBoth(t *testing.T) {
	body := []string{
		"chr1\t10\t.\tA\tG\t.\tPASS\t.",
		"chr1\t10\t.\tA\tT\t.\tPASS\t.",
		"chr1\t10\t.\tACGT\tA\t.\tPASS\t.",
		"chr1\t10\t.\tACGT\tACGTACGT\t.\tPASS\t.",
	}
	input := vcfDoc(nil, body...)
	_, _, n := runNorm(t, input, NormOptions{RmDup: RmDupBoth})
	// SNPs collapse to 1, indels collapse to 1.
	if n != 2 {
		t.Fatalf("rm-dup both n=%d, want 2", n)
	}
}

func TestRmDupAll(t *testing.T) {
	body := []string{
		"chr1\t10\t.\tA\tG\t.\tPASS\t.",
		"chr1\t10\t.\tA\tT\t.\tPASS\t.",
		"chr1\t10\t.\tACGT\tA\t.\tPASS\t.",
	}
	input := vcfDoc(nil, body...)
	_, _, n := runNorm(t, input, NormOptions{RmDup: RmDupAll})
	if n != 1 {
		t.Fatalf("rm-dup all n=%d, want 1", n)
	}
}

func TestNormFileEntry(t *testing.T) {
	dir := t.TempDir()
	path := filepath.Join(dir, "in.vcf")
	body := []string{"chr1\t10\t.\tA\tG\t.\tPASS\t."}
	if err := os.WriteFile(path, []byte(vcfDoc(nil, body...)), 0o644); err != nil {
		t.Fatalf("write: %v", err)
	}
	var out, stderr bytes.Buffer
	n, err := NormFile(path, &out, NormOptions{}, &stderr)
	if err != nil {
		t.Fatalf("NormFile: %v", err)
	}
	if n != 1 {
		t.Fatalf("n=%d", n)
	}
}

func TestNormFileMissing(t *testing.T) {
	var out, stderr bytes.Buffer
	if _, err := NormFile("/nonexistent.vcf", &out, NormOptions{}, &stderr); err == nil {
		t.Fatalf("expected error on missing file")
	}
}

func TestSplitStrandsRoundtrip(t *testing.T) {
	cases := []string{"0/1", "1|2", "0", "./0", "0/1/2"}
	for _, c := range cases {
		got := joinStrands(splitStrands(c))
		if got != c {
			t.Errorf("roundtrip %q -> %q", c, got)
		}
	}
}

func TestAtomizePassThrough(t *testing.T) {
	// A pure SNP should pass through atomize unchanged.
	body := []string{"chr1\t10\t.\tA\tG\t.\tPASS\t."}
	input := vcfDoc(nil, body...)
	out, _, n := runNorm(t, input, NormOptions{Atomize: true})
	if n != 1 {
		t.Fatalf("n=%d", n)
	}
	if !strings.Contains(out, "chr1\t10\t.\tA\tG") {
		t.Fatalf("SNP should be unchanged:\n%s", out)
	}
}

func TestSplitMultiallelicIndelsOnly(t *testing.T) {
	// SNP multiallelic with -indels mode should NOT split.
	body := []string{"chr1\t10\t.\tA\tG,C\t.\tPASS\t."}
	input := vcfDoc(nil, body...)
	opts := NormOptions{
		Multiallelics: MultiallelicMode{Active: true, Split: true, Indels: true},
	}
	_, _, n := runNorm(t, input, opts)
	if n != 1 {
		t.Fatalf("indels-only split SNP n=%d, want 1", n)
	}
}

func TestJoinMultiallelicNoMatch(t *testing.T) {
	// Two records at different positions cannot be joined.
	body := []string{
		"chr1\t10\t.\tA\tG\t.\tPASS\t.",
		"chr1\t11\t.\tA\tC\t.\tPASS\t.",
	}
	input := vcfDoc(nil, body...)
	opts := NormOptions{
		Multiallelics: MultiallelicMode{Active: true, Split: false, Snps: true},
	}
	_, _, n := runNorm(t, input, opts)
	if n != 2 {
		t.Fatalf("join across positions should leave 2: got %d", n)
	}
}

func TestNormBCFInput(t *testing.T) {
	// Build a small BCF stream in memory using the existing writer, then
	// run Norm on it. This exercises the BCF detection branch in Norm()
	// and the readBCFAll helper.
	hdrSrc := vcfDoc(nil)
	r := vcf.NewReader(strings.NewReader(hdrSrc))
	hdr, err := r.ReadHeader()
	if err != nil {
		t.Fatalf("vcf ReadHeader: %v", err)
	}
	var buf bytes.Buffer
	bw, err := newTestBCFWriter(&buf, hdr)
	if err != nil {
		t.Fatalf("new bcf writer: %v", err)
	}
	bw.MustWriteVariants(t, []*vcf.Variant{
		{Chrom: "chr1", Pos: 10, ID: ".", Ref: "A", Alt: []string{"G"}, Qual: -1, Filter: []string{"PASS"}, Info: map[string]string{}},
		{Chrom: "chr1", Pos: 20, ID: ".", Ref: "A", Alt: []string{"T"}, Qual: -1, Filter: []string{"PASS"}, Info: map[string]string{}},
	})

	var out, stderr bytes.Buffer
	n, err := Norm(&buf, &out, NormOptions{}, &stderr)
	if err != nil {
		t.Fatalf("Norm on BCF: %v\nstderr: %s", err, stderr.String())
	}
	if n != 2 {
		t.Fatalf("bcf n=%d, want 2", n)
	}
	if !strings.Contains(out.String(), "chr1\t10") || !strings.Contains(out.String(), "chr1\t20") {
		t.Fatalf("expected both records:\n%s", out.String())
	}
}

func TestNormEmitGzipOutput(t *testing.T) {
	input := vcfDoc(nil, "chr1\t10\t.\tA\tG\t.\tPASS\t.")
	var out, stderr bytes.Buffer
	n, err := Norm(strings.NewReader(input), &out, NormOptions{
		OutputFormat:  OutputVCFGz,
		CompressLevel: -1,
	}, &stderr)
	if err != nil {
		t.Fatalf("Norm: %v", err)
	}
	if n != 1 {
		t.Fatalf("n=%d", n)
	}
	// A gzip stream begins with 0x1f 0x8b.
	b := out.Bytes()
	if len(b) < 2 || b[0] != 0x1f || b[1] != 0x8b {
		t.Fatalf("output not gzip: % x", b[:minInt(8, len(b))])
	}
}

func minInt(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func TestCheckRefMissingContig(t *testing.T) {
	ref := writeRefFasta(t, map[string]string{"chr1": "AAAA"})
	body := []string{"chrX\t1\t.\tA\tG\t.\tPASS\t."}
	input := vcfDoc(nil, body...)

	// error mode → error.
	var out, stderr bytes.Buffer
	if _, err := Norm(strings.NewReader(input), &out, NormOptions{FastaRef: ref, CheckRef: CheckRefError}, &stderr); err == nil {
		t.Fatalf("expected error for missing contig")
	}
	// skip mode → no error, record dropped.
	out.Reset()
	stderr.Reset()
	n, err := Norm(strings.NewReader(input), &out, NormOptions{FastaRef: ref, CheckRef: CheckRefSkip}, &stderr)
	if err != nil {
		t.Fatalf("skip err: %v", err)
	}
	if n != 0 {
		t.Fatalf("skip n=%d, want 0", n)
	}
}

// TestUnitMergeAlleles exercises the pure allele-merge helper (a port of
// htslib merge_alleles) without invoking any binary: identical SNPs collapse,
// distinct SNPs append, and differing-length REFs are reconciled to a common
// padded reference with the correct allele-index map.
func TestUnitMergeAlleles(t *testing.T) {
	cases := []struct {
		name    string
		a       []string // incoming line: REF, ALT...
		b       []string // accumulator:   REF, ALT...
		wantAls []string
		wantMap []int
		wantOK  bool
	}{
		{
			name:    "same-snp",
			a:       []string{"A", "C"},
			b:       []string{"A", "C"},
			wantAls: []string{"A", "C"},
			wantMap: []int{0, 1},
			wantOK:  true,
		},
		{
			name:    "distinct-snp",
			a:       []string{"A", "G"},
			b:       []string{"A", "C"},
			wantAls: []string{"A", "C", "G"},
			wantMap: []int{0, 2},
			wantOK:  true,
		},
		{
			name:    "longer-incoming-ref-pads-accumulator",
			a:       []string{"ATG", "A"}, // incoming REF longer
			b:       []string{"AT", "A"},  // accumulator REF shorter -> padded to ATG
			wantAls: []string{"ATG", "AG", "A"},
			wantMap: []int{0, 2},
			wantOK:  true,
		},
		{
			name:   "incompatible-ref-prefix",
			a:      []string{"C", "T"},
			b:      []string{"A", "G"},
			wantOK: false,
		},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			gotAls, gotMap, ok := mergeAlleles(tc.a, tc.b)
			if ok != tc.wantOK {
				t.Fatalf("ok=%v want %v", ok, tc.wantOK)
			}
			if !ok {
				return
			}
			if strings.Join(gotAls, ",") != strings.Join(tc.wantAls, ",") {
				t.Errorf("alleles=%v want %v", gotAls, tc.wantAls)
			}
			if len(gotMap) != len(tc.wantMap) {
				t.Fatalf("map len=%d want %d (%v)", len(gotMap), len(tc.wantMap), gotMap)
			}
			for i := range gotMap {
				if gotMap[i] != tc.wantMap[i] {
					t.Errorf("map[%d]=%d want %d", i, gotMap[i], tc.wantMap[i])
				}
			}
		})
	}
}

// TestUnitMaxQual checks the QUAL-maximum helper handles missing (-1) values.
func TestUnitMaxQual(t *testing.T) {
	mk := func(qs ...float64) []*vcf.Variant {
		out := make([]*vcf.Variant, len(qs))
		for i, q := range qs {
			out[i] = &vcf.Variant{Qual: q}
		}
		return out
	}
	if got := maxQual(mk(-1, -1)); got != -1 {
		t.Errorf("all-missing: got %v want -1", got)
	}
	if got := maxQual(mk(10, -1, 30, 20)); got != 30 {
		t.Errorf("max: got %v want 30", got)
	}
	if got := maxQual(mk(-1, 5)); got != 5 {
		t.Errorf("single: got %v want 5", got)
	}
}
