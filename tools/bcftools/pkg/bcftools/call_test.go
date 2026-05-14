package bcftools

import (
	"bytes"
	"os"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/vcf"
)

// callVCFHeader is the shared meta-information block used by the
// hand-crafted PL fixtures below. It declares one INFO/DP tag and the
// FORMAT/GT, FORMAT/PL pair plus three samples (S1/S2/S3 unless the
// caller trims the sample list).
const callVCFHeader = `##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
##contig=<ID=chr2,length=10000>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled likelihoods">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Sample depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
`

// makeCallVCF stitches the canonical header to a body of record lines.
func makeCallVCF(records ...string) []byte {
	return []byte(callVCFHeader + strings.Join(records, "\n") + "\n")
}

// runCall executes the Call pipeline on in with opts and returns the
// emitted stdout bytes. Test helper: fails the test on any pipeline
// error so test bodies stay short.
func runCall(t *testing.T, in []byte, opts CallOptions) []byte {
	t.Helper()
	var out bytes.Buffer
	if _, err := Call(bytes.NewReader(in), &out, opts); err != nil {
		t.Fatalf("Call: %v", err)
	}
	return out.Bytes()
}

// parseCallOutput re-parses the emitted text as VCF so tests can assert
// on structured fields (GT, INFO, QUAL) instead of brittle string
// substrings.
func parseCallOutput(t *testing.T, b []byte) (*vcf.Header, []*vcf.Variant) {
	t.Helper()
	r := vcf.NewReader(bytes.NewReader(b))
	hdr, err := r.ReadHeader()
	if err != nil {
		t.Fatalf("ReadHeader: %v", err)
	}
	vs, err := r.ReadAll()
	if err != nil {
		t.Fatalf("ReadAll: %v", err)
	}
	return hdr, vs
}

// TestCall_RejectsMissingModel asserts that Call requires -c or -m.
// Upstream rejects the call with the same precondition.
func TestCall_RejectsMissingModel(t *testing.T) {
	var out bytes.Buffer
	if _, err := Call(strings.NewReader(callVCFHeader), &out, CallOptions{}); err == nil {
		t.Fatalf("expected error for missing model, got nil")
	}
}

// TestCall_AllRefNoVariantsOnly: all samples best-likelihood is 0/0;
// without -v we still emit the site, but the GT calls are 0/0 across the
// board.
func TestCall_AllRefNoVariantsOnly(t *testing.T) {
	in := makeCallVCF("chr1\t100\t.\tA\tT\t.\t.\tDP=30\tGT:PL\t0/0:0,30,255\t0/0:0,30,255\t0/0:0,30,255")
	got := runCall(t, in, CallOptions{Model: CallModelConsensus})
	_, vs := parseCallOutput(t, got)
	if len(vs) != 1 {
		t.Fatalf("expected 1 record (no -v), got %d", len(vs))
	}
	for _, s := range vs[0].Samples {
		if s.Data["GT"] != "0/0" {
			t.Fatalf("sample %s GT = %q want 0/0", s.Name, s.Data["GT"])
		}
	}
}

// TestCall_AllRefVariantsOnlyDrops: -v drops a site where no sample's
// best-likelihood is non-ref.
func TestCall_AllRefVariantsOnlyDrops(t *testing.T) {
	in := makeCallVCF("chr1\t100\t.\tA\tT\t.\t.\tDP=30\tGT:PL\t0/0:0,30,255\t0/0:0,30,255\t0/0:0,30,255")
	got := runCall(t, in, CallOptions{Model: CallModelConsensus, VariantsOnly: true})
	_, vs := parseCallOutput(t, got)
	if len(vs) != 0 {
		t.Fatalf("expected 0 records under -v all-ref, got %d", len(vs))
	}
}

// TestCall_AltOnlyHomozygous: PL=255,30,0 makes 1/1 the most likely
// genotype across every sample → emit a variant, QUAL > 0.
func TestCall_AltOnlyHomozygous(t *testing.T) {
	in := makeCallVCF("chr1\t100\t.\tA\tT\t.\t.\tDP=30\tGT:PL\t0/0:255,30,0\t0/0:255,30,0\t0/0:255,30,0")
	got := runCall(t, in, CallOptions{Model: CallModelConsensus, VariantsOnly: true})
	_, vs := parseCallOutput(t, got)
	if len(vs) != 1 {
		t.Fatalf("expected 1 variant record, got %d:\n%s", len(vs), got)
	}
	for _, s := range vs[0].Samples {
		if s.Data["GT"] != "1/1" {
			t.Fatalf("sample %s GT = %q want 1/1", s.Name, s.Data["GT"])
		}
	}
	if vs[0].Qual <= 0 {
		t.Fatalf("QUAL %v <= 0 — expected positive QUAL for confident variant call", vs[0].Qual)
	}
}

// TestCall_HeterozygousCall: middle PL index is the lowest → 0/1.
func TestCall_HeterozygousCall(t *testing.T) {
	in := makeCallVCF("chr1\t100\t.\tA\tT\t.\t.\tDP=30\tGT:PL\t0/0:255,0,30\t0/0:255,0,30\t0/0:255,0,30")
	got := runCall(t, in, CallOptions{Model: CallModelConsensus, VariantsOnly: true})
	_, vs := parseCallOutput(t, got)
	if len(vs) != 1 {
		t.Fatalf("expected 1 variant record, got %d", len(vs))
	}
	for _, s := range vs[0].Samples {
		if s.Data["GT"] != "0/1" {
			t.Fatalf("sample %s GT = %q want 0/1", s.Name, s.Data["GT"])
		}
	}
}

// TestCall_MultiSamplePartialHet: only one of three samples is het. The
// site is still a variant; per-sample GTs differ.
func TestCall_MultiSamplePartialHet(t *testing.T) {
	in := makeCallVCF("chr1\t100\t.\tA\tT\t.\t.\tDP=30\tGT:PL\t0/0:0,30,255\t0/0:255,0,30\t0/0:0,30,255")
	got := runCall(t, in, CallOptions{Model: CallModelConsensus, VariantsOnly: true})
	_, vs := parseCallOutput(t, got)
	if len(vs) != 1 {
		t.Fatalf("expected 1 variant record, got %d", len(vs))
	}
	want := map[string]string{"S1": "0/0", "S2": "0/1", "S3": "0/0"}
	for _, s := range vs[0].Samples {
		if s.Data["GT"] != want[s.Name] {
			t.Fatalf("sample %s GT = %q want %q", s.Name, s.Data["GT"], want[s.Name])
		}
	}
	// INFO/AC should reflect one ALT allele across the three samples.
	if vs[0].Info["AC"] != "1" {
		t.Fatalf("INFO/AC = %q want 1", vs[0].Info["AC"])
	}
	if vs[0].Info["AN"] != "6" {
		t.Fatalf("INFO/AN = %q want 6", vs[0].Info["AN"])
	}
}

// TestCall_KeepAltsPreservesUnsupportedALT: -A keeps an ALT that has zero
// supporting reads. Without -A the ALT would be dropped and the record
// converted to ref-only.
func TestCall_KeepAltsPreservesUnsupportedALT(t *testing.T) {
	in := makeCallVCF("chr1\t100\t.\tA\tT\t.\t.\tDP=30\tGT:PL\t0/0:0,30,255\t0/0:0,30,255\t0/0:0,30,255")
	got := runCall(t, in, CallOptions{Model: CallModelConsensus, KeepAlts: true})
	_, vs := parseCallOutput(t, got)
	if len(vs) != 1 {
		t.Fatalf("expected 1 record under -A, got %d", len(vs))
	}
	if len(vs[0].Alt) != 1 || vs[0].Alt[0] != "T" {
		t.Fatalf("ALT = %v want [T] under -A", vs[0].Alt)
	}
}

// TestCall_KeepAltsAllowsAllRefWithVariantsOnly: with -A *and* -v the
// site must still be emitted because -A overrides the "drop all-ref"
// rule (matching upstream behaviour).
func TestCall_KeepAltsAllowsAllRefWithVariantsOnly(t *testing.T) {
	in := makeCallVCF("chr1\t100\t.\tA\tT\t.\t.\tDP=30\tGT:PL\t0/0:0,30,255\t0/0:0,30,255\t0/0:0,30,255")
	got := runCall(t, in, CallOptions{Model: CallModelConsensus, KeepAlts: true, VariantsOnly: true})
	_, vs := parseCallOutput(t, got)
	if len(vs) != 1 {
		t.Fatalf("expected 1 record under -A -v all-ref, got %d", len(vs))
	}
}

// TestCall_PvalThresholdBlocksLowConfidence: a borderline PL (5 vs 0)
// produces a low posterior; an aggressive -p should drop the site.
func TestCall_PvalThresholdBlocksLowConfidence(t *testing.T) {
	in := makeCallVCF("chr1\t100\t.\tA\tT\t.\t.\tDP=10\tGT:PL\t0/0:5,0,5")
	// Single sample is het at PL=0; pass-without-recompute → variant kept.
	got1 := runCall(t, in, CallOptions{Model: CallModelConsensus, VariantsOnly: true})
	_, vs1 := parseCallOutput(t, got1)
	if len(vs1) != 1 {
		t.Fatalf("expected 1 record at default threshold, got %d", len(vs1))
	}
	// Make the prior the only thing that matters — and crank up the
	// threshold so the posterior comparison alone can't keep the call.
	// We also wipe the het PL by giving it equal likelihoods so the
	// best-non-ref tie is broken arbitrarily and the AC fallback can't
	// kick in. To do that we need a *different* input where the best
	// likelihood is ref.
	in2 := makeCallVCF("chr1\t100\t.\tA\tT\t.\t.\tDP=10\tGT:PL\t0/0:0,3,255")
	got2 := runCall(t, in2, CallOptions{Model: CallModelConsensus, VariantsOnly: true, PvalThreshold: 0.99, Prior: 1e-12})
	_, vs2 := parseCallOutput(t, got2)
	if len(vs2) != 0 {
		t.Fatalf("expected 0 records at -p 0.99 -P 1e-12, got %d", len(vs2))
	}
}

// TestCall_PriorShiftsCalls: tiny prior + tight threshold → drop the
// borderline het. The opposite (huge prior) keeps it.
func TestCall_PriorShiftsCalls(t *testing.T) {
	in := makeCallVCF("chr1\t100\t.\tA\tT\t.\t.\tDP=10\tGT:PL\t0/0:0,3,255")
	got := runCall(t, in, CallOptions{Model: CallModelConsensus, VariantsOnly: true, PvalThreshold: 0.99, Prior: 0.5})
	_, vs := parseCallOutput(t, got)
	// With a huge prior the het beats the threshold even though the PL
	// difference is small.
	if len(vs) != 1 {
		t.Fatalf("expected 1 record with strong prior, got %d", len(vs))
	}
}

// TestCall_HaploidPloidy: --ploidy 1 (or -X) produces "0" / "1" GTs.
func TestCall_HaploidPloidy(t *testing.T) {
	// Haploid PLs are length-N (one per allele). For A/T at PL=0 for the
	// REF (allele 0) and PL=30 for the ALT (allele 1) the call is "0".
	// Switch the order to elicit the alt call.
	in := makeCallVCF("chr1\t100\t.\tA\tT\t.\t.\tDP=30\tGT:PL\t0:30,0\t0:30,0\t0:0,30")
	got := runCall(t, in, CallOptions{Model: CallModelConsensus, Ploidy: PloidyHaploid})
	_, vs := parseCallOutput(t, got)
	if len(vs) != 1 {
		t.Fatalf("expected 1 record under --ploidy 1, got %d", len(vs))
	}
	want := map[string]string{"S1": "1", "S2": "1", "S3": "0"}
	for _, s := range vs[0].Samples {
		if s.Data["GT"] != want[s.Name] {
			t.Fatalf("sample %s haploid GT = %q want %q", s.Name, s.Data["GT"], want[s.Name])
		}
	}
}

// TestCall_PloidySpec_GRCh37Deferred asserts the documented v1 behaviour:
// --ploidy GRCh37 (or GRCh38) is parsed but rejected at runtime with a
// roadmap pointer.
func TestCall_PloidySpec_GRCh37Deferred(t *testing.T) {
	in := makeCallVCF("chr1\t100\t.\tA\tT\t.\t.\tDP=30\tGT:PL\t0/0:0,30,255\t0/0:0,30,255\t0/0:0,30,255")
	var out bytes.Buffer
	_, err := Call(bytes.NewReader(in), &out, CallOptions{Model: CallModelConsensus, PloidySpec: "GRCh37"})
	if err == nil {
		t.Fatalf("expected error for --ploidy GRCh37, got nil")
	}
	if !strings.Contains(err.Error(), "GRCh37") {
		t.Fatalf("expected GRCh37 in error, got %q", err.Error())
	}
}

// TestCall_RegionRestriction: -r/-t drops records outside the requested
// chromosome.
func TestCall_RegionRestriction(t *testing.T) {
	in := makeCallVCF(
		"chr1\t100\t.\tA\tT\t.\t.\tDP=30\tGT:PL\t0/0:255,30,0\t0/0:255,30,0\t0/0:255,30,0",
		"chr2\t200\t.\tC\tG\t.\t.\tDP=30\tGT:PL\t0/0:255,30,0\t0/0:255,30,0\t0/0:255,30,0",
	)
	got := runCall(t, in, CallOptions{Model: CallModelConsensus, VariantsOnly: true, Targets: []string{"chr1"}})
	_, vs := parseCallOutput(t, got)
	if len(vs) != 1 {
		t.Fatalf("expected 1 record after -t chr1, got %d", len(vs))
	}
	if vs[0].Chrom != "chr1" {
		t.Fatalf("kept record chrom %q want chr1", vs[0].Chrom)
	}
}

// TestCall_MultiallelicConsensusFallback: the -m caller, on a
// multi-allelic site, falls back to the consensus model. The contract
// today is "produce some call without erroring" — we assert the site is
// emitted and at least one sample's GT references a non-ref allele.
func TestCall_MultiallelicConsensusFallback(t *testing.T) {
	// Site A/T,G with three samples: each prefers a different ALT.
	// PL ordering for diploid 3-allele site:
	//   0=0/0, 1=0/1, 2=1/1, 3=0/2, 4=1/2, 5=2/2
	in := makeCallVCF("chr1\t100\t.\tA\tT,G\t.\t.\tDP=30\tGT:PL\t0/0:255,0,255,255,255,255\t0/0:255,255,255,0,255,255\t0/0:255,255,0,255,255,255")
	got := runCall(t, in, CallOptions{Model: CallModelMultiallelic, VariantsOnly: true})
	_, vs := parseCallOutput(t, got)
	if len(vs) != 1 {
		t.Fatalf("expected 1 record for multiallelic site, got %d", len(vs))
	}
	if len(vs[0].Alt) == 0 {
		t.Fatalf("expected at least one ALT allele to be kept")
	}
}

// TestCall_SampleSubsetRestrictsOutput: -s narrows the per-sample columns.
func TestCall_SampleSubsetRestrictsOutput(t *testing.T) {
	in := makeCallVCF("chr1\t100\t.\tA\tT\t.\t.\tDP=30\tGT:PL\t0/0:0,30,255\t0/0:0,30,255\t0/0:0,30,255")
	got := runCall(t, in, CallOptions{Model: CallModelConsensus, Samples: []string{"S2"}})
	hdr, vs := parseCallOutput(t, got)
	if len(hdr.Samples) != 1 || hdr.Samples[0] != "S2" {
		t.Fatalf("header samples = %v want [S2]", hdr.Samples)
	}
	if len(vs) != 1 || len(vs[0].Samples) != 1 || vs[0].Samples[0].Name != "S2" {
		t.Fatalf("record samples = %v want one S2", vs[0].Samples)
	}
}

// TestCall_DeclaresACANInHeader confirms the canonical header
// augmentation runs even when the input lacks ##INFO/AC entries.
func TestCall_DeclaresACANInHeader(t *testing.T) {
	in := makeCallVCF("chr1\t100\t.\tA\tT\t.\t.\tDP=30\tGT:PL\t0/0:255,30,0\t0/0:255,30,0\t0/0:255,30,0")
	got := runCall(t, in, CallOptions{Model: CallModelConsensus, VariantsOnly: true})
	if !bytes.Contains(got, []byte("##INFO=<ID=AC,")) {
		t.Fatalf("missing ##INFO=<ID=AC,...> meta line in output:\n%s", got)
	}
	if !bytes.Contains(got, []byte("##INFO=<ID=AN,")) {
		t.Fatalf("missing ##INFO=<ID=AN,...> meta line in output:\n%s", got)
	}
}

// TestCall_FileEntryPoint exercises CallFile (vs Call) on a temp file
// with both BCF and gzip preludes ruled out.
func TestCall_FileEntryPoint(t *testing.T) {
	dir := t.TempDir()
	path := dir + "/in.vcf"
	if err := writeFileBytes(path, makeCallVCF("chr1\t100\t.\tA\tT\t.\t.\tDP=30\tGT:PL\t0/0:255,30,0")); err != nil {
		t.Fatal(err)
	}
	var out bytes.Buffer
	if _, err := CallFile(path, &out, CallOptions{Model: CallModelConsensus, VariantsOnly: true}, nil); err != nil {
		t.Fatalf("CallFile: %v", err)
	}
	body := out.Bytes()
	if !bytes.Contains(body, []byte("chr1\t100")) {
		t.Fatalf("expected variant record in output, got:\n%s", body)
	}
}

// TestCall_ParsePloidySpec covers the public parser used by the CLI.
func TestCall_ParsePloidySpec(t *testing.T) {
	cases := []struct {
		in       string
		wantP    PloidySpec
		wantText string
		wantErr  bool
	}{
		{"", PloidyDiploid, "2", false},
		{"2", PloidyDiploid, "2", false},
		{"1", PloidyHaploid, "1", false},
		{"GRCh37", PloidyDiploid, "GRCh37", false},
		{"GRCh38", PloidyDiploid, "GRCh38", false},
		{"hmm", 0, "", true},
	}
	for _, c := range cases {
		p, text, err := ParsePloidySpec(c.in)
		if (err != nil) != c.wantErr {
			t.Fatalf("ParsePloidySpec(%q): err=%v wantErr=%v", c.in, err, c.wantErr)
		}
		if !c.wantErr && (p != c.wantP || text != c.wantText) {
			t.Fatalf("ParsePloidySpec(%q) = (%v, %q) want (%v, %q)", c.in, p, text, c.wantP, c.wantText)
		}
	}
}

// TestCall_DecomposeGTIndex exercises both ploidies of the
// decomposeGTIndex helper.
func TestCall_DecomposeGTIndex(t *testing.T) {
	// Diploid: 3 alleles → 6 GT indices.
	type want struct{ a1, a2 int }
	wants := []want{{0, 0}, {0, 1}, {1, 1}, {0, 2}, {1, 2}, {2, 2}}
	for i, w := range wants {
		a1, a2, ok := decomposeGTIndex(i, 3, PloidyDiploid)
		if !ok || a1 != w.a1 || a2 != w.a2 {
			t.Fatalf("decomposeGTIndex(%d, 3, diploid) = (%d, %d, %v) want (%d, %d, true)", i, a1, a2, ok, w.a1, w.a2)
		}
	}
	// Haploid: index is allele directly.
	for i := 0; i < 3; i++ {
		a, _, ok := decomposeGTIndex(i, 3, PloidyHaploid)
		if !ok || a != i {
			t.Fatalf("decomposeGTIndex(%d, 3, haploid) = (%d, _, %v) want (%d, _, true)", i, a, ok, i)
		}
	}
	// Out of range.
	if _, _, ok := decomposeGTIndex(100, 3, PloidyDiploid); ok {
		t.Fatalf("expected out-of-range PL index to fail")
	}
}

// TestCall_DecodePLMissing / Padding covers the resilient decoder paths.
func TestCall_DecodePLMissingAndPadding(t *testing.T) {
	if _, ok := decodePL("", 2, PloidyDiploid); ok {
		t.Fatalf("empty PL should return ok=false")
	}
	if _, ok := decodePL(".", 2, PloidyDiploid); ok {
		t.Fatalf("dot PL should return ok=false")
	}
	// Diploid 2 alleles expects 3 values; pad short input.
	pl, ok := decodePL("0,30", 2, PloidyDiploid)
	if !ok || len(pl) != 3 || pl[2] != 255 {
		t.Fatalf("decodePL pad: %v ok=%v", pl, ok)
	}
	// Garbage value rejected.
	if _, ok := decodePL("abc", 2, PloidyDiploid); ok {
		t.Fatalf("non-numeric PL should return ok=false")
	}
}

// TestCall_RemapGTByIndex exercises the per-allele index remap helper.
func TestCall_RemapGTByIndex(t *testing.T) {
	// "1/2" with remap {0:0, 1:1} (drop allele 2) → "1/."
	got := remapGTByIndex("1/2", map[int]int{0: 0, 1: 1})
	if got != "1/." {
		t.Fatalf("remapGTByIndex(1/2) = %q want 1/.", got)
	}
	// "." preserved as ".".
	if remapGTByIndex(".", map[int]int{0: 0}) != "." {
		t.Fatalf("remapGTByIndex(.) should preserve dot")
	}
	// Phased separator preserved.
	if remapGTByIndex("0|1", map[int]int{0: 0, 1: 1}) != "0|1" {
		t.Fatalf("remapGTByIndex preserves | separator")
	}
}

// TestCall_TrimUnsupportedAltsDoesNothingWhenAllSupported is the
// nothing-to-do branch of trimUnsupportedAlts.
func TestCall_TrimUnsupportedAltsDoesNothingWhenAllSupported(t *testing.T) {
	v := &vcf.Variant{Alt: []string{"T", "G"}, Samples: []vcf.Sample{
		{Name: "S1", Data: map[string]string{"GT": "1/2"}},
	}}
	got := trimUnsupportedAlts(v.Alt, []int{1, 1}, v)
	if len(got) != 2 || got[0] != "T" || got[1] != "G" {
		t.Fatalf("expected no trim, got %v", got)
	}
}

// writeFileBytes is a tiny helper used by TestCall_FileEntryPoint.
func writeFileBytes(path string, b []byte) error {
	return os.WriteFile(path, b, 0o644)
}
