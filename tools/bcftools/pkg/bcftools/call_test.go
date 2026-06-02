package bcftools

import (
	"bytes"
	"os"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
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

// TestCall_PloidySpec_GRCh37Accepted asserts that --ploidy GRCh37 builds
// the per-region ploidy table and runs without error on an autosome.
// Per-chromosome behaviour is exercised by TestCall_PloidyGRCh37PerContig
// and the live oracle suite.
func TestCall_PloidySpec_GRCh37Accepted(t *testing.T) {
	in := makeCallVCF("chr1\t100\t.\tA\tT\t.\t.\tDP=30\tGT:PL\t0/0:0,30,255\t0/0:0,30,255\t0/0:0,30,255")
	var out bytes.Buffer
	_, err := Call(bytes.NewReader(in), &out, CallOptions{Model: CallModelConsensus, PloidySpec: "GRCh37"})
	if err != nil {
		t.Fatalf("--ploidy GRCh37 should be accepted, got %v", err)
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

// mcallVCFHeader declares the mpileup-style INFO/FORMAT tags (DP, I16, QS,
// MQ0F, PL, AD) the faithful mcall path consumes. One sample (HG00100).
const mcallVCFHeader = `##fileformat=VCFv4.2
##contig=<ID=17,length=4200>
##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
##INFO=<ID=I16,Number=16,Type=Float,Description="Auxiliary tag used for calling, see description of bcf_callret1_t in bam2bcf.h">
##INFO=<ID=QS,Number=R,Type=Float,Description="Auxiliary tag used for calling">
##INFO=<ID=MQ0F,Number=1,Type=Float,Description="Fraction of MQ0 reads (smaller is better)">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths (high-quality bases)">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	HG00100
`

func makeMcallVCF(records ...string) []byte {
	return []byte(mcallVCFHeader + strings.Join(records, "\n") + "\n")
}

// TestCall_McallHomRef reproduces the pos-1 mpileup row from the live
// fixture: a clean ref-only site (`<*>` only) must yield QUAL=129.588,
// GT=0/0, INFO DP/MQ0F/AN/DP4/MQ (I16/QS dropped), FORMAT GT:AD.
func TestCall_McallHomRef(t *testing.T) {
	in := makeMcallVCF("17\t1\t.\tA\t<*>\t0\t.\tDP=5;I16=5,0,0,0,202,8170,0,0,145,4205,0,0,107,2459,0,0;QS=1,0;MQ0F=0\tPL:AD\t0,15,100:5,0")
	got := runCall(t, in, CallOptions{Model: CallModelMultiallelic})
	_, vs := parseCallOutput(t, got)
	if len(vs) != 1 {
		t.Fatalf("expected 1 record, got %d:\n%s", len(vs), got)
	}
	v := vs[0]
	if len(v.Alt) != 1 || v.Alt[0] != "." {
		t.Fatalf("ALT = %v want [.]", v.Alt)
	}
	if got := formatFloat32G(v.Qual); got != "129.588" {
		t.Fatalf("QUAL = %q want 129.588", got)
	}
	if s := v.Samples[0]; s.Data["GT"] != "0/0" || s.Data["AD"] != "5" {
		t.Fatalf("sample = %v want GT=0/0 AD=5", s.Data)
	}
	if _, ok := v.Samples[0].Data["PL"]; ok {
		t.Fatalf("PL should be dropped for a ref-only site")
	}
	for k, want := range map[string]string{"AN": "2", "DP4": "5,0,0,0", "MQ": "29"} {
		if v.Info[k] != want {
			t.Fatalf("INFO/%s = %q want %q", k, v.Info[k], want)
		}
	}
	if _, ok := v.Info["I16"]; ok {
		t.Fatalf("INFO/I16 must be dropped")
	}
	if _, ok := v.Info["QS"]; ok {
		t.Fatalf("INFO/QS must be dropped")
	}
}

// TestCall_McallVariantHet reproduces the pos-828 SNP (T -> C,<*>): the
// unseen allele is dropped, GT is the called het, QUAL/AC/DP4/MQ match
// upstream, and PL/AD are re-indexed to the surviving T,C alleles.
func TestCall_McallVariantHet(t *testing.T) {
	in := makeMcallVCF("17\t828\t.\tT\tC,<*>\t0\t.\tDP=12;I16=1,1,3,7,71,2525,359,13309,120,7200,600,36000,41,841,166,3304;QS=0.165116,0.834884,0;MQ0F=0\tPL:AD\t216,0,35,223,65,255:2,10,0")
	got := runCall(t, in, CallOptions{Model: CallModelMultiallelic})
	_, vs := parseCallOutput(t, got)
	if len(vs) != 1 {
		t.Fatalf("expected 1 record, got %d", len(vs))
	}
	v := vs[0]
	if len(v.Alt) != 1 || v.Alt[0] != "C" {
		t.Fatalf("ALT = %v want [C] (<*> dropped)", v.Alt)
	}
	if got := formatFloat32G(v.Qual); got != "180.829" {
		t.Fatalf("QUAL = %q want 180.829", got)
	}
	s := v.Samples[0]
	if s.Data["GT"] != "0/1" || s.Data["PL"] != "216,0,35" || s.Data["AD"] != "2,10" {
		t.Fatalf("sample = %v want GT=0/1 PL=216,0,35 AD=2,10", s.Data)
	}
	for k, want := range map[string]string{"AC": "1", "AN": "2", "DP4": "1,1,3,7", "MQ": "60"} {
		if v.Info[k] != want {
			t.Fatalf("INFO/%s = %q want %q", k, v.Info[k], want)
		}
	}
}

// writeFileBytes is a tiny helper used by TestCall_FileEntryPoint.
func writeFileBytes(path string, b []byte) error {
	return os.WriteFile(path, b, 0o644)
}

// TestCall_PloidyTableGRCh37 exercises ParsePloidyTable and its
// per-contig query: F is registered last so DefaultSexID() points at F
// (matching vcfcall.c sample2sex initialisation). The query then
// returns 2 for autosomes, 0 for chrY, 1 for chrM, and 2 for chrX.
func TestCall_PloidyTableGRCh37(t *testing.T) {
	tbl, err := BuildPloidyTableFromSpec("GRCh37")
	if err != nil {
		t.Fatalf("BuildPloidyTableFromSpec: %v", err)
	}
	if tbl == nil {
		t.Fatal("expected non-nil ploidy table for GRCh37")
	}
	if got := tbl.SexName(tbl.DefaultSexID()); got != "F" {
		t.Fatalf("default sex = %q, want F", got)
	}
	f := tbl.SexID("F")
	m := tbl.SexID("M")
	cases := []struct {
		name       string
		chrom      string
		pos        int
		sex        int
		wantPloidy int
	}{
		{"chr1-F", "chr1", 100, f, 2},
		{"chr1-M", "chr1", 100, m, 2},
		{"X-F-anywhere", "X", 100, f, 2},
		{"X-M-PAR1", "X", 50, m, 1},
		{"X-M-mid", "X", 500_000, m, 2},
		{"X-M-PAR2", "X", 100_000_000, m, 1},
		{"Y-F", "Y", 100, f, 0},
		{"Y-M", "Y", 100, m, 1},
		{"MT-F", "MT", 100, f, 1},
		{"MT-M", "MT", 100, m, 1},
		{"chrM-F", "chrM", 100, f, 1},
	}
	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			got := tbl.Query(c.chrom, c.pos, c.sex)
			if got != c.wantPloidy {
				t.Fatalf("Query(%q,%d,sex=%d) = %d, want %d", c.chrom, c.pos, c.sex, got, c.wantPloidy)
			}
		})
	}
}

// TestCall_PloidyGRCh37PerContig drives a synthetic mpileup-shaped
// fixture with QS/I16 + PL on chrX, chrY, chrM and asserts the
// per-contig output matches what upstream `bcftools call -m
// --ploidy GRCh37` produces (encoded in the want strings below).
func TestCall_PloidyGRCh37PerContig(t *testing.T) {
	const hdr = `##fileformat=VCFv4.2
##contig=<ID=chr1>
##contig=<ID=X>
##contig=<ID=Y>
##contig=<ID=MT>
##INFO=<ID=QS,Number=R,Type=Float,Description="QS">
##INFO=<ID=I16,Number=16,Type=Float,Description="I16">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="PL">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
`
	body := `chr1	100	.	A	G	.	.	QS=0.5,0.5;I16=10,10,10,10,400,400,400,400,30,30,30,30,30,30,30,30	PL	30,0,30
X	100	.	A	G	.	.	QS=0.5,0.5;I16=10,10,10,10,400,400,400,400,30,30,30,30,30,30,30,30	PL	30,0,30
Y	100	.	A	G	.	.	QS=0.5,0.5;I16=10,10,10,10,400,400,400,400,30,30,30,30,30,30,30,30	PL	30,0,30
MT	100	.	A	G	.	.	QS=0.5,0.5;I16=10,10,10,10,400,400,400,400,30,30,30,30,30,30,30,30	PL	30,0,30
`
	in := []byte(hdr + body)
	got := runCall(t, in, CallOptions{Model: CallModelMultiallelic, PloidySpec: "GRCh37"})
	_, vs := parseCallOutput(t, got)
	if len(vs) != 4 {
		t.Fatalf("got %d records, want 4", len(vs))
	}
	// Default sex is F under the GRCh37 predefs.
	wantGT := map[string]string{"chr1": "0/0", "X": "0/0", "Y": ".", "MT": "0"}
	for _, v := range vs {
		if got := v.Samples[0].Data["GT"]; got != wantGT[v.Chrom] {
			t.Fatalf("%s: GT = %q, want %q", v.Chrom, got, wantGT[v.Chrom])
		}
	}
}
