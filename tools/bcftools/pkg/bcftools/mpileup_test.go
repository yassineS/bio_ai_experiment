package bcftools

import (
	"bytes"
	"os"
	"os/exec"
	"path/filepath"
	"strconv"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/errmod"
	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// samRecordStub helps build minimal sam.Records for filter tests.
type samRecordStub struct {
	Flag  uint16
	MapQ  uint8
	RName string
	Pos   int32
}

func (s *samRecordStub) toRecord() *sam.Record {
	return &sam.Record{Flag: s.Flag, MapQ: s.MapQ, RName: s.RName, Pos: s.Pos}
}

// mkPile builds a pileup column of n reads carrying the given uppercase
// base at the given quality, on the forward strand, with MAPQ 60. The
// neighbour qualities are set equal to qual so the delta_baseQ cap is a
// no-op.
func mkPile(n int, base byte, qual uint8) []pileupBase {
	out := make([]pileupBase, 0, n)
	for i := 0; i < n; i++ {
		out = append(out, pileupBase{
			base4:   seqNt16Int[baseToNt16(base)],
			rawQual: qual,
			prevQ:   int(qual),
			nextQ:   int(qual),
			mapq:    60,
			qlen:    10,
			qpos:    5,
		})
	}
	return out
}

// glfPL runs glfgen + combine on a single sample and returns the per-
// sample PL grid plus the resolved alleles.
func glfPL(t *testing.T, pile []pileupBase, ref byte) (bcfCall, []int) {
	t.Helper()
	em := errmod.Init(1.0 - mpileupTheta)
	ref4 := seqNt16Int[baseToNt16(ref)]
	var cr bcfCallret
	bcfCallGlfgen(pile, ref4, MpileupOptions{MinBQ: 1, MaxBQ: 60, DeltaBQ: 30}, em, &cr)
	call := bcfCallCombine([]bcfCallret{cr}, ref4)
	return call, call.pl[0]
}

// TestGlfgenHomRef: an all-REF pile yields PL[0/0]=0 and rising PL for
// the het and hom-alt genotypes against the appended <*> allele.
func TestGlfgenHomRef(t *testing.T) {
	call, pl := glfPL(t, mkPile(10, 'A', 30), 'A')
	// REF + <*> => 2 alleles => 3 PL values.
	if call.nAlleles != 2 {
		t.Fatalf("nAlleles = %d, want 2 (REF + <*>)", call.nAlleles)
	}
	if call.unseen != 1 {
		t.Errorf("unseen index = %d, want 1", call.unseen)
	}
	if len(pl) != 3 {
		t.Fatalf("PL grid len = %d, want 3", len(pl))
	}
	if pl[0] != 0 {
		t.Errorf("hom-ref PL[0/0] = %d, want 0", pl[0])
	}
	if pl[2] <= pl[1] || pl[1] <= 0 {
		t.Errorf("hom-ref PL = %v, want rising 0 < pl[1] < pl[2]", pl)
	}
}

// TestGlfgenHetSite: an even split of REF and ALT bases makes the
// heterozygous genotype the most likely (PL of the 0/1 cell is 0).
func TestGlfgenHetSite(t *testing.T) {
	pile := append(mkPile(5, 'A', 30), mkPile(5, 'C', 30)...)
	call, pl := glfPL(t, pile, 'A')
	// REF=A, ALT=C, <*> => 3 alleles => 6 PL values.
	if call.nAlleles != 3 {
		t.Fatalf("nAlleles = %d, want 3 (REF, C, <*>)", call.nAlleles)
	}
	if len(pl) != 6 {
		t.Fatalf("PL grid len = %d, want 6", len(pl))
	}
	// PL ordering for alleles (0=REF,1=C,2=<*>): index 0=0/0, 1=0/1,
	// 2=1/1, 3=0/2, 4=1/2, 5=2/2.
	if pl[1] != 0 {
		t.Errorf("het PL[0/1] = %d, want 0 (het most likely)", pl[1])
	}
	if pl[0] <= 0 || pl[2] <= 0 {
		t.Errorf("het: hom-ref/hom-alt PL should be > 0, got %v", pl)
	}
}

// TestGlfgenEmptyPile: a zero-coverage column returns 0 bases.
func TestGlfgenEmptyPile(t *testing.T) {
	em := errmod.Init(1.0 - mpileupTheta)
	var cr bcfCallret
	n := bcfCallGlfgen(nil, 0, MpileupOptions{MinBQ: 1, MaxBQ: 60, DeltaBQ: 30}, em, &cr)
	if n != 0 {
		t.Errorf("empty pile glfgen returned %d, want 0", n)
	}
}

// TestCombineAlleleOrdering: ALT alleles are ordered by descending
// coverage-normalised QS sum, with <*> always last.
func TestCombineAlleleOrdering(t *testing.T) {
	// REF=A. Pile: 2xA, 6xC, 4xG. Expect order REF(A), C, G, <*>.
	pile := mkPile(2, 'A', 30)
	pile = append(pile, mkPile(6, 'C', 30)...)
	pile = append(pile, mkPile(4, 'G', 30)...)
	em := errmod.Init(1.0 - mpileupTheta)
	var cr bcfCallret
	bcfCallGlfgen(pile, seqNt16Int[baseToNt16('A')], MpileupOptions{MinBQ: 1, MaxBQ: 60, DeltaBQ: 30}, em, &cr)
	call := bcfCallCombine([]bcfCallret{cr}, seqNt16Int[baseToNt16('A')])
	if call.nAlleles != 4 {
		t.Fatalf("nAlleles = %d, want 4", call.nAlleles)
	}
	// alleles[0]=A(0), [1]=C(1), [2]=G(2), [3]=<*>.
	if call.alleles[0] != 0 || call.alleles[1] != 1 || call.alleles[2] != 2 {
		t.Errorf("allele order = %v, want [A,C,G,...]", call.alleles[:3])
	}
	if call.unseen != 3 {
		t.Errorf("unseen index = %d, want 3", call.unseen)
	}
	// PL grid is upper triangle of a 4x4 => 10 values.
	if len(call.pl[0]) != 10 {
		t.Errorf("PL grid len = %d, want 10", len(call.pl[0]))
	}
}

// TestBcfCall2bcfRecord checks the emitted vcf.Variant: REF, ALT incl.
// <*>, QUAL=0, the INFO ordering (DP/I16/QS first, the bias tags, then
// MQ0F last), and FORMAT/PL.
func TestBcfCall2bcfRecord(t *testing.T) {
	pile := append(mkPile(5, 'A', 30), mkPile(3, 'C', 30)...)
	em := errmod.Init(1.0 - mpileupTheta)
	var cr bcfCallret
	bcfCallGlfgen(pile, seqNt16Int[baseToNt16('A')], MpileupOptions{MinBQ: 1, MaxBQ: 60, DeltaBQ: 30}, em, &cr)
	call := bcfCallCombine([]bcfCallret{cr}, seqNt16Int[baseToNt16('A')])
	v := bcfCall2bcf("chr1", 42, 'A', &call, 0)
	if v.Chrom != "chr1" || v.Pos != 42 || v.Ref != "A" {
		t.Errorf("record locus = %s:%d %s, want chr1:42 A", v.Chrom, v.Pos, v.Ref)
	}
	if v.Qual != 0 {
		t.Errorf("QUAL = %v, want 0", v.Qual)
	}
	if len(v.Alt) != 2 || v.Alt[0] != "C" || v.Alt[1] != "<*>" {
		t.Errorf("ALT = %v, want [C <*>]", v.Alt)
	}
	// INFO begins with DP/I16/QS and ends with MQ0F; the bias tags
	// (present here because the site has a real ALT) sit in between.
	wantPrefix := []string{"DP", "I16", "QS"}
	for i, k := range wantPrefix {
		if i >= len(v.InfoOrder) || v.InfoOrder[i] != k {
			t.Errorf("InfoOrder = %v, want prefix %v", v.InfoOrder, wantPrefix)
			break
		}
	}
	if n := len(v.InfoOrder); n == 0 || v.InfoOrder[n-1] != "MQ0F" {
		t.Errorf("InfoOrder = %v, want MQ0F last", v.InfoOrder)
	}
	if v.Info["DP"] != "8" {
		t.Errorf("INFO/DP = %q, want 8", v.Info["DP"])
	}
	i16 := strings.Split(v.Info["I16"], ",")
	if len(i16) != 16 {
		t.Errorf("INFO/I16 has %d values, want 16", len(i16))
	}
	// QS has one value per allele (REF, C, <*>) => 3.
	if qs := strings.Split(v.Info["QS"], ","); len(qs) != 3 {
		t.Errorf("INFO/QS = %q, want 3 values", v.Info["QS"])
	}
	if len(v.Format) != 1 || v.Format[0] != "PL" {
		t.Errorf("FORMAT = %v, want [PL]", v.Format)
	}
	pl := strings.Split(v.Samples[0].Data["PL"], ",")
	if len(pl) != 6 { // 3 alleles => 6 PL values.
		t.Errorf("FORMAT/PL = %q, want 6 values", v.Samples[0].Data["PL"])
	}
}

// TestGlfgenDeltaBQCap verifies the neighbour-quality cap: a high-Q base
// flanked by low-Q neighbours is downgraded to neighbour_q + DeltaBQ.
func TestGlfgenDeltaBQCap(t *testing.T) {
	// One base Q=60 with neighbour qualities 2. With DeltaBQ=30 the
	// effective quality is min(60, 2+30) = 32.
	pile := []pileupBase{{
		base4: seqNt16Int[baseToNt16('A')], rawQual: 60,
		prevQ: 2, nextQ: 2, mapq: 60, qlen: 10, qpos: 5,
	}}
	em := errmod.Init(1.0 - mpileupTheta)
	var cr bcfCallret
	bcfCallGlfgen(pile, seqNt16Int[baseToNt16('A')], MpileupOptions{MinBQ: 1, MaxBQ: 60, DeltaBQ: 30}, em, &cr)
	// The capped quality 32 lands in QS[0] (the A allele).
	if cr.qs[0] != 32 {
		t.Errorf("delta_baseQ cap: QS[A] = %v, want 32 (min(60,2+30))", cr.qs[0])
	}
}

// TestMpileupGoldenStructure runs mpileup against the upstream
// mpileup.1.bam fixture and validates the output structurally against
// the upstream golden mpileup/mpileup.3.out.
//
// Full byte-for-byte parity is NOT achievable in slice 2: the upstream
// golden was produced with `-B --ff 0x14` and bakes in (a) the --ff
// flag-filter (out of scope here) and, for ALT sites, (b) the bias
// annotations (slice 4). We therefore assert the structure: one record
// per covered position, the <*> allele always present, QUAL=0, the
// INFO tag set, and a PL grid whose length matches the allele count.
func TestMpileupGoldenStructure(t *testing.T) {
	bam := referenceFixture(t, "mpileup/mpileup.1.bam")
	ref := referenceFixture(t, "mpileup/mpileup.ref.fa")

	var buf bytes.Buffer
	opts := MpileupOptions{
		Inputs:   []string{bam},
		FastaRef: ref,
		Regions:  []string{"17:1050-1060"},
		NoBAQ:    true,
	}
	if err := MpileupFile(opts, &buf); err != nil {
		t.Fatalf("MpileupFile: %v", err)
	}
	var data []string
	for _, ln := range strings.Split(buf.String(), "\n") {
		if ln == "" || strings.HasPrefix(ln, "#") {
			continue
		}
		data = append(data, ln)
	}
	// 17:1050-1060 inclusive = 11 covered positions.
	if len(data) != 11 {
		t.Fatalf("got %d records for 17:1050-1060, want 11", len(data))
	}
	for _, ln := range data {
		f := strings.Split(ln, "\t")
		if len(f) != 10 {
			t.Fatalf("record %q has %d fields, want 10", ln, len(f))
		}
		if f[5] != "0" {
			t.Errorf("%s:%s QUAL=%q, want 0", f[0], f[1], f[5])
		}
		alt := strings.Split(f[4], ",")
		if alt[len(alt)-1] != "<*>" {
			t.Errorf("%s:%s ALT=%q, want trailing <*>", f[0], f[1], f[4])
		}
		for _, tag := range []string{"DP=", "I16=", "QS=", "MQ0F="} {
			if !strings.Contains(f[7], tag) {
				t.Errorf("%s:%s INFO=%q missing %s", f[0], f[1], f[7], tag)
			}
		}
		// PL grid length must be n_alleles*(n_alleles+1)/2. The default
		// FORMAT is PL:AD; PL is the first subfield.
		nAll := len(alt) + 1 // ALT count + REF
		wantPL := nAll * (nAll + 1) / 2
		plField := strings.Split(f[9], ":")[0]
		if pl := strings.Split(plField, ","); len(pl) != wantPL {
			t.Errorf("%s:%s PL has %d values, want %d for %d alleles",
				f[0], f[1], len(pl), wantPL, nAll)
		}
	}
}

// mpileupStarRecords runs mpileup over the given region of one upstream
// fixture BAM and returns the data records whose ALT is exactly the
// `<*>` unseen allele (i.e. no called ALT, hence no slice-4 bias INFO
// tags), each as a tab-joined string.
func mpileupStarRecords(t *testing.T, bam, ref, region string, noBAQ bool) []string {
	t.Helper()
	var buf bytes.Buffer
	opts := MpileupOptions{
		Inputs:   []string{bam},
		FastaRef: ref,
		Regions:  []string{region},
		NoBAQ:    noBAQ,
		// upstream golden was generated with `-a -AD`; mirror it so the
		// FORMAT column is just PL (no AD subfield) for the comparison.
		Annotate: "-AD",
	}
	if err := MpileupFile(opts, &buf); err != nil {
		t.Fatalf("MpileupFile %s: %v", region, err)
	}
	var out []string
	for _, ln := range strings.Split(buf.String(), "\n") {
		if ln == "" || strings.HasPrefix(ln, "#") {
			continue
		}
		f := strings.Split(ln, "\t")
		if len(f) >= 5 && f[4] == "<*>" {
			out = append(out, ln)
		}
	}
	return out
}

// upstreamStarRecords runs the live `bcftools mpileup -a -AD -r <region>`
// binary over the same BAM/FASTA and keeps records on chrom in the
// inclusive [beg,end] window whose ALT is `<*>`. With `-a -AD` the FORMAT
// column is PL-only, matching mpileupStarRecords, so the lines can be
// compared byte-for-byte. The upstream binary needs a BAM index for `-r`;
// a temp indexed copy is built if the fixture lacks a sidecar.
func upstreamStarRecords(t *testing.T, bin, bam, ref, region, chrom string, beg, end int) []string {
	t.Helper()
	inBAM := bam
	if !fileExists(bam+".bai") && !fileExists(bam+".csi") {
		inBAM = indexedBAMCopy(t, bin, bam)
	}
	cmd := exec.Command(bin, "mpileup", "--no-version", "-a", "-AD",
		"-f", ref, "-r", region, "-Ov", inBAM)
	var out, errBuf bytes.Buffer
	cmd.Stdout = &out
	cmd.Stderr = &errBuf
	if err := cmd.Run(); err != nil {
		t.Fatalf("upstream bcftools mpileup -r %s: %v\n%s", region, err, errBuf.String())
	}
	var recs []string
	for _, ln := range strings.Split(out.String(), "\n") {
		if ln == "" || strings.HasPrefix(ln, "#") {
			continue
		}
		f := strings.Split(ln, "\t")
		if len(f) < 10 || f[0] != chrom || f[4] != "<*>" {
			continue
		}
		pos, err := strconv.Atoi(f[1])
		if err != nil || pos < beg || pos > end {
			continue
		}
		recs = append(recs, ln)
	}
	return recs
}

// TestMpileupBAQGoldens is the slice-3 byte-for-byte parity check. With
// BAQ wired, mpileup over the upstream mpileup.3.bam fixture must match the
// live upstream `bcftools mpileup` output on every `<*>`-only record (those
// carry no slice-4 bias annotations).
//
// Region 17:1-1116 is used: it is free of overlapping read pairs, so the
// not-yet-ported MPLP_SMART_OVERLAPS quality merging cannot perturb the
// comparison. The default (BAQ-on) output is compared to the upstream
// binary; the test also runs -B (BAQ-off) and asserts it diverges on at
// least one column, proving the MPLP_REALN_PARTIAL heuristic does trigger
// BAQ on the indel-bearing columns in this region (and that BAQ then
// matches upstream byte-for-byte, since the compare above is the BAQ-on
// path).
func TestMpileupBAQGoldens(t *testing.T) {
	bam := referenceFixture(t, "mpileup/mpileup.3.bam")
	ref := referenceFixture(t, "mpileup/mpileup.ref.fa")
	bin := upstreamBcftools(t)

	got := mpileupStarRecords(t, bam, ref, "17:1-1116", false)
	want := upstreamStarRecords(t, bin, bam, ref, "17:1-1116", "17", 1, 1116)
	if len(got) != len(want) {
		t.Fatalf("record count: got %d, want %d", len(got), len(want))
	}
	for i := range got {
		if got[i] != want[i] {
			t.Errorf("record %d mismatch (BAQ-on path):\n got:  %s\n want: %s", i, got[i], want[i])
		}
	}

	// -B disables BAQ; it must differ from the BAQ-on path on at least
	// one column, otherwise BAQ is not actually being exercised here.
	noBAQ := mpileupStarRecords(t, bam, ref, "17:1-1116", true)
	if len(noBAQ) != len(got) {
		t.Fatalf("-B record count: got %d, want %d", len(noBAQ), len(got))
	}
	diffs := 0
	for i := range got {
		if noBAQ[i] != got[i] {
			diffs++
		}
	}
	if diffs == 0 {
		t.Error("-B output identical to default: BAQ was never applied, so this test would not actually exercise BAQ")
	}
	t.Logf("BAQ active: %d/%d `<*>`-only columns differ between -B and the default", diffs, len(got))
}

// TestMpileupBCFRoundTrip verifies that -O b output is well-formed BCF
// that round-trips through the project's BCF reader.
func TestMpileupBCFRoundTrip(t *testing.T) {
	dir := t.TempDir()
	famPath := filepath.Join(dir, "ref.fa")
	if err := os.WriteFile(famPath, []byte(">chr1\n"+strings.Repeat("A", 50)+"\n"), 0o644); err != nil {
		t.Fatalf("write fasta: %v", err)
	}
	if err := os.WriteFile(famPath+".fai", []byte("chr1\t50\t6\t50\t51\n"), 0o644); err != nil {
		t.Fatalf("write fai: %v", err)
	}
	samPath := filepath.Join(dir, "in.sam")
	samText := strings.Join([]string{
		"@HD\tVN:1.6\tSO:coordinate",
		"@SQ\tSN:chr1\tLN:50",
		"@RG\tID:rg1\tSM:s1",
		"r1\t0\tchr1\t1\t60\t8M\t*\t0\t0\tAAAACAAA\t????????\tRG:Z:rg1",
		"r2\t0\tchr1\t1\t60\t8M\t*\t0\t0\tAAAACAAA\t????????\tRG:Z:rg1",
		"",
	}, "\n")
	if err := os.WriteFile(samPath, []byte(samText), 0o644); err != nil {
		t.Fatalf("write sam: %v", err)
	}
	var buf bytes.Buffer
	opts := MpileupOptions{
		Inputs:       []string{samPath},
		FastaRef:     famPath,
		OutputFormat: OutputBCF,
	}
	if err := MpileupFile(opts, &buf); err != nil {
		t.Fatalf("MpileupFile -O b: %v", err)
	}
	if buf.Len() == 0 {
		t.Fatal("BCF output is empty")
	}
	// The BCF stream is BGZF-wrapped; a well-formed file starts with the
	// gzip magic 0x1f 0x8b.
	out := buf.Bytes()
	if len(out) < 2 || out[0] != 0x1f || out[1] != 0x8b {
		t.Errorf("BCF output missing BGZF/gzip magic, got % x", out[:min2(4, len(out))])
	}
}

func min2(a, b int) int {
	if a < b {
		return a
	}
	return b
}

// TestParseMpileupRegionSpec covers the chr / chr:beg / chr:beg-end
// shapes plus open-ended chr:beg-.
func TestParseMpileupRegionSpec(t *testing.T) {
	chromLen := map[string]int{"chr1": 1000, "chr2": 500}
	cases := []struct {
		in              string
		chrom           string
		beg, end        int
		wantErr         bool
		acceptOpenEnded bool
	}{
		{in: "chr1", chrom: "chr1", beg: 1, end: 1000},
		{in: "chr1:100", chrom: "chr1", beg: 100, end: 100},
		{in: "chr1:100-200", chrom: "chr1", beg: 100, end: 200},
		{in: "chr1:50-", chrom: "chr1", beg: 50, end: 1000, acceptOpenEnded: true},
		{in: "chrX", chrom: "chrX", beg: 1, end: 1 << 30, acceptOpenEnded: true},
		{in: "chr1:bad", wantErr: true},
		{in: "chr1:1-bad", wantErr: true},
	}
	for _, tc := range cases {
		chrom, beg, end, err := parseMpileupRegionSpec(tc.in, chromLen)
		if tc.wantErr {
			if err == nil {
				t.Errorf("%q: expected error", tc.in)
			}
			continue
		}
		if err != nil {
			t.Errorf("%q: unexpected error %v", tc.in, err)
			continue
		}
		if chrom != tc.chrom || beg != tc.beg {
			t.Errorf("%q: got (%s, %d, %d) want (%s, %d, %d)", tc.in, chrom, beg, end, tc.chrom, tc.beg, tc.end)
		}
		if !tc.acceptOpenEnded && end != tc.end {
			t.Errorf("%q: end got %d want %d", tc.in, end, tc.end)
		}
	}
}

// TestRegionContains is a tiny smoke test on regionContains.
func TestRegionContains(t *testing.T) {
	if !regionContains(nil, "chr1", 100) {
		t.Error("nil windows should always contain")
	}
	w := map[string][][2]int{
		"chr1": {{100, 200}, {300, 400}},
	}
	if !regionContains(w, "chr1", 150) {
		t.Error("chr1:150 should be in [100,200]")
	}
	if regionContains(w, "chr1", 250) {
		t.Error("chr1:250 should NOT be in either window")
	}
	if regionContains(w, "chr2", 150) {
		t.Error("chr2 not in window map")
	}
}

// TestValidateMpileupOptions covers the default substitutions (-d, -Q,
// --max-bq) and that previously deferred flags are now accepted.
func TestValidateMpileupOptions(t *testing.T) {
	// RedoBAQ is wired in slice 3: no longer rejected.
	opts := MpileupOptions{RedoBAQ: true}
	if err := validateMpileupOptions(&opts); err != nil {
		t.Errorf("RedoBAQ should be accepted (slice 3), got %v", err)
	}
	// BCF output is now supported (slice 2): no rejection.
	opts = MpileupOptions{OutputFormat: OutputBCF}
	if err := validateMpileupOptions(&opts); err != nil {
		t.Errorf("OutputBCF should be accepted, got %v", err)
	}
	opts = MpileupOptions{}
	if err := validateMpileupOptions(&opts); err != nil {
		t.Errorf("default opts should be OK, got %v", err)
	}
	if opts.MaxDepth != DefaultMpileupMaxDepth {
		t.Errorf("default MaxDepth = %d, want %d", opts.MaxDepth, DefaultMpileupMaxDepth)
	}
	if opts.MinBQ != DefaultMpileupMinBQ {
		t.Errorf("default MinBQ = %d, want %d", opts.MinBQ, DefaultMpileupMinBQ)
	}
	if opts.MaxBQ != 60 {
		t.Errorf("default MaxBQ = %d, want 60", opts.MaxBQ)
	}
}

// TestMpileupEndToEndSAM exercises the full pipeline: a hand-built SAM
// + a tiny FASTA, run through MpileupFile, and assert the VCF output.
//
// Fixture:
//
//	REF: 100bp of A on chr1 (positions 1..100).
//	Reads: 4x 10bp, all aligned to chr1:1..10. The 5th base is mutated
//	in 2 of the 4 reads (REF=A, ALT=C at position 5) → expect a
//	heterozygous-looking site at chr1:5.
//
// Asserts: a single VCF record at chr1:5 with REF=A, ALT=C, INFO/DP=4
// (one base per read at that column), and a PL triple that picks 0/1
// or a low PL[1].
//
// The fixture runs with NoBAQ (the variant-detection intent): a 10 bp
// read is far too short for BAQ to align confidently, so BAQ would
// collapse the lone mismatch's quality to 0 and the ALT would vanish.
// BAQ behaviour itself is exercised against upstream goldens in
// TestMpileupBAQGoldens.
func TestMpileupEndToEndSAM(t *testing.T) {
	dir := t.TempDir()
	famPath := filepath.Join(dir, "ref.fa")
	famSeq := strings.Repeat("A", 100)
	famContent := ">chr1\n" + famSeq + "\n"
	if err := os.WriteFile(famPath, []byte(famContent), 0o644); err != nil {
		t.Fatalf("write fasta: %v", err)
	}
	// .fai sidecar so OpenRandomAccess works (FAI: name, length,
	// offset, linebases, linewidth).
	faiPath := famPath + ".fai"
	if err := os.WriteFile(faiPath, []byte("chr1\t100\t6\t100\t101\n"), 0o644); err != nil {
		t.Fatalf("write fai: %v", err)
	}

	// 4 reads in SAM text format. Read 3 and 4 have a C at column 5
	// (0-based offset 4); reads 1 and 2 have an A.
	samPath := filepath.Join(dir, "in.sam")
	// SAM column 11 (QUAL) is ASCII Phred+33: '?' = 30.
	sam := strings.Join([]string{
		"@HD\tVN:1.6\tSO:coordinate",
		"@SQ\tSN:chr1\tLN:100",
		"@RG\tID:rg1\tSM:sample1",
		"r1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tAAAAAAAAAA\t??????????\tRG:Z:rg1",
		"r2\t0\tchr1\t1\t60\t10M\t*\t0\t0\tAAAAAAAAAA\t??????????\tRG:Z:rg1",
		"r3\t0\tchr1\t1\t60\t10M\t*\t0\t0\tAAAACAAAAA\t??????????\tRG:Z:rg1",
		"r4\t0\tchr1\t1\t60\t10M\t*\t0\t0\tAAAACAAAAA\t??????????\tRG:Z:rg1",
		"",
	}, "\n")
	if err := os.WriteFile(samPath, []byte(sam), 0o644); err != nil {
		t.Fatalf("write sam: %v", err)
	}

	var buf bytes.Buffer
	opts := MpileupOptions{
		Inputs:   []string{samPath},
		FastaRef: famPath,
		NoBAQ:    true,
	}
	if err := MpileupFile(opts, &buf); err != nil {
		t.Fatalf("MpileupFile: %v", err)
	}

	out := buf.String()
	if !strings.Contains(out, "##fileformat=VCFv4.2") {
		t.Errorf("missing VCF header in output:\n%s", out)
	}
	if !strings.Contains(out, "##contig=<ID=chr1,length=100>") {
		t.Errorf("missing contig line in output:\n%s", out)
	}
	if !strings.Contains(out, "##INFO=<ID=DP") || !strings.Contains(out, "##INFO=<ID=I16") {
		t.Errorf("missing DP / I16 INFO lines:\n%s", out)
	}
	if !strings.Contains(out, "##FORMAT=<ID=PL") {
		t.Errorf("missing PL FORMAT line:\n%s", out)
	}
	if !strings.Contains(out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1") {
		t.Errorf("missing #CHROM header line with sample name:\n%s", out)
	}

	// Find the data lines (no leading '#').
	var data []string
	for _, ln := range strings.Split(out, "\n") {
		if ln == "" || strings.HasPrefix(ln, "#") {
			continue
		}
		data = append(data, ln)
	}
	if len(data) == 0 {
		t.Fatalf("expected ≥1 VCF record, got 0:\n%s", out)
	}
	// Look for the chr1:5 line.
	var found string
	for _, ln := range data {
		fields := strings.Split(ln, "\t")
		if len(fields) >= 2 && fields[0] == "chr1" && fields[1] == "5" {
			found = ln
			break
		}
	}
	if found == "" {
		t.Fatalf("expected a record at chr1:5, got:\n%s", strings.Join(data, "\n"))
	}
	f := strings.Split(found, "\t")
	if f[3] != "A" {
		t.Errorf("REF at chr1:5 = %q, want A", f[3])
	}
	if !strings.Contains(f[4], "C") {
		t.Errorf("ALT at chr1:5 = %q, want to include C", f[4])
	}
	if !strings.Contains(f[7], "DP=4") {
		t.Errorf("INFO at chr1:5 = %q, want DP=4", f[7])
	}
	if !strings.Contains(f[7], "I16=") {
		t.Errorf("INFO at chr1:5 = %q, want I16= tag", f[7])
	}
	// FORMAT is PL:AD by default (B2B_FMT_AD on in DefaultMpileupFmtFlag).
	if f[8] != "PL:AD" {
		t.Errorf("FORMAT at chr1:5 = %q, want PL:AD", f[8])
	}
	// chr1:5 is a heterozygous-looking site: REF=A, ALT=C, plus the
	// appended <*> unseen allele => 3 alleles => 6 PL values.
	sample := strings.Split(f[9], ":")
	if len(sample) != 2 {
		t.Fatalf("sample at chr1:5 = %q, want PL:AD subfields", f[9])
	}
	plParts := strings.Split(sample[0], ",")
	if len(plParts) != 6 {
		t.Errorf("PL at chr1:5 = %q, want 6 comma-separated values", sample[0])
	}
	adParts := strings.Split(sample[1], ",")
	if len(adParts) != 3 {
		t.Errorf("AD at chr1:5 = %q, want 3 comma-separated values", sample[1])
	}
	// Every covered position (chr1:1..10) emits a record now, not just
	// the variant site.
	if len(data) != 10 {
		t.Errorf("got %d records, want 10 (one per covered position)", len(data))
	}
	// Non-variant positions still carry the <*> allele.
	for _, ln := range data {
		ff := strings.Split(ln, "\t")
		if len(ff) >= 5 && !strings.Contains(ff[4], "<*>") {
			t.Errorf("record %s:%s ALT=%q missing <*>", ff[0], ff[1], ff[4])
		}
	}
}

// TestMpileupFileMissingFasta exercises the "missing required flag"
// branch.
func TestMpileupFileMissingFasta(t *testing.T) {
	var buf bytes.Buffer
	err := MpileupFile(MpileupOptions{Inputs: []string{"x.bam"}}, &buf)
	if err == nil {
		t.Error("expected error when -f is empty")
	}
}

// TestMpileupKeepRecordFilters covers the per-read filters: unmapped,
// secondary, supplementary, dup, qcfail, orphan-pair (without -A),
// and MAPQ floor.
func TestMpileupKeepRecordFilters(t *testing.T) {
	// helpers
	mk := func(flag uint16, mq uint8) *samRecordStub {
		return &samRecordStub{Flag: flag, MapQ: mq, RName: "chr1", Pos: 1}
	}
	checkDrops := []struct {
		name string
		rec  *samRecordStub
		want bool
	}{
		{"unmapped", mk(0x4, 60), false},
		{"secondary", mk(0x100, 60), false},
		{"qcfail", mk(0x200, 60), false},
		{"duplicate", mk(0x400, 60), false},
		// Supplementary alignments pass the default mask in upstream
		// mpileup (0x704 = UNMAP|SECONDARY|QCFAIL|DUP, no SUPP). The
		// earlier behaviour was a regression caught in PR #111 review.
		{"supplementary keeps", mk(0x800, 60), true},
		{"mate-unmapped without -A", mk(0x1|0x8, 60), false},
		{"not proper-pair without -A", mk(0x1, 60), false},
		{"low MAPQ", mk(0x2|0x1, 9), false},
		{"OK proper-pair", mk(0x2|0x1, 60), true},
		{"OK single", mk(0, 60), true},
	}
	opts := MpileupOptions{MinMQ: 10}
	// Apply the same default RflagSkipAnySet mask that
	// validateMpileupOptions does (mirroring mpileup.c:1392): UNMAP |
	// SECONDARY | QCFAIL | DUP. The salvaged mpileupKeepRecord no longer
	// hardcodes this default (so --ff can REPLACE it); callers that build
	// MpileupOptions by hand must go through validateMpileupOptions or set
	// the mask themselves.
	opts.RflagSkipAnySet = uint16(sam.FlagUnmapped) | uint16(sam.FlagSecondary) |
		uint16(sam.FlagQCFail) | uint16(sam.FlagDuplicate)
	for _, tc := range checkDrops {
		got := mpileupKeepRecord(tc.rec.toRecord(), opts)
		if got != tc.want {
			t.Errorf("%s: keep=%v want %v", tc.name, got, tc.want)
		}
	}
	// -A flips orphan and not-proper-pair to keep.
	optsA := MpileupOptions{CountOrphans: true, RflagSkipAnySet: opts.RflagSkipAnySet}
	if !mpileupKeepRecord(mk(0x1|0x8, 60).toRecord(), optsA) {
		t.Error("-A should keep mate-unmapped paired reads")
	}
	if !mpileupKeepRecord(mk(0x1, 60).toRecord(), optsA) {
		t.Error("-A should keep not-proper-pair reads")
	}
}

// TestParseFormatFlag covers the parser for the `-a/--annotate` token
// list (mpileup.c:1141 parse_format_flag). It checks token recognition
// (bare names, FORMAT/* and INFO/* prefixes), exclusion with "-", and
// the upstream default's compatibility with `-AD`.
func TestParseFormatFlag(t *testing.T) {
	cases := []struct {
		name  string
		seed  uint32
		input string
		want  uint32
	}{
		{"empty leaves seed", DefaultMpileupFmtFlag, "", DefaultMpileupFmtFlag},
		{"bare AD on zero", 0, "AD", B2BFmtAD},
		{"format-prefix AD", 0, "FORMAT/AD", B2BFmtAD},
		{"info-prefix AD", 0, "INFO/AD", B2BInfoAD},
		{"add DP,AD,SP", 0, "DP,AD,SP", B2BFmtDP | B2BFmtAD | B2BFmtSP},
		{"strip AD from default", DefaultMpileupFmtFlag, "-AD", DefaultMpileupFmtFlag &^ B2BFmtAD},
		{"toggle ADF,ADR", 0, "ADF,ADR", B2BFmtADF | B2BFmtADR},
		{"info NM/NMBZ", 0, "INFO/NM,INFO/NMBZ", B2BInfoNM | B2BInfoNMBZ},
		{"case insensitive", 0, "ad", B2BFmtAD},
	}
	for _, tc := range cases {
		t.Run(tc.name, func(t *testing.T) {
			flag := tc.seed
			if err := parseFormatFlag(&flag, tc.input); err != nil {
				t.Fatalf("parseFormatFlag(%q): %v", tc.input, err)
			}
			if flag != tc.want {
				t.Errorf("parseFormatFlag(%q) = %#x, want %#x", tc.input, flag, tc.want)
			}
		})
	}

	// Unknown tag must error out.
	var f uint32
	if err := parseFormatFlag(&f, "NOT_A_TAG"); err == nil {
		t.Error("parseFormatFlag(\"NOT_A_TAG\") returned nil, want error")
	}
}
