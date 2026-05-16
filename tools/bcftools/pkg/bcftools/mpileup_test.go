package bcftools

import (
	"bytes"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/sam"
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

// TestMpileupDiploidPL_HomRef checks that an all-REF pile yields PL[0]=0
// and PL[1], PL[2] >> 0. A 10-base pile of REF=A with quality 30 gives
// strong support for hom-ref.
func TestMpileupDiploidPL_HomRef(t *testing.T) {
	bases := make([]mpileupBase, 0, 10)
	for i := 0; i < 10; i++ {
		bases = append(bases, mpileupBase{base: 'A', qual: 30})
	}
	pl := mpileupDiploidPL(bases, 'A', 'C')
	if pl[0] != 0 {
		t.Errorf("hom-ref PL[0/0] = %d, want 0", pl[0])
	}
	if pl[1] <= 20 {
		t.Errorf("hom-ref PL[0/1] = %d, want > 20", pl[1])
	}
	if pl[2] <= pl[1] {
		t.Errorf("hom-ref PL[1/1] = %d, want > PL[0/1]=%d", pl[2], pl[1])
	}
}

// TestMpileupDiploidPL_HomAlt is the mirror: a 10-base pile of all-C
// bases when REF=A, ALT=C should yield PL[2]=0 and PL[0] >> 0.
func TestMpileupDiploidPL_HomAlt(t *testing.T) {
	bases := make([]mpileupBase, 0, 10)
	for i := 0; i < 10; i++ {
		bases = append(bases, mpileupBase{base: 'C', qual: 30})
	}
	pl := mpileupDiploidPL(bases, 'A', 'C')
	if pl[2] != 0 {
		t.Errorf("hom-alt PL[1/1] = %d, want 0", pl[2])
	}
	if pl[0] <= 20 {
		t.Errorf("hom-alt PL[0/0] = %d, want > 20", pl[0])
	}
}

// TestMpileupDiploidPL_Het verifies that a balanced REF/ALT pile picks
// the heterozygous genotype 0/1 (PL[1]=0) and that both homozygous
// hypotheses are penalised symmetrically.
func TestMpileupDiploidPL_Het(t *testing.T) {
	bases := make([]mpileupBase, 0, 10)
	for i := 0; i < 5; i++ {
		bases = append(bases, mpileupBase{base: 'A', qual: 30})
	}
	for i := 0; i < 5; i++ {
		bases = append(bases, mpileupBase{base: 'C', qual: 30})
	}
	pl := mpileupDiploidPL(bases, 'A', 'C')
	if pl[1] != 0 {
		t.Errorf("het PL[0/1] = %d, want 0", pl[1])
	}
	if pl[0] <= 0 || pl[2] <= 0 {
		t.Errorf("het PL ends should be > 0, got (%d,%d)", pl[0], pl[2])
	}
	// Symmetric: hom-ref vs hom-alt penalty should be equal-ish (the
	// rounding can shift them by 1 phred unit).
	delta := pl[0] - pl[2]
	if delta < -1 || delta > 1 {
		t.Errorf("het PL should be symmetric, got [%d,%d,%d]", pl[0], pl[1], pl[2])
	}
}

// TestMpileupDiploidPL_NoCoverage exercises the zero-base edge case
// (every PL value should be 0 — uninformative).
func TestMpileupDiploidPL_NoCoverage(t *testing.T) {
	pl := mpileupDiploidPL(nil, 'A', 'C')
	if pl != [3]int{0, 0, 0} {
		t.Errorf("zero-cov PL = %v, want [0,0,0]", pl)
	}
}

// TestMpileupDiploidPL_LowQuality: high-coverage but with low qualities
// (Q=2) should give muddled likelihoods (small numerical gap).
func TestMpileupDiploidPL_LowQuality(t *testing.T) {
	bases := make([]mpileupBase, 0, 10)
	for i := 0; i < 5; i++ {
		bases = append(bases, mpileupBase{base: 'A', qual: 2})
	}
	for i := 0; i < 5; i++ {
		bases = append(bases, mpileupBase{base: 'C', qual: 2})
	}
	pl := mpileupDiploidPL(bases, 'A', 'C')
	if pl[1] != 0 {
		t.Errorf("low-Q het PL[0/1] = %d, want 0", pl[1])
	}
	if pl[0] > 5 {
		t.Errorf("low-Q het: PL[0/0] = %d should be small, got > 5", pl[0])
	}
}

// TestMpileupI16 verifies the per-strand / per-quality count layout
// returned by mpileupI16.
func TestMpileupI16(t *testing.T) {
	bases := []mpileupBase{
		{base: 'A', qual: 20, mapq: 30, isReverse: false}, // ref forward
		{base: 'A', qual: 30, mapq: 40, isReverse: true},  // ref reverse
		{base: 'C', qual: 25, mapq: 35, isReverse: false}, // alt forward
		{base: 'C', qual: 35, mapq: 45, isReverse: true},  // alt reverse
	}
	got := mpileupI16(bases, 'A')
	if got[0] != 1 || got[1] != 1 || got[2] != 1 || got[3] != 1 {
		t.Errorf("strand counts wrong: %v", got[:4])
	}
	if got[4] != 50 {
		t.Errorf("sum baseQ ref = %v, want 50 (20+30)", got[4])
	}
	if got[5] != 1300 { // 20*20 + 30*30 = 400+900
		t.Errorf("sum baseQ^2 ref = %v, want 1300", got[5])
	}
	if got[6] != 60 {
		t.Errorf("sum baseQ non-ref = %v, want 60 (25+35)", got[6])
	}
	if got[8] != 70 {
		t.Errorf("sum mapq ref = %v, want 70 (30+40)", got[8])
	}
}

// TestChooseALTs verifies ALT-choice ordering (descending count, then
// lexicographic tie-break).
func TestChooseALTs(t *testing.T) {
	perSample := [][]mpileupBase{
		{
			{base: 'A'}, {base: 'A'}, {base: 'C'},
			{base: 'C'}, {base: 'C'}, {base: 'G'},
		},
	}
	got := chooseALTs(perSample, 'A')
	// v1 caps at 1 ALT to keep biallelic PL formatting valid; "C"
	// is the higher-count non-REF base. Multi-allelic ALT grid lands
	// with the MAQ port (PR #111 reviewer-caught regression).
	if string(got) != "C" {
		t.Errorf("chooseALTs got %q want \"C\"", string(got))
	}
	// Pure-reference site.
	pure := [][]mpileupBase{{{base: 'A'}, {base: 'A'}}}
	if a := chooseALTs(pure, 'A'); len(a) != 0 {
		t.Errorf("pure ref: ALTs %v want []", a)
	}
	// Cap at 3.
	multi := [][]mpileupBase{{
		{base: 'C'}, {base: 'G'}, {base: 'T'},
		{base: 'N'},
	}}
	if a := chooseALTs(multi, 'A'); len(a) > 3 {
		t.Errorf("cap: ALTs %v exceeds 3", a)
	}
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

// TestValidateMpileupOptions covers the v1 hard-rejections and default
// substitutions (-d, -Q, --max-bq).
func TestValidateMpileupOptions(t *testing.T) {
	opts := MpileupOptions{RedoBAQ: true}
	if err := validateMpileupOptions(&opts); err == nil {
		t.Error("RedoBAQ should be rejected")
	}
	opts = MpileupOptions{OutputFormat: OutputBCF}
	if err := validateMpileupOptions(&opts); err == nil {
		t.Error("OutputBCF should be rejected")
	}
	opts = MpileupOptions{OutputFormat: OutputBCFUncompressed}
	if err := validateMpileupOptions(&opts); err == nil {
		t.Error("OutputBCFUncompressed should be rejected")
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
	if f[8] != "PL" {
		t.Errorf("FORMAT at chr1:5 = %q, want PL", f[8])
	}
	plParts := strings.Split(f[9], ",")
	if len(plParts) != 3 {
		t.Errorf("PL at chr1:5 = %q, want 3 comma-separated values", f[9])
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
	for _, tc := range checkDrops {
		got := mpileupKeepRecord(tc.rec.toRecord(), opts)
		if got != tc.want {
			t.Errorf("%s: keep=%v want %v", tc.name, got, tc.want)
		}
	}
	// -A flips orphan and not-proper-pair to keep.
	optsA := MpileupOptions{CountOrphans: true}
	if !mpileupKeepRecord(mk(0x1|0x8, 60).toRecord(), optsA) {
		t.Error("-A should keep mate-unmapped paired reads")
	}
	if !mpileupKeepRecord(mk(0x1, 60).toRecord(), optsA) {
		t.Error("-A should keep not-proper-pair reads")
	}
}
