package bcftools

import (
	"reflect"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/sam"
)

// samRecordForTest builds a minimal mapped sam.Record from a CIGAR, 1-based
// position and sequence, for the CIGAR-walk unit tests.
func samRecordForTest(t *testing.T, cigar string, pos1 int, seq string) *sam.Record {
	t.Helper()
	c, err := sam.ParseCigar(cigar)
	if err != nil {
		t.Fatalf("ParseCigar(%q): %v", cigar, err)
	}
	return &sam.Record{
		QName: "r",
		RName: "chr1",
		Pos:   int64(pos1),
		MapQ:  60,
		Cigar: c,
		Seq:   seq,
	}
}

// Binary-free unit tests for the pure vrfs helpers. These pass with the
// upstream submodules unpopulated (no bcftools/samtools binary needed): they
// exercise the sites/aln-list parsers, the VAF binning, the per-site ref/alt
// counting from a synthetic pileup column, the profile mean/var aggregation,
// the variance-profile rescaling, the score function and the C-style float
// formatting — all the logic the oracle test then validates end-to-end.

func TestUnitVrfsNN2Bin(t *testing.T) {
	tests := []struct {
		nbins, nref, nalt, want int
	}{
		{20, 0, 0, -1}, // no reads
		{20, 5, 0, 0},  // no alt -> bin 0
		{20, 4, 1, 3},  // floor(19*1/5)=3
		{20, 7, 3, 5},  // floor(19*3/10)=5
		{20, 0, 1, 19}, // all alt -> top bin
		{10, 4, 6, 5},  // floor(9*6/10)=5
		{25, 1, 1, 12}, // floor(24*1/2)=12
		{20, 8, 0, 0},  // ref-only
		{20, 4, 6, 11}, // floor(19*6/10)=11
	}
	for _, tt := range tests {
		if got := nn2bin(tt.nbins, tt.nref, tt.nalt); got != tt.want {
			t.Errorf("nn2bin(%d,%d,%d)=%d, want %d", tt.nbins, tt.nref, tt.nalt, got, tt.want)
		}
	}
}

func TestUnitVrfsParseSites(t *testing.T) {
	in := `# comment
chr1	11	G	A
chr1 11 G T
chr2	5	AT	A
chr2	7	A	ACG

chr3	9	ACGT	TGCA
`
	sites, err := parseVrfsSites(strings.NewReader(in), 20)
	if err != nil {
		t.Fatalf("parseVrfsSites: %v", err)
	}
	want := []struct {
		chrom    string
		pos0     int
		ref, alt string
		isIndel  bool
		altClass int
	}{
		{"chr1", 10, "G", "A", false, 0},
		{"chr1", 10, "G", "T", false, 3},
		{"chr2", 4, "AT", "A", true, 4},  // deletion
		{"chr2", 6, "A", "ACG", true, 4}, // insertion
		{"chr3", 8, "A", "T", false, 3},  // equal-length MNV -> SNV (first base)
	}
	if len(sites) != len(want) {
		t.Fatalf("got %d sites, want %d", len(sites), len(want))
	}
	for i, w := range want {
		s := sites[i]
		if s.chrom != w.chrom || s.pos0 != w.pos0 || s.ref != w.ref || s.alt != w.alt {
			t.Errorf("site %d = {%s %d %s %s}, want {%s %d %s %s}",
				i, s.chrom, s.pos0, s.ref, s.alt, w.chrom, w.pos0, w.ref, w.alt)
		}
		if s.isIndel() != w.isIndel {
			t.Errorf("site %d isIndel=%v, want %v", i, s.isIndel(), w.isIndel)
		}
		if s.altClass() != w.altClass {
			t.Errorf("site %d altClass=%d, want %d", i, s.altClass(), w.altClass)
		}
		if len(s.dist) != 20 {
			t.Errorf("site %d dist len=%d, want 20", i, len(s.dist))
		}
	}
}

func TestUnitVrfsParseAlnList(t *testing.T) {
	in := "  /path/a.bam \n\n/path/b.cram\n   \n/path/c.bam\n"
	got, err := parseAlnList(strings.NewReader(in))
	if err != nil {
		t.Fatalf("parseAlnList: %v", err)
	}
	want := []string{"/path/a.bam", "/path/b.cram", "/path/c.bam"}
	if !reflect.DeepEqual(got, want) {
		t.Errorf("parseAlnList = %v, want %v", got, want)
	}
}

func TestUnitVrfsBaseClass(t *testing.T) {
	cases := map[byte]int{'A': 0, 'a': 0, 'C': 1, 'c': 1, 'G': 2, 'g': 2, 'T': 3, 't': 3, 'N': -1, 'x': -1}
	for b, want := range cases {
		if got := baseClass(b); got != want {
			t.Errorf("baseClass(%q)=%d, want %d", b, got, want)
		}
	}
}

// TestUnitVrfsAccumulatorFlush builds a synthetic per-(sample,site) pileup
// directly via the accumulator and checks the histogram increments match the
// nn2bin + min-depth logic of vrfs.c's count loop.
func TestUnitVrfsAccumulatorFlush(t *testing.T) {
	// One SNV site chr1:11 G->A with nbins 20.
	siteA := &vrfsSite{chrom: "chr1", pos0: 10, ref: "G", alt: "A", dist: make([]uint32, 20)}
	acc := newVrfsAccumulator()

	// Sample S1: 4 ref, 1 alt-A -> ntot=5, nalt[A]=1, bin nn2bin(20,4,1)=3.
	c1 := &vrfsPosCounts{}
	c1.site[0] = siteA // alt class A = 0
	c1.ntot = 5
	c1.nalt[0] = 1
	acc.cells["S1\x00chr1\x0010"] = c1

	// Sample S2: 8 ref, 0 alt-A -> bin 0.
	c2 := &vrfsPosCounts{}
	c2.site[0] = siteA
	c2.ntot = 8
	acc.cells["S2\x00chr1\x0010"] = c2

	// Sample S3: depth below min -> skipped.
	c3 := &vrfsPosCounts{}
	c3.site[0] = siteA
	c3.ntot = 2
	c3.nalt[0] = 1
	acc.cells["S3\x00chr1\x0010"] = c3

	acc.flush(20, 5)

	if siteA.nval != 2 {
		t.Fatalf("nval=%d, want 2 (S3 below min-depth)", siteA.nval)
	}
	want := make([]uint32, 20)
	want[0] = 1 // S2
	want[3] = 1 // S1
	if !reflect.DeepEqual(siteA.dist, want) {
		t.Errorf("dist=%v, want %v", siteA.dist, want)
	}
}

// TestUnitVrfsComputeProfileHC checks the mean/var2 aggregation against a
// hand-computed example in the default "hc" recalc mode.
func TestUnitVrfsComputeProfileHC(t *testing.T) {
	// Three sites, each with a single value in one bin, mirroring the live
	// fixture (sites in bins 3, 4, 2 respectively, after the bin-0 ref vote).
	mk := func(b0, bk, k int) *vrfsSite {
		s := &vrfsSite{dist: make([]uint32, 20), nval: b0 + bk}
		s.dist[0] = uint32(b0)
		s.dist[k] = uint32(bk)
		return s
	}
	sites := []*vrfsSite{mk(1, 1, 3), mk(1, 1, 4), mk(1, 1, 2)}
	prof := computeVrfsProfile(sites, vrfsConfig{nbins: 20, recalc: "hc"})

	if prof.nval != 3 {
		t.Fatalf("nval=%d, want 3", prof.nval)
	}
	// Each site normalises to max=1, so bin0 contributes 1 at every site -> mean[0]=1.
	if prof.mean[0] != 1.0 {
		t.Errorf("mean[0]=%v, want 1", prof.mean[0])
	}
	// bins 2,3,4 each get one site contributing 1/3.
	for _, k := range []int{2, 3, 4} {
		if got := prof.mean[k]; got < 0.3333332 || got > 0.3333334 {
			t.Errorf("mean[%d]=%v, want ~0.333333", k, got)
		}
	}
	// hc mode uses the hard-coded var2 table verbatim at nbins=20.
	if prof.var2[1] != vrfsHardcodedVar2[1] {
		t.Errorf("var2[1]=%v, want %v (hc table)", prof.var2[1], vrfsHardcodedVar2[1])
	}
}

// TestUnitVrfsComputeProfileEmpty checks that an all-empty profile yields the
// "nan" MEAN line: the mean is a positive 0.0/0.0 NaN (division by zero nval),
// which glibc printf renders as "nan" (verified against the live upstream
// +vrfs oracle). This previously asserted "-nan", encoding an earlier
// hard-coded NaN sign rather than upstream's actual output.
func TestUnitVrfsComputeProfileEmpty(t *testing.T) {
	sites := []*vrfsSite{{dist: make([]uint32, 20)}} // nval==0, skipped
	prof := computeVrfsProfile(sites, vrfsConfig{nbins: 20, recalc: "hc"})
	if prof.nval != 0 {
		t.Fatalf("nval=%d, want 0", prof.nval)
	}
	if got := vrfsFormatE(prof.mean[0]); got != "nan" {
		t.Errorf("empty mean format = %q, want %q", got, "nan")
	}
}

func TestUnitVrfsRescaleVar2(t *testing.T) {
	// Same length -> identical copy.
	got := vrfsRescaleVar2([]float64{1, 2, 3}, 3)
	if !reflect.DeepEqual(got, []float64{1, 2, 3}) {
		t.Errorf("rescale same-len = %v, want [1 2 3]", got)
	}
	// Endpoints are preserved on interpolation.
	out := vrfsRescaleVar2([]float64{0, 10}, 5)
	if out[0] != 0 || out[len(out)-1] != 10 {
		t.Errorf("rescale endpoints = %v, want [0 ... 10]", out)
	}
	// Linear interpolation: 5 bins over [0,10] -> 0,2.5,5,7.5,10.
	wantMid := []float64{0, 2.5, 5, 7.5, 10}
	for i := range wantMid {
		if d := out[i] - wantMid[i]; d > 1e-9 || d < -1e-9 {
			t.Errorf("rescale[%d]=%v, want %v", i, out[i], wantMid[i])
		}
	}
}

func TestUnitVrfsFormatE(t *testing.T) {
	cases := map[float64]string{
		1.713053e+02: "1.713053e+02",
		0:            "0.000000e+00",
		3.333333e-01: "3.333333e-01",
		7.443527e-10: "7.443527e-10",
	}
	for v, want := range cases {
		if got := vrfsFormatE(v); got != want {
			t.Errorf("vrfsFormatE(%g)=%q, want %q", v, got, want)
		}
	}
}

// TestUnitVrfsReadColumns exercises the CIGAR walk that maps a read onto site
// positions, including the deletion (is_del reads the post-gap base) and the
// indel-on-prior-base cases that match htslib's bam_pileup1_t semantics.
func TestUnitVrfsReadColumns(t *testing.T) {
	mkSites := func(pos ...int) map[int][]*vrfsSite {
		m := map[int][]*vrfsSite{}
		for _, p := range pos {
			m[p] = []*vrfsSite{{}}
		}
		return m
	}

	// Plain 20M read at pos5 (0-based 4): site at 0-based 10 -> base seq[6].
	rec := samRecordForTest(t, "20M", 5, "ACGTACGTACGTACGTACGT")
	cols := vrfsReadColumns(rec, mkSites(10))
	if len(cols) != 1 || cols[0].pos0 != 10 || cols[0].indel {
		t.Fatalf("plain M cols = %+v", cols)
	}
	// seq index 10-4=6 -> 'G'
	if cols[0].base != 'G' {
		t.Errorf("plain M base = %q, want 'G'", cols[0].base)
	}

	// 6M2D12M read: at the last M base (0-based 9) indel=true; the deleted
	// positions (0-based 10,11) read the post-deletion base (query index 6).
	delRec := samRecordForTest(t, "6M2D12M", 5, "ACGTACACGTACGTACGT")
	dcols := vrfsReadColumns(delRec, mkSites(9, 10, 11))
	byPos := map[int]vrfsReadColumn{}
	for _, c := range dcols {
		byPos[c.pos0] = c
	}
	if !byPos[9].indel {
		t.Errorf("col at 0-based 9 should have indel=true (D follows)")
	}
	if byPos[10].indel || byPos[11].indel {
		t.Errorf("is_del columns should have indel=false")
	}
	if byPos[10].base != 'A' || byPos[11].base != 'A' {
		t.Errorf("is_del base = %q/%q, want post-deletion 'A'", byPos[10].base, byPos[11].base)
	}

	// 6M2I12M read: the last M base before the insertion carries indel=true.
	insRec := samRecordForTest(t, "6M2I12M", 5, "ACGTACTTGTACGTACGTAC")
	icols := vrfsReadColumns(insRec, mkSites(9))
	if len(icols) != 1 || !icols[0].indel {
		t.Fatalf("insertion-prior col = %+v, want indel=true", icols)
	}
}
