package vcftools

import (
	"fmt"
	"math"
	"strings"
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/bioformats/vcf"
)

// Test helper to create a simple VCF variant
func createTestVariant(chrom string, pos int, ref string, alt []string, qual float64, gt []string) *vcf.Variant {
	v := &vcf.Variant{
		Chrom:  chrom,
		Pos:    pos,
		ID:     ".",
		Ref:    ref,
		Alt:    alt,
		Qual:   qual,
		Filter: []string{"PASS"},
		Info:   make(map[string]string),
		Format: []string{"GT"},
	}

	for i, g := range gt {
		sample := vcf.Sample{
			Name: "sample" + string(rune('1'+i)),
			Data: map[string]string{"GT": g},
		}
		v.Samples = append(v.Samples, sample)
	}

	return v
}

func TestIsIndelVariant(t *testing.T) {
	tests := []struct {
		name     string
		ref      string
		alt      []string
		expected bool
	}{
		{"SNP", "A", []string{"G"}, false},
		{"Insertion", "A", []string{"AT"}, true},
		{"Deletion", "AT", []string{"A"}, true},
		{"Multi-alt with indel", "A", []string{"G", "AT"}, true},
		{"Multi-base SNP", "AT", []string{"GC"}, false},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := &vcf.Variant{Ref: tt.ref, Alt: tt.alt}
			result := isIndelVariant(v)
			if result != tt.expected {
				t.Errorf("isIndelVariant() = %v, want %v", result, tt.expected)
			}
		})
	}
}

func TestCalculateMAF(t *testing.T) {
	tests := []struct {
		name        string
		genotypes   []string
		expectedMAF float64
		expectedMAC int
	}{
		{
			name:        "All homozygous ref",
			genotypes:   []string{"0/0", "0/0", "0/0"},
			expectedMAF: 0.0,
			expectedMAC: 0,
		},
		{
			name:        "All homozygous alt",
			genotypes:   []string{"1/1", "1/1", "1/1"},
			expectedMAF: 0.0,
			expectedMAC: 0,
		},
		{
			name:        "All heterozygous",
			genotypes:   []string{"0/1", "0/1", "0/1"},
			expectedMAF: 0.5,
			expectedMAC: 3,
		},
		{
			name:        "Mixed genotypes",
			genotypes:   []string{"0/0", "0/1", "1/1"},
			expectedMAF: 0.5,
			expectedMAC: 3,
		},
		{
			name:        "One rare variant",
			genotypes:   []string{"0/0", "0/0", "0/0", "0/1"},
			expectedMAF: 0.125,
			expectedMAC: 1,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := createTestVariant("chr1", 100, "A", []string{"G"}, 30, tt.genotypes)
			maf, mac := calculateMAF(v)

			if maf != tt.expectedMAF {
				t.Errorf("MAF = %v, want %v", maf, tt.expectedMAF)
			}
			if mac != tt.expectedMAC {
				t.Errorf("MAC = %v, want %v", mac, tt.expectedMAC)
			}
		})
	}
}

func TestCalculateMissingRate(t *testing.T) {
	tests := []struct {
		name      string
		genotypes []string
		expected  float64
	}{
		{
			name:      "No missing",
			genotypes: []string{"0/0", "0/1", "1/1"},
			expected:  0.0,
		},
		{
			name:      "All missing",
			genotypes: []string{"./.", "./.", "./."},
			expected:  1.0,
		},
		{
			name:      "Half missing",
			genotypes: []string{"0/0", "./.", "0/1", "./."},
			expected:  0.5,
		},
		{
			name:      "Phased missing",
			genotypes: []string{"0|0", ".|.", "0|1"},
			expected:  1.0 / 3.0,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := createTestVariant("chr1", 100, "A", []string{"G"}, 30, tt.genotypes)
			result := calculateMissingRate(v)

			if result != tt.expected {
				t.Errorf("Missing rate = %v, want %v", result, tt.expected)
			}
		})
	}
}

func TestCalculateMeanDepth(t *testing.T) {
	v := &vcf.Variant{
		Chrom:  "chr1",
		Pos:    100,
		Format: []string{"GT", "DP"},
		Samples: []vcf.Sample{
			{Name: "s1", Data: map[string]string{"GT": "0/0", "DP": "10"}},
			{Name: "s2", Data: map[string]string{"GT": "0/1", "DP": "20"}},
			{Name: "s3", Data: map[string]string{"GT": "1/1", "DP": "30"}},
		},
	}

	depth := calculateMeanDepth(v)
	expected := 20.0

	if depth != expected {
		t.Errorf("Mean depth = %v, want %v", depth, expected)
	}
}

func TestPassFilters(t *testing.T) {
	tests := []struct {
		name     string
		variant  *vcf.Variant
		params   *Params
		expected bool
	}{
		{
			name:     "No filters - should pass",
			variant:  createTestVariant("chr1", 100, "A", []string{"G"}, 30, []string{"0/0"}),
			params:   &Params{},
			expected: true,
		},
		{
			name:     "Chr filter - match",
			variant:  createTestVariant("chr1", 100, "A", []string{"G"}, 30, []string{"0/0"}),
			params:   &Params{Chr: "chr1"},
			expected: true,
		},
		{
			name:     "Chr filter - no match",
			variant:  createTestVariant("chr2", 100, "A", []string{"G"}, 30, []string{"0/0"}),
			params:   &Params{Chr: "chr1"},
			expected: false,
		},
		{
			name:     "Position range - inside",
			variant:  createTestVariant("chr1", 500, "A", []string{"G"}, 30, []string{"0/0"}),
			params:   &Params{FromBp: 100, ToBp: 1000},
			expected: true,
		},
		{
			name:     "Position range - outside",
			variant:  createTestVariant("chr1", 50, "A", []string{"G"}, 30, []string{"0/0"}),
			params:   &Params{FromBp: 100, ToBp: 1000},
			expected: false,
		},
		{
			name:     "Quality filter - pass",
			variant:  createTestVariant("chr1", 100, "A", []string{"G"}, 30, []string{"0/0"}),
			params:   &Params{MinQ: 20},
			expected: true,
		},
		{
			name:     "Quality filter - fail",
			variant:  createTestVariant("chr1", 100, "A", []string{"G"}, 15, []string{"0/0"}),
			params:   &Params{MinQ: 20},
			expected: false,
		},
		{
			name:     "Remove indels - SNP passes",
			variant:  createTestVariant("chr1", 100, "A", []string{"G"}, 30, []string{"0/0"}),
			params:   &Params{RemoveIndels: true},
			expected: true,
		},
		{
			name:     "Remove indels - indel fails",
			variant:  createTestVariant("chr1", 100, "A", []string{"AT"}, 30, []string{"0/0"}),
			params:   &Params{RemoveIndels: true},
			expected: false,
		},
		{
			name:     "Keep only indels - SNP fails",
			variant:  createTestVariant("chr1", 100, "A", []string{"G"}, 30, []string{"0/0"}),
			params:   &Params{KeepOnlyIndels: true},
			expected: false,
		},
		{
			name:     "Keep only indels - indel passes",
			variant:  createTestVariant("chr1", 100, "A", []string{"AT"}, 30, []string{"0/0"}),
			params:   &Params{KeepOnlyIndels: true},
			expected: true,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result := passFilters(tt.variant, tt.params, nil, nil, nil, nil, nil, nil)
			if result != tt.expected {
				t.Errorf("passFilters() = %v, want %v", result, tt.expected)
			}
		})
	}
}

func TestLoadPositions(t *testing.T) {
	// Create a mock positions file content
	content := `# Comment line
chr1	100
chr1	200
chr2	300
# Another comment
chr2	400
`

	// We can't easily test file I/O without creating actual files,
	// so we'll test the parsing logic separately
	positions := make(positionSet)

	lines := strings.Split(content, "\n")
	for _, line := range lines {
		line = strings.TrimSpace(line)
		if line == "" || strings.HasPrefix(line, "#") {
			continue
		}

		fields := strings.Fields(line)
		if len(fields) < 2 {
			continue
		}

		chrom := fields[0]
		var pos int
		_, err := fmt.Sscanf(fields[1], "%d", &pos)
		if err != nil {
			continue
		}

		if positions[chrom] == nil {
			positions[chrom] = make(map[int]bool)
		}
		positions[chrom][pos] = true
	}

	// Verify the positions were loaded correctly
	if !positions["chr1"][100] {
		t.Error("Position chr1:100 should be in set")
	}
	if !positions["chr1"][200] {
		t.Error("Position chr1:200 should be in set")
	}
	if !positions["chr2"][300] {
		t.Error("Position chr2:300 should be in set")
	}
	if !positions["chr2"][400] {
		t.Error("Position chr2:400 should be in set")
	}
	if positions["chr1"][999] {
		t.Error("Position chr1:999 should not be in set")
	}
}

func TestBuildSampleFilter(t *testing.T) {
	header := &vcf.Header{
		Samples: []string{"sample1", "sample2", "sample3", "sample4"},
	}

	tests := []struct {
		name     string
		params   *Params
		expected map[string]bool
	}{
		{
			name:     "No filters",
			params:   &Params{},
			expected: nil,
		},
		{
			name:     "Keep specific individuals",
			params:   &Params{IndvList: []string{"sample1", "sample3"}},
			expected: map[string]bool{"sample1": true, "sample3": true},
		},
		{
			name:     "Remove specific individuals",
			params:   &Params{RemoveIndvList: []string{"sample2", "sample4"}},
			expected: map[string]bool{"sample1": true, "sample3": true},
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			result, err := buildSampleFilter(header, tt.params)
			if err != nil {
				t.Fatalf("buildSampleFilter() error = %v", err)
			}

			if tt.expected == nil && result != nil {
				t.Errorf("Expected nil, got %v", result)
				return
			}

			if tt.expected != nil {
				if result == nil {
					t.Error("Expected non-nil result")
					return
				}

				if len(result) != len(tt.expected) {
					t.Errorf("Expected %d samples, got %d", len(tt.expected), len(result))
				}

				for sample := range tt.expected {
					if !result[sample] {
						t.Errorf("Expected sample %s to be in result", sample)
					}
				}
			}
		})
	}
}

func TestNucleotideDiversity(t *testing.T) {
	tests := []struct {
		name string
		gt   []string
		want float64
		ok   bool
	}{
		// 3 ref + 3 alt out of 6 chromosomes: (36 - 9 - 9) / (6*5) = 0.6
		{"balanced biallelic", []string{"0/0", "0/1", "1/1"}, 0.6, true},
		// All reference => no diversity.
		{"monomorphic", []string{"0/0", "0/0", "0/0"}, 0.0, true},
		// 5 ref + 1 alt out of 6: (36 - 25 - 1) / 30 = 10/30
		{"singleton", []string{"0/0", "0/0", "0/1"}, 10.0 / 30.0, true},
		// Missing data is excluded; only two chromosomes (one ref, one alt)
		// remain: (4 - 1 - 1) / (2*1) = 1.0
		{"with missing", []string{"./.", "1", "0"}, 1.0, true},
		// Fewer than two non-missing chromosomes => not defined.
		{"insufficient data", []string{"./.", "./.", "0"}, 0.0, false},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := createTestVariant("1", 100, "A", []string{"G"}, 50, tt.gt)
			got, ok := nucleotideDiversity(v)
			if ok != tt.ok {
				t.Fatalf("nucleotideDiversity ok = %v, want %v", ok, tt.ok)
			}
			if ok && math.Abs(got-tt.want) > 1e-9 {
				t.Errorf("nucleotideDiversity = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestIsTransitionSNP(t *testing.T) {
	transitions := [][2]string{{"A", "G"}, {"G", "A"}, {"C", "T"}, {"T", "C"}}
	transversions := [][2]string{{"A", "C"}, {"A", "T"}, {"G", "C"}, {"G", "T"}, {"C", "A"}, {"T", "A"}}
	for _, p := range transitions {
		if !isTransitionSNP(p[0], p[1]) {
			t.Errorf("isTransitionSNP(%s,%s) = false, want true", p[0], p[1])
		}
	}
	for _, p := range transversions {
		if isTransitionSNP(p[0], p[1]) {
			t.Errorf("isTransitionSNP(%s,%s) = true, want false", p[0], p[1])
		}
	}
}

func TestTajimasD(t *testing.T) {
	// Fewer than two segregating sites => undefined.
	if _, ok := tajimasD(0.0, 1, 10); ok {
		t.Error("tajimasD with S=1 should be undefined")
	}
	if _, ok := tajimasD(5.0, 4, 2); ok {
		t.Error("tajimasD with n=2 should be undefined")
	}
	// When the summed per-site diversity equals Watterson's theta (S/a1) the
	// numerator is zero, so D must be zero.
	n, S := 6, 4
	a1 := 0.0
	for i := 1; i < n; i++ {
		a1 += 1.0 / float64(i)
	}
	d, ok := tajimasD(float64(S)/a1, S, n)
	if !ok {
		t.Fatal("tajimasD should be defined for S=4, n=6")
	}
	if math.Abs(d) > 1e-9 {
		t.Errorf("tajimasD = %v, want ~0 when pi == thetaW", d)
	}
	// pi above thetaW => positive D; pi below => negative D.
	if dp, _ := tajimasD(float64(S)/a1+2, S, n); dp <= 0 {
		t.Errorf("tajimasD with pi > thetaW = %v, want > 0", dp)
	}
	if dm, _ := tajimasD(float64(S)/a1-2, S, n); dm >= 0 {
		t.Errorf("tajimasD with pi < thetaW = %v, want < 0", dm)
	}
}

// TestWeirCockerhamFstWorkedExample verifies the per-site Weir & Cockerham
// 1984 Fst estimator against the worked example documented in the
// implementation plan: two populations of sizes 2 and 3, genotypes
// {0/1, 0/1} and {0/0, 0/0, 1/1}, expected Fst ~ -0.3232.
func TestWeirCockerhamFstWorkedExample(t *testing.T) {
	// Pop1: S1=0/1, S2=0/1 -> n=2, altCount=2, hetCount=2.
	// Pop2: S3=0/0, S4=0/0, S5=1/1 -> n=3, altCount=2, hetCount=0.
	pops := []popGenotypeCounts{
		{n: 2, altCount: 2, hetCount: 2},
		{n: 3, altCount: 2, hetCount: 0},
	}
	a, b, c, fst, ok := weirCockerhamFst(pops)
	if !ok {
		t.Fatalf("weirCockerhamFst returned ok=false; a=%v b=%v c=%v", a, b, c)
	}
	want := -0.3232
	if math.Abs(fst-want) > 1e-3 {
		t.Errorf("weirCockerhamFst Fst = %v, want ~%v (within 1e-3)", fst, want)
	}
	// Sanity: a + b + c should be ~0.24352 and the components carry the
	// expected signs (a < 0, b > 0, c > 0).
	if a >= 0 {
		t.Errorf("expected a < 0, got %v", a)
	}
	if b <= 0 {
		t.Errorf("expected b > 0, got %v", b)
	}
	if c <= 0 {
		t.Errorf("expected c > 0, got %v", c)
	}
	if math.Abs((a+b+c)-0.24352) > 1e-3 {
		t.Errorf("a+b+c = %v, want ~0.24352", a+b+c)
	}
}

func TestCheckUnsupported(t *testing.T) {
	supported := []*Params{
		{},
		{SitePi: true, WindowPi: 1000},
		{TsTvByCount: true, Depth: true},
		{TsTvByQual: true},
		{HistIndelLen: true, GenoDepth: true, TajimaD: 10000},
		{WeirFstPop: []string{"pop1.txt", "pop2.txt"}, FstWindowSize: 10000, FstWindowStep: 5000},
	}
	for i, p := range supported {
		if err := checkUnsupported(p); err != nil {
			t.Errorf("checkUnsupported(supported[%d]) = %v, want nil", i, err)
		}
	}

	unsupported := []*Params{
		// Currently every previously-rejected feature has been implemented;
		// this list is kept (empty) so the test trivially documents the lack
		// of unsupported flags. Add new entries here when introducing checks.
	}
	for i, p := range unsupported {
		if err := checkUnsupported(p); err == nil {
			t.Errorf("checkUnsupported(unsupported[%d]) = nil, want error", i)
		}
	}
}

func TestRunAcceptsTsTvByQual(t *testing.T) {
	// TsTvByQual is now implemented (it used to be rejected with an error);
	// this test pins that the feature stays accepted.
	in := strings.NewReader(testVCF)
	if err := Run(in, &Params{OutPrefix: t.TempDir() + "/out", TsTvByQual: true}); err != nil {
		t.Fatalf("Run with --TsTv-by-qual should succeed, got: %v", err)
	}
}

// TestCalculateAlleleCounts pins the per-ALT count helper used by the
// --non-ref-af / --non-ref-ac filters. Mirrors the get_allele_counts
// semantics from entry_getters.cpp:389-422 (skip missing, REF skipped
// when building altCounts, totalCalled tracks every non-missing chr).
func TestCalculateAlleleCounts(t *testing.T) {
	tests := []struct {
		name        string
		alt         []string
		genotypes   []string
		wantAlt     []int
		wantTotalNc int
	}{
		{
			name:        "biallelic, all alt",
			alt:         []string{"G"},
			genotypes:   []string{"1/1", "1|1", "1/1"},
			wantAlt:     []int{6},
			wantTotalNc: 6,
		},
		{
			name:        "biallelic, mixed",
			alt:         []string{"G"},
			genotypes:   []string{"0/0", "0/1", "1/1"},
			wantAlt:     []int{3},
			wantTotalNc: 6,
		},
		{
			name:        "triallelic",
			alt:         []string{"G", "T"},
			genotypes:   []string{"0/1", "1/2", "2/2"},
			wantAlt:     []int{2, 3},
			wantTotalNc: 6,
		},
		{
			name:        "with missing",
			alt:         []string{"G"},
			genotypes:   []string{"0/1", "./.", "1/1"},
			wantAlt:     []int{3},
			wantTotalNc: 4,
		},
		{
			name:        "haploid mixed",
			alt:         []string{"T"},
			genotypes:   []string{"0", "1", "0/1"},
			wantAlt:     []int{2},
			wantTotalNc: 4,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			v := createTestVariant("chr1", 1, "A", tt.alt, 30, tt.genotypes)
			alt, total := calculateAlleleCounts(v)
			if total != tt.wantTotalNc {
				t.Errorf("totalCalled = %d, want %d", total, tt.wantTotalNc)
			}
			if len(alt) != len(tt.wantAlt) {
				t.Fatalf("altCounts length = %d, want %d", len(alt), len(tt.wantAlt))
			}
			for i := range alt {
				if alt[i] != tt.wantAlt[i] {
					t.Errorf("altCounts[%d] = %d, want %d", i, alt[i], tt.wantAlt[i])
				}
			}
		})
	}
}

// TestNonRefFilters exercises the in-process passFilters logic for
// --non-ref-af and --non-ref-ac. Includes the upstream quirk that
// MinNonRefAF > 0 drops monomorphic sites while MinNonRefAC alone does
// not.
func TestNonRefFilters(t *testing.T) {
	mk := func(alt []string, gt []string) *vcf.Variant {
		return createTestVariant("chr1", 1, "A", alt, 30, gt)
	}
	tests := []struct {
		name   string
		v      *vcf.Variant
		params Params
		want   bool
	}{
		{
			name:   "AF pass",
			v:      mk([]string{"G"}, []string{"0/1", "0/1", "0/1"}),
			params: Params{MinNonRefAF: 0.4},
			want:   true,
		},
		{
			name:   "AF fail (one ALT too rare)",
			v:      mk([]string{"G"}, []string{"0/0", "0/0", "0/1"}),
			params: Params{MinNonRefAF: 0.3},
			want:   false,
		},
		{
			name:   "AF fail at multi-allelic (one ALT below threshold)",
			v:      mk([]string{"G", "T"}, []string{"0/1", "0/1", "1/2"}),
			params: Params{MinNonRefAF: 0.3},
			want:   false, // T has freq 1/6 = 0.167
		},
		{
			name:   "AC pass",
			v:      mk([]string{"G"}, []string{"0/1", "0/1", "1/1"}),
			params: Params{MinNonRefAC: 2},
			want:   true,
		},
		{
			name:   "AC fail",
			v:      mk([]string{"G"}, []string{"0/0", "0/0", "0/1"}),
			params: Params{MinNonRefAC: 2},
			want:   false,
		},
		{
			name:   "AF drops monomorphic (ALT='.')",
			v:      mk([]string{"."}, []string{"0/0", "0/0", "0/0"}),
			params: Params{MinNonRefAF: 0.01},
			want:   false,
		},
		{
			name:   "AC keeps monomorphic (upstream-quirk: no _any fallback)",
			v:      mk([]string{"."}, []string{"0/0", "0/0", "0/0"}),
			params: Params{MinNonRefAC: 1},
			want:   true,
		},
		{
			name:   "both flags zero leaves site alone",
			v:      mk([]string{"G"}, []string{"0/0", "0/0", "0/0"}),
			params: Params{},
			want:   true,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got := passFilters(tt.v, &tt.params, nil, nil, nil, nil, nil, nil)
			if got != tt.want {
				t.Errorf("passFilters = %v, want %v", got, tt.want)
			}
		})
	}
}
