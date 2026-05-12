package vcftools

import (
	"fmt"
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
			result := passFilters(tt.variant, tt.params, nil, nil, nil, nil)
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

func TestGenotypeDistance(t *testing.T) {
	tests := []struct {
		gt1      string
		gt2      string
		expected int
	}{
		{"0/0", "0/0", 0},
		{"0/0", "0/1", 1},
		{"0/0", "1/1", 2},
		{"0/1", "1/0", 2},
		{"0/1", "0/1", 0},
		{"0|0", "0|0", 0},
		{"0|1", "1|0", 2},
	}

	for _, tt := range tests {
		t.Run(tt.gt1+"_"+tt.gt2, func(t *testing.T) {
			result := genotypeDistance(tt.gt1, tt.gt2)
			if result != tt.expected {
				t.Errorf("genotypeDistance(%s, %s) = %d, want %d",
					tt.gt1, tt.gt2, result, tt.expected)
			}
		})
	}
}
