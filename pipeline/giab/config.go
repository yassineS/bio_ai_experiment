package giab

import (
	"encoding/json"
	"fmt"
	"os"
	"sort"
)

// DocPointer is the canonical pointer printed alongside every SKIP so an
// operator knows where to find the data-acquisition recipe.
const DocPointer = "see docs/GIAB_CONCORDANCE.md"

// Stratification names a GA4GH/GIAB stratification BED used to break the
// hap.py/vcfeval comparison down by region difficulty (e.g. CMRG, low-mappability,
// segmental-duplications, alldifficultregions).
type Stratification struct {
	Name string `json:"name"`
	BED  string `json:"bed"`
}

// Config is the full description of one concordance run. All paths are absolute
// or relative to the working directory. Any field may be empty; the harness
// checks each prerequisite and SKIPs the dependent stage with a clear reason
// rather than failing, so a partially populated config still runs the parts it
// can (e.g. the ours-vs-upstream comparison without a truth set).
type Config struct {
	// Sample is the GIAB sample label (HG001..HG007 / NA12878 etc.).
	Sample string `json:"sample"`
	// Build is the reference build label (GRCh37 / GRCh38 / T2T), informational.
	Build string `json:"build"`

	// Reference is the indexed reference FASTA (.fa with a .fai).
	Reference string `json:"reference"`
	// ReadsBAM is the aligned-reads BAM the call sets are produced from.
	ReadsBAM string `json:"reads_bam"`

	// TruthVCF is the GIAB benchmark truth VCF (v4.2.1).
	TruthVCF string `json:"truth_vcf"`
	// HighConfBED is the GIAB high-confidence region BED.
	HighConfBED string `json:"high_conf_bed"`

	// Stratifications is the optional list of GA4GH stratification BEDs.
	Stratifications []Stratification `json:"stratifications,omitempty"`

	// OurBcftools and UpstreamBcftools are the two bcftools binaries to drive.
	// If empty, the harness resolves OUR binary by building it and UPSTREAM via
	// pipeline/internal/upstream.
	OurBcftools      string `json:"our_bcftools,omitempty"`
	UpstreamBcftools string `json:"upstream_bcftools,omitempty"`

	// HappyBin / VcfevalBin point at the benchmarking engines. If both are
	// empty the harness probes PATH (hap.py, then rtg/vcfeval); if neither is
	// found the biological-concordance stage SKIPs.
	HappyBin   string `json:"happy_bin,omitempty"`
	VcfevalBin string `json:"vcfeval_bin,omitempty"`
	// SDFTemplate is the RTG SDF (reference) directory vcfeval requires. Only
	// used when the engine is vcfeval.
	SDFTemplate string `json:"sdf_template,omitempty"`

	// QualULP overrides the ULP/Phred tolerance (defaults to DefaultQualULP).
	QualULP float64 `json:"qual_ulp,omitempty"`

	// OutDir is where giab_concordance.{json,md} are written.
	OutDir string `json:"out_dir,omitempty"`

	// MaxDiffs caps the number of differing sites embedded in the report
	// (0 = the package default; negative = unlimited).
	MaxDiffs int `json:"max_diffs,omitempty"`
}

// LoadConfig reads a JSON config file.
func LoadConfig(path string) (*Config, error) {
	b, err := os.ReadFile(path)
	if err != nil {
		return nil, err
	}
	var c Config
	if err := json.Unmarshal(b, &c); err != nil {
		return nil, fmt.Errorf("parsing config %s: %w", path, err)
	}
	return &c, nil
}

// effectiveQualULP returns the configured ULP tolerance or the package default.
func (c *Config) effectiveQualULP() float64 {
	if c.QualULP > 0 {
		return c.QualULP
	}
	return DefaultQualULP
}

// effectiveMaxDiffs returns the configured diff cap or the default (200).
func (c *Config) effectiveMaxDiffs() int {
	if c.MaxDiffs < 0 {
		return 0 // unlimited sentinel handled by caller
	}
	if c.MaxDiffs == 0 {
		return 200
	}
	return c.MaxDiffs
}

// missingFile reports whether a path is empty or does not name an existing
// regular file/dir. Empty path returns true (treated as "not provided").
func missingFile(path string) bool {
	if path == "" {
		return true
	}
	_, err := os.Stat(path)
	return err != nil
}

// sortedStratNames returns the stratification names in stable order, for
// deterministic reports.
func (c *Config) sortedStratNames() []string {
	names := make([]string, 0, len(c.Stratifications))
	for _, s := range c.Stratifications {
		names = append(names, s.Name)
	}
	sort.Strings(names)
	return names
}
