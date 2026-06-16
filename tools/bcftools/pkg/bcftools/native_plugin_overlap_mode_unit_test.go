// Binary-free unit tests for the overlap-mode predicate (region_target.go) and
// the gtisec arbitrary-ploidy genotype key (native_plugin_gtisec.go). These pin
// the pure helpers without any upstream binary or reference_code submodule, so
// they pass with submodules UNPOPULATED.
package bcftools

import (
	"testing"

	"github.com/yassineS/bio_ai_experiment/pkg/htsgo/vcf"
)

// TestUnitParseOverlapOption pins the MODE-string parsing (pos|0, record|1,
// variant|2, case-insensitive on the words), mirroring bcftools'
// parse_overlap_option.
func TestUnitParseOverlapOption(t *testing.T) {
	ok := []struct {
		in   string
		want int
	}{
		{"pos", overlapPos}, {"POS", overlapPos}, {"0", overlapPos},
		{"record", overlapRecord}, {"Record", overlapRecord}, {"1", overlapRecord},
		{"variant", overlapVariant}, {"VARIANT", overlapVariant}, {"2", overlapVariant},
	}
	for _, tt := range ok {
		got, err := parseOverlapOption(tt.in)
		if err != nil || got != tt.want {
			t.Errorf("parseOverlapOption(%q) = (%d,%v), want (%d,nil)", tt.in, got, err, tt.want)
		}
	}
	for _, bad := range []string{"", "3", "pos2", "rec", "x", " pos"} {
		if _, err := parseOverlapOption(bad); err == nil {
			t.Errorf("parseOverlapOption(%q) expected error", bad)
		}
	}
}

// TestUnitVariantLeadingOffset pins the off computation of htslib
// _set_variant_boundaries: the longest leading run of bases common to REF and
// every ALT, capped at rlen.
func TestUnitVariantLeadingOffset(t *testing.T) {
	tests := []struct {
		name string
		ref  string
		alts []string
		want int
	}{
		{"snp no common prefix", "A", []string{"T"}, 0},
		{"deletion shares first base", "CGT", []string{"C"}, 1},               // rlen=3, off=1
		{"insertion shares first base", "G", []string{"GAAA"}, 1},             // rlen=1, off capped to 1
		{"mnp two common", "ACGT", []string{"ACTT"}, 2},                       // common "AC"
		{"multi-alt takes shortest prefix", "ACG", []string{"ACT", "AGG"}, 1}, // min(2,1)=1
		{"alt equals ref prefix exhausted", "AC", []string{"AC"}, 2},          // off=min(rlen,2)=2
		{"empty alt", "ACGT", []string{""}, 0},
		{"no alts", "ACGT", nil, 4}, // off starts at rlen, no alt reduces it -> capped at rlen
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			rlen := len(tt.ref)
			if got := variantLeadingOffset(tt.ref, tt.alts, rlen); got != tt.want {
				t.Errorf("variantLeadingOffset(%q,%v) = %d, want %d", tt.ref, tt.alts, got, tt.want)
			}
		})
	}
}

// TestUnitOverlapBoundaries pins the [beg,end] interval each mode derives for a
// record, replicating htslib synced_bcf_reader.c.
func TestUnitOverlapBoundaries(t *testing.T) {
	// Deletion CGT>C at POS=150: rlen=3, off=1.
	del := &vcf.Variant{Chrom: "chr1", Pos: 150, Ref: "CGT", Alt: []string{"C"}}
	// Insertion G>GAAA at POS=400: rlen=1, off capped to 1 -> beg=401,end=400.
	ins := &vcf.Variant{Chrom: "chr2", Pos: 400, Ref: "G", Alt: []string{"GAAA"}}
	// SNP A>T at POS=100: rlen=1, off=0.
	snp := &vcf.Variant{Chrom: "chr1", Pos: 100, Ref: "A", Alt: []string{"T"}}

	tests := []struct {
		name     string
		v        *vcf.Variant
		mode     int
		beg, end int
	}{
		{"del pos", del, overlapPos, 150, 150},
		{"del record", del, overlapRecord, 150, 152},
		{"del variant", del, overlapVariant, 151, 152},
		{"ins pos", ins, overlapPos, 400, 400},
		{"ins record", ins, overlapRecord, 400, 400},
		{"ins variant", ins, overlapVariant, 401, 400}, // inverted (empty) span
		{"snp pos", snp, overlapPos, 100, 100},
		{"snp record", snp, overlapRecord, 100, 100},
		{"snp variant", snp, overlapVariant, 100, 100},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			beg, end := overlapBoundaries(tt.v, tt.mode)
			if beg != tt.beg || end != tt.end {
				t.Errorf("overlapBoundaries(%s, mode %d) = (%d,%d), want (%d,%d)",
					tt.name, tt.mode, beg, end, tt.beg, tt.end)
			}
		})
	}
}

// TestUnitIntervalOverlapsAny pins the unified interval-overlap predicate
// (region.start <= end && region.end >= beg), including the chrom mismatch and
// the inverted-interval (insertion variant span) edge.
func TestUnitIntervalOverlapsAny(t *testing.T) {
	regs := []region{{chrom: "chr1", beg: 150, end: 150}}
	tests := []struct {
		name     string
		chrom    string
		beg, end int
		want     bool
	}{
		{"point hits", "chr1", 150, 150, true},
		{"record span reaches", "chr1", 150, 152, true},
		{"variant span misses 150", "chr1", 151, 152, false},
		{"wrong chrom", "chr2", 150, 150, false},
		// Inverted insertion span [401,400] against a window bracketing it.
		{"inverted span bracketed", "chr1", 151, 150, true},
	}
	regsBracket := []region{{chrom: "chr1", beg: 100, end: 200}}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			use := regs
			if tt.name == "inverted span bracketed" {
				use = regsBracket
			}
			if got := intervalOverlapsAny(tt.chrom, tt.beg, tt.end, use); got != tt.want {
				t.Errorf("intervalOverlapsAny(%s,%d,%d) = %v, want %v", tt.chrom, tt.beg, tt.end, got, tt.want)
			}
		})
	}
}

// TestUnitRegionTargetParseOverlapDefaults pins that parseRegionTargetArgs seeds
// the upstream defaults (regions=record, targets=pos) and that the option
// overrides them, only when the plugin exposes the matching family.
func TestUnitRegionTargetParseOverlapDefaults(t *testing.T) {
	// No overlap option: defaults.
	_, f, err := parseRegionTargetArgs([]string{"-r", "chr1"}, allRegionTargetCaps)
	if err != nil {
		t.Fatalf("parse: %v", err)
	}
	if f.regionsOverlap != overlapRecord || f.targetsOverlap != overlapPos {
		t.Fatalf("default modes = (r %d, t %d), want (%d,%d)", f.regionsOverlap, f.targetsOverlap, overlapRecord, overlapPos)
	}
	// Explicit override — only honored for plugins that expose the overlap
	// option (overlapRegionTargetCaps); allRegionTargetCaps passes it through.
	_, f, err = parseRegionTargetArgs([]string{"-r", "chr1", "--regions-overlap", "variant", "--targets-overlap", "record"}, overlapRegionTargetCaps)
	if err != nil {
		t.Fatalf("parse: %v", err)
	}
	if f.regionsOverlap != overlapVariant || f.targetsOverlap != overlapRecord {
		t.Fatalf("override modes = (r %d, t %d), want (%d,%d)", f.regionsOverlap, f.targetsOverlap, overlapVariant, overlapRecord)
	}
	// A plugin without the targets family leaves --targets-overlap for its own
	// parser (passed through in remaining).
	rem, _, err := parseRegionTargetArgs([]string{"--targets-overlap", "record"}, regionsOnlyCaps)
	if err != nil {
		t.Fatalf("parse: %v", err)
	}
	if len(rem) != 2 || rem[0] != "--targets-overlap" {
		t.Fatalf("expected --targets-overlap passed through for regions-only caps, got %v", rem)
	}
	// An unparseable mode is an error (when the plugin exposes the option, so it
	// is consumed and validated rather than passed through).
	if _, _, err := parseRegionTargetArgs([]string{"--regions-overlap", "nope"}, overlapRegionTargetCaps); err == nil {
		t.Fatal("expected error for bad --regions-overlap value")
	}
}

// TestUnitKeepOverlapModes pins keep() end to end under each region overlap
// mode against the discriminating deletion.
func TestUnitKeepOverlapModes(t *testing.T) {
	del := &vcf.Variant{Chrom: "chr1", Pos: 150, Ref: "CGT", Alt: []string{"C"}}
	mk := func(beg, end, mode int) regionTargetFilter {
		return regionTargetFilter{
			regions:        []region{{chrom: "chr1", beg: beg, end: end}},
			hasRegions:     true,
			regionsOverlap: mode,
		}
	}
	tests := []struct {
		name     string
		f        regionTargetFilter
		wantKeep bool
	}{
		{"w150 pos keeps", mk(150, 150, overlapPos), true},
		{"w150 record keeps", mk(150, 150, overlapRecord), true},
		{"w150 variant drops", mk(150, 150, overlapVariant), false},
		{"w151 pos drops", mk(151, 151, overlapPos), false},
		{"w151 record keeps", mk(151, 151, overlapRecord), true},
		{"w151 variant keeps", mk(151, 151, overlapVariant), true},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if got := tt.f.keep(del); got != tt.wantKeep {
				t.Errorf("keep = %v, want %v", got, tt.wantKeep)
			}
		})
	}
}

// TestUnitGtisecGTKey pins the arbitrary-ploidy genotype key: same unordered
// multiset => same key (across ploidies), any missing allele => missing, and
// haploid stays distinct from a homozygous diploid.
func TestUnitGtisecGTKey(t *testing.T) {
	key := func(s string) (string, bool) {
		gt, ok := parseGT(s)
		if !ok {
			t.Fatalf("parseGT(%q) failed", s)
		}
		return gtisecGTKey(gt)
	}
	// Order-independence at ploidy 3.
	k1, m1 := key("0/1/2")
	k2, m2 := key("2/1/0")
	if m1 || m2 || k1 != k2 {
		t.Fatalf("0/1/2 and 2/1/0 must share a key: %q vs %q (missing %v,%v)", k1, k2, m1, m2)
	}
	// Distinct multisets at ploidy 3.
	if a, _ := key("0/0/1"); a == k1 {
		t.Fatalf("0/0/1 must differ from 0/1/2")
	}
	if a, _ := key("0/0/1"); a == func() string { k, _ := key("0/1/1"); return k }() {
		t.Fatalf("0/0/1 must differ from 0/1/1")
	}
	// Haploid 1 distinct from diploid 1/1.
	if h, _ := key("1"); h == func() string { k, _ := key("1/1"); return k }() {
		t.Fatal("haploid 1 must differ from diploid 1/1")
	}
	// Tetraploid order-independence.
	t1, _ := key("0/0/1/1")
	t2, _ := key("1/0/1/0")
	if t1 != t2 {
		t.Fatalf("0/0/1/1 and 1/0/1/0 must share a key: %q vs %q", t1, t2)
	}
	// Any missing allele => missing.
	for _, s := range []string{"./.", "0/.", "./1/2", "0/1/."} {
		if _, missing := key(s); !missing {
			t.Errorf("%q should be treated as missing", s)
		}
	}
	// Fully present => not missing.
	if _, missing := key("0/1/2"); missing {
		t.Error("0/1/2 should not be missing")
	}
}
